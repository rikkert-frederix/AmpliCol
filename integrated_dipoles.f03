module integrated_dipoles
  ! Inventory of integrated subtraction histories.  The inventory is built
  ! from the local dipoles themselves; it never independently enumerates
  ! splittings.  This is the central consistency guarantee between D and
  ! I/P/K at leading colour.
  use handling_processes
  use cs_dipole_mappings, only: cs_dipole_topology
  use cs_integrated_kernels
  use cs_massive_integrated_kernels
  use pdf_wrap, only: evaluate_pdf_flavour
  use common, only: alpha_dipole
  implicit none
  private

  type, public :: integrated_history
     integer :: real_group=0
     integer :: real_process=0
     integer :: real_copy=0
     integer :: local_dipole=0
     integer :: topology=0
     integer :: born_emitter=0
     integer :: born_spectator=0
     integer :: incoming_leg=0
     integer :: real_incoming_flavour=0
     integer :: born_incoming_flavour=0
     real(kind=8) :: emitter_mass=0d0
     real(kind=8) :: unresolved_mass=0d0
     real(kind=8) :: parent_mass=0d0
     real(kind=8) :: spectator_mass=0d0
     integer, allocatable :: real_flavours(:)
     integer, allocatable :: real_colour_order(:)
     integer, allocatable :: born_flavours(:)
     integer, allocatable :: born_colour_order(:)
  end type integrated_history

  type(integrated_history), allocatable, save, public :: integrated_history_list(:)
  integer, save, public :: n_integrated_histories=0
  integer, save, public :: integrated_dimensional_scheme=cs_scheme_hv
  integer, save, public :: integrated_n_active_flavours=5

  public :: initialise_integrated_dipoles, history_matches_born
  public :: integrated_endpoint, integrated_beam

contains

  subroutine initialise_integrated_dipoles(nborn_groups,scheme_name,n_active_flavours)
    integer, intent(in) :: nborn_groups
    integer, intent(in) :: n_active_flavours
    character(len=*), intent(in) :: scheme_name
    type(integrated_history), allocatable :: candidates(:),unique_histories(:)
    integer :: igroup,iproc,icopy,idip,max_histories,ncandidate,i
    logical :: matched

    if (trim(scheme_name).eq.'hv') then
       integrated_dimensional_scheme=cs_scheme_hv
    elseif (trim(scheme_name).eq.'fdh') then
       integrated_dimensional_scheme=cs_scheme_fdh
    else
       write(*,*) 'ERROR: unsupported integrated-dipole dimensional scheme: ',trim(scheme_name)
       stop 1
    endif
    if (n_active_flavours.lt.1 .or. n_active_flavours.gt.5) then
       write(*,*) 'ERROR: integrated subtraction requires 1 <= nf <= 5:',n_active_flavours
       stop 1
    endif
    integrated_n_active_flavours=n_active_flavours

    max_histories=0
    do igroup=nborn_groups+1,ngroups
       if (.not.allocated(pgl(igroup)%dpl)) cycle
       do iproc=1,pgl(igroup)%nproc
          max_histories=max_histories+pgl(igroup)%iden_iproc(iproc)*pgl(igroup)%dpl(iproc)%ndip
       enddo
    enddo
    if (max_histories.eq.0) then
       write(*,*) 'ERROR: no local dipoles are available for integrated subtraction'
       stop 1
    endif

    allocate(candidates(max_histories))
    ncandidate=0
    do igroup=nborn_groups+1,ngroups
       if (.not.pgl(igroup)%is_subtracted_real) cycle
       do iproc=1,pgl(igroup)%nproc
          do icopy=1,pgl(igroup)%iden_iproc(iproc)
             do idip=1,pgl(igroup)%dpl(iproc)%ndip
                call ensure_supported_integrated_history(igroup,iproc,idip)
                call fill_history(candidates(ncandidate+1),igroup,iproc,icopy,idip)
                call canonicalize_history_born(candidates(ncandidate+1),nborn_groups,matched)
                if (.not.matched) then
                   write(*,*) 'ERROR: integrated dipole has no corresponding Born process/order'
                   write(*,*) ' real group/process/copy/dipole:',igroup,iproc,icopy,idip
                   write(*,*) ' reduced process:',candidates(ncandidate+1)%born_flavours
                   write(*,*) ' reduced colour order:',candidates(ncandidate+1)%born_colour_order
                   call print_born_flavour_matches(candidates(ncandidate+1)%born_flavours,nborn_groups)
                   stop 1
                endif
                if (is_duplicate_history(candidates(ncandidate+1),candidates,ncandidate)) cycle
                ncandidate=ncandidate+1
             enddo
          enddo
       enddo
    enddo

    allocate(unique_histories(ncandidate))
    do i=1,ncandidate
       unique_histories(i)=candidates(i)
    enddo
    if (allocated(integrated_history_list)) deallocate(integrated_history_list)
    call move_alloc(unique_histories,integrated_history_list)
    deallocate(candidates)
    n_integrated_histories=ncandidate
    call validate_integrated_history_poles(nborn_groups)

    write(*,'(a,i0,a,a)') 'Integrated subtraction registry: ',n_integrated_histories,&
         ' locally matched histories, scheme=',trim(scheme_name)
    write(99,'(a,i0,a,a)') 'Integrated subtraction registry: ',n_integrated_histories,&
         ' locally matched histories, scheme=',trim(scheme_name)
  end subroutine initialise_integrated_dipoles

  subroutine validate_integrated_history_poles(nborn_groups)
    ! Every massless Born emitter must reconstruct its universal I-operator
    ! poles from physical-parent histories after summing all physical
    ! flavours, ordered splittings, and leading-colour neighbours.  The
    ! auxiliary U(1) histories are colour-flow corrections rather than an
    ! additional physical gluon splitting and are checked through their
    ! explicit normalization in history_i_primitive.
    integer, intent(in) :: nborn_groups
    real(kind=8), parameter :: tolerance=1d-10
    real(kind=8) :: primitive(-2:0),pole_sum(-2:-1),expected(-2:-1),weight
    integer :: igroup,iproc,icopy,emitter,ih,fi,fj,fp,parton,nf_available
    integer :: mapped_emitter,mapped_spectator
    integer :: leg_map(maxval(pgl(1:nborn_groups)%next))
    logical :: invalid,massive_emitter_histories

    do igroup=1,nborn_groups
       do iproc=1,pgl(igroup)%nproc
          do icopy=1,pgl(igroup)%iden_iproc(iproc)
             do emitter=1,pgl(igroup)%next
                fp=pgl(igroup)%iden_processes(emitter,icopy,iproc)
                call parton_kind(fp,parton)
                if (parton.eq.0) cycle
                pole_sum=0d0
                massive_emitter_histories=.false.
                do ih=1,n_integrated_histories
                   if (.not.history_matches_copy(integrated_history_list(ih),igroup,iproc,icopy,&
                        leg_map(1:pgl(igroup)%next))) cycle
                   mapped_emitter=leg_map(integrated_history_list(ih)%born_emitter)
                   if (mapped_emitter.ne.emitter) cycle
                   if (history_has_mass(integrated_history_list(ih))) then
                      massive_emitter_histories=.true.
                      cycle
                   endif
                   fi=integrated_history_list(ih)%real_flavours(&
                        pgl(integrated_history_list(ih)%real_group)%dpl(&
                        integrated_history_list(ih)%real_process)%dl(&
                        integrated_history_list(ih)%local_dipole)%dip_ijk(1))
                   fj=integrated_history_list(ih)%real_flavours(&
                        pgl(integrated_history_list(ih)%real_group)%dpl(&
                        integrated_history_list(ih)%real_process)%dl(&
                        integrated_history_list(ih)%local_dipole)%dip_ijk(2))
                   fp=integrated_history_list(ih)%born_flavours(emitter)
                   if (fp.eq.99) cycle
                   call history_i_primitive(fi,fj,fp,integrated_history_list(ih)%topology,&
                        primitive,parton)
                   weight=pgl(integrated_history_list(ih)%real_group)%dpl(&
                        integrated_history_list(ih)%real_process)%dl(&
                        integrated_history_list(ih)%local_dipole)%lc_weight
                   pole_sum=pole_sum+weight*primitive(-2:-1)
                enddo
                if (massive_emitter_histories) cycle
                call parton_kind(pgl(igroup)%iden_processes(emitter,icopy,iproc),parton)
                if (parton.eq.cs_parton_q) then
                   expected=[cs_cf_lc,1.5d0*cs_cf_lc]
                   invalid=any(abs(pole_sum-expected).gt.tolerance*max(1d0,abs(expected)))
                else
                   ! A process list may deliberately omit real states with
                   ! three or more quark lines.  Require the universal gg
                   ! poles and an integer number (zero through five) of
                   ! complete massless q-qbar flavour sectors.
                   nf_available=nint(((11d0/6d0)*cs_ca-pole_sum(-1))/&
                        ((2d0/3d0)*cs_tr))
                   expected=[cs_ca,(11d0/6d0)*cs_ca-&
                        (2d0/3d0)*cs_tr*dble(nf_available)]
                   invalid=nf_available.lt.0 .or. &
                        nf_available.gt.integrated_n_active_flavours .or. &
                        any(abs(pole_sum-expected).gt.tolerance*max(1d0,abs(expected)))
                endif
                if (invalid) then
                   write(*,*) 'ERROR: integrated histories do not reconstruct universal poles'
                   write(*,*) ' Born group/process/copy/emitter:',igroup,iproc,icopy,emitter
                   write(*,*) ' reconstructed double/single poles:',pole_sum
                   write(*,*) ' expected double/single poles:',expected
                   do ih=1,n_integrated_histories
                      if (.not.history_matches_copy(integrated_history_list(ih),igroup,iproc,icopy,&
                           leg_map(1:pgl(igroup)%next))) cycle
                      mapped_emitter=leg_map(integrated_history_list(ih)%born_emitter)
                      if (mapped_emitter.ne.emitter) cycle
                      mapped_spectator=leg_map(integrated_history_list(ih)%born_spectator)
                      fi=integrated_history_list(ih)%real_flavours(&
                           pgl(integrated_history_list(ih)%real_group)%dpl(&
                           integrated_history_list(ih)%real_process)%dl(&
                           integrated_history_list(ih)%local_dipole)%dip_ijk(1))
                      fj=integrated_history_list(ih)%real_flavours(&
                           pgl(integrated_history_list(ih)%real_group)%dpl(&
                           integrated_history_list(ih)%real_process)%dl(&
                           integrated_history_list(ih)%local_dipole)%dip_ijk(2))
                      fp=integrated_history_list(ih)%born_flavours(emitter)
                      if (fp.eq.99) cycle
                      call history_i_primitive(fi,fj,fp,integrated_history_list(ih)%topology,&
                           primitive,parton)
                      weight=pgl(integrated_history_list(ih)%real_group)%dpl(&
                           integrated_history_list(ih)%real_process)%dl(&
                           integrated_history_list(ih)%local_dipole)%lc_weight
                      write(*,*) ' history/topology/fi/fj/fp/emitter/spectator/weight/poles:',ih,&
                           integrated_history_list(ih)%topology,fi,fj,fp,mapped_emitter,&
                           mapped_spectator,weight,primitive(-2:-1)
                   enddo
                   stop 1
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine validate_integrated_history_poles

  subroutine fill_history(history,igroup,iproc,icopy,idip)
    type(integrated_history), intent(out) :: history
    integer, intent(in) :: igroup,iproc,icopy,idip
    type(dipole) :: dip
    integer :: emitter,spectator

    dip=pgl(igroup)%dpl(iproc)%dl(idip)
    history%real_group=igroup
    history%real_process=iproc
    history%real_copy=icopy
    history%local_dipole=idip
    history%topology=cs_dipole_topology(dip%dip_ijk)
    history%born_emitter=dip%dip_r_ijk(1)
    history%born_spectator=dip%dip_r_ijk(2)
    history%real_flavours=pgl(igroup)%iden_processes(:,icopy,iproc)
    history%real_colour_order=pgl(igroup)%color_orders(:,iproc)
    history%emitter_mass=abs(phys_model%get_mass(dip%dip_ijk_f(1)))
    history%unresolved_mass=abs(phys_model%get_mass(dip%dip_ijk_f(2)))
    history%parent_mass=abs(phys_model%get_mass(dip%dip_r_ijk_f(1)))
    history%spectator_mass=abs(phys_model%get_mass(dip%dip_r_ijk_f(2)))
    allocate(history%born_flavours(size(dip%process_r)))
    call make_history_reduced_process(pgl(igroup)%processes(:,iproc),history%real_flavours,&
         dip%process_r,history%born_flavours)
    history%born_colour_order=dip%reduced_color_order
    emitter=dip%dip_ijk(1)
    if (emitter.le.2) then
       history%incoming_leg=emitter
       history%real_incoming_flavour=history%real_flavours(emitter)
       history%born_incoming_flavour=history%born_flavours(dip%dip_r_ijk(1))
    elseif (dip%dip_ijk(3).le.2) then
       ! A final--initial dipole changes the momentum fraction of its
       ! initial-state spectator, so its non-endpoint alpha contribution
       ! belongs to that beam even though the emitter is final state.
       spectator=dip%dip_ijk(3)
       history%incoming_leg=spectator
       history%real_incoming_flavour=history%real_flavours(spectator)
       history%born_incoming_flavour=history%born_flavours(dip%dip_r_ijk(2))
    endif
  end subroutine fill_history

  pure subroutine make_history_reduced_process(base_real,copy_real,base_reduced,copy_reduced)
    integer, intent(in) :: base_real(:),copy_real(:),base_reduced(:)
    integer, intent(out) :: copy_reduced(:)
    integer :: flavour_map(6)
    integer :: i,base_flavour,copy_flavour

    flavour_map=0
    do i=1,size(base_real)
       if (abs(base_real(i)).lt.1 .or. abs(base_real(i)).gt.6) cycle
       base_flavour=abs(base_real(i))
       copy_flavour=abs(copy_real(i))
       if (copy_flavour.lt.1 .or. copy_flavour.gt.6) cycle
       if (flavour_map(base_flavour).eq.0) flavour_map(base_flavour)=copy_flavour
    enddo

    copy_reduced=base_reduced
    do i=1,size(copy_reduced)
       if (abs(copy_reduced(i)).lt.1 .or. abs(copy_reduced(i)).gt.6) cycle
       base_flavour=abs(copy_reduced(i))
       if (flavour_map(base_flavour).eq.0) cycle
       copy_reduced(i)=sign(flavour_map(base_flavour),copy_reduced(i))
    enddo
  end subroutine make_history_reduced_process

  subroutine ensure_supported_integrated_history(igroup,iproc,idip)
    ! The unresolved leg and both incoming legs remain massless.  A massive
    ! emitter is supported only for final-state Q->Qg, while a massive
    ! spectator is supported in FF and IF histories.
    integer, intent(in) :: igroup,iproc,idip
    type(dipole) :: dip
    real(kind=8) :: emitter_mass,unresolved_mass,parent_mass,spectator_mass,mass_tolerance
    integer :: topology,fi,fj,fp
    logical :: supported

    dip=pgl(igroup)%dpl(iproc)%dl(idip)
    fi=dip%dip_ijk_f(1)
    fj=dip%dip_ijk_f(2)
    fp=dip%dip_r_ijk_f(1)
    emitter_mass=abs(phys_model%get_mass(dip%dip_ijk_f(1)))
    unresolved_mass=abs(phys_model%get_mass(dip%dip_ijk_f(2)))
    parent_mass=abs(phys_model%get_mass(dip%dip_r_ijk_f(1)))
    spectator_mass=abs(phys_model%get_mass(dip%dip_r_ijk_f(2)))
    mass_tolerance=100d0*epsilon(1d0)*max(1d0,emitter_mass,unresolved_mass,parent_mass,spectator_mass)
    if (emitter_mass.le.mass_tolerance .and. unresolved_mass.le.mass_tolerance .and. &
         parent_mass.le.mass_tolerance .and. spectator_mass.le.mass_tolerance) return

    topology=cs_dipole_topology(dip%dip_ijk)
    supported=.true.
    if (unresolved_mass.gt.mass_tolerance) supported=.false.
    if (emitter_mass.gt.mass_tolerance .or. parent_mass.gt.mass_tolerance) then
       if (topology.ne.1 .and. topology.ne.2) supported=.false.
       if (abs(fi).gt.6 .or. abs(fp).gt.6 .or. abs(fi).ne.abs(fp)) supported=.false.
       if (fj.ne.21 .and. fj.ne.99) supported=.false.
       if (abs(emitter_mass-parent_mass).gt.mass_tolerance) supported=.false.
    endif
    if (spectator_mass.gt.mass_tolerance .and. topology.ne.1 .and. topology.ne.3) supported=.false.
    if (supported) return

    write(*,*) 'ERROR: unsupported massive integrated dipole history'
    write(*,*) ' real group/process/dipole/topology:',igroup,iproc,idip,topology
    write(*,*) ' real emitter/unresolved masses:',emitter_mass,unresolved_mass
    write(*,*) ' reduced parent/spectator flavours:',dip%dip_r_ijk_f
    write(*,*) ' reduced parent/spectator masses:',parent_mass,spectator_mass
    write(99,*) 'ERROR: unsupported massive integrated dipole history'
    write(99,*) ' real group/process/dipole/topology:',igroup,iproc,idip,topology
    write(99,*) ' real emitter/unresolved masses:',emitter_mass,unresolved_mass
    write(99,*) ' reduced parent/spectator flavours:',dip%dip_r_ijk_f
    write(99,*) ' reduced parent/spectator masses:',parent_mass,spectator_mass
    stop 1
  end subroutine ensure_supported_integrated_history

  subroutine canonicalize_history_born(history,nborn_groups,matched)
    ! Reduced histories can place a collapsed U(1)/gluon parent in a
    ! different array slot from the canonical physical Born process.  Find
    ! the closest equivalent Born layout once and translate all reduced-leg
    ! indices to it.  Subsequent matching is deliberately strict so distinct
    ! colour-ordered partial amplitudes are never conflated.
    type(integrated_history), intent(inout) :: history
    integer, intent(in) :: nborn_groups
    logical, intent(out) :: matched
    integer :: igroup,iproc,ip,i,score,best_score
    integer :: leg_map(size(history%born_flavours))
    integer :: best_map(size(history%born_flavours))
    integer :: best_process(size(history%born_flavours))
    integer :: best_order(size(history%born_colour_order))
    integer :: old_emitter,old_spectator,old_incoming,parent_flavour
    logical :: trial_matched

    matched=.false.
    best_score=-1
    best_map=0
    do igroup=1,nborn_groups
       if (pgl(igroup)%next.ne.size(history%born_flavours)) cycle
       do iproc=1,pgl(igroup)%nproc
          do ip=1,pgl(igroup)%iden_iproc(iproc)
             call match_coloured_legs(history%born_flavours,history%born_colour_order,&
                  pgl(igroup)%iden_processes(:,ip,iproc),pgl(igroup)%color_orders(:,iproc),&
                  leg_map,trial_matched)
             if (.not.trial_matched) cycle
             score=count(leg_map.eq.[(i,i=1,size(leg_map))])
             if (score.gt.best_score) then
                best_score=score
                best_map=leg_map
                best_process=pgl(igroup)%iden_processes(:,ip,iproc)
                best_order=pgl(igroup)%color_orders(:,iproc)
             endif
          enddo
       enddo
    enddo
    if (best_score.lt.0) return

    old_emitter=history%born_emitter
    old_spectator=history%born_spectator
    old_incoming=history%incoming_leg
    parent_flavour=history%born_flavours(old_emitter)
    history%born_flavours=best_process
    history%born_emitter=best_map(old_emitter)
    history%born_spectator=best_map(old_spectator)
    history%born_flavours(history%born_emitter)=parent_flavour
    history%born_colour_order=best_order
    if (old_incoming.gt.0) then
       history%incoming_leg=best_map(old_incoming)
       history%born_incoming_flavour=history%born_flavours(history%incoming_leg)
    endif
    matched=.true.
  end subroutine canonicalize_history_born

  subroutine print_born_flavour_matches(process,nborn_groups)
    integer, intent(in) :: process(:),nborn_groups
    integer :: igroup,iproc,ip
    do igroup=1,nborn_groups
       do iproc=1,pgl(igroup)%nproc
          do ip=1,pgl(igroup)%iden_iproc(iproc)
             if (.not.same_flavour_pattern(process,pgl(igroup)%iden_processes(:,ip,iproc))) cycle
             write(*,*) ' flavour-compatible Born group/process/copy:',igroup,iproc,ip
             write(*,*) ' Born process:',pgl(igroup)%iden_processes(:,ip,iproc)
             write(*,*) ' Born colour order:',pgl(igroup)%color_orders(:,iproc)
          enddo
       enddo
    enddo
  end subroutine print_born_flavour_matches

  logical function history_matches_born(history,igroup,iproc)
    type(integrated_history), intent(in) :: history
    integer, intent(in) :: igroup,iproc
    integer :: ip
    history_matches_born=.false.
    if (pgl(igroup)%next.ne.size(history%born_flavours)) return
    do ip=1,pgl(igroup)%iden_iproc(iproc)
       if (.not.same_flavour_vector(history%born_flavours,&
            pgl(igroup)%iden_processes(:,ip,iproc))) cycle
       if (.not.same_coloured_order(history%born_flavours,history%born_colour_order,&
            pgl(igroup)%iden_processes(:,ip,iproc),pgl(igroup)%color_orders(:,iproc))) cycle
       history_matches_born=.true.
       return
    enddo
  end function history_matches_born

  logical function history_matches_copy(history,igroup,iproc,icopy,leg_map)
    type(integrated_history), intent(in) :: history
    integer, intent(in) :: igroup,iproc,icopy
    integer, intent(out), optional :: leg_map(:)
    integer :: i

    history_matches_copy=.false.
    if (pgl(igroup)%next.ne.size(history%born_flavours)) return
    if (present(leg_map)) then
       if (size(leg_map).ne.size(history%born_flavours)) then
          write(*,*) 'ERROR: integrated-history leg map has incompatible size'
          stop 1
       endif
    endif
    if (.not.same_flavour_vector(history%born_flavours,&
         pgl(igroup)%iden_processes(:,icopy,iproc))) return
    if (.not.same_coloured_order(history%born_flavours,history%born_colour_order,&
         pgl(igroup)%iden_processes(:,icopy,iproc),pgl(igroup)%color_orders(:,iproc))) return
    history_matches_copy=.true.
    if (present(leg_map)) leg_map=[(i,i=1,size(leg_map))]
  end function history_matches_copy

  subroutine integrated_endpoint(igroup,iproc,born_copy,p,mu_ren,alpha_s,coeff,coeff_copy)
    ! Laurent coefficients of I(eps), multiplied by the corresponding
    ! PDF-weighted Born contributions.  Every entry is inherited from an
    ! actual local-dipole history.
    integer, intent(in) :: igroup,iproc
    real(kind=8), intent(in) :: born_copy(:),p(0:,:)
    real(kind=8), intent(in) :: mu_ren,alpha_s
    real(kind=8), intent(out) :: coeff(-2:0)
    real(kind=8), intent(out),optional :: coeff_copy(-2:,:)
    real(kind=8) :: primitive(-2:0),expanded(-2:0),sij,ell,weight
    real(kind=8) :: endpoint_alpha
    real(kind=8) :: copy_contribution(-2:0)
    integer :: ih,icopy,fi,fj,fp,parton,info
    integer :: emitter,spectator,leg_map(pgl(igroup)%next)
    logical :: shifted(size(born_copy),pgl(igroup)%next)
    logical :: massive_history

    coeff=0d0
    if (present(coeff_copy)) then
       if (size(coeff_copy,2).ne.size(born_copy)) then
          write(*,*) 'ERROR: integrated endpoint copy array has incompatible size'
          stop 1
       endif
       coeff_copy=0d0
    endif
    shifted=.false.
    do icopy=1,size(born_copy)
       if (born_copy(icopy).eq.0d0) cycle
       do ih=1,n_integrated_histories
          if (.not.history_matches_copy(integrated_history_list(ih),igroup,iproc,icopy,leg_map)) cycle
          fi=integrated_history_list(ih)%real_flavours(&
               pgl(integrated_history_list(ih)%real_group)%dpl(&
               integrated_history_list(ih)%real_process)%dl(&
               integrated_history_list(ih)%local_dipole)%dip_ijk(1))
          fj=integrated_history_list(ih)%real_flavours(&
               pgl(integrated_history_list(ih)%real_group)%dpl(&
               integrated_history_list(ih)%real_process)%dl(&
               integrated_history_list(ih)%local_dipole)%dip_ijk(2))
          fp=integrated_history_list(ih)%born_flavours(integrated_history_list(ih)%born_emitter)
          emitter=leg_map(integrated_history_list(ih)%born_emitter)
          spectator=leg_map(integrated_history_list(ih)%born_spectator)
          massive_history=history_has_mass(integrated_history_list(ih))
          if (massive_history) then
             call massive_history_endpoint(integrated_history_list(ih),fi,fj,fp,p,&
                  emitter,spectator,mu_ren,expanded,parton,info)
             if (info.ne.0) then
                write(*,*) 'ERROR: massive integrated endpoint failed'
                write(*,*) ' history/topology/fi/fj/fp/status:',ih,&
                     integrated_history_list(ih)%topology,fi,fj,fp,info
                stop 1
             endif
             if (parton.eq.0) cycle
          else
             call history_i_primitive(fi,fj,fp,integrated_history_list(ih)%topology,&
                  primitive,parton)
             if (parton.eq.0) cycle
          endif
          sij=abs(2d0*minkowski_dot(p(:,emitter),p(:,spectator)))
          if (sij.le.tiny(1d0)) cycle
          if (.not.massive_history) then
             ell=log(4d0*cs_pi*mu_ren*mu_ren/sij)-0.5772156649015328606d0
             expanded(-2)=primitive(-2)
             expanded(-1)=primitive(-1)+ell*primitive(-2)
             expanded(0)=primitive(0)+ell*primitive(-1)+&
                  (0.5d0*ell*ell-cs_pi**2/12d0)*primitive(-2)
          endif
          weight=pgl(integrated_history_list(ih)%real_group)%dpl(&
               integrated_history_list(ih)%real_process)%dl(&
               integrated_history_list(ih)%local_dipole)%lc_weight
          if (.not.massive_history) then
             endpoint_alpha=0d0
             select case (integrated_history_list(ih)%topology)
             case (1)
                call cs_ff_alpha_endpoint(alpha_dipole(1),primitive,endpoint_alpha,info)
             case (2)
                ! The FI delta baseline differs from the FF insertion by the
                ! negative single-pole coefficient.
                expanded(0)=expanded(0)-primitive(-1)
                call cs_fi_alpha_endpoint(alpha_dipole(2),primitive,endpoint_alpha,info)
             case (3,4)
                info=0
             case default
                info=-4
             end select
             if (info.ne.0) then
                write(*,*) 'ERROR: invalid integrated endpoint topology/alpha:',&
                     integrated_history_list(ih)%topology,info
                stop 1
             endif
             expanded(0)=expanded(0)+endpoint_alpha
          endif
          copy_contribution=born_copy(icopy)*weight*expanded*alpha_s/(2d0*cs_pi)
          coeff=coeff+copy_contribution
          if (present(coeff_copy)) coeff_copy(:,icopy)=coeff_copy(:,icopy)+copy_contribution
          if (.not.massive_history .and. integrated_dimensional_scheme.eq.cs_scheme_fdh .and. &
               phys_model%get_mass(pgl(igroup)%iden_processes(emitter,icopy,iproc)).eq.0d0 .and. &
               .not.shifted(icopy,emitter)) then
             copy_contribution=0d0
             copy_contribution(0)=born_copy(icopy)*cs_fdh_endpoint_shift(parton)*alpha_s/(2d0*cs_pi)
             coeff=coeff+copy_contribution
             if (present(coeff_copy)) coeff_copy(:,icopy)=coeff_copy(:,icopy)+copy_contribution
             shifted(icopy,emitter)=.true.
          endif
       enddo
    enddo
  end subroutine integrated_endpoint

  subroutine massive_history_endpoint(history,fi,fj,fp,p,emitter,spectator,mu_ren,coeff,parton,info)
    type(integrated_history), intent(in) :: history
    integer, intent(in) :: fi,fj,fp,emitter,spectator
    real(kind=8), intent(in) :: p(0:,:),mu_ren
    real(kind=8), intent(out) :: coeff(-2:0)
    integer, intent(out) :: parton,info
    real(kind=8) :: sij,q2,ell
    integer :: split,pa,pb

    coeff=0d0
    parton=0
    info=0
    sij=abs(2d0*minkowski_dot(p(:,emitter),p(:,spectator)))
    if (sij.le.tiny(1d0)) then
       info=-10
       return
    endif
    call parton_kind(fp,parton)

    select case (history%topology)
    case (1)
       split=0
       if (parton.eq.cs_parton_q .and. (fj.eq.21 .or. fj.eq.99)) then
          split=cs_massive_split_qg
       elseif (parton.eq.cs_parton_g .and. (fi.eq.21 .or. fi.eq.99) .and. &
            (fj.eq.21 .or. fj.eq.99)) then
          split=cs_massive_split_gg
       elseif (parton.eq.cs_parton_g .and. abs(fi).le.6 .and. fi.eq.-fj) then
          split=cs_massive_split_qqbar
       endif
       if (split.eq.0) then
          info=-11
          return
       endif
       q2=history%parent_mass**2+history%spectator_mass**2+sij
       ! The explicit massive final--final formulae use Q_ik^2 as their
       ! dimensional scale, rather than 2 p_i.p_k.
       ell=log(4d0*cs_pi*mu_ren*mu_ren/q2)-0.5772156649015328606d0
       call cs_massive_ff_endpoint(parton,split,history%parent_mass,&
            history%spectator_mass,q2,ell,alpha_dipole(1),&
            integrated_dimensional_scheme,coeff,info)
       if (info.eq.0 .and. fp.eq.99 .and. split.eq.cs_massive_split_qqbar) &
            coeff=coeff/(cs_ca*cs_ca)
    case (2)
       if (parton.ne.cs_parton_q .or. history%parent_mass.le.0d0) then
          info=-12
          return
       endif
       ell=log(4d0*cs_pi*mu_ren*mu_ren/sij)-0.5772156649015328606d0
       call cs_massive_fi_endpoint(history%parent_mass,sij,ell,&
            alpha_dipole(2),coeff,info)
    case (3)
       call parton_kind(fi,pa)
       call parton_kind(fp,pb)
       if (pa.eq.0 .or. pb.eq.0 .or. history%spectator_mass.le.0d0) then
          info=-13
          return
       endif
       ell=log(4d0*cs_pi*mu_ren*mu_ren/sij)-0.5772156649015328606d0
       call cs_massive_if_endpoint(pa,pb,history%spectator_mass,sij,ell,&
            integrated_dimensional_scheme,integrated_n_active_flavours,coeff,info)
    case default
       info=-14
    end select
  end subroutine massive_history_endpoint

  subroutine integrated_beam(igroup,iproc,beam,z,hard_copy,xbj,mu_ren,mu_fac,alpha_s,&
       pterm,kterm,pterm_copy,kterm_copy)
    integer, intent(in) :: igroup,iproc,beam
    real(kind=8), intent(in) :: z,hard_copy(:),xbj(2),mu_ren,mu_fac,alpha_s
    real(kind=8), intent(out) :: pterm,kterm
    real(kind=8), intent(out),optional :: pterm_copy(:),kterm_copy(:)
    type(cs_distribution) :: pk,kk,tilde_kernel,alpha_kernel
    type(cs_convolution_kernel) :: massive_kernel
    real(kind=8) :: fa,fb,gz,g1,other_pdf,regularised_p,regularised_k
    real(kind=8) :: p_delta,sij,colour_log,colour_log_endpoint
    real(kind=8) :: history_weight,kernel_factor,primitive(-2:0)
    real(kind=8) :: fi_regular,fi_subtracted,szone,sx
    integer :: ih,icopy,a,b,other,info,pa,pb,fi,fj,fp,parton,spectator,emitter
    integer :: leg_map(pgl(igroup)%next)
    logical :: seen(n_integrated_histories)
    logical :: massive_history,massive_if_history

    pterm=0d0
    kterm=0d0
    if (present(pterm_copy)) then
       if (size(pterm_copy).ne.size(hard_copy)) then
          write(*,*) 'ERROR: integrated P copy array has incompatible size'
          stop 1
       endif
       pterm_copy=0d0
    endif
    if (present(kterm_copy)) then
       if (size(kterm_copy).ne.size(hard_copy)) then
          write(*,*) 'ERROR: integrated K copy array has incompatible size'
          stop 1
       endif
       kterm_copy=0d0
    endif
    if (z.le.0d0 .or. z.ge.1d0) return
    other=3-beam
    do icopy=1,size(hard_copy)
       if (hard_copy(icopy).eq.0d0) cycle
       b=pgl(igroup)%iden_processes(beam,icopy,iproc)
       call evaluate_pdf_flavour(pgl(igroup)%iden_processes(other,icopy,iproc),&
            xbj(other),mu_fac,other_pdf)
       if (other_pdf.eq.0d0) cycle
       seen=.false.
       do ih=1,n_integrated_histories
          if (seen(ih)) cycle
          if (integrated_history_list(ih)%incoming_leg.le.0) cycle
          if (.not.history_matches_copy(integrated_history_list(ih),igroup,iproc,icopy,leg_map)) cycle
          if (leg_map(integrated_history_list(ih)%incoming_leg).ne.beam) cycle
          a=integrated_history_list(ih)%real_incoming_flavour
          call parton_kind(a,pa)
          call parton_kind(b,pb)
          if (pa.eq.0 .or. pb.eq.0) cycle
          call evaluate_pdf_flavour(a,xbj(beam)/z,mu_fac,fa)
          call evaluate_pdf_flavour(b,xbj(beam),mu_fac,fb)
          if (z.ge.xbj(beam)) then
             gz=fa/z
          else
             gz=0d0
          endif
          g1=fb
          p_delta=0d0
          massive_history=history_has_mass(integrated_history_list(ih))
          massive_if_history=.false.
          if (massive_history) then
             fi=integrated_history_list(ih)%real_flavours(&
                  pgl(integrated_history_list(ih)%real_group)%dpl(&
                  integrated_history_list(ih)%real_process)%dl(&
                  integrated_history_list(ih)%local_dipole)%dip_ijk(1))
             fj=integrated_history_list(ih)%real_flavours(&
                  pgl(integrated_history_list(ih)%real_group)%dpl(&
                  integrated_history_list(ih)%real_process)%dl(&
                  integrated_history_list(ih)%local_dipole)%dip_ijk(2))
             fp=integrated_history_list(ih)%born_flavours(&
                  integrated_history_list(ih)%born_emitter)
             emitter=leg_map(integrated_history_list(ih)%born_emitter)
             spectator=leg_map(integrated_history_list(ih)%born_spectator)
             select case (integrated_history_list(ih)%topology)
             case (2)
                szone=abs(2d0*minkowski_dot(pgl(igroup)%ps(1)%p(:,emitter),&
                     pgl(igroup)%ps(1)%p(:,beam)))
                ! The real incoming momentum is p_a/z while the reduced
                ! Born momentum p_a is held fixed by the convolution.
                sx=szone/z
                call cs_massive_fi_convolution(integrated_history_list(ih)%parent_mass,&
                     sx,szone,z,alpha_dipole(2),massive_kernel,info)
             case (3)
                szone=abs(2d0*minkowski_dot(pgl(igroup)%ps(1)%p(:,beam),&
                     pgl(igroup)%ps(1)%p(:,spectator)))
                sx=szone/z
                call cs_massive_if_convolution(pa,pb,&
                     integrated_history_list(ih)%spectator_mass,sx,szone,z,&
                     mu_ren,mu_fac,alpha_dipole(3),&
                     integrated_n_active_flavours,massive_kernel,info)
                if (info.eq.0) then
                   call cs_ap_distribution(pa,pb,z,integrated_n_active_flavours,&
                        pk,info)
                   ! The auxiliary U(1) parent already carries one 1/Nc
                   ! through the local g->q qbar splitting trace.  Relative
                   ! to the matched physical-gluon Born contribution, the
                   ! initial-state convolution therefore needs only the
                   ! remaining 1/Nc factor.
                   if (info.eq.0 .and. fp.eq.99) &
                        call cs_scale_distribution(pk,cs_u1_initial_factor)
                endif
                if (info.eq.0) then
                   regularised_p=pk%regular*gz+pk%plus_one*(gz-g1)/(1d0-z)
                   p_delta=pk%delta*g1
                   massive_if_history=.true.
                endif
             case default
                info=-20
             end select
             if (info.ne.0) then
                write(*,*) 'ERROR: massive integrated convolution failed'
                write(*,*) ' history/topology/fi/fj/fp/z/status:',ih,&
                     integrated_history_list(ih)%topology,fi,fj,fp,z,info
                stop 1
             endif
             regularised_k=cs_apply_convolution(massive_kernel,gz,g1)
             if (fp.eq.99) regularised_k=regularised_k*cs_u1_initial_factor
             if (.not.massive_if_history) regularised_p=0d0
          elseif (integrated_history_list(ih)%topology.eq.2) then
             ! Final--initial histories have no universal P/Kbar term.  Their
             ! complete finite distribution, including the alpha=1 baseline,
             ! is a convolution on the initial spectator beam.
             fi=integrated_history_list(ih)%real_flavours(&
                  pgl(integrated_history_list(ih)%real_group)%dpl(&
                  integrated_history_list(ih)%real_process)%dl(&
                  integrated_history_list(ih)%local_dipole)%dip_ijk(1))
             fj=integrated_history_list(ih)%real_flavours(&
                  pgl(integrated_history_list(ih)%real_group)%dpl(&
                  integrated_history_list(ih)%real_process)%dl(&
                  integrated_history_list(ih)%local_dipole)%dip_ijk(2))
             fp=integrated_history_list(ih)%born_flavours(&
                  integrated_history_list(ih)%born_emitter)
             call history_i_primitive(fi,fj,fp,integrated_history_list(ih)%topology,&
                  primitive,parton)
             if (parton.eq.0) cycle
             call cs_fi_distribution(primitive,z,alpha_dipole(2),&
                  fi_regular,fi_subtracted,info)
             if (info.ne.0) cycle
             regularised_p=0d0
             regularised_k=fi_regular*gz+fi_subtracted*(gz-g1)
          else
             call cs_ap_distribution(pa,pb,z,integrated_n_active_flavours,pk,info)
             if (info.ne.0) cycle
             call cs_kbar_distribution(pa,pb,z,integrated_n_active_flavours,kk,info)
             if (info.ne.0) cycle
             select case (integrated_history_list(ih)%topology)
             case (3)
                call cs_if_tilde_distribution(pa,pb,z,integrated_n_active_flavours,&
                     tilde_kernel,info)
                if (info.ne.0) cycle
                call cs_if_alpha_distribution(pa,pb,z,integrated_n_active_flavours,alpha_dipole(3),&
                     alpha_kernel,info)
             case (4)
                call cs_ii_tilde_distribution(pa,pb,z,integrated_n_active_flavours,&
                     tilde_kernel,info)
                if (info.ne.0) cycle
                call cs_ii_alpha_distribution(pa,pb,z,integrated_n_active_flavours,alpha_dipole(4),&
                     alpha_kernel,info)
             case default
                tilde_kernel=cs_distribution()
                alpha_kernel=cs_distribution()
                info=0
             end select
             if (info.ne.0) cycle
             kk%regular=kk%regular+tilde_kernel%regular+alpha_kernel%regular
             kk%plus_one=kk%plus_one+tilde_kernel%plus_one+alpha_kernel%plus_one
             kk%plus_log=kk%plus_log+tilde_kernel%plus_log+alpha_kernel%plus_log
             kk%plus_log_one=kk%plus_log_one+tilde_kernel%plus_log_one+&
                  alpha_kernel%plus_log_one
             kk%delta=kk%delta+tilde_kernel%delta+alpha_kernel%delta
             fp=integrated_history_list(ih)%born_flavours(&
                  integrated_history_list(ih)%born_emitter)
             kernel_factor=1d0
             if (fp.eq.99) kernel_factor=cs_u1_initial_factor
             call cs_scale_distribution(pk,kernel_factor)
             call cs_scale_distribution(kk,kernel_factor)
             regularised_p=pk%regular*gz+pk%plus_one*(gz-g1)/(1d0-z)
             p_delta=pk%delta*g1
             regularised_k=kk%regular*gz+kk%plus_one*(gz-g1)/(1d0-z)+&
                  kk%plus_log*log((1d0-z)/z)*(gz-g1)/(1d0-z)+&
                  kk%plus_log_one*log(1d0-z)*(gz-g1)/(1d0-z)+kk%delta*g1
          endif
          spectator=leg_map(integrated_history_list(ih)%born_spectator)
          sij=abs(2d0*minkowski_dot(pgl(igroup)%ps(1)%p(:,beam),&
               pgl(igroup)%ps(1)%p(:,spectator)))
          colour_log=0d0
          colour_log_endpoint=0d0
          if (sij.gt.tiny(1d0)) then
             colour_log=log(mu_fac*mu_fac/(z*sij))
             colour_log_endpoint=log(mu_fac*mu_fac/sij)
          endif
          if (massive_if_history) then
             ! finiteif returns the complete mass-factorised convolution.
             ! Expose the same universal P contribution as the massless
             ! path and leave the compensating remainder in K, without
             ! changing their sum.
             regularised_k=regularised_k+colour_log*regularised_p+&
                  colour_log_endpoint*p_delta
          endif
          history_weight=pgl(integrated_history_list(ih)%real_group)%dpl(&
               integrated_history_list(ih)%real_process)%dl(&
               integrated_history_list(ih)%local_dipole)%lc_weight
          pterm=pterm-hard_copy(icopy)*other_pdf*history_weight*&
               (colour_log*regularised_p+colour_log_endpoint*p_delta)
          kterm=kterm+hard_copy(icopy)*other_pdf*history_weight*regularised_k
          if (present(pterm_copy)) pterm_copy(icopy)=pterm_copy(icopy)-&
               hard_copy(icopy)*other_pdf*history_weight*&
               (colour_log*regularised_p+colour_log_endpoint*p_delta)
          if (present(kterm_copy)) kterm_copy(icopy)=kterm_copy(icopy)+&
               hard_copy(icopy)*other_pdf*history_weight*regularised_k
          seen(ih)=.true.
       enddo
    enddo
    pterm=pterm*alpha_s/(2d0*cs_pi)
    kterm=kterm*alpha_s/(2d0*cs_pi)
    if (present(pterm_copy)) pterm_copy=pterm_copy*alpha_s/(2d0*cs_pi)
    if (present(kterm_copy)) kterm_copy=kterm_copy*alpha_s/(2d0*cs_pi)
  end subroutine integrated_beam

  subroutine history_i_primitive(fi,fj,fp,topology,coeff,parton)
    integer, intent(in) :: fi,fj,fp,topology
    real(kind=8), intent(out) :: coeff(-2:0)
    integer, intent(out) :: parton
    coeff=0d0
    parton=0
    if (abs(fp).le.6 .and. (fj.eq.21 .or. fj.eq.99)) then
       call cs_i_qg(coeff)
       parton=cs_parton_q
    elseif ((fp.eq.21 .or. fp.eq.99) .and. &
         (fi.eq.21 .or. fi.eq.99) .and. (fj.eq.21 .or. fj.eq.99)) then
       ! Every LC history is tied to one ordered colour neighbour.  The
       ! complementary soft pole belongs to the history on the other side of
       ! the gluon, independently of whether that neighbour is initial or
       ! final state.
       call cs_i_gg_ordered(coeff)
       parton=cs_parton_g
    elseif ((fp.eq.21 .or. fp.eq.99) .and. abs(fi).le.6 .and. abs(fi).eq.abs(fj) .and. &
         ((topology.le.2 .and. fi.eq.-fj) .or. topology.ge.3)) then
       call cs_i_qqbar(coeff)
       ! Relative to a matched physical-gluon Born contribution, an
       ! auxiliary U(1) history carries 1/Nc from its splitting trace and
       ! another 1/Nc from the reduced amplitude's LC colour factor.
       if (fp.eq.99) coeff=coeff/(cs_ca*cs_ca)
       parton=cs_parton_g
    endif
  end subroutine history_i_primitive

  subroutine parton_kind(flavour,kind)
    integer, intent(in) :: flavour
    integer, intent(out) :: kind
    if (abs(flavour).le.6) then
       kind=cs_parton_q
    elseif (flavour.eq.21 .or. flavour.eq.99) then
       kind=cs_parton_g
    else
       kind=0
    endif
  end subroutine parton_kind

  pure logical function history_has_mass(history)
    type(integrated_history), intent(in) :: history
    real(kind=8) :: tolerance
    tolerance=100d0*epsilon(1d0)*max(1d0,history%emitter_mass,&
         history%unresolved_mass,history%parent_mass,history%spectator_mass)
    history_has_mass=history%emitter_mass.gt.tolerance .or. &
         history%unresolved_mass.gt.tolerance .or. &
         history%parent_mass.gt.tolerance .or. &
         history%spectator_mass.gt.tolerance
  end function history_has_mass

  pure real(kind=8) function minkowski_dot(a,b)
    real(kind=8), intent(in) :: a(0:3),b(0:3)
    minkowski_dot=a(0)*b(0)-dot_product(a(1:3),b(1:3))
  end function minkowski_dot

  logical function is_duplicate_history(candidate,histories,nhistory)
    type(integrated_history), intent(in) :: candidate
    integer, intent(in) :: nhistory
    type(integrated_history), intent(in) :: histories(:)
    integer :: i

    is_duplicate_history=.false.
    do i=1,nhistory
       ! Distinct real colour flows can reduce to the same integrated
       ! insertion.  Its physical identity is the reduced colour string,
       ! emitter, spectator, topology, and splitting flavour—not the labels
       ! of the higher-multiplicity real process.
       if (.not.same_flavour_vector(candidate%born_flavours,histories(i)%born_flavours)) cycle
       if (.not.same_coloured_order(candidate%born_flavours,candidate%born_colour_order,&
            histories(i)%born_flavours,histories(i)%born_colour_order)) cycle
       if (candidate%born_emitter.ne.histories(i)%born_emitter) cycle
       if (candidate%born_spectator.ne.histories(i)%born_spectator) cycle
       ! Physical and auxiliary-U(1) parents have the same external quantum
       ! numbers for matching but different splitting normalisations.
       if (candidate%born_flavours(candidate%born_emitter).ne.&
            histories(i)%born_flavours(histories(i)%born_emitter)) cycle
       if (.not.same_integrated_splitting(candidate,histories(i))) cycle
       is_duplicate_history=.true.
       return
    enddo
  end function is_duplicate_history

  logical function same_integrated_splitting(first,second)
    type(integrated_history), intent(in) :: first,second
    integer :: first_ijk(3),second_ijk(3)
    integer :: fi1,fj1,fp1,fi2,fj2,fp2

    same_integrated_splitting=.false.
    if (first%topology.ne.second%topology) return
    first_ijk=pgl(first%real_group)%dpl(first%real_process)%dl(first%local_dipole)%dip_ijk
    second_ijk=pgl(second%real_group)%dpl(second%real_process)%dl(second%local_dipole)%dip_ijk
    fi1=first%real_flavours(first_ijk(1))
    fj1=first%real_flavours(first_ijk(2))
    fp1=first%born_flavours(first%born_emitter)
    fi2=second%real_flavours(second_ijk(1))
    fj2=second%real_flavours(second_ijk(2))
    fp2=second%born_flavours(second%born_emitter)
    if (fp1.ne.fp2) return

    if ((fp1.eq.21 .or. fp1.eq.99) .and. first%topology.le.2 .and. &
         abs(fi1).le.6 .and. fi1.eq.-fj1 .and. &
         abs(fi2).le.6 .and. fi2.eq.-fj2) then
       ! The two orientations of a final-state q-qbar pair have the same
       ! integrated kernel.
       same_integrated_splitting=abs(fi1).eq.abs(fi2)
    else
       same_integrated_splitting=fi1.eq.fi2 .and. fj1.eq.fj2
    endif
  end function same_integrated_splitting

  pure logical function same_flavour_vector(first,second)
    integer, intent(in) :: first(:),second(:)
    integer :: i,a,b
    same_flavour_vector=.false.
    if (size(first).ne.size(second)) return
    do i=1,size(first)
       a=first(i)
       b=second(i)
       if (a.eq.99) a=21
       if (b.eq.99) b=21
       if (a.ne.b) return
    enddo
    same_flavour_vector=.true.
  end function same_flavour_vector

  pure logical function same_flavour_pattern(template,actual)
    ! The process reader uses light-quark labels as placeholders and expands
    ! an amplitude into identical flavour copies.  Integrated histories are
    ! built from the placeholder process, so matching literal PDG values
    ! would retain only its first (u-like) copy.  Match the equality and sign
    ! pattern under a one-to-one relabelling instead.
    integer, intent(in) :: template(:),actual(:)
    integer :: mapped
    logical :: ok

    call map_history_flavour(template,actual,21,mapped,ok)
    same_flavour_pattern=ok
  end function same_flavour_pattern

  pure subroutine map_history_flavour(template,actual,template_flavour,mapped_flavour,ok)
    integer, intent(in) :: template(:),actual(:),template_flavour
    integer, intent(out) :: mapped_flavour
    logical, intent(out) :: ok
    integer :: flavour_map(6),reverse_map(6)
    integer :: i,t,a,at,aa

    ok=.false.
    mapped_flavour=0
    if (size(template).ne.size(actual)) return
    flavour_map=0
    reverse_map=0
    do i=1,size(template)
       t=template(i)
       a=actual(i)
       if (t.eq.99) t=21
       if (a.eq.99) a=21
       if (abs(t).ge.1 .and. abs(t).le.6) then
          if (.not.(abs(a).ge.1 .and. abs(a).le.6)) return
          if (sign(1,t).ne.sign(1,a)) return
          at=abs(t)
          aa=abs(a)
          if (flavour_map(at).eq.0) then
             if (reverse_map(aa).ne.0 .and. reverse_map(aa).ne.at) return
             flavour_map(at)=aa
             reverse_map(aa)=at
          elseif (flavour_map(at).ne.aa) then
             return
          endif
       elseif (t.ne.a) then
          return
       endif
    enddo

    if (template_flavour.eq.99 .or. template_flavour.eq.21) then
       mapped_flavour=template_flavour
    elseif (abs(template_flavour).ge.1 .and. abs(template_flavour).le.6) then
       at=abs(template_flavour)
       if (flavour_map(at).eq.0) return
       mapped_flavour=sign(flavour_map(at),template_flavour)
    else
       mapped_flavour=template_flavour
    endif
    ok=.true.
  end subroutine map_history_flavour

  logical function same_coloured_order(process1,order1,process2,order2)
    integer, intent(in) :: process1(:),order1(:),process2(:),order2(:)
    integer :: coloured1(size(order1)),coloured2(size(order2))
    integer :: n1,n2,i,j,shift

    same_coloured_order=.false.
    n1=0
    do i=1,size(order1)
       if (phys_model%is_singlet(process1(order1(i))) .and. process1(order1(i)).ne.99) cycle
       n1=n1+1
       coloured1(n1)=order1(i)
    enddo
    n2=0
    do i=1,size(order2)
       if (phys_model%is_singlet(process2(order2(i))) .and. process2(order2(i)).ne.99) cycle
       n2=n2+1
       coloured2(n2)=order2(i)
    enddo
    if (n1.ne.n2) return
    if (n1.eq.0) then
       same_coloured_order=.true.
       return
    endif
    do shift=0,n1-1
       do i=1,n1
          j=mod(i-1+shift,n1)+1
          if (coloured1(i).ne.coloured2(j)) exit
       enddo
       if (i.gt.n1) then
          same_coloured_order=.true.
          return
       endif
    enddo
  end function same_coloured_order

  subroutine match_coloured_legs(process1,order1,process2,order2,leg_map,matched)
    ! Match equivalent colour strings even when identical final-state gluons
    ! occupy different array slots.  Reduced U(1) parents are identified
    ! with the corresponding physical gluon, while incoming and final-state
    ! legs are never interchanged.  leg_map translates indices in process1
    ! to the matched indices in process2.
    integer, intent(in) :: process1(:),order1(:),process2(:),order2(:)
    integer, intent(out) :: leg_map(:)
    logical, intent(out) :: matched
    integer :: coloured1(size(order1)),coloured2(size(order2))
    integer :: trial_map(size(leg_map)),best_map(size(leg_map))
    integer :: n1,n2,i,j,shift,a,b,score,best_score
    logical :: valid

    matched=.false.
    leg_map=0
    if (size(process1).ne.size(process2)) return
    if (size(order1).ne.size(order2)) return
    if (size(leg_map).ne.size(process1)) return
    n1=0
    do i=1,size(order1)
       if (phys_model%is_singlet(process1(order1(i))) .and. process1(order1(i)).ne.99) cycle
       n1=n1+1
       coloured1(n1)=order1(i)
    enddo
    n2=0
    do i=1,size(order2)
       if (phys_model%is_singlet(process2(order2(i))) .and. process2(order2(i)).ne.99) cycle
       n2=n2+1
       coloured2(n2)=order2(i)
    enddo
    if (n1.ne.n2) return
    if (n1.eq.0) then
       matched=.true.
       return
    endif
    best_score=-1
    best_map=0
    do shift=0,n1-1
       trial_map=0
       score=0
       valid=.true.
       do i=1,n1
          j=mod(i-1+shift,n1)+1
          if ((coloured1(i).le.2) .neqv. (coloured2(j).le.2)) then
             valid=.false.
             exit
          endif
          a=process1(coloured1(i))
          b=process2(coloured2(j))
          if (a.eq.99) a=21
          if (b.eq.99) b=21
          if (a.ne.b) then
             valid=.false.
             exit
          endif
          trial_map(coloured1(i))=coloured2(j)
          if (coloured1(i).eq.coloured2(j)) score=score+1
       enddo
       if (valid .and. score.gt.best_score) then
          best_score=score
          best_map=trial_map
       endif
    enddo
    if (best_score.ge.0) then
       leg_map=best_map
       matched=.true.
    endif
  end subroutine match_coloured_legs

end module integrated_dipoles
