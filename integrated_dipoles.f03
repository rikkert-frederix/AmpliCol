module integrated_dipoles
  ! Inventory of integrated subtraction histories.  The inventory is built
  ! from the local dipoles themselves; it never independently enumerates
  ! splittings.  This is the central consistency guarantee between D and
  ! I/P/K at leading colour.
  use handling_processes
  use cs_dipole_mappings, only: cs_dipole_topology
  use cs_integrated_kernels
  use cs_massive_integrated_kernels
  use pdf_wrap, only: evaluate_pdf_table,pdf_table_flavour
  use run_parameters, only: alpha_dipole
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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

  type(integrated_history), allocatable, save :: integrated_history_list(:)
  integer, save :: n_integrated_histories=0
  integer, save :: integrated_dimensional_scheme=cs_scheme_hv
  integer, save :: integrated_n_active_flavours=5
  real(kind=8), parameter :: integrated_value_limit=0.25d0*huge(1d0)**0.25d0
  integer(kind=8), parameter :: max_integrated_histories=20000000_8

  public :: initialise_integrated_dipoles, history_matches_born
  public :: integrated_endpoint, integrated_beam

contains

  pure logical function light_quark_code(flavour)
    integer, intent(in) :: flavour
    light_quark_code=(flavour.ge.1 .and. flavour.le.6) .or. &
         (flavour.le.-1 .and. flavour.ge.-6)
  end function light_quark_code

  pure integer function light_quark_magnitude(flavour)
    integer, intent(in) :: flavour
    light_quark_magnitude=0
    if (flavour.ge.1 .and. flavour.le.6) then
       light_quark_magnitude=flavour
    elseif (flavour.le.-1 .and. flavour.ge.-6) then
       light_quark_magnitude=-flavour
    endif
  end function light_quark_magnitude

  pure logical function same_light_quark_flavour(first,second)
    integer, intent(in) :: first,second
    same_light_quark_flavour=.false.
    if (.not.light_quark_code(first) .or. .not.light_quark_code(second)) return
    same_light_quark_flavour=first.eq.second .or. first.eq.-second
  end function same_light_quark_flavour

  pure logical function light_quark_antiquark_pair(first,second)
    integer, intent(in) :: first,second
    light_quark_antiquark_pair=.false.
    if (.not.light_quark_code(first) .or. .not.light_quark_code(second)) return
    light_quark_antiquark_pair=first.eq.-second
  end function light_quark_antiquark_pair

  logical function valid_born_process_state(igroup,iproc,ncopy,require_momenta)
    integer, intent(in) :: igroup,iproc
    integer, intent(out) :: ncopy
    logical, intent(in), optional :: require_momenta
    logical :: need_momenta

    valid_born_process_state=.false.
    ncopy=0
    need_momenta=.false.
    if (present(require_momenta)) need_momenta=require_momenta
    if (.not.allocated(pgl)) return
    if (igroup.lt.1) return
    if (igroup.gt.size(pgl)) return
    if (pgl(igroup)%next.lt.1) return
    if (pgl(igroup)%nproc.lt.1) return
    if (iproc.lt.1) return
    if (iproc.gt.pgl(igroup)%nproc) return
    if (.not.allocated(pgl(igroup)%iden_iproc)) return
    if (size(pgl(igroup)%iden_iproc).lt.pgl(igroup)%nproc) return
    if (.not.allocated(pgl(igroup)%iden_processes)) return
    if (size(pgl(igroup)%iden_processes,1).ne.pgl(igroup)%next) return
    if (size(pgl(igroup)%iden_processes,3).lt.pgl(igroup)%nproc) return
    if (.not.allocated(pgl(igroup)%color_orders)) return
    if (size(pgl(igroup)%color_orders,1).ne.pgl(igroup)%next) return
    if (size(pgl(igroup)%color_orders,2).lt.pgl(igroup)%nproc) return
    ncopy=pgl(igroup)%iden_iproc(iproc)
    if (ncopy.lt.1) return
    if (ncopy.gt.size(pgl(igroup)%iden_processes,2)) return
    if (need_momenta) then
       if (.not.allocated(pgl(igroup)%ps)) return
       if (size(pgl(igroup)%ps).lt.1) return
       if (.not.allocated(pgl(igroup)%ps(1)%p)) return
       if (size(pgl(igroup)%ps(1)%p,1).ne.4) return
       if (size(pgl(igroup)%ps(1)%p,2).lt.pgl(igroup)%next) return
    endif
    valid_born_process_state=.true.
  end function valid_born_process_state

  logical function integrated_registry_is_ready()
    integrated_registry_is_ready=.false.
    if (n_integrated_histories.lt.1) return
    if (.not.allocated(integrated_history_list)) return
    if (size(integrated_history_list).ne.n_integrated_histories) return
    integrated_registry_is_ready=.true.
  end function integrated_registry_is_ready

  logical function valid_local_dipole_history(igroup,iproc,idip)
    integer, intent(in) :: igroup,iproc,idip
    integer :: ncopy,reduced_size

    valid_local_dipole_history=.false.
    if (.not.valid_born_process_state(igroup,iproc,ncopy)) return
    if (.not.allocated(pgl(igroup)%processes)) return
    if (size(pgl(igroup)%processes,1).ne.pgl(igroup)%next) return
    if (size(pgl(igroup)%processes,2).lt.pgl(igroup)%nproc) return
    if (.not.allocated(pgl(igroup)%dpl)) return
    if (iproc.gt.size(pgl(igroup)%dpl)) return
    if (.not.allocated(pgl(igroup)%dpl(iproc)%dl)) return
    if (idip.lt.1 .or. idip.gt.pgl(igroup)%dpl(iproc)%ndip) return
    if (idip.gt.size(pgl(igroup)%dpl(iproc)%dl)) return
    associate(local_dipole=>pgl(igroup)%dpl(iproc)%dl(idip))
      if (pgl(igroup)%next.lt.3) return
      if (any(local_dipole%dip_ijk.lt.1) .or. &
           any(local_dipole%dip_ijk.gt.pgl(igroup)%next)) return
      if (count(local_dipole%dip_ijk.eq.local_dipole%dip_ijk(1)).ne.1 .or. &
           count(local_dipole%dip_ijk.eq.local_dipole%dip_ijk(2)).ne.1 .or. &
           count(local_dipole%dip_ijk.eq.local_dipole%dip_ijk(3)).ne.1) return
      reduced_size=pgl(igroup)%next-1
      if (any(local_dipole%dip_r_ijk.lt.1) .or. &
           any(local_dipole%dip_r_ijk.gt.reduced_size)) return
      if (local_dipole%dip_r_ijk(1).eq.local_dipole%dip_r_ijk(2)) return
      if (.not.allocated(local_dipole%process_r)) return
      if (.not.allocated(local_dipole%reduced_color_order)) return
      if (size(local_dipole%process_r).ne.reduced_size) return
      if (.not.valid_colour_order(local_dipole%process_r,&
           local_dipole%reduced_color_order)) return
      if (.not.integrated_value_is_safe(local_dipole%lc_weight)) return
      if (cs_dipole_topology(local_dipole%dip_ijk).lt.1 .or. &
           cs_dipole_topology(local_dipole%dip_ijk).gt.4) return
    end associate
    valid_local_dipole_history=.true.
  end function valid_local_dipole_history

  subroutine initialise_integrated_dipoles(nborn_groups,scheme_name,n_active_flavours)
    integer, intent(in) :: nborn_groups
    integer, intent(in) :: n_active_flavours
    character(len=*), intent(in) :: scheme_name
    type(integrated_history), allocatable :: candidates(:),unique_histories(:)
    integer :: igroup,iproc,icopy,idip,max_histories,ncandidate,i,ncopy
    integer :: allocation_status
    integer(kind=8) :: max_histories64,additional_histories
    character(len=256) :: allocation_message
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

    if (.not.allocated(pgl)) then
       write(*,*) 'ERROR: integrated subtraction requires initialized process groups'
       stop 1
    endif
    if (ngroups.ne.size(pgl) .or. nborn_groups.lt.1 .or. &
         nborn_groups.ge.ngroups) then
       write(*,*) 'ERROR: invalid Born/real group partition for integrated subtraction:',&
            nborn_groups,ngroups,size(pgl)
       stop 1
    endif
    do igroup=1,nborn_groups
       if (pgl(igroup)%nproc.lt.1) then
          write(*,*) 'ERROR: Born group has no processes for integrated subtraction:',igroup
          stop 1
       endif
       do iproc=1,pgl(igroup)%nproc
          if (.not.valid_born_process_state(igroup,iproc,ncopy)) then
             write(*,*) 'ERROR: malformed Born process metadata for integrated subtraction:',&
                  igroup,iproc
             stop 1
          endif
       enddo
    enddo

    max_histories64=0_8
    do igroup=nborn_groups+1,ngroups
       if (.not.pgl(igroup)%is_subtracted_real) cycle
       if (pgl(igroup)%nproc.lt.1) then
          write(*,*) 'ERROR: subtracted-real group has no processes:',igroup
          stop 1
       endif
       if (.not.allocated(pgl(igroup)%processes)) then
          write(*,*) 'ERROR: subtracted-real group has no process table:',igroup
          stop 1
       endif
       if (size(pgl(igroup)%processes,1).ne.pgl(igroup)%next .or. &
            size(pgl(igroup)%processes,2).lt.pgl(igroup)%nproc) then
          write(*,*) 'ERROR: malformed subtracted-real process table:',igroup
          stop 1
       endif
       if (.not.allocated(pgl(igroup)%dpl)) then
          write(*,*) 'ERROR: subtracted-real group has no local dipoles:',igroup
          stop 1
       endif
       if (size(pgl(igroup)%dpl).lt.pgl(igroup)%nproc) then
          write(*,*) 'ERROR: malformed local-dipole table:',igroup
          stop 1
       endif
       do iproc=1,pgl(igroup)%nproc
          if (.not.valid_born_process_state(igroup,iproc,ncopy)) then
             write(*,*) 'ERROR: malformed real-process metadata for integrated subtraction:',&
                  igroup,iproc
             stop 1
          endif
          if (pgl(igroup)%dpl(iproc)%ndip.lt.0) then
             write(*,*) 'ERROR: negative local-dipole count:',igroup,iproc
             stop 1
          endif
          if (pgl(igroup)%dpl(iproc)%ndip.gt.0) then
             if (.not.allocated(pgl(igroup)%dpl(iproc)%dl)) then
                write(*,*) 'ERROR: missing local-dipole list:',igroup,iproc
                stop 1
             endif
             if (size(pgl(igroup)%dpl(iproc)%dl).lt.pgl(igroup)%dpl(iproc)%ndip) then
                write(*,*) 'ERROR: truncated local-dipole list:',igroup,iproc
                stop 1
             endif
             do idip=1,pgl(igroup)%dpl(iproc)%ndip
                if (.not.valid_local_dipole_history(igroup,iproc,idip)) then
                   write(*,*) 'ERROR: malformed local dipole for integrated subtraction:',&
                        igroup,iproc,idip
                   stop 1
                endif
             enddo
          endif
          additional_histories=int(ncopy,kind=8)*&
               int(pgl(igroup)%dpl(iproc)%ndip,kind=8)
          if (additional_histories.gt.max_integrated_histories-max_histories64) then
             write(*,*) 'ERROR: integrated-history registry exceeds supported size'
             stop 1
          endif
          max_histories64=max_histories64+additional_histories
       enddo
    enddo
    if (max_histories64.eq.0_8) then
       write(*,*) 'ERROR: no local dipoles are available for integrated subtraction'
       stop 1
    endif
    max_histories=int(max_histories64)

    allocation_message=''
    allocate(candidates(max_histories),stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write(*,*) 'ERROR: could not allocate integrated-history candidates: ',&
            trim(allocation_message)
       stop 1
    endif
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

    allocation_message=''
    allocate(unique_histories(ncandidate),stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write(*,*) 'ERROR: could not allocate integrated-history registry: ',&
            trim(allocation_message)
       stop 1
    endif
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
    integer :: incomplete_flavour_emitters,unsupported_four_line_emitters
    integer :: first_group,first_process,first_copy,first_emitter,first_nf,quark_lines
    integer :: mapped_emitter,mapped_spectator
    integer :: leg_map(maxval(pgl(1:nborn_groups)%next))
    logical :: invalid,massive_emitter_histories

    incomplete_flavour_emitters=0
    unsupported_four_line_emitters=0
    first_group=0
    first_process=0
    first_copy=0
    first_emitter=0
    first_nf=0
    do igroup=1,nborn_groups
       do iproc=1,pgl(igroup)%nproc
          do icopy=1,pgl(igroup)%iden_iproc(iproc)
             do emitter=1,pgl(igroup)%next
                nf_available=integrated_n_active_flavours
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
                if (parton.eq.cs_parton_g .and. nf_available.lt.integrated_n_active_flavours) then
                   incomplete_flavour_emitters=incomplete_flavour_emitters+1
                   quark_lines=count_quark_lines(pgl(igroup)%iden_processes(:,icopy,iproc))
                   if (quark_lines.ge.3) unsupported_four_line_emitters=unsupported_four_line_emitters+1
                   if (incomplete_flavour_emitters.eq.1) then
                      first_group=igroup
                      first_process=iproc
                      first_copy=icopy
                      first_emitter=emitter
                      first_nf=nf_available
                   endif
                endif
             enddo
          enddo
       enddo
    enddo
    if (incomplete_flavour_emitters.gt.0) then
       write(*,'(a,i0,a,i0,a)') 'WARNING: ',incomplete_flavour_emitters,&
            ' Born gluon emitter(s) have incomplete g -> q qbar flavour coverage; the first has ',first_nf,' sectors.'
       write(*,'(a,i0,a)') ' The configured integrated kernels use nf=',integrated_n_active_flavours,&
            '; this process list is therefore not a complete full-nf NLO real contribution.'
       write(*,'(a,4(1x,i0))') ' First affected Born group/process/copy/emitter:',&
            first_group,first_process,first_copy,first_emitter
       if (unsupported_four_line_emitters.gt.0) then
          write(*,'(a,i0,a)') ' Of these, ',unsupported_four_line_emitters,&
               ' emitter(s) belong to a three-quark-line Born state: their real q-qbar sectors require four lines.'
          write(*,'(a)') ' Four-quark-line amplitudes are not implemented, so --include_3qqbar cannot complete those NLO processes.'
       endif
       if (unsupported_four_line_emitters.lt.incomplete_flavour_emitters) &
            write(*,'(a)') ' For lower-line Born states, regenerate with --include_3qqbar when those sectors are physical.'
       write(99,'(a,i0,a,i0,a)') 'WARNING: ',incomplete_flavour_emitters,&
            ' Born gluon emitter(s) have incomplete g -> q qbar flavour coverage; the first has ',first_nf,' sectors.'
       write(99,'(a,i0,a)') ' The configured integrated kernels use nf=',integrated_n_active_flavours,&
            '; this process list is therefore not a complete full-nf NLO real contribution.'
       write(99,'(a,4(1x,i0))') ' First affected Born group/process/copy/emitter:',&
            first_group,first_process,first_copy,first_emitter
       if (unsupported_four_line_emitters.gt.0) then
          write(99,'(a,i0,a)') ' Of these, ',unsupported_four_line_emitters,&
               ' emitter(s) belong to a three-quark-line Born state: their real q-qbar sectors require four lines.'
          write(99,'(a)') ' Four-quark-line amplitudes are not implemented, so --include_3qqbar cannot complete'//&
               ' those NLO processes.'
       endif
       if (unsupported_four_line_emitters.lt.incomplete_flavour_emitters) &
            write(99,'(a)') ' For lower-line Born states, regenerate with --include_3qqbar when those sectors are physical.'
    endif
  end subroutine validate_integrated_history_poles

  integer function count_quark_lines(flavours)
    integer, intent(in) :: flavours(:)
    integer :: i,nfermions
    nfermions=0
    do i=1,size(flavours)
       if (phys_model%is_quark(flavours(i)) .or. phys_model%is_antiquark(flavours(i))) &
            nfermions=nfermions+1
    enddo
    count_quark_lines=nfermions/2
  end function count_quark_lines

  subroutine fill_history(history,igroup,iproc,icopy,idip)
    type(integrated_history), intent(out) :: history
    integer, intent(in) :: igroup,iproc,icopy,idip
    integer :: emitter,spectator

    associate(dip=>pgl(igroup)%dpl(iproc)%dl(idip))
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
    end associate
  end subroutine fill_history

  pure subroutine make_history_reduced_process(base_real,copy_real,base_reduced,copy_reduced)
    integer, intent(in) :: base_real(:),copy_real(:),base_reduced(:)
    integer, intent(out) :: copy_reduced(:)
    integer :: flavour_map(6)
    integer :: i,base_flavour,copy_flavour

    copy_reduced=0
    if (size(copy_real).ne.size(base_real) .or. &
         size(copy_reduced).ne.size(base_reduced)) return
    flavour_map=0
    do i=1,size(base_real)
       if (.not.light_quark_code(base_real(i))) cycle
       if (.not.light_quark_code(copy_real(i))) cycle
       base_flavour=light_quark_magnitude(base_real(i))
       copy_flavour=light_quark_magnitude(copy_real(i))
       if (flavour_map(base_flavour).eq.0) flavour_map(base_flavour)=copy_flavour
    enddo

    copy_reduced=base_reduced
    do i=1,size(copy_reduced)
       if (.not.light_quark_code(copy_reduced(i))) cycle
       base_flavour=light_quark_magnitude(copy_reduced(i))
       if (flavour_map(base_flavour).eq.0) cycle
       copy_reduced(i)=sign(flavour_map(base_flavour),copy_reduced(i))
    enddo
  end subroutine make_history_reduced_process

  subroutine ensure_supported_integrated_history(igroup,iproc,idip)
    ! The unresolved leg and both incoming legs remain massless.  A massive
    ! emitter is supported only for final-state Q->Qg, while a massive
    ! spectator is supported in FF and IF histories.
    integer, intent(in) :: igroup,iproc,idip
    real(kind=8) :: emitter_mass,unresolved_mass,parent_mass,spectator_mass,mass_tolerance
    integer :: topology,fi,fj,fp
    logical :: supported

    associate(dip=>pgl(igroup)%dpl(iproc)%dl(idip))
      fi=dip%dip_ijk_f(1)
      fj=dip%dip_ijk_f(2)
      fp=dip%dip_r_ijk_f(1)
      emitter_mass=abs(phys_model%get_mass(dip%dip_ijk_f(1)))
      unresolved_mass=abs(phys_model%get_mass(dip%dip_ijk_f(2)))
      parent_mass=abs(phys_model%get_mass(dip%dip_r_ijk_f(1)))
      spectator_mass=abs(phys_model%get_mass(dip%dip_r_ijk_f(2)))
      mass_tolerance=100d0*epsilon(1d0)*max(tiny(1d0),emitter_mass,&
           unresolved_mass,parent_mass,spectator_mass)
      if (emitter_mass.eq.0d0 .and. unresolved_mass.eq.0d0 .and. &
           parent_mass.eq.0d0 .and. spectator_mass.eq.0d0) return

      topology=cs_dipole_topology(dip%dip_ijk)
      supported=.true.
      if (unresolved_mass.gt.0d0) supported=.false.
      if (emitter_mass.gt.0d0 .or. parent_mass.gt.0d0) then
         if (topology.ne.1 .and. topology.ne.2) supported=.false.
         if (.not.same_light_quark_flavour(fi,fp)) supported=.false.
         if (fj.ne.21 .and. fj.ne.99) supported=.false.
         if (abs(emitter_mass-parent_mass).gt.mass_tolerance) supported=.false.
      endif
      if (spectator_mass.gt.0d0 .and. topology.ne.1 .and. topology.ne.3) supported=.false.
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
    end associate
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
    if (.not.allocated(history%born_flavours) .or. &
         .not.allocated(history%born_colour_order)) return
    if (.not.allocated(pgl)) return
    if (igroup.lt.1 .or. igroup.gt.size(pgl)) return
    if (pgl(igroup)%nproc.lt.1 .or. iproc.lt.1 .or. &
         iproc.gt.pgl(igroup)%nproc) return
    if (.not.allocated(pgl(igroup)%iden_iproc) .or. &
         .not.allocated(pgl(igroup)%iden_processes) .or. &
         .not.allocated(pgl(igroup)%color_orders)) return
    if (size(pgl(igroup)%iden_iproc).lt.pgl(igroup)%nproc .or. &
         size(pgl(igroup)%iden_processes,1).ne.pgl(igroup)%next .or. &
         size(pgl(igroup)%iden_processes,3).lt.pgl(igroup)%nproc .or. &
         size(pgl(igroup)%color_orders,1).ne.pgl(igroup)%next .or. &
         size(pgl(igroup)%color_orders,2).lt.pgl(igroup)%nproc) return
    if (pgl(igroup)%next.ne.size(history%born_flavours)) return
    if (pgl(igroup)%iden_iproc(iproc).lt.1 .or. &
         pgl(igroup)%iden_iproc(iproc).gt.&
         size(pgl(igroup)%iden_processes,2)) return
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
    integer :: i,ncopy

    history_matches_copy=.false.
    if (present(leg_map)) leg_map=0
    if (.not.allocated(history%born_flavours)) return
    if (.not.allocated(history%born_colour_order)) return
    if (.not.valid_born_process_state(igroup,iproc,ncopy)) return
    if (icopy.lt.1 .or. icopy.gt.ncopy) return
    if (pgl(igroup)%next.ne.size(history%born_flavours)) return
    if (size(history%born_colour_order).ne.size(history%born_flavours)) return
    if (history%born_emitter.lt.1 .or. &
         history%born_emitter.gt.size(history%born_flavours)) return
    if (history%born_spectator.lt.1 .or. &
         history%born_spectator.gt.size(history%born_flavours)) return
    if (history%incoming_leg.lt.0 .or. &
         history%incoming_leg.gt.size(history%born_flavours)) return
    if (present(leg_map)) then
       if (size(leg_map).ne.size(history%born_flavours)) return
    endif
    if (.not.same_flavour_vector(history%born_flavours,&
         pgl(igroup)%iden_processes(:,icopy,iproc))) return
    if (.not.same_coloured_order(history%born_flavours,history%born_colour_order,&
         pgl(igroup)%iden_processes(:,icopy,iproc),pgl(igroup)%color_orders(:,iproc))) return
    history_matches_copy=.true.
    if (present(leg_map)) leg_map=[(i,i=1,size(leg_map))]
  end function history_matches_copy

  subroutine history_splitting_data(history,fi,fj,fp,weight,valid)
    type(integrated_history), intent(in) :: history
    integer, intent(out) :: fi,fj,fp
    real(kind=8), intent(out) :: weight
    logical, intent(out) :: valid
    integer :: ijk(3),real_group,real_process,local_dipole

    fi=0
    fj=0
    fp=0
    weight=0d0
    valid=.false.
    if (.not.allocated(history%real_flavours)) return
    if (.not.allocated(history%real_colour_order)) return
    if (.not.allocated(history%born_flavours)) return
    if (.not.allocated(history%born_colour_order)) return
    if (size(history%real_flavours).lt.1) return
    if (.not.valid_colour_order(history%real_flavours,history%real_colour_order)) return
    if (.not.valid_colour_order(history%born_flavours,history%born_colour_order)) return
    if (history%born_emitter.lt.1 .or. &
         history%born_emitter.gt.size(history%born_flavours)) return
    if (history%born_spectator.lt.1 .or. &
         history%born_spectator.gt.size(history%born_flavours)) return
    if (history%incoming_leg.lt.0 .or. &
         history%incoming_leg.gt.size(history%born_flavours)) return
    if (history%topology.lt.1 .or. history%topology.gt.4) return
    if (.not.all(ieee_is_finite([history%emitter_mass,history%unresolved_mass,&
         history%parent_mass,history%spectator_mass]))) return
    if (history%emitter_mass.lt.0d0 .or. history%unresolved_mass.lt.0d0 .or. &
         history%parent_mass.lt.0d0 .or. history%spectator_mass.lt.0d0) return
    if (.not.allocated(pgl)) return
    real_group=history%real_group
    if (real_group.lt.1 .or. real_group.gt.size(pgl)) return
    real_process=history%real_process
    if (real_process.lt.1 .or. real_process.gt.pgl(real_group)%nproc) return
    if (.not.allocated(pgl(real_group)%dpl)) return
    if (real_process.gt.size(pgl(real_group)%dpl)) return
    if (.not.allocated(pgl(real_group)%dpl(real_process)%dl)) return
    local_dipole=history%local_dipole
    if (local_dipole.lt.1) return
    if (local_dipole.gt.pgl(real_group)%dpl(real_process)%ndip) return
    if (local_dipole.gt.size(pgl(real_group)%dpl(real_process)%dl)) return
    ijk=pgl(real_group)%dpl(real_process)%dl(local_dipole)%dip_ijk
    if (any(ijk.lt.1) .or. any(ijk.gt.size(history%real_flavours))) return
    if (cs_dipole_topology(ijk).ne.history%topology) return
    fi=history%real_flavours(ijk(1))
    fj=history%real_flavours(ijk(2))
    fp=history%born_flavours(history%born_emitter)
    weight=pgl(real_group)%dpl(real_process)%dl(local_dipole)%lc_weight
    if (.not.integrated_value_is_safe(weight)) return
    valid=.true.
  end subroutine history_splitting_data

  subroutine integrated_endpoint(igroup,iproc,born_copy,p,mu_ren,alpha_s,coeff,coeff_copy,status)
    ! Laurent coefficients of I(eps), multiplied by the corresponding
    ! PDF-weighted Born contributions.  Every entry is inherited from an
    ! actual local-dipole history.
    integer, intent(in) :: igroup,iproc
    real(kind=8), intent(in) :: born_copy(:),p(0:,:)
    real(kind=8), intent(in) :: mu_ren,alpha_s
    real(kind=8), intent(out) :: coeff(-2:0)
    real(kind=8), intent(out),optional :: coeff_copy(-2:,:)
    integer, intent(out),optional :: status
    real(kind=8) :: primitive(-2:0),expanded(-2:0),sij,ell,weight
    real(kind=8) :: endpoint_alpha
    real(kind=8) :: copy_contribution(-2:0),alpha_factor,updated_value
    integer :: ih,icopy,fi,fj,fp,parton,info,ipole,ncopy,allocation_status
    integer :: emitter,spectator
    integer, allocatable :: leg_map(:)
    logical, allocatable :: shifted(:,:)
    logical :: massive_history,arithmetic_ok,history_valid
    character(len=256) :: allocation_message

    coeff=0d0
    if (present(status)) status=0
    if (present(coeff_copy)) coeff_copy=0d0
    if (.not.all(integrated_value_is_safe(born_copy)) .or. &
         .not.all(integrated_value_is_safe(p)) .or. &
         .not.integrated_value_is_safe(mu_ren) .or. &
         .not.integrated_value_is_safe(alpha_s)) then
       if (present(status)) then
          status=-20
          return
       endif
       write(*,*) 'ERROR: invalid input to integrated endpoint'
       stop 1
    endif
    if (mu_ren.le.0d0 .or. alpha_s.lt.0d0 .or. size(p,1).ne.4) then
       if (present(status)) then
          status=-20
          return
       endif
       write(*,*) 'ERROR: invalid input to integrated endpoint'
       stop 1
    endif
    if (.not.valid_born_process_state(igroup,iproc,ncopy)) then
       call endpoint_failure(-4,'invalid Born group/process state')
       return
    endif
    if (size(born_copy).ne.ncopy .or. size(p,2).lt.pgl(igroup)%next) then
       call endpoint_failure(-4,'incompatible integrated endpoint input shape')
       return
    endif
    if (.not.integrated_registry_is_ready()) then
       call endpoint_failure(-4,'integrated-history registry is not initialized')
       return
    endif
    if (present(coeff_copy)) then
       if (size(coeff_copy,1).ne.3 .or. size(coeff_copy,2).ne.size(born_copy)) then
          write(*,*) 'ERROR: integrated endpoint copy array has incompatible size'
          stop 1
       endif
    endif
    allocation_message=''
    allocate(leg_map(pgl(igroup)%next),shifted(size(born_copy),pgl(igroup)%next),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       call endpoint_failure(-4,'could not allocate integrated endpoint workspace: '//&
            trim(allocation_message))
       return
    endif
    shifted=.false.
    do icopy=1,size(born_copy)
       if (born_copy(icopy).eq.0d0) cycle
       do ih=1,n_integrated_histories
          if (.not.history_matches_copy(integrated_history_list(ih),igroup,iproc,icopy,leg_map)) cycle
          call history_splitting_data(integrated_history_list(ih),fi,fj,fp,weight,history_valid)
          if (.not.history_valid) then
             call endpoint_failure(-4,'invalid integrated endpoint history')
             return
          endif
          emitter=leg_map(integrated_history_list(ih)%born_emitter)
          spectator=leg_map(integrated_history_list(ih)%born_spectator)
          massive_history=history_has_mass(integrated_history_list(ih))
          if (massive_history) then
             call massive_history_endpoint(integrated_history_list(ih),fi,fj,fp,p,&
                  emitter,spectator,mu_ren,expanded,parton,info)
             if (info.ne.0) then
                if (present(status)) then
                   status=info
                   coeff=0d0
                   if (present(coeff_copy)) coeff_copy=0d0
                   return
                endif
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
          if (.not.ieee_is_finite(sij)) then
             if (present(status)) then
                status=-20
                coeff=0d0
                if (present(coeff_copy)) coeff_copy=0d0
                return
             endif
             write(*,*) 'ERROR: non-finite invariant in integrated endpoint:',ih,sij
             stop 1
          endif
          if (sij.le.0d0) then
             if (present(status)) then
                status=-20
                coeff=0d0
                if (present(coeff_copy)) coeff_copy=0d0
                return
             endif
             write(*,*) 'ERROR: invalid invariant in integrated endpoint:',ih,sij
             stop 1
          endif
          if (.not.massive_history) then
             ell=log(4d0*cs_pi)+2d0*log(mu_ren)-log(sij)-0.5772156649015328606d0
             expanded(-2)=primitive(-2)
             expanded(-1)=primitive(-1)+ell*primitive(-2)
             expanded(0)=primitive(0)+ell*primitive(-1)+&
                  (0.5d0*ell*ell-cs_pi**2/12d0)*primitive(-2)
          endif
          if (.not.all(integrated_value_is_safe(expanded))) then
             if (present(status)) then
                status=-20
                coeff=0d0
                if (present(coeff_copy)) coeff_copy=0d0
                return
             endif
             write(*,*) 'ERROR: non-finite integrated endpoint coefficients:',ih
             stop 1
          endif
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
                if (present(status)) then
                   status=info
                   coeff=0d0
                   if (present(coeff_copy)) coeff_copy=0d0
                   return
                endif
                write(*,*) 'ERROR: invalid integrated endpoint topology/alpha:',&
                     integrated_history_list(ih)%topology,info
                stop 1
             endif
             expanded(0)=expanded(0)+endpoint_alpha
          endif
          if (.not.all(integrated_value_is_safe(expanded))) then
             call endpoint_failure(-20,'unsafe alpha-restricted endpoint coefficients')
             return
          endif
          call safe_quotient(alpha_s,2d0*cs_pi,alpha_factor,arithmetic_ok)
          if (.not.arithmetic_ok) then
             call endpoint_failure(-20,'unsafe integrated endpoint coupling factor')
             return
          endif
          do ipole=-2,0
             call safe_product_four(born_copy(icopy),weight,expanded(ipole),alpha_factor,&
                  copy_contribution(ipole),arithmetic_ok)
             if (.not.arithmetic_ok) then
                call endpoint_failure(-20,'unsafe integrated endpoint contribution')
                return
             endif
          enddo
          do ipole=-2,0
             call safe_sum(coeff(ipole),copy_contribution(ipole),updated_value,arithmetic_ok)
             if (.not.arithmetic_ok) then
                call endpoint_failure(-20,'unsafe integrated endpoint accumulation')
                return
             endif
             coeff(ipole)=updated_value
             if (present(coeff_copy)) then
                call safe_sum(coeff_copy(ipole,icopy),copy_contribution(ipole),&
                     updated_value,arithmetic_ok)
                if (.not.arithmetic_ok) then
                   call endpoint_failure(-20,'unsafe integrated endpoint copy accumulation')
                   return
                endif
                coeff_copy(ipole,icopy)=updated_value
             endif
          enddo
          if (.not.massive_history .and. integrated_dimensional_scheme.eq.cs_scheme_fdh .and. &
               phys_model%get_mass(pgl(igroup)%iden_processes(emitter,icopy,iproc)).eq.0d0 .and. &
               .not.shifted(icopy,emitter)) then
             copy_contribution=0d0
             call safe_product_three(born_copy(icopy),cs_fdh_endpoint_shift(parton),alpha_factor,&
                  copy_contribution(0),arithmetic_ok)
             if (.not.arithmetic_ok) then
                call endpoint_failure(-20,'unsafe FDH endpoint shift')
                return
             endif
             call safe_sum(coeff(0),copy_contribution(0),updated_value,arithmetic_ok)
             if (.not.arithmetic_ok) then
                call endpoint_failure(-20,'unsafe FDH endpoint accumulation')
                return
             endif
             coeff(0)=updated_value
             if (present(coeff_copy)) then
                call safe_sum(coeff_copy(0,icopy),copy_contribution(0),updated_value,arithmetic_ok)
                if (.not.arithmetic_ok) then
                   call endpoint_failure(-20,'unsafe FDH endpoint copy accumulation')
                   return
                endif
                coeff_copy(0,icopy)=updated_value
             endif
             shifted(icopy,emitter)=.true.
          endif
       enddo
    enddo
    if (.not.all(integrated_value_is_safe(coeff))) then
       coeff=0d0
       if (present(coeff_copy)) coeff_copy=0d0
       if (present(status)) then
          status=-20
          return
       endif
       write(*,*) 'ERROR: non-finite integrated endpoint total'
       stop 1
    endif
    if (present(coeff_copy)) then
       if (.not.all(integrated_value_is_safe(coeff_copy))) then
          coeff=0d0
          coeff_copy=0d0
          if (present(status)) then
             status=-20
             return
          endif
          write(*,*) 'ERROR: non-finite integrated endpoint copy total'
          stop 1
       endif
    endif
  contains
    subroutine endpoint_failure(code,message)
      integer,intent(in) :: code
      character(len=*),intent(in) :: message
      coeff=0d0
      if (present(coeff_copy)) coeff_copy=0d0
      if (present(status)) then
         status=code
      else
         write(*,*) 'ERROR: '//trim(message),code
         stop 1
      endif
    end subroutine endpoint_failure
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
    if (.not.ieee_is_finite(mu_ren) .or. .not.all(ieee_is_finite(p))) then
       info=-20
       return
    endif
    if (mu_ren.le.0d0) then
       info=-20
       return
    endif
    sij=abs(2d0*minkowski_dot(p(:,emitter),p(:,spectator)))
    if (.not.ieee_is_finite(sij)) then
       info=-20
       return
    endif
    if (sij.le.0d0) then
       info=-20
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
       elseif (parton.eq.cs_parton_g .and. &
            light_quark_antiquark_pair(fi,fj)) then
          split=cs_massive_split_qqbar
       endif
       if (split.eq.0) then
          info=-11
          return
       endif
       q2=history%parent_mass**2+history%spectator_mass**2+sij
       if (.not.ieee_is_finite(q2)) then
          info=-20
          return
       endif
       if (q2.le.0d0) then
          info=-20
          return
       endif
       ! The explicit massive final--final formulae use Q_ik^2 as their
       ! dimensional scale, rather than 2 p_i.p_k.
       ell=log(4d0*cs_pi)+2d0*log(mu_ren)-log(q2)-0.5772156649015328606d0
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
       ell=log(4d0*cs_pi)+2d0*log(mu_ren)-log(sij)-0.5772156649015328606d0
       call cs_massive_fi_endpoint(history%parent_mass,sij,ell,&
            alpha_dipole(2),coeff,info)
    case (3)
       call parton_kind(fi,pa)
       call parton_kind(fp,pb)
       if (pa.eq.0 .or. pb.eq.0 .or. history%spectator_mass.le.0d0) then
          info=-13
          return
       endif
       ell=log(4d0*cs_pi)+2d0*log(mu_ren)-log(sij)-0.5772156649015328606d0
       call cs_massive_if_endpoint(pa,pb,history%spectator_mass,sij,ell,&
            integrated_dimensional_scheme,integrated_n_active_flavours,coeff,info)
    case default
       info=-14
    end select
    if (info.eq.0 .and. .not.all(ieee_is_finite(coeff))) then
       coeff=0d0
       info=-20
    endif
  end subroutine massive_history_endpoint

  subroutine integrated_beam(igroup,iproc,beam,z,hard_copy,xbj,mu_ren,mu_fac,alpha_s,&
       pterm,kterm,pterm_copy,kterm_copy,status)
    integer, intent(in) :: igroup,iproc,beam
    real(kind=8), intent(in) :: z,hard_copy(:),xbj(2),mu_ren,mu_fac,alpha_s
    real(kind=8), intent(out) :: pterm,kterm
    real(kind=8), intent(out),optional :: pterm_copy(:),kterm_copy(:)
    integer, intent(out),optional :: status
    type(cs_distribution) :: pk,kk,tilde_kernel,alpha_kernel
    type(cs_convolution_kernel) :: massive_kernel
    real(kind=8) :: fa,fb,gz,g1,other_pdf,regularised_p,regularised_k
    real(kind=8) :: pdf_beam(-6:7),pdf_other(-6:7),pdf_convolved(-6:7)
    real(kind=8) :: p_delta,sij,colour_log,colour_log_endpoint
    real(kind=8) :: history_weight,kernel_factor,primitive(-2:0)
    real(kind=8) :: fi_regular,fi_subtracted,szone,sx
    real(kind=8) :: colour_combination,p_contribution,k_contribution
    real(kind=8) :: alpha_factor,updated_value
    integer :: ih,icopy,a,b,other,info,pa,pb,fi,fj,fp,parton,spectator,emitter,pdf_info
    integer :: ncopy,allocation_status
    integer, allocatable :: leg_map(:)
    logical :: massive_history,massive_if_history,arithmetic_ok,history_valid
    character(len=256) :: allocation_message

    pterm=0d0
    kterm=0d0
    if (present(status)) status=0
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
    if (.not.integrated_value_is_safe(z) .or. &
         .not.all(integrated_value_is_safe(hard_copy)) .or. &
         .not.all(integrated_value_is_safe(xbj)) .or. &
         .not.integrated_value_is_safe(mu_ren) .or. &
         .not.integrated_value_is_safe(mu_fac) .or. &
         .not.integrated_value_is_safe(alpha_s)) then
       call beam_failure(-20,'invalid integrated-convolution input')
       return
    endif
    if (any(xbj.le.0d0) .or. any(xbj.gt.1d0) .or. mu_ren.le.0d0 .or. &
         mu_fac.le.0d0 .or. alpha_s.lt.0d0) then
       call beam_failure(-20,'invalid integrated-convolution input')
       return
    endif
    if (z.le.0d0 .or. z.ge.1d0) return
    if (beam.lt.1 .or. beam.gt.2) then
       call beam_failure(-4,'invalid beam index in integrated convolution')
       return
    endif
    if (.not.valid_born_process_state(igroup,iproc,ncopy,.true.)) then
       call beam_failure(-4,'invalid Born group/process state in integrated convolution')
       return
    endif
    if (size(hard_copy).ne.ncopy) then
       call beam_failure(-4,'incompatible integrated-convolution copy count')
       return
    endif
    if (.not.integrated_registry_is_ready()) then
       call beam_failure(-4,'integrated-history registry is not initialized')
       return
    endif
    if (.not.all(integrated_value_is_safe(pgl(igroup)%ps(1)%p))) then
       call beam_failure(-20,'invalid Born momentum in integrated convolution')
       return
    endif
    allocation_message=''
    allocate(leg_map(pgl(igroup)%next),stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       call beam_failure(-4,'could not allocate integrated-convolution workspace: '//&
            trim(allocation_message))
       return
    endif
    other=3-beam
    call evaluate_pdf_table(xbj(beam),mu_fac,pdf_beam,pdf_info)
    if (pdf_info.ne.0) then
       call beam_failure(pdf_info,'incoming-beam PDF evaluation failed')
       return
    endif
    call evaluate_pdf_table(xbj(other),mu_fac,pdf_other,pdf_info)
    if (pdf_info.ne.0) then
       call beam_failure(pdf_info,'other-beam PDF evaluation failed')
       return
    endif
    pdf_convolved=0d0
    if (z.ge.xbj(beam)) then
       call evaluate_pdf_table(xbj(beam)/z,mu_fac,pdf_convolved,pdf_info)
       if (pdf_info.ne.0) then
          call beam_failure(pdf_info,'convolved PDF evaluation failed')
          return
       endif
    endif
    if (.not.all(integrated_value_is_safe(pdf_beam)) .or. &
         .not.all(integrated_value_is_safe(pdf_other)) .or. &
         .not.all(integrated_value_is_safe(pdf_convolved))) then
       call beam_failure(-20,'non-finite PDF table in integrated convolution')
       return
    endif
    do icopy=1,size(hard_copy)
       if (hard_copy(icopy).eq.0d0) cycle
       b=pgl(igroup)%iden_processes(beam,icopy,iproc)
       call parton_kind(b,pb)
       if (pb.eq.0) cycle
       other_pdf=pdf_table_flavour(pdf_other,&
            pgl(igroup)%iden_processes(other,icopy,iproc))
       if (other_pdf.eq.0d0) cycle
       fb=pdf_table_flavour(pdf_beam,b)
       g1=fb
       do ih=1,n_integrated_histories
          if (integrated_history_list(ih)%incoming_leg.le.0) cycle
          if (.not.history_matches_copy(integrated_history_list(ih),igroup,iproc,icopy,leg_map)) cycle
          if (leg_map(integrated_history_list(ih)%incoming_leg).ne.beam) cycle
          call history_splitting_data(integrated_history_list(ih),fi,fj,fp,&
               history_weight,history_valid)
          if (.not.history_valid) then
             call beam_failure(-4,'invalid integrated-convolution history')
             return
          endif
          a=integrated_history_list(ih)%real_incoming_flavour
          call parton_kind(a,pa)
          if (pa.eq.0 .or. pb.eq.0) cycle
          fa=pdf_table_flavour(pdf_convolved,a)
          call safe_quotient(fa,z,gz,arithmetic_ok)
          if (.not.arithmetic_ok) then
             call beam_failure(-20,'unsafe convolved PDF quotient')
             return
          endif
          p_delta=0d0
          massive_history=history_has_mass(integrated_history_list(ih))
          massive_if_history=.false.
          if (massive_history) then
             emitter=leg_map(integrated_history_list(ih)%born_emitter)
             spectator=leg_map(integrated_history_list(ih)%born_spectator)
             select case (integrated_history_list(ih)%topology)
             case (2)
                szone=abs(2d0*minkowski_dot(pgl(igroup)%ps(1)%p(:,emitter),&
                     pgl(igroup)%ps(1)%p(:,beam)))
                if (.not.integrated_value_is_safe(szone)) then
                   call beam_failure(-20,'invalid FI invariant in integrated convolution')
                   return
                endif
                if (szone.le.0d0) then
                   call beam_failure(-20,'invalid FI invariant in integrated convolution')
                   return
                endif
                ! The real incoming momentum is p_a/z while the reduced
                ! Born momentum p_a is held fixed by the convolution.
                call safe_quotient(szone,z,sx,arithmetic_ok)
                if (.not.arithmetic_ok) then
                   call beam_failure(-20,'unsafe FI real invariant in integrated convolution')
                   return
                endif
                call cs_massive_fi_convolution(integrated_history_list(ih)%parent_mass,&
                     sx,szone,z,alpha_dipole(2),massive_kernel,info)
             case (3)
                szone=abs(2d0*minkowski_dot(pgl(igroup)%ps(1)%p(:,beam),&
                     pgl(igroup)%ps(1)%p(:,spectator)))
                if (.not.integrated_value_is_safe(szone)) then
                   call beam_failure(-20,'invalid IF invariant in integrated convolution')
                   return
                endif
                if (szone.le.0d0) then
                   call beam_failure(-20,'invalid IF invariant in integrated convolution')
                   return
                endif
                call safe_quotient(szone,z,sx,arithmetic_ok)
                if (.not.arithmetic_ok) then
                   call beam_failure(-20,'unsafe IF real invariant in integrated convolution')
                   return
                endif
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
                   call apply_p_distribution(pk,z,gz,g1,regularised_p,p_delta,arithmetic_ok)
                   if (.not.arithmetic_ok) info=-20
                   massive_if_history=.true.
                endif
             case default
                info=-20
             end select
             if (info.ne.0) then
                call beam_failure(info,'massive integrated convolution failed')
                return
             endif
             regularised_k=cs_apply_convolution(massive_kernel,gz,g1,info)
             if (info.eq.0 .and. fp.eq.99) then
                call safe_product(regularised_k,cs_u1_initial_factor,&
                     updated_value,arithmetic_ok)
                if (arithmetic_ok) then
                   regularised_k=updated_value
                else
                   info=-20
                endif
             endif
             if (info.ne.0) then
                call beam_failure(info,'unsafe massive integrated convolution value')
                return
             endif
             if (.not.massive_if_history) regularised_p=0d0
          elseif (integrated_history_list(ih)%topology.eq.2) then
             ! Final--initial histories have no universal P/Kbar term.  Their
             ! complete finite distribution, including the alpha=1 baseline,
             ! is a convolution on the initial spectator beam.
             call history_i_primitive(fi,fj,fp,integrated_history_list(ih)%topology,&
                  primitive,parton)
             if (parton.eq.0) cycle
             call cs_fi_distribution(primitive,z,alpha_dipole(2),&
                  fi_regular,fi_subtracted,info)
             if (info.ne.0) then
                call beam_failure(info,'final-initial integrated distribution failed')
                return
             endif
             regularised_p=0d0
             call apply_fi_distribution(fi_regular,fi_subtracted,gz,g1,&
                  regularised_k,arithmetic_ok)
             if (.not.arithmetic_ok) then
                call beam_failure(-20,'unsafe final-initial integrated distribution')
                return
             endif
          else
             call cs_ap_distribution(pa,pb,z,integrated_n_active_flavours,pk,info)
             if (info.ne.0) then
                call beam_failure(info,'integrated AP distribution failed')
                return
             endif
             call cs_kbar_distribution(pa,pb,z,integrated_n_active_flavours,kk,info)
             if (info.ne.0) then
                call beam_failure(info,'integrated Kbar distribution failed')
                return
             endif
             select case (integrated_history_list(ih)%topology)
             case (3)
                   call cs_if_tilde_distribution(pa,pb,z,integrated_n_active_flavours,&
                         tilde_kernel,info)
                if (info.ne.0) then
                   call beam_failure(info,'integrated IF tilde distribution failed')
                   return
                endif
                call cs_if_alpha_distribution(pa,pb,z,integrated_n_active_flavours,alpha_dipole(3),&
                     alpha_kernel,info)
             case (4)
                   call cs_ii_tilde_distribution(pa,pb,z,integrated_n_active_flavours,&
                         tilde_kernel,info)
                if (info.ne.0) then
                   call beam_failure(info,'integrated II tilde distribution failed')
                   return
                endif
                call cs_ii_alpha_distribution(pa,pb,z,integrated_n_active_flavours,alpha_dipole(4),&
                     alpha_kernel,info)
             case default
                tilde_kernel=cs_distribution()
                alpha_kernel=cs_distribution()
                info=0
             end select
             if (info.ne.0) then
                call beam_failure(info,'integrated alpha distribution failed')
                return
             endif
             kk%regular=kk%regular+tilde_kernel%regular+alpha_kernel%regular
             kk%plus_one=kk%plus_one+tilde_kernel%plus_one+alpha_kernel%plus_one
             kk%plus_log=kk%plus_log+tilde_kernel%plus_log+alpha_kernel%plus_log
             kk%plus_log_one=kk%plus_log_one+tilde_kernel%plus_log_one+&
                  alpha_kernel%plus_log_one
             kk%delta=kk%delta+tilde_kernel%delta+alpha_kernel%delta
             kernel_factor=1d0
             if (fp.eq.99) kernel_factor=cs_u1_initial_factor
             call cs_scale_distribution(pk,kernel_factor)
             call cs_scale_distribution(kk,kernel_factor)
             call apply_p_distribution(pk,z,gz,g1,regularised_p,p_delta,arithmetic_ok)
             if (.not.arithmetic_ok) then
                call beam_failure(-20,'unsafe integrated P distribution')
                return
             endif
             call apply_k_distribution(kk,z,gz,g1,regularised_k,arithmetic_ok)
             if (.not.arithmetic_ok) then
                call beam_failure(-20,'unsafe integrated K distribution')
                return
             endif
          endif
          spectator=leg_map(integrated_history_list(ih)%born_spectator)
          sij=abs(2d0*minkowski_dot(pgl(igroup)%ps(1)%p(:,beam),&
               pgl(igroup)%ps(1)%p(:,spectator)))
          if (.not.integrated_value_is_safe(sij)) then
             call beam_failure(-20,'invalid invariant in integrated convolution')
             return
          endif
          if (sij.le.0d0) then
             call beam_failure(-20,'invalid invariant in integrated convolution')
             return
          endif
          colour_log=2d0*log(mu_fac)-log(z)-log(sij)
          colour_log_endpoint=2d0*log(mu_fac)-log(sij)
          if (massive_if_history) then
             ! finiteif returns the complete mass-factorised convolution.
             ! Expose the same universal P contribution as the massless
             ! path and leave the compensating remainder in K, without
             ! changing their sum.
             call safe_linear_pair(colour_log,regularised_p,colour_log_endpoint,p_delta,&
                  colour_combination,arithmetic_ok)
             if (arithmetic_ok) call safe_sum(regularised_k,colour_combination,&
                  updated_value,arithmetic_ok)
             if (.not.arithmetic_ok) then
                call beam_failure(-20,'unsafe massive IF P/K separation')
                return
             endif
             regularised_k=updated_value
          endif
          if (.not.integrated_value_is_safe(regularised_p) .or. &
               .not.integrated_value_is_safe(regularised_k) .or. &
               .not.integrated_value_is_safe(p_delta) .or. &
               .not.integrated_value_is_safe(colour_log) .or. &
               .not.integrated_value_is_safe(colour_log_endpoint)) then
             call beam_failure(-20,'non-finite integrated convolution kernel')
             return
          endif
          call safe_linear_pair(colour_log,regularised_p,colour_log_endpoint,p_delta,&
               colour_combination,arithmetic_ok)
          if (arithmetic_ok) call safe_product_four(hard_copy(icopy),other_pdf,&
               history_weight,colour_combination,p_contribution,arithmetic_ok)
          if (arithmetic_ok) call safe_product_four(hard_copy(icopy),other_pdf,&
               history_weight,regularised_k,k_contribution,arithmetic_ok)
          if (.not.arithmetic_ok) then
             call beam_failure(-20,'unsafe integrated convolution contribution')
             return
          endif
          call safe_sum(pterm,-p_contribution,updated_value,arithmetic_ok)
          if (arithmetic_ok) pterm=updated_value
          if (arithmetic_ok) call safe_sum(kterm,k_contribution,updated_value,arithmetic_ok)
          if (arithmetic_ok) kterm=updated_value
          if (arithmetic_ok .and. present(pterm_copy)) then
             call safe_sum(pterm_copy(icopy),-p_contribution,updated_value,arithmetic_ok)
             if (arithmetic_ok) pterm_copy(icopy)=updated_value
          endif
          if (arithmetic_ok .and. present(kterm_copy)) then
             call safe_sum(kterm_copy(icopy),k_contribution,updated_value,arithmetic_ok)
             if (arithmetic_ok) kterm_copy(icopy)=updated_value
          endif
          if (.not.arithmetic_ok) then
             call beam_failure(-20,'unsafe integrated convolution accumulation')
             return
          endif
       enddo
    enddo
    call safe_quotient(alpha_s,2d0*cs_pi,alpha_factor,arithmetic_ok)
    if (.not.arithmetic_ok) then
       call beam_failure(-20,'unsafe integrated convolution coupling factor')
       return
    endif
    call safe_product(pterm,alpha_factor,updated_value,arithmetic_ok)
    if (arithmetic_ok) pterm=updated_value
    if (arithmetic_ok) call safe_product(kterm,alpha_factor,updated_value,arithmetic_ok)
    if (arithmetic_ok) kterm=updated_value
    if (.not.arithmetic_ok) then
       call beam_failure(-20,'non-finite integrated convolution total')
       return
    endif
    if (present(pterm_copy)) then
       do icopy=1,size(pterm_copy)
          call safe_product(pterm_copy(icopy),alpha_factor,updated_value,arithmetic_ok)
          if (.not.arithmetic_ok) then
             call beam_failure(-20,'unsafe integrated P-copy total')
             return
          endif
          pterm_copy(icopy)=updated_value
       enddo
    endif
    if (present(kterm_copy)) then
       do icopy=1,size(kterm_copy)
          call safe_product(kterm_copy(icopy),alpha_factor,updated_value,arithmetic_ok)
          if (.not.arithmetic_ok) then
             call beam_failure(-20,'unsafe integrated K-copy total')
             return
          endif
          kterm_copy(icopy)=updated_value
       enddo
    endif
  contains
    subroutine beam_failure(code,message)
      integer,intent(in) :: code
      character(len=*),intent(in) :: message
      pterm=0d0
      kterm=0d0
      if (present(pterm_copy)) pterm_copy=0d0
      if (present(kterm_copy)) kterm_copy=0d0
      if (present(status)) then
         status=code
      else
         write(*,*) 'ERROR: '//trim(message),code
         stop 1
      endif
    end subroutine beam_failure
  end subroutine integrated_beam

  subroutine history_i_primitive(fi,fj,fp,topology,coeff,parton)
    integer, intent(in) :: fi,fj,fp,topology
    real(kind=8), intent(out) :: coeff(-2:0)
    integer, intent(out) :: parton
    coeff=0d0
    parton=0
    if (light_quark_code(fp) .and. (fj.eq.21 .or. fj.eq.99)) then
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
    elseif ((fp.eq.21 .or. fp.eq.99) .and. &
         same_light_quark_flavour(fi,fj) .and. &
         ((topology.le.2 .and. light_quark_antiquark_pair(fi,fj)) .or. &
         topology.ge.3)) then
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
    if (light_quark_code(flavour)) then
       kind=cs_parton_q
    elseif (flavour.eq.21 .or. flavour.eq.99) then
       kind=cs_parton_g
    else
       kind=0
    endif
  end subroutine parton_kind

  elemental pure logical function integrated_value_is_safe(value)
    real(kind=8),intent(in) :: value
    integrated_value_is_safe=.false.
    if (.not.ieee_is_finite(value)) return
    if (abs(value).gt.integrated_value_limit) return
    if (value.ne.0d0 .and. abs(value).lt.tiny(1d0)) return
    integrated_value_is_safe=.true.
  end function integrated_value_is_safe

  pure subroutine safe_product(first,second,result,valid)
    real(kind=8),intent(in) :: first,second
    real(kind=8),intent(out) :: result
    logical,intent(out) :: valid
    result=0d0
    valid=.false.
    if (.not.integrated_value_is_safe(first) .or. &
         .not.integrated_value_is_safe(second)) return
    if (first.eq.0d0 .or. second.eq.0d0) then
       valid=.true.
       return
    endif
    if (abs(first).gt.integrated_value_limit/abs(second)) return
    if (abs(second).lt.1d0) then
       if (abs(first).lt.tiny(1d0)/abs(second)) return
    elseif (abs(first).lt.1d0) then
       if (abs(second).lt.tiny(1d0)/abs(first)) return
    endif
    result=first*second
    valid=integrated_value_is_safe(result)
    if (.not.valid) result=0d0
  end subroutine safe_product

  pure subroutine safe_sum(first,second,result,valid)
    real(kind=8),intent(in) :: first,second
    real(kind=8),intent(out) :: result
    logical,intent(out) :: valid
    result=0d0
    valid=.false.
    if (.not.integrated_value_is_safe(first) .or. &
         .not.integrated_value_is_safe(second)) return
    if ((first.ge.0d0 .and. second.ge.0d0) .or. &
         (first.lt.0d0 .and. second.lt.0d0)) then
       if (abs(first).gt.integrated_value_limit-abs(second)) return
    endif
    result=first+second
    valid=integrated_value_is_safe(result)
    if (.not.valid) result=0d0
  end subroutine safe_sum

  pure subroutine safe_product_three(first,second,third,result,valid)
    real(kind=8),intent(in) :: first,second,third
    real(kind=8),intent(out) :: result
    logical,intent(out) :: valid
    real(kind=8) :: partial
    call safe_product(first,second,partial,valid)
    if (valid) call safe_product(partial,third,result,valid)
    if (.not.valid) result=0d0
  end subroutine safe_product_three

  pure subroutine safe_product_four(first,second,third,fourth,result,valid)
    real(kind=8),intent(in) :: first,second,third,fourth
    real(kind=8),intent(out) :: result
    logical,intent(out) :: valid
    real(kind=8) :: partial
    call safe_product_three(first,second,third,partial,valid)
    if (valid) call safe_product(partial,fourth,result,valid)
    if (.not.valid) result=0d0
  end subroutine safe_product_four

  pure subroutine safe_quotient(numerator,denominator,result,valid)
    real(kind=8),intent(in) :: numerator,denominator
    real(kind=8),intent(out) :: result
    logical,intent(out) :: valid
    result=0d0
    valid=.false.
    if (.not.integrated_value_is_safe(numerator) .or. &
         .not.integrated_value_is_safe(denominator)) return
    if (denominator.eq.0d0) return
    if (numerator.eq.0d0) then
       valid=.true.
       return
    endif
    if (abs(numerator).gt.integrated_value_limit*abs(denominator)) return
    if (abs(denominator).gt.1d0) then
       if (abs(numerator).lt.tiny(1d0)*abs(denominator)) return
    endif
    result=numerator/denominator
    valid=integrated_value_is_safe(result)
    if (.not.valid) result=0d0
  end subroutine safe_quotient

  pure subroutine safe_accumulate_product(total,first,second,valid)
    real(kind=8),intent(inout) :: total
    real(kind=8),intent(in) :: first,second
    logical,intent(inout) :: valid
    real(kind=8) :: contribution,updated
    if (.not.valid) return
    call safe_product(first,second,contribution,valid)
    if (valid) call safe_sum(total,contribution,updated,valid)
    if (valid) total=updated
  end subroutine safe_accumulate_product

  pure subroutine safe_accumulate_product_three(total,first,second,third,valid)
    real(kind=8),intent(inout) :: total
    real(kind=8),intent(in) :: first,second,third
    logical,intent(inout) :: valid
    real(kind=8) :: contribution,updated
    if (.not.valid) return
    call safe_product_three(first,second,third,contribution,valid)
    if (valid) call safe_sum(total,contribution,updated,valid)
    if (valid) total=updated
  end subroutine safe_accumulate_product_three

  pure subroutine safe_linear_pair(first,second,third,fourth,result,valid)
    real(kind=8),intent(in) :: first,second,third,fourth
    real(kind=8),intent(out) :: result
    logical,intent(out) :: valid
    result=0d0
    valid=.true.
    call safe_accumulate_product(result,first,second,valid)
    call safe_accumulate_product(result,third,fourth,valid)
    if (.not.valid) result=0d0
  end subroutine safe_linear_pair

  pure subroutine apply_p_distribution(kernel,z,gz,g1,value,delta,valid)
    type(cs_distribution),intent(in) :: kernel
    real(kind=8),intent(in) :: z,gz,g1
    real(kind=8),intent(out) :: value,delta
    logical,intent(out) :: valid
    real(kind=8) :: difference,difference_quotient
    value=0d0
    delta=0d0
    difference=0d0
    difference_quotient=0d0
    valid=.true.
    call safe_sum(gz,-g1,difference,valid)
    if (valid) call safe_quotient(difference,1d0-z,difference_quotient,valid)
    call safe_accumulate_product(value,kernel%regular,gz,valid)
    call safe_accumulate_product(value,kernel%plus_one,difference_quotient,valid)
    if (valid) call safe_product(kernel%delta,g1,delta,valid)
    if (.not.valid) then
       value=0d0
       delta=0d0
    endif
  end subroutine apply_p_distribution

  pure subroutine apply_k_distribution(kernel,z,gz,g1,value,valid)
    type(cs_distribution),intent(in) :: kernel
    real(kind=8),intent(in) :: z,gz,g1
    real(kind=8),intent(out) :: value
    logical,intent(out) :: valid
    real(kind=8) :: difference,difference_quotient,log_ratio,log_one_minus
    value=0d0
    difference=0d0
    difference_quotient=0d0
    log_ratio=0d0
    log_one_minus=0d0
    valid=.true.
    if (.not.integrated_value_is_safe(z)) valid=.false.
    if (valid) then
       if (z.le.0d0 .or. z.ge.1d0) valid=.false.
    endif
    if (valid) call safe_sum(gz,-g1,difference,valid)
    if (valid) call safe_quotient(difference,1d0-z,difference_quotient,valid)
    if (valid) then
       log_one_minus=log(1d0-z)
       log_ratio=log_one_minus-log(z)
       valid=integrated_value_is_safe(log_one_minus) .and. &
            integrated_value_is_safe(log_ratio)
    endif
    call safe_accumulate_product(value,kernel%regular,gz,valid)
    call safe_accumulate_product(value,kernel%plus_one,difference_quotient,valid)
    call safe_accumulate_product_three(value,kernel%plus_log,log_ratio,&
         difference_quotient,valid)
    call safe_accumulate_product_three(value,kernel%plus_log_one,log_one_minus,&
         difference_quotient,valid)
    call safe_accumulate_product(value,kernel%delta,g1,valid)
    if (.not.valid) value=0d0
  end subroutine apply_k_distribution

  pure subroutine apply_fi_distribution(regular,subtracted,gz,g1,value,valid)
    real(kind=8),intent(in) :: regular,subtracted,gz,g1
    real(kind=8),intent(out) :: value
    logical,intent(out) :: valid
    real(kind=8) :: difference
    value=0d0
    difference=0d0
    valid=.true.
    call safe_sum(gz,-g1,difference,valid)
    call safe_accumulate_product(value,regular,gz,valid)
    call safe_accumulate_product(value,subtracted,difference,valid)
    if (.not.valid) value=0d0
  end subroutine apply_fi_distribution

  pure logical function history_has_mass(history)
    type(integrated_history), intent(in) :: history
    history_has_mass=history%emitter_mass.gt.0d0 .or. &
         history%unresolved_mass.gt.0d0 .or. &
         history%parent_mass.gt.0d0 .or. &
         history%spectator_mass.gt.0d0
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
       if (.not.same_integrated_splitting(candidate,histories(i))) cycle
       is_duplicate_history=.true.
       return
    enddo
  end function is_duplicate_history

  logical function same_integrated_splitting(first,second)
    type(integrated_history), intent(in) :: first,second
    integer :: fi1,fj1,fp1,fi2,fj2,fp2
    real(kind=8) :: weight1,weight2
    logical :: valid1,valid2

    same_integrated_splitting=.false.
    if (first%topology.ne.second%topology) return
    call history_splitting_data(first,fi1,fj1,fp1,weight1,valid1)
    if (.not.valid1) return
    call history_splitting_data(second,fi2,fj2,fp2,weight2,valid2)
    if (.not.valid2) return
    if (fp1.ne.fp2) return

    if ((fp1.eq.21 .or. fp1.eq.99) .and. first%topology.le.2 .and. &
         light_quark_antiquark_pair(fi1,fj1) .and. &
         light_quark_antiquark_pair(fi2,fj2)) then
       ! The two orientations of a final-state q-qbar pair have the same
       ! integrated kernel.
       same_integrated_splitting=same_light_quark_flavour(fi1,fi2)
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
       if (light_quark_code(t)) then
          if (.not.light_quark_code(a)) return
          if (sign(1,t).ne.sign(1,a)) return
          at=light_quark_magnitude(t)
          aa=light_quark_magnitude(a)
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
    elseif (light_quark_code(template_flavour)) then
       at=light_quark_magnitude(template_flavour)
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
    if (.not.valid_colour_order(process1,order1)) return
    if (.not.valid_colour_order(process2,order2)) return
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
    if (.not.valid_colour_order(process1,order1)) return
    if (.not.valid_colour_order(process2,order2)) return
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

  pure logical function valid_colour_order(process,order)
    integer, intent(in) :: process(:),order(:)
    integer :: label

    valid_colour_order=.false.
    if (size(process).lt.1 .or. size(order).ne.size(process)) return
    if (any(order.lt.1) .or. any(order.gt.size(process))) return
    do label=1,size(process)
       if (count(order.eq.label).ne.1) return
    enddo
    valid_colour_order=.true.
  end function valid_colour_order

end module integrated_dipoles
