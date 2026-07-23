module integrated_dipoles
  ! Inventory of integrated subtraction histories.  The inventory is built
  ! from the local dipoles themselves; it never independently enumerates
  ! splittings.  This is the central consistency guarantee between D and
  ! I/P/K at leading colour.
  use handling_processes
  use cs_dipole_mappings, only: cs_dipole_topology
  use cs_integrated_kernels
  use pdf_wrap, only: evaluate_pdf_flavour
  use common, only: alpha_dipole
  implicit none
  private

  type, public :: integrated_history
     integer :: real_group=0
     integer :: real_process=0
     integer :: local_dipole=0
     integer :: topology=0
     integer :: born_emitter=0
     integer :: born_spectator=0
     integer :: incoming_leg=0
     integer :: real_incoming_flavour=0
     integer :: born_incoming_flavour=0
     integer, allocatable :: real_flavours(:)
     integer, allocatable :: real_colour_order(:)
     integer, allocatable :: born_flavours(:)
     integer, allocatable :: born_colour_order(:)
  end type integrated_history

  type(integrated_history), allocatable, save, public :: integrated_history_list(:)
  integer, save, public :: n_integrated_histories=0
  integer, save, public :: integrated_dimensional_scheme=cs_scheme_hv

  public :: initialise_integrated_dipoles, history_matches_born
  public :: integrated_endpoint, integrated_beam

contains

  subroutine initialise_integrated_dipoles(nborn_groups,scheme_name)
    integer, intent(in) :: nborn_groups
    character(len=*), intent(in) :: scheme_name
    type(integrated_history), allocatable :: candidates(:),unique_histories(:)
    integer :: igroup,iproc,idip,max_histories,ncandidate,i,imatch

    if (trim(scheme_name).eq.'hv') then
       integrated_dimensional_scheme=cs_scheme_hv
    elseif (trim(scheme_name).eq.'fdh') then
       integrated_dimensional_scheme=cs_scheme_fdh
    else
       write(*,*) 'ERROR: unsupported integrated-dipole dimensional scheme: ',trim(scheme_name)
       stop 1
    endif

    max_histories=0
    do igroup=nborn_groups+1,ngroups
       if (.not.allocated(pgl(igroup)%dpl)) cycle
       do iproc=1,pgl(igroup)%nproc
          max_histories=max_histories+pgl(igroup)%dpl(iproc)%ndip
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
          do idip=1,pgl(igroup)%dpl(iproc)%ndip
             call ensure_massless_integrated_history(igroup,iproc,idip)
             imatch=find_born_match(pgl(igroup)%dpl(iproc)%dl(idip)%process_r,&
                  pgl(igroup)%dpl(iproc)%dl(idip)%reduced_color_order,nborn_groups)
             if (imatch.eq.0) then
                write(*,*) 'ERROR: integrated dipole has no corresponding Born process/order'
                write(*,*) ' real group/process/dipole:',igroup,iproc,idip
                write(*,*) ' reduced process:',pgl(igroup)%dpl(iproc)%dl(idip)%process_r
                write(*,*) ' reduced colour order:',pgl(igroup)%dpl(iproc)%dl(idip)%reduced_color_order
                call print_born_flavour_matches(pgl(igroup)%dpl(iproc)%dl(idip)%process_r,nborn_groups)
                stop 1
             endif
             if (is_duplicate_history(igroup,iproc,idip,candidates,ncandidate)) cycle
             ncandidate=ncandidate+1
             call fill_history(candidates(ncandidate),igroup,iproc,idip)
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

    write(*,'(a,i0,a,a)') 'Integrated subtraction registry: ',n_integrated_histories,&
         ' locally matched histories, scheme=',trim(scheme_name)
    write(99,'(a,i0,a,a)') 'Integrated subtraction registry: ',n_integrated_histories,&
         ' locally matched histories, scheme=',trim(scheme_name)
  end subroutine initialise_integrated_dipoles

  subroutine fill_history(history,igroup,iproc,idip)
    type(integrated_history), intent(out) :: history
    integer, intent(in) :: igroup,iproc,idip
    type(dipole) :: dip
    integer :: emitter,spectator

    dip=pgl(igroup)%dpl(iproc)%dl(idip)
    history%real_group=igroup
    history%real_process=iproc
    history%local_dipole=idip
    history%topology=cs_dipole_topology(dip%dip_ijk)
    history%born_emitter=dip%dip_r_ijk(1)
    history%born_spectator=dip%dip_r_ijk(2)
    history%real_flavours=pgl(igroup)%processes(:,iproc)
    history%real_colour_order=pgl(igroup)%color_orders(:,iproc)
    history%born_flavours=dip%process_r
    history%born_colour_order=dip%reduced_color_order
    emitter=dip%dip_ijk(1)
    if (emitter.le.2) then
       history%incoming_leg=emitter
       history%real_incoming_flavour=pgl(igroup)%processes(emitter,iproc)
       history%born_incoming_flavour=dip%process_r(dip%dip_r_ijk(1))
    elseif (dip%dip_ijk(3).le.2) then
       ! A final--initial dipole changes the momentum fraction of its
       ! initial-state spectator, so its non-endpoint alpha contribution
       ! belongs to that beam even though the emitter is final state.
       spectator=dip%dip_ijk(3)
       history%incoming_leg=spectator
       history%real_incoming_flavour=pgl(igroup)%processes(spectator,iproc)
       history%born_incoming_flavour=dip%process_r(dip%dip_r_ijk(2))
    endif
  end subroutine fill_history

  subroutine ensure_massless_integrated_history(igroup,iproc,idip)
    ! The kernels in cs_integrated_kernels are the massless Catani--Seymour
    ! formulae.  Local mappings support massive parents/spectators as well,
    ! but silently attaching the massless integrated formula to such a
    ! history gives a finite, physically incorrect result.
    integer, intent(in) :: igroup,iproc,idip
    type(dipole) :: dip
    real(kind=8) :: parent_mass,spectator_mass,mass_tolerance
    integer :: topology

    dip=pgl(igroup)%dpl(iproc)%dl(idip)
    parent_mass=abs(phys_model%get_mass(dip%dip_r_ijk_f(1)))
    spectator_mass=abs(phys_model%get_mass(dip%dip_r_ijk_f(2)))
    mass_tolerance=100d0*epsilon(1d0)*max(1d0,parent_mass,spectator_mass)
    if (parent_mass.le.mass_tolerance .and. spectator_mass.le.mass_tolerance) return

    topology=cs_dipole_topology(dip%dip_ijk)
    write(*,*) 'ERROR: massive integrated dipole history is not supported'
    write(*,*) ' real group/process/dipole/topology:',igroup,iproc,idip,topology
    write(*,*) ' reduced parent/spectator flavours:',dip%dip_r_ijk_f
    write(*,*) ' reduced parent/spectator masses:',parent_mass,spectator_mass
    write(99,*) 'ERROR: massive integrated dipole history is not supported'
    write(99,*) ' real group/process/dipole/topology:',igroup,iproc,idip,topology
    write(99,*) ' reduced parent/spectator flavours:',dip%dip_r_ijk_f
    write(99,*) ' reduced parent/spectator masses:',parent_mass,spectator_mass
    stop 1
  end subroutine ensure_massless_integrated_history

  integer function find_born_match(process,order,nborn_groups) result(nmatch)
    integer, intent(in) :: process(:),order(:),nborn_groups
    integer :: igroup,iproc,ip
    logical :: found
    nmatch=0
    found=.false.
    do igroup=1,nborn_groups
       if (pgl(igroup)%next.ne.size(process)) cycle
       do iproc=1,pgl(igroup)%nproc
          do ip=1,pgl(igroup)%iden_iproc(iproc)
             if (.not.same_flavour_pattern(process,pgl(igroup)%iden_processes(:,ip,iproc))) cycle
             if (.not.same_coloured_order(process,order,pgl(igroup)%iden_processes(:,ip,iproc),&
                  pgl(igroup)%color_orders(:,iproc))) cycle
             if (.not.found) then
                nmatch=1
                found=.true.
             endif
          enddo
       enddo
    enddo
  end function find_born_match

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
       if (.not.same_flavour_pattern(history%born_flavours,pgl(igroup)%iden_processes(:,ip,iproc))) cycle
       if (.not.same_coloured_order(history%born_flavours,history%born_colour_order,&
            pgl(igroup)%iden_processes(:,ip,iproc),pgl(igroup)%color_orders(:,iproc))) cycle
       history_matches_born=.true.
       return
    enddo
  end function history_matches_born

  logical function history_matches_copy(history,igroup,iproc,icopy)
    type(integrated_history), intent(in) :: history
    integer, intent(in) :: igroup,iproc,icopy
    history_matches_copy=.false.
    if (pgl(igroup)%next.ne.size(history%born_flavours)) return
    if (.not.same_flavour_pattern(history%born_flavours,&
         pgl(igroup)%iden_processes(:,icopy,iproc))) return
    history_matches_copy=same_coloured_order(history%born_flavours,history%born_colour_order,&
         pgl(igroup)%iden_processes(:,icopy,iproc),pgl(igroup)%color_orders(:,iproc))
  end function history_matches_copy

  subroutine integrated_endpoint(igroup,iproc,born_copy,p,mu,alpha_s,coeff)
    ! Laurent coefficients of I(eps), multiplied by the corresponding
    ! PDF-weighted Born contributions.  Every entry is inherited from an
    ! actual local-dipole history.
    integer, intent(in) :: igroup,iproc
    real(kind=8), intent(in) :: born_copy(:),p(0:,:)
    real(kind=8), intent(in) :: mu,alpha_s
    real(kind=8), intent(out) :: coeff(-2:0)
    real(kind=8) :: primitive(-2:0),expanded(-2:0),sij,ell,weight
    real(kind=8) :: endpoint_alpha
    integer :: ih,icopy,fi,fj,fp,parton,info
    integer :: emitter
    logical :: shifted(size(born_copy),pgl(igroup)%next)

    coeff=0d0
    shifted=.false.
    do icopy=1,size(born_copy)
       if (born_copy(icopy).eq.0d0) cycle
       do ih=1,n_integrated_histories
          if (.not.history_matches_copy(integrated_history_list(ih),igroup,iproc,icopy)) cycle
          fi=integrated_history_list(ih)%real_flavours(&
               pgl(integrated_history_list(ih)%real_group)%dpl(&
               integrated_history_list(ih)%real_process)%dl(&
               integrated_history_list(ih)%local_dipole)%dip_ijk(1))
          fj=integrated_history_list(ih)%real_flavours(&
               pgl(integrated_history_list(ih)%real_group)%dpl(&
               integrated_history_list(ih)%real_process)%dl(&
               integrated_history_list(ih)%local_dipole)%dip_ijk(2))
          fp=integrated_history_list(ih)%born_flavours(integrated_history_list(ih)%born_emitter)
          call history_i_primitive(fi,fj,fp,integrated_history_list(ih)%topology,&
               primitive,parton)
          if (parton.eq.0) cycle
          sij=abs(2d0*minkowski_dot(p(:,integrated_history_list(ih)%born_emitter),&
               p(:,integrated_history_list(ih)%born_spectator)))
          if (sij.le.tiny(1d0)) cycle
          ell=log(4d0*cs_pi*mu*mu/sij)-0.5772156649015328606d0
          expanded(-2)=primitive(-2)
          expanded(-1)=primitive(-1)+ell*primitive(-2)
          expanded(0)=primitive(0)+ell*primitive(-1)+&
               (0.5d0*ell*ell-cs_pi**2/12d0)*primitive(-2)
          weight=pgl(integrated_history_list(ih)%real_group)%dpl(&
               integrated_history_list(ih)%real_process)%dl(&
               integrated_history_list(ih)%local_dipole)%lc_weight
          endpoint_alpha=0d0
          select case (integrated_history_list(ih)%topology)
          case (1)
             call cs_ff_alpha_endpoint(alpha_dipole(1),primitive,endpoint_alpha,info)
          case (2)
             call cs_fi_alpha_endpoint(alpha_dipole(2),primitive,endpoint_alpha,info)
          case (3,4)
             ! Initial--final and initial--initial restrictions have no
             ! alpha-dependent endpoint.  Their complete finite changes
             ! are applied as beam convolutions below.
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
          coeff=coeff+born_copy(icopy)*weight*expanded*alpha_s/(2d0*cs_pi)
          emitter=integrated_history_list(ih)%born_emitter
          if (integrated_dimensional_scheme.eq.cs_scheme_fdh .and. &
               phys_model%get_mass(integrated_history_list(ih)%born_flavours(emitter)).eq.0d0 .and. &
               .not.shifted(icopy,emitter)) then
             coeff(0)=coeff(0)+born_copy(icopy)*cs_fdh_endpoint_shift(parton)*&
                  alpha_s/(2d0*cs_pi)
             shifted(icopy,emitter)=.true.
          endif
       enddo
    enddo
  end subroutine integrated_endpoint

  subroutine integrated_beam(igroup,iproc,beam,z,hard_copy,xbj,mu_fac,alpha_s,pterm,kterm)
    integer, intent(in) :: igroup,iproc,beam
    real(kind=8), intent(in) :: z,hard_copy(:),xbj(2),mu_fac,alpha_s
    real(kind=8), intent(out) :: pterm,kterm
    type(cs_distribution) :: pk,kk,alpha_kernel
    real(kind=8) :: fa,fb,gz,g1,other_pdf,regularised_p,regularised_k
    real(kind=8) :: sij,colour_log,history_weight,primitive(-2:0)
    real(kind=8) :: fi_regular,fi_subtracted
    integer :: ih,icopy,a,b,other,info,pa,pb,fi,fj,fp,parton
    logical :: flavour_map_ok
    logical :: seen(n_integrated_histories)

    pterm=0d0
    kterm=0d0
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
          if (integrated_history_list(ih)%incoming_leg.ne.beam) cycle
          if (.not.history_matches_copy(integrated_history_list(ih),igroup,iproc,icopy)) cycle
          call map_history_flavour(integrated_history_list(ih)%born_flavours,&
               pgl(igroup)%iden_processes(:,icopy,iproc),&
               integrated_history_list(ih)%real_incoming_flavour,a,flavour_map_ok)
          if (.not.flavour_map_ok) cycle
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
          if (integrated_history_list(ih)%topology.eq.2) then
             ! Final--initial histories have no universal P/Kbar term.  A
             ! restricted local dipole does, however, produce a finite
             ! convolution on the initial spectator beam.
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
             call cs_fi_alpha_terms(primitive,z,alpha_dipole(2),&
                  fi_regular,fi_subtracted,info)
             if (info.ne.0) cycle
             regularised_p=0d0
             regularised_k=fi_regular*gz+fi_subtracted*(gz-g1)
          else
             call cs_ap_distribution(pa,pb,z,5,pk,info)
             if (info.ne.0) cycle
             call cs_kbar_distribution(pa,pb,z,5,kk,info)
             if (info.ne.0) cycle
             regularised_p=pk%regular*gz+pk%plus_one*(gz-g1)/(1d0-z)+pk%delta*g1
             regularised_k=kk%regular*gz+kk%plus_one*(gz-g1)/(1d0-z)+&
                  kk%plus_log*log((1d0-z)/z)*(gz-g1)/(1d0-z)+kk%delta*g1
             select case (integrated_history_list(ih)%topology)
             case (3)
                call cs_if_alpha_distribution(pa,pb,z,alpha_dipole(3),&
                     alpha_kernel,info)
             case (4)
                call cs_ii_alpha_distribution(pa,pb,z,alpha_dipole(4),&
                     alpha_kernel,info)
             case default
                alpha_kernel=cs_distribution()
                info=0
             end select
             if (info.ne.0) cycle
             regularised_k=regularised_k+alpha_kernel%regular*gz
          endif
          sij=abs(2d0*minkowski_dot(pgl(igroup)%ps(1)%p(:,beam),&
               pgl(igroup)%ps(1)%p(:,integrated_history_list(ih)%born_spectator)))
          colour_log=0d0
          if (sij.gt.tiny(1d0)) colour_log=log(mu_fac*mu_fac/(z*sij))
          history_weight=pgl(integrated_history_list(ih)%real_group)%dpl(&
               integrated_history_list(ih)%real_process)%dl(&
               integrated_history_list(ih)%local_dipole)%lc_weight
          pterm=pterm-hard_copy(icopy)*other_pdf*history_weight*colour_log*regularised_p
          kterm=kterm+hard_copy(icopy)*other_pdf*history_weight*regularised_k
          seen(ih)=.true.
       enddo
    enddo
    pterm=pterm*alpha_s/(2d0*cs_pi)
    kterm=kterm*alpha_s/(2d0*cs_pi)
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
       if (topology.eq.1 .or. topology.eq.2) then
          ! FF/FI local g -> gg histories contain one ordered side of the
          ! symmetric kernel.  Attaching the full primitive to both the
          ! (i,j) and (j,i) histories doubles their integrated alpha change.
          call cs_i_gg_ordered(coeff)
       else
          call cs_i_gg(coeff)
       endif
       parton=cs_parton_g
    elseif ((fp.eq.21 .or. fp.eq.99) .and. abs(fi).le.6 .and. fi.eq.-fj) then
       call cs_i_qqbar(coeff)
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

  pure real(kind=8) function minkowski_dot(a,b)
    real(kind=8), intent(in) :: a(0:3),b(0:3)
    minkowski_dot=a(0)*b(0)-dot_product(a(1:3),b(1:3))
  end function minkowski_dot

  logical function is_duplicate_history(igroup,iproc,idip,histories,nhistory)
    integer, intent(in) :: igroup,iproc,idip,nhistory
    type(integrated_history), intent(in) :: histories(:)
    integer :: i
    type(dipole) :: dip

    dip=pgl(igroup)%dpl(iproc)%dl(idip)
    is_duplicate_history=.false.
    do i=1,nhistory
       if (.not.same_flavour_pattern(pgl(igroup)%processes(:,iproc),histories(i)%real_flavours)) cycle
       ! Colour-singlet permutations are distinct phase-space channels, not
       ! distinct subtraction histories.  Comparing the full order retained
       ! one copy of an otherwise identical integrated dipole for every
       ! singlet sampling order.
       if (.not.same_coloured_order(pgl(igroup)%processes(:,iproc),&
            pgl(igroup)%color_orders(:,iproc),histories(i)%real_flavours,&
            histories(i)%real_colour_order)) cycle
       if (.not.all(dip%dip_ijk.eq.&
            pgl(histories(i)%real_group)%dpl(histories(i)%real_process)%dl(histories(i)%local_dipole)%dip_ijk)) cycle
       if (.not.same_flavour_pattern(dip%process_r,histories(i)%born_flavours)) cycle
       if (.not.same_coloured_order(dip%process_r,dip%reduced_color_order,&
            histories(i)%born_flavours,histories(i)%born_colour_order)) cycle
       is_duplicate_history=.true.
       return
    enddo
  end function is_duplicate_history

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

end module integrated_dipoles
