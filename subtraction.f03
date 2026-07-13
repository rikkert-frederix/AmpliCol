module subtraction
  use handling_processes
  use particles
  use amplitude_QCD_mod, only: use_real_gluons
  use common, only: alpha_dipole
  implicit none
  integer :: n
  integer,parameter :: alpha_status_excluded=-100
  real(kind=8),parameter :: alpha_cut_tolerance=1d-12
  private
  public :: initialise_subtraction, test_limits_integrand, generate_limit_phase_space_point
  public :: print_limit_failure_fractions, compute_the_amps, square_the_amps
  public :: evaluate_real_dipoles, evaluate_recycled_born, compute_the_dipole_amps, square_the_dipole_amps
  public :: initialise_recycling_history_weights
contains
  subroutine initialise_subtraction(igroup,iamp)
    implicit none
    integer,intent(in) :: iamp,igroup
    integer :: ipart,is_dipole,ipart_l,ipart_r,idip,ichan,ichannel,nchannels,ncandidates,max_candidates
    integer,allocatable :: candidate_dipoles(:,:)
    logical,allocatable :: candidate_reverse(:)
    n=pgl(igroup)%next
    if (.not.allocated(pgl(igroup)%dpl)) allocate(pgl(igroup)%dpl(pgl(igroup)%nproc))
    if (allocated(pgl(igroup)%dpl(iamp)%dl)) then
       call finalize_dipole_set(pgl(igroup)%dpl(iamp))
    endif
    ! A matrix element can be integrated by several phase-space channels when
    ! colour singlets are permuted.  Its dipoles must cover the union of the
    ! colour-adjacent limits of all those channels, not just igroup's order.
    nchannels=pgl(igroup)%multichan%unique_channelgroup_list(0, &
         pgl(igroup)%multichan%map_proc_to_channelgroup(iamp))
    max_candidates=2*(pgl(igroup)%next-2)*nchannels
    allocate(candidate_dipoles(3,max_candidates),candidate_reverse(max_candidates))
    ncandidates=0
    do ichannel=1,nchannels
       ichan=pgl(igroup)%multichan%unique_channel_list( &
            pgl(igroup)%multichan%unique_channelgroup_list(ichannel, &
            pgl(igroup)%multichan%map_proc_to_channelgroup(iamp)))
       do ipart=3,pgl(igroup)%next
          call is_valid_dipole(ipart,pgl(igroup)%processes(:,iamp),pgl(ichan)%phase_space_orders(:), &
               is_dipole,ipart_l,ipart_r)
          if (btest(is_dipole,0)) then
             if (dipole_has_valid_born(igroup,iamp,ipart_l,ipart,ipart_r,.true.)) then
                call add_dipole_candidate(ipart_l,ipart,ipart_r,.true., &
                     candidate_dipoles,candidate_reverse,ncandidates)
             endif
          endif
          if (btest(is_dipole,1)) then
             if (dipole_has_valid_born(igroup,iamp,ipart_r,ipart,ipart_l,.false.)) then
                call add_dipole_candidate(ipart_r,ipart,ipart_l,.false., &
                     candidate_dipoles,candidate_reverse,ncandidates)
             endif
          endif
       enddo
    enddo
    pgl(igroup)%dpl(iamp)%ndip=ncandidates
    allocate(pgl(igroup)%dpl(iamp)%dl(pgl(igroup)%dpl(iamp)%ndip))

    do idip=1,pgl(igroup)%dpl(iamp)%ndip
       call fill_dipole(pgl(igroup)%dpl(iamp)%dl(idip),pgl(igroup)%processes(1:n,iamp), &
            candidate_dipoles(1,idip),candidate_dipoles(2,idip),candidate_dipoles(3,idip),candidate_reverse(idip))
    enddo
    do idip=1,pgl(igroup)%dpl(iamp)%ndip
       allocate(pgl(igroup)%dpl(iamp)%dl(idip)%reduced_color_order(n-1))
       call build_reduced_color_order(pgl(igroup)%color_orders(1:n,iamp), &
            pgl(igroup)%dpl(iamp)%dl(idip)%dip_ijk(2),pgl(igroup)%dpl(iamp)%dl(idip)%process_r, &
            pgl(igroup)%dpl(iamp)%dl(idip)%reduced_color_order)
       pgl(igroup)%dpl(iamp)%dl(idip)%col_fac=lc_colour_factor(pgl(igroup)%dpl(iamp)%dl(idip)%process_r)
       call pgl(igroup)%dpl(iamp)%dl(idip)%amp%init(1,n-1,1,pgl(igroup)%dpl(iamp)%dl(idip)%process_r,&
            pgl(igroup)%spin(0:3,pgl(igroup)%dpl(iamp)%dl(idip)%dip_map(1:n-1)), &
            pgl(igroup)%dpl(iamp)%dl(idip)%reduced_color_order,&
            phys_model)
       allocate(pgl(igroup)%dpl(iamp)%dl(idip)%p_mapped(0:3,n-1))
       call initialise_rho_lookup(pgl(igroup)%dpl(iamp)%dl(idip))
    enddo
!!$    call print_dipoles(pgl(igroup)%processes(:,iamp),pgl(igroup)%color_orders(:,iamp),pgl(igroup)%dpl(iamp)%dl)
  end subroutine initialise_subtraction

  subroutine initialise_recycling_history_weights()
    ! A physical Born subprocess can be reached from several real crossing
    ! channels and flavour replicas.  Partition those histories according to
    ! their process-file multiplicity within each reduced process and coloured
    ! ordering.  The remaining local CS-projection partition is evaluated
    ! point by point in evaluate_recycled_born.
    integer :: ig,ip,idip,iiden,jg,jp,jdip,jiden
    integer,allocatable :: process_r(:),process_r_other(:)
    logical :: history_matches,target_supported
    real(kind=8) :: history_norm
    do ig=1,ngroups
       do ip=1,pgl(ig)%nproc
          do idip=1,pgl(ig)%dpl(ip)%ndip
             associate(dip => pgl(ig)%dpl(ip)%dl(idip))
             if (dip%dip_r_ijk_f(1).eq.99) cycle
             target_supported=dipole_has_unit_recycling_measure(dip)
             if (target_supported) allocate(dip%recycling_history_weight(pgl(ig)%iden_iproc(ip)))
             allocate(process_r(size(dip%process_r)))
             do iiden=1,pgl(ig)%iden_iproc(ip)
                call reduced_process_for_dipole(dip,pgl(ig)%iden_processes(:,iiden,ip),process_r)
                history_norm=0d0
                do jg=1,ngroups
                   do jp=1,pgl(jg)%nproc
                      if (jg.ne.representative_recycling_group(jg,jp)) cycle
                      do jiden=1,pgl(jg)%iden_iproc(jp)
                         if (pgl(jg)%idenCOandMAPfactor(jiden,jp).eq.0d0) cycle
                         history_matches=.false.
                         do jdip=1,pgl(jg)%dpl(jp)%ndip
                            associate(dip_other => pgl(jg)%dpl(jp)%dl(jdip))
                            if (.not.dipole_has_unit_recycling_measure(dip_other)) cycle
                            allocate(process_r_other(size(dip_other%process_r)))
                            call reduced_process_for_dipole(dip_other,&
                                 pgl(jg)%iden_processes(:,jiden,jp),process_r_other)
                            history_matches=same_recycled_born_key(dip,process_r,dip_other,process_r_other)
                            deallocate(process_r_other)
                            if (history_matches) exit
                            end associate
                         enddo
                         if (history_matches) history_norm=history_norm+ &
                              abs(pgl(jg)%recycling_history_multiplicity(jiden,jp))
                      enddo
                   enddo
                enddo
                if (history_norm.eq.0d0) then
                   write (*,*) 'ERROR: no massless recycling history for reduced Born process'
                   write (*,*) process_r
                   write (*,*) 'Massive radiation-volume normalization is not implemented for Born recycling.'
                   stop 1
                endif
                if (target_supported) dip%recycling_history_weight(iiden)=1d0/history_norm
             enddo
             deallocate(process_r)
             end associate
          enddo
       enddo
    enddo
  end subroutine initialise_recycling_history_weights

  integer function representative_recycling_group(igroup,iproc)
    ! The same physical subprocess is present in every phase-space channel
    ! listed for its multichannel group.  Count it only in the lowest-numbered
    ! explicit channel; compute_multichannel_weight combines the other copies.
    integer,intent(in) :: igroup,iproc
    integer :: ichannelgroup,ichannel,nchannels,iunique,igroup_other
    logical :: contains_self
    ichannelgroup=pgl(igroup)%multichan%map_proc_to_channelgroup(iproc)
    nchannels=pgl(igroup)%multichan%unique_channelgroup_list(0,ichannelgroup)
    if (nchannels.lt.1) then
       write (*,*) 'ERROR: empty recycling multichannel group',igroup,iproc
       stop 1
    endif
    representative_recycling_group=huge(0)
    contains_self=.false.
    do ichannel=1,nchannels
       iunique=pgl(igroup)%multichan%unique_channelgroup_list(ichannel,ichannelgroup)
       igroup_other=pgl(igroup)%multichan%unique_channel_list(iunique)
       if (igroup_other.lt.1 .or. igroup_other.gt.ngroups) then
          write (*,*) 'ERROR: invalid recycling multichannel group',igroup_other
          stop 1
       endif
       representative_recycling_group=min(representative_recycling_group,igroup_other)
       contains_self=contains_self .or. igroup_other.eq.igroup
    enddo
    if (.not.contains_self) then
       write (*,*) 'ERROR: recycling multichannel group omits source group',igroup,iproc
       stop 1
    endif
  end function representative_recycling_group

  logical function same_recycled_born_key(dip_a,process_a,dip_b,process_b)
    ! Colour-singlet permutations are phase-space multichannel copies, not
    ! distinct Born histories.  Preserve the order of coloured legs while
    ! removing every singlet from the comparison.
    type(dipole),intent(in) :: dip_a,dip_b
    integer,intent(in) :: process_a(:),process_b(:)
    integer :: i,leg,ncol_a,ncol_b
    integer,dimension(size(process_a)) :: coloured_order_a,coloured_order_b
    same_recycled_born_key=.false.
    if (size(process_a).ne.size(process_b)) return
    if (.not.all(process_a.eq.process_b)) return
    if (size(dip_a%reduced_color_order).ne.size(dip_b%reduced_color_order)) return
    ncol_a=0
    do i=1,size(dip_a%reduced_color_order)
       leg=dip_a%reduced_color_order(i)
       if (phys_model%is_singlet(process_a(leg))) cycle
       ncol_a=ncol_a+1
       coloured_order_a(ncol_a)=leg
    enddo
    ncol_b=0
    do i=1,size(dip_b%reduced_color_order)
       leg=dip_b%reduced_color_order(i)
       if (phys_model%is_singlet(process_b(leg))) cycle
       ncol_b=ncol_b+1
       coloured_order_b(ncol_b)=leg
    enddo
    if (ncol_a.ne.ncol_b) return
    same_recycled_born_key=all(coloured_order_a(1:ncol_a).eq.coloured_order_b(1:ncol_b))
  end function same_recycled_born_key

  logical function dipole_has_unit_recycling_measure(dip)
    ! The current push-back weight normalises the massless CS radiation cube.
    ! Massive emitters or spectators have mass-dependent radiation bounds and
    ! therefore cannot be used as independent Born-recycling histories.
    type(dipole),intent(in) :: dip
    dipole_has_unit_recycling_measure=.false.
    if (dip%dip_r_ijk_f(1).eq.99) return
    if (phys_model%get_mass(dip%dip_ijk_f(1)).ne.0d0) return
    if (phys_model%get_mass(dip%dip_ijk_f(2)).ne.0d0) return
    if (phys_model%get_mass(dip%dip_ijk_f(3)).ne.0d0) return
    if (phys_model%get_mass(dip%dip_r_ijk_f(1)).ne.0d0) return
    dipole_has_unit_recycling_measure=.true.
  end function dipole_has_unit_recycling_measure

  subroutine add_dipole_candidate(dip_i,dip_j,dip_k,reverse,candidates,reverses,ncandidates)
    implicit none
    integer,intent(in) :: dip_i,dip_j,dip_k
    logical,intent(in) :: reverse
    integer,intent(inout) :: candidates(:,:),ncandidates
    logical,intent(inout) :: reverses(:)
    integer :: i
    do i=1,ncandidates
       if (all(candidates(:,i).eq.[dip_i,dip_j,dip_k])) return
    enddo
    ncandidates=ncandidates+1
    candidates(:,ncandidates)=[dip_i,dip_j,dip_k]
    reverses(ncandidates)=reverse
  end subroutine add_dipole_candidate

  subroutine initialise_rho_lookup(dip)
    implicit none
    type(dipole),intent(inout) :: dip
    integer :: ih1,ih2,ij,nmatch
    logical :: is_vector_parent

    ij=dip%dip_r_ijk(1)
    is_vector_parent=phys_model%is_gluon(dip%dip_r_ijk_f(1))
    nmatch=0
    do ih1=1,dip%amp%n_amps
       do ih2=1,dip%amp%n_amps
          if (.not.same_rho_spins(dip%amp%spins(:,1,ih1),dip%amp%spins(:,1,ih2),ij, &
               is_vector_parent)) cycle
          nmatch=nmatch+1
       enddo
    enddo
    allocate(dip%rho_lookup_ih1(nmatch),dip%rho_lookup_ih2(nmatch))
    nmatch=0
    do ih1=1,dip%amp%n_amps
       do ih2=1,dip%amp%n_amps
          if (.not.same_rho_spins(dip%amp%spins(:,1,ih1),dip%amp%spins(:,1,ih2),ij, &
               is_vector_parent)) cycle
          nmatch=nmatch+1
          dip%rho_lookup_ih1(nmatch)=ih1
          dip%rho_lookup_ih2(nmatch)=ih2
       enddo
    enddo
    dip%rho_lookup_upper=.not.is_vector_parent
    dip%rho_hermitian_checked=.not.is_vector_parent
  end subroutine initialise_rho_lookup

  logical function same_rho_spins(spins1,spins2,ij,is_vector_parent)
    implicit none
    integer,intent(in) :: spins1(:),spins2(:),ij
    logical,intent(in) :: is_vector_parent
    integer :: i
    same_rho_spins=.false.
    if (.not.is_vector_parent .and. spins1(ij).ne.spins2(ij)) return
    do i=1,size(spins1)
       if (i.eq.ij) cycle
       if (spins1(i).ne.spins2(i)) return
    enddo
    same_rho_spins=.true.
  end function same_rho_spins

  logical function dipole_has_valid_born(igroup,iamp,dip_i,dip_j,dip_k,reverse)
    use amplitude_QCD_mod, only: amplitude_QCD
    implicit none
    integer,intent(in) :: igroup,iamp,dip_i,dip_j,dip_k
    logical,intent(in) :: reverse
    type(dipole) :: trial_dipole
    type(amplitude_QCD) :: trial_amplitude
    logical :: valid

    call fill_dipole(trial_dipole,pgl(igroup)%processes(1:n,iamp),dip_i,dip_j,dip_k,reverse)
    allocate(trial_dipole%reduced_color_order(n-1))
    call build_reduced_color_order(pgl(igroup)%color_orders(1:n,iamp),trial_dipole%dip_ijk(2), &
         trial_dipole%process_r,trial_dipole%reduced_color_order)
    call trial_amplitude%init(1,n-1,1,trial_dipole%process_r, &
         pgl(igroup)%spin(0:3,trial_dipole%dip_map(1:n-1)),trial_dipole%reduced_color_order,phys_model,valid)
    dipole_has_valid_born=valid
  end function dipole_has_valid_born
  subroutine print_dipoles(process,order,dips)
    implicit none
    integer,dimension(*),intent(in) :: process,order
    type(dipole),dimension(:),intent(in) :: dips
    integer :: idip
    write (*,*) 'process',process(1:n)
    write (*,*) 'color-order',order(1:n)
    write (*,*) 'Need',size(dips),'dipoles'
    do idip=1,size(dips)
       write (*,*) '------------------'
       write (*,*) 'dipole',idip
       write (*,*) 'i,j,k',dips(idip)%dip_ijk
       write (*,*) 'i,j,k',dips(idip)%dip_ijk_f
       write (*,*) 'process reduced',dips(idip)%process_r
       write (*,*) 'color order reduced',dips(idip)%reduced_color_order
    enddo
    write (*,*) '------------------'
    write (*,*) ''
    write (*,*) ''
  end subroutine print_dipoles
  subroutine fill_dipole(dip,process,dip_i,dip_j,dip_k,reverse)
    implicit none
    integer,dimension(*),intent(in) :: process
    type(dipole),intent(inout) :: dip
    integer,intent(in) :: dip_i,dip_j,dip_k
    logical,intent(in) :: reverse
    integer :: ipart,i
    dip%dip_ijk(1:3)=[dip_i,dip_j,dip_k]
    dip%dip_ijk_f(1:3)=[process(dip_i),process(dip_j),process(dip_k)]
    if (dip_i.gt.2) dip%dipole_type=ibset(dip%dipole_type,0)
    if (dip_k.gt.2) dip%dipole_type=ibset(dip%dipole_type,1)
    ! reduced process and dipole info
    if (dip_j .lt. dip_i) then
       dip%dip_r_ijk(1)=dip_i-1
    else
       dip%dip_r_ijk(1)=dip_i
    endif
    if (dip_j .lt. dip_k) then
       dip%dip_r_ijk(2)=dip_k-1
    else
       dip%dip_r_ijk(2)=dip_k
    endif
    if (phys_model%is_gluon(process(dip_j))) then
       dip%dip_r_ijk_f(1)=dip%dip_ijk_f(1)
    elseif (phys_model%is_gluon(process(dip_i))) then
       if (btest(dip%dipole_type,0)) then
          write (*,*) 'error in dipoles: emitter is a final-state gluon and '// &
               'emitted is a quark'
          write (*,*) dip%dip_ijk
          write (*,*) dip%dipole_type
          write (*,*) dip%dip_ijk_f
          stop 1
       endif
       dip%dip_r_ijk_f(1)=phys_model%get_antipart(dip%dip_ijk_f(2))
    else
       dip%dip_r_ijk_f(1)=combined_gluon_type(dip_i,process(dip_i),dip_j,process(dip_j),reverse)
    endif
    if (dip%dip_r_ijk_f(1).eq.21 .or. dip%dip_r_ijk_f(1).eq.99) dip%lc_weight=0.5d0
    dip%dip_r_ijk_f(2)=dip%dip_ijk_f(3)
    allocate(dip%process_r(n-1))
    allocate(dip%dip_map(n-1))
    i=0
    do ipart=1,n
       if (ipart.eq.dip_j) cycle
       i=i+1
       dip%dip_map(i)=ipart
       if (ipart.eq.dip_i) then
          dip%process_r(i)=dip%dip_r_ijk_f(1)
       elseif(ipart.eq.dip_k) then
          dip%process_r(i)=dip%dip_r_ijk_f(2)
       else
          dip%process_r(i)=process(ipart)
       endif
    enddo
  end subroutine fill_dipole
  integer function lc_colour_factor(process)
    implicit none
    integer,dimension(:),intent(in) :: process
    integer :: i,ifac
    real(kind=8) :: fac
    fac=0d0
    do i=1,size(process)
       if (process(i).eq.21) then
          fac=fac+1d0
       elseif (phys_model%is_quark(process(i)) .or. phys_model%is_antiquark(process(i))) then
          fac=fac+0.5d0
       endif
    enddo
    ifac=nint(fac)
    if (dble(ifac).ne.fac) then
       write (*,*) 'There is some issue with the reduced LC colour factor computation: ',ifac,fac
       stop 1
    endif
    lc_colour_factor=3**ifac
  end function lc_colour_factor

  integer function combined_gluon_type(dip_i,part_i,dip_j,part_j,reverse)
    implicit none
    integer,intent(in) :: dip_i,part_i,dip_j,part_j
    logical,intent(in) :: reverse
    if (((dip_i.le.2 .and. phys_model%is_quark(part_i)) .or. &
         (dip_i.gt.2 .and. phys_model%is_antiquark(part_i))) .and. &
        ((dip_j.le.2 .and. phys_model%is_antiquark(part_j)) .or. &
         (dip_j.gt.2 .and. phys_model%is_quark(part_j)))) then
       if (reverse) then
          combined_gluon_type=21
       else
          combined_gluon_type=99
       endif
    elseif (((dip_i.le.2 .and. phys_model%is_antiquark(part_i)) .or. &
             (dip_i.gt.2 .and. phys_model%is_quark(part_i))) .and. &
            ((dip_j.le.2 .and. phys_model%is_quark(part_j)) .or. &
             (dip_j.gt.2 .and. phys_model%is_antiquark(part_j)))) then
       if (reverse) then
          combined_gluon_type=99
       else
          combined_gluon_type=21
       endif
    else
       write (*,*) 'ERROR: cannot infer combined gluon type from dipole pair'
       write (*,*) part_i,part_j
       stop 1
    endif
  end function combined_gluon_type
  subroutine build_reduced_color_order(parent_order,removed_pos,process,order)
    implicit none
    integer,dimension(:),intent(in) :: parent_order
    integer,intent(in) :: removed_pos
    integer,dimension(:),intent(in) :: process
    integer,dimension(:),intent(out) :: order
    integer :: i,ipos,insert
    logical :: found_quark,found_antiquark
    ipos=0
    do i=1,n
       if (parent_order(i).eq.removed_pos) cycle
       ipos=ipos+1
       if (parent_order(i).gt.removed_pos) then
          order(ipos)=parent_order(i)-1
       else
          order(ipos)=parent_order(i)
       endif
    enddo
    found_quark=.false.
    do i=1,n-1
       if ((order(i).le.2 .and. phys_model%is_antiquark(process(order(i)))) .or. &
            (order(i).gt.2 .and. phys_model%is_quark(process(order(i))))) then
          ! found quark to start colour order with
          order=[order(i:),order(:i-1)]
          found_quark=.true.
          exit
       endif
    enddo
    if (.not.found_quark) return
    found_antiquark=.false.
    do i=n-1,1,-1
       if ((order(i).le.2 .and. phys_model%is_quark(process(order(i)))) .or. &
            (order(i).gt.2 .and. phys_model%is_antiquark(process(order(i))))) then
          ! found the last antiquark
          insert=i
          found_antiquark=.true.
          exit
       endif
    enddo
    if (.not.found_antiquark) return
    i=n-1
    do
       if (phys_model%is_singlet(process(order(i)))) then
          if (i.gt.insert) then
             order=[order(1:insert-1),order(i),order(insert:i-1),order(i+1:n-1)]
             insert=insert+1
             i=i+1
          elseif (i.lt.insert) then
             order=[order(1:i-1),order(i+1:insert-1),order(i),order(insert:n-1)]
             insert=insert-1 ! new position
          endif
       endif
       i=i-1
       if (i.eq.0) exit
    enddo
  end subroutine build_reduced_color_order
  subroutine is_valid_dipole(ipart,process,order,is_dipole,ipart_l,ipart_r)
    ! Checks whether the nearest coloured neighbours of ipart in the colour
    ! order form a valid dipole that could have radiated particle ipart.
    implicit none
    integer,intent(in) :: ipart
    integer,dimension(n),intent(in) :: process,order
    integer,intent(out) :: is_dipole,ipart_l,ipart_r
    integer :: i,ipart_l_pos,ipart_r_pos
    is_dipole=0
    do i=1,n
       if (order(i).eq.ipart) exit
    enddo
    if (phys_model%get_mass(process(ipart)).ne.0d0) return
    if (phys_model%get_colour_rep(process(ipart)).eq.1) return
    ! Colour singlets do not interrupt a colour line.  Find the nearest
    ! coloured neighbours rather than only the immediately adjacent legs.
    ipart_l_pos=i
    do
       ipart_l_pos=mod(n+ipart_l_pos-2,n)+1
       if (phys_model%get_colour_rep(process(order(ipart_l_pos))).ne.1) exit
       if (ipart_l_pos.eq.i) return
    enddo
    ipart_r_pos=i
    do
       ipart_r_pos=mod(ipart_r_pos,n)+1
       if (phys_model%get_colour_rep(process(order(ipart_r_pos))).ne.1) exit
       if (ipart_r_pos.eq.i) return
    enddo
    ipart_l=order(ipart_l_pos)
    ipart_r=order(ipart_r_pos)

    if (phys_model%is_gluon(process(ipart))) then
       is_dipole=3 ! both left and right can be emitters
    else ! must be a quark
       if (ipart_l.gt.2) then
          if (phys_model%get_antipart(process(ipart)).eq.process(ipart_l)) is_dipole=ibset(is_dipole,0)
       else
          if (process(ipart).eq.process(ipart_l)) is_dipole=ibset(is_dipole,0)
          if (phys_model%is_gluon(process(ipart_l))) is_dipole=ibset(is_dipole,0)
       endif
       if (ipart_r.gt.2) then
          if (phys_model%get_antipart(process(ipart)).eq.process(ipart_r)) is_dipole=ibset(is_dipole,1)
       else
          if (process(ipart).eq.process(ipart_r)) is_dipole=ibset(is_dipole,1)
          if (phys_model%is_gluon(process(ipart_r))) is_dipole=ibset(is_dipole,1)
       endif
    endif
  end subroutine is_valid_dipole
  subroutine test_limits_integrand(ichan,iint,limit_point,soft_fail,soft_tested,&
       collinear_fail,collinear_tested,use_amplitude_library)
    use phase_space_module
    use cuts, only: pass_real_subtracted_cuts
    implicit none
    integer,intent(in) :: ichan
    integer,intent(in) :: iint
    integer,intent(in) :: limit_point
    integer,intent(inout) :: soft_fail(:),soft_tested(:)
    integer,intent(inout) :: collinear_fail(:,:),collinear_tested(:,:)
    logical,intent(in) :: use_amplitude_library
    integer,parameter :: nsteps=11
    real(kind=8),parameter :: limit_tolerance=1d-2,nonsingular_growth_tolerance=5d-2
    integer :: i,j,k,status,nselected,nmatched,nalpha_selected
    real(kind=8) :: lambda,mass(pgl(ichan)%next-2),amp2_dip,ratio
    real(kind=8) :: ratios(nsteps),residuals(nsteps)
    real(kind=8) :: lambdas(nsteps),amp_values(nsteps),dip_values(nsteps)
    integer :: mapping_status(nsteps)
    logical :: sequence_ok,no_dipoles,soft_no_dipoles,valid_values(nsteps),mapping_failed,is_gluon,real_pass
    real(kind=8),dimension(0:3,pgl(ichan)%next) :: p_save

    p_save=pgl(ichan)%ps(1)%p

    do i=3,pgl(ichan)%next
       mass(i-2)=phys_model%get_mass(pgl(ichan)%processes(i,iint))
    enddo

    do i=3,pgl(ichan)%next
       if (phys_model%get_mass(pgl(ichan)%processes(i,iint)).gt.0d0) cycle
       is_gluon=phys_model%is_gluon(pgl(ichan)%processes(i,iint))
       ratios=0d0
       residuals=-1d0
       lambdas=0d0
       amp_values=0d0
       dip_values=0d0
       mapping_status=0
       valid_values=.false.
       mapping_failed=.false.
       soft_no_dipoles=.false.
       do k = 0, nsteps-1
          lambda = 10.0_dp**(-real(k,kind=8)/2d0)
          lambdas(k+1)=lambda
          call soft_deform_event(pgl(ichan)%next-2, mass, p_save, i, lambda, pgl(ichan)%ps(1)%p, status)
          if (status.ne.0) then
             mapping_status(k+1)=status
             call write_limit_failure('soft',ichan,iint,limit_point,i,0,nsteps,lambdas,amp_values,dip_values,ratios,&
                  residuals,valid_values,mapping_status)
             soft_fail(i)=soft_fail(i)+1
             mapping_failed=.true.
             exit
          endif
          call compute_the_amps(iint,ichan,use_amplitude_library)
          call square_the_amps(iint,ichan)
          amp_values(k+1)=pgl(ichan)%amp2(1)
          real_pass=pass_real_subtracted_cuts(pgl(ichan),iint)
          if (.not.real_pass) amp_values(k+1)=0d0
          if (is_gluon) then
             call compute_the_dipole_amps(iint,ichan)
             call square_the_dipole_amps(iint,ichan,amp2_dip,i,nselected=nselected,nmatched=nmatched,&
                  nalpha_selected=nalpha_selected)
             if (nmatched.eq.0) then
                soft_no_dipoles=.true.
                valid_values(k+1)=finite_nonsingular_value(amp_values(k+1))
             elseif (nselected.eq.0) then
                if (.not.real_pass) then
                   ! Both measurement functions vanish, so this point cannot
                   ! contribute to the cut real-subtracted observable.
                   ratios(k+1)=1d0
                   residuals(k+1)=0d0
                   valid_values(k+1)=.true.
                elseif (nalpha_selected.eq.0) then
                   mapping_status(k+1)=alpha_status_excluded
                else
                   mapping_status(k+1)=-101
                endif
                cycle
             else
                dip_values(k+1)=amp2_dip
                valid_values(k+1)=finite_limit_values(amp_values(k+1),amp2_dip,ratio,residuals(k+1))
                if (valid_values(k+1)) then
                   ratios(k+1)=ratio
                endif
             endif
          else
             valid_values(k+1)=finite_nonsingular_value(amp_values(k+1))
          endif
       enddo
       soft_tested(i)=soft_tested(i)+1
       if (mapping_failed) cycle
       if (is_gluon .and. soft_no_dipoles) then
          call assess_nonsingular_limit_sequence(amp_values,valid_values,nsteps,nonsingular_growth_tolerance,sequence_ok)
       elseif (is_gluon) then
          call assess_limit_sequence(ratios,residuals,nsteps,limit_tolerance,sequence_ok)
       else
          call assess_integrable_soft_limit_sequence(amp_values,lambdas,valid_values,nsteps,&
               nonsingular_growth_tolerance,sequence_ok)
       endif
       if (.not.sequence_ok) then
          soft_fail(i)=soft_fail(i)+1
          if (is_gluon .and. soft_no_dipoles) then
             call write_limit_failure('nonsingular soft',ichan,iint,limit_point,i,0,nsteps,lambdas,amp_values,dip_values,&
                  ratios,residuals,valid_values,mapping_status)
          elseif (is_gluon) then
             call write_limit_failure('soft',ichan,iint,limit_point,i,0,nsteps,lambdas,amp_values,dip_values,ratios,&
                  residuals,valid_values,mapping_status)
          else
             call write_limit_failure('integrable soft',ichan,iint,limit_point,i,0,nsteps,lambdas,amp_values,dip_values,ratios,&
                  residuals,valid_values,mapping_status)
          endif
       endif
    enddo

    do i=1,pgl(ichan)%next-1
       do j=max(3,i+1),pgl(ichan)%next
          if (phys_model%get_mass(pgl(ichan)%processes(i,iint)).gt.0d0 .or. &
               phys_model%get_mass(pgl(ichan)%processes(j,iint)).gt.0d0) cycle
          ratios=0d0
          residuals=-1d0
          lambdas=0d0
          amp_values=0d0
          dip_values=0d0
          mapping_status=0
          valid_values=.false.
          no_dipoles=.false.
          mapping_failed=.false.
          do k = 0, nsteps-1
             lambda = 10.0_dp**(-real(k, dp)/2d0)
             lambdas(k+1)=lambda
             call collinear_deform_event(pgl(ichan)%next-2, mass, p_save, i, j, lambda, pgl(ichan)%ps(1)%p, status)
             if (status.ne.0) then
                mapping_status(k+1)=status
                call write_limit_failure('collinear',ichan,iint,limit_point,i,j,nsteps,lambdas,amp_values,dip_values,ratios,&
                     residuals,valid_values,mapping_status)
                collinear_fail(i,j)=collinear_fail(i,j)+1
                mapping_failed=.true.
                exit
             endif
             call compute_the_amps(iint,ichan,use_amplitude_library)
             call square_the_amps(iint,ichan)
             amp_values(k+1)=pgl(ichan)%amp2(1)
             real_pass=pass_real_subtracted_cuts(pgl(ichan),iint)
             if (.not.real_pass) amp_values(k+1)=0d0
             if (.not.no_dipoles) then
                call compute_the_dipole_amps(iint,ichan)
                call square_the_dipole_amps(iint,ichan,amp2_dip,icol1=i,icol2=j,&
                     nselected=nselected,nmatched=nmatched,nalpha_selected=nalpha_selected)
                if (nmatched.eq.0) no_dipoles=.true.
                if (nmatched.gt.0 .and. nselected.eq.0) then
                   if (.not.real_pass) then
                      ratios(k+1)=1d0
                      residuals(k+1)=0d0
                      valid_values(k+1)=.true.
                      cycle
                   elseif (nalpha_selected.eq.0) then
                      mapping_status(k+1)=alpha_status_excluded
                   else
                      mapping_status(k+1)=-101
                   endif
                   cycle
                endif
             endif
             if (no_dipoles) then
                valid_values(k+1)=finite_nonsingular_value(amp_values(k+1))
             else
                dip_values(k+1)=amp2_dip
                valid_values(k+1)=finite_limit_values(amp_values(k+1),amp2_dip,ratio,residuals(k+1))
                if (valid_values(k+1)) then
                   ratios(k+1)=ratio
                endif
             endif
          enddo
          collinear_tested(i,j)=collinear_tested(i,j)+1
          if (mapping_failed) cycle
          if (no_dipoles) then
             call assess_nonsingular_limit_sequence(amp_values,valid_values,nsteps,nonsingular_growth_tolerance,sequence_ok)
          else
             call assess_limit_sequence(ratios,residuals,nsteps,limit_tolerance,sequence_ok)
          endif
          if (.not.sequence_ok) then
             collinear_fail(i,j)=collinear_fail(i,j)+1
             if (no_dipoles) then
                call write_limit_failure('nonsingular collinear',ichan,iint,limit_point,i,j,nsteps,lambdas,amp_values,dip_values,&
                     ratios,residuals,valid_values,mapping_status)
             else
                call write_limit_failure('collinear',ichan,iint,limit_point,i,j,nsteps,lambdas,amp_values,dip_values,ratios,&
                     residuals,valid_values,mapping_status)
             endif
          endif
       enddo
    enddo
  end subroutine test_limits_integrand

  subroutine generate_limit_phase_space_point(ichan)
    use phase_space_module
    use cuts
    implicit none
    integer,intent(in) :: ichan
    integer :: ix
    real(kind=8),external :: ran2
    do
       do ix=1,size(pgl(ichan)%ps(1)%x)
          pgl(ichan)%ps(1)%x(ix)=ran2()
       enddo
       call pgl(ichan)%phase_space%generate_momenta(pgl(ichan)%ps(1))
       if (pgl(ichan)%ps(1)%jac.gt.0d0 .and. pass_real_subtracted_cuts(pgl(ichan),1)) then
          return
       endif
   enddo
  end subroutine generate_limit_phase_space_point

  logical function finite_limit_values(amp2,amp2_dip,ratio,residual)
    use, intrinsic :: ieee_arithmetic
    implicit none
    real(kind=8),intent(in) :: amp2,amp2_dip
    real(kind=8),intent(out) :: ratio,residual
    real(kind=8) :: scale,norm_amp2,norm_dip
    finite_limit_values=.false.
    ratio=0d0
    residual=-1d0
    if (.not.ieee_is_finite(amp2) .or. .not.ieee_is_finite(amp2_dip)) return
    scale=max(abs(amp2),abs(amp2_dip))
    if (scale.eq.0d0 .or. .not.ieee_is_finite(scale)) return
    norm_amp2=amp2/scale
    norm_dip=amp2_dip/scale
    residual=abs(norm_amp2-norm_dip)
    if (abs(norm_dip).lt.sqrt(tiny(1d0)) .or. .not.ieee_is_finite(residual)) then
       residual=-1d0
       return
    endif
    ratio=norm_amp2/norm_dip
    if (.not.ieee_is_finite(ratio)) return
    if (ratio.le.0d0) then
       residual=-1d0
       return
    endif
    finite_limit_values=.true.
  end function finite_limit_values

  logical function finite_nonsingular_value(amp2)
    use, intrinsic :: ieee_arithmetic
    implicit none
    real(kind=8),intent(in) :: amp2
    finite_nonsingular_value=ieee_is_finite(amp2)
  end function finite_nonsingular_value

  subroutine write_limit_failure(limit_name,ichan,iint,limit_point,ileg1,ileg2,nsteps,lambdas,amp_values,dip_values,&
       ratios,residuals,valid_values,mapping_status)
    implicit none
    character(len=*),intent(in) :: limit_name
    integer,intent(in) :: ichan,iint,limit_point,ileg1,ileg2,nsteps
    real(kind=8),intent(in) :: lambdas(nsteps),amp_values(nsteps),dip_values(nsteps)
    real(kind=8),intent(in) :: ratios(nsteps),residuals(nsteps)
    logical,intent(in) :: valid_values(nsteps)
    integer,intent(in) :: mapping_status(nsteps)
    integer :: k
    write (100,'(a)') '------------------------------------------------------------'
    if (index(limit_name,'soft').ne.0) then
       write (100,*) 'FAILED ',trim(limit_name),' limit: point ',limit_point,' channel ',ichan,' iint ',iint,' leg ',ileg1
    else
       write (100,*) 'FAILED ',trim(limit_name),' limit: point ',limit_point,' channel ',ichan,' iint ',iint,' legs ',ileg1,ileg2
    endif
    write (100,'(a)') '# step lambda matrix_element dipole ratio residual valid mapping_status'
    write (100,'(a)') '# mapping_status -100: all matching real dipoles excluded by alpha'
    write (100,'(a)') '# mapping_status -101: all alpha-active matching dipoles fail mapped cuts'
    do k=1,nsteps
       write (100,'(i4,1x,5(es24.16,1x),l1,1x,i0)') k-1,lambdas(k),amp_values(k),&
            dip_values(k),ratios(k),residuals(k),valid_values(k),mapping_status(k)
    enddo
    write (100,'(a)') ''
  end subroutine write_limit_failure

  subroutine print_limit_failure_fractions(soft_fail,soft_tested,collinear_fail,collinear_tested,&
       failure_threshold,all_ok)
    implicit none
    integer,intent(in) :: soft_fail(:,:,:),soft_tested(:,:,:)
    integer,intent(in) :: collinear_fail(:,:,:,:),collinear_tested(:,:,:,:)
    integer,intent(in) :: failure_threshold
    logical,intent(out) :: all_ok
    integer :: ichan,i,j,iint
    real(kind=8) :: fraction
    character(len=4) :: result
    all_ok=.true.
    write (*,'(a)') 'Limit failure fractions:'
    write (99,'(a)') 'Limit failure fractions:'
    do ichan=1,ngroups
       do iint=1,pgl(ichan)%nproc
          do i=3,pgl(ichan)%next
             if (phys_model%get_mass(pgl(ichan)%processes(i,iint)).gt.0d0) then
                write (*,'(2x,"channel ",i0," integral ",i0," soft leg ",i0,": SKIP (massive)")') ichan,iint,i
                write (99,'(2x,"channel ",i0," integral ",i0," soft leg ",i0,": SKIP (massive)")') ichan,iint,i
                cycle
             endif
             fraction=100d0*dble(soft_fail(ichan,iint,i))/dble(soft_tested(ichan,iint,i))
             if (soft_fail(ichan,iint,i).ge.failure_threshold) then
                result='FAIL'
                all_ok=.false.
             else
                result='PASS'
             endif
             if (phys_model%is_gluon(pgl(ichan)%processes(i,iint))) then
                write (*,'(2x,"channel ",i0," integral ",i0," soft leg ",i0,": ",i0,"/",i0, &
                     " failed (",f5.1,"%)",t86,a4)') ichan,iint,i,soft_fail(ichan,iint,i), &
                     soft_tested(ichan,iint,i),fraction,result
                write (99,'(2x,"channel ",i0," integral ",i0," soft leg ",i0,": ",i0,"/",i0, &
                     " failed (",f5.1,"%)",t86,a4)') ichan,iint,i,soft_fail(ichan,iint,i), &
                     soft_tested(ichan,iint,i),fraction,result
             else
                write (*,'(2x,"channel ",i0," integral ",i0," soft leg ",i0," (integrable): ",i0,"/",i0, &
                     " failed (",f5.1,"%)",t86,a4)') ichan,iint,i,soft_fail(ichan,iint,i), &
                     soft_tested(ichan,iint,i),fraction,result
                write (99,'(2x,"channel ",i0," integral ",i0," soft leg ",i0," (integrable): ",i0,"/",i0, &
                     " failed (",f5.1,"%)",t86,a4)') ichan,iint,i,soft_fail(ichan,iint,i), &
                     soft_tested(ichan,iint,i),fraction,result
             endif
          enddo
          do i=1,pgl(ichan)%next-1
             do j=max(3,i+1),pgl(ichan)%next
                if (phys_model%get_mass(pgl(ichan)%processes(i,iint)).gt.0d0 .or. &
                     phys_model%get_mass(pgl(ichan)%processes(j,iint)).gt.0d0) then
                   write (*,'(2x,"channel ",i0," integral ",i0," collinear legs ",i0,"/",i0, &
                        ": SKIP (massive)")') ichan,iint,i,j
                   write (99,'(2x,"channel ",i0," integral ",i0," collinear legs ",i0,"/",i0, &
                        ": SKIP (massive)")') ichan,iint,i,j
                   cycle
                endif
                fraction=100d0*dble(collinear_fail(ichan,iint,i,j))/&
                     dble(collinear_tested(ichan,iint,i,j))
                if (collinear_fail(ichan,iint,i,j).ge.failure_threshold) then
                   result='FAIL'
                   all_ok=.false.
                else
                   result='PASS'
                endif
                if (has_collinear_dipoles(iint,ichan,i,j)) then
                   write (*,'(2x,"channel ",i0," integral ",i0," collinear legs ",i0,"/",i0,": ", &
                        i0,"/",i0," failed (",f5.1,"%)",t86,a4)') ichan,iint,i,j, &
                        collinear_fail(ichan,iint,i,j),collinear_tested(ichan,iint,i,j),fraction,result
                   write (99,'(2x,"channel ",i0," integral ",i0," collinear legs ",i0,"/",i0,": ", &
                        i0,"/",i0," failed (",f5.1,"%)",t86,a4)') ichan,iint,i,j, &
                        collinear_fail(ichan,iint,i,j),collinear_tested(ichan,iint,i,j),fraction,result
                else
                   write (*,'(2x,"channel ",i0," integral ",i0," collinear legs ",i0,"/",i0," (finite): ", &
                        i0,"/",i0," failed (",f5.1,"%)",t86,a4)') ichan,iint,i,j, &
                        collinear_fail(ichan,iint,i,j),collinear_tested(ichan,iint,i,j),fraction,result
                   write (99,'(2x,"channel ",i0," integral ",i0," collinear legs ",i0,"/",i0," (finite): ", &
                        i0,"/",i0," failed (",f5.1,"%)",t86,a4)') ichan,iint,i,j, &
                        collinear_fail(ichan,iint,i,j),collinear_tested(ichan,iint,i,j),fraction,result
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine print_limit_failure_fractions

  subroutine assess_limit_sequence(ratios,residuals,nsteps,tolerance,passed)
    implicit none
    integer,intent(in) :: nsteps
    real(kind=8),intent(in) :: ratios(nsteps),residuals(nsteps),tolerance
    logical,intent(out) :: passed
    integer :: k
    real(kind=8) :: window_error,window_change
    passed=.false.
    do k=1,nsteps-2
       if (any(residuals(k:k+2).lt.0d0) .or. any(residuals(k:k+2).gt.1d0)) cycle
       window_error=maxval(residuals(k:k+2))
       window_change=maxval(abs(ratios(k+1:k+2)-ratios(k:k+1)))
       if (window_error.le.tolerance .and. window_change.le.2d0*tolerance) then
          passed=.true.
       elseif (residuals(k+2).le.tolerance .and. residuals(k).le.3d0*tolerance .and. &
            residuals(k).gt.residuals(k+1) .and. residuals(k+1).gt.residuals(k+2) .and. &
            window_change.le.3d0*tolerance) then
          ! At the smallest usable lambda the asymptotic window may contain
          ! only two points before roundoff spoils the next point.  Accept a
          ! monotonic approach when its last point is already within the
          ! strict tolerance.
          passed=.true.
       endif
    enddo
    do k=2,nsteps-1
       if (any(residuals(k-1:k+1).lt.0d0) .or. any(residuals(k-1:k+1).gt.1d0)) cycle
       ! A resolved local minimum is sufficient when smaller lambdas are
       ! already dominated by roundoff.  Its neighbours must still be close
       ! enough to one to exclude an accidental crossing.
       if (residuals(k).le.tolerance .and. residuals(k).le.residuals(k-1) .and. &
            residuals(k).le.residuals(k+1) .and. maxval(residuals(k-1:k+1)).le.5d0*tolerance) then
          passed=.true.
       endif
    enddo
    if (residuals(nsteps-1).ge.0d0 .and. residuals(nsteps).ge.0d0 .and. &
         residuals(nsteps-1).le.3d0*tolerance .and. residuals(nsteps).le.tolerance .and. &
         residuals(nsteps-1).gt.residuals(nsteps)) then
       ! At the end of the scan only two usable points may remain before the
       ! next deformation would be numerically unreliable.
       passed=.true.
    endif
  end subroutine assess_limit_sequence

  subroutine assess_nonsingular_limit_sequence(amp_values,valid_values,nsteps,growth_tolerance,passed)
    implicit none
    integer,intent(in) :: nsteps
    real(kind=8),intent(in) :: amp_values(nsteps),growth_tolerance
    logical,intent(in) :: valid_values(nsteps)
    logical,intent(out) :: passed
    integer :: k
    real(kind=8) :: previous_scale
    passed=.false.
    do k=1,nsteps-2
       if (.not.all(valid_values(k:k+2))) cycle
       previous_scale=max(abs(amp_values(k)),abs(amp_values(k+1)),tiny(1d0))
       ! A nonsingular matrix element must level off as lambda decreases.
       ! A small tolerance avoids rejecting harmless phase-space variation.
       if (abs(amp_values(k+2)).le.(1d0+growth_tolerance)*previous_scale) then
          passed=.true.
       endif
    enddo
  end subroutine assess_nonsingular_limit_sequence

  subroutine assess_integrable_soft_limit_sequence(amp_values,lambdas,valid_values,nsteps,growth_tolerance,passed)
    implicit none
    integer,intent(in) :: nsteps
    real(kind=8),intent(in) :: amp_values(nsteps),lambdas(nsteps),growth_tolerance
    logical,intent(in) :: valid_values(nsteps)
    logical,intent(out) :: passed
    integer :: k
    real(kind=8) :: scaled_values(nsteps),previous_scale
    scaled_values=lambdas*abs(amp_values)
    passed=.false.
    do k=1,nsteps-2
       if (.not.all(valid_values(k:k+2))) cycle
       previous_scale=max(scaled_values(k),scaled_values(k+1),tiny(1d0))
       ! Soft-fermion amplitudes may grow as 1/lambda while remaining
       ! integrable.  Test the phase-space-weighted matrix element instead.
       if (scaled_values(k+2).le.(1d0+growth_tolerance)*previous_scale) then
          passed=.true.
       endif
    enddo
  end subroutine assess_integrable_soft_limit_sequence

  subroutine compute_the_dipole_amps(iint,ichan)
    use cs_dipole_mappings
    use common, only: alpha_dipole
    use cuts, only: pass_mapped_dipole_cuts
    implicit none
    integer,intent(in) :: iint,ichan
    integer :: idip,info,cut_info,topology
    real(kind=8) :: cut_variable
    real(kind=8),dimension(0:3,pgl(ichan)%next-1) :: ps_mapped
    real(kind=8),dimension(pgl(ichan)%next) :: mass_real
    integer,dimension(pgl(ichan)%next-1) :: hel_mapped
    integer :: ipart
    do ipart=1,pgl(ichan)%next
       mass_real(ipart)=phys_model%get_mass(pgl(ichan)%processes(ipart,iint))
    enddo
    do idip=1,pgl(ichan)%dpl(iint)%ndip
       pgl(ichan)%dpl(iint)%dl(idip)%active=.false.
       pgl(ichan)%dpl(iint)%dl(idip)%alpha_active=.false.
       pgl(ichan)%dpl(iint)%dl(idip)%passes_cuts=.false.
       call cs_dipole_cut_variable(pgl(ichan)%ps(1)%p,pgl(ichan)%dpl(iint)%dl(idip)%dip_ijk,&
            mass_real,phys_model%get_mass(pgl(ichan)%dpl(iint)%dl(idip)%dip_r_ijk_f(1)),cut_variable,cut_info)
       if (cut_info.ne.0) then
          write (*,*) 'error computing cs alpha variable',cut_info
          stop 1
       endif
       topology=cs_dipole_topology(pgl(ichan)%dpl(iint)%dl(idip)%dip_ijk)
       if (cut_variable.gt.alpha_dipole(topology)+alpha_cut_tolerance) cycle
       pgl(ichan)%dpl(iint)%dl(idip)%alpha_active=.true.
       call cs_map(pgl(ichan)%ps(1)%p,pgl(ichan)%dpl(iint)%dl(idip)%dip_ijk,ps_mapped,info, &
            mass_real=mass_real, &
            mass_parent=phys_model%get_mass(pgl(ichan)%dpl(iint)%dl(idip)%dip_r_ijk_f(1)))
       pgl(ichan)%dpl(iint)%dl(idip)%p_mapped=ps_mapped
       pgl(ichan)%dpl(iint)%dl(idip)%p_mapped_ij(0:3)= &
            ps_mapped(0:3,pgl(ichan)%dpl(iint)%dl(idip)%dip_r_ijk(1))
       if (info.ne.0) then
          write (*,*) 'error in cs momentum mapping',info
          stop 1
       endif
       pgl(ichan)%dpl(iint)%dl(idip)%passes_cuts=pass_mapped_dipole_cuts(ps_mapped,&
            pgl(ichan)%dpl(iint)%dl(idip)%process_r)
       if (.not.pgl(ichan)%dpl(iint)%dl(idip)%passes_cuts) cycle
       pgl(ichan)%dpl(iint)%dl(idip)%active=.true.
       hel_mapped=pgl(ichan)%hel(pgl(ichan)%dpl(iint)%dl(idip)%dip_map(1:pgl(ichan)%next-1))
       call pgl(ichan)%dpl(iint)%dl(idip)%amp%evaluate(pgl(ichan)%next-1,ps_mapped,&
            hel_mapped,read_proc_from_file,phys_model)
    enddo
  end subroutine compute_the_dipole_amps

  subroutine evaluate_real_dipoles(iint,ichan,amp2_dip)
    ! Sum all alpha-active local dipoles for one real-emission subprocess.
    implicit none
    integer,intent(in) :: iint,ichan
    real(kind=8),intent(out) :: amp2_dip

    call compute_the_dipole_amps(iint,ichan)
    call square_the_dipole_amps(iint,ichan,amp2_dip)
  end subroutine evaluate_real_dipoles

  subroutine evaluate_recycled_born(iint,ichan,vol,multichannel_weight,born_value)
    ! Evaluate the tree Born contribution through the mapped points that were
    ! already evaluated for the CS dipoles.  The real phase-space Jacobian
    ! (including its flux) is kept outside the mapping; for initial-state
    ! maps the Bjorken and flux ratios cancel, leaving precisely 1/J_CS.
    use cs_dipole_mappings, only: cs_born_pushback_weight
    use scales
    use math_functions, only: factorial8
    use common
    implicit none
    integer,intent(in) :: iint,ichan
    real(kind=8),intent(in) :: vol,multichannel_weight
    real(kind=8),intent(out) :: born_value
    integer :: idip,info,ip,mapped_initial
    real(kind=8) :: pushback,born2,scale_born,as_born,coupling,luminosity,xborn_initial
    real(kind=8),dimension(pgl(ichan)%next) :: mass_real
    real(kind=8),parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real(kind=8),external :: alphaspdf

    born_value=0d0
    do ip=1,pgl(ichan)%next
       mass_real(ip)=phys_model%get_mass(pgl(ichan)%processes(ip,iint))
    enddo
    do idip=1,pgl(ichan)%dpl(iint)%ndip
       associate(dip => pgl(ichan)%dpl(iint)%dl(idip))
       if (.not.dipole_has_unit_recycling_measure(dip) .or. .not.dip%active) cycle
       call cs_born_pushback_weight(pgl(ichan)%ps(1)%p,dip%p_mapped,dip%dip_ijk,mass_real, &
            phys_model%get_mass(dip%dip_r_ijk_f(1)),pushback,info)
       if (info.ne.0) then
          write (*,*) 'error computing CS Born push-back weight',info
          stop 1
       endif
       if (dip%dip_ijk(1).le.2 .or. dip%dip_ijk(3).le.2) then
          if (dip%dip_ijk(1).le.2) then
             mapped_initial=dip%dip_r_ijk(1)
          else
             mapped_initial=dip%dip_r_ijk(2)
          endif
          xborn_initial=2d0*dip%p_mapped(0,mapped_initial)/sqrts
          if (xborn_initial.lt.1d0) pushback=pushback/(1d0-xborn_initial)
       endif
       born2=dipole_born_square(dip)
       call set_scale(scale_choice,size(dip%process_r),dip%p_mapped,dip%process_r,scale_born)
       if (use_lhapdf) then
          as_born=alphaspdf(scale_born)
       else
          as_born=alphas_Q(scale_born,2,alphas_MZ)
       endif
       coupling=1d0
       if (dip%amp%n_sing(1).lt.size(dip%process_r)-2) then
          coupling=coupling*(4d0*pi*as_born)**(size(dip%process_r)-2-dip%amp%n_sing(1))
       endif
       if (dip%amp%n_sing(1).ge.1) coupling=coupling*(2d0*4d0*pi*alphaEW)**dip%amp%n_sing(1)
       luminosity=reduced_process_luminosity(ichan,iint,dip,scale_born)
       born_value=born_value+pushback*born2*coupling*luminosity
       end associate
    enddo
    born_value=born_value*vol*pgl(ichan)%ps(1)%jac*conv*multichannel_weight
  end subroutine evaluate_recycled_born

  real(kind=8) function dipole_born_square(dip)
    type(dipole),intent(in) :: dip
    integer :: ih
    dipole_born_square=0d0
    if (use_real_gluons .and. all(dip%amp%n_qqbar(1:1).eq.0)) then
       do ih=1,dip%amp%n_amps
          dipole_born_square=dipole_born_square+dip%amp%amps_r(ih)*dip%col_fac*dip%amp%amps_r(ih)
       enddo
    else
       do ih=1,dip%amp%n_amps
          dipole_born_square=dipole_born_square+dble(dip%amp%amps(ih)*dip%col_fac*dconjg(dip%amp%amps(ih)))
       enddo
    endif
  end function dipole_born_square

  real(kind=8) function reduced_process_luminosity(ichan,iint,dip,scale)
    integer,intent(in) :: ichan,iint
    type(dipole),intent(in) :: dip
    real(kind=8),intent(in) :: scale
    integer :: ip,nlocal
    integer,dimension(size(dip%process_r)) :: process_r
    real(kind=8),dimension(2) :: xborn
    reduced_process_luminosity=0d0
    xborn(1)=2d0*dip%p_mapped(0,1)/sqrts
    xborn(2)=2d0*dip%p_mapped(0,2)/sqrts
    if (any(xborn.le.0d0) .or. any(xborn.ge.1d0+1d-12)) return
    xborn=min(1d0,xborn)
    do ip=1,pgl(ichan)%iden_iproc(iint)
       call reduced_process_for_dipole(dip,pgl(ichan)%iden_processes(:,ip,iint),process_r)
       if (.not.allocated(dip%recycling_history_weight)) then
          write (*,*) 'ERROR: recycling history weights were not initialised'
          stop 1
       endif
       nlocal=count_local_recycling_projections(ichan,iint,ip,dip,process_r)
       if (nlocal.eq.0) then
          write (*,*) 'ERROR: no local recycling projection for active dipole'
          stop 1
       endif
       reduced_process_luminosity=reduced_process_luminosity+dip%recycling_history_weight(ip)* &
            pgl(ichan)%idenCOandMAPfactor(ip,iint)* &
            pdf_for_pdg(process_r(1),xborn(1),scale)*pdf_for_pdg(process_r(2),xborn(2),scale)/ &
            (dble(process_normalisation(process_r))*dble(nlocal))
    enddo
  end function reduced_process_luminosity

  integer function count_local_recycling_projections(ichan,iint,iiden,dip,process_r)
    integer,intent(in) :: ichan,iint,iiden
    type(dipole),intent(in) :: dip
    integer,intent(in) :: process_r(:)
    integer :: idip
    integer,dimension(size(process_r)) :: process_r_other
    count_local_recycling_projections=0
    do idip=1,pgl(ichan)%dpl(iint)%ndip
       associate(dip_other => pgl(ichan)%dpl(iint)%dl(idip))
       if (.not.dipole_has_unit_recycling_measure(dip_other) .or. .not.dip_other%active) cycle
       call reduced_process_for_dipole(dip_other,pgl(ichan)%iden_processes(:,iiden,iint),process_r_other)
       if (same_recycled_born_key(dip,process_r,dip_other,process_r_other)) &
            count_local_recycling_projections=count_local_recycling_projections+1
       end associate
    enddo
  end function count_local_recycling_projections

  subroutine reduced_process_for_dipole(dip,process,process_r)
    type(dipole),intent(in) :: dip
    integer,intent(in) :: process(:)
    integer,intent(out) :: process_r(:)
    integer :: old,new,parent
    if (size(process_r).ne.size(process)-1) then
       write (*,*) 'invalid reduced-process storage'
       stop 1
    endif
    if (phys_model%is_gluon(process(dip%dip_ijk(2)))) then
       parent=process(dip%dip_ijk(1))
    elseif (phys_model%is_gluon(process(dip%dip_ijk(1)))) then
       parent=phys_model%get_antipart(process(dip%dip_ijk(2)))
    else
       ! The stored canonical parent is sufficient to fix whether the
       ! colour line is the ordinary or U1 gluon for flavour replicas.
       parent=dip%dip_r_ijk_f(1)
    endif
    new=0
    do old=1,size(process)
       if (old.eq.dip%dip_ijk(2)) cycle
       new=new+1
       process_r(new)=process(old)
       if (old.eq.dip%dip_ijk(1)) process_r(new)=parent
    enddo
  end subroutine reduced_process_for_dipole

  integer(kind=8) function process_normalisation(process)
    use math_functions, only: factorial8
    integer,intent(in) :: process(:)
    integer :: i,j,nkind
    integer,dimension(size(process),2) :: kinds
    process_normalisation=1_8
    do i=1,2
       if (process(i).eq.21) then
          process_normalisation=process_normalisation*16_8
       elseif (abs(process(i)).ge.1 .and. abs(process(i)).le.6) then
          process_normalisation=process_normalisation*6_8
       else
          process_normalisation=process_normalisation*2_8
       endif
    enddo
    nkind=0
    do i=3,size(process)
       do j=1,nkind
          if (kinds(j,1).eq.process(i)) then
             kinds(j,2)=kinds(j,2)+1
             exit
          endif
       enddo
       if (j.eq.nkind+1) then
          nkind=nkind+1
          kinds(nkind,1)=process(i)
          kinds(nkind,2)=1
       endif
    enddo
    do i=1,nkind
       process_normalisation=process_normalisation*factorial8(kinds(i,2))
    enddo
  end function process_normalisation

  real(kind=8) function pdf_for_pdg(pdg,x,scale)
    use common
    integer,intent(in) :: pdg
    real(kind=8),intent(in) :: x,scale
    real(kind=8),dimension(-6:7) :: pdf
    logical,dimension(-6:7) :: all_flavours
    integer :: index
    pdf_for_pdg=0d0
    if (.not.include_pdf) then
       pdf_for_pdg=1d0
       return
    endif
    if (x.le.0d0 .or. x.gt.1d0) return
    if (pdg.eq.21) then
       index=0
    elseif (pdg.eq.22) then
       index=7
    elseif (pdg.ge.-6 .and. pdg.le.6) then
       index=pdg
    else
       return
    endif
    if (use_lhapdf) then
       call evolvePDF(x,scale,pdf(-6))
       pdf_for_pdg=pdf(index)/x
    else
       ! The internal driver has the same flavour layout.  Ask for all
       ! flavours here because mapped initial states vary per dipole.
       all_flavours=.true.
       call PDF_eval(1,all_flavours,x,scale,pdf(-6))
       pdf_for_pdg=pdf(index)
    endif
  end function pdf_for_pdg

  subroutine square_the_dipole_amps(iint,ichan,amp2_dip,iunres,icol1,icol2,nselected,nmatched,nalpha_selected)
    use cs_lc_spin_dipoles
    use FeynmanRules
    implicit none
    integer,intent(in) :: iint,ichan
    integer,intent(in),optional :: iunres,icol1,icol2
    integer,intent(out),optional :: nselected
    integer,intent(out),optional :: nmatched
    integer,intent(out),optional :: nalpha_selected
    integer :: idip
    real(kind=8) :: amp2_dip,dip
    real(kind=8),parameter :: pi=3.14159265358979323846d0
    real(kind=8),dimension(pgl(ichan)%next) :: mass_real
    integer :: ij
    integer :: ipart
    logical :: use_collinear
    complex(kind=8),dimension(2,2) :: rho
    complex(kind=8),dimension(0:3,2) :: eps_parent
    use_collinear=present(icol1).or.present(icol2)
    if (use_collinear .and. .not.(present(icol1).and.present(icol2))) then
       write (*,*) 'collinear dipole selection needs both collinear legs'
       stop 1
    endif
    amp2_dip=0d0
    do ipart=1,pgl(ichan)%next
       mass_real(ipart)=phys_model%get_mass(pgl(ichan)%processes(ipart,iint))
    enddo
    if (present(nselected)) nselected=0
    if (present(nmatched)) nmatched=0
    if (present(nalpha_selected)) nalpha_selected=0
    do idip=1,pgl(ichan)%dpl(iint)%ndip
       if (use_collinear) then
          if (.not.collinear_dipole_matches(iint,ichan,idip,icol1,icol2)) cycle
       elseif (present(iunres)) then
          ! For a raw soft diagnostic, isolate the dipoles where the
          ! tested leg is the CS unresolved leg j.
          if (pgl(ichan)%dpl(iint)%dl(idip)%dip_ijk(2).ne.iunres) cycle
       endif
       if (present(nmatched)) nmatched=nmatched+1
       if (pgl(ichan)%dpl(iint)%dl(idip)%alpha_active) then
          if (present(nalpha_selected)) nalpha_selected=nalpha_selected+1
       endif
       if (.not.pgl(ichan)%dpl(iint)%dl(idip)%active) cycle
       if (present(nselected)) nselected=nselected+1
       ij=pgl(ichan)%dpl(iint)%dl(idip)%dip_r_ijk(1)
       call create_rho(iint,ichan,idip,rho)
       if (ij.gt.2) then
          call ext_gluon_cmplx(pgl(ichan)%dpl(iint)%dl(idip)%p_mapped_ij,-1, 1, eps_parent(0:3,1))
          call ext_gluon_cmplx(pgl(ichan)%dpl(iint)%dl(idip)%p_mapped_ij, 1, 1, eps_parent(0:3,2))
       else
          call ext_gluon_cmplx(-pgl(ichan)%dpl(iint)%dl(idip)%p_mapped_ij,-1, 1, eps_parent(0:3,1))
          call ext_gluon_cmplx(-pgl(ichan)%dpl(iint)%dl(idip)%p_mapped_ij, 1, 1, eps_parent(0:3,2))
       endif
       ! HELAS returns the external gluon wavefunction epsilon^*.  The CS
       ! helicity projection expects the physical polarization epsilon.
       eps_parent=conjg(eps_parent)
       call cs_lc_dipole_spinrho(pgl(ichan)%ps(1)%p,pgl(ichan)%processes(:,iint), &
            pgl(ichan)%dpl(iint)%dl(idip)%process_r,pgl(ichan)%dpl(iint)%dl(idip)%dip_ijk,1d0/(4d0*pi), &
            rho,eps_parent,dip,lc_weight=pgl(ichan)%dpl(iint)%dl(idip)%lc_weight, &
            mass_real=mass_real, &
            mass_parent=phys_model%get_mass(pgl(ichan)%dpl(iint)%dl(idip)%dip_r_ijk_f(1)))
       amp2_dip=amp2_dip+dip
    enddo
  end subroutine square_the_dipole_amps

  logical function collinear_dipole_matches(iint,ichan,idip,icol1,icol2)
    implicit none
    integer,intent(in) :: iint,ichan,idip,icol1,icol2
    integer :: dip_i,dip_j,iini,ifin
    dip_i=pgl(ichan)%dpl(iint)%dl(idip)%dip_ijk(1)
    dip_j=pgl(ichan)%dpl(iint)%dl(idip)%dip_ijk(2)
    collinear_dipole_matches=.false.
    if (icol1.le.2 .and. icol2.le.2) return
    if (icol1.le.2 .or. icol2.le.2) then
       if (icol1.le.2) then
          iini=icol1
          ifin=icol2
       else
          iini=icol2
          ifin=icol1
       endif
       ! For an initial-final collinear pair, the initial leg must be the
       ! emitter.  A dipole with the initial leg as spectator belongs to a
       ! different singular configuration and must not be selected here.
       collinear_dipole_matches=(dip_i.eq.iini .and. dip_j.eq.ifin)
    else
       collinear_dipole_matches=((dip_i.eq.icol1 .and. dip_j.eq.icol2) .or. &
            (dip_i.eq.icol2 .and. dip_j.eq.icol1))
    endif
  end function collinear_dipole_matches

  logical function has_collinear_dipoles(iint,ichan,icol1,icol2)
    implicit none
    integer,intent(in) :: iint,ichan,icol1,icol2
    integer :: idip
    has_collinear_dipoles=.false.
    do idip=1,pgl(ichan)%dpl(iint)%ndip
       if (collinear_dipole_matches(iint,ichan,idip,icol1,icol2)) then
          has_collinear_dipoles=.true.
          return
       endif
    enddo
  end function has_collinear_dipoles

  subroutine create_rho(iint,ichan,idip,rho)
    implicit none
    integer,intent(in) :: iint,ichan,idip
    complex(kind=8),dimension(2,2),intent(out) :: rho
    integer :: ilook,ih1,ih2,ij,a,b
    complex(kind=8) :: term
    real(kind=8) :: hermitian_scale
    associate(dip => pgl(ichan)%dpl(iint)%dl(idip))
    rho=(0d0,0d0)
    ij=dip%dip_r_ijk(1)
    do ilook=1,size(dip%rho_lookup_ih1)
       ih1=dip%rho_lookup_ih1(ilook)
       ih2=dip%rho_lookup_ih2(ilook)
       a=(dip%amp%spins(ij,1,ih1)+3)/2
       b=(dip%amp%spins(ij,1,ih2)+3)/2
       term=dip%amp%amps(ih1)*dip%col_fac*dconjg(dip%amp%amps(ih2))
       rho(a,b)=rho(a,b)+term
    enddo
    if (dip%rho_lookup_upper) rho(2,1)=dconjg(rho(1,2))
    if (.not.dip%rho_hermitian_checked .and. phys_model%is_gluon(dip%dip_r_ijk_f(1))) then
       hermitian_scale=max(1d0,maxval(abs(rho)))
       dip%rho_hermitian_checked=(abs(rho(1,2)-dconjg(rho(2,1))).le.1d-10*hermitian_scale .and. &
            abs(aimag(rho(1,1))).le.1d-10*hermitian_scale .and. &
            abs(aimag(rho(2,2))).le.1d-10*hermitian_scale)
       if (dip%rho_hermitian_checked) then
          call compact_rho_lookup(dip)
          dip%rho_lookup_upper=.true.
       endif
    endif
    end associate
  end subroutine create_rho

  subroutine compact_rho_lookup(dip)
    implicit none
    type(dipole),intent(inout) :: dip
    integer,allocatable :: keep_ih1(:),keep_ih2(:)
    integer :: i,nkeep,ij,a,b
    ij=dip%dip_r_ijk(1)
    nkeep=0
    do i=1,size(dip%rho_lookup_ih1)
       a=(dip%amp%spins(ij,1,dip%rho_lookup_ih1(i))+3)/2
       b=(dip%amp%spins(ij,1,dip%rho_lookup_ih2(i))+3)/2
       if (a.le.b) nkeep=nkeep+1
    enddo
    allocate(keep_ih1(nkeep),keep_ih2(nkeep))
    nkeep=0
    do i=1,size(dip%rho_lookup_ih1)
       a=(dip%amp%spins(ij,1,dip%rho_lookup_ih1(i))+3)/2
       b=(dip%amp%spins(ij,1,dip%rho_lookup_ih2(i))+3)/2
       if (a.gt.b) cycle
       nkeep=nkeep+1
       keep_ih1(nkeep)=dip%rho_lookup_ih1(i)
       keep_ih2(nkeep)=dip%rho_lookup_ih2(i)
    enddo
    call move_alloc(keep_ih1,dip%rho_lookup_ih1)
    call move_alloc(keep_ih2,dip%rho_lookup_ih2)
  end subroutine compact_rho_lookup

  subroutine compute_the_amps(iint,ichan,use_amplitude_library)
    use amp_lib
    implicit none
    integer,intent(in) :: iint,ichan
    logical,intent(in) :: use_amplitude_library
    if (.not. use_amplitude_library) then
       call pgl(ichan)%amps(iint)%evaluate(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
            pgl(ichan)%hel,read_proc_from_file,phys_model)
    else
       call evaluate_amp(ichan,iint,pgl(ichan)%ps(1)%p,pgl(ichan)%amps(iint)%amps)
    endif
  end subroutine compute_the_amps

  subroutine square_the_amps(iint,ichan)
    implicit none
    integer,intent(in) :: iint,ichan
    integer :: iproc,ih
    iproc=0
    pgl(ichan)%amp2=0d0
    if (keep_processes_separate) then
       if (use_real_gluons .and. all(pgl(ichan)%amps(iint)%n_qqbar(1:1).eq.0)) then
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih)=pgl(ichan)%amps(iint)%amps_r(ih)*&
                  pgl(ichan)%col_fac(iint)*pgl(ichan)%amps(iint)%amps_r(ih)*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
          enddo
       else
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih)=dble(pgl(ichan)%amps(iint)%amps(ih)*&
                  pgl(ichan)%col_fac(iint)*dconjg(pgl(ichan)%amps(iint)%amps(ih)))*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
          enddo
       endif
    else
       if (use_real_gluons .and. all(pgl(ichan)%amps(iint)%n_qqbar(1:1).eq.0)) then
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih)=pgl(ichan)%amps(iint)%amps_r(ih)*&
                  pgl(ichan)%col_fac(iproc)*pgl(ichan)%amps(iint)%amps_r(ih)*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
          enddo
       else
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih)=dble(pgl(ichan)%amps(iint)%amps(ih)*&
                  pgl(ichan)%col_fac(iproc)*dconjg(pgl(ichan)%amps(iint)%amps(ih)))*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
          enddo
       endif
    endif
  end subroutine square_the_amps
end module subtraction
