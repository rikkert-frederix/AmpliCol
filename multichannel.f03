module multichannel
  use handling_processes
  use simple_integrator_mod
  use phase_space_gen23_mod, only: phase_space_gen23,&
       gen23_momentum_cache,build_gen23_momentum_cache
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!  use mint_module
contains
  subroutine compute_family_multichannel_weight(ifamily,ps,weight)
    ! Evaluate the fixed-prior hierarchical mixture owned by one family.
    ! Every submap uses the family's frozen reference grid in the partition,
    ! while the selected forward contribution retains the shared live-grid
    ! Jacobian supplied by the integrator.
    implicit none
    integer,intent(in) :: ifamily
    type(psv),intent(in) :: ps
    real(kind=8),intent(out) :: weight
    type(psv) :: ps_local
    type(gen23_momentum_cache),dimension(:),allocatable :: caches
    logical,dimension(:),allocatable :: cache_ready
    real(kind=8),dimension(0:3,pgl(ifamily)%next) :: p_target
    real(kind=8),dimension(pgl(ifamily)%nsubmaps) :: densities
    real(kind=8) :: reference_wgt,total_density
    integer :: imap,a,selected,permutation_id

    weight=0d0
    densities=0d0
    selected=pgl(ifamily)%selected_map
    if (selected.lt.1 .or. selected.gt.pgl(ifamily)%nsubmaps) then
       write (*,*) 'No valid selected phase map in integration family',&
            ifamily,selected
       stop 1
    endif
    if (.not.use_colour_singlet_multichannel) then
       weight=1d0
       return
    endif
    if (pgl(ifamily)%nsubmaps.eq.1) then
       weight=1d0
       return
    endif
    p_target=ps%p
    allocate(caches(nphase_permutations))
    allocate(cache_ready(nphase_permutations))
    cache_ready=.false.
    do imap=1,pgl(ifamily)%nsubmaps
       if (.not.pgl(ifamily)%phase_maps(imap)%phase_space%can_invert_momenta) then
          weight=1d0
          return
       endif
       ps_local=ps
       do a=1,pgl(ifamily)%next
          ps_local%p(:,a)=p_target(:,&
               pgl(ifamily)%phase_maps(imap)%permutation(a))
       enddo
       if (imap.eq.selected) then
          ps_local%jac=ps%jac
          ps_local%x(1:pgl(ifamily)%ndim)=ps%x(1:pgl(ifamily)%ndim)
       else
          permutation_id=pgl(ifamily)%phase_maps(imap)%permutation_id
          select type (phase_map=>&
               pgl(ifamily)%phase_maps(imap)%phase_space)
          type is (phase_space_gen23)
             if (.not.cache_ready(permutation_id)) then
                call build_gen23_momentum_cache(ps_local%p,&
                     caches(permutation_id))
                cache_ready(permutation_id)=.true.
             endif
             call phase_map%compute_x_from_cache(ps_local,&
                  caches(permutation_id))
          class default
             write (*,*) 'A multi-map integration family requires gen23',&
                  ifamily,imap
             stop 1
          end select
       endif
       ! Partner inversions outside their support contribute zero silently.
       if (ps_local%jac.le.0d0 .or. .not.ieee_is_finite(ps_local%jac)) cycle
       if (any(.not.ieee_is_finite(&
            ps_local%x(1:pgl(ifamily)%ndim)))) cycle
       if (any(ps_local%x(1:pgl(ifamily)%ndim).le.0d0) .or.&
            any(ps_local%x(1:pgl(ifamily)%ndim).ge.1d0)) cycle
       call simple_integrator%compute_reference_wgt_from_x(&
            ifamily,ps_local%x,reference_wgt)
       if (reference_wgt.le.0d0 .or.&
            .not.ieee_is_finite(reference_wgt)) cycle
       densities(imap)=1d0/(ps_local%jac*reference_wgt)
       if (.not.ieee_is_finite(densities(imap)) .or.&
            densities(imap).le.0d0) densities(imap)=0d0
    enddo
    if (densities(selected).le.0d0) then
       write (*,*) 'Unexpected failure of selected family phase map',&
            ifamily,selected,ps%jac
       write (99,*) 'Unexpected failure of selected family phase map',&
            ifamily,selected,ps%jac
       stop 1
    endif
    total_density=sum(densities)
    if (total_density.le.0d0 .or. .not.ieee_is_finite(total_density)) then
       write (*,*) 'Invalid family phase-map density sum',ifamily,total_density
       stop 1
    endif
    weight=dble(pgl(ifamily)%nsubmaps)*densities(selected)/total_density
    deallocate(caches,cache_ready)
  end subroutine compute_family_multichannel_weight

  subroutine setup_optimised_multichannel_weight_computation(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer,dimension(pgl%nproc*pgl%multichan%max_channels) :: all_unique_chans
    integer,dimension(pgl%nproc*ngroups) :: all_unique_chans_inv
    integer,dimension(0:pgl%multichan%max_channels,pgl%nproc) :: all_unique_channelgroups
    integer,dimension(pgl%nproc) :: map_proc_to_group
    integer :: iproc,ichan,nchans,nchans_group,i
    logical :: found
    nchans=0
    nchans_group=0
    do iproc=1,pgl%nproc
       do ichan=1,pgl%multichan%number_of_channels(iproc)
          if (all(all_unique_chans(1:nchans).ne.pgl%multichan%channels(ichan,iproc))) then
             nchans=nchans+1
             all_unique_chans(nchans)=pgl%multichan%channels(ichan,iproc)
             all_unique_chans_inv(pgl%multichan%channels(ichan,iproc))=nchans
          endif
       enddo
    enddo
    do iproc=1,pgl%nproc
       found=.false.
       do i=1,nchans_group
          if (all_unique_channelgroups(0,i).ne.pgl%multichan%number_of_channels(iproc))cycle
          ! Stored channel groups use compact positions in all_unique_chans,
          ! whereas pgl%multichan%channels still contains global group IDs.
          if (all(all_unique_channelgroups(1:pgl%multichan%number_of_channels(iproc),i).eq. &
               all_unique_chans_inv(pgl%multichan%channels(1:pgl%multichan%number_of_channels(iproc),iproc)))) then
             map_proc_to_group(iproc)=i
             found=.true.
             exit
          endif
       enddo
       if (.not.found) then
          nchans_group=nchans_group+1
          all_unique_channelgroups(0,nchans_group)= &
               pgl%multichan%number_of_channels(iproc)
          all_unique_channelgroups(1:pgl%multichan%number_of_channels(iproc),nchans_group)= &
               all_unique_chans_inv(pgl%multichan%channels(1:pgl%multichan%number_of_channels(iproc),iproc))
          map_proc_to_group(iproc)=nchans_group
       endif
    enddo
    pgl%multichan%n_unique_channels=nchans
    pgl%multichan%n_unique_channelgroups=nchans_group
    allocate(pgl%multichan%unique_channel_list(nchans))
    pgl%multichan%unique_channel_list=all_unique_chans(1:nchans)
    allocate(pgl%multichan%unique_channelgroup_list(0:pgl%multichan%max_channels,nchans_group))
    pgl%multichan%unique_channelgroup_list(0:pgl%multichan%max_channels,1:nchans_group)= &
         all_unique_channelgroups(0:pgl%multichan%max_channels,1:nchans_group)
    allocate(pgl%multichan%map_proc_to_channelgroup(pgl%nproc))
    pgl%multichan%map_proc_to_channelgroup(1:pgl%nproc)=map_proc_to_group(1:pgl%nproc)
    ! Keep the uncompressed channel list: phase-space permutations make a
    ! density specific to an individual matrix-element row.
  end subroutine setup_optimised_multichannel_weight_computation

  subroutine determine_multi_channel_size(part,n_ps)
    use math_functions
    implicit none
    integer,dimension(1:next),intent(in) :: part
    integer,intent(out) :: n_ps
    integer :: i
    if (.not. use_colour_singlet_multichannel) then
       n_ps=1
       return
    endif
    n_ps=0
    do i=1,next
       if (phys_model%is_singlet(part(i))) n_ps=n_ps+1
    enddo
    n_ps=factorial(n_ps)
  end subroutine determine_multi_channel_size

end module multichannel
