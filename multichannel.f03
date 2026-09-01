module multichannel
  use handling_processes
  use simple_integrator_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!  use mint_module
contains
  subroutine compute_multichannel_weight(ichan,iint,ps,weight,adaptation_class)
    ! Computes the multichannel weight 'weight' when there are
    ! 'chans(0)' channels (that are listed in the array 'chans(1:)') and
    ! the current channel is 'ichan'. The momenta 'p' have been
    ! generate with phase-space jacobian 'jac' using the random
    ! variables 'x' within the channel 'ichan'.
    !
    ! weight = 1/( J_{ichan}*[ sum_{i=1}^{chans(0)} 1/J_i ] )
    !
    ! with J_i the combined Jacobian coming from MINT and the
    ! phase-space.
    implicit none
    integer,intent(in) :: ichan,iint
    integer,intent(in),optional :: adaptation_class
    type(psv),intent(in) :: ps
    type(psv) :: ps_local
    real(kind=8),dimension(:,:),allocatable :: p_target
    real(kind=8),dimension(:),intent(out) :: weight
    integer :: i,iproc,k,a,iproc_first,iproc_last,self_count,iadapt
    integer :: label
    real(kind=8) :: vol_ichan,vol,log_denominator,log_ratio
    weight=0d0
    iadapt=1
    if (present(adaptation_class)) iadapt=adaptation_class
    if (iadapt.lt.1) then
       write (*,*) 'Invalid adaptation class in multichannel weight',iadapt
       stop 1
    endif
    if (.not.allocated(pgl)) then
       write (*,*) 'Cannot compute a multichannel weight before process initialisation'
       stop 1
    endif
    if (ichan.lt.1 .or. ichan.gt.ngroups) then
       write (*,*) 'Invalid current channel in multichannel weight',ichan,ngroups
       stop 1
    endif
    if (pgl(ichan)%next.lt.4 .or. pgl(ichan)%nproc.lt.1) then
       write (*,*) 'Invalid current process dimensions in multichannel weight',&
            pgl(ichan)%next,pgl(ichan)%nproc
       stop 1
    endif
    if (size(weight).ne.pgl(ichan)%nproc) then
       write (*,*) 'Invalid output extent in multichannel weight',size(weight),pgl(ichan)%nproc
       stop 1
    endif
    if (.not.allocated(pgl(ichan)%phase_space_permutations)) then
       write (*,*) 'Missing phase-space permutations in multichannel weight'
       stop 1
    endif
    if (size(pgl(ichan)%phase_space_permutations,1).lt.pgl(ichan)%next .or. &
         size(pgl(ichan)%phase_space_permutations,2).lt.pgl(ichan)%nproc) then
       write (*,*) 'Incomplete phase-space permutations in multichannel weight'
       stop 1
    endif
    if (.not.allocated(pgl(ichan)%multichan%number_of_channels)) then
       write (*,*) 'Missing channel counts in multichannel weight'
       stop 1
    endif
    if (size(pgl(ichan)%multichan%number_of_channels).lt.pgl(ichan)%nproc) then
       write (*,*) 'Incomplete channel counts in multichannel weight'
       stop 1
    endif
    if (.not.allocated(pgl(ichan)%multichan%channels)) then
       write (*,*) 'Missing channel list in multichannel weight'
       stop 1
    endif
    if (.not.allocated(pgl(ichan)%multichan%channel_permutations)) then
       write (*,*) 'Missing channel permutations in multichannel weight'
       stop 1
    endif
    if (size(pgl(ichan)%multichan%channels,2).lt.pgl(ichan)%nproc .or. &
         size(pgl(ichan)%multichan%channel_permutations,1).lt.pgl(ichan)%next .or. &
         size(pgl(ichan)%multichan%channel_permutations,3).lt.pgl(ichan)%nproc) then
       write (*,*) 'Incomplete channel metadata in multichannel weight'
       stop 1
    endif
    if (keep_processes_separate) then
       if (iint.lt.1 .or. iint.gt.pgl(ichan)%nproc) then
          write (*,*) 'Invalid process in multichannel weight',iint,pgl(ichan)%nproc
          stop 1
       endif
       iproc_first=iint
       iproc_last=iint
    else
       iproc_first=1
       iproc_last=pgl(ichan)%nproc
       do iproc=2,pgl(ichan)%nproc
          if (any(pgl(ichan)%phase_space_permutations(:,iproc).ne.&
               pgl(ichan)%phase_space_permutations(:,1))) then
             write (*,*) 'keep_processes_separate=false is incompatible with different phase-space permutations'
             stop 1
          endif
       enddo
    endif
    do iproc=iproc_first,iproc_last
       if (pgl(ichan)%multichan%number_of_channels(iproc).lt.1 .or. &
            pgl(ichan)%multichan%number_of_channels(iproc).gt.&
            size(pgl(ichan)%multichan%channels,1) .or. &
            pgl(ichan)%multichan%number_of_channels(iproc).gt.&
            size(pgl(ichan)%multichan%channel_permutations,2)) then
          write (*,*) 'Invalid channel count in multichannel weight',iproc,&
               pgl(ichan)%multichan%number_of_channels(iproc)
          stop 1
       endif
       do label=1,pgl(ichan)%next
          if (count(pgl(ichan)%phase_space_permutations(1:pgl(ichan)%next,iproc).eq.label).ne.1) then
             write (*,*) 'Invalid current phase-space permutation in multichannel weight',iproc
             stop 1
          endif
       enddo
       do k=1,pgl(ichan)%multichan%number_of_channels(iproc)
          i=pgl(ichan)%multichan%channels(k,iproc)
          if (i.lt.1 .or. i.gt.ngroups) then
             write (*,*) 'Out-of-range partner channel in multichannel weight',i,ngroups
             stop 1
          endif
          if (pgl(i)%next.ne.pgl(ichan)%next .or. &
               pgl(i)%ndim.ne.pgl(ichan)%ndim) then
             write (*,*) 'Incompatible partner phase-space dimensions in multichannel weight',ichan,i
             stop 1
          endif
          if (.not.allocated(pgl(i)%phase_space)) then
             write (*,*) 'Missing partner phase-space generator in multichannel weight',i
             stop 1
          endif
          do label=1,pgl(ichan)%next
             if (count(pgl(ichan)%multichan%channel_permutations(&
                  1:pgl(ichan)%next,k,iproc).eq.label).ne.1) then
                write (*,*) 'Invalid partner permutation in multichannel weight',iproc,k
                stop 1
             endif
          enddo
       enddo
    enddo
    if (.not. use_colour_singlet_multichannel) then
       do iproc=iproc_first,iproc_last
          weight(iproc)=1d0/dble(pgl(ichan)%multichan%number_of_channels(iproc))
       enddo
       return
    endif
    if (.not.pgl(ichan)%phase_space%can_invert_momenta) then
       ! HAAG and the pT-based generator do not implement the inverse map.
       ! Use the same uniform partition in every partner channel rather than
       ! attempting a point-dependent multichannel weight.
       do iproc=iproc_first,iproc_last
          weight(iproc)=1d0/dble(pgl(ichan)%multichan%number_of_channels(iproc))
       enddo
       return
    endif
    if (.not.allocated(ps%p) .or. .not.allocated(ps%x)) then
       weight=huge(1d0)
       return
    endif
    if (lbound(ps%p,1).ne.0 .or. ubound(ps%p,1).lt.3 .or. &
         size(ps%p,2).lt.pgl(ichan)%next .or. size(ps%x).lt.pgl(ichan)%ndim) then
       weight=huge(1d0)
       return
    endif
    if (.not.all(ieee_is_finite(ps%p(0:3,1:pgl(ichan)%next))) .or. &
         .not.all(ieee_is_finite(ps%x)) .or. .not.ieee_is_finite(ps%jac)) then
       weight=huge(1d0)
       return
    endif
    if (ps%jac.le.0d0) then
       weight=huge(1d0)
       return
    endif
    allocate(p_target(0:3,pgl(ichan)%next))
    p_target=ps%p(0:3,1:pgl(ichan)%next)
    call simple_integrator%compute_wgt_from_x(ichan,ps%x(1:pgl(ichan)%ndim),vol_ichan,&
         adaptation_class=iadapt)
    if (.not.ieee_is_finite(vol_ichan)) then
       weight=huge(1d0)
       return
    endif
    if (vol_ichan.le.0d0) then
       weight=huge(1d0)
       return
    endif
    do iproc=iproc_first,iproc_last
       ! The current density contributes log(1)=0.  Accumulate all partner
       ! density ratios with log-sum-exp so products and ratios cannot overflow
       ! before the final reciprocal is formed.
       log_denominator=0d0
       self_count=0
       do k=1,pgl(ichan)%multichan%number_of_channels(iproc)
          if (pgl(ichan)%multichan%channels(k,iproc).eq.ichan .and. all(&
               pgl(ichan)%multichan%channel_permutations(:,k,iproc).eq.&
               pgl(ichan)%phase_space_permutations(:,iproc))) &
               self_count=self_count+1
       enddo
       if (self_count.ne.1) then
          write (*,*) 'Could not identify exactly one current multichannel density',ichan,iproc,self_count
          stop 1
       endif
       do k=1,pgl(ichan)%multichan%number_of_channels(iproc)
          i=pgl(ichan)%multichan%channels(k,iproc)
          if (i.eq.ichan .and. all(&
               pgl(ichan)%multichan%channel_permutations(:,k,iproc).eq.&
               pgl(ichan)%phase_space_permutations(:,iproc))) then
             ! The current channel is already represented by log(1)=0 above.
             cycle
          endif
          ps_local=ps
          do a=1,pgl(ichan)%next
             ps_local%p(:,a)=p_target(:,&
                  pgl(ichan)%multichan%channel_permutations(a,k,iproc))
          enddo
          call pgl(i)%phase_space%compute_x_from_momenta(ps_local)
          ! A point outside a partner map's support has zero density in that
          ! map.  It must be omitted from the MIS denominator; replacing the
          ! whole point by a local 1/N weight would not form a partition of
          ! unity across the partner channels.
          if (.not.ieee_is_finite(ps_local%jac)) cycle
          if (ps_local%jac.le.0d0) cycle
          if (.not.allocated(ps_local%x)) cycle
          if (size(ps_local%x).lt.pgl(i)%ndim) cycle
          if (any(.not.ieee_is_finite(ps_local%x(1:pgl(i)%ndim)))) cycle
          call simple_integrator%compute_wgt_from_x(i,ps_local%x(1:pgl(i)%ndim),vol,&
               adaptation_class=iadapt)
          if (.not.ieee_is_finite(vol)) cycle
          if (vol.le.0d0) cycle
          log_ratio=log(ps%jac)+log(vol_ichan)-log(ps_local%jac)-log(vol)
          if (.not.ieee_is_finite(log_ratio)) cycle
          call add_log_density(log_denominator,log_ratio)
       enddo
       weight(iproc)=exp(-log_denominator)
       if (.not.ieee_is_finite(weight(iproc))) then
          weight=huge(1d0)
          return
       endif
       if (weight(iproc).lt.0d0 .or. weight(iproc).gt.1d0) then
          weight=huge(1d0)
          return
       endif
    enddo
  contains
    subroutine add_log_density(log_total,log_value)
      real(kind=8),intent(inout) :: log_total
      real(kind=8),intent(in) :: log_value
      if (log_value.gt.log_total) then
         log_total=log_value+log(1d0+exp(log_total-log_value))
      else
         log_total=log_total+log(1d0+exp(log_value-log_total))
      endif
    end subroutine add_log_density
  end subroutine compute_multichannel_weight

  subroutine setup_optimised_multichannel_weight_computation(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer,dimension(:),allocatable :: all_unique_chans,all_unique_chans_inv,map_proc_to_group
    integer,dimension(:,:),allocatable :: all_unique_channelgroups
    integer :: iproc,ichan,nchans,nchans_group,i,ios,channel_capacity
    integer(kind=8) :: channel_capacity8
    logical :: found
    if (pgl%nproc.lt.1 .or. pgl%multichan%max_channels.lt.1 .or. ngroups.lt.1) then
       write (*,*) 'Invalid dimensions while optimising multichannel metadata',&
            pgl%nproc,pgl%multichan%max_channels,ngroups
       stop 1
    endif
    if (.not.allocated(pgl%multichan%number_of_channels)) then
       write (*,*) 'Missing channel counts while optimising multichannel metadata'
       stop 1
    endif
    if (.not.allocated(pgl%multichan%channels)) then
       write (*,*) 'Missing channel list while optimising multichannel metadata'
       stop 1
    endif
    if (size(pgl%multichan%number_of_channels).lt.pgl%nproc .or. &
         size(pgl%multichan%channels,1).lt.pgl%multichan%max_channels .or. &
         size(pgl%multichan%channels,2).lt.pgl%nproc) then
       write (*,*) 'Incomplete channel data while optimising multichannel metadata'
       stop 1
    endif
    channel_capacity8=int(pgl%nproc,kind=8)*int(pgl%multichan%max_channels,kind=8)
    if (channel_capacity8.gt.200000000_8 .or. channel_capacity8.gt.huge(channel_capacity)) then
       write (*,*) 'Multichannel metadata is too large to optimise',channel_capacity8
       stop 1
    endif
    channel_capacity=int(channel_capacity8)
    allocate(all_unique_chans(channel_capacity),all_unique_chans_inv(ngroups),&
         all_unique_channelgroups(0:pgl%multichan%max_channels,pgl%nproc),&
         map_proc_to_group(pgl%nproc),stat=ios)
    if (ios.ne.0) then
       write (*,*) 'Could not allocate optimised multichannel metadata',ios
       stop 1
    endif
    all_unique_chans=0
    all_unique_chans_inv=0
    all_unique_channelgroups=0
    map_proc_to_group=0
    nchans=0
    nchans_group=0
    do iproc=1,pgl%nproc
       if (pgl%multichan%number_of_channels(iproc).lt.1 .or. &
            pgl%multichan%number_of_channels(iproc).gt.pgl%multichan%max_channels) then
          write (*,*) 'Invalid process channel count while optimising metadata',iproc,&
               pgl%multichan%number_of_channels(iproc)
          stop 1
       endif
       do ichan=1,pgl%multichan%number_of_channels(iproc)
          if (pgl%multichan%channels(ichan,iproc).lt.1 .or. &
               pgl%multichan%channels(ichan,iproc).gt.ngroups) then
             write (*,*) 'Out-of-range channel while optimising metadata',&
                  pgl%multichan%channels(ichan,iproc),ngroups
             stop 1
          endif
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
    if (allocated(pgl%multichan%unique_channel_list)) deallocate(pgl%multichan%unique_channel_list)
    if (allocated(pgl%multichan%unique_channelgroup_list)) &
         deallocate(pgl%multichan%unique_channelgroup_list)
    if (allocated(pgl%multichan%map_proc_to_channelgroup)) &
         deallocate(pgl%multichan%map_proc_to_channelgroup)
    allocate(pgl%multichan%unique_channel_list(nchans),stat=ios)
    if (ios.ne.0) then
       write (*,*) 'Could not allocate unique multichannel list',ios
       stop 1
    endif
    pgl%multichan%unique_channel_list=all_unique_chans(1:nchans)
    allocate(pgl%multichan%unique_channelgroup_list(0:pgl%multichan%max_channels,nchans_group),stat=ios)
    if (ios.ne.0) then
       write (*,*) 'Could not allocate unique multichannel groups',ios
       stop 1
    endif
    pgl%multichan%unique_channelgroup_list(0:pgl%multichan%max_channels,1:nchans_group)= &
         all_unique_channelgroups(0:pgl%multichan%max_channels,1:nchans_group)
    allocate(pgl%multichan%map_proc_to_channelgroup(pgl%nproc),stat=ios)
    if (ios.ne.0) then
       write (*,*) 'Could not allocate multichannel process map',ios
       stop 1
    endif
    pgl%multichan%map_proc_to_channelgroup(1:pgl%nproc)=map_proc_to_group(1:pgl%nproc)
    ! Keep the uncompressed channel list: phase-space permutations make a
    ! density specific to an individual matrix-element row.
  end subroutine setup_optimised_multichannel_weight_computation

  subroutine determine_multi_channel_size(part,n_ps)
    use math_functions
    implicit none
    integer,dimension(:),intent(in) :: part
    integer,intent(out) :: n_ps
    integer :: i
    if (.not. use_colour_singlet_multichannel) then
       n_ps=1
       return
    endif
    n_ps=0
    do i=1,size(part)
       if (phys_model%is_singlet(part(i))) n_ps=n_ps+1
    enddo
    n_ps=factorial(n_ps)
  end subroutine determine_multi_channel_size

end module multichannel
