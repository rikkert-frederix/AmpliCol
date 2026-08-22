module multichannel
  use handling_processes
  use simple_integrator_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!  use mint_module
contains
  subroutine compute_multichannel_weight(ichan,iint,ps,weight)
    ! Computes the multichannel weight 'weight' when there are
    ! 'chans(0)' channels (that are listed in the array 'chans(1:)') and
    ! the current channel is 'ichan'. The momenta 'p' have been
    ! generate with phase-space jacobian 'jac' using the random
    ! variables 'x' within the channel 'ichan'.
    !
    ! weight = 1/( Jref_{ichan}*[ sum_{i=1}^{chans(0)} 1/Jref_i ] )
    !
    ! with Jref_i the phase-space Jacobian times the fixed reference-grid
    ! Jacobian.  The live grid Jacobian remains in the generated integrand.
    implicit none
    integer,intent(in) :: ichan,iint
    type(psv),intent(in) :: ps
    type(psv) :: ps_local
    real(kind=8),dimension(0:3,pgl(ichan)%next) :: p_target
    real(kind=8),dimension(pgl(ichan)%nproc),intent(out) :: weight
    integer :: i,iproc,k,a,iproc_first,iproc_last,self_count
    real(kind=8) :: vol_ichan,vol,denominator,density_ratio
    weight=0d0
    if (keep_processes_separate) then
       iproc_first=iint
       iproc_last=iint
    else
       iproc_first=1
       iproc_last=pgl(ichan)%nproc
       do iproc=2,pgl(ichan)%nproc
          if (any(pgl(ichan)%phase_space_permutations(:,iproc).ne.&
               pgl(ichan)%phase_space_permutations(:,1))) then
             write (*,*) '--combine_subprocesses is incompatible with different phase-space permutations'
             stop 1
          endif
       enddo
    endif
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
    p_target=ps%p
    ! The generated point retains its live-grid Jacobian in the integrand, but
    ! the multichannel partition is built exclusively from the fixed reference
    ! grids.  During warm-up the reference follows the live grid.
    call simple_integrator%compute_reference_wgt_from_x(ichan,ps%x,vol_ichan)
    if (ps%jac.le.0d0 .or. .not.ieee_is_finite(ps%jac) .or. &
         vol_ichan.le.0d0 .or. .not.ieee_is_finite(vol_ichan)) then
       write (*,*) 'Unexpected failure of the current multichannel map',&
            ichan,ps%jac,vol_ichan
       write (99,*) 'Unexpected failure of the current multichannel map',&
            ichan,ps%jac,vol_ichan
       stop 1
    endif
    do iproc=iproc_first,iproc_last
       denominator=0d0
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
             denominator=denominator+1d0
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
          if (ps_local%jac.le.0d0 .or. .not.ieee_is_finite(ps_local%jac)) cycle
          if (any(.not.ieee_is_finite(ps_local%x(1:pgl(i)%ndim)))) cycle
          if (any(ps_local%x(1:pgl(i)%ndim).le.0d0) .or. &
               any(ps_local%x(1:pgl(i)%ndim).ge.1d0)) cycle
          call simple_integrator%compute_reference_wgt_from_x(i,ps_local%x,vol)
          if (vol.le.0d0 .or. .not.ieee_is_finite(vol)) cycle
          density_ratio=ps%jac*vol_ichan/(ps_local%jac*vol)
          if (density_ratio.le.0d0 .or. &
               .not.ieee_is_finite(density_ratio)) cycle
          denominator=denominator+density_ratio
       enddo
       if (denominator.gt.0d0 .and. ieee_is_finite(denominator)) then
          weight(iproc)=1d0/denominator
       else
          write (*,*) 'Unexpected invalid multichannel denominator',&
               ichan,iproc,denominator
          write (99,*) 'Unexpected invalid multichannel denominator',&
               ichan,iproc,denominator
          stop 1
       endif
    enddo
  end subroutine compute_multichannel_weight

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
