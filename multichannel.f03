module multichannel
  use handling_processes
  use mint_module
contains
  subroutine compute_multichannel_weight(ichan,x,p,jac,weight)
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
    integer,intent(in) :: ichan
    real(kind=8),dimension(0:3,next),intent(in) :: p
    real(kind=8),dimension(ndim),intent(in) :: x
    real(kind=8),intent(in) :: jac
    real(kind=8),dimension(pgl(ichan)%multichan%n_unique_channels) :: factors
    real(kind=8),dimension(pgl(ichan)%multichan%n_unique_channelgroups) :: weight_factors
    real(kind=8),dimension(pgl(ichan)%nproc),intent(out) :: weight
    integer :: i,j,iproc,ii
    real(kind=8) :: vol_ichan,vol
    if (.not. use_colour_singlet_multichannel) then
       weight(1:pgl(ichan)%nproc)=1d0/dble(pgl(ichan)%multichan%number_of_channels(1:pgl(ichan)%nproc))
       return
    endif
    call mint_get_jacobian_from_x(ichan,x,vol_ichan)
    do j=1,pgl(ichan)%multichan%n_unique_channels
       i=pgl(ichan)%multichan%unique_channel_list(j)
       if (i.eq.ichan) then
          ii=j
          cycle
       endif
       call pgl(i)%phase_space%compute_x_from_momenta(p)
       if (pgl(i)%phase_space%jac.lt.0d0) then
          ! The x's could not be correctly computed from the momenta
          write (*,*) 'WARNING: multi-channel weight not included'
          weight(1:pgl(ichan)%nproc)=1d0/dble(pgl(ichan)%multichan%number_of_channels(1:pgl(ichan)%nproc))
          return
       endif
       call mint_get_jacobian_from_x(i,pgl(i)%phase_space%x,vol)
       factors(j)=pgl(i)%phase_space%jac*vol
    enddo
    do i=1,pgl(ichan)%multichan%n_unique_channelgroups
       weight_factors(i)=1d0
       do j=1,pgl(ichan)%multichan%unique_channelgroup_list(0,i)
          if (pgl(ichan)%multichan%unique_channelgroup_list(j,i).eq.ii) cycle
          weight_factors(i)=weight_factors(i)+jac*vol_ichan/factors(pgl(ichan)%multichan%unique_channelgroup_list(j,i))
       enddo
    enddo
    weight(1:pgl(ichan)%nproc)=1d0/weight_factors(pgl(ichan)%multichan%map_proc_to_channelgroup(1:pgl(ichan)%nproc))
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
          if (all(all_unique_chans_inv(all_unique_channelgroups(1:pgl%multichan%number_of_channels(iproc),i)).eq. &
               pgl%multichan%channels(1:pgl%multichan%number_of_channels(iproc),iproc))) then
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
    deallocate(pgl%multichan%channels)
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
       if (is_singlet(part(i))) n_ps=n_ps+1
    enddo
    n_ps=factorial(n_ps)
  end subroutine determine_multi_channel_size
  


end module multichannel
