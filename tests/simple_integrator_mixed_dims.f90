program simple_integrator_mixed_dims_test
  use simple_integrator_mod
  implicit none
  integer, parameter :: nchannel=2,naux=2
  type(integrator) :: integ
  integer :: ichan,iint
  integer :: ndim(nchannel),ndim_extra(nchannel),nintegral(nchannel)
  real(kind=8) :: f(1),f_abs(1),aux(naux,1)
  real(kind=8),allocatable :: result(:,:),unc(:,:),aux_result(:,:),aux_unc(:,:)
  logical :: to_write(1),done

  ndim=(/1,2/)
  ndim_extra=(/0,1/)
  nintegral=1
  call integ%init(nchannel,ndim,ndim_extra,nintegral,100,1,accuracy=0.9d0,naux=naux)
  done=.false.
  do while (.not.done)
     call integ%get_points(1,ichan,iint)
     if (size(integ%x,1).ne.ndim(ichan)+ndim_extra(ichan)) then
        write(*,*) 'FAIL: mixed channel dimension'
        stop 1
     endif
     f=1d0
     f_abs=1d0
     aux(:,1)=(/2d0,3d0/)
     call integ%fill_points(1,f_abs,f,to_write,done,f_aux=aux)
  enddo
  call integ%get_channel_results(result,unc)
  call integ%get_channel_aux_results(aux_result,aux_unc)
  call assert_all_close(result(2,:),1d0,'primary result')
  call assert_all_close(aux_result(1,:),2d0,'first auxiliary result')
  call assert_all_close(aux_result(2,:),3d0,'second auxiliary result')
  write(*,'(a)') 'simple integrator mixed-dimension test: PASS'

contains

  subroutine assert_all_close(values,expected,label)
    real(kind=8),intent(in) :: values(:),expected
    character(len=*),intent(in) :: label
    if (any(abs(values-expected).gt.1d-12)) then
       write(*,*) 'FAIL: ',trim(label),values,expected
       stop 1
    endif
  end subroutine assert_all_close

end program simple_integrator_mixed_dims_test
