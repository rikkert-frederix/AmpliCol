module math_functions
contains
  integer function factorial(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    integer,save :: ifact_save=-1
    integer,dimension(:),allocatable,save :: factorial_save
    if (ifact.gt.ifact_save) then
       if (allocated(factorial_save)) deallocate(factorial_save)
       allocate(factorial_save(0:ifact))
       factorial_save(0)=1
       do i=1,ifact
          factorial_save(i)=factorial_save(i-1)*i
       enddo
       ifact_save=ifact
    endif
    factorial=factorial_save(ifact)
  end function factorial

  integer(kind=8) function factorial8(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    integer,save :: ifact_save=-1
    integer(kind=8),dimension(:),allocatable,save :: factorial_save
    if (ifact.gt.ifact_save) then
       if (allocated(factorial_save)) deallocate(factorial_save)
       allocate(factorial_save(0:ifact))
       factorial_save(0)=1
       do i=1,ifact
          factorial_save(i)=factorial_save(i-1)*i
       enddo
       ifact_save=ifact
    endif
    factorial8=factorial_save(ifact)
  end function factorial8

  real(kind=8) function factorial_dble(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    integer,save :: ifact_save=-1
    real(kind=8),dimension(:),allocatable,save :: factorial_save
    if (ifact.gt.ifact_save) then
       if (allocated(factorial_save)) deallocate(factorial_save)
       allocate(factorial_save(0:ifact))
       factorial_save(0)=1d0
       do i=1,ifact
          factorial_save(i)=factorial_save(i-1)*dble(i)
       enddo
       ifact_save=ifact
    endif
    factorial_dble=factorial_save(ifact)
  end function factorial_dble
end module math_functions
