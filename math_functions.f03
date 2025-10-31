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
  subroutine get_next_iperm(ip,ips_in,ips,n)
    ! Given a permutation ips_in, find the next one and return it through ips.
    ! For example for ip=3 (length of permutation list), n=4 (elements to be
    ! considered in the permutation) this gives
    !
    !    ips_in        ips
    !-------------------------
    !    1,2,3   -->   1,2,4
    !    1,2,4   -->   1,3,2
    !    1,3,2   -->   1,3,4
    !    1,3,4   -->   1,4,2
    !    1,4,2   -->   1,4,3
    !    1,4,3   -->   2,1,3
    !    2,1,3   -->   2,1,4
    !    2,1,4   -->   2,3,1
    !    2,3,1   -->   2,3,4
    !    2,3,4   -->   2,4,1
    !    2,4,1   -->   2,4,3
    !    2,4,3   -->   3,1,2
    !    3,1,2   -->   3,1,4
    !    3,1,4   -->   3,2,1
    !    3,2,1   -->   3,2,4
    !    3,2,4   -->   3,4,1
    !    3,4,1   -->   3,4,2
    !    3,4,2   -->   4,1,2
    !    4,1,2   -->   4,1,3
    !    4,1,3   -->   4,2,1
    !    4,2,1   -->   4,2,3
    !    4,2,3   -->   4,3,1
    !    4,3,1   -->   4,3,2
    !    4,3,2   -->   XXXXX
    !
    ! Note that when giving non-sensical inputs (e.g., the last one in the
    ! list above), the code either goes into an infinite loop, or returns some
    ! bogus result. There is no check on the consistency of the input.
    implicit none
    integer :: ip,n,i_up,i,j
    integer,dimension(ip) :: ips,ips_in
    logical :: found
    found=.false.
    ips(1:ip)=ips_in(1:ip)
    do i_up=ip,1,-1
       do while (ips(i_up).lt.n)
          ips(i_up)=ips(i_up)+1
          if (any(ips(1:i_up-1).eq.ips(i_up))) cycle
          found=.true.
          exit
       enddo
       if (found) exit
    enddo
    do i=i_up+1,ip
       do j=1,n
          if (any(ips(1:i).eq.j)) then
             continue
          else
             ips(i)=j
             exit
          endif
       enddo
    enddo
  end subroutine get_next_iperm

  integer(kind=4) function ifindloc(a,n,i)
    ! returns the location of the value 'i' in the array 'a' (of size 'n'). If
    ! the array 'a' does not contain 'i', the value 'n+1' is returned.
    implicit none
    integer,intent(in) :: i,n
    integer,dimension(n),intent(in) :: a
    do ifindloc=1,n
       if (a(ifindloc).eq.i) return
    enddo
  end function ifindloc

end module math_functions
