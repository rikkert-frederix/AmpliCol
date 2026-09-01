module math_functions
  implicit none
contains
  integer function factorial(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    if (ifact.lt.0) then
       write (*,*) 'ERROR: factorial argument cannot be negative:',ifact
       stop 1
    endif
    factorial=1
    do i=2,ifact
       if (factorial.gt.huge(factorial)/i) then
          write (*,*) 'ERROR: factorial does not fit in default integer:',ifact
          stop 1
       endif
       factorial=factorial*i
    enddo
  end function factorial

  integer(kind=8) function factorial8(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    if (ifact.lt.0) then
       write (*,*) 'ERROR: factorial8 argument cannot be negative:',ifact
       stop 1
    endif
    factorial8=1_8
    do i=2,ifact
       if (factorial8.gt.huge(factorial8)/int(i,kind=8)) then
          write (*,*) 'ERROR: factorial does not fit in 64-bit integer:',ifact
          stop 1
       endif
       factorial8=factorial8*int(i,kind=8)
    enddo
  end function factorial8

  real(kind=8) function factorial_dble(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    if (ifact.lt.0) then
       write (*,*) 'ERROR: floating-point factorial argument cannot be negative:',ifact
       stop 1
    endif
    factorial_dble=1d0
    do i=2,ifact
       if (factorial_dble.gt.huge(factorial_dble)/dble(i)) then
          write (*,*) 'ERROR: factorial does not fit in double precision:',ifact
          stop 1
       endif
       factorial_dble=factorial_dble*dble(i)
    enddo
  end function factorial_dble

  integer function checked_multiply(first,second,context)
    implicit none
    integer,intent(in) :: first,second
    character(len=*),intent(in) :: context
    if (first.lt.0 .or. second.lt.0) then
       write (*,*) 'ERROR: negative integer factor in ',trim(context),first,second
       stop 1
    endif
    if (first.ne.0) then
       if (second.gt.huge(checked_multiply)/first) then
          write (*,*) 'ERROR: integer overflow in ',trim(context),first,second
          stop 1
       endif
    endif
    checked_multiply=first*second
  end function checked_multiply

  integer(kind=8) function checked_multiply8(first,second,context)
    implicit none
    integer(kind=8),intent(in) :: first,second
    character(len=*),intent(in) :: context
    if (first.lt.0_8 .or. second.lt.0_8) then
       write (*,*) 'ERROR: negative 64-bit integer factor in ',trim(context),first,second
       stop 1
    endif
    if (first.ne.0_8) then
       if (second.gt.huge(checked_multiply8)/first) then
          write (*,*) 'ERROR: 64-bit integer overflow in ',trim(context),first,second
          stop 1
       endif
    endif
    checked_multiply8=first*second
  end function checked_multiply8

  integer function checked_add(first,second,context)
    implicit none
    integer,intent(in) :: first,second
    character(len=*),intent(in) :: context
    if (first.lt.0 .or. second.lt.0) then
       write (*,*) 'ERROR: negative integer addend in ',trim(context),first,second
       stop 1
    endif
    if (second.gt.huge(checked_add)-first) then
       write (*,*) 'ERROR: integer overflow in ',trim(context),first,second
       stop 1
    endif
    checked_add=first+second
  end function checked_add

  integer(kind=8) function checked_add8(first,second,context)
    implicit none
    integer(kind=8),intent(in) :: first,second
    character(len=*),intent(in) :: context
    if (first.lt.0_8 .or. second.lt.0_8) then
       write (*,*) 'ERROR: negative 64-bit integer addend in ',trim(context),first,second
       stop 1
    endif
    if (second.gt.huge(checked_add8)-first) then
       write (*,*) 'ERROR: 64-bit integer overflow in ',trim(context),first,second
       stop 1
    endif
    checked_add8=first+second
  end function checked_add8

  integer function permutation_count(n,k)
    implicit none
    integer,intent(in) :: n,k
    integer :: i
    if (n.lt.0 .or. k.lt.0 .or. k.gt.n) then
       write (*,*) 'ERROR: invalid permutation-count arguments:',n,k
       stop 1
    endif
    permutation_count=1
    if (k.eq.0) return
    do i=n-k+1,n
       permutation_count=checked_multiply(permutation_count,i,'permutation count')
    enddo
  end function permutation_count

  integer function checked_integer_power(base,exponent,context)
    implicit none
    integer,intent(in) :: base,exponent
    character(len=*),intent(in) :: context
    integer :: remaining,factor
    if (base.lt.0 .or. exponent.lt.0) then
       write (*,*) 'ERROR: invalid integer power in ',trim(context),base,exponent
       stop 1
    endif
    if (exponent.eq.0) then
       checked_integer_power=1
       return
    endif
    if (base.eq.0 .or. base.eq.1) then
       checked_integer_power=base
       return
    endif
    checked_integer_power=1
    factor=base
    remaining=exponent
    do while (remaining.gt.0)
       if (btest(remaining,0)) checked_integer_power=&
            checked_multiply(checked_integer_power,factor,context)
       remaining=remaining/2
       if (remaining.gt.0) factor=checked_multiply(factor,factor,context)
    enddo
  end function checked_integer_power
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
    ! The final permutation has no successor and is rejected explicitly.
    implicit none
    integer :: ip,n,i_up,i,j
    integer,dimension(ip) :: ips,ips_in
    logical :: found
    if (ip.lt.1 .or. n.lt.1 .or. ip.gt.n) then
       write (*,*) 'ERROR: invalid partial-permutation dimensions:',ip,n
       stop 1
    endif
    if (any(ips_in.lt.1) .or. any(ips_in.gt.n)) then
       write (*,*) 'ERROR: partial permutation contains an out-of-range entry:',ips_in
       stop 1
    endif
    do i=1,ip-1
       if (any(ips_in(i+1:ip).eq.ips_in(i))) then
          write (*,*) 'ERROR: partial permutation contains a duplicate entry:',ips_in
          stop 1
       endif
    enddo
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
    if (.not.found) then
       write (*,*) 'ERROR: requested the successor of the final partial permutation:',ips_in
       stop 1
    endif
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
