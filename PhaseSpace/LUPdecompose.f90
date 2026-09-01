module LUPdecomposition
! Implementation of the LUP decomposition. Can be used to efficiently
! 1. invert a matrix,
! 2. solve a linear system or
! 3. compute a determinant.
! Based on the C code that can be found here https://en.wikipedia.org/wiki/LU_decomposition
! Written January 2022.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: LUPdecompose,LUPsolve,LUPinvert,LUPdeterminant
contains
  subroutine LUPdecompose(a,n,tol,p,success)
! INPUT: A - array of pointers to rows of a square matrix having dimension N
!        Tol - small tolerance number to detect failure when the matrix is near degenerate
! OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
!        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
!        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
!        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    

    implicit none
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(0:n),intent(out) :: p
    real(kind=8),dimension(n,n),intent(inout) :: a
    logical,intent(out) :: success
    real(kind=8),intent(in) :: tol 
    integer(kind=4) :: i,j,k,imax
    real(kind=8) :: maxA,absA,multiplier,updated_entry
    real(kind=8),dimension(n) :: tmp
    success=.false.
    p=0
    if (n.le.0) return
    if (.not.ieee_is_finite(tol) .or. tol.lt.0d0 .or. &
         .not.all(ieee_is_finite(a))) return
    do i=1,n
       p(i)=i
    enddo
    p(0)=n
    do i=1,n
       maxA=0d0
       imax=i
       do k=i,n
          absA=abs(a(k,i))
          if (absA.gt.maxA) then
             maxA=absA
             imax=k
          endif
       enddo
       if (maxA.le.tol) then
          write (99,*) 'LUP decomposition failure: matrix is degenerate'
          return
       endif
       if (imax .ne. i ) then
          ! pivoting P
          j=p(i)
          p(i)=p(imax)
          p(imax)=j
          
          ! pivoting rows of a
          tmp(1:n)=a(i,1:n)
          a(i,1:n)=a(imax,1:n)
          a(imax,1:n)=tmp(1:n)
          
          ! counting pivots starting from n (for determinant)
          p(0)=p(0)+1
       endif
       do j=i+1,n
          a(j,i)=a(j,i)/a(i,i)
          if (.not.ieee_is_finite(a(j,i))) return
          do k=i+1,n
             if (.not.safe_lup_product(a(j,i),a(i,k),multiplier)) return
             if (.not.safe_lup_difference(a(j,k),multiplier,updated_entry)) return
             a(j,k)=updated_entry
          enddo
       enddo
    enddo
    success=.true.
  end subroutine LUPdecompose

  subroutine LUPSolve(a,p,b,n,x,success)
! INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
! OUTPUT: x - solution vector of A*x=b
    implicit none
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(0:n),intent(in) :: p
    real(kind=8),dimension(n,n),intent(in) :: a
    real(kind=8),dimension(n),intent(in) :: b
    real(kind=8),dimension(n),intent(out) :: x
    logical,intent(out),optional :: success
    integer(kind=4) i,k
    real(kind=8) :: term,updated_value

    x=0d0
    if (present(success)) success=.false.
    if (.not.valid_lup_inputs(a,p,n) .or. .not.all(ieee_is_finite(b))) return
    do i=1,n
       x(i)=b(p(i))
       do k=1,i-1
          if (.not.safe_lup_product(a(i,k),x(k),term)) then
             x=0d0
             return
          endif
          if (.not.safe_lup_difference(x(i),term,updated_value)) then
             x=0d0
             return
          endif
          x(i)=updated_value
       enddo
    enddo
    do i=n,1,-1
       do k=i+1,n
          if (.not.safe_lup_product(a(i,k),x(k),term)) then
             x=0d0
             return
          endif
          if (.not.safe_lup_difference(x(i),term,updated_value)) then
             x=0d0
             return
          endif
          x(i)=updated_value
       enddo
       if (.not.safe_lup_ratio(x(i),a(i,i),updated_value)) then
          x=0d0
          return
       endif
       x(i)=updated_value
    enddo
    if (present(success)) success=.true.
  end subroutine LUPSolve
  
  subroutine LUPinvert(a,p,n,ia,success)
! INPUT: A,P filled in LUPDecompose; N - dimension
! OUTPUT: IA is the inverse of the initial matrix
    implicit none
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(0:n),intent(in) :: p
    real(kind=8),dimension(n,n),intent(in) :: a
    real(kind=8),dimension(n,n),intent(out) :: ia
    logical,intent(out),optional :: success
    integer(kind=4) :: j
    real(kind=8),dimension(n) :: rhs,column
    logical :: column_success
    ia=0d0
    if (present(success)) success=.false.
    if (.not.valid_lup_inputs(a,p,n)) return
    do j=1,n
       rhs=0d0
       rhs(j)=1d0
       call LUPsolve(a,p,rhs,n,column,column_success)
       if (.not.column_success) then
          ia=0d0
          return
       endif
       ia(:,j)=column
    enddo
    if (present(success)) success=.true.
  end subroutine LUPinvert

  subroutine LUPdeterminant(a,p,n,det,success)
! INPUT: A,P filled in LUPDecompose; N - dimension. 
! OUTPUT: det returns the determinant of the initial matrix
    implicit none
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(0:n),intent(in) :: p
    real(kind=8),dimension(n,n),intent(in) :: a
    real(kind=8),intent(out) :: det
    logical,intent(out),optional :: success
    integer(kind=4) i
    real(kind=8) :: updated_det
    det=1d0
    if (present(success)) success=.false.
    if (.not.valid_lup_inputs(a,p,n)) then
       det=0d0
       return
    endif
    do i=1,n
       if (.not.safe_lup_product(det,a(i,i),updated_det)) then
          det=0d0
          return
       endif
       det=updated_det
    enddo
    if (mod(p(0)-n,2).ne.0) det=-det
    if (present(success)) success=.true.
  end subroutine LUPdeterminant

  logical function valid_lup_inputs(a,p,n) result(valid)
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(0:n),intent(in) :: p
    real(kind=8),dimension(n,n),intent(in) :: a
    integer :: i
    valid=.false.
    if (n.le.0 .or. .not.all(ieee_is_finite(a))) return
    if (p(0).lt.n .or. p(0).gt.2*n-1) return
    do i=1,n
       if (count(p(1:n).eq.i).ne.1) return
       if (a(i,i).eq.0d0) return
    enddo
    valid=.true.
  end function valid_lup_inputs

  logical function safe_lup_product(first,second,value) result(valid)
    real(kind=8),intent(in) :: first,second
    real(kind=8),intent(out) :: value
    valid=.false.
    value=0d0
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    if (first.eq.0d0 .or. second.eq.0d0) then
       valid=.true.
       return
    endif
    if (abs(second).gt.1d0) then
       if (abs(first).gt.huge(1d0)/abs(second)) return
    elseif (abs(first).gt.1d0) then
       if (abs(second).gt.huge(1d0)/abs(first)) return
    endif
    value=first*second
    valid=ieee_is_finite(value)
    if (.not.valid) value=0d0
  end function safe_lup_product

  logical function safe_lup_difference(first,second,value) result(valid)
    real(kind=8),intent(in) :: first,second
    real(kind=8),intent(out) :: value
    valid=.false.
    value=0d0
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    if (second.gt.0d0) then
       if (first.lt.-huge(1d0)+second) return
    elseif (second.lt.0d0) then
       if (first.gt.huge(1d0)+second) return
    endif
    value=first-second
    valid=ieee_is_finite(value)
    if (.not.valid) value=0d0
  end function safe_lup_difference

  logical function safe_lup_ratio(numerator,denominator,value) result(valid)
    real(kind=8),intent(in) :: numerator,denominator
    real(kind=8),intent(out) :: value
    valid=.false.
    value=0d0
    if (.not.ieee_is_finite(numerator) .or. .not.ieee_is_finite(denominator)) return
    if (denominator.eq.0d0) return
    if (abs(denominator).lt.1d0) then
       if (abs(numerator).gt.huge(1d0)*abs(denominator)) return
    endif
    value=numerator/denominator
    valid=ieee_is_finite(value)
    if (.not.valid) value=0d0
  end function safe_lup_ratio
end module LUPdecomposition
