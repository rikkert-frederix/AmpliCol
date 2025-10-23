module LUPdecomposition
! Implementation of the LUP decomposition. Can be used to efficiently
! 1. invert a matrix,
! 2. solve a linear system or
! 3. compute a determinant.
! Based on the C code that can be found here https://en.wikipedia.org/wiki/LU_decomposition
! Written January 2022.
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
    real(kind=8) :: maxA,absA
    real(kind=8),dimension(n) :: tmp
    success=.true.
    do i=1,n
       p(i)=i
    enddo
    p(0)=n
    do i=1,n
       maxA=0d0
       imax=i
       do k=i,n
          absA=abs(a(i,k))
          if (absA.gt.maxA) then
             maxA=absA
             imax=k
          endif
       enddo
       if (maxA .lt. tol) then
          write (99,*) 'LUP decomposition failure: matrix is degenerate'
          success=.false.
          return
       endif
       if (imax .ne. i ) then
          ! pivoting P
          j=p(i)
          p(i)=p(imax)
          p(imax)=j
          
          ! pivoting rows of a
          tmp(1:n)=a(1:n,i)
          a(1:n,i)=a(1:n,imax)
          a(1:n,imax)=tmp(1:n)
          
          ! counting pivots starting from n (for determinant)
          p(0)=p(0)+1
       endif
       do j=i+1,n
          A(i,j)=A(i,j)/A(i,i)
          do k=i+1,n
             A(k,j)=A(k,j)-A(i,j)*A(k,i)
          enddo
       enddo
    enddo
  end subroutine LUPdecompose

  subroutine LUPSolve(a,p,b,n,x)
! INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
! OUTPUT: x - solution vector of A*x=b
    implicit none
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(0:n),intent(in) :: p
    real(kind=8),dimension(n,n),intent(in) :: a
    real(kind=8),dimension(n),intent(in) :: b
    real(kind=8),dimension(n),intent(out) :: x
    integer(kind=4) i,k

    do i=1,n
       x(i)=b(p(i))
       do k=1,i
          x(i)=x(i)-a(k,i)*x(k)
       enddo
    enddo
    do i=n,1,-1
       do k=i+1,n
          x(i)=x(i)-a(k,i)*x(k)
       enddo
       x(i)=x(i)/a(i,i)
    enddo
  end subroutine LUPSolve
  
  subroutine LUPinvert(a,p,n,ia)
! INPUT: A,P filled in LUPDecompose; N - dimension
! OUTPUT: IA is the inverse of the initial matrix
    implicit none
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(0:n),intent(in) :: p
    real(kind=8),dimension(n,n),intent(in) :: a
    real(kind=8),dimension(n,n),intent(out) :: ia
    integer(kind=4) i,j,k
    do j=1,n
       do i=1,n
          if (p(i).eq.j) then
             ia(j,i)=1d0
          else
             ia(j,i)=0d0
          endif
          do k=1,i-1
             ia(j,i)=ia(j,i)-ia(k,i)*ia(j,k)
          enddo
       enddo
       do i=n,1,-1
          do k=i+1,n
             ia(j,i)=ia(j,i)-a(k,i)*ia(j,k)
          enddo
          ia(j,i)=ia(j,i)/ia(i,i)
       enddo
    enddo
  end subroutine LUPinvert

  subroutine LUPdeterminant(a,p,n,det)
! INPUT: A,P filled in LUPDecompose; N - dimension. 
! OUTPUT: det returns the determinant of the initial matrix
    implicit none
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(0:n),intent(in) :: p
    real(kind=8),dimension(n,n),intent(in) :: a
    real(kind=8),intent(out) :: det
    integer(kind=4) i
    det=1d0
    do i=1,n
       det=det*a(i,i)
    enddo
    if (mod(p(0)-n,2).ne.0) det=-det
  end subroutine LUPdeterminant
end module LUPdecomposition
