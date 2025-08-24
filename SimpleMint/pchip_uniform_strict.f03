module pchip_uniform_strict
  use iso_fortran_env, only: real64
  implicit none
contains

  ! Compute PCHIP (Fritsch–Carlson) slopes for a UNIFORM grid.
  ! Input y(:) MUST be strictly increasing (y(i+1) > y(i) for all i).
  subroutine pchip_slopes_uniform_strict(y, m)
    real(real64), intent(in)  :: y(:)   ! can have any lower/upper bounds
    real(real64), intent(out) :: m(lbound(y,1):ubound(y,1))

    integer :: lb, ub, n, i
    real(real64), allocatable :: d(:)  ! secant slopes (uniform spacing => delta)
    real(real64) :: di, dim1, m0

    lb = lbound(y,1)
    ub = ubound(y,1)
    n  = ub - lb + 1
    if (n < 2) error stop "pchip: need at least two samples"

    allocate(d(lb:ub-1))
    do i = lb, ub-1
       d(i) = y(i+1) - y(i)
       if (d(i) <= 0.0_real64) then
          write (*,*) d(i),i,y(i+1),y(i)
          error stop "pchip: input must be strictly increasing"
       endif
    end do

    if (n == 2) then
       ! Only one interval: derivative equals the secant on both ends
       m(lb) = d(lb)
       m(ub) = d(lb)
       deallocate(d)
       return
    end if

    ! Interior slopes (Fritsch–Carlson; uniform grid => simple harmonic mean)
    do i = lb+1, ub-1
       dim1 = d(i-1); di = d(i)
       ! Both positive by strict-increasing check → harmonic mean is positive
       m(i) = 2.0_real64 * dim1 * di / (dim1 + di)
    end do

    ! Endpoint slopes (PCHIP one-sided formulas, uniform spacing)
    ! Left
    m0 = 0.5_real64 * (3.0_real64*d(lb) - d(lb+1))
    if (m0 < 0.0_real64) then
       m(lb) = 0.0_real64  ! would create a tiny flat start, but only if data near-flat
    elseif (m0 > 3.0_real64*d(lb)) then
       m(lb) = 3.0_real64*d(lb)
    else
       m(lb) = m0
    end if
    ! Right
    m0 = 0.5_real64 * (3.0_real64*d(ub-1) - d(ub-2))
    if (m0 < 0.0_real64) then
       m(ub) = 0.0_real64
    elseif (m0 > 3.0_real64*d(ub-1)) then
       m(ub) = 3.0_real64*d(ub-1)
    else
       m(ub) = m0
    end if

    deallocate(d)
  end subroutine pchip_slopes_uniform_strict


  ! Evaluate the cubic Hermite on a larger UNIFORM grid 0..size_new
  ! Original y(:) is on integer grid lb..ub (spacing 1). Endpoints are matched exactly.
  subroutine upsample_pchip_uniform(y, m, arr_new)
    real(real64), intent(in)  :: y(:)
    real(real64), intent(in)  :: m(lbound(y,1):ubound(y,1))
    real(real64), intent(out) :: arr_new(:)  ! expected bounds 0..size_new

    integer :: lb, ub, n_orig, j, size_new, k, lbn, ubn
    real(real64) :: pos, t
    real(real64) :: y0, y1, m0, m1
    real(real64) :: t2, t3, h00, h10, h01, h11

    lb = lbound(y,1)
    ub = ubound(y,1)
    n_orig = ub - lb

    lbn = lbound(arr_new,1)
    ubn = ubound(arr_new,1)
    size_new = ubn - lbn
    if (size_new < 1) then
       if (size_new == 0) then
          arr_new(lbn) = y(lb)
          return
       else
          error stop "upsample_pchip_uniform: size_new must be >= 0"
       end if
    end if

    ! Exact endpoints
    arr_new(lbn) = y(lb)
    arr_new(ubn) = y(ub)

    ! Map linearly from new integer grid to original continuous index
    do j = lbn+1, ubn-1
       pos = real(j - lbn, real64) / real(size_new, real64) * real(n_orig, real64) + real(lb, real64)

       k = int(floor(pos))
       if (k < lb) k = lb
       if (k > ub-1) k = ub-1
       t = pos - real(k, real64)

       y0 = y(k);     y1 = y(k+1)
       m0 = m(k);     m1 = m(k+1)

       t2 = t*t; t3 = t2*t
       h00 =  2.0_real64*t3 - 3.0_real64*t2 + 1.0_real64
       h10 =        t3 - 2.0_real64*t2 + t
       h01 = -2.0_real64*t3 + 3.0_real64*t2
       h11 =        t3 -       t2

       arr_new(j) = h00*y0 + h10*m0 + h01*y1 + h11*m1
    end do
  end subroutine upsample_pchip_uniform


  ! Convenience wrapper:
  ! Given arr(0:isize) strictly increasing with arr(0)=0 and arr(isize)=1,
  ! produce arr_new(0:size_new) (allocated) with strictly monotone cubic interpolation.
  subroutine resize_arr_pchip_strict(arr, size_new, arr_new)
    real(real64), intent(in)  :: arr(:)     ! e.g., bounds 0:isize
    integer,      intent(in)  :: size_new
    real(real64), intent(out) :: arr_new(0:size_new)

    real(real64), allocatable :: m(:)
    integer :: lb, ub

    lb = lbound(arr,1)
    ub = ubound(arr,1)

    if (arr(lb) /= 0.0_real64 .or. arr(ub) /= 1.0_real64) then
       error stop "resize_arr_pchip_strict: endpoints must be 0 and 1"
    end if

    allocate(m(lb:ub))
    call pchip_slopes_uniform_strict(arr, m)

    call upsample_pchip_uniform(arr, m, arr_new)

    ! No clamping needed; PCHIP is shape-preserving and stays within [0,1]
    ! If you want to guard against rare floating noise at 1 ulp, you could:
    ! arr_new(0)        = 0.0_real64
    ! arr_new(size_new) = 1.0_real64
  end subroutine resize_arr_pchip_strict

end module pchip_uniform_strict

