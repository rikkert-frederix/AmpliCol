program native_smoke
  use, intrinsic :: iso_c_binding
  use rusticol
  implicit none

  type(rusticol_runtime) :: runtime
  character(len=4096) :: process_dir
  real(c_double), target :: momenta(16)
  real(c_double), allocatable, target :: total(:)
  real(c_double), allocatable, target :: resolved(:, :, :)
  real(c_double) :: resolved_total
  integer(c_int) :: status

  if (command_argument_count() /= 1) then
    write(*, '(A)') "usage: native_smoke PROCESS_OUTPUT"
    stop 2
  end if
  call get_command_argument(1, process_dir)
  call runtime%load(trim(process_dir), ierr=status)
  if (status /= RUSTICOL_STATUS_OK) then
    write(*, '(A)') rusticol_last_error()
    stop 1
  end if

  momenta = [ &
      500.0_c_double, 0.0_c_double, 0.0_c_double, 500.0_c_double, &
      500.0_c_double, 0.0_c_double, 0.0_c_double, -500.0_c_double, &
      504.15762567199999_c_double, -304.10842628649999_c_double, &
      208.76026523528103_c_double, 331.35611794513767_c_double, &
      495.84237432800001_c_double, 304.10842628649999_c_double, &
      -208.76026523528103_c_double, -331.35611794513767_c_double]

  call runtime%evaluate(momenta, 1_c_size_t, total, ierr=status)
  if (status /= RUSTICOL_STATUS_OK) stop 1
  call runtime%evaluate_resolved(momenta, 1_c_size_t, resolved, ierr=status)
  if (status /= RUSTICOL_STATUS_OK) stop 1
  resolved_total = sum(resolved(:, :, 1))
  if (abs(total(1) - resolved_total) > 1.0e-12_c_double * abs(total(1))) then
    write(*, '(A)') "resolved sum does not reproduce the compatibility total"
    stop 1
  end if
  write(*, '(ES24.16,1X,ES24.16)') total(1), resolved_total
end program native_smoke
