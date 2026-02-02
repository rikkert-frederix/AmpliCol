module amp_lib
  implicit none
contains
  subroutine evaluate_amp(ichan,iint,p,amps)
    implicit none
    integer :: ichan,iint
    real(kind=8),dimension(*) :: p
    complex(kind=8),dimension(*) :: amps
    write (*,*) 'This subroutine should not be used'
    write (*,*) 'evaluate_amp() from the dummy.f file'
    stop 1
  end subroutine evaluate_amp
end module amp_lib
