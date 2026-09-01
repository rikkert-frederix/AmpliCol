program mg_checks_regression
  use mg_checks, only: perform_check
  use handling_processes, only: pgl,ngroups
  use run_parameters, only: reset_run_parameters,alphaEW
  implicit none
  character(len=32) :: mode

  mode='pass'
  if (command_argument_count().ge.1) call get_command_argument(1,mode)
  call reset_run_parameters()
  ngroups=1
  allocate(pgl(1))
  pgl(1)%next=4
  pgl(1)%nproc=1
  allocate(pgl(1)%amp2(1),pgl(1)%iden(1),pgl(1)%amps(1))
  allocate(pgl(1)%amps(1)%n_sing(1))
  pgl(1)%amp2=1d0
  pgl(1)%iden=1_8
  pgl(1)%amps(1)%n_sing=0

  select case (trim(mode))
  case ('pass')
     call perform_check(1,1)
     write(*,'(a)') 'MadGraph comparison regression: PASS'
  case ('zero')
     pgl(1)%amp2=0d0
     call perform_check(1,1)
     write(*,'(a)') 'MadGraph zero comparison regression: PASS'
  case ('mismatch','bad-reference','bad-momentum','bad-count','truncated','trailing')
     call perform_check(1,1)
     error stop 'invalid MadGraph reference unexpectedly succeeded'
  case ('coupling-overflow')
     pgl(1)%amps(1)%n_sing=2
     alphaEW=huge(1d0)/100d0
     call perform_check(1,1)
     error stop 'overflowing coupling normalization unexpectedly succeeded'
  case ('overflow')
     pgl(1)%amp2=huge(1d0)
     call perform_check(1,1)
     error stop 'overflowing MadGraph normalization unexpectedly succeeded'
  case ('negative-code')
     pgl(1)%amp2=-1d0
     call perform_check(1,1)
     error stop 'negative matrix element unexpectedly succeeded'
  case ('bad-next')
     pgl(1)%next=21
     call perform_check(1,1)
     error stop 'invalid process dimension unexpectedly succeeded'
  case default
     error stop 'unknown MadGraph regression mode'
  end select
end program mg_checks_regression
