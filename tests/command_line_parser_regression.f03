program command_line_parser_regression
  use argument_parser
  implicit none
  character(len=80) :: filename,library,tag,timing
  character(len=256) :: input_file
  integer :: ncalls0,itmax,PS_choice,timing_sample
  integer(kind=8) :: seed
  logical :: combine_subprocesses

  call parse_argument(filename,input_file,ncalls0,itmax,PS_choice,seed,library,tag,&
       combine_subprocesses,timing,timing_sample)
  write (*,'(l1)') combine_subprocesses
end program command_line_parser_regression
