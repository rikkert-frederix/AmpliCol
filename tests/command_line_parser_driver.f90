program command_line_parser_driver
  use argument_parser
  implicit none
  character(len=80) :: filename,real_filename,library,tag,timing,dim_reg_scheme
  character(len=256) :: input_file,tail_replay_file
  integer :: ncalls0,itmax,ps_choice,me_points,timing_sample
  integer(kind=8) :: seed
  logical :: read_momenta,limit_test,has_real_process,replay_tail
  real(kind=8) :: accuracy,migration_tail_fraction_limit

  call parse_argument(filename,real_filename,input_file,ncalls0,itmax,ps_choice,seed,library,tag,&
       read_momenta,me_points,limit_test,timing,timing_sample,accuracy,dim_reg_scheme,&
       has_real_process,tail_replay_file,replay_tail,migration_tail_fraction_limit)

  write (*,'(a,1x,5(i0,1x),l1,1x,l1)') 'PARSER_OK',ncalls0,itmax,ps_choice,me_points,&
       timing_sample,read_momenta,replay_tail
end program command_line_parser_driver
