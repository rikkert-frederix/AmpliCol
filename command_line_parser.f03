module argument_parser
  implicit none
contains
  subroutine parse_argument(filename,ncalls0,itmax,PS_choice,seed,library,tag,read_momenta,me_points)
    integer :: i
    character(len=256) :: arg
    character(len=256) :: input_file,tmp
    logical :: verbose, show_help
    character(len=80) :: filename,library,tag
    integer :: ncalls0,itmax,PS_choice
    integer(kind=8) :: seed
    logical :: read_momenta 
    integer :: me_points

    ! Default values:
    show_help=.false.
    filename='processes.txt'
    ncalls0=10000
    PS_choice=1
    seed=0
    itmax=48
    library='none'
    tag=''
    read_momenta=.false.

    do i = 1, command_argument_count()
       call get_command_argument(i, arg)
       arg = trim(arg)
       if (index(arg,"--help").eq.1 .or. index(arg,"-h").eq.1) then
          show_help = .true.
       elseif (index(arg, "--process=").eq.1 .or. index(arg, "-p=").eq.1) then
          filename = arg(index(arg, "=")+1:)
       elseif (index(arg, "--nevents=").eq.1 .or. index(arg, "-n=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) ncalls0
       elseif (index(arg, "--seed=").eq.1 .or. index(arg, "-s=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) seed
       elseif (index(arg, "--phasespace=").eq.1 .or. index(arg, "-ps=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) PS_choice
       elseif (index(arg, "--itmax=").eq.1 .or. index(arg, "-i=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) itmax
       elseif (index(arg, "--library=").eq.1 .or. index(arg, "-l=").eq.1) then
          library = arg(index(arg, "=")+1:)
       elseif (index(arg, "--tag=").eq.1 .or. index(arg, "-t=").eq.1) then
          tag = trim(arg(index(arg, "=")+1:))//'_'
       elseif (index(arg, "--me_test=").eq.1 .or. index(arg, "-mt=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) me_points
          read_momenta=.true.
       else
          write (*,*) 'Unknown argument: ',arg
          stop 1
       endif
    end do

    if (show_help) then
       write (*,'(a)') ""
       write (*,'(a)') "Usage: 'matrix_integrate_QCD <arguments>'. Possible arguments are"
       write (*,'(a)') ""
       write (*,'(a)') "  --help,           -h      : Show this message."
       write (*,'(a)') "  --process=[X],    -p=[X]  : Process specified in file [X] (default is './processes.txt')."
       write (*,'(a)') "  --nevents=[X],    -n=[X]  : Number of unweighted events to generate (default is 10000)."
       write (*,'(a)') "  --phasespace=[X], -ps=[X] : Phase-space parametrisation to use "// &
            "-- 1=gen23 (default), 2=HAAG, 3=pT-based, 4=t-channel."
       write (*,'(a)') "  --seed=[X],       -s=[X]  : The random number seed to use (default is read from randinit file)."
       write (*,'(a)') "  --itmax=[X],      -i=[X]  : The maximum number of iterations to use (detault is 48)."
       write (*,'(a)') "  --library=[X],    -l=[X]  : To create or use a library for the amplitudes, "// &
            "set [X] to 'create' or 'use', respectively. (To use a library, re-compile code with 'make "// &
            "matrix_integrate_library' after a library has been created). Default is 'none'."
       write (*,'(a)') "  --tag=[X],        -t=[X]  : Event file (and log file) names will be prepended with with a tag '[X]_'."
       write (*,'(a)') "  --me_test=[X],        -mt=[X]  : Perform ME level test against MG with [X] "// &
            "points tested (single PS kinematics)"
       write (*,'(a)') ""
       stop
    end if
  end subroutine parse_argument
end module argument_parser
