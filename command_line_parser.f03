module argument_parser
  implicit none
contains
  subroutine parse_argument(filename,ncalls0,itmax,PS_choice,seed,library,tag,read_momenta,me_points,&
       limit_test,timing,timing_sample,accuracy,alpha_dipole)
    integer :: i
    character(len=256) :: arg
    character(len=256) :: input_file,tmp
    logical :: verbose, show_help
    character(len=80) :: filename,library,tag,timing
    integer :: ncalls0,itmax,PS_choice
    integer(kind=8) :: seed
    logical :: read_momenta,limit_test
    integer :: me_points,timing_sample
    real(kind=8) :: accuracy,alpha_dipole(4)

    ! Default values:
    show_help=.false.
    filename='processes.txt'
    ncalls0=10000
    PS_choice=1
    seed=0
    itmax=128
    library='none'
    tag=''
    read_momenta=.false.
    limit_test=.false.
    timing='basic'
    timing_sample=100
    accuracy=0d0
    alpha_dipole=(/1d0,1d0,1d0,1d0/)

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
       elseif (arg.eq."--limit_test" .or. arg.eq."-lt") then
          limit_test=.true.
       elseif (index(arg, "--timing=").eq.1) then
          timing = arg(index(arg, "=")+1:)
       elseif (index(arg, "--timing-sample=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) timing_sample
       elseif (index(arg, "--accuracy=").eq.1 .or. index(arg, "-a=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) accuracy
          if (accuracy.le.0d0 .or. accuracy.ge.1d0) then
             write (*,*) 'Accuracy must be between 0 and 1: ',accuracy
             stop 1
          endif
       elseif (index(arg, "--alpha=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call parse_alpha(tmp,alpha_dipole)
       else
          write (*,*) 'Unknown argument: ',arg
          stop 1
       endif
    end do

    if (show_help) then
       write (*,'(a)') ""
       write (*,'(a)') "Usage: 'amplicol_generate <arguments>'. Possible arguments are"
       write (*,'(a)') ""
       write (*,'(a)') "  --help,           -h      : Show this message."
       write (*,'(a)') "  --process=[X],    -p=[X]  : Process specified in file [X] (default is './processes.txt')."
       write (*,'(a)') "  --nevents=[X],    -n=[X]  : Number of unweighted events to generate (default is 10000)."
       write (*,'(a)') "  --phasespace=[X], -ps=[X] : Phase-space parametrisation to use "// &
            "-- 1=gen23 (default), 2=HAAG, 3=pT-based, 4=t-channel."
       write (*,'(a)') "  --seed=[X],       -s=[X]  : The random number seed to use (default is read from randinit file)."
       write (*,'(a)') "  --itmax=[X],      -i=[X]  : The maximum number of iterations to use (detault is 128)."
       write (*,'(a)') "  --library=[X],    -l=[X]  : To create or use a library for the amplitudes, "// &
            "set [X] to 'create' or 'use', respectively. (To use a library, re-compile code with 'make "// &
            "amplicol_generate_library' after a library has been created). Default is 'none'."
       write (*,'(a)') "  --tag=[X],        -t=[X]  : Event file (and log file) names will be prepended with with a tag '[X]_'."
       write (*,'(a)') "  --me_test=[X],    -mt=[X] : Perform ME level test against MG "//& 
            "with [X] points tested (single PS kinematics)"
       write (*,'(a)') "  --limit_test,     -lt      : Test soft and collinear CS limits and exit."
       write (*,'(a)') "  --timing=[X]              : Timing mode: none, basic (default), or detailed."
       write (*,'(a)') "  --timing-sample=[X]       : In detailed timing, sample point timers every [X] points. Default is 100."
       write (*,'(a)') "  --accuracy=[X],   -a=[X]  : Disable event generation and integrate until "//&
            "the relative error is below [X] (0 < X < 1)."
       write (*,'(a)') "  --alpha=[X]               : Real-dipole restriction; one value or four comma-separated values FF,FI,IF,II."
       write (*,'(a)') ""
       stop
    end if
  end subroutine parse_argument

  subroutine parse_alpha(text,alpha)
    implicit none
    character(len=*),intent(in) :: text
    real(kind=8),intent(out) :: alpha(4)
    real(kind=8) :: value
    integer :: i,ios,ncomma

    ncomma=0
    do i=1,len_trim(text)
       if (text(i:i).eq.',') ncomma=ncomma+1
    enddo
    if (ncomma.eq.0) then
       read(text,*,iostat=ios) value
       if (ios.ne.0) then
          write (*,*) 'Invalid --alpha value: ',trim(text)
          stop 1
       endif
       alpha=value
    elseif (ncomma.eq.3) then
       read(text,*,iostat=ios) alpha
       if (ios.ne.0) then
          write (*,*) 'Invalid --alpha list: ',trim(text)
          stop 1
       endif
    else
       write (*,*) '--alpha requires one value or four comma-separated values: ',trim(text)
       stop 1
    endif
    do i=1,4
       if (.not.(alpha(i).gt.0d0 .and. alpha(i).le.1d0)) then
          write (*,*) '--alpha values must satisfy 0 < alpha <= 1: ',alpha(i)
          stop 1
       endif
    enddo
  end subroutine parse_alpha
end module argument_parser
