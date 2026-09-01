module argument_parser
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer(kind=8),parameter :: maximum_ranmar_seed=30081_8*30081_8
contains
  subroutine parse_argument(filename,real_filename,input_file,ncalls0,itmax,PS_choice,seed,library,tag,read_momenta,me_points,&
       limit_test,timing,timing_sample,accuracy,dim_reg_scheme,has_real_process,&
       tail_replay_file,replay_tail,migration_tail_fraction_limit)
    integer :: i,arg_length,arg_status
    character(len=256) :: arg
    character(len=*) :: filename,real_filename
    character(len=256) :: input_file,tmp,tail_replay_file
    logical :: show_help
    character(len=80) :: library,tag,timing,dim_reg_scheme
    integer :: ncalls0,itmax,PS_choice
    integer(kind=8) :: seed
    logical :: read_momenta,limit_test,has_real_process,replay_tail
    integer :: me_points,timing_sample
    real(kind=8) :: accuracy,migration_tail_fraction_limit

    ! Default values:
    show_help=.false.
    filename='processes.txt'
    real_filename=''
    input_file='run_card.dat'
    ncalls0=10000
    PS_choice=1
    seed=0
    itmax=128
    library='none'
    tag=''
    read_momenta=.false.
    me_points=0
    limit_test=.false.
    has_real_process=.false.
    tail_replay_file=''
    replay_tail=.false.
    migration_tail_fraction_limit=0.2d0
    timing='basic'
    timing_sample=100
    accuracy=0d0
    dim_reg_scheme='hv'

    do i = 1, command_argument_count()
       arg_length=-1
       arg_status=-1
       call get_command_argument(i,length=arg_length,status=arg_status)
       if (arg_status.ne.0) then
          write (*,*) 'Could not query command-line argument at position',i
          stop 1
       endif
       if (arg_length.gt.len(arg)) then
          write (*,*) 'Command-line argument is too long at position',i,arg_length
          stop 1
       endif
       arg=''
       arg_status=-1
       call get_command_argument(i,arg,status=arg_status)
       if (arg_status.ne.0) then
          write (*,*) 'Could not read command-line argument at position',i
          stop 1
       endif
       arg = trim(arg)
       if (arg.eq."--help" .or. arg.eq."-h") then
          show_help = .true.
       elseif (index(arg, "--process=").eq.1 .or. index(arg, "-p=").eq.1) then
          call read_text_option(arg(index(arg,"=")+1:),'process file',filename)
       elseif (index(arg, "--real-process=").eq.1) then
          call read_text_option(arg(index(arg,"=")+1:),'real-process file',real_filename)
          has_real_process=.true.
       elseif (index(arg, "--input=").eq.1 .or. index(arg, "--card=").eq.1) then
          call read_text_option(arg(index(arg,"=")+1:),'input card',input_file)
       elseif (index(arg, "--nevents=").eq.1 .or. index(arg, "-n=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call read_integer_option(tmp,'nevents',ncalls0)
          if (ncalls0.lt.1) then
             write (*,*) 'Number of events must be at least 1: ',ncalls0
             stop 1
          endif
       elseif (index(arg, "--seed=").eq.1 .or. index(arg, "-s=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call read_integer8_option(tmp,'seed',seed)
          if (seed.lt.0_8 .or. seed.gt.maximum_ranmar_seed) then
             write (*,*) 'Seed must be between 0 and ',maximum_ranmar_seed,': ',seed
             stop 1
          endif
       elseif (index(arg, "--phasespace=").eq.1 .or. index(arg, "-ps=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call read_integer_option(tmp,'phasespace',PS_choice)
          if (PS_choice.lt.1 .or. PS_choice.gt.4) then
             write (*,*) 'Phase-space mode must be 1, 2, 3, or 4: ',PS_choice
             stop 1
          endif
       elseif (index(arg, "--itmax=").eq.1 .or. index(arg, "-i=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call read_integer_option(tmp,'itmax',itmax)
          if (itmax.lt.1) then
             write (*,*) 'Maximum iterations must be at least 1: ',itmax
             stop 1
          endif
       elseif (index(arg, "--library=").eq.1 .or. index(arg, "-l=").eq.1) then
          call read_text_option(arg(index(arg,"=")+1:),'library mode',library)
       elseif (index(arg, "--tag=").eq.1 .or. index(arg, "-t=").eq.1) then
          call read_tag_option(arg(index(arg,"=")+1:),tag)
       elseif (index(arg, "--me_test=").eq.1 .or. index(arg, "-mt=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call read_integer_option(tmp,'me_test',me_points)
          if (me_points.lt.1) then
             write (*,*) 'Number of matrix-element test points must be at least 1: ',me_points
             stop 1
          endif
          read_momenta=.true.
       elseif (arg.eq."--limit_test" .or. arg.eq."-lt") then
          limit_test=.true.
       elseif (arg.eq."--subtracted-real") then
          write (*,*) '--subtracted-real has been replaced by --real-process=FILE'
          stop 1
       elseif (index(arg, "--timing=").eq.1) then
          call read_text_option(arg(index(arg,"=")+1:),'timing mode',timing)
       elseif (index(arg, "--timing-sample=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call read_integer_option(tmp,'timing-sample',timing_sample)
          if (timing_sample.lt.1) then
             write (*,*) 'Timing sample must be at least 1: ',timing_sample
             stop 1
          endif
       elseif (index(arg, "--accuracy=").eq.1 .or. index(arg, "-a=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call read_real_option(tmp,'accuracy',accuracy)
          if (.not.ieee_is_finite(accuracy)) then
             write (*,*) 'Accuracy must be between 0 and 1: ',accuracy
             stop 1
          endif
          if (accuracy.le.0d0 .or. accuracy.ge.1d0) then
             write (*,*) 'Accuracy must be between 0 and 1: ',accuracy
             stop 1
          endif
       elseif (index(arg, "--dim-reg=").eq.1) then
          call read_text_option(arg(index(arg,"=")+1:),'dimensional-regularization scheme',dim_reg_scheme)
       elseif (index(arg, "--tail-replay=").eq.1) then
          call read_text_option(arg(index(arg,"=")+1:),'tail replay file',tail_replay_file)
          replay_tail=.true.
       elseif (index(arg, "--migration-tail-fraction=").eq.1) then
          tmp=arg(index(arg, "=")+1:)
          call read_real_option(tmp,'migration-tail-fraction',migration_tail_fraction_limit)
          if (.not.ieee_is_finite(migration_tail_fraction_limit)) then
             write(*,*) 'Migration tail fraction must satisfy 0 <= value <= 1: ',migration_tail_fraction_limit
             stop 1
          endif
          if (migration_tail_fraction_limit.lt.0d0 .or. migration_tail_fraction_limit.gt.1d0) then
             write(*,*) 'Migration tail fraction must satisfy 0 <= value <= 1: ',migration_tail_fraction_limit
             stop 1
          endif
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
       write (*,'(a)') "  --process=[X],    -p=[X]  : Born process specified in file [X] (default is './processes.txt')."
       write (*,'(a)') "  --real-process=[X]         : Real-emission process file; integrates B + R-sum(D) in one run."
       write (*,'(a)') "  --input=[X], --card=[X]   : Physics/run input card (default is './run_card.dat')."
       write (*,'(a)') "  --nevents=[X],    -n=[X]  : Number of unweighted events to generate; ignored with"// &
            " --accuracy (default is 10000)."
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
            "the absolute-envelope relative error is below [X] (0 < X < 1)."
       write (*,'(a)') "  --dim-reg=[hv|fdh]        : Dimensional scheme for integrated dipoles (default: hv)."
       write (*,'(a)') "  --tail-replay=[FILE]      : Re-evaluate the saved maximum-weight point in FILE and exit."
       write (*,'(a)') "  --migration-tail-fraction=[X] : Require the largest migration point to contribute"// &
            " at most X of its iteration variance before accuracy-based stopping (default 0.2; 0 disables)."
       write (*,'(a)') ""
       stop
    end if
  end subroutine parse_argument

  subroutine read_integer_option(raw,name,value)
    implicit none
    character(len=*),intent(in) :: raw,name
    integer,intent(out) :: value
    character(len=:),allocatable :: token
    integer :: i,ios
    token=trim(raw)
    if (len(token).eq.0) call invalid_numeric_option(name,raw)
    do i=1,len(token)
       if (i.eq.1 .and. (token(i:i).eq.'+' .or. token(i:i).eq.'-')) cycle
       if (token(i:i).lt.'0' .or. token(i:i).gt.'9') call invalid_numeric_option(name,raw)
    enddo
    if (token.eq.'+' .or. token.eq.'-') call invalid_numeric_option(name,raw)
    read(token,*,iostat=ios) value
    if (ios.ne.0) call invalid_numeric_option(name,raw)
  end subroutine read_integer_option

  subroutine read_integer8_option(raw,name,value)
    implicit none
    character(len=*),intent(in) :: raw,name
    integer(kind=8),intent(out) :: value
    character(len=:),allocatable :: token
    integer :: i,ios
    token=trim(raw)
    if (len(token).eq.0) call invalid_numeric_option(name,raw)
    do i=1,len(token)
       if (i.eq.1 .and. (token(i:i).eq.'+' .or. token(i:i).eq.'-')) cycle
       if (token(i:i).lt.'0' .or. token(i:i).gt.'9') call invalid_numeric_option(name,raw)
    enddo
    if (token.eq.'+' .or. token.eq.'-') call invalid_numeric_option(name,raw)
    read(token,*,iostat=ios) value
    if (ios.ne.0) call invalid_numeric_option(name,raw)
  end subroutine read_integer8_option

  subroutine read_real_option(raw,name,value)
    implicit none
    character(len=*),intent(in) :: raw,name
    real(kind=8),intent(out) :: value
    character(len=:),allocatable :: token
    character(len=*),parameter :: allowed='0123456789+-.eEdD'
    integer :: i,ios
    token=trim(raw)
    if (len(token).eq.0) call invalid_numeric_option(name,raw)
    do i=1,len(token)
       if (index(allowed,token(i:i)).eq.0) call invalid_numeric_option(name,raw)
    enddo
    value=0d0
    read(token,*,iostat=ios) value
    if (ios.ne.0) call invalid_numeric_option(name,raw)
    if (.not.ieee_is_finite(value)) call invalid_numeric_option(name,raw)
  end subroutine read_real_option

  subroutine invalid_numeric_option(name,raw)
    implicit none
    character(len=*),intent(in) :: name,raw
    write (*,*) 'Invalid numeric value for --',trim(name),': ',trim(raw)
    stop 1
  end subroutine invalid_numeric_option

  subroutine read_text_option(raw,name,value)
    implicit none
    character(len=*),intent(in) :: raw,name
    character(len=*),intent(out) :: value
    integer :: n
    n=len_trim(raw)
    if (n.lt.1) then
       write (*,*) 'A non-empty value is required for --',trim(name)
       stop 1
    endif
    if (n.gt.len(value)) then
       write (*,*) 'Value for --',trim(name),' is too long; maximum length is',len(value)
       stop 1
    endif
    value=raw(1:n)
  end subroutine read_text_option

  subroutine read_tag_option(raw,value)
    implicit none
    character(len=*),intent(in) :: raw
    character(len=*),intent(out) :: value
    integer :: n
    n=len_trim(raw)
    if (n.lt.1) then
       write (*,*) 'A non-empty value is required for --tag'
       stop 1
    endif
    if (n+1.gt.len(value)) then
       write (*,*) 'Tag is too long; maximum length is',len(value)-1
       stop 1
    endif
    if (index(raw(1:n),'/').ne.0 .or. index(raw(1:n),'\').ne.0) then
       write (*,*) 'Tag cannot contain a path separator: ',raw(1:n)
       stop 1
    endif
    value=raw(1:n)//'_'
  end subroutine read_tag_option

end module argument_parser
