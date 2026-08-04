module argument_parser
  implicit none
contains
  subroutine parse_argument(filename,ncalls0,itmax,PS_choice,seed,library,tag,read_momenta,me_points,&
       amplicol_probe_points,amplicol_fixed_probe_points,amplicol_momenta_probe_points,&
       amplicol_probe_quiet,timing,timing_sample)
    integer :: i
    character(len=256) :: arg
    character(len=256) :: input_file,tmp
    logical :: verbose, show_help
    character(len=80) :: filename,library,tag,timing
    integer :: ncalls0,itmax,PS_choice
    integer(kind=8) :: seed
    logical :: read_momenta,amplicol_probe_quiet
    integer :: me_points,amplicol_probe_points,amplicol_fixed_probe_points
    integer :: amplicol_momenta_probe_points,timing_sample

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
    amplicol_probe_points=0
    amplicol_fixed_probe_points=0
    amplicol_momenta_probe_points=0
    amplicol_probe_quiet=.false.
    timing='basic'
    timing_sample=100

    call get_environment_variable("AMPICOL_PROBE_QUIET", tmp, status=i)
    if (i.eq.0) call parse_boolish(tmp,amplicol_probe_quiet)

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
       elseif (index(arg, "--amplicol_probe=").eq.1 .or. index(arg, "--amplicol-probe=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) amplicol_probe_points
       elseif (index(arg, "--amplicol_fixed_probe=").eq.1 .or. &
            index(arg, "--amplicol-fixed-probe=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) amplicol_fixed_probe_points
       elseif (index(arg, "--amplicol_momenta_probe=").eq.1 .or. &
            index(arg, "--amplicol-momenta-probe=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) amplicol_momenta_probe_points
       elseif (arg.eq."--amplicol_probe_quiet" .or. &
            arg.eq."--amplicol-probe-quiet") then
          amplicol_probe_quiet=.true.
       elseif (index(arg, "--amplicol_probe_quiet=").eq.1 .or. &
            index(arg, "--amplicol-probe-quiet=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          call parse_boolish(tmp,amplicol_probe_quiet)
       elseif (index(arg, "--timing=").eq.1) then
          timing = arg(index(arg, "=")+1:)
       elseif (index(arg, "--timing-sample=").eq.1) then
          tmp = arg(index(arg, "=")+1:)
          read(tmp,*) timing_sample
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
       write (*,'(a)') "  --amplicol_probe=[X]      : Print [X] direct AmpliCol ME values and kinematics "//&
            "without running MadGraph."
       write (*,'(a)') "  --amplicol_fixed_probe=[X]: Print [X] direct AmpliCol ME values at a fixed "//&
            "2-to-1 on-shell kinematic point. Currently supports q q~ > Z."
       write (*,'(a)') "  --amplicol_momenta_probe=[X]: Print direct AmpliCol ME values using momenta "//&
            "from Utilities/ME_checks/momenta_<group>_<integral>.txt, bypassing integration."
       write (*,'(a)') "  --amplicol_probe_quiet   : Suppress per-point direct-probe ME/momentum stdout. "//&
            "Can also be enabled with AMPICOL_PROBE_QUIET=1."
       write (*,'(a)') "  --timing=[X]              : Timing mode: none, basic (default), detailed, or numeric "//&
            "detailed timing sample."
       write (*,'(a)') "  --timing-sample=[X]       : In detailed timing, sample point timers every [X] points. Default is 100."
       write (*,'(a)') ""
       stop
    end if
  end subroutine parse_argument

  subroutine parse_boolish(value,result)
    implicit none
    character(len=*),intent(in) :: value
    logical,intent(out) :: result
    character(len=256) :: local
    local=trim(adjustl(value))
    result=.true.
    if (local.eq."" .or. local.eq."0" .or. local.eq."false" .or. local.eq."FALSE" .or. &
         local.eq."False" .or. local.eq."no" .or. local.eq."NO" .or. local.eq."No" .or. &
         local.eq."off" .or. local.eq."OFF" .or. local.eq."Off") then
       result=.false.
    endif
  end subroutine parse_boolish
end module argument_parser
