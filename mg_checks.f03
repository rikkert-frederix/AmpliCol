module mg_checks
  use common
  use amplitude_QCD_mod
  use argument_parser
  use handling_processes
  integer :: k,nord
  real(kind=8),dimension(:),allocatable :: mg_check
  real(kind=8),dimension(100) :: me_code
  logical :: match
  logical :: printed_first_me_point=.false.
  logical :: printed_first_amplicol_probe_point=.false.
  logical :: me_test_done=.false.
  logical :: amplicol_probe_done=.false.
  real(kind=8),dimension(:,:),allocatable :: p_read
  real(kind=8),parameter  :: alpha_check=0.118d0
  real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
  integer :: me_points
  integer :: amplicol_probe_points=0
  integer :: amplicol_fixed_probe_points=0
  integer :: amplicol_momenta_probe_points=0
  logical :: amplicol_probe_quiet=.false.

contains

  subroutine get_madgraph_results(n,ichan,iint,me,nlines)
   implicit none
   integer :: i,k,n,io
   integer,intent(in) :: ichan,iint
   real(kind=8),dimension(:),allocatable :: me
   real(kind=8),dimension(0:3) :: dum
   character(len=50) :: filename
   integer :: nlines
   write(filename, '("Utilities/ME_checks/momenta_", I0, "_",I0,".txt")') ichan,iint
   open(20,file=trim(filename),status="old",iostat=io)
   if (io.ne.0) then
      write(*,*) 'Could not open MadGraph check file: ',trim(filename)
      stop 1
   endif
   do i=1,n
        read(20,*,iostat=io) dum
        if (io.ne.0) then
           write(*,*) 'Could not read momentum line from MadGraph check file: ',trim(filename)
           stop 1
        endif
   enddo
   read(20,*,iostat=io) nlines
   if (io.ne.0) then
      write(*,*) 'Could not read number of MadGraph matrix elements from: ',trim(filename)
      stop 1
   endif
   if (allocated(me)) deallocate(me)
   allocate(me(nlines))
   do i=1,nlines
      read(20,*,iostat=io) me(i)
      if (io.ne.0) then
         write(*,*) 'Could not read MadGraph matrix element from: ',trim(filename)
         stop 1
      endif
   enddo
   close(20)
  end subroutine get_madgraph_results

  subroutine read_in_momenta(n,igroup,iamp,p_in)
   implicit none
   integer :: i,k,n,io
   integer,intent(in) :: iamp,igroup
   real(kind=8),dimension(n,0:3),intent(out) :: p_in
   character(len=50) :: filename
   write(filename, '("Utilities/ME_checks/momenta_", I0, "_",I0,".txt")') igroup,iamp
   open(40,file=trim(filename),status="old",iostat=io)
   if (io.ne.0) then
      write(*,*) 'Could not open MadGraph momentum file: ',trim(filename)
      stop 1
   endif
   do i=1,n
        read(40,*,iostat=io) p_in(i,:)
        if (io.ne.0) then
           write(*,*) 'Could not read momentum line from MadGraph momentum file: ',trim(filename)
           stop 1
        endif
   enddo
   close(40)
  end subroutine read_in_momenta

  subroutine run_madgraph_check(n,igroup,iamp,list)
   implicit none
   integer,intent(in) :: iamp
   integer,intent(in) :: n,igroup
   integer,dimension(n),intent(in) :: list
   integer :: status
   character(len=100) :: command
   integer :: i
   write(command, '(A,I0)') "./Utilities/ME_checks/mg_wrapper.sh ", igroup
   write(command(len_trim(command)+1:), '(1X,I0)') iamp
   do i = 1, n
        write(command(len_trim(command)+1:), '(1X,I0)') list(i)
   end do
   call execute_command_line(command, exitstat=status)
    if (status /= 0) then
        print *, "MadGraph check script failed with status:", status
        stop 1
    else
        print *, "MG script finished successfully fori igroup, iamp: ",igroup,iamp
    end if
  end subroutine run_madgraph_check

  subroutine perform_check(iint,ichan)
    implicit none
    integer :: i,nord,iint,ichan,matched_index
    if (.not.allocated(mg_check)) allocate(mg_check(1000))
    me_code = pgl(ichan)%amp2(:)*(4*pi*alpha_check)**(pgl(ichan)%next-2-pgl(ichan)%amps(iint)%n_sing(1))&
               *(2d0*4d0*pi*alphaEW)**pgl(ichan)%amps(iint)%n_sing(1)/dble(pgl(ichan)%iden(iint))
    call get_madgraph_results(pgl(ichan)%next,ichan,iint,mg_check,nord)
    match=.false.
    matched_index=0
    do i=1,nord
       if (abs((mg_check(i)-me_code(1))/me_code(1)).lt.1d-4) then
          match=.true.
          if (matched_index.eq.0) matched_index=i
       endif
    enddo
    if (.not.printed_first_me_point) then
       write(*,*) ' '
       write(*,*) 'ME-check first phase-space point'
       write(*,*) '  group, integral:',ichan,iint
       write(*,'(a)') '     i      pdg                      E                     px                     py                     pz'
       do i=1,pgl(ichan)%next
          write(*,'(i6,1x,i8,4(1x,es31.23))') i,pgl(ichan)%processes(i,iint),pgl(ichan)%ps(1)%p(0:3,i)
       enddo
       write(*,'(a,1x,es31.23)') '  AmpliCol matrix element:',dble(me_code(1))
       write(*,'(a)') '  Matching MadGraph matrix element:'
       if (matched_index.gt.0) then
          write(*,'(i6,1x,es31.23)') matched_index,dble(mg_check(matched_index))
       else
          write(*,'(a)') '    none'
       endif
       write(*,*) ' '
       printed_first_me_point=.true.
    endif
    if (.not.match) then
       write(*,*) 'ERROR: disagreement with MG in ME-level check!',me_code(1)
       write(*,*) mg_check
       !stop 4
    endif
    if (pgl(ichan)%passed(iint).ge.me_points) then
            write(*,*) ' '
            write(*,*) '***** Passed all the', me_points,' ME-level tests. Stop the ME evaluation test.'
            me_test_done=.true.
    endif
  end subroutine perform_check

  subroutine perform_amplicol_probe(iint,ichan)
    implicit none
    integer :: i,iint,ichan,env_status
    character(len=16) :: debug_helicities
    if (amplicol_probe_points.le.0) return
    me_code = pgl(ichan)%amp2(:)*(4*pi*alpha_check)**(pgl(ichan)%next-2-pgl(ichan)%amps(iint)%n_sing(1))&
               *(2d0*4d0*pi*alphaEW)**pgl(ichan)%amps(iint)%n_sing(1)/dble(pgl(ichan)%iden(iint))
    call get_environment_variable('AMPICOL_PROBE_HELICITY_DEBUG',debug_helicities,status=env_status)
    if (env_status.eq.0 .and. trim(debug_helicities).eq.'1') then
       write(*,'(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,es31.23,1x,i0,1x,i0)') &
            'AMPICOL_PROBE_DEBUG_HEADER',pgl(ichan)%passed(iint),ichan,iint,&
            pgl(ichan)%amps(iint)%n_amps,dble(pgl(ichan)%col_fac(1)),&
            pgl(ichan)%iden(iint),pgl(ichan)%amps(iint)%n_sing(1)
       write(*,'(a,1x,i0,1x,i0,1x,i0,1x,es31.23,1x,es31.23,1x,es31.23)') &
            'AMPICOL_PROBE_DEBUG_NORMALIZATION',pgl(ichan)%passed(iint),&
            next,pgl(ichan)%next,alpha_check,alphaEW,&
            (4*pi*alpha_check)**(pgl(ichan)%next-2-pgl(ichan)%amps(iint)%n_sing(1))&
            *(2d0*4d0*pi*alphaEW)**pgl(ichan)%amps(iint)%n_sing(1)
       write(*,'(a,1x,i0,1x,i0,99(1x,i0))') &
            'AMPICOL_PROBE_DEBUG_IPROC_START',pgl(ichan)%passed(iint),&
            size(pgl(ichan)%amps(iint)%iproc_start),&
            pgl(ichan)%amps(iint)%iproc_start
       do i=1,size(pgl(ichan)%amp2)
          write(*,'(a,1x,i0,1x,i0,1x,es31.23,1x,es31.23)') &
               'AMPICOL_PROBE_DEBUG_PROC_AMP2',pgl(ichan)%passed(iint),i,&
               pgl(ichan)%amp2(i),me_code(i)
       enddo
       do i=1,pgl(ichan)%nproc
          write(*,'(a,1x,i0,1x,i0,99(1x,i0))') &
               'AMPICOL_PROBE_DEBUG_PROC_ROW',pgl(ichan)%passed(iint),i,&
               pgl(ichan)%processes(1:pgl(ichan)%next,i)
       enddo
       do i=1,pgl(ichan)%amps(iint)%n_amps
          write(*,'(a,1x,i0,1x,i0,1x,es31.23,1x,i0,99(1x,i0))') &
               'AMPICOL_PROBE_DEBUG_HEL',pgl(ichan)%passed(iint),i,&
               pgl(ichan)%amp2_hel(i),pgl(ichan)%hel_fac(i,iint),&
               pgl(ichan)%amps(iint)%spins(1:pgl(ichan)%next,1,i)
       enddo
    endif
    if (.not.amplicol_probe_quiet) then
       write(*,'(a,1x,i0,1x,i0,1x,i0,1x,es31.23)') 'AMPICOL_PROBE_VALUE',&
            pgl(ichan)%passed(iint),ichan,iint,dble(me_code(1))
       do i=1,pgl(ichan)%next
          write(*,'(a,1x,i0,1x,i0,1x,i0,4(1x,es31.23))') 'AMPICOL_PROBE_MOM',&
               pgl(ichan)%passed(iint),i,pgl(ichan)%processes(i,iint),pgl(ichan)%ps(1)%p(0:3,i)
       enddo
    endif
    if ((.not.amplicol_probe_quiet) .and. (.not.printed_first_amplicol_probe_point)) then
       write(*,*) ' '
       write(*,*) 'AmpliCol probe first phase-space point'
       write(*,*) '  group, integral:',ichan,iint
       write(*,'(a)') '     i      pdg                      E                     px                     py                     pz'
       do i=1,pgl(ichan)%next
          write(*,'(i6,1x,i8,4(1x,es31.23))') i,pgl(ichan)%processes(i,iint),pgl(ichan)%ps(1)%p(0:3,i)
       enddo
       write(*,'(a,1x,es31.23)') '  AmpliCol matrix element:',dble(me_code(1))
       write(*,*) ' '
       printed_first_amplicol_probe_point=.true.
    endif
    if (pgl(ichan)%passed(iint).ge.amplicol_probe_points) then
       write(*,*) ' '
       write(*,*) '***** Completed all the', amplicol_probe_points,' direct AmpliCol probe points. Stop the probe.'
       amplicol_probe_done=.true.
    endif
  end subroutine perform_amplicol_probe


end module mg_checks
