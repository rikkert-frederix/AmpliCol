module mg_checks
  use common
  use amplitude_QCD_mod
  use argument_parser
  use handling_processes
  integer :: k,nord
  real(kind=4),dimension(:),allocatable :: mg_check
  real(kind=4),dimension(100) :: me_code
  logical :: match
  real(kind=8),dimension(:,:),allocatable :: p_read
  real(kind=8),parameter  :: alpha_check=0.118d0
  real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
  integer :: me_points

contains

  subroutine get_madgraph_results(n,ichan,iint,me,nlines)
   implicit none
   integer :: i,k,n
   integer,intent(in) :: ichan,iint
   real(kind=4),dimension(:),allocatable :: me
   real(kind=8),dimension(0:3) :: dum
   character(len=50) :: filename
   integer :: nlines
   write(filename, '("Utilities/ME_checks/momenta_", I0, "_",I0,".txt")') ichan,iint
   open(20,file=trim(filename),status="old")
   do i=1,n
        read(20,*) dum
   enddo
   read(20,*) nlines
   if (allocated(me)) deallocate(me)
   allocate(me(nlines))
   do i=1,nlines
      read(20,*) me(i)
   enddo
   close(20)
  end subroutine get_madgraph_results

  subroutine read_in_momenta(n,igroup,iamp,p_in)
   implicit none
   integer :: i,k,n
   integer,intent(in) :: iamp,igroup
   real(kind=8),dimension(n,0:3),intent(out) :: p_in
   character(len=50) :: filename
   write(filename, '("Utilities/ME_checks/momenta_", I0, "_",I0,".txt")') igroup,iamp
   open(40,file=trim(filename),status="old")
   do i=1,n
        read(40,*) p_in(i,:)
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
        print *, "Bash script failed with status:", status
    else
        print *, "MG script finished successfully fori igroup, iamp: ",igroup,iamp
    end if
  end subroutine run_madgraph_check

  subroutine perform_check(iint,ichan)
    implicit none
    integer :: i,nord,iint,ichan
    if (.not.allocated(mg_check)) allocate(mg_check(1000))
    ! Coupling-order sectors are scaled before squaring.  Applying the old
    ! n_sing-based monomial here would multiply the couplings a second time
    ! and is not meaningful once QCD/EW interference sectors are retained.
    me_code = pgl(ichan)%amp2(:)/dble(pgl(ichan)%iden(iint))
    call get_madgraph_results(pgl(ichan)%next,ichan,iint,mg_check,nord)
    match=.false.
    do i=1,nord
    if (abs((mg_check(i)-me_code(1))/me_code(1)).lt.1d-4) then
            match=.true.
    endif
    enddo
    if (.not.match) then
       write(*,*) 'ERROR: disagreement with MG in ME-level check!',me_code(1)
       write(*,*) mg_check
       !stop 4
    endif
    if (pgl(ichan)%passed(iint).gt.me_points) then
            write(*,*) ' '
            write(*,*) '***** Passed all the', me_points,' ME-level tests. Stop the ME evaluation test.'
    endif
  end subroutine perform_check


end module mg_checks
