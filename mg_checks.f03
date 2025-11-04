module mg_checks
  use common
  integer :: k,nord
  real(kind=4),dimension(:),allocatable :: mg_check
  real(kind=4),dimension(100) :: me_code
  logical :: match
  real(kind=8),dimension(:,:),allocatable :: p_read
  real(kind=8),parameter  :: alpha_check=0.118d0

contains

  subroutine get_madgraph_results(n,ichan,iint,me,nlines)
   implicit none
   integer :: i,k,n
   integer,intent(in) :: ichan,iint
   real(kind=4),dimension(:),allocatable :: me
   real(kind=8),dimension(0:3) :: dum
   character(len=50) :: filename
   integer :: nlines
   write(filename, '("momenta_", I0, "_",I0,".txt")') ichan,iint
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
   write(filename, '("momenta_", I0, "_",I0,".txt")') igroup,iamp
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

   write(command, '(A,I0)') "./myscript.sh ", igroup
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



end module mg_checks
