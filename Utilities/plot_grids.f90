program plot_grids
  implicit none
  character*150 :: filename
  integer,parameter :: nintervals=32,maxchannels=20,ndimmax=60
  double precision, dimension(0:nintervals,ndimmax) :: xgrid
  integer :: next,ndim,i,j
  
  read(*,'(i2,a)') next,filename
  ndim=3*(next-2)-2

  call read_grids_from_file
  
  filename='grids.data'
  open(unit=12,file=filename,status='unknown')
  do i=1,ndim
     do j=0,nintervals
        write (12,*) j,xgrid(j,i)
     enddo
     write (12,*) ''
     write (12,*) ''
     write (12,*) ''
     write (12,*) ''
  enddo
  close(12)
  filename='grids.gnuplot'
  open(unit=12,file=filename,status='unknown')
  write (12,*) 'set lmargin 10'
  write (12,*) 'set rmargin 0'
  write (12,*) 'set terminal postscript portrait enhanced color "Helvetica" 15'
  write (12,*) 'set key font ",15"'
  write (12,*) 'set key samplen "2"'
  write (12,*) 'set output "grids.ps"'
  do i=1,ndim
     write (12,*) 'plot "grids.data" index',i-1
  enddo
  write (12,*) '!ps2pdf "grids.ps" &> /dev/null'
  close(12)
  
  contains
    subroutine read_grids_from_file
! Read the MINT integration grids from file
    implicit none
    integer :: i,j,k,kchan
    character(len=3) :: dummy
    open (unit=12, file=filename,status='old')
    do j=0,nintervals
       read (12,*) dummy,(xgrid(j,i),i=1,ndim)
    enddo
    close (12)
  end subroutine read_grids_from_file


  
end program plot_grids
