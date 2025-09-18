program check_wgts

  integer :: ifile,nPSpoints,i,nover
  logical :: done
  character(len=256) :: argv,filename
  integer,parameter :: max_PSpoints=10000000
  real(kind=8),dimension(2,max_PSpoints) :: wgts
  real(kind=8),dimension(max_PSpoints) :: excess
  real(kind=8) :: evt_wgt,evt_overwgt,evt_ex

  CALL GET_COMMAND_ARGUMENT(1, argv)
  read(argv,'(a)') filename
  ifile=11
  open(unit=ifile,file=filename,status='OLD')

  nPSpoints=0
  do
     call read_event(ifile,done)
     if (done) exit
     if (nPSpoints.eq.max_PSpoints) then
        write (*,*) 'WARNING: more than',max_PSpoints, &
             'events in event file. Taking only first',max_PSpoints,'events'
        exit
     endif
     nPSpoints=nPSpoints+1
     wgts(1,nPSpoints)=evt_wgt
     wgts(2,nPSpoints)=evt_overwgt
     excess(nPSpoints)=evt_ex
  enddo
  write (*,*) 'sum of weights, number of weights, average weight                :', &
       sum(wgts(1,1:nPSpoints)),nPSpoints,sum(wgts(1,1:nPSpoints))/dble(nPSpoints)
  write (*,*) 'sum of weights, number of weights, average weight (from overwgts):', &
       sum(wgts(2,1:nPSpoints)),nPSpoints,sum(wgts(2,1:nPSpoints))/dble(nPSpoints)
  write (14,*) '    event ,   overweight'
  nover=0
  do i=1,nPSpoints
     if(excess(i).ne.0d0) then
        nover=nover+1
        write (14,'(i10,1x,a,1x,e18.10)') i,',',excess(i)
     endif
  enddo
  write (*,'(a,1x,i6,1x,e10.4,1x,a,1x,e10.4,1x,e10.4)') &
       'number of overweight (fraction), cross section overweight (fraction)',&
       nover,dble(nover)/dble(nPSpoints),',',&
       sum(excess(1:nPSpoints))*sum(wgts(1,1:nPSpoints))/dble(nPSpoints)**2,sum(excess(1:nPSpoints))/dble(nPSpoints)
contains
  subroutine read_event(iunit,done)
    implicit none
    integer :: i,iunit,next
    logical :: done
    character(len=100) :: dummy
    character(len=15) :: dummy2
    real(kind=8) :: dum,wgt,amp2,weight
    done=.false.
    do
       read(iunit,'(a)',end=99,err=99) dummy
       if (index(dummy,"<event>").ne.0) exit
    enddo
    read (iunit,*,err=99,end=99) next,evt_wgt,evt_overwgt,evt_ex!,wgt,amp2,weight
    read (iunit,*) dum ! helicity
    read (iunit,*) dum ! color
    do i=1,next
       read (iunit,*,err=99,end=99) dum
    enddo
    read (iunit,*,err=99,end=99) dummy
    return
99  done=.true.
  end subroutine read_event

  recursive subroutine quicksort(a, left, right)
    real(kind=8), intent(inout) :: a(:)
    integer, intent(in) :: left, right
    integer :: i, j
    real(kind=8) :: pivot, tmp
    if (left >= right) return
    pivot = a((left+right)/2)
    i = left
    j = right
    do
      do while (a(i) < pivot)
        i = i + 1
      end do
      do while (a(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = a(i); a(i) = a(j); a(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call quicksort(a, left, j)
    if (i < right) call quicksort(a, i, right)
  end subroutine quicksort

  subroutine sort_doubles(isize,a)
    real(kind=8), intent(inout) :: a(isize)
    call quicksort(a, 1, isize)
  end subroutine sort_doubles
  
end program
