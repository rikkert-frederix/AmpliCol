module handling_events
  use common
  use handling_processes
  integer :: iproc_picked,iproc_iden_picked
  integer,dimension(2) :: hel_picked
  real(kind=8) :: evt_sign
contains
  subroutine write_event(iunit,pgl,wgt)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer :: i,iunit
    real(kind=8) :: wgt
    real(kind=8),external :: ran2
    write (iunit,*) '<event>'
    write (iunit,'(i4,e18.10)') pgl%next,sign(wgt,evt_sign)
    if (keep_processes_separate) then
       write (iunit,'(100i3)') pgl%amps(iproc_picked)%spins(1:pgl%next,hel_picked(1),hel_picked(2))
    else
       write (iunit,'(100i3)') pgl%amps(1)%spins(1:pgl%next,hel_picked(1),hel_picked(2))
    endif
!!$    if (.not.read_proc_from_file) iproc_picked=1
    write (iunit,'(100i3)') pgl%color_orders(1:pgl%next,iproc_picked)
    ! Since some of the symmetry factors (in particular for gg->qqbar+ng)
    ! compensate for reducing the number of integration channels assuming
    ! symmetric initial states, we need to randomly flip all z-components in
    ! those cases. Easiest to always do this if the two incoming particles are
    ! identical.
    if (pgl%processes(1,iproc_picked).ne.pgl%processes(2,iproc_picked) .or. ran2().lt.0.5d0) then
       ! do not flip
       do i=1,pgl%next
          write (iunit,*) pgl%iden_processes(i,iproc_iden_picked,iproc_picked),&
               pgl%phase_space%p(1:3,i),pgl%phase_space%p(0,i)
       enddo
    else
       ! do flip
       do i=1,pgl%next
          if (i.le.2) then
             write (iunit,*) pgl%iden_processes(i,iproc_iden_picked,iproc_picked),&
                  pgl%phase_space%p(1:2,3-i),-pgl%phase_space%p(3,3-i),pgl%phase_space%p(0,3-i)
          else
             write (iunit,*) pgl%iden_processes(i,iproc_iden_picked,iproc_picked),&
                  pgl%phase_space%p(1:2,i),-pgl%phase_space%p(3,i),pgl%phase_space%p(0,i)
          endif
       enddo
    endif
    write (iunit,*) '</event>'
  end subroutine write_event

  subroutine event_update_wgt(iunit,ounit,wgt)
    implicit none
    character(len=1024) :: string
    integer :: next,iunit,ounit
    real(kind=8),dimension(3) :: wgt
    logical,save :: firsttime=.true.
    if (firsttime) then
       do
          read(iunit,'(a)') string
          if (index(string,"<event>").ne.0) then
             if (wgt(1).ne.0d0) write(ounit,'(a)') trim(string)
             exit
          endif
          write(ounit,'(a)') trim(string)
       enddo
       firsttime=.false.
    else
       read(iunit,'(a)') string
       if (wgt(1).ne.0d0) write(ounit,'(a)') trim(string)
    endif
    read(iunit,*) next
    if (wgt(1).ne.0d0) write(ounit,'(i4,1x,e18.10,1x,e18.10,1x,e18.10)') next,wgt(1:3)
    do
       read(iunit,'(a)') string
       if (wgt(1).ne.0d0) write(ounit,'(a)') trim(string)
       if (index(string,"</event>").ne.0) exit
    enddo
  end subroutine event_update_wgt
  
  subroutine unwgt_process(pgl,iint)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: iint
    integer :: i,iproc
    real(kind=8) :: random,accum,target
    real(kind=8),external :: ran2
    target=0d0
    do iproc=1,pgl%nproc
       if (keep_processes_separate .and. iproc.ne.iint) cycle
       do i=1,pgl%iden_iproc(iproc)
          target=target+abs(pgl%val_procs(i,iproc))
       enddo
    enddo
    random=ran2()*target
    if (keep_processes_separate) then
       iproc=iint
    else
       iproc=1
    endif
    i=1
    accum=abs(pgl%val_procs(i,iproc))
    do
       if (accum.gt.random) then
          exit
       else
          i=i+1
          if (i.gt.pgl%iden_iproc(iproc)) then
             i=1
             iproc=iproc+1
          endif
          accum=accum+abs(pgl%val_procs(i,iproc))
       endif
    enddo
    iproc_picked=iproc
    iproc_iden_picked=i
    if (iproc_iden_picked.gt.pgl%iden_iproc(iproc)) then
       write (*,*) "Could not unweight process",iproc,iproc_iden_picked,pgl%iden_iproc(iproc)
       stop 1
    endif
    if (pgl%val_procs(iproc_iden_picked,iproc_picked).lt.0d0) then
       evt_sign=-1d0
    else
       evt_sign=+1d0
    endif
    if (keep_processes_separate .and. iproc_picked.ne.iint) then
       write (*,*) 'Could not unweight process correctly (keep_processes_separate=true)'
       stop 1
    endif
  end subroutine unwgt_process

  subroutine unwgt_helicity(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i
    real(kind=8) :: random
    real(kind=8),external :: ran2
    if (keep_processes_separate) then
       random=ran2()*pgl%amp2(1)
       i=pgl%amps(iproc_picked)%iproc_start(1)
    else
       random=ran2()*pgl%amp2(iproc_picked)
       i=pgl%amps(1)%iproc_start(iproc_picked)
    endif
    do
       if (pgl%amp2_hel(i).gt.random) then
          exit
       else
          i=i+1
          pgl%amp2_hel(i)=pgl%amp2_hel(i)+pgl%amp2_hel(i-1)
       endif
    enddo
    hel_picked(2)=i
    if (keep_processes_separate) then
       if ( hel_picked(2).lt.pgl%amps(iproc_picked)%iproc_start(1) .or. &
            hel_picked(2).ge.pgl%amps(iproc_picked)%iproc_start(1+1)) then
          write (*,*) 'Could not unweight helicity',hel_picked,iproc_picked, &
               pgl%amps(iproc_picked)%iproc_start(1),pgl%amps(iproc_picked)%iproc_start(1+1)
          stop 1
       endif
       if (pgl%hel_fac(hel_picked(2),iproc_picked).gt.1) then
          hel_picked(1)=1+int(ran2()*pgl%hel_fac(hel_picked(2),iproc_picked))
       else
          hel_picked(1)=1
       endif
    else
       if ( hel_picked(2).lt.pgl%amps(1)%iproc_start(iproc_picked) .or. &
            hel_picked(2).ge.pgl%amps(1)%iproc_start(iproc_picked+1)) then
          write (*,*) 'Could not unweight helicity',hel_picked,iproc_picked,pgl%amps(1)%iproc_start(iproc_picked),&
               pgl%amps(1)%iproc_start(iproc_picked+1)
          stop 1
       endif
       if (pgl%hel_fac(hel_picked(2),1).gt.1) then
          hel_picked(1)=1+int(ran2()*pgl%hel_fac(hel_picked(2),1))
       else
          hel_picked(1)=1
       endif
    endif
  end subroutine unwgt_helicity

  subroutine write_unique_in_file(pgl_unique,unique_map,unique_map_value)
    implicit none
    type(phase_space_order_group),allocatable :: pgl_unique
    real(kind=8),dimension(pgl_unique%nproc) :: unique_map_value
    integer,dimension(pgl_unique%nproc) :: unique_map
    integer :: iproc
    write(11,*) pgl_unique%next,pgl_unique%nproc
    do iproc=1,pgl_unique%nproc
       write(11,*) unique_map(iproc),unique_map_value(iproc),pgl_unique%processes(1:pgl_unique%next,iproc)
    enddo
  end subroutine write_unique_in_file
end module handling_events
