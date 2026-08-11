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
    integer :: NUP,IDPRUP
    integer,dimension(pgl%next) :: IDUP,ISTUP
    integer,dimension(2,pgl%next) :: MOTHUP,ICOLUP
    real(kind=8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP
    real(kind=8),dimension(5,pgl%next) :: PUP
    real(kind=8),dimension(pgl%next) :: VTIMUP,SPINUP
    ! determine LHEF info
    NUP=pgl%next
    IDPRUP=1
    XWGTUP=sign(wgt,evt_sign)
    SCALUP=scale_shower
    AQEDUP=alphaEW
    AQCDUP=alphaS
    do i=1,pgl%next
       IDUP(i)=pgl%iden_processes(i,iproc_iden_picked,iproc_picked)
       if(i.le.2) then
          ISTUP(i)=-1
          MOTHUP(1:2,i)=0
       else
          ISTUP(i)=+1
          MOTHUP(1:2,i)=[1,2]
       endif
       PUP(1:3,i)=pgl%ps(1)%p(1:3,i)
       PUP(4,i)=pgl%ps(1)%p(0,i)
       PUP(5,i)=phys_model%get_mass(IDUP(i))
       VTIMUP(i)=0d0
       if (keep_processes_separate) then
          SPINUP(i)=dble(pgl%amps(iproc_picked)%spins(i,hel_picked(1),hel_picked(2)))
       else
          SPINUP(i)=dble(pgl%amps(1)%spins(i,hel_picked(1),hel_picked(2)))
       endif
    enddo
    call get_col_info(pgl,ICOLUP)
    ! write event to file
    write (iunit,'(a)') '<event>'
    write(iunit,503)NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
    do i=1,NUP
       write(iunit,504)IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I), &
            ICOLUP(1,I),ICOLUP(2,I),&
            PUP(1,I),PUP(2,I),PUP(3,I),PUP(4,I),PUP(5,I),&
            VTIMUP(I),SPINUP(I)
    enddo
    write (iunit,506) '#color',pgl%color_orders(1:pgl%next,iproc_picked)
    write (iunit,'(a)') '</event>'
503 format(1x,i2,1x,i6,4(1x,e14.8))
504 format(1x,i8,1x,i2,4(1x,i4),5(1x,e24.17),2(1x,e10.4))
506 format(a,100i3)
  end subroutine write_event

  subroutine event_update_wgt(iunit,ounit,wgt)
    implicit none
    character(len=1024) :: string
    integer :: next,iunit,ounit
    real(kind=8),dimension(3) :: wgt
    logical,save :: firsttime=.true.
    integer :: NUP,IDPRUP
    integer :: LPRUP
    real(kind=8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,XSECUP,XERRUP,XMAXUP
    if (firsttime) then
       do
          read(iunit,'(a)') string
          if (index(string,"<init>").ne.0) then
             ! Update the cross section
             write(ounit,'(a)') trim(string)
             read(iunit,'(a)') string
             write(ounit,'(a)') trim(string)
             read(iunit,'(a)') string
             XSECUP=simple_integrator%res(2)
             XERRUP=simple_integrator%unc(2)
             XMAXUP=simple_integrator%res(1)
             LPRUP=1
             write(ounit,502)XSECUP,XERRUP,XMAXUP,LPRUP
             read(iunit,'(a)') string
          endif
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
    read(iunit,503)NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
    XWGTUP=wgt(1)
    if (wgt(1).ne.0d0) &
         write(ounit,503)NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
    do
       read(iunit,'(a)') string
       if (index(string,"</event>").ne.0) exit
       if (wgt(1).ne.0d0) write(ounit,'(a)') trim(string)
    enddo
    if (wgt(1).ne.0d0) then
       write (ounit,505) '#overwgt',wgt(1:3)
       write (ounit,'(a)') '</event>'
    endif
502 format(3(1x,e14.8),1x,i6)
503 format(1x,i2,1x,i6,4(1x,e14.8))
505 format(a,3(1x,e14.8))
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

  subroutine write_unique_in_file(pgl_unique,unique_map,unique_map_value,nevents)
    implicit none
    type(phase_space_order_group),allocatable :: pgl_unique
    real(kind=8),dimension(pgl_unique%nproc) :: unique_map_value
    integer,dimension(pgl_unique%nproc) :: unique_map
    integer :: iproc,nevents
    integer :: IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP,LPRUP
    real(kind=8) :: EBMUP(2),XSECUP,XERRUP,XMAXUP
    integer(kind=8) iseed
    common /to_seed/iseed
    IDBMUP(1:2)=2212     ! two protons
    EBMUP(1:2)=sqrts/2d0 ! half of collision energy
    PDFGUP(1:2)=-1
    PDFSUP(1:2)=pdf_lhaid
    IDWTUP=-3
    NPRUP=1
    XSECUP=0d0
    XERRUP=0d0
    XMAXUP=0d0
    LPRUP=1
    write(11,'(a)') '<LesHouchesEvents version="3.0">'
    write(11,'(a)') '<header>'
    write(11,*) pgl_unique%next,pgl_unique%nproc
    do iproc=1,pgl_unique%nproc
       write(11,*) unique_map(iproc),unique_map_value(iproc),pgl_unique%processes(1:pgl_unique%next,iproc)
    enddo
    write(11,'(a,1x,i12,1x,a)') '<nevents>',nevents,'</nevents>'
    write(11,'(a,1x,i12,1x,a)') '<seed>   ',iseed,  '</seed>'
    write(11,'(a)') '</header>'
    write(11,'(a)') '<init>'
    write(11,501)IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),PDFGUP(1)&
         &,PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
    write(11,502)XSECUP,XERRUP,XMAXUP,LPRUP
    write(11,'(a)') "<generator name='AmpliCol' version='1.0'>please cite arXiv:2601.19483</generator>"
    write(11,'(a)') '</init>'
501 format(2(1x,i6),2(1x,e14.8),2(1x,i2),2(1x,i8),1x,i2,1x,i3)
502 format(3(1x,e14.8),1x,i6)
  end subroutine write_unique_in_file


  subroutine get_col_info(pgl,ICOLUP)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,dimension(2,pgl%next),intent(out) :: ICOLUP
    integer,dimension(pgl%next) :: ipdg,ord
    integer :: i,label
    ICOLUP=0
    ! Take anti-particles for the initial state
    do i=1,2
       ipdg(i)=phys_model%get_antipart(pgl%iden_processes(i,iproc_iden_picked,iproc_picked))
    enddo
    ipdg(3:pgl%next)=pgl%iden_processes(3:pgl%next,iproc_iden_picked,iproc_picked)
    ord(1:pgl%next)=pgl%color_orders(1:pgl%next,iproc_picked)
    i=1
    label=501
    do
       if (phys_model%is_quark(iPDG(ord(i))) .or. phys_model%is_gluon(iPDG(ord(i)))) then
          if (ICOLUP(1,ord(i)).ne.0) exit ! already filled. We have made the full loop and are done.
          ! found a colour line.
          ICOLUP(1,ord(i))=label
          ! find the corresponding anti-colour line. This should be
          ! the next particle in the colour ordering (as long as
          ! that's not a colour singlet.
          do
             i=mod(i,pgl%next)+1
             if (phys_model%is_antiquark(iPDG(ord(i))) .or. phys_model%is_gluon(iPDG(ord(i)))) then
                ICOLUP(2,ord(i))=label
                label=label+1
                exit
             elseif (phys_model%is_quark(iPDG(ord(i)))) then
                write (*,*) 'Cannot have a quark after another one'
                stop 1
             endif
          enddo
       else
          i=mod(i,pgl%next)+1
       endif
    enddo
    ! Flip initial states
    ICOLUP(1:2,1:2)=ICOLUP(2:1:-1,1:2)
  end subroutine get_col_info
end module handling_events
