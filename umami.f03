module umami
  implicit none
  private
  type channel_group
     integer :: n,nvals,n_sing
     integer,dimension(:),allocatable :: namps,iproc_start
     real(kind=8),dimension(:),allocatable :: col_fac,iden
     real(kind=8),dimension(:,:),allocatable :: hel_fac
  end type channel_group
  integer :: nchans
  type(channel_group),dimension(:),allocatable :: chans
  logical :: keep_processes_separate
  real(kind=8) :: alphaEW
  public :: umami_initialize,umami_matrix_element,umami_free
contains

  subroutine umami_initialize()
    implicit none
    call read_amplitude_lib_light()
  end subroutine umami_initialize

  subroutine umami_free()
    implicit none
    integer :: ichan
    do ichan=1,nchans
       if (allocated(chans(ichan)%namps)) deallocate(chans(ichan)%namps)
       if (allocated(chans(ichan)%iproc_start)) deallocate(chans(ichan)%iproc_start)
       if (allocated(chans(ichan)%col_fac)) deallocate(chans(ichan)%col_fac)
       if (allocated(chans(ichan)%iden)) deallocate(chans(ichan)%iden)
       if (allocated(chans(ichan)%hel_fac)) deallocate(chans(ichan)%hel_fac)
    enddo
    if (allocated(chans)) deallocate(chans)
  end subroutine umami_free
  
  subroutine umami_matrix_element(ichan,iint,p,alphas,rnd,amp2,hels)
    use amp_lib
    implicit none
    integer,intent(in) :: ichan,iint
    real(kind=8),dimension(0:3,chans(ichan)%n),intent(in) :: p
    real(kind=8),intent(in) :: alphas
    real(kind=8),dimension(chans(ichan)%nvals),intent(in) :: rnd
    real(kind=8),dimension(chans(ichan)%nvals),intent(out) :: amp2
    integer,dimension(chans(ichan)%nvals),intent(out) :: hels

    integer :: iproc,ih
    complex(kind=8),dimension(chans(ichan)%namps(iint)) :: amps
    real(kind=8),dimension(chans(ichan)%namps(iint)) :: amp2_hel
    real(kind=8) :: coupling
    real(kind=8), parameter :: pi=3.14159265358979323846d0
    
    call evaluate_amp(ichan,iint,p,amps)

    if (keep_processes_separate) then
       do ih=1,chans(ichan)%namps(iint)
          amp2_hel(ih)=dble(amps(ih)*chans(ichan)%col_fac(iint)*dconjg(amps(ih)))*chans(ichan)%hel_fac(ih,iint)
       enddo
       amp2(1)=sum(amp2_hel)
    else
       amp2=0d0
       iproc=0
       do ih=1,chans(ichan)%namps(1)
          do while (chans(ichan)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
          amp2_hel(ih)=dble(amps(ih)*chans(ichan)%col_fac(iproc)*dconjg(amps(ih)))*chans(ichan)%hel_fac(ih,1)
          amp2(iproc)=amp2(iproc)+amp2_hel(ih)
       enddo
    endif
    
    ! Pick a random helicity for each process
    call unwgt_helicities(ichan,iint,amp2_hel,rnd,hels)

    coupling=1d0
    ! The strong coupling
    if (chans(ichan)%n_sing.lt.chans(ichan)%n-2) then
       coupling=coupling*(4*pi*alphas)**(chans(ichan)%n-2-chans(ichan)%n_sing)
    endif
    ! The EW coupling
    if (chans(ichan)%n_sing.ge.1) then
       coupling=coupling*(2d0*4d0*pi*alphaEW)**chans(ichan)%n_sing
    endif

    ! multiply the matrix elements by the QCD and EW couplings
    amp2=amp2*coupling

    ! Include identical particle symmetry factors
    if (keep_processes_separate) then
       amp2(1)=amp2(1)/chans(ichan)%iden(iint)
    else
       amp2(1:chans(ichan)%nvals)=amp2(1:chans(ichan)%nvals)/chans(ichan)%iden(1:chans(ichan)%nvals)
    endif
    
  end subroutine umami_matrix_element

  subroutine unwgt_helicities(ichan,iint,amp2_hel,rnd,hels)
    implicit none
    integer,intent(in) :: ichan,iint
    real(kind=8),dimension((chans(ichan)%namps(iint))),intent(in) :: amp2_hel
    real(kind=8),dimension(chans(ichan)%nvals),intent(in) :: rnd
    integer,dimension((chans(ichan)%namps(iint))),intent(out) :: hels
    integer :: iproc,isize
    if (keep_processes_separate) then
       call pick_one(chans(ichan)%namps(iint),amp2_hel(1),rnd(1),hels(1))
    else
       do iproc=1,chans(ichan)%nvals
          isize=chans(ichan)%iproc_start(iproc+1)-chans(ichan)%iproc_start(iproc)
          call pick_one(isize,amp2_hel(chans(ichan)%iproc_start(iproc)),rnd(iproc),hels(iproc))
       enddo
    end if
  end subroutine unwgt_helicities

  subroutine pick_one(isize,vals,rnd,ipicked)
    implicit none
    integer,intent(in) :: isize
    real(kind=8),dimension(isize),intent(in) :: vals
    real(kind=8),intent(in) :: rnd
    integer,intent(out) :: ipicked
    real(kind=8) :: accum,target
    target=rnd*sum(vals)
    ipicked=1
    accum=vals(1)
    do
       if (accum.gt.target) then
          exit
       else
          ipicked=ipicked+1
          accum=accum+vals(ipicked)
       endif
    enddo
    if (ipicked.gt.isize) then
       write (*,*) 'Inconsistency in "pick_one" subroutine',isize,ipicked
       stop 1
    endif
  end subroutine pick_one

  subroutine read_amplitude_lib_light()
    implicit none
    character(len=170) :: filename
    integer :: idim,ichan,max_hel,max_iint
    integer,dimension(:),allocatable :: iproc_start
    filename='Library/amplitudes_light.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='old')
!!$    open(unit=14,file=filename,status='old')
    read(14)keep_processes_separate
    read(14)nchans
    read(14)alphaEW
    allocate(chans(nchans))
    do ichan=1,nchans
       read(14) chans(ichan)%n
       read(14) chans(ichan)%nvals
       read(14) chans(ichan)%n_sing
       read(14) max_iint
       read(14) max_hel

       if (chans(ichan)%nvals.eq.1) then
          allocate(chans(ichan)%namps(max_iint))
       else
          allocate(chans(ichan)%namps(1))
       endif
       read(14) chans(ichan)%namps

       if (.not.keep_processes_separate) then
          allocate(chans(ichan)%iproc_start(chans(ichan)%nvals+1))
          read(14) idim,chans(ichan)%iproc_start
          if (idim.ne.chans(ichan)%nvals+1) then
             write (*,*) 'iproc_start has wrong size',ichan,idim,chans(ichan)%nvals+1
             stop 1
          endif
       endif
       
       allocate(chans(ichan)%col_fac(max_iint))
       read(14) idim,chans(ichan)%col_fac
       if (idim.ne.max_iint) then
          write (*,*) 'col_fac has wrong size',ichan,idim,max_iint
          stop 1
       endif
       
       allocate(chans(ichan)%iden(max_iint))
       read(14) idim,chans(ichan)%iden
       if (idim.ne.max_iint) then
          write (*,*) 'iden has wrong size',ichan,idim,max_iint
          stop 1
       endif

       if (keep_processes_separate) then
          allocate(chans(ichan)%hel_fac(max_hel,max_iint))
       else
          allocate(chans(ichan)%hel_fac(max_hel,1))
       endif
       read(14) idim,chans(ichan)%hel_fac

    enddo
    close(14)
  end subroutine read_amplitude_lib_light

end module umami
