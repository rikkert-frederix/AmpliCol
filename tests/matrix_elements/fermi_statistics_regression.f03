program fermi_statistics_regression
  use amplitude_QCD_mod
  use particles
  implicit none
  type(physics_model) :: model
  integer,parameter :: dp=kind(1d0)
  ! Stock MadGraph 3.2.0 full-colour SMATRIX ratios for d d~ -> 4l and
  ! d d~ -> 4l+2g at the points below.  MadGraph includes the 1/(2! 2!)
  ! identical-lepton factor, so multiply that result by four to compare
  ! labelled phase-space points.  The 4l+2g result has two colour orders.
  real(kind=dp),parameter :: madgraph_fc_ratio_4l=115.97215136285963d0
  real(kind=dp),parameter :: madgraph_fc_ratio_4l2g=1.1834035986208604d0
  integer,dimension(6) :: process_distinct,process_identical
  integer,dimension(8) :: process_distinct_gg,process_identical_gg
  real(kind=dp),dimension(0:3,6) :: p4l
  real(kind=dp),dimension(0:3,8) :: p4l2g
  real(kind=dp),dimension(3) :: col_distinct,col_identical,col_distinct_gg,col_identical_gg
  integer :: ncolour_orders
  character(len=64) :: option

  open(unit=99,file='/dev/null',status='unknown',action='write')
  call model%init_part(173d0,0d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call model%init_vert()
  if (command_argument_count().gt.0) then
     call get_command_argument(1,option)
     if (trim(option).eq.'--emit-library') then
        call emit_identical_library()
        stop
     endif
     write (*,*) 'Unknown option: ',trim(option)
     stop 1
  endif

  process_distinct=[1,-1,-11,11,-13,13]
  process_identical=[1,-1,-11,11,-11,11]
  process_distinct_gg=[1,-1,-11,11,-13,13,21,21]
  process_identical_gg=[1,-1,-11,11,-11,11,21,21]
  call fill_four_lepton_momenta(p4l)
  call fill_six_body_momenta(p4l2g)

  call check_identical_exchange(process_identical,p4l,3,5,'four-lepton positron exchange')
  call check_six_lepton_exchange()

  call colour_squared(process_distinct,p4l,col_distinct)
  call colour_squared(process_identical,p4l,col_identical)
  write (*,'(a,es24.16)') 'FC_DISTINCT_4L=',col_distinct(3)
  write (*,'(a,es24.16)') 'FC_IDENTICAL_4L=',col_identical(3)
  write (*,'(a,es24.16)') 'FC_IDENTICAL_TO_DISTINCT_4L=',col_identical(3)/col_distinct(3)
  if (abs(col_identical(3)/col_distinct(3)/madgraph_fc_ratio_4l-1d0).gt.1d-7) then
     write (*,*) 'Full-colour four-lepton result disagrees with MadGraph:',&
          col_identical(3)/col_distinct(3),madgraph_fc_ratio_4l
     stop 1
  endif

  call colour_squared(process_distinct_gg,p4l2g,col_distinct_gg,ncolour_orders)
  call colour_squared(process_identical_gg,p4l2g,col_identical_gg)
  if (ncolour_orders.le.1) then
     write (*,*) 'The 4l+2g regression does not have multiple colour orders'
     stop 1
  endif
  write (*,'(a,3es24.16)') 'LC_NLC_FC_DISTINCT_4L2G=',col_distinct_gg
  write (*,'(a,3es24.16)') 'LC_NLC_FC_IDENTICAL_4L2G=',col_identical_gg
  write (*,'(a,es24.16)') 'FC_IDENTICAL_TO_DISTINCT_4L2G=',&
       col_identical_gg(3)/col_distinct_gg(3)
  if (abs(col_identical_gg(3)/col_distinct_gg(3)/madgraph_fc_ratio_4l2g-1d0).gt.1d-7) then
     write (*,*) 'Full-colour 4l+2g result disagrees with MadGraph:',&
          col_identical_gg(3)/col_distinct_gg(3),madgraph_fc_ratio_4l2g
     stop 1
  endif
  if (abs(col_distinct_gg(3)/col_distinct_gg(1)-1d0).lt.0.05d0) then
     write (*,*) 'The full-colour 4l+2g check is too close to leading colour'
     stop 1
  endif
  write (*,'(a)') 'Fermi-statistics regression passed'

contains

  subroutine check_identical_exchange(process,p,leg1,leg2,label)
    implicit none
    integer,dimension(:),intent(in) :: process
    real(kind=dp),dimension(0:,:),intent(in) :: p
    integer,intent(in) :: leg1,leg2
    character(len=*),intent(in) :: label
    type(amplitude_QCD) :: amp
    integer :: n,i,j,match
    integer,dimension(:,:),allocatable :: part,orders,spin,spins_before
    integer,dimension(:),allocatable :: hel,target_spin
    real(kind=dp),dimension(:,:),allocatable :: p_swap
    complex(kind=dp),dimension(:),allocatable :: amps_before,amps_after
    real(kind=dp) :: max_amp,max_residual

    n=size(process)
    allocate(part(n,1),orders(n,1),spin(0:3,n),hel(n))
    allocate(target_spin(n),p_swap(0:3,n))
    part(:,1)=process
    orders(1,1)=2
    do i=3,n
       orders(i-1,1)=i
    enddo
    orders(n,1)=1
    spin=0
    spin(0,:)=2
    spin(1,:)=-1
    spin(2,:)=1
    hel=0
    call amp%init(1,n,1,part,spin,orders,model)
    call amp%evaluate(n,p,hel,.false.,model)
    allocate(amps_before(amp%n_amps),amps_after(amp%n_amps))
    allocate(spins_before(n,amp%n_amps))
    amps_before=amp%amps
    spins_before=amp%spins(:,1,:)

    p_swap=p
    p_swap(:,leg1)=p(:,leg2)
    p_swap(:,leg2)=p(:,leg1)
    call amp%evaluate(n,p_swap,hel,.false.,model)
    amps_after=amp%amps

    max_amp=max(maxval(abs(amps_before)),maxval(abs(amps_after)))
    max_residual=0d0
    do i=1,amp%n_amps
       target_spin=spins_before(:,i)
       target_spin(leg1)=spins_before(leg2,i)
       target_spin(leg2)=spins_before(leg1,i)
       match=0
       do j=1,amp%n_amps
          if (all(spins_before(:,j).eq.target_spin)) then
             match=j
             exit
          endif
       enddo
       if (match.eq.0) then
          write (*,*) 'Missing exchanged helicity in ',trim(label),target_spin
          stop 1
       endif
       max_residual=max(max_residual,abs(amps_before(i)+amps_after(match)))
    enddo
    if (max_residual.gt.1d-10*max(1d-30,max_amp)) then
       write (*,*) 'Identical-lepton antisymmetry failure in ',trim(label)
       write (*,*) 'maximum amplitude/residual:',max_amp,max_residual
       stop 1
    endif
    write (*,'(a,a,a,es12.4)') 'Passed ',trim(label),'; residual=',max_residual
  end subroutine check_identical_exchange

  subroutine check_six_lepton_exchange()
    implicit none
    integer,dimension(8) :: process
    real(kind=dp),dimension(0:3,8) :: p
    process=[1,-1,-11,-11,-11,11,11,11]
    call fill_six_body_momenta(p)
    call check_identical_exchange(process,p,3,4,'six-lepton positron exchange')
  end subroutine check_six_lepton_exchange

  subroutine emit_identical_library()
    implicit none
    type(amplitude_QCD) :: amp
    integer,parameter :: n=6
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    real(kind=dp),dimension(0:3,n) :: p

    part(:,1)=[1,-1,-11,11,-11,11]
    orders(:,1)=[2,3,4,5,6,1]
    hel=[-1,1,1,-1,1,-1]
    spin=0
    spin(0,:)=1
    spin(1,:)=hel
    call fill_four_lepton_momenta(p)
    call amp%init(1,n,1,part,spin,orders,model)
    call amp%evaluate(n,p,hel,.false.,model)
    if (maxval(abs(amp%amps)).lt.1d-30) then
       write (*,*) 'Cannot create regression library from a zero amplitude'
       stop 1
    endif
    call amp%create_library(n,hel,1,1,model,p)
    write (*,'(a)') 'Generated identical-lepton regression library'
  end subroutine emit_identical_library

  subroutine colour_squared(process,p,matrix2,ncolour_orders)
    implicit none
    integer,dimension(:),intent(in) :: process
    real(kind=dp),dimension(0:,:),intent(in) :: p
    real(kind=dp),dimension(3),intent(out) :: matrix2
    integer,intent(out),optional :: ncolour_orders
    type(amplitude_QCD) :: amp
    integer :: n,i,ih,iacc,irow,ival,ic,icol,ioff
    integer,dimension(:,:),allocatable :: part,orders,spin
    integer,dimension(:),allocatable :: hel
    complex(kind=dp) :: amp2_colour,amp_colour

    n=size(process)
    allocate(part(n,1),orders(n,1),spin(0:3,n),hel(n))
    part(:,1)=process
    do i=1,n
       orders(i,1)=i
    enddo
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    call amp%init(2,n,1,part,spin,orders,model)
    call amp%init_col(n,20)
    if (present(ncolour_orders)) ncolour_orders=amp%nColOrd
    matrix2=0d0
    ioff=amp%iproc_start(amp%nprocs)-1
    do ih=0,2**n-1
       do i=1,n
          if (btest(ih,i-1)) then
             hel(i)=1
          else
             hel(i)=-1
          endif
       enddo
       call amp%evaluate(n,p,hel,.false.,model)
       do iacc=1,3
          do irow=1,amp%nColOrd
             amp_colour=(0d0,0d0)
             do ival=1,amp%n_col_vals(iacc)
                amp2_colour=(0d0,0d0)
                do ic=amp%row_index(irow-1,ival,iacc)+1,amp%row_index(irow,ival,iacc)
                   icol=amp%col_index(amp%i_col_i(ival,iacc)+ic)
                   amp2_colour=amp2_colour+amp%amps(ioff+icol)
                enddo
                amp_colour=amp_colour+amp2_colour*amp%diff_col_vals(ival,iacc)
             enddo
             matrix2(iacc)=matrix2(iacc)+dble(amp_colour*conjg(amp%amps(ioff+irow)))
          enddo
       enddo
    enddo
  end subroutine colour_squared

  subroutine fill_four_lepton_momenta(p)
    implicit none
    real(kind=dp),dimension(0:3,6),intent(out) :: p
    p(:,1)=[500.0000000d0,   0.00000000d0,  0.0000000000d0,   500.0000000d0]
    p(:,2)=[500.0000000d0,   0.00000000d0,  0.0000000000d0,  -500.0000000d0]
    p(:,3)=[88.55133305d0,  -22.1006902d0,  40.080353191d0,  -75.80543095d0]
    p(:,4)=[328.3294192d0,  -103.849611d0, -301.93375538d0,   76.49492138d0]
    p(:,5)=[152.3581094d0,  -105.880959d0, -97.709638326d0,   49.54838522d0]
    p(:,6)=[430.7611382d0,   231.831261d0,  359.56304052d0,  -50.23787565d0]
  end subroutine fill_four_lepton_momenta

  subroutine fill_six_body_momenta(p)
    implicit none
    real(kind=dp),dimension(0:3,8),intent(out) :: p
    p=0d0
    p(:,1)=[468d0,0d0,0d0,468d0]
    p(:,2)=[468d0,0d0,0d0,-468d0]
    p(:,3)=[133d0,39.9d0,53.2d0,133d0*sqrt(0.75d0)]
    p(:,4)=[133d0,-39.9d0,-53.2d0,-133d0*sqrt(0.75d0)]
    p(:,5)=[156d0,93.6d0,-31.2d0,156d0*sqrt(0.60d0)]
    p(:,6)=[156d0,-93.6d0,31.2d0,-156d0*sqrt(0.60d0)]
    p(:,7)=[179d0,-17.9d0,125.3d0,179d0*sqrt(0.50d0)]
    p(:,8)=[179d0,17.9d0,-125.3d0,-179d0*sqrt(0.50d0)]
  end subroutine fill_six_body_momenta

end program fermi_statistics_regression
