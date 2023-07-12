! gfortran -O3 -o test mint_module.f90 test.f90 MC_integer.f ranmar.f HwU.f
! ./test

program TEST
  use mint_module
  implicit none
  integer :: j 
  
! relevant input parameters for integration
  ncalls0=-100    ! Number of events to generate. (If negative, start
                   ! from a small number of points and double it each
                   ! iteration. If positive, this is the number of
                   ! points per iteration as well).

  ndim=1           ! Number of dimensions of the integration.

  itmax=12        ! Number of iterations. (If ncalls0 < 0, the
                   ! integration is aborted if accuracy (next line)
                   ! has been reached.

  accuracy=0.001d0 ! Accuracy of the integration. (Ignored if ncalls0 > 0).

! Not so relevant parameters: only used in special cases.
  fixed_order=.false.
  nlo_ps=.true.
  n_ord_virt=1
  nchans=1
  iconfig=1
  ichan=1
  ifold_energy=1
  ifold_yij=1
  ifold_phi=1
  ifold(1:ndimmax)=1
  iconfigs(1:maxchannels)=1
  min_virt_fraction_mint=1d0
  virt_fraction=1d0
  wgt_mult=1d0
  average_virtual(0:n_ave_virt,maxchannels)=0d0
  virt_wgt_mint(0:n_ave_virt)=0d0
  born_wgt_mint(0:n_ave_virt)=0d0
  virtual_fraction(1:maxchannels)=1d0
  ans(1:nintegrals,0:maxchannels)=0d0
  unc(1:nintegrals,0:maxchannels)=0d0
  only_virt=.false.

! integration steps:
  ! imode=0  (Setting up grids)
  ! imode=-1 (same as imode=0, but starting from existing grids)
  ! imode=1  (computing bounding envelope)
  ! imode=2  (event generation)
  write (*,*) 'Give imode:'
  read (*,*) imode
  
  if (imode.le.1) then
     call mint(integrand)
  else
     call read_grids_from_file
     call gen(integrand,0,-1) ! initialise counters
     do j=1,abs(ncalls0)
        call gen(integrand,1,2) ! generate an unweighted event
        call write_event(j)
     enddo
     call gen(integrand,3,-1) ! print counters
  endif
     
contains
  function integrand(x,vol,ifirst,f1)
    implicit none
    double precision integrand
    integer :: ifirst
    double precision :: vol
    double precision, dimension(ndimmax) :: x
    double precision, dimension(nintegrals) :: f1
    double precision,save :: val
    ! some point-by-point initialisation
    new_point=.true.
    pass_cuts_check=.true.
    f1(1:nintegrals)=0d0

    ! function to integrate is sin(x)
    if (ifirst .eq.0) then
       val=sin(x(1))*vol
    elseif (ifirst.eq.2) then
       continue ! use previously computed integrand
    endif

    ! pass the result to the mint module
    f1(1)=abs(val)
    f1(2)=val
    
  end function integrand

  
  subroutine write_event(j)
    implicit none
    integer :: j
    write (*,*) 'found event',j
  end subroutine write_event

end program TEST

