program resonance_phase_space_regression
  use phase_space_base
  use phase_space_gen23_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use, intrinsic :: ieee_exceptions, only: ieee_set_flag,ieee_all
  implicit none
  character(len=256) :: log_line
  integer :: io_status

  open(unit=99,status='scratch',action='readwrite')
  call check_breit_wigner_quantile()
  call check_massive_quark_external_legs()
  call check_root_breit_wigner_with_pdfs()
  call check_nested_and_disjoint_inverse()
  call check_alternative_maps_across_scales()
  call check_singular_boundary_and_nonfinite_rejection()
  call ieee_set_flag(ieee_all,.false.)
  rewind(99)
  do
     read(99,'(a)',iostat=io_status) log_line
     if (io_status.ne.0) exit
     if (index(log_line,'LUP decomposition').gt.0 .or. &
          index(log_line,'Warning: gram4').gt.0) then
        write (*,*) 'Inverse-map failure produced determinant warning spam:',&
             trim(log_line)
        stop 1
     endif
  enddo
  close(99)
  write (*,*) 'resonance phase-space regression: PASS'

contains

  subroutine check_breit_wigner_quantile()
    type(phase_space_gen23) :: phase_space
    type(psv) :: point,inverse_point
    integer,parameter :: n=5
    integer,dimension(n) :: order
    integer,dimension(1) :: resonance_pdgs,resonance_masks
    real(kind=8),dimension(1) :: resonance_masses,resonance_widths
    real(kind=8),dimension(n) :: masses,pt_cut,rap_cut
    real(kind=8),dimension(n,n) :: dr_cut,sqrt_s_min
    real(kind=8),dimension(5) :: original_x
    real(kind=8),dimension(0:3) :: pair_momentum
    real(kind=8) :: pair_mass2,expected_mass2,theta_min,theta_max,theta
    real(kind=8) :: forward_jac

    order=[1,2,3,4,5]
    masses=0d0
    pt_cut=0d0
    rap_cut=0d0
    dr_cut=0d0
    sqrt_s_min=0d0
    resonance_pdgs=[23]
    resonance_masks=[ibset(ibset(0,3),4)]
    resonance_masses=[91.188d0]
    resonance_widths=[2.441404d0]
    call phase_space%configure_resonances(1,resonance_pdgs,resonance_masks,&
         resonance_masses,resonance_widths)
    call phase_space%init(500d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
         sqrt_s_min,.false.,.false.)
    if (phase_space%ndim.ne.5 .or. phase_space%ndim_extra.ne.0) then
       write (*,*) 'Unexpected resonance-map dimensionality',phase_space%ndim,&
            phase_space%ndim_extra
       stop 1
    endif
    allocate(point%x(phase_space%ndim))
    allocate(point%p(0:3,n))
    original_x=[0.37d0,0.21d0,0.63d0,0.44d0,0.72d0]
    point%x=original_x
    call phase_space%generate_momenta(point)
    if (point%jac.le.0d0) then
       write (*,*) 'Forward Breit-Wigner map failed',point%jac
       stop 1
    endif
    pair_momentum=point%p(:,4)+point%p(:,5)
    pair_mass2=minkowski_square(pair_momentum)
    theta_min=atan((0d0-resonance_masses(1)**2)/&
         (resonance_masses(1)*resonance_widths(1)))
    theta_max=atan((500d0**2-resonance_masses(1)**2)/&
         (resonance_masses(1)*resonance_widths(1)))
    theta=theta_min+(theta_max-theta_min)*original_x(1)
    expected_mass2=resonance_masses(1)**2+&
         resonance_masses(1)*resonance_widths(1)*tan(theta)
    if (abs(pair_mass2-expected_mass2).gt.1d-9*max(1d0,abs(expected_mass2))) then
       write (*,*) 'Breit-Wigner quantile mismatch',pair_mass2,expected_mass2
       stop 1
    endif
    forward_jac=point%jac
    allocate(inverse_point%x(phase_space%ndim))
    allocate(inverse_point%p(0:3,n))
    inverse_point%p=point%p
    call phase_space%compute_x_from_momenta(inverse_point)
    call check_inverse(original_x,forward_jac,inverse_point)
    call phase_space%cleanup()
  end subroutine check_breit_wigner_quantile

  subroutine check_massive_quark_external_legs()
    type(phase_space_gen23) :: phase_space
    type(psv) :: point,inverse_point
    integer,parameter :: n=5
    integer,dimension(n) :: order
    integer,dimension(1) :: resonance_pdgs,resonance_masks
    real(kind=8),dimension(1) :: resonance_masses,resonance_widths
    real(kind=8),dimension(n) :: masses,pt_cut,rap_cut
    real(kind=8),dimension(n,n) :: dr_cut,sqrt_s_min
    real(kind=8),dimension(5) :: original_x
    real(kind=8) :: forward_jac,mass_shell
    integer :: i

    order=[1,2,3,4,5]
    masses=[0d0,0d0,4.7d0,1.42d0,1.42d0]
    pt_cut=0d0
    rap_cut=0d0
    dr_cut=0d0
    sqrt_s_min=0d0
    resonance_pdgs=[23]
    resonance_masks=[ibset(ibset(0,3),4)]
    resonance_masses=[91.188d0]
    resonance_widths=[2.441404d0]
    call phase_space%configure_resonances(1,resonance_pdgs,resonance_masks,&
         resonance_masses,resonance_widths)
    call phase_space%init(500d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
         sqrt_s_min,.false.,.false.)
    allocate(point%x(phase_space%ndim),point%p(0:3,n))
    original_x=[0.29d0,0.67d0,0.38d0,0.81d0,0.46d0]
    point%x=original_x
    call phase_space%generate_momenta(point)
    if (point%jac.le.0d0) then
       write (*,*) 'Massive-quark resonance map failed',point%jac
       stop 1
    endif
    do i=3,n
       mass_shell=minkowski_square(point%p(:,i))
       if (abs(mass_shell-masses(i)**2).gt.1d-9*max(1d0,masses(i)**2)) then
          write (*,*) 'Massive external leg is off shell',i,mass_shell,masses(i)**2
          stop 1
       endif
    enddo
    forward_jac=point%jac
    allocate(inverse_point%x(phase_space%ndim),inverse_point%p(0:3,n))
    inverse_point%p=point%p
    call phase_space%compute_x_from_momenta(inverse_point)
    call check_inverse(original_x,forward_jac,inverse_point)
    call phase_space%cleanup()
  end subroutine check_massive_quark_external_legs

  subroutine check_root_breit_wigner_with_pdfs()
    type(phase_space_gen23) :: phase_space
    type(psv) :: point,inverse_point
    integer,parameter :: n=4
    integer,dimension(n) :: order
    integer,dimension(1) :: resonance_pdgs,resonance_masks
    real(kind=8),dimension(1) :: resonance_masses,resonance_widths
    real(kind=8),dimension(n) :: masses,pt_cut,rap_cut
    real(kind=8),dimension(n,n) :: dr_cut,sqrt_s_min
    real(kind=8),dimension(4) :: original_x
    real(kind=8),dimension(0:3) :: total_momentum
    real(kind=8) :: invariant_mass2,expected_mass2,theta_min,theta_max,theta
    real(kind=8) :: forward_jac

    order=[1,2,3,4]
    masses=0d0
    pt_cut=0d0
    rap_cut=0d0
    dr_cut=0d0
    sqrt_s_min=0d0
    resonance_pdgs=[23]
    resonance_masks=[ibset(ibset(0,2),3)]
    resonance_masses=[91.188d0]
    resonance_widths=[2.441404d0]
    call phase_space%configure_resonances(1,resonance_pdgs,resonance_masks,&
         resonance_masses,resonance_widths)
    call phase_space%init(700d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
         sqrt_s_min,.false.,.true.)
    allocate(point%x(phase_space%ndim))
    allocate(point%p(0:3,n))
    original_x=[0.41d0,0.56d0,0.33d0,0.77d0]
    point%x=original_x
    call phase_space%generate_momenta(point)
    if (point%jac.le.0d0) then
       write (*,*) 'Root Breit-Wigner map failed',point%jac
       stop 1
    endif
    total_momentum=point%p(:,3)+point%p(:,4)
    invariant_mass2=minkowski_square(total_momentum)
    theta_min=atan((0d0-resonance_masses(1)**2)/&
         (resonance_masses(1)*resonance_widths(1)))
    theta_max=atan((700d0**2-resonance_masses(1)**2)/&
         (resonance_masses(1)*resonance_widths(1)))
    theta=theta_min+(theta_max-theta_min)*original_x(1)
    expected_mass2=resonance_masses(1)**2+&
         resonance_masses(1)*resonance_widths(1)*tan(theta)
    if (abs(invariant_mass2-expected_mass2).gt.&
         1d-9*max(1d0,abs(expected_mass2))) then
       write (*,*) 'Root Breit-Wigner quantile mismatch',invariant_mass2,expected_mass2
       stop 1
    endif
    forward_jac=point%jac
    allocate(inverse_point%x(phase_space%ndim))
    allocate(inverse_point%p(0:3,n))
    inverse_point%p=point%p
    call phase_space%compute_x_from_momenta(inverse_point)
    call check_inverse(original_x,forward_jac,inverse_point)
    call phase_space%cleanup()
  end subroutine check_root_breit_wigner_with_pdfs

  subroutine check_nested_and_disjoint_inverse()
    type(phase_space_gen23) :: phase_space
    type(psv) :: point,inverse_point
    integer,parameter :: n=8,nres=4
    integer,dimension(n) :: order
    integer,dimension(nres) :: resonance_pdgs,resonance_masks
    real(kind=8),dimension(nres) :: resonance_masses,resonance_widths
    real(kind=8),dimension(n) :: masses,pt_cut,rap_cut
    real(kind=8),dimension(n,n) :: dr_cut,sqrt_s_min
    real(kind=8),dimension(14) :: original_x
    real(kind=8) :: forward_jac
    integer :: i

    order=[1,2,3,4,5,6,7,8]
    masses=0d0
    pt_cut=0d0
    rap_cut=0d0
    dr_cut=0d0
    sqrt_s_min=0d0
    resonance_pdgs=[24,6,23,23]
    resonance_masks=[&
         ibset(ibset(0,3),4),&
         ibset(ibset(ibset(0,2),3),4),&
         ibset(ibset(0,5),6),&
         ibset(ibset(ibset(0,5),6),7)]
    ! Two disjoint branches each contain a nested pole: W=(4,5) inside
    ! top=(3,4,5), and Z=(6,7) inside an outer neutral current=(6,7,8).
    resonance_masses=[80.419002445756163d0,173d0,91.188d0,250d0]
    resonance_widths=[2.0476d0,1.4915d0,2.441404d0,8d0]
    call phase_space%configure_resonances(nres,resonance_pdgs,resonance_masks,&
         resonance_masses,resonance_widths)
    call phase_space%init(700d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
         sqrt_s_min,.false.,.false.)
    do i=1,size(original_x)
       original_x(i)=0.08d0+0.83d0*dble(i)/dble(size(original_x)+1)
    enddo
    allocate(point%x(phase_space%ndim))
    allocate(point%p(0:3,n))
    point%x=original_x
    call phase_space%generate_momenta(point)
    if (point%jac.le.0d0) then
       write (*,*) 'Nested/disjoint forward map failed',point%jac
       stop 1
    endif
    forward_jac=point%jac
    allocate(inverse_point%x(phase_space%ndim))
    allocate(inverse_point%p(0:3,n))
    inverse_point%p=point%p
    call phase_space%compute_x_from_momenta(inverse_point)
    call check_inverse(original_x,forward_jac,inverse_point)
    call phase_space%cleanup()
  end subroutine check_nested_and_disjoint_inverse

  subroutine check_alternative_maps_across_scales()
    type(phase_space_gen23) :: generating_map,alternative_map
    type(psv) :: point,alternative_point
    integer,parameter :: n=6
    integer,dimension(n) :: generating_order,alternative_order
    real(kind=8),dimension(n) :: masses,pt_cut,rap_cut
    real(kind=8),dimension(n,n) :: dr_cut,sqrt_s_min
    real(kind=8),dimension(3) :: energy_scales
    real(kind=8),allocatable,dimension(:) :: original_x
    integer :: iscale,i,attempt,nforward
    real(kind=8) :: last_alternative_jac
    logical :: found

    generating_order=[1,3,4,5,6,2]
    alternative_order=[1,4,3,5,6,2]
    energy_scales=[7d0,700d0,7d5]
    masses=0d0
    pt_cut=0d0
    rap_cut=0d0
    dr_cut=0d0
    sqrt_s_min=0d0
    do iscale=1,size(energy_scales)
       call generating_map%init(energy_scales(iscale),n,masses,&
            generating_order,pt_cut,rap_cut,dr_cut,sqrt_s_min,.false.,.false.)
       call alternative_map%init(energy_scales(iscale),n,masses,&
            alternative_order,pt_cut,rap_cut,dr_cut,sqrt_s_min,.false.,.false.)
       allocate(point%x(generating_map%ndim+generating_map%ndim_extra))
       allocate(point%p(0:3,n))
       allocate(original_x(size(point%x)))
       allocate(alternative_point%x(alternative_map%ndim+&
            alternative_map%ndim_extra))
       allocate(alternative_point%p(0:3,n))
       found=.false.
       nforward=0
       last_alternative_jac=-999d0
       do attempt=1,400
          do i=1,size(original_x)
             original_x(i)=0.05d0+0.9d0*&
                  modulo(dble(37*attempt+53*i),997d0)/997d0
          enddo
          point%x=original_x
          call generating_map%generate_momenta(point)
          if (point%jac.le.0d0 .or. .not.ieee_is_finite(point%jac)) cycle
          nforward=nforward+1
          alternative_point%x=0.5d0
          alternative_point%p=point%p
          call alternative_map%compute_x_from_momenta(alternative_point)
          last_alternative_jac=alternative_point%jac
          if (alternative_point%jac.le.0d0 .or. &
               .not.ieee_is_finite(alternative_point%jac)) cycle
          if (any(.not.ieee_is_finite(&
               alternative_point%x(1:alternative_map%ndim)))) cycle
          found=.true.
          exit
       enddo
       if (.not.found) then
          write (*,*) 'Could not find a valid alternative-map inversion at scale',&
               energy_scales(iscale),nforward,last_alternative_jac
          stop 1
       endif
       deallocate(original_x)
       deallocate(point%x,point%p)
       deallocate(alternative_point%x,alternative_point%p)
       call generating_map%cleanup()
       call alternative_map%cleanup()
    enddo
  end subroutine check_alternative_maps_across_scales

  subroutine check_singular_boundary_and_nonfinite_rejection()
    type(phase_space_gen23) :: phase_space
    type(psv) :: point,inverse_point
    integer,parameter :: n=6
    integer,dimension(n) :: order
    real(kind=8),dimension(n) :: masses,pt_cut,rap_cut
    real(kind=8),dimension(n,n) :: dr_cut,sqrt_s_min
    integer :: i,attempt
    logical :: found

    order=[1,3,4,5,6,2]
    masses=0d0
    pt_cut=0d0
    rap_cut=0d0
    dr_cut=0d0
    sqrt_s_min=0d0
    call phase_space%init(700d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
         sqrt_s_min,.false.,.false.)
    allocate(point%x(phase_space%ndim+phase_space%ndim_extra))
    allocate(point%p(0:3,n))
    found=.false.
    do attempt=1,400
       do i=1,size(point%x)
          point%x(i)=0.05d0+0.9d0*&
               modulo(dble(37*attempt+53*i),997d0)/997d0
       enddo
       call phase_space%generate_momenta(point)
       if (point%jac.gt.0d0 .and. ieee_is_finite(point%jac)) then
          found=.true.
          exit
       endif
    enddo
    if (.not.found) then
       write (*,*) 'Forward point for inverse rejection tests failed',point%jac
       stop 1
    endif
    allocate(inverse_point%x(phase_space%ndim+phase_space%ndim_extra))
    allocate(inverse_point%p(0:3,n))

    ! Collapse one external leg into its neighbour while preserving total
    ! momentum.  This is a kinematic boundary/singular configuration and must
    ! simply have zero density in an alternative map.
    inverse_point%x=0.5d0
    inverse_point%p=point%p
    inverse_point%p(:,4)=inverse_point%p(:,4)+inverse_point%p(:,3)
    inverse_point%p(:,3)=0d0
    call phase_space%compute_x_from_momenta(inverse_point)
    if (inverse_point%jac.ge.0d0 .or. .not.ieee_is_finite(inverse_point%jac)) then
       write (*,*) 'Singular inverse point was not rejected cleanly',&
            inverse_point%jac
       stop 1
    endif

    inverse_point%x=0.5d0
    inverse_point%p=point%p
    inverse_point%p(0,3)=ieee_value(0d0,ieee_quiet_nan)
    call phase_space%compute_x_from_momenta(inverse_point)
    if (inverse_point%jac.ge.0d0 .or. .not.ieee_is_finite(inverse_point%jac)) then
       write (*,*) 'Non-finite inverse point was not rejected cleanly',&
            inverse_point%jac
       stop 1
    endif
    call phase_space%cleanup()

    call check_exact_random_boundary()
  end subroutine check_singular_boundary_and_nonfinite_rejection

  subroutine check_exact_random_boundary()
    type(phase_space_gen23) :: phase_space
    type(psv) :: inverse_point
    integer,parameter :: n=4
    integer,dimension(n) :: order
    real(kind=8),dimension(n) :: masses,pt_cut,rap_cut
    real(kind=8),dimension(n,n) :: dr_cut,sqrt_s_min

    order=[1,3,2,4]
    masses=0d0
    pt_cut=0d0
    rap_cut=0d0
    dr_cut=0d0
    sqrt_s_min=0d0
    call phase_space%init(700d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
         sqrt_s_min,.false.,.false.)
    allocate(inverse_point%x(phase_space%ndim+phase_space%ndim_extra))
    allocate(inverse_point%p(0:3,n))
    inverse_point%x=0.5d0
    inverse_point%p=0d0
    inverse_point%p(:,1)=[350d0,0d0,0d0,350d0]
    inverse_point%p(:,2)=[350d0,0d0,0d0,-350d0]
    inverse_point%p(:,3)=inverse_point%p(:,1)
    inverse_point%p(:,4)=inverse_point%p(:,2)
    call phase_space%compute_x_from_momenta(inverse_point)
    if (inverse_point%jac.ge.0d0 .or. .not.ieee_is_finite(inverse_point%jac)) then
       write (*,*) 'Exact random-variable boundary was not rejected cleanly',&
            inverse_point%jac
       stop 1
    endif
    call phase_space%cleanup()
  end subroutine check_exact_random_boundary

  subroutine check_inverse(expected_x,forward_jac,inverse_point)
    real(kind=8),dimension(:),intent(in) :: expected_x
    real(kind=8),intent(in) :: forward_jac
    type(psv),intent(in) :: inverse_point
    if (inverse_point%jac.le.0d0) then
       write (*,*) 'Inverse resonance map failed',inverse_point%jac
       stop 1
    endif
    if (maxval(abs(inverse_point%x-expected_x)).gt.2d-9) then
       write (*,*) 'Forward/inverse random variables disagree'
       write (*,*) expected_x
       write (*,*) inverse_point%x
       stop 1
    endif
    if (abs(inverse_point%jac-forward_jac).gt.&
         2d-9*max(1d0,abs(forward_jac))) then
       write (*,*) 'Forward/inverse Jacobians disagree',forward_jac,inverse_point%jac
       stop 1
    endif
  end subroutine check_inverse

  real(kind=8) function minkowski_square(momentum)
    real(kind=8),dimension(0:3),intent(in) :: momentum
    minkowski_square=momentum(0)**2-sum(momentum(1:3)**2)
  end function minkowski_square

end program resonance_phase_space_regression
