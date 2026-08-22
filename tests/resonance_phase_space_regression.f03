program resonance_phase_space_regression
  use phase_space_base
  use phase_space_gen23_mod
  implicit none

  open(unit=99,status='scratch',action='write')
  call check_breit_wigner_quantile()
  call check_massive_quark_external_legs()
  call check_root_breit_wigner_with_pdfs()
  call check_nested_and_disjoint_inverse()
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
