program phase_space_gen23_regression
  use phase_space_gen23_mod, only: phase_space_gen23,&
       set_use_soft_bounds_as_actual_limits
  use phase_space_base, only: psv
  implicit none

  integer,parameter :: n=4
  integer,parameter :: all_final_mask=12
  integer,parameter :: initial_final_mask=9
  integer,parameter :: complementary_mask=6
  type(phase_space_gen23) :: generator
  real(kind=8) :: masses(n),pt_cut(n),rap_cut(n)
  real(kind=8) :: dr_cut(n,n),sqrt_s_min(n,n)
  integer :: order(n)

  open(unit=99,status='scratch',action='write')
  call set_use_soft_bounds_as_actual_limits(.true.)
  order=[1,3,2,4]
  rap_cut=0d0
  dr_cut=0d0

  ! For two massless final particles, M_final^2=s_34.  A pair cut of
  ! sqrt(s_34)>=10 must therefore set the lower bound to 100, not 200.
  masses=0d0
  pt_cut=0d0
  sqrt_s_min=0d0
  sqrt_s_min(3,4)=10d0
  sqrt_s_min(4,3)=10d0
  call generator%init(1000d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
       sqrt_s_min,.false.,.true.)
  call assert_close(generator%invm_min(all_final_mask,1),100d0,&
       'all-final hard lower bound')
  call assert_close(generator%invm_min(all_final_mask,2),100d0,&
       'all-final soft lower bound')
  call generator%cleanup()

  ! The final cut uses (p_1+p_4)^2>=Q^2 while the phase-space mask stores
  ! t_14=(p_4-p_1)^2=2*m_4^2-(p_1+p_4)^2.  For m_4=30, Q=100 and pT=20,
  ! the strongest upper bound is min(-20^2,2*30^2-100^2)=-8200.
  masses=0d0
  masses(4)=30d0
  pt_cut=0d0
  pt_cut(4)=20d0
  sqrt_s_min=0d0
  sqrt_s_min(1,4)=100d0
  call generator%init(1000d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
       sqrt_s_min,.false.,.true.)
  call assert_close(generator%invm_max(initial_final_mask,1),-8200d0,&
       'massive initial-final hard upper bound')
  call assert_close(generator%invm_max(initial_final_mask,2),-8200d0,&
       'massive initial-final soft upper bound')
  call assert_close(generator%invm_max(complementary_mask,1),-8200d0,&
       'complementary massive initial-final bound')
  call generator%cleanup()

  ! In a two-body final state, {1,4} and {2,3} are complementary masks for
  ! the same transfer.  Both physical pair cuts constrain it, so retain the
  ! stronger upper bound regardless of mask traversal order.
  sqrt_s_min(2,3)=120d0
  call generator%init(1000d0,n,masses,order,pt_cut,rap_cut,dr_cut,&
       sqrt_s_min,.false.,.true.)
  call assert_close(generator%invm_max(initial_final_mask,1),-14400d0,&
       'strongest complementary initial-final bound')
  call assert_close(generator%invm_max(complementary_mask,1),-14400d0,&
       'strongest reverse complementary bound')
  call generator%cleanup()

  call test_forward_inverse_closure()
  close(99)

contains

  subroutine assert_close(actual,expected,label)
    real(kind=8),intent(in) :: actual,expected
    character(len=*),intent(in) :: label
    real(kind=8) :: tolerance
    tolerance=1d-12*max(1d0,abs(expected))
    if (abs(actual-expected).gt.tolerance) then
       write (*,*) 'FAIL:',trim(label),'expected',expected,'got',actual
       stop 1
    endif
  end subroutine assert_close

  subroutine test_forward_inverse_closure()
    integer,parameter :: nclosure=5
    type(phase_space_gen23) :: closure_generator
    type(psv) :: point
    real(kind=8) :: closure_masses(nclosure),closure_pt(nclosure)
    real(kind=8) :: closure_rap(nclosure),closure_dr(nclosure,nclosure)
    real(kind=8) :: closure_smin(nclosure,nclosure),forward_jacobian
    integer :: closure_order(nclosure),i

    closure_masses=0d0
    closure_pt=0d0
    closure_rap=0d0
    closure_dr=0d0
    closure_smin=0d0
    closure_order=[1,3,4,2,5]
    call closure_generator%init(1000d0,nclosure,closure_masses,closure_order,&
         closure_pt,closure_rap,closure_dr,closure_smin,.false.,.false.)
    allocate(point%x(closure_generator%ndim+closure_generator%ndim_extra))
    allocate(point%p(0:3,nclosure))
    do i=1,size(point%x)
       point%x(i)=dble(mod(37*i+11,89))/89d0
    enddo
    call closure_generator%generate_momenta(point)
    if (point%jac.le.0d0) then
       write (*,*) 'FAIL: unable to generate BK closure point',point%jac
       stop 1
    endif
    forward_jacobian=point%jac
    call closure_generator%compute_x_from_momenta(point)
    if (point%jac.le.0d0) then
       write (*,*) 'FAIL: unable to invert BK closure point',point%jac
       stop 1
    endif
    call assert_close(point%jac,forward_jacobian,'BK forward/inverse Jacobian')
    call closure_generator%cleanup()
  end subroutine test_forward_inverse_closure

end program phase_space_gen23_regression
