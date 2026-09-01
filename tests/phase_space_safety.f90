program phase_space_safety
  use phase_space_base, only: phase_space_type,psv,generated_momenta_are_valid,&
       safe_phase_space_product,safe_phase_space_ratio,safe_phase_space_scaled_ratio,&
       stable_phase_space_power_map,stable_rotate_from_z_axis,stable_rotate_to_z_axis
  use phase_space_gen23_mod, only: phase_space_gen23
  use phase_space_haag_mod, only: phase_space_haag
  use phase_space_genpt_mod, only: phase_space_genpt
  use LUPdecomposition, only: LUPdecompose,LUPsolve,LUPinvert,LUPdeterminant
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none

  type(phase_space_gen23) :: gen23
  type(phase_space_haag) :: haag
  type(phase_space_genpt) :: genpt
  type(psv) :: blank_point

  open(unit=99,status='scratch',action='write')
  call require(blank_point%jac.eq.0d0 .and. all(blank_point%xbjrk.eq.0d0),&
       'fresh phase-space point has undefined scalar state')
  call require_clean_generator(gen23,'fresh gen23')
  call require_clean_generator(haag,'fresh HAAG')
  call require_clean_generator(genpt,'fresh pT phase space')
  call exercise_numeric_helpers()
  call exercise_lup_solver()
  call exercise_gen23(gen23,4,[1,3,2,4],.true.,'gen23 2->2 t-channel')
  call exercise_gen23(gen23,5,[1,3,4,2,5],.false.,'gen23 2->3 Gram channel')
  call exercise_haag(haag,4,[1,3,2,4],.false.,'HAAG 2->2')
  call exercise_haag(haag,4,[1,3,2,4],.true.,'HAAG 2->2 with PDFs')
  call exercise_haag(haag,5,[1,3,4,2,5],.false.,'HAAG 2->3')
  call exercise_genpt(genpt,5,[1,3,4,2,5],'pT phase space 2->3')
  close(99)
  write (*,'(a)') 'phase-space safety regression passed'

contains

  subroutine exercise_numeric_helpers()
    real(kind=8) :: value,var,jac_factor
    real(kind=8) :: momentum(0:3),axis(0:3),rotated(0:3),round_trip(0:3)
    logical :: valid

    valid=safe_phase_space_product(2d0,3d0,value)
    call require(valid .and. value.eq.6d0,'checked phase-space product changed a finite result')
    valid=safe_phase_space_product(huge(1d0),2d0,value)
    call require(.not.valid,'checked phase-space product accepted overflow')
    valid=safe_phase_space_ratio(1d0,tiny(1d0),value)
    call require(valid .and. ieee_is_finite(value),&
         'checked phase-space ratio rejected a representable quotient')
    valid=safe_phase_space_ratio(10d0,tiny(1d0),value)
    call require(.not.valid,'checked phase-space ratio accepted overflow')
    valid=safe_phase_space_scaled_ratio(1d-300,1d300,1d-300,value)
    call require(valid .and. ieee_is_finite(value),&
         'checked scaled ratio rejected a representable divide-first overflow')
    call require_close_scalar(value,1d300,2d-15,&
         'checked scaled ratio changed a divide-first overflow result')
    valid=safe_phase_space_scaled_ratio(1d300,1d-300,1d300,value)
    call require(valid .and. ieee_is_finite(value) .and. value.gt.0d0,&
         'checked scaled ratio rejected a representable divide-first underflow')
    call require(abs(value/1d-300-1d0).le.2d-15,&
         'checked scaled ratio changed a divide-first underflow result')
    valid=safe_phase_space_scaled_ratio(huge(1d0),2d0,1d0,value)
    call require(.not.valid,'checked scaled ratio accepted overflow')

    call stable_phase_space_power_map(1d0,-1.5d0,tiny(1d0),1d102,&
         var,jac_factor,valid)
    call require(valid .and. ieee_is_finite(var) .and. &
         ieee_is_finite(jac_factor),'stable power map rejected a representable endpoint')
    call require_close_scalar(var,1d102,2d-13,&
         'stable power map changed its upper endpoint')
    call stable_phase_space_power_map(1d0,-1.5d0,tiny(1d0),1d103,&
         var,jac_factor,valid)
    call require(.not.valid,'stable power map accepted an unrepresentable Jacobian')

    momentum=[5d0,1d0,-2d0,3d0]
    axis=[10d0,1d-8,-2d-8,4d0]
    call stable_rotate_from_z_axis(momentum,axis,rotated,valid)
    call require(valid,'stable forward rotation rejected a near-axis direction')
    call stable_rotate_to_z_axis(rotated,axis,round_trip,valid)
    call require(valid,'stable inverse rotation rejected a near-axis direction')
    call require(maxval(abs(round_trip-momentum)).le.3d-15*maxval(abs(momentum)),&
         'near-axis forward and inverse rotations are inconsistent')
    axis=[1d-300,1d-300,-2d-300,3d-300]
    call stable_rotate_from_z_axis(momentum,axis,rotated,valid)
    call require(valid,'stable rotation underflowed a uniformly tiny axis')
    call stable_rotate_to_z_axis(rotated,axis,round_trip,valid)
    call require(valid .and. &
         maxval(abs(round_trip-momentum)).le.3d-15*maxval(abs(momentum)),&
         'tiny-axis forward and inverse rotations are inconsistent')
  end subroutine exercise_numeric_helpers

  subroutine exercise_lup_solver()
    real(kind=8) :: original(3,3),lu(3,3),inverse(3,3),identity(3,3)
    real(kind=8) :: solution(3),expected(3),rhs(3),determinant
    integer :: permutation(0:3),i
    logical :: success

    original=reshape([4d0,0d0,2d0,1d0,3d0,0d0,2d0,-1d0,5d0],[3,3])
    expected=[1d0,-2d0,0.5d0]
    rhs=matmul(original,expected)
    lu=original
    call LUPdecompose(lu,3,128d0*epsilon(1d0),permutation,success)
    call require(success,'LUP decomposition rejected a nonsingular matrix')
    call LUPsolve(lu,permutation,rhs,3,solution,success)
    call require(success,'LUP solve failed after a successful decomposition')
    call require(maxval(abs(solution-expected)).le.2d-14,&
         'LUP solve returned the wrong solution')
    call LUPinvert(lu,permutation,3,inverse,success)
    call require(success,'LUP inverse failed after a successful decomposition')
    identity=matmul(original,inverse)
    do i=1,3
       identity(i,i)=identity(i,i)-1d0
    enddo
    call require(maxval(abs(identity)).le.3d-14,'LUP inverse is incorrect')
    call LUPdeterminant(lu,permutation,3,determinant,success)
    call require(success,'LUP determinant failed for a representable result')
    call require_close_scalar(determinant,46d0,2d-14,'LUP determinant is incorrect')

    original(2,:)=original(1,:)
    lu=original
    call LUPdecompose(lu,3,128d0*epsilon(1d0),permutation,success)
    call require(.not.success,'LUP decomposition accepted a singular matrix')
    lu=0d0
    lu(1,1)=huge(1d0)
    lu(2,2)=2d0
    lu(3,3)=1d0
    permutation=[3,1,2,3]
    call LUPdeterminant(lu,permutation,3,determinant,success)
    call require(.not.success,'LUP determinant overflow was not reported')
  end subroutine exercise_lup_solver

  subroutine setup_inputs(n,masses,ptcut,rapcut,drcut,smin)
    integer,intent(in) :: n
    real(kind=8),intent(out) :: masses(n),ptcut(n),rapcut(n),drcut(n,n),smin(n,n)
    integer :: i,j

    masses=0d0
    ptcut=0d0
    rapcut=0d0
    drcut=0d0
    smin=0d0
    ptcut(3:n)=10d0
    do i=3,n-1
       do j=i+1,n
          drcut(i,j)=0.4d0
          drcut(j,i)=drcut(i,j)
       enddo
    enddo
  end subroutine setup_inputs

  subroutine exercise_gen23(generator,n,order,t_channel,label)
    type(phase_space_gen23),intent(inout) :: generator
    integer,intent(in) :: n,order(n)
    logical,intent(in) :: t_channel
    character(len=*),intent(in) :: label
    real(kind=8) :: masses(n),ptcut(n),rapcut(n),drcut(n,n),smin(n,n)
    real(kind=8),allocatable :: input_x(:),reference_p(:,:)
    real(kind=8) :: reference_jac,nan_value
    type(psv) :: point,repeated,inverse,regenerated,short_point

    call setup_inputs(n,masses,ptcut,rapcut,drcut,smin)
    call generator%init(1000d0,n,masses,order,ptcut,rapcut,drcut,smin, &
         t_channel,.false.,flat=.true.)
    call allocate_point(point,n,generator%ndim+generator%ndim_extra)
    call find_valid_point(generator,masses,point,label)
    allocate(input_x(size(point%x)),reference_p(0:3,n))
    input_x=point%x
    reference_p=point%p
    reference_jac=point%jac

    call allocate_point(repeated,n,size(input_x))
    repeated%x=input_x
    call generator%generate_momenta(repeated)
    call require(repeated%jac.gt.0d0,trim(label)//' repeat rejected')
    call require_close(repeated%p,reference_p,2d-12,trim(label)//' is not deterministic')
    call require_close_scalar(repeated%jac,reference_jac,2d-12,trim(label)//' Jacobian changed')

    call allocate_point(inverse,n,size(input_x))
    inverse%p=reference_p
    inverse%xbjrk=point%xbjrk
    call generator%compute_x_from_momenta(inverse)
    call require(inverse%jac.gt.0d0 .and. ieee_is_finite(inverse%jac), &
         trim(label)//' inverse rejected a generated point')
    call require(all(inverse%x(1:generator%ndim).ge.0d0) .and. &
         all(inverse%x(1:generator%ndim).le.1d0),trim(label)//' inverse coordinate outside [0,1]')
    call require_close_scalar(inverse%jac,reference_jac,2d-7, &
         trim(label)//' forward/inverse Jacobian mismatch')

    call allocate_point(regenerated,n,size(input_x))
    regenerated%x=inverse%x
    call generator%generate_momenta(regenerated)
    call require(regenerated%jac.gt.0d0,trim(label)//' inverse point cannot regenerate')
    call require(generated_momenta_are_valid(regenerated%p,masses,regenerated%xbjrk,.false.), &
         trim(label)//' inverse point regenerated invalid momenta')

    nan_value=ieee_value(0d0,ieee_quiet_nan)
    repeated%x=input_x
    repeated%x(1)=nan_value
    call generator%generate_momenta(repeated)
    call require(repeated%jac.eq.-47d0,trim(label)//' accepted a NaN coordinate')
    call require(all(repeated%p.eq.0d0),trim(label)//' exposed output after a NaN coordinate')

    call allocate_point(short_point,n,generator%ndim)
    short_point%x=0.5d0
    call generator%generate_momenta(short_point)
    if (generator%ndim_extra.gt.0) then
       call require(short_point%jac.eq.-47d0,trim(label)//' accepted a short extra-coordinate vector')
    endif
    call generator%cleanup()
    call require_clean_generator(generator,trim(label)//' cleanup')
  end subroutine exercise_gen23

  subroutine exercise_haag(generator,n,order,include_pdf,label)
    type(phase_space_haag),intent(inout) :: generator
    integer,intent(in) :: n,order(n)
    logical,intent(in) :: include_pdf
    character(len=*),intent(in) :: label
    real(kind=8) :: masses(n),ptcut(n),rapcut(n),drcut(n,n),smin(n,n)
    real(kind=8),allocatable :: input_x(:),reference_p(:,:)
    real(kind=8) :: reference_jac,nan_value
    type(psv) :: point,repeated,inverse

    call setup_inputs(n,masses,ptcut,rapcut,drcut,smin)
    call generator%init(1000d0,n,masses,order,ptcut,rapcut,drcut,smin, &
         .false.,include_pdf,flat=.true.)
    call require(generator%ndim_extra.gt.0,trim(label)//' did not declare discrete flat choices')
    call allocate_point(point,n,generator%ndim+generator%ndim_extra)
    call find_valid_point(generator,masses,point,label)
    allocate(input_x(size(point%x)),reference_p(0:3,n))
    input_x=point%x
    reference_p=point%p
    reference_jac=point%jac

    call allocate_point(repeated,n,size(input_x))
    repeated%x=input_x
    call generator%generate_momenta(repeated)
    call require(repeated%jac.gt.0d0,trim(label)//' repeat rejected')
    call require_close(repeated%p,reference_p,2d-12,trim(label)//' is not deterministic')
    call require_close_scalar(repeated%jac,reference_jac,2d-12,trim(label)//' Jacobian changed')

    call allocate_point(inverse,n,size(input_x))
    inverse%p=reference_p
    inverse%x=0.5d0
    call generator%compute_x_from_momenta(inverse)
    call require(inverse%jac.eq.-14d0,trim(label)//' unsupported inverse did not return a status')
    call require(all(inverse%x.eq.0d0),trim(label)//' unsupported inverse retained coordinates')

    nan_value=ieee_value(0d0,ieee_quiet_nan)
    repeated%x=input_x
    repeated%x(generator%ndim+1)=nan_value
    call generator%generate_momenta(repeated)
    call require(repeated%jac.eq.-13d0,trim(label)//' accepted a NaN flat coordinate')
    call require(all(repeated%p.eq.0d0),trim(label)//' exposed output after a NaN flat coordinate')
    call generator%cleanup()
    call require_clean_generator(generator,trim(label)//' cleanup')
  end subroutine exercise_haag

  subroutine exercise_genpt(generator,n,order,label)
    type(phase_space_genpt),intent(inout) :: generator
    integer,intent(in) :: n,order(n)
    character(len=*),intent(in) :: label
    real(kind=8),parameter :: low_scale=1d-20,high_scale=1d20
    real(kind=8) :: masses(n),ptcut(n),rapcut(n),drcut(n,n),smin(n,n),nan_value
    real(kind=8),allocatable :: input_x(:),reference_p(:,:)
    type(psv) :: point,repeated,inverse,scaled

    call setup_inputs(n,masses,ptcut,rapcut,drcut,smin)
    rapcut(3:n)=5d0
    call generator%init(1000d0,n,masses,order,ptcut,rapcut,drcut,smin, &
         .false.,.true.,flat=.true.)
    call allocate_point(point,n,generator%ndim+generator%ndim_extra)
    call find_valid_point(generator,masses,point,label)
    allocate(input_x(size(point%x)),reference_p(0:3,n))
    input_x=point%x
    reference_p=point%p

    call allocate_point(repeated,n,size(input_x))
    repeated%x=input_x
    call generator%generate_momenta(repeated)
    call require(repeated%jac.gt.0d0,trim(label)//' repeat rejected')
    call require_close(repeated%p,reference_p,2d-12,trim(label)//' is not deterministic')

    call allocate_point(inverse,n,size(input_x))
    inverse%p=reference_p
    inverse%x=0.5d0
    call generator%compute_x_from_momenta(inverse)
    call require(inverse%jac.eq.-14d0,trim(label)//' unsupported inverse did not return a status')
    call require(all(inverse%x.eq.0d0),trim(label)//' unsupported inverse retained coordinates')

    nan_value=ieee_value(0d0,ieee_quiet_nan)
    repeated%x=input_x
    repeated%x(1)=nan_value
    call generator%generate_momenta(repeated)
    call require(repeated%jac.eq.-10d0,trim(label)//' accepted a NaN coordinate')
    call require(all(repeated%p.eq.0d0),trim(label)//' exposed output after a NaN coordinate')

    call generator%cleanup()
    ptcut=low_scale*ptcut
    call generator%init(low_scale*1000d0,n,masses,order,ptcut,rapcut,drcut,smin, &
         .false.,.true.,flat=.true.)
    call allocate_point(scaled,n,size(input_x))
    scaled%x=input_x
    call generator%generate_momenta(scaled)
    call require(scaled%jac.gt.0d0,trim(label)//' low-scale point rejected')
    call require_close(scaled%p/low_scale,reference_p,2d-10, &
         trim(label)//' is not covariant at low scale')

    call generator%cleanup()
    ptcut=(high_scale/low_scale)*ptcut
    call generator%init(high_scale*1000d0,n,masses,order,ptcut,rapcut,drcut,smin, &
         .false.,.true.,flat=.true.)
    scaled%x=input_x
    call generator%generate_momenta(scaled)
    call require(scaled%jac.gt.0d0,trim(label)//' high-scale point rejected')
    call require_close(scaled%p/high_scale,reference_p,2d-10, &
         trim(label)//' is not covariant at high scale')
    call generator%cleanup()
    call require_clean_generator(generator,trim(label)//' cleanup')
  end subroutine exercise_genpt

  subroutine require_clean_generator(generator,label)
    class(phase_space_type),intent(in) :: generator
    character(len=*),intent(in) :: label
    call require(generator%jac.eq.0d0 .and. all(generator%xbjrk.eq.0d0),&
         trim(label)//' retained point scalars')
    call require(generator%s0.eq.0d0 .and. generator%tot_mass.eq.0d0 .and. &
         generator%sqrtshat.eq.0d0 .and. generator%sqrts.eq.0d0,&
         trim(label)//' retained energy state')
    call require(generator%ndim.eq.0 .and. generator%next.eq.0 .and. &
         generator%ndim_extra.eq.0,trim(label)//' retained dimensions')
    call require(.not.generator%t_channel .and. .not.generator%can_invert_momenta,&
         trim(label)//' retained mode state')
  end subroutine require_clean_generator

  subroutine allocate_point(point,n,nrandom)
    type(psv),intent(inout) :: point
    integer,intent(in) :: n,nrandom
    if (allocated(point%p)) deallocate(point%p)
    if (allocated(point%x)) deallocate(point%x)
    allocate(point%p(0:3,n),point%x(nrandom))
    point%p=0d0
    point%x=0d0
    point%xbjrk=1d0
    point%jac=-1d0
  end subroutine allocate_point

  subroutine find_valid_point(generator,masses,point,label)
    class(phase_space_type),intent(inout) :: generator
    real(kind=8),intent(in) :: masses(:)
    type(psv),intent(inout) :: point
    character(len=*),intent(in) :: label
    integer :: attempt,i

    do attempt=1,20000
       do i=1,size(point%x)
          point%x(i)=0.05d0+0.9d0*mod(0.7548776662466927d0*dble(attempt*(i+3))+ &
               0.5698402909980532d0*dble(i),1d0)
       enddo
       call generator%generate_momenta(point)
       if (point%jac.gt.0d0) exit
    enddo
    call require(point%jac.gt.0d0,trim(label)//' found no valid deterministic point')
    call require(ieee_is_finite(point%jac),trim(label)//' returned a non-finite Jacobian')
    call require(generated_momenta_are_valid(point%p,masses,point%xbjrk,.false.), &
         trim(label)//' returned invalid momenta')
  end subroutine find_valid_point

  subroutine require_close(actual,expected,tolerance,label)
    real(kind=8),intent(in) :: actual(0:,:),expected(0:,:),tolerance
    character(len=*),intent(in) :: label
    real(kind=8) :: scale
    scale=max(1d0,maxval(abs(expected)))
    call require(maxval(abs(actual-expected)).le.tolerance*scale,label)
  end subroutine require_close

  subroutine require_close_scalar(actual,expected,tolerance,label)
    real(kind=8),intent(in) :: actual,expected,tolerance
    character(len=*),intent(in) :: label
    call require(abs(actual-expected).le.tolerance*max(1d0,abs(expected)),label)
  end subroutine require_close_scalar

  subroutine require(condition,message)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: message
    if (.not.condition) then
       write (*,'(a)') 'FAIL: '//trim(message)
       stop 1
    endif
  end subroutine require

end program phase_space_safety
