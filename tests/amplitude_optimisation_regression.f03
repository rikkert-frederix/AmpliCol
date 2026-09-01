program amplitude_optimisation_regression
  use amplitude_QCD_mod
  use particles
  use run_parameters, only: reset_run_parameters
  implicit none
  integer,parameter :: dp=kind(1d0)
  type(physics_model) :: model
  character(len=64) :: test_mode

  call reset_run_parameters()
  call model%init_part()
  call model%init_vert()
  open(unit=99,file='/dev/null',status='unknown',action='write')
  call get_command_argument(1,test_mode)
  if (len_trim(test_mode).gt.0) then
     call check_invalid_filter_state(trim(test_mode))
     write (*,*) 'Malformed helicity-filter state was not rejected:',trim(test_mode)
     stop
  endif

  call check_process([21,21,21,21],[1,2,3,4],0d0,0d0,&
       'four-gluon parity pairs',.true.)
  call check_process([2,-1,21,24],[2,3,4,1],0d0,model%get_mass(24),&
       'chiral W plus gluon',.false.)

  write (*,'(a)') 'Amplitude optimisation regression passed'

contains

  subroutine check_invalid_filter_state(mode)
    implicit none
    character(len=*),intent(in) :: mode
    type(amplitude_QCD) :: amp
    real(kind=dp) :: samples(3,2)
    integer :: include_hel(3)

    samples=1d0
    select case (mode)
    case ('missing-filter-offsets')
       amp%n_amps=3
       amp%nprocs=1
    case ('nonmonotone-filter-offsets')
       amp%n_amps=3
       amp%nprocs=3
       allocate(amp%iproc_start(4))
       amp%iproc_start=[1,3,2,4]
    case default
       write (*,*) 'Unknown amplitude-optimisation regression mode:',trim(mode)
       stop 1
    end select
    call build_helicity_filter(amp,samples,include_hel,.true.)
  end subroutine check_invalid_filter_state

  subroutine check_process(process,order,mass3,mass4,label,require_degeneracy)
    implicit none
    integer,intent(in) :: process(4),order(4)
    real(kind=dp),intent(in) :: mass3,mass4
    character(len=*),intent(in) :: label
    logical,intent(in) :: require_degeneracy
    type(amplitude_QCD) :: amp
    integer :: spin(0:3,4),helicity(4),processes(4,1),orders(4,1)
    integer :: isample,ihel,jhel,member,nhel,nhel_original,nmatches,matched,representative
    integer,allocatable :: include_hel(:),filter_map(:),original_spins(:,:),original_perm(:,:)
    real(kind=dp) :: p(0:3,4),sample_scale,full_scale,filtered_scale
    real(kind=dp),allocatable :: samples(:,:),full_weights(:),filtered_weights(:)
    complex(kind=dp),allocatable :: full_amplitudes(:)
    logical,allocatable :: seen(:)

    processes(:,1)=process
    orders(:,1)=order
    call setup_spin_states(process,spin)
    helicity=spin(1,:)
    call amp%init(1,4,1,processes,spin,orders,model)
    nhel_original=amp%n_amps
    allocate(samples(nhel_original,10),include_hel(nhel_original),filter_map(nhel_original))

    do isample=1,10
       call fill_two_body_point(800d0+37d0*dble(isample),mass3,mass4,&
            -0.82d0+0.15d0*dble(isample),0.19d0*dble(isample),p)
       call amp%evaluate(4,p,helicity,.false.,model)
       samples(:,isample)=abs(amp%amps)**2
       sample_scale=maxval(samples(:,isample))
       if (sample_scale.le.tiny(1d0)) then
          write (*,*) trim(label),': all helicities vanished at warm-up point',isample
          stop 1
       endif
       samples(:,isample)=samples(:,isample)/sample_scale
       call amp%record_optimisation_sample(isample,10)
    enddo

    call fill_two_body_point(1375d0,mass3,mass4,0.347d0,0.913d0,p)
    call amp%evaluate(4,p,helicity,.false.,model)
    allocate(full_amplitudes(nhel_original),full_weights(nhel_original),&
         original_spins(4,nhel_original),original_perm(size(amp%perm,1),nhel_original))
    full_amplitudes=amp%amps
    full_weights=abs(full_amplitudes)**2
    original_spins=amp%spins(:,1,1:nhel_original)
    original_perm=amp%perm

    call amp%optimise_evaluation(4)
    call amp%evaluate(4,p,helicity,.false.,model)
    do ihel=1,nhel_original
       call assert_complex_close(amp%amps(ihel),full_amplitudes(ihel),&
            trim(label)//' current sharing')
    enddo

    call build_helicity_filter(amp,samples,include_hel,.true.)
    if (require_degeneracy .and. maxval(include_hel).le.1) then
       write (*,*) trim(label),': expected an equivalent-helicity pair'
       stop 1
    endif
    filter_map=include_hel
    nhel=nhel_original
    call amp%filter_helicity(4,nhel,include_hel)
    if (nhel.ne.amp%n_amps) then
       write (*,*) trim(label),': filtered helicity count is inconsistent'
       stop 1
    endif
    if (.not.allocated(amp%include_amp)) then
       write (*,*) trim(label),': filtering discarded serialized inclusion metadata'
       stop 1
    endif
    if (size(amp%include_amp).lt.nhel) then
       write (*,*) trim(label),': filtered inclusion metadata is too short'
       stop 1
    endif
    if (.not.all(amp%include_amp(1:nhel))) then
       write (*,*) trim(label),': retained amplitudes are not marked for inclusion'
       stop 1
    endif
    call amp%evaluate(4,p,helicity,.false.,model)
    allocate(filtered_weights(nhel),seen(nhel_original))
    filtered_weights=abs(amp%amps)**2
    seen=.false.

    do ihel=1,nhel
       representative=0
       nmatches=0
       do jhel=1,nhel_original
          if (filter_map(jhel).le.0) cycle
          nmatches=nmatches+1
          if (nmatches.eq.ihel) representative=jhel
       enddo
       if (representative.eq.0 .or. &
            any(amp%perm(:,ihel).ne.original_perm(:,representative))) then
          write (*,*) trim(label),': retained colour permutation is inconsistent',ihel
          stop 1
       endif
       if (include_hel(ihel).ne.filter_map(representative)) then
          write (*,*) trim(label),': compacted helicity multiplicity is inconsistent',ihel
          stop 1
       endif
       if (include_hel(ihel).lt.1 .or. include_hel(ihel).gt.size(amp%spins,2)) then
          write (*,*) trim(label),': invalid retained-helicity multiplicity',&
               ihel,include_hel(ihel)
          stop 1
       endif
       do member=1,include_hel(ihel)
          nmatches=0
          matched=0
          do jhel=1,nhel_original
             if (all(amp%spins(:,member,ihel).eq.original_spins(:,jhel))) then
                nmatches=nmatches+1
                matched=jhel
             endif
          enddo
          if (nmatches.ne.1) then
             write (*,*) trim(label),': helicity labels were lost or duplicated',&
                  ihel,member,nmatches
             stop 1
          endif
          if (seen(matched)) then
             write (*,*) trim(label),': helicity label was duplicated',ihel,member
             stop 1
          endif
          if ((member.eq.1 .and. matched.ne.representative) .or. &
               (member.gt.1 .and. filter_map(matched).ne.-representative)) then
             write (*,*) trim(label),': helicity label was assigned to the wrong group',&
                  ihel,member,matched,representative
             stop 1
          endif
          seen(matched)=.true.
          call assert_real_close(filtered_weights(ihel),full_weights(matched),&
               trim(label)//' grouped helicity weight')
       enddo
    enddo

    full_scale=maxval(full_weights)
    do ihel=1,nhel_original
       if (full_weights(ihel).gt.100d0*helicity_zero_tolerance*full_scale .and. &
            .not.seen(ihel)) then
          write (*,*) trim(label),': non-zero helicity was removed',&
               original_spins(:,ihel),full_weights(ihel)
          stop 1
       endif
    enddo
    filtered_scale=sum(filtered_weights*dble(include_hel(1:nhel)))
    call assert_real_close(filtered_scale,sum(full_weights),&
         trim(label)//' helicity sum')
  end subroutine check_process

  subroutine setup_spin_states(process,spin)
    implicit none
    integer,intent(in) :: process(4)
    integer,intent(out) :: spin(0:3,4)
    integer :: ipart
    spin=0
    do ipart=1,4
       spin(0,ipart)=model%get_spin(process(ipart))
       select case (spin(0,ipart))
       case (1)
          spin(1,ipart)=0
       case (2)
          spin(1:2,ipart)=[-1,1]
       case (3)
          spin(1:3,ipart)=[-1,0,1]
       case default
          write (*,*) 'Unsupported spin multiplicity in optimisation test',&
               process(ipart),spin(0,ipart)
          stop 1
       end select
    enddo
  end subroutine setup_spin_states

  subroutine fill_two_body_point(sqrts,mass3,mass4,costheta,phi,p)
    implicit none
    real(kind=dp),intent(in) :: sqrts,mass3,mass4,costheta,phi
    real(kind=dp),intent(out) :: p(0:3,4)
    real(kind=dp) :: shat,lambda,q,e3,e4,sintheta,px,py,pz
    shat=sqrts*sqrts
    lambda=(shat-(mass3+mass4)**2)*(shat-(mass3-mass4)**2)
    if (lambda.le.0d0 .or. abs(costheta).ge.1d0) then
       write (*,*) 'Invalid two-body point in amplitude optimisation test'
       stop 1
    endif
    q=sqrt(lambda)/(2d0*sqrts)
    e3=(shat+mass3*mass3-mass4*mass4)/(2d0*sqrts)
    e4=(shat+mass4*mass4-mass3*mass3)/(2d0*sqrts)
    sintheta=sqrt(1d0-costheta*costheta)
    px=q*sintheta*cos(phi)
    py=q*sintheta*sin(phi)
    pz=q*costheta
    p(:,1)=[0.5d0*sqrts,0d0,0d0,0.5d0*sqrts]
    p(:,2)=[0.5d0*sqrts,0d0,0d0,-0.5d0*sqrts]
    p(:,3)=[e3,px,py,pz]
    p(:,4)=[e4,-px,-py,-pz]
  end subroutine fill_two_body_point

  subroutine assert_real_close(value,reference,label)
    implicit none
    real(kind=dp),intent(in) :: value,reference
    character(len=*),intent(in) :: label
    real(kind=dp) :: scale
    scale=max(abs(value),abs(reference),tiny(1d0))
    if (abs(value-reference).gt.2d-10*scale) then
       write (*,*) trim(label),' mismatch:',value,reference
       stop 1
    endif
  end subroutine assert_real_close

  subroutine assert_complex_close(value,reference,label)
    implicit none
    complex(kind=dp),intent(in) :: value,reference
    character(len=*),intent(in) :: label
    real(kind=dp) :: scale
    scale=max(abs(value),abs(reference),tiny(1d0))
    if (abs(value-reference).gt.2d-10*scale) then
       write (*,*) trim(label),' mismatch:',value,reference
       stop 1
    endif
  end subroutine assert_complex_close

end program amplitude_optimisation_regression
