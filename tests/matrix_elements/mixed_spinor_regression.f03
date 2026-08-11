program mixed_spinor_regression
  use amplitude_QCD_mod
  use particles
  implicit none

  integer,parameter :: dp=kind(1d0),n=4,nhelicities=24
  real(kind=dp),parameter :: pi=3.14159265358979323846d0
  real(kind=dp),parameter :: alpha_s=0.118d0
  real(kind=dp),parameter :: alpha_ew=1d0/132.507d0
  real(kind=dp),parameter :: madgraph_averaged=0.10696191393778821d0
  real(kind=dp),parameter :: amplitude_rel_tol=5d-12
  real(kind=dp),parameter :: forbidden_abs_tol=1d-12
  real(kind=dp),parameter :: madgraph_rel_tol=5d-12
  type(physics_model) :: model
  real(kind=dp),dimension(0:3,n) :: momenta
  character(len=64) :: option

  open(unit=99,file='/dev/null',status='unknown',action='write')
  call model%init_part(173d0,1.4915d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call model%init_vert()
  call fill_momenta(momenta)
  if (command_argument_count().gt.0) then
     call get_command_argument(1,option)
     if (trim(option).eq.'--emit-library') then
        call emit_mixed_spinor_libraries()
        stop
     endif
     write (*,*) 'Unknown option: ',trim(option)
     stop 1
  endif

  ! The two orders are the direct-amplitude versions of the process-list
  ! orders: colour singlets precede the coloured leg which closes the current.
  call check_process([5,21,6,-24],[3,2,4,1],[3,1,2,4],&
       [-6,5,21,-24],1,&
       'b g > t W-')
  call check_process([-5,21,-6,24],[1,2,4,3],[1,3,2,4],&
       [-5,6,21,24],-1,&
       'b~ g > t~ W+')
  write (*,'(a)') 'Mixed Weyl/Dirac t-b-W regression passed'

contains

  subroutine check_process(process,order,canonical_mapping,&
       expected_canonical_process,forbidden_bottom_helicity,label)
    implicit none
    integer,dimension(n),intent(in) :: process,order,canonical_mapping
    integer,dimension(n),intent(in) :: expected_canonical_process
    integer,intent(in) :: forbidden_bottom_helicity
    character(len=*),intent(in) :: label
    type(amplitude_QCD) :: mixed_amp,dirac_amp,canonical_amp
    integer,dimension(n,1) :: part,orders,canonical_part,canonical_orders
    integer,dimension(0:3,n) :: all_spins,one_spin
    integer,dimension(n) :: helicity,dummy_helicity,canonical_helicity
    real(kind=dp),dimension(0:3,n) :: canonical_momenta
    logical,dimension(nhelicities) :: seen
    integer :: i,index
    real(kind=dp) :: difference,scale,forbidden_max,dirac_squared
    real(kind=dp) :: mixed_sum,dirac_sum,canonical_sum,normalized

    part(:,1)=process
    orders(:,1)=order
    call setup_all_spins(all_spins)
    dummy_helicity=0
    call mixed_amp%init(1,n,1,part,all_spins,orders,model)
    call mixed_amp%evaluate(n,momenta,dummy_helicity,.false.,model)
    if (mixed_amp%n_amps.ne.nhelicities) then
       write (*,*) trim(label),': expected helicity amplitudes:',&
            nhelicities,'found:',mixed_amp%n_amps
       stop 1
    endif

    ! imode=3 preserves this physical order.  The -9 spin sentinel both
    ! suppresses Weyl tagging during init and takes the requested helicity
    ! from evaluate, giving a pointwise all-Dirac reference.
    one_spin=0
    one_spin(0,:)=1
    one_spin(1,:)=-9
    call dirac_amp%init(3,n,1,part,one_spin,orders,model)
    if (dirac_amp%n_amps.ne.1 .or. dirac_amp%nColOrd.ne.1) then
       write (*,*) trim(label),': unexpected all-Dirac colour basis:',&
            dirac_amp%n_amps,dirac_amp%nColOrd
       stop 1
    endif

    ! Reproduce the reweighter's all-outgoing sort and crossing before using
    ! imode=2.  Direct physical labels are not a canonical imode=2 input.
    call canonicalize_event(process,canonical_mapping,canonical_part(:,1),&
         canonical_momenta)
    if (any(canonical_part(:,1).ne.expected_canonical_process)) then
       write (*,*) trim(label),': unexpected canonical process:',&
            canonical_part(:,1)
       stop 1
    endif
    ! imode=2 derives its own canonical colour order; this argument is unused.
    canonical_orders=0
    call canonical_amp%init(2,n,1,canonical_part,one_spin,&
         canonical_orders,model)
    if (canonical_amp%n_amps.ne.1 .or. canonical_amp%nColOrd.ne.1) then
       write (*,*) trim(label),': unexpected canonical colour basis:',&
            canonical_amp%n_amps,canonical_amp%nColOrd
       stop 1
    endif

    seen=.false.
    mixed_sum=0d0
    dirac_sum=0d0
    canonical_sum=0d0
    forbidden_max=0d0
    do i=1,mixed_amp%n_amps
       helicity=mixed_amp%spins(:,1,i)
       index=helicity_index(helicity)
       if (seen(index)) then
          write (*,*) trim(label),': duplicate helicity:',helicity
          stop 1
       endif
       seen(index)=.true.

       call dirac_amp%evaluate(n,momenta,helicity,.false.,model)
       canonical_helicity=helicity(canonical_mapping)
       call canonical_amp%evaluate(n,canonical_momenta,canonical_helicity,&
            .false.,model)
       dirac_squared=abs(dirac_amp%amps(1))**2
       difference=abs(dirac_squared-abs(canonical_amp%amps(1))**2)
       scale=max(1d0,dirac_squared,abs(canonical_amp%amps(1))**2)
       if (difference.gt.amplitude_rel_tol*scale) then
          write (*,*) trim(label),': fixed-order and canonical amplitudes disagree'
          write (*,*) 'helicity:',helicity
          write (*,*) 'fixed-order squared:',dirac_squared
          write (*,*) 'canonical squared  :',abs(canonical_amp%amps(1))**2
          stop 1
       endif

       difference=abs(mixed_amp%amps(i)-dirac_amp%amps(1))
       scale=max(1d0,abs(mixed_amp%amps(i)),abs(dirac_amp%amps(1)))
       if (difference.gt.amplitude_rel_tol*scale) then
          write (*,*) trim(label),': mixed and Dirac amplitudes disagree'
          write (*,*) 'helicity:',helicity
          write (*,*) 'mixed  :',mixed_amp%amps(i)
          write (*,*) 'Dirac  :',dirac_amp%amps(1)
          stop 1
       endif

       mixed_sum=mixed_sum+abs(mixed_amp%amps(i))**2
       dirac_sum=dirac_sum+dirac_squared
       canonical_sum=canonical_sum+abs(canonical_amp%amps(1))**2
       if (helicity(1).eq.forbidden_bottom_helicity) then
          forbidden_max=max(forbidden_max,abs(mixed_amp%amps(i)))
       endif
    enddo
    if (.not.all(seen)) then
       write (*,*) trim(label),': incomplete helicity coverage'
       stop 1
    endif
    if (forbidden_max.gt.forbidden_abs_tol) then
       write (*,*) trim(label),': forbidden bottom helicity is non-zero:',&
            forbidden_max
       stop 1
    endif

    scale=max(1d0,mixed_sum,dirac_sum,canonical_sum)
    if (abs(mixed_sum-dirac_sum).gt.amplitude_rel_tol*scale .or.&
         abs(dirac_sum-canonical_sum).gt.amplitude_rel_tol*scale) then
       write (*,*) trim(label),': helicity sums disagree:',mixed_sum,&
            dirac_sum,canonical_sum
       stop 1
    endif

    ! Stock MadGraph 3.6.0 SMATRIX for this point, with MB=YMB=0.  SMATRIX
    ! includes the initial helicity and colour averages.  The AmpliCol
    ! amplitudes are coupling-stripped; 12 is their corresponding initial
    ! spin/colour normalization for this single-colour-structure process.
    normalized=mixed_sum*(4d0*pi*alpha_s)*(2d0*4d0*pi*alpha_ew)/12d0
    if (abs(normalized/madgraph_averaged-1d0).gt.madgraph_rel_tol) then
       write (*,*) trim(label),': result disagrees with MadGraph:',&
            normalized,madgraph_averaged
       stop 1
    endif
    write (*,'(a,a,a,es24.16,a,es12.4)') 'MIXED_SPINOR[',trim(label),&
         ']=',normalized,' forbidden-max=',forbidden_max
  end subroutine check_process

  subroutine canonicalize_event(process,mapping,canonical_process,&
       canonical_momenta)
    implicit none
    integer,dimension(n),intent(in) :: process,mapping
    integer,dimension(n),intent(out) :: canonical_process
    real(kind=dp),dimension(0:3,n),intent(out) :: canonical_momenta
    integer,dimension(n) :: crossed_process
    real(kind=dp),dimension(0:3,n) :: crossed_momenta
    integer :: i

    crossed_process=process
    crossed_momenta=momenta
    do i=1,2
       crossed_process(i)=model%get_antipart(crossed_process(i))
       crossed_momenta(:,i)=-crossed_momenta(:,i)
    enddo
    do i=1,n
       canonical_process(i)=crossed_process(mapping(i))
       if (i.le.2) then
          canonical_process(i)=model%get_antipart(canonical_process(i))
          canonical_momenta(:,i)=-crossed_momenta(:,mapping(i))
       else
          canonical_momenta(:,i)=crossed_momenta(:,mapping(i))
       endif
    enddo
  end subroutine canonicalize_event

  subroutine emit_mixed_spinor_libraries()
    implicit none

    call emit_process_library([5,21,6,-24],[3,2,4,1],&
         [-1,-1,-1,-1],1,'b g > t W-')
    call emit_process_library([-5,21,-6,24],[1,2,4,3],&
         [1,1,1,1],2,'b~ g > t~ W+')
    write (*,'(a)') 'Generated mixed-spinor regression libraries'
  end subroutine emit_mixed_spinor_libraries

  subroutine emit_process_library(process,order,helicity,igroup,label)
    implicit none
    integer,dimension(n),intent(in) :: process,order,helicity
    integer,intent(in) :: igroup
    character(len=*),intent(in) :: label
    type(amplitude_QCD) :: amp
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin

    part(:,1)=process
    orders(:,1)=order
    spin=0
    spin(0,:)=1
    spin(1,:)=helicity
    call amp%init(1,n,1,part,spin,orders,model)
    call amp%evaluate(n,momenta,helicity,.false.,model)
    if (amp%n_amps.ne.1 .or. abs(amp%amps(1)).lt.1d-12) then
       write (*,*) trim(label),': cannot create nonzero regression library:',&
            amp%n_amps,amp%amps
       stop 1
    endif
    call amp%create_library(n,helicity,igroup,1,model,momenta)
  end subroutine emit_process_library

  subroutine setup_all_spins(spin)
    implicit none
    integer,dimension(0:3,n),intent(out) :: spin
    spin=0
    spin(0,:)=[2,2,2,3]
    spin(1:2,1)=[-1,1]
    spin(1:2,2)=[-1,1]
    spin(1:2,3)=[-1,1]
    spin(1:3,4)=[-1,0,1]
  end subroutine setup_all_spins

  integer function helicity_index(helicity)
    implicit none
    integer,dimension(n),intent(in) :: helicity
    integer :: w_index
    if (any(abs(helicity(1:3)).ne.1)) then
       write (*,*) 'Invalid fermion/gluon helicity:',helicity
       stop 1
    endif
    select case (helicity(4))
    case (-1)
       w_index=0
    case (0)
       w_index=1
    case (1)
       w_index=2
    case default
       write (*,*) 'Invalid W helicity:',helicity
       stop 1
    end select
    helicity_index=1
    if (helicity(1).eq.1) helicity_index=helicity_index+12
    if (helicity(2).eq.1) helicity_index=helicity_index+6
    if (helicity(3).eq.1) helicity_index=helicity_index+3
    helicity_index=helicity_index+w_index
  end function helicity_index

  subroutine fill_momenta(p)
    implicit none
    real(kind=dp),dimension(0:3,n),intent(out) :: p
    p(:,1)=[500d0,0d0,0d0,500d0]
    p(:,2)=[500d0,0d0,0d0,-500d0]
    p(:,3)=[511.73089202281471d0,144.48029459598254d0,&
         192.64039279464342d0,417.07868488793480d0]
    p(:,4)=[488.26910797718529d0,-144.48029459598254d0,&
         -192.64039279464342d0,-417.07868488793480d0]
  end subroutine fill_momenta

end program mixed_spinor_regression
