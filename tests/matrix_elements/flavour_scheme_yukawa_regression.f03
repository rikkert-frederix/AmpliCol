program flavour_scheme_yukawa_regression
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amplitude_QCD_mod
  use particles
  use run_parameters
  implicit none

  integer,parameter :: dp=kind(1d0)
  real(kind=dp) :: bottom_result,charm_result,top_decay_result

  open(unit=99,status='scratch',action='write')
  call evaluate_yukawa_process(4,5,bottom_mass,bottom_result)
  call evaluate_yukawa_process(3,4,charm_mass,charm_result)
  call evaluate_massive_bottom_top_decay(top_decay_result)
  close(99)
  write (*,'(a,3es24.16)') 'FLAVOUR_SCHEME_AMPLITUDE_CHECK=',&
       bottom_result,charm_result,top_decay_result
  write (*,'(a)') 'Flavour-scheme matrix-element regression passed'

contains

  subroutine evaluate_yukawa_process(scheme,pdg,mass,result)
    implicit none
    integer,intent(in) :: scheme,pdg
    real(kind=dp),intent(in) :: mass
    real(kind=dp),intent(out) :: result
    integer,parameter :: n=3
    type(physics_model) :: model
    type(amplitude_QCD) :: amplitude
    integer,dimension(n,1) :: process,order
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: helicity
    real(kind=dp),dimension(0:3,n) :: momenta
    real(kind=dp) :: energy,pz,yukawa,expected,tolerance

    call set_flavour_scheme(scheme)
    call model%init_part()
    call model%init_vert()
    if (abs(model%get_mass(pdg)-mass).gt.1d-14) then
       write (*,*) 'Unexpected effective quark mass',scheme,pdg,model%get_mass(pdg),mass
       stop 1
    endif

    ! The amplitude colour-order convention crosses the two incoming legs:
    ! an incoming anti-quark opens the fermion line and the quark closes it.
    process(:,1)=[-pdg,pdg,25]
    order(:,1)=[1,3,2]
    spin=0
    spin(0,:)=[2,2,1]
    spin(1:2,1)=[-1,1]
    spin(1:2,2)=[-1,1]
    spin(1,3)=0
    helicity=0

    energy=model%get_mass(25)/2d0
    pz=sqrt(energy**2-mass**2)
    momenta(:,1)=[energy,0d0,0d0,pz]
    momenta(:,2)=[energy,0d0,0d0,-pz]
    momenta(:,3)=[2d0*energy,0d0,0d0,0d0]

    call amplitude%init(1,n,1,process,spin,order,model)
    call amplitude%evaluate(n,momenta,helicity,.false.,model)
    result=sum(abs(amplitude%amps)**2)
    yukawa=mass/(2d0*model%get_mass(24)*sw)
    ! AmpliCol's colour-ordered external spinors carry the corresponding
    ! 1/sqrt(2) normalization per non-zero helicity pair.
    expected=yukawa**2*(model%get_mass(25)**2-4d0*mass**2)
    tolerance=2d-11*max(1d0,expected)
    if (abs(result-expected).gt.tolerance) then
       write (*,*) 'Higgs-Yukawa matrix element mismatch',scheme,pdg
       write (*,*) 'computed, expected:',result,expected
       stop 1
    endif
  end subroutine evaluate_yukawa_process

  subroutine evaluate_massive_bottom_top_decay(result)
    implicit none
    real(kind=dp),intent(out) :: result
    integer,parameter :: n=4
    type(physics_model) :: model
    type(amplitude_QCD) :: amplitude
    integer,dimension(n,1) :: process,order
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: helicity
    real(kind=dp),dimension(0:3,n) :: momenta
    real(kind=dp) :: mb,mw,mt,sqrts_local,qin,qout,eb,eg,et,ew
    real(kind=dp),parameter :: cosine=0.3d0

    call set_flavour_scheme(4)
    call model%init_part()
    call model%init_vert()
    mb=model%get_mass(5)
    mw=model%get_mass(24)
    mt=model%get_mass(6)
    sqrts_local=1000d0
    qin=(sqrts_local**2-mb**2)/(2d0*sqrts_local)
    eb=(sqrts_local**2+mb**2)/(2d0*sqrts_local)
    eg=qin
    qout=sqrt((sqrts_local**2-(mt+mw)**2)*&
         (sqrts_local**2-(mt-mw)**2))/(2d0*sqrts_local)
    et=(sqrts_local**2+mt**2-mw**2)/(2d0*sqrts_local)
    ew=(sqrts_local**2+mw**2-mt**2)/(2d0*sqrts_local)

    ! b g -> t W- is a crossed t -> b W+ current.  In FS4 both fermion
    ! wavefunctions are massive, unlike the existing five-flavour check.
    process(:,1)=[5,21,6,-24]
    order(:,1)=[3,2,4,1]
    spin=0
    spin(0,:)=[2,2,2,3]
    spin(1:2,1)=[-1,1]
    spin(1:2,2)=[-1,1]
    spin(1:2,3)=[-1,1]
    spin(1:3,4)=[-1,0,1]
    helicity=0
    momenta(:,1)=[eb,0d0,0d0,qin]
    momenta(:,2)=[eg,0d0,0d0,-qin]
    momenta(:,3)=[et,qout*sqrt(1d0-cosine**2),0d0,qout*cosine]
    momenta(:,4)=[ew,-qout*sqrt(1d0-cosine**2),0d0,-qout*cosine]

    call amplitude%init(1,n,1,process,spin,order,model)
    call amplitude%evaluate(n,momenta,helicity,.false.,model)
    result=sum(abs(amplitude%amps)**2)
    if (.not.ieee_is_finite(result) .or. result.le.1d-12) then
       write (*,*) 'Massive-bottom t-b-W amplitude vanished in FS4',result
       stop 1
    endif
  end subroutine evaluate_massive_bottom_top_decay

end program flavour_scheme_yukawa_regression
