module run_parameters
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(kind=8),parameter,private :: pi=3.1415926535897932384626433832795d0
  real(kind=8),parameter,private :: default_z_mass=91.188d0
  real(kind=8),parameter,private :: default_w_mass=80.419002445756163d0
  real(kind=8),parameter,private :: default_alphaEW=0.007546771114d0
  real(kind=8),parameter,private :: default_sw=&
       sqrt(1d0-(default_w_mass/default_z_mass)**2)
  real(kind=8),parameter,private :: default_heft_vev=&
       2d0*default_w_mass*default_sw/sqrt(4d0*pi*default_alphaEW)

  ! Particle properties used by the built-in Standard Model.
  real(kind=8) :: top_mass=173d0
  real(kind=8) :: top_width=1.4915d0
  real(kind=8) :: z_mass=default_z_mass
  real(kind=8) :: z_width=2.441404d0
  real(kind=8) :: w_mass=default_w_mass
  real(kind=8) :: w_width=2.0476d0
  real(kind=8) :: higgs_mass=125d0
  real(kind=8) :: higgs_width=0.0063823389999999999d0

  ! Couplings.  The on-shell weak-mixing sine and electroweak vev are derived
  ! from mW, mZ and alphaEW whenever an input card is read.  PROTECTED keeps
  ! callers from making these dependent quantities inconsistent with them.
  real(kind=8),protected :: sw=default_sw
  real(kind=8) :: alphaS_MZ=0.119d0
  real(kind=8) :: alphaEW=default_alphaEW
  real(kind=8) :: heft_kappa=1d0
  real(kind=8),protected :: heft_vev=default_heft_vev

  private :: update_derived_electroweak_parameters

  ! Collider, scale and PDF setup.
  real(kind=8) :: sqrts=14000d0
  integer :: scale_choice=2
  logical :: include_pdf=.true.
  logical :: use_lhapdf=.true.
  character(len=128) :: lhapdfset='NNPDF23_nlo_as_0119_qed'
  integer :: lhapdf_member=0
  integer :: pdf_lhaid=244800
  character(len=256) :: internal_pdf_grid='NNPDF23nlo_as_0119_qed_mem0.grid'

  ! By default an unstable particle requested in the physical final state is
  ! treated as on shell throughout every subprocess.  Set this switch to true
  ! only to retain its configured width.
  logical :: ignore_final_state_width_fix=.false.

  ! Cuts.  A negative value disables the corresponding cut.
  real(kind=8) :: pTj_min=30d0
  real(kind=8) :: DRjj_min=0.4d0
  real(kind=8) :: etaj_max=6d0
  real(kind=8) :: sqrt_sjj_min=-1d0

  real(kind=8) :: pTa_min=30d0
  real(kind=8) :: DRaa_min=0.4d0
  real(kind=8) :: etaa_max=6d0
  real(kind=8) :: sqrt_saa_min=-1d0

  real(kind=8) :: pTl_min=20d0
  real(kind=8) :: DRll_min=0.4d0
  real(kind=8) :: etal_max=2.5d0
  real(kind=8) :: sqrt_sll_min=-1d0

  real(kind=8) :: DRja_min=0.4d0
  real(kind=8) :: sqrt_sja_min=-1d0
  real(kind=8) :: DRjl_min=0.4d0
  real(kind=8) :: sqrt_sjl_min=-1d0
  real(kind=8) :: DRla_min=0.4d0
  real(kind=8) :: sqrt_sla_min=-1d0

  character(len=256) :: active_run_card=''

  namelist /amplicol/ top_mass,top_width,z_mass,z_width,w_mass,w_width,&
       higgs_mass,higgs_width,alphaS_MZ,alphaEW,heft_kappa,&
       sqrts,scale_choice,&
       include_pdf,use_lhapdf,lhapdfset,lhapdf_member,pdf_lhaid,&
       internal_pdf_grid,ignore_final_state_width_fix,&
       pTj_min,DRjj_min,etaj_max,sqrt_sjj_min,&
       pTa_min,DRaa_min,etaa_max,sqrt_saa_min,&
       pTl_min,DRll_min,etal_max,sqrt_sll_min,&
       DRja_min,sqrt_sja_min,DRjl_min,sqrt_sjl_min,DRla_min,sqrt_sla_min

contains

  real(kind=8) function heft_coupling(alpha_s)
    implicit none
    real(kind=8),intent(in) :: alpha_s
    heft_coupling=heft_kappa*alpha_s/(3d0*pi*heft_vev)
  end function heft_coupling

  subroutine update_derived_electroweak_parameters()
    implicit none
    sw=sqrt(1d0-(w_mass/z_mass)**2)
    heft_vev=2d0*w_mass*sw/sqrt(4d0*pi*alphaEW)
  end subroutine update_derived_electroweak_parameters

  subroutine reset_run_parameters()
    implicit none
    top_mass=173d0
    top_width=1.4915d0
    z_mass=default_z_mass
    z_width=2.441404d0
    w_mass=default_w_mass
    w_width=2.0476d0
    higgs_mass=125d0
    higgs_width=0.0063823389999999999d0
    alphaS_MZ=0.119d0
    alphaEW=default_alphaEW
    heft_kappa=1d0
    sqrts=14000d0
    scale_choice=2
    include_pdf=.true.
    use_lhapdf=.true.
    lhapdfset='NNPDF23_nlo_as_0119_qed'
    lhapdf_member=0
    pdf_lhaid=244800
    internal_pdf_grid='NNPDF23nlo_as_0119_qed_mem0.grid'
    ignore_final_state_width_fix=.false.
    pTj_min=30d0
    DRjj_min=0.4d0
    etaj_max=6d0
    sqrt_sjj_min=-1d0
    pTa_min=30d0
    DRaa_min=0.4d0
    etaa_max=6d0
    sqrt_saa_min=-1d0
    pTl_min=20d0
    DRll_min=0.4d0
    etal_max=2.5d0
    sqrt_sll_min=-1d0
    DRja_min=0.4d0
    sqrt_sja_min=-1d0
    DRjl_min=0.4d0
    sqrt_sjl_min=-1d0
    DRla_min=0.4d0
    sqrt_sla_min=-1d0
    active_run_card=''
    call update_derived_electroweak_parameters()
  end subroutine reset_run_parameters

  subroutine read_run_parameters(filename)
    implicit none
    character(len=*),intent(in) :: filename
    integer :: iunit,ios
    character(len=512) :: message

    call reset_run_parameters()
    open(newunit=iunit,file=trim(filename),status='old',action='read',&
         iostat=ios,iomsg=message)
    if (ios.ne.0) then
       write (*,*) 'Could not open AmpliCol input file ',trim(filename)
       write (*,*) trim(message)
       stop 1
    endif
    read(iunit,nml=amplicol,iostat=ios,iomsg=message)
    close(iunit)
    if (ios.ne.0) then
       write (*,*) 'Could not read &amplicol namelist from ',trim(filename)
       write (*,*) trim(message)
       stop 1
    endif
    active_run_card=trim(filename)
    call validate_run_parameters()
  end subroutine read_run_parameters

  subroutine validate_run_parameters()
    implicit none
    if (any(.not.ieee_is_finite([top_mass,top_width,z_mass,z_width,&
         w_mass,w_width,higgs_mass,higgs_width,alphaS_MZ,alphaEW,&
         heft_kappa,sqrts,&
         pTj_min,DRjj_min,etaj_max,sqrt_sjj_min,pTa_min,DRaa_min,etaa_max,&
         sqrt_saa_min,pTl_min,DRll_min,etal_max,sqrt_sll_min,DRja_min,&
         sqrt_sja_min,DRjl_min,sqrt_sjl_min,DRla_min,sqrt_sla_min]))) then
       write (*,*) 'All real-valued input parameters must be finite'
       stop 1
    endif
    if (top_mass.lt.0d0 .or. z_mass.le.0d0 .or. w_mass.le.0d0 .or.&
         higgs_mass.le.0d0) then
       write (*,*) 'Particle masses in the input file must be non-negative (and boson masses positive)'
       stop 1
    endif
    if (top_width.lt.0d0 .or. z_width.lt.0d0 .or. w_width.lt.0d0 .or.&
         higgs_width.lt.0d0) then
       write (*,*) 'Particle widths in the input file must be non-negative'
       stop 1
    endif
    if (w_mass.ge.z_mass) then
       write (*,*) 'w_mass must be smaller than z_mass to derive sin(theta_W)',&
            w_mass,z_mass
       stop 1
    endif
    if (alphaS_MZ.le.0d0 .or. alphaEW.le.0d0) then
       write (*,*) 'alphaS_MZ and alphaEW must be positive',alphaS_MZ,alphaEW
       stop 1
    endif
    call update_derived_electroweak_parameters()
    if (.not.ieee_is_finite(sw) .or. .not.ieee_is_finite(heft_vev) .or.&
         sw.le.0d0 .or. sw.ge.1d0 .or. heft_vev.le.0d0) then
       write (*,*) 'Could not derive finite positive electroweak parameters',&
            sw,heft_vev
       stop 1
    endif
    if (sqrts.le.0d0) then
       write (*,*) 'sqrts must be positive',sqrts
       stop 1
    endif
    if (scale_choice.lt.0 .or. scale_choice.gt.5) then
       write (*,*) 'scale_choice must be between 0 and 5',scale_choice
       stop 1
    endif
    if (lhapdf_member.lt.0) then
       write (*,*) 'lhapdf_member must be non-negative',lhapdf_member
       stop 1
    endif
    if (len_trim(lhapdfset).eq.0 .and. use_lhapdf) then
       write (*,*) 'lhapdfset must not be empty when use_lhapdf is true'
       stop 1
    endif
    if (len_trim(internal_pdf_grid).eq.0 .and. include_pdf .and. (.not.use_lhapdf)) then
       write (*,*) 'internal_pdf_grid must not be empty when the internal PDF is used'
       stop 1
    endif
  end subroutine validate_run_parameters

  subroutine write_run_parameters(iunit)
    implicit none
    integer,intent(in) :: iunit
    write(iunit,'(a)') 'AmpliCol input parameters from '//trim(active_run_card)
    write(iunit,nml=amplicol)
  end subroutine write_run_parameters

end module run_parameters
