module phase_space_base
  implicit none
  type :: psv
     real(kind=8),dimension(:,:),allocatable,public :: p
     real(kind=8),dimension(:),allocatable,public :: x
     real(kind=8),public :: jac,xbjrk(2)
   contains
     final :: finalize_psv
  end type psv
  type,abstract :: phase_space_type
     ! If adding variables here, make sure to also update the
     ! 'cleanup' subroutines for all phase-space parametrisation
     ! module.
     real(kind=8),dimension(:,:),allocatable,public :: p
     real(kind=8),public :: jac,xbjrk(2)
     real(kind=8) :: s0,tot_mass,sqrtshat,sqrts
     real(kind=8),dimension(:),allocatable :: masses,ptcut,ycut,drcut
     integer(kind=4) :: ndim,next,ndim_extra
     integer(kind=4),dimension(:,:),allocatable :: sets
     real(kind=8),dimension(:,:),allocatable :: pp,sqrt_s_min
     real(kind=8),dimension(:),allocatable :: x,invm,invm_min,invm_max,ETmin
     integer(kind=4),dimension(:),allocatable :: order
     integer(kind=4),dimension(:),allocatable :: resonance_pdgs,resonance_masks
     real(kind=8),dimension(:),allocatable :: resonance_masses,resonance_widths
     integer(kind=4) :: nresonances=0
     logical :: t_channel
     logical :: can_invert_momenta=.false.
   contains
     procedure(phase_space_interface_init),deferred :: init
     procedure(phase_space_interface_generate_momenta),deferred :: generate_momenta
     procedure(phase_space_interface_compute_x_from_momenta),deferred :: compute_x_from_momenta
     procedure(phase_space_interface_cleanup),deferred :: cleanup
     procedure :: configure_resonances => phase_space_configure_resonances
  end type phase_space_type
  
  ! Declare the abstract interface for the procedures
  abstract interface
     subroutine phase_space_interface_init(this,sqrts,n,m,o,pt_cut,rap_cut,dr_cut,sqrt_s_min,t_chan,include_pdf,flat)
       import :: phase_space_type
       class(phase_space_type),intent(inout) :: this
       real(kind=8),intent(in) :: sqrts
       integer(kind=4),intent(in) :: n
       integer(kind=4),dimension(n),intent(in) :: o
       real(kind=8),dimension(n,n),intent(in) :: dr_cut,sqrt_s_min
       real(kind=8),dimension(n),intent(in) :: m,pt_cut,rap_cut
       logical,intent(in) :: t_chan
       logical,intent(in) :: include_pdf
       logical,intent(in),optional :: flat
     end subroutine phase_space_interface_init
     subroutine phase_space_interface_generate_momenta(this,ps)
       import :: phase_space_type,psv
       class(phase_space_type),intent(inout) :: this
       type(psv),intent(inout) :: ps
     end subroutine phase_space_interface_generate_momenta
     subroutine phase_space_interface_compute_x_from_momenta(this,ps)
       import :: phase_space_type,psv
       class(phase_space_type),intent(inout) :: this
       type(psv),intent(inout) :: ps
     end subroutine phase_space_interface_compute_x_from_momenta
     subroutine phase_space_interface_cleanup(this)
       import :: phase_space_type
       class(phase_space_type),intent(inout) :: this
     end subroutine phase_space_interface_cleanup
  end interface
contains
  subroutine phase_space_configure_resonances(this,nresonances,pdgs,masks,masses,widths)
    class(phase_space_type),intent(inout) :: this
    integer(kind=4),intent(in) :: nresonances
    integer(kind=4),dimension(nresonances),intent(in) :: pdgs,masks
    real(kind=8),dimension(nresonances),intent(in) :: masses,widths
    if (allocated(this%resonance_pdgs)) deallocate(this%resonance_pdgs)
    if (allocated(this%resonance_masks)) deallocate(this%resonance_masks)
    if (allocated(this%resonance_masses)) deallocate(this%resonance_masses)
    if (allocated(this%resonance_widths)) deallocate(this%resonance_widths)
    this%nresonances=nresonances
    if (nresonances.eq.0) return
    if (any(masses.le.0d0) .or. any(widths.le.0d0)) then
       write (*,*) 'A mapped resonance must have positive mass and width',pdgs,masses,widths
       stop 1
    endif
    allocate(this%resonance_pdgs(nresonances))
    allocate(this%resonance_masks(nresonances))
    allocate(this%resonance_masses(nresonances))
    allocate(this%resonance_widths(nresonances))
    this%resonance_pdgs=pdgs
    this%resonance_masks=masks
    this%resonance_masses=masses
    this%resonance_widths=widths
  end subroutine phase_space_configure_resonances

  subroutine finalize_psv(this)
    type(psv),intent(inout) :: this
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%x)) deallocate(this%x)
  end subroutine finalize_psv
end module phase_space_base
