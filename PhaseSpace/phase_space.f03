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
     real(kind=8),dimension(:),allocatable :: x,invm
     real(kind=8),dimension(:,:),allocatable :: invm_min,invm_max,ETmin
     integer(kind=4),dimension(:),allocatable :: order
     logical :: t_channel
   contains
     procedure(phase_space_interface_init),deferred :: init
     procedure(phase_space_interface_generate_momenta),deferred :: generate_momenta
     procedure(phase_space_interface_compute_x_from_momenta),deferred :: compute_x_from_momenta
     procedure(phase_space_interface_cleanup),deferred :: cleanup
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
  subroutine finalize_psv(this)
    type(psv),intent(inout) :: this
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%x)) deallocate(this%x)
  end subroutine finalize_psv
end module phase_space_base
