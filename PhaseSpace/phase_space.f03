module phase_space_base
  implicit none
  ! define an abstract base for the phase_space types
  type,abstract :: phase_space_type
     real(kind=8),dimension(:,:),allocatable,public :: p
     real(kind=8),public :: jac,xbjrk(2)
     real(kind=8) :: s0,tot_mass,sqrtshat,sqrts,drcut,ptcut,ycut
     real(kind=8),dimension(:),allocatable :: masses
     integer(kind=4) :: ndim,next
     integer(kind=4),dimension(:,:),allocatable :: sets
     real(kind=8),dimension(:,:),allocatable :: pp
     real(kind=8),dimension(:),allocatable :: x,invm,invm_min,invm_max,ETmin
     integer(kind=4),dimension(:),allocatable :: order
     logical :: t_channel
   contains
     procedure(phase_space_interface_init),deferred :: init
     procedure(phase_space_interface_generate_momenta),deferred :: generate_momenta
  end type phase_space_type
  
  ! Declare the abstract interface for the procedures
  abstract interface
     subroutine phase_space_interface_init(this,sqrts,n,m,o,s_cut,pt_cut,rap_cut,dr_cut,sqrt_s_min,t_chan,include_pdf,part)
       import :: phase_space_type
       class(phase_space_type),intent(inout) :: this
       real(kind=8),intent(in) :: sqrts
       integer(kind=4),intent(in) :: n
       integer(kind=4),dimension(n),intent(in) :: o,part
       real(kind=8),intent(in) :: s_cut(2),pt_cut,dr_cut,rap_cut,sqrt_s_min
       real(kind=8),dimension(n),intent(in) :: m
       logical,intent(in) :: t_chan
       logical,intent(in) :: include_pdf
     end subroutine phase_space_interface_init
     subroutine phase_space_interface_generate_momenta(this,xx)
       import :: phase_space_type
       class(phase_space_type),intent(inout) :: this
       real(kind=8),dimension(99),intent(in) :: xx
     end subroutine phase_space_interface_generate_momenta
  end interface
end module phase_space_base
