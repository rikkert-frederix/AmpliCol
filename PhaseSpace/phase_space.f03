module phase_space_base
  implicit none
  integer,parameter,public :: transform_breit_wigner=1
  integer,parameter,public :: transform_massless_pole=2
  integer,parameter,public :: transform_massive_power=3
  integer,parameter,public :: transform_flat_contact=4
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
     integer(kind=4),dimension(:),allocatable :: topology_pdgs,topology_masks,&
          topology_left_masks,topology_kinds,topology_parameters
     real(kind=8),dimension(:),allocatable :: topology_masses,topology_widths
     integer(kind=4) :: ntopology_nodes=0
     logical :: t_channel
     logical :: can_invert_momenta=.false.
   contains
     procedure(phase_space_interface_init),deferred :: init
     procedure(phase_space_interface_generate_momenta),deferred :: generate_momenta
     procedure(phase_space_interface_compute_x_from_momenta),deferred :: compute_x_from_momenta
     procedure(phase_space_interface_cleanup),deferred :: cleanup
     procedure :: configure_resonances => phase_space_configure_resonances
     procedure :: configure_topology => phase_space_configure_topology
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
    integer :: i
    if (any(masses.le.0d0) .or. any(widths.le.0d0)) then
       write (*,*) 'A mapped resonance must have positive mass and width',pdgs,masses,widths
       stop 1
    endif
    call this%configure_topology(nresonances,pdgs,masks,&
         [(0, i=1,nresonances)],masses,widths)
  end subroutine phase_space_configure_resonances

  subroutine phase_space_configure_topology(this,nnodes,pdgs,masks,left_masks,&
       masses,widths,kinds,parameters)
    class(phase_space_type),intent(inout) :: this
    integer(kind=4),intent(in) :: nnodes
    integer(kind=4),dimension(nnodes),intent(in) :: pdgs,masks,left_masks
    real(kind=8),dimension(nnodes),intent(in) :: masses,widths
    integer(kind=4),dimension(nnodes),intent(in),optional :: kinds,parameters
    integer :: i
    if (allocated(this%topology_pdgs)) deallocate(this%topology_pdgs)
    if (allocated(this%topology_masks)) deallocate(this%topology_masks)
    if (allocated(this%topology_left_masks)) deallocate(this%topology_left_masks)
    if (allocated(this%topology_kinds)) deallocate(this%topology_kinds)
    if (allocated(this%topology_parameters)) deallocate(this%topology_parameters)
    if (allocated(this%topology_masses)) deallocate(this%topology_masses)
    if (allocated(this%topology_widths)) deallocate(this%topology_widths)
    this%ntopology_nodes=nnodes
    if (nnodes.eq.0) return
    if (any(masses.lt.0d0) .or. any(widths.lt.0d0)) then
       write (*,*) 'A topology node cannot have a negative mass or width',&
            pdgs,masses,widths
       stop 1
    endif
    allocate(this%topology_pdgs(nnodes))
    allocate(this%topology_masks(nnodes))
    allocate(this%topology_left_masks(nnodes))
    allocate(this%topology_kinds(nnodes))
    allocate(this%topology_parameters(nnodes))
    allocate(this%topology_masses(nnodes))
    allocate(this%topology_widths(nnodes))
    this%topology_pdgs=pdgs
    this%topology_masks=masks
    this%topology_left_masks=left_masks
    if (present(kinds)) then
       if (any(kinds.lt.transform_breit_wigner) .or.&
            any(kinds.gt.transform_flat_contact)) then
          write (*,*) 'Unknown topology transform kind',kinds
          stop 1
       endif
       this%topology_kinds=kinds
    else
       do i=1,nnodes
          if (pdgs(i).eq.0 .or. pdgs(i).eq.-21 .or. pdgs(i).eq.-23 .or.&
               abs(pdgs(i)).eq.26 .or. (pdgs(i).ge.125 .and. pdgs(i).le.127)) then
             this%topology_kinds(i)=transform_flat_contact
          elseif (masses(i).gt.0d0 .and. widths(i).gt.0d0) then
             this%topology_kinds(i)=transform_breit_wigner
          elseif (masses(i).eq.0d0) then
             this%topology_kinds(i)=transform_massless_pole
          else
             this%topology_kinds(i)=transform_massive_power
          endif
       enddo
    endif
    if (present(parameters)) then
       this%topology_parameters=parameters
    else
       this%topology_parameters=abs(pdgs)
    endif
    this%topology_masses=masses
    this%topology_widths=widths
  end subroutine phase_space_configure_topology

  subroutine finalize_psv(this)
    type(psv),intent(inout) :: this
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%x)) deallocate(this%x)
  end subroutine finalize_psv
end module phase_space_base
