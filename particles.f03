module particles
  implicit none
  type particle
     integer :: type,anti_type,spin
     real(kind=8) :: mass,width
  end type particle
  type physics_model
     type(particle),dimension(:),allocatable :: particle_list
     integer :: npart
   contains
     procedure,public :: init_part,get_mass,get_width,get_spin,get_antipart
  end type physics_model
contains
  subroutine init_part(this,tmass,twidth)
    implicit none
    class(physics_model) :: this
    integer :: i
    real(kind=8) :: tmass,twidth
    this%npart=9 ! gluon, 6 quarks, tensor, and the photon
    allocate(this%particle_list(this%npart))

    ! 5 massless quarks
    do i=1,5
       this%particle_list(i)%type=i
       this%particle_list(i)%mass=0d0
       this%particle_list(i)%width=0d0
       this%particle_list(i)%spin=2 ! two spin states
       this%particle_list(i)%anti_type=-i
    enddo

    ! top quark
    this%particle_list(6)%type=6
    this%particle_list(6)%mass=tmass
    this%particle_list(6)%width=twidth
    this%particle_list(6)%spin=2 ! two spin states
    this%particle_list(6)%anti_type=-6

    ! gluon
    this%particle_list(7)%type=21
    this%particle_list(7)%mass=0d0
    this%particle_list(7)%width=0d0
    this%particle_list(7)%spin=2 ! two spin states
    this%particle_list(7)%anti_type=21

    ! tensor (non-propagator auxiliary particle to decompose 4-gluon interaction)
    this%particle_list(8)%type=-21
    this%particle_list(8)%mass=0d0
    this%particle_list(8)%width=0d0
    this%particle_list(8)%spin=-1 ! ill-defined
    this%particle_list(8)%anti_type=-21

    ! photon
    this%particle_list(9)%type=22
    this%particle_list(9)%mass=0d0
    this%particle_list(9)%width=0d0
    this%particle_list(9)%spin=2 ! two spin states
    this%particle_list(9)%anti_type=22
    
  end subroutine init_part
  integer function get_antipart(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    do i=1,this%npart
       if (this%particle_list(i)%type.eq.ipdg) then
          get_antipart=this%particle_list(i)%anti_type
          return
       elseif (this%particle_list(i)%anti_type.eq.ipdg) then
          get_antipart=this%particle_list(i)%type
          return
       endif
    enddo
    write (*,*) 'Particle not in model (mass)',ipdg
    stop 1
  end function get_antipart
  real(kind=8) function get_mass(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    do i=1,this%npart
       if (this%particle_list(i)%type.eq.ipdg .or. this%particle_list(i)%anti_type.eq.ipdg) then
          get_mass=this%particle_list(i)%mass
          return
       endif
    enddo
    write (*,*) 'Particle not in model (mass)',ipdg
    stop 1
  end function get_mass
  real(kind=8) function get_width(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    do i=1,this%npart
       if (this%particle_list(i)%type.eq.ipdg .or. this%particle_list(i)%anti_type.eq.ipdg) then
          get_width=this%particle_list(i)%width
          return
       endif
    enddo
    write (*,*) 'Particle not in model (width)',ipdg
    stop 1
  end function get_width
  integer function get_spin(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    do i=1,this%npart
       if (this%particle_list(i)%type.eq.ipdg .or. this%particle_list(i)%anti_type.eq.ipdg) then
          get_spin=this%particle_list(i)%spin
          if (get_spin.lt.0) then
             write (*,*) 'Spin ill-defined for particle',ipdg
             stop 1
          endif
          return
       endif
    enddo
    write (*,*) 'Particle not in model (spin)',ipdg
    stop 1
  end function get_spin
  logical function is_quark(iPDG)
    integer :: iPDG
    if (iPDG.ge.1 .and. iPDG.le.6) then
       is_quark=.true.
    else
       is_quark=.false.
    endif
  end function is_quark
  logical function is_antiquark(iPDG)
    integer :: iPDG
    if (iPDG.le.-1 .and. iPDG.ge.-6) then
       is_antiquark=.true.
    else
       is_antiquark=.false.
    endif
  end function is_antiquark
  logical function is_gluon(iPDG)
    integer :: iPDG
    if (iPDG.eq.21) then
       is_gluon=.true.
    else
       is_gluon=.false.
    endif
  end function is_gluon
  logical function is_tensor(i)
    implicit none
    integer :: i
    if (i.eq.-21) then
       is_tensor=.true.
    else
       is_tensor=.false.
    endif
  end function is_tensor
  logical function is_singlet(i)
    implicit none
    integer :: i
    if (abs(i).ge.22) then
       is_singlet=.true.
    else
       is_singlet=.false.
    endif
  end function is_singlet
end module particles
