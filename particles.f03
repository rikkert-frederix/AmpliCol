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
  subroutine init_part(this,tmass,twidth,zmass,zwidth,wmass,wwidth,hmass,hwidth)
    implicit none
    class(physics_model) :: this
    integer :: i
    real(kind=8) :: tmass,twidth
    real(kind=8) :: zmass,zwidth
    real(kind=8) :: wmass,wwidth
    real(kind=8) :: hmass,hwidth

    this%npart=23 ! gluon, 6 quarks, tensor, photon, higgs, z-boson, w-boson + 10 tensors for EW bosons
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

    ! z-boson
    this%particle_list(10)%type=23
    this%particle_list(10)%mass=zmass
    this%particle_list(10)%width=zwidth
    this%particle_list(10)%spin=3 ! three spin states
    this%particle_list(10)%anti_type=23

    ! w-boson
    this%particle_list(11)%type=24
    this%particle_list(11)%mass=wmass
    this%particle_list(11)%width=wwidth
    this%particle_list(11)%spin=3 ! three spin states
    this%particle_list(11)%anti_type=-24

    this%particle_list(12)%type=-24
    this%particle_list(12)%mass=wmass
    this%particle_list(12)%width=wwidth
    this%particle_list(12)%spin=3 ! three spin states
    this%particle_list(12)%anti_type=24

    ! w+w+ tensor T1
    this%particle_list(13)%type=-101
    this%particle_list(13)%mass=0d0
    this%particle_list(13)%width=0d0
    this%particle_list(13)%spin=-1 
    this%particle_list(13)%anti_type=-101

    ! w+w- tensor T2
    this%particle_list(14)%type=-102
    this%particle_list(14)%mass=0d0
    this%particle_list(14)%width=0d0
    this%particle_list(14)%spin=-1
    this%particle_list(14)%anti_type=-102

    ! w+z tensor T3
    this%particle_list(15)%type=-103
    this%particle_list(15)%mass=0d0
    this%particle_list(15)%width=0d0
    this%particle_list(15)%spin=-1
    this%particle_list(15)%anti_type=-103

    ! w+a tensor T4
    this%particle_list(16)%type=-104
    this%particle_list(16)%mass=0d0
    this%particle_list(16)%width=0d0
    this%particle_list(16)%spin=-1
    this%particle_list(16)%anti_type=-104

    ! w-w- tensor T5
    this%particle_list(17)%type=-105
    this%particle_list(17)%mass=0d0
    this%particle_list(17)%width=0d0
    this%particle_list(17)%spin=-1
    this%particle_list(17)%anti_type=-105

    ! w-z tensor T6
    this%particle_list(18)%type=-106
    this%particle_list(18)%mass=0d0
    this%particle_list(18)%width=0d0
    this%particle_list(18)%spin=-1
    this%particle_list(18)%anti_type=-106

    ! w-a tensor T7
    this%particle_list(19)%type=-107
    this%particle_list(19)%mass=0d0
    this%particle_list(19)%width=0d0
    this%particle_list(19)%spin=-1
    this%particle_list(19)%anti_type=-107

    ! zz tensor T8
    this%particle_list(20)%type=-108
    this%particle_list(20)%mass=0d0
    this%particle_list(20)%width=0d0
    this%particle_list(20)%spin=-1
    this%particle_list(20)%anti_type=-108

    ! za tensor T9
    this%particle_list(21)%type=-109
    this%particle_list(21)%mass=0d0
    this%particle_list(21)%width=0d0
    this%particle_list(21)%spin=-1
    this%particle_list(21)%anti_type=-109

    ! aa tensor T10
    this%particle_list(22)%type=-110
    this%particle_list(22)%mass=0d0
    this%particle_list(22)%width=0d0
    this%particle_list(22)%spin=-1
    this%particle_list(22)%anti_type=-110

    ! higgs
    this%particle_list(23)%type=25
    this%particle_list(23)%mass=hmass
    this%particle_list(23)%width=hwidth
    this%particle_list(23)%spin=1
    this%particle_list(23)%anti_type=25
    
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
  logical function is_scalar(i)
    implicit none
    integer :: i
    if (i.eq.25) then
       is_scalar=.true.
    else
       is_scalar=.false.
    endif
  end function is_scalar
  logical function is_tensor_v(i)
    implicit none
    integer :: i
    if (is_tensor_t1(i).or.&
        is_tensor_t2(i).or.&
        is_tensor_t3(i).or.&
        is_tensor_t4(i).or.&
        is_tensor_t5(i).or.&
        is_tensor_t6(i).or.&
        is_tensor_t7(i).or.&
        is_tensor_t8(i).or.&
        is_tensor_t9(i).or.&
        is_tensor_t10(i)) then
       is_tensor_v=.true.
    else
       is_tensor_v=.false.
    endif
  end function is_tensor_v
  logical function is_tensor_t1(i)
    implicit none
    integer :: i
    if (i.eq.-101) then
       is_tensor_t1=.true.
    else
       is_tensor_t1=.false.
    endif
  end function is_tensor_t1
  logical function is_tensor_t2(i)
    implicit none
    integer :: i
    if (i.eq.-102) then
       is_tensor_t2=.true.
    else
       is_tensor_t2=.false.
    endif
  end function is_tensor_t2
  logical function is_tensor_t3(i)
    implicit none
    integer :: i
    if (i.eq.-103) then
       is_tensor_t3=.true.
    else
       is_tensor_t3=.false.
    endif
  end function is_tensor_t3
  logical function is_tensor_t4(i)
    implicit none
    integer :: i
    if (i.eq.-104) then
       is_tensor_t4=.true.
    else
       is_tensor_t4=.false.
    endif
  end function is_tensor_t4
  logical function is_tensor_t5(i)
    implicit none
    integer :: i
    if (i.eq.-105) then
       is_tensor_t5=.true.
    else
       is_tensor_t5=.false.
    endif
  end function is_tensor_t5
  logical function is_tensor_t6(i)
    implicit none
    integer :: i
    if (i.eq.-106) then
       is_tensor_t6=.true.
    else
       is_tensor_t6=.false.
    endif
  end function is_tensor_t6
  logical function is_tensor_t7(i)
    implicit none
    integer :: i
    if (i.eq.-107) then
       is_tensor_t7=.true.
    else
       is_tensor_t7=.false.
    endif
  end function is_tensor_t7
  logical function is_tensor_t8(i)
    implicit none
    integer :: i
    if (i.eq.-108) then
       is_tensor_t8=.true.
    else
       is_tensor_t8=.false.
    endif
  end function is_tensor_t8
  logical function is_tensor_t9(i)
    implicit none
    integer :: i
    if (i.eq.-109) then
       is_tensor_t9=.true.
    else
       is_tensor_t9=.false.
    endif
  end function is_tensor_t9
  logical function is_tensor_t10(i)
    implicit none
    integer :: i
    if (i.eq.-110) then
       is_tensor_t10=.true.
    else
       is_tensor_t10=.false.
    endif
  end function is_tensor_t10
  logical function is_singlet(i)
    implicit none
    integer :: i
    if (abs(i).ge.22.and.abs(i).le.25) then
       is_singlet=.true.
    else
       is_singlet=.false.
    endif
  end function is_singlet
  logical function is_singlet_a(i)
    implicit none
    integer :: i
    if (abs(i).eq.22) then
       is_singlet_a=.true.
    else
       is_singlet_a=.false.
    endif
  end function is_singlet_a
  logical function is_singlet_z(i)
    implicit none
    integer :: i
    if (abs(i).eq.23) then
       is_singlet_z=.true.
    else
       is_singlet_z=.false.
    endif
  end function is_singlet_z
  logical function is_singlet_w(i)
    implicit none
    integer :: i
    if (abs(i).eq.24) then
       is_singlet_w=.true.
    else
       is_singlet_w=.false.
    endif
  end function is_singlet_w
end module particles
