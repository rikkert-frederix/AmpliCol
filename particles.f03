module particles
  implicit none
  real(kind=8),parameter :: sw = 0.47143025548407230d0 
  type particle
     integer :: type,anti_type,spin,dim
     real(kind=8) :: mass,width
  end type particle
  type vertex
     integer :: type
     integer,dimension(3) :: particles
     real(kind=8),dimension(2) :: coupl
  end type vertex
  type physics_model
     type(particle),dimension(:),allocatable :: particle_list
     type(vertex),dimension(:),allocatable :: vertex_list
     integer :: npart,nint
   contains
     procedure,public :: init_part,get_mass,get_width,get_spin,get_antipart,init_vert,get_dim,get_inter_dim
  end type physics_model
contains
  subroutine init_part(this,tmass,twidth,zmass,zwidth,wmass,wwidth)
    implicit none
    class(physics_model) :: this
    integer :: i,l
    real(kind=8) :: tmass,twidth
    real(kind=8) :: zmass,zwidth
    real(kind=8) :: wmass,wwidth
    l=0
    this%npart=13! gluon, 6 quarks, tensor, photon, Z-boson and W-boson, etc.
    allocate(this%particle_list(this%npart))

    ! 5 massless quarks
    do i=1,5
       l=l+1
       this%particle_list(l)%type=i
       this%particle_list(l)%mass=0d0
       this%particle_list(l)%width=0d0
       this%particle_list(l)%spin=2 ! two spin states
       this%particle_list(l)%anti_type=-i
       this%particle_list(l)%dim=4
    enddo

    ! top quark
    l=l+1
    this%particle_list(l)%type=6
    this%particle_list(l)%mass=tmass
    this%particle_list(l)%width=twidth
    this%particle_list(l)%spin=2 ! two spin states
    this%particle_list(l)%anti_type=-6
    this%particle_list(l)%dim=4
       
    ! gluon
    l=l+1
    this%particle_list(l)%type=21
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=2 ! two spin states
    this%particle_list(l)%anti_type=21
    this%particle_list(l)%dim=4

    ! tensor (non-propagator auxiliary particle to decompose 4-gluon interaction)
    l=l+1
    this%particle_list(l)%type=-21
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=-1 ! ill-defined
    this%particle_list(l)%anti_type=-21
    this%particle_list(l)%dim=6

    ! photon
    l=l+1
    this%particle_list(l)%type=22
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=2 ! two spin states
    this%particle_list(l)%anti_type=22
    this%particle_list(l)%dim=4
    
    ! Z-boson
    l=l+1
    this%particle_list(l)%type=23
    this%particle_list(l)%mass=zmass
    this%particle_list(l)%width=zwidth
    this%particle_list(l)%spin=3 ! three spin states
    this%particle_list(l)%anti_type=23
    this%particle_list(l)%dim=4
    
    ! Ztensor (non-propagator auxiliary particle to decompose 4-Wboson interaction)
    l=l+1
    this%particle_list(l)%type=-23
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=-1 ! ill defined
    this%particle_list(l)%anti_type=-23
    this%particle_list(l)%dim=6
    
    ! W-boson
    l=l+1
    this%particle_list(l)%type=24
    this%particle_list(l)%mass=wmass
    this%particle_list(l)%width=wwidth
    this%particle_list(l)%spin=3 ! three spin states
    this%particle_list(l)%anti_type=-24
    this%particle_list(l)%dim=4

    ! Wtensor (non-propagator auxiliary particle to decompose 4-boson interactions)
    l=l+1
    this%particle_list(l)%type=26
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=-1 ! ill-defined
    this%particle_list(l)%anti_type=-26
    this%particle_list(l)%dim=6

    write (*,*) l,'particles loaded'
    
  end subroutine init_part

  subroutine init_vert(this)
    implicit none
    class(physics_model) :: this
    integer :: i,l
    real(kind=8) :: fact,gw,Vf,Af
    l=0
    this%nint = 118 ! number of vertices
    allocate(this%vertex_list(this%nint))
    ! gluon-gluon to gluon vertex
    l=l+1
    this%vertex_list(l)%type=0
    this%vertex_list(l)%particles(1)=21
    this%vertex_list(l)%particles(2)=21
    this%vertex_list(l)%particles(3)=21
    this%vertex_list(l)%coupl=[1d0,0d0]
    ! gluon-gluon to tensor vertex
    l=l+1
    this%vertex_list(l)%type=1
    this%vertex_list(l)%particles(1)=21
    this%vertex_list(l)%particles(2)=21
    this%vertex_list(l)%particles(3)=-21
    this%vertex_list(l)%coupl=[1d0,0d0]
    ! tensor-gluon to gluon vertex
    l=l+1
    this%vertex_list(l)%type=2
    this%vertex_list(l)%particles(1)=-21
    this%vertex_list(l)%particles(2)=21
    this%vertex_list(l)%particles(3)=21
    this%vertex_list(l)%coupl=[1d0,0d0]
    ! gluon-tensor to gluon vertex
    l=l+1
    this%vertex_list(l)%type=3
    this%vertex_list(l)%particles(1)=21
    this%vertex_list(l)%particles(2)=-21
    this%vertex_list(l)%particles(3)=21
    this%vertex_list(l)%coupl=[1d0,0d0]
    ! gluon-quark to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=4
       this%vertex_list(l)%particles(1)=21
       this%vertex_list(l)%particles(2)=i
       this%vertex_list(l)%particles(3)=i
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! gluon-antiquark to antiquark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=5
       this%vertex_list(l)%particles(1)=21
       this%vertex_list(l)%particles(2)=-i
       this%vertex_list(l)%particles(3)=-i
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! quark-gluon to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=6
       this%vertex_list(l)%particles(1)=i
       this%vertex_list(l)%particles(2)=21
       this%vertex_list(l)%particles(3)=i
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! antiquark-gluon to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=7
       this%vertex_list(l)%particles(1)=-i
       this%vertex_list(l)%particles(2)=21
       this%vertex_list(l)%particles(3)=-i
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! quark-antiquark to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=8
       this%vertex_list(l)%particles(1)=i
       this%vertex_list(l)%particles(2)=-i
       this%vertex_list(l)%particles(3)=21
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! antiquark-quark to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=9
       this%vertex_list(l)%particles(1)=-i
       this%vertex_list(l)%particles(2)=i
       this%vertex_list(l)%particles(3)=21
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! quark-photon to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=10
       this%vertex_list(l)%particles(1)=i
       this%vertex_list(l)%particles(2)=22
       this%vertex_list(l)%particles(3)=i
       if (mod(i,2).eq.0) then
          this%vertex_list(l)%coupl=[ 2d0/3d0, 2d0/3d0]
       else
          this%vertex_list(l)%coupl=[-1d0/3d0,-1d0/3d0]
       endif
    enddo
    ! antiquark-photon to antiquark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=11
       this%vertex_list(l)%particles(1)=-i
       this%vertex_list(l)%particles(2)=22
       this%vertex_list(l)%particles(3)=-i
       if (mod(i,2).eq.0) then
          this%vertex_list(l)%coupl=[ 2d0/3d0, 2d0/3d0]
       else
          this%vertex_list(l)%coupl=[-1d0/3d0,-1d0/3d0]
       endif
    enddo
    ! quark-Zboson to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=10
       this%vertex_list(l)%particles(1)=i
       this%vertex_list(l)%particles(2)=23
       this%vertex_list(l)%particles(3)=i
       gw=1d0/sw
       fact=1d0/(2d0*sqrt(1d0-sw**2))
       if (mod(i,2).eq.0) then
          Vf=0.5d0-4d0*sw**2/3d0
          Af=0.5d0
       else
          Vf=-0.5d0+2d0*sw**2/3d0
          Af=-0.5d0
       endif
       this%vertex_list(l)%coupl=[Vf+Af,Vf-Af]*gw*fact
    enddo
    ! antiquark-Zboson to antiquark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=11
       this%vertex_list(l)%particles(1)=-i
       this%vertex_list(l)%particles(2)=23
       this%vertex_list(l)%particles(3)=-i
       gw=1d0/sw
       fact=1d0/(2d0*sqrt(1d0-sw**2))
       if (mod(i,2).eq.0) then
          Vf=0.5d0-4d0*sw**2/3d0
          Af=0.5d0
       else
          Vf=-0.5d0+2d0*sw**2/3d0
          Af=-0.5d0
       endif
       this%vertex_list(l)%coupl=[Vf+Af,Vf-Af]*gw*fact
    enddo
    ! quark-Wboson to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=10
       this%vertex_list(l)%particles(1)=i
       if (mod(i,2).eq.0) then
          this%vertex_list(l)%particles(2)=-24
          this%vertex_list(l)%particles(3)=i-1
       else
          this%vertex_list(l)%particles(2)=24
          this%vertex_list(l)%particles(3)=i+1
       endif
       gw=1d0/sw
       fact=1d0/(sqrt(2d0))
       this%vertex_list(l)%coupl=[gw*fact,0d0]
    enddo
    ! antiquark-Wboson to antiquark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=11
       this%vertex_list(l)%particles(1)=-i
       if (mod(i,2).eq.0) then
          this%vertex_list(l)%particles(2)=+24
          this%vertex_list(l)%particles(3)=-i+1
       else
          this%vertex_list(l)%particles(2)=-24
          this%vertex_list(l)%particles(3)=-i-1
       endif
       gw=1d0/sw
       fact=1d0/(sqrt(2d0))
       this%vertex_list(l)%coupl=[gw*fact,0d0]
    enddo
    ! Wboson-Wboson to Z-boson
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=23
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0]
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=23
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0]
    ! Wboson-Wboson to photon
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=22
    this%vertex_list(l)%coupl=[-1d0,0d0]
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=22
    this%vertex_list(l)%coupl=[1d0,0d0]
    ! Wboson-photon to Wboson
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=22
    this%vertex_list(l)%particles(3)=24
    this%vertex_list(l)%coupl=[1d0,0d0]
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=22
    this%vertex_list(l)%particles(3)=-24
    this%vertex_list(l)%coupl=[-1d0,0d0]
    ! photon-Wboson to Wboson
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=22
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=24
    this%vertex_list(l)%coupl=[-1d0,0d0]
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=22
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=-24
    this%vertex_list(l)%coupl=[1d0,0d0]
    ! Wboson-Zboson to Wboson
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=24
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0]
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=-24
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0]
    ! Zboson-Wboson to Wboson
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=24
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0]
    l=l+1
    this%vertex_list(l)%type=12
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=-24
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0]


    ! Wboson-Wboson to Ztensor
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=-23
    gw=1d0/sw
    this%vertex_list(l)%coupl=[gw,0d0] !!
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=-23
    gw=1d0/sw
    this%vertex_list(l)%coupl=[-gw,0d0] !!
    ! Ztensor-Wboson to Wboson
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=-23
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=24
    gw=1d0/sw
    this%vertex_list(l)%coupl=[gw,0d0] !!
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=-23
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=-24
    gw=1d0/sw
    this%vertex_list(l)%coupl=[-gw,0d0] !!
    ! Wboson-Ztensor to Wboson
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=-23
    this%vertex_list(l)%particles(3)=24
    gw=1d0/sw
    this%vertex_list(l)%coupl=[-gw,0d0] !!
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=-23
    this%vertex_list(l)%particles(3)=-24
    gw=1d0/sw
    this%vertex_list(l)%coupl=[gw,0d0] !!


    ! Wboson-photon to Wtensor
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=22
    this%vertex_list(l)%particles(3)=26
    this%vertex_list(l)%coupl=[1d0,0d0]
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=22
    this%vertex_list(l)%particles(3)=-26
    this%vertex_list(l)%coupl=[1d0,0d0]
    ! photon-Wboson to Wtensor
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=22
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=26
    this%vertex_list(l)%coupl=[-1d0,0d0]
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=22
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=-26
    this%vertex_list(l)%coupl=[-1d0,0d0]
    ! Wboson-Zboson to Wtensor
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=26
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0] !!
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=-26
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0] !!
    ! Zboson-Wboson to Wtensor
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=26
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0] !!
    l=l+1
    this%vertex_list(l)%type=13
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=-26
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0] !!


    ! Wtensor-photon to Wboson
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=26
    this%vertex_list(l)%particles(2)=22
    this%vertex_list(l)%particles(3)=24
    this%vertex_list(l)%coupl=[1d0,0d0] 
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=-26
    this%vertex_list(l)%particles(2)=22
    this%vertex_list(l)%particles(3)=-24
    this%vertex_list(l)%coupl=[1d0,0d0]
    ! Wtensor-Wboson to photon
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=26
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=22
    this%vertex_list(l)%coupl=[-1d0,0d0] 
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=-26
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=22
    this%vertex_list(l)%coupl=[-1d0,0d0]
    ! Wtensor-Zboson to Wboson
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=26
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=24
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0] 
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=-26
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=-24
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0]
    ! Wtensor-Wboson to Zboson
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=26
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=23
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0] !!
    l=l+1
    this%vertex_list(l)%type=14
    this%vertex_list(l)%particles(1)=-26
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=23
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0] !!
    ! photon-Wtensor to Wboson
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=22
    this%vertex_list(l)%particles(2)=26
    this%vertex_list(l)%particles(3)=24
    this%vertex_list(l)%coupl=[-1d0,0d0] 
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=22
    this%vertex_list(l)%particles(2)=-26
    this%vertex_list(l)%particles(3)=-24
    this%vertex_list(l)%coupl=[-1d0,0d0]
    ! Wboson-Wtensor to photon
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=-26
    this%vertex_list(l)%particles(3)=22
    this%vertex_list(l)%coupl=[+1d0,0d0] 
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=26
    this%vertex_list(l)%particles(3)=22
    this%vertex_list(l)%coupl=[+1d0,0d0]
    ! Zboson-Wtensor to Wboson
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=26
    this%vertex_list(l)%particles(3)=24
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0] 
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=-26
    this%vertex_list(l)%particles(3)=-24
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[-gw*fact,0d0]
    ! Wboson-Wtensor to Zboson
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=-26
    this%vertex_list(l)%particles(3)=23
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0]  !!
    l=l+1
    this%vertex_list(l)%type=15
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=26
    this%vertex_list(l)%particles(3)=23
    gw=1d0/sw
    fact=sqrt(1d0-sw**2)
    this%vertex_list(l)%coupl=[gw*fact,0d0] !!
    
    write (*,*) l,'interactions loaded'
  end subroutine init_vert
    
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
  integer function get_dim(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    do i=1,this%npart
       if (this%particle_list(i)%type.eq.ipdg .or. this%particle_list(i)%anti_type.eq.ipdg) then
          get_dim=this%particle_list(i)%dim
          return
       endif
    enddo
    write (*,*) 'Particle not in model (dim)',ipdg
    stop 1
  end function get_dim
  integer function get_inter_dim(this,itype)
    implicit none
    class(physics_model) :: this
    integer :: i,itype
    do i=1,this%nint
       if (this%vertex_list(i)%type.eq.itype) then
          get_inter_dim=this%get_dim(this%vertex_list(i)%particles(3))
          return
       endif
    enddo
    write (*,*) 'Interaction not in model (dim)',itype
    stop 1
  end function get_inter_dim
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
  logical function is_tensor_g(iPDG)
    implicit none
    integer :: iPDG
    if (iPDG.eq.-21) then
       is_tensor_g=.true.
    else
       is_tensor_g=.false.
    endif
  end function is_tensor_g
  logical function is_tensor_z(iPDG)
    implicit none
    integer :: iPDG
    if (iPDG.eq.-23) then
       is_tensor_z=.true.
    else
       is_tensor_z=.false.
    endif
  end function is_tensor_z
  logical function is_tensor_w(iPDG)
    implicit none
    integer :: iPDG
    if (abs(iPDG).eq.26) then
       is_tensor_w=.true.
    else
       is_tensor_w=.false.
    endif
  end function is_tensor_w
  logical function is_tensor6(iPDG)
    implicit none
    integer :: iPDG
    is_tensor6=is_tensor_g(iPDG) .or. is_tensor_z(iPDG) .or. is_tensor_w(iPDG)
  end function is_tensor6
  logical function is_tensor(iPDG)
    implicit none
    integer :: iPDG
    is_tensor=is_tensor_g(iPDG) .or. is_tensor_z(iPDG) .or. is_tensor_w(iPDG)
  end function is_tensor
  logical function is_singlet(iPDG)
    implicit none
    integer :: iPDG
    if (abs(iPDG).lt.6 .or. iPDG.eq.21) then
       is_singlet=.false.
    else
       is_singlet=.true.
    endif
  end function is_singlet
  logical function is_photon(iPDG)
    implicit none
    integer :: iPDG
    if (iPDG.eq.22) then
       is_photon=.true.
    else
       is_photon=.false.
    endif
  end function is_photon
  logical function is_massiveboson(iPDG)
    implicit none
    integer :: iPDG
    if (iPDG.eq.23 .or. abs(iPDG).eq.24) then
       is_massiveboson=.true.
    else
       is_massiveboson=.false.
    endif
  end function is_massiveboson
  logical function is_higgs(iPDG)
    implicit none
    integer :: iPDG
    if (iPDG.eq.25) then
       is_higgs=.true.
    else
       is_higgs=.false.
    endif
  end function is_higgs
  logical function is_jet(iPDG)
    implicit none
    integer :: iPDG
    if (is_quark(iPDG).or.is_antiquark(iPDG).or.is_gluon(iPDG)) then
       is_jet=.true.
    else
       is_jet=.false.
    endif
  end function is_jet
end module particles
