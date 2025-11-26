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
     procedure,public :: init_part,get_mass,get_width,get_spin&
          &,get_antipart,init_vert,get_dim,get_inter_dim,is_quark&
          &,is_antiquark,is_lepton,is_antilepton,is_lepton_any,&
          &is_gluon,is_tensor_g,is_tensor_z,is_tensor_w&
          &,is_tensor6,is_tensor,is_singlet,is_photon,is_massiveboson&
          &,is_higgs,is_jet,is_higgsor
  end type physics_model
contains
  subroutine init_part(this,tmass,twidth,zmass,zwidth,wmass,wwidth,hmass,hwidth)
    implicit none
    class(physics_model) :: this
    integer :: i,l
    real(kind=8) :: tmass,twidth
    real(kind=8) :: zmass,zwidth
    real(kind=8) :: wmass,wwidth
    real(kind=8) :: hmass,hwidth
    
    l=0
    this%npart=24! gluon, 6 quarks, tensor, photon, Z-boson and W-boson, H-boson,etc.,6 leptons
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

    ! gluon-U1
    l=l+1
    this%particle_list(l)%type=99
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=-1 ! ill-defined
    this%particle_list(l)%anti_type=99
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

    ! Higgs-boson
    l=l+1
    this%particle_list(l)%type=25
    this%particle_list(l)%mass=hmass
    this%particle_list(l)%width=hwidth
    this%particle_list(l)%spin=1 ! one spin states
    this%particle_list(l)%anti_type=25
    this%particle_list(l)%dim=1

    ! Higgs"or"A (non-propagator scalar auxiliary particle to decompose 4-boson interactions)
    l=l+1
    this%particle_list(l)%type=125
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=-1 ! ill-defined
    this%particle_list(l)%anti_type=125
    this%particle_list(l)%dim=1
    ! Higgs"or"B (non-propagator scalar auxiliary particle to decompose 4-boson interactions)
    l=l+1
    this%particle_list(l)%type=126
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=-1 ! ill-defined
    this%particle_list(l)%anti_type=126
    this%particle_list(l)%dim=1
    ! Higgs"or"C (non-propagator scalar auxiliary particle to decompose 4-boson interactions)
    l=l+1
    this%particle_list(l)%type=127
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=-1 ! ill-defined
    this%particle_list(l)%anti_type=127
    this%particle_list(l)%dim=1

    ! Wtensor (non-propagator auxiliary particle to decompose 4-boson interactions)
    l=l+1
    this%particle_list(l)%type=26
    this%particle_list(l)%mass=0d0
    this%particle_list(l)%width=0d0
    this%particle_list(l)%spin=-1 ! ill-defined
    this%particle_list(l)%anti_type=-26
    this%particle_list(l)%dim=6

    ! charged leptons
    do i=1,3
       l=l+1
       this%particle_list(l)%type=11+(2*i-2)
       this%particle_list(l)%mass=0d0
       this%particle_list(l)%width=0d0
       this%particle_list(l)%spin=2 ! two spin states
       this%particle_list(l)%anti_type=-(11+(2*i-2))
       this%particle_list(l)%dim=4
    enddo

    ! neutral leptons
    do i=1,3
       l=l+1
       this%particle_list(l)%type=12+(2*i-2)
       this%particle_list(l)%mass=0d0
       this%particle_list(l)%width=0d0
       this%particle_list(l)%spin=2 ! two spin states
       this%particle_list(l)%anti_type=-(12+(2*i-2))
       this%particle_list(l)%dim=4
    enddo

    write (99,*) l,'particles loaded'

  end subroutine init_part

  subroutine init_vert(this)
    implicit none
    class(physics_model) :: this
    integer :: i,l
    real(kind=8) :: fact,gw,Vf,Af
    l=0
    this%nint = 222 ! number of vertices
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
    ! antiquark-quark to gluon vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=9
       this%vertex_list(l)%particles(1)=-i
       this%vertex_list(l)%particles(2)=i
       this%vertex_list(l)%particles(3)=21
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! quark-antiquark to gluonU1 vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=8
       this%vertex_list(l)%particles(1)=i
       this%vertex_list(l)%particles(2)=-i
       this%vertex_list(l)%particles(3)=99
       this%vertex_list(l)%coupl=[1d0/3d0,0d0]
    enddo
    ! gluonU1-quark to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=4
       this%vertex_list(l)%particles(1)=99
       this%vertex_list(l)%particles(2)=i
       this%vertex_list(l)%particles(3)=i
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! gluonU1-antiquark to antiquark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=5
       this%vertex_list(l)%particles(1)=99
       this%vertex_list(l)%particles(2)=-i
       this%vertex_list(l)%particles(3)=-i
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! quark-gluonU1 to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=6
       this%vertex_list(l)%particles(1)=i
       this%vertex_list(l)%particles(2)=99
       this%vertex_list(l)%particles(3)=i
       this%vertex_list(l)%coupl=[1d0,0d0]
    enddo
    ! antiquark-gluonU1 to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=7
       this%vertex_list(l)%particles(1)=-i
       this%vertex_list(l)%particles(2)=99
       this%vertex_list(l)%particles(3)=-i
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

    ! quark-Higgs to quark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=16
       this%vertex_list(l)%particles(1)=i
       this%vertex_list(l)%particles(2)=25
       this%vertex_list(l)%particles(3)=i
       this%vertex_list(l)%coupl=[this%get_mass(i)/(this%get_mass(24)*2d0*sw),0d0]
    enddo
    ! antiquark-Higgs to antiquark vertices
    do i=1,6
       l=l+1
       this%vertex_list(l)%type=16
       this%vertex_list(l)%particles(1)=-i
       this%vertex_list(l)%particles(2)=25
       this%vertex_list(l)%particles(3)=-i
       this%vertex_list(l)%coupl=[this%get_mass(i)/(this%get_mass(24)*2d0*sw),0d0]
    enddo
    ! Wboson-Wboson to Higgs
    l=l+1
    this%vertex_list(l)%type=17
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=25
    this%vertex_list(l)%coupl=[this%get_mass(24)/sw,0d0]  !!
    l=l+1
    this%vertex_list(l)%type=17
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=25
    this%vertex_list(l)%coupl=[this%get_mass(24)/sw,0d0]  !!
    ! Higgs-Wboson to Wboson
    l=l+1
    this%vertex_list(l)%type=18
    this%vertex_list(l)%particles(1)=25
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=24
    this%vertex_list(l)%coupl=[this%get_mass(24)/sw,0d0]  !!
    l=l+1
    this%vertex_list(l)%type=18
    this%vertex_list(l)%particles(1)=25
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=-24
    this%vertex_list(l)%coupl=[this%get_mass(24)/sw,0d0]  !!
    ! Wboson-Higgs to Wboson
    l=l+1
    this%vertex_list(l)%type=19
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=25
    this%vertex_list(l)%particles(3)=24
    this%vertex_list(l)%coupl=[this%get_mass(24)/sw,0d0]  !!
    l=l+1
    this%vertex_list(l)%type=19
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=25
    this%vertex_list(l)%particles(3)=-24
    this%vertex_list(l)%coupl=[this%get_mass(24)/sw,0d0]  !!
    ! Zboson-Zboson to Higgs
    l=l+1
    this%vertex_list(l)%type=17
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=25
    this%vertex_list(l)%coupl=[this%get_mass(23)/(sw*dsqrt(1d0-sw**2)),0d0]  !!
    ! Higgs-Zboson to Zboson
    l=l+1
    this%vertex_list(l)%type=18
    this%vertex_list(l)%particles(1)=25
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=23
    this%vertex_list(l)%coupl=[this%get_mass(23)/(sw*dsqrt(1d0-sw**2)),0d0]  !!
    ! Zboson-Higgs to Zboson
    l=l+1
    this%vertex_list(l)%type=19
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=25
    this%vertex_list(l)%particles(3)=23
    this%vertex_list(l)%coupl=[this%get_mass(23)/(sw*dsqrt(1d0-sw**2)),0d0]  !!
    ! Higgs-Higgs to Higgs
    l=l+1
    this%vertex_list(l)%type=20
    this%vertex_list(l)%particles(1)=25
    this%vertex_list(l)%particles(2)=25
    this%vertex_list(l)%particles(3)=25
    this%vertex_list(l)%coupl=[(-3d0/2d0)/sw*(this%get_mass(25)**2/this%get_mass(24)),0d0]  !!

    ! Wboson-Wboson to HiggsorC
    l=l+1
    this%vertex_list(l)%type=17
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=127
    this%vertex_list(l)%coupl=[1d0/2d0/sw**2,0d0]  !!
    l=l+1
    this%vertex_list(l)%type=17
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=127
    this%vertex_list(l)%coupl=[1d0/2d0/sw**2,0d0]  !!
    !! Zboson-Zboson to HiggsorC
    l=l+1
    this%vertex_list(l)%type=17
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=127
    this%vertex_list(l)%coupl=[1d0/2d0/sw**2/(1d0-sw**2),0d0]  !!
    ! HiggsorA-Higgs to Higgs
    l=l+1
    this%vertex_list(l)%type=20
    this%vertex_list(l)%particles(1)=125
    this%vertex_list(l)%particles(2)=25
    this%vertex_list(l)%particles(3)=25
    this%vertex_list(l)%coupl=[1d0,0d0]  !!
    ! Higgs-HiggsorA to Higgs
    l=l+1
    this%vertex_list(l)%type=20
    this%vertex_list(l)%particles(1)=25
    this%vertex_list(l)%particles(2)=125
    this%vertex_list(l)%particles(3)=25
    this%vertex_list(l)%coupl=[1d0,0d0]  !!
    ! Higgs-Higgs to HiggsorA
    l=l+1
    this%vertex_list(l)%type=20
    this%vertex_list(l)%particles(1)=25
    this%vertex_list(l)%particles(2)=25
    this%vertex_list(l)%particles(3)=125
    this%vertex_list(l)%coupl=[(-3d0/4d0)/sw**2*this%get_mass(25)**2/this%get_mass(24)**2,0d0]  !!
    ! HiggsorC-Higgs to Higgs
    l=l+1
    this%vertex_list(l)%type=20
    this%vertex_list(l)%particles(1)=127
    this%vertex_list(l)%particles(2)=25
    this%vertex_list(l)%particles(3)=25
    this%vertex_list(l)%coupl=[1d0,0d0]  !!
    ! Higgs-HiggsorC to Higgs
    l=l+1
    this%vertex_list(l)%type=20
    this%vertex_list(l)%particles(1)=25
    this%vertex_list(l)%particles(2)=127
    this%vertex_list(l)%particles(3)=25
    this%vertex_list(l)%coupl=[1d0,0d0]  !!
    ! Higgs-Higgs to HiggsorB
    l=l+1
    this%vertex_list(l)%type=20
    this%vertex_list(l)%particles(1)=25
    this%vertex_list(l)%particles(2)=25
    this%vertex_list(l)%particles(3)=126
    this%vertex_list(l)%coupl=[1d0,-10d0]  !!

    ! HiggsorB-Zboson to Zboson
    l=l+1
    this%vertex_list(l)%type=18
    this%vertex_list(l)%particles(1)=126
    this%vertex_list(l)%particles(2)=23
    this%vertex_list(l)%particles(3)=23
    this%vertex_list(l)%coupl=[-1d0/2d0/sw**2/(1d0-sw**2),0d0]  !!
    ! Zboson-HiggsorB to Zboson
    l=l+1
    this%vertex_list(l)%type=19
    this%vertex_list(l)%particles(1)=23
    this%vertex_list(l)%particles(2)=126
    this%vertex_list(l)%particles(3)=23
    this%vertex_list(l)%coupl=[-1d0/2d0/sw**2/(1d0-sw**2),0d0]  !!

    ! HiggsorB-Wboson to Wboson
    l=l+1
    this%vertex_list(l)%type=18
    this%vertex_list(l)%particles(1)=126
    this%vertex_list(l)%particles(2)=24
    this%vertex_list(l)%particles(3)=24
    this%vertex_list(l)%coupl=[-1d0/2d0/sw**2,0d0]  !!
    l=l+1
    this%vertex_list(l)%type=18
    this%vertex_list(l)%particles(1)=126
    this%vertex_list(l)%particles(2)=-24
    this%vertex_list(l)%particles(3)=-24
    this%vertex_list(l)%coupl=[-1d0/2d0/sw**2,0d0]  !!
    ! Wboson-HiggsorB to Wboson
    l=l+1
    this%vertex_list(l)%type=19
    this%vertex_list(l)%particles(1)=24
    this%vertex_list(l)%particles(2)=126
    this%vertex_list(l)%particles(3)=24
    this%vertex_list(l)%coupl=[-1d0/2d0/sw**2,0d0]  !!
    l=l+1
    this%vertex_list(l)%type=19
    this%vertex_list(l)%particles(1)=-24
    this%vertex_list(l)%particles(2)=126
    this%vertex_list(l)%particles(3)=-24
    this%vertex_list(l)%coupl=[-1d0/2d0/sw**2,0d0]  !!

    ! lepton-alepton to photon
    do i=1,3
    l=l+1
    this%vertex_list(l)%type=21
    this%vertex_list(l)%particles(1)=11+(2*i-2)
    this%vertex_list(l)%particles(2)=-(11+(2*i-2))
    this%vertex_list(l)%particles(3)=22
    this%vertex_list(l)%coupl=[ -1d0, -1d0]
    !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo
    ! alepton-lepton to photon
    do i=1,3
    l=l+1
    this%vertex_list(l)%type=22
    this%vertex_list(l)%particles(1)=-(11+(2*i-2))
    this%vertex_list(l)%particles(2)=11+(2*i-2)
    this%vertex_list(l)%particles(3)=22
    this%vertex_list(l)%coupl=[ -1d0, -1d0]
    !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo
    ! lepton-alepton to Zboson
    do i=1,3
    l=l+1
    this%vertex_list(l)%type=21
    this%vertex_list(l)%particles(1)=(11+(2*i-2))
    this%vertex_list(l)%particles(2)=-(11+(2*i-2))
    this%vertex_list(l)%particles(3)=23
    gw=1d0/sw
    fact=1d0/(2d0*sqrt(1d0-sw**2))
    Vf=-0.5d0+2d0*sw**2
    Af=-0.5d0
    this%vertex_list(l)%coupl=[Vf+Af,Vf-Af]*gw*fact
    !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo
    ! alepton-lepton to Zboson
    do i=1,3
    l=l+1
    this%vertex_list(l)%type=22
    this%vertex_list(l)%particles(1)=-(11+(2*i-2))
    this%vertex_list(l)%particles(2)=(11+(2*i-2))
    this%vertex_list(l)%particles(3)=23
    gw=1d0/sw
    fact=1d0/(2d0*sqrt(1d0-sw**2))
    Vf=-0.5d0+2d0*sw**2
    Af=-0.5d0
    this%vertex_list(l)%coupl=[Vf+Af,Vf-Af]*gw*fact
    !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo

    ! charged lepton-lepton to Wboson
    do i=1,3
    l=l+1
    this%vertex_list(l)%type=21
    this%vertex_list(l)%particles(1)=11+(2*i-2)
    this%vertex_list(l)%particles(2)=-(12+(2*i-2))
    this%vertex_list(l)%particles(3)=-24
    gw=1d0/sw
    fact=1d0/(sqrt(2d0))
    this%vertex_list(l)%coupl=[gw*fact,0d0]
    enddo
    ! charged lepton-lepton to Wboson
    do i=1,3
    l=l+1
    this%vertex_list(l)%type=22
    this%vertex_list(l)%particles(1)=-(11+(2*i-2))
    this%vertex_list(l)%particles(2)=(12+(2*i-2))
    this%vertex_list(l)%particles(3)=24
    gw=1d0/sw
    fact=1d0/(sqrt(2d0))
    this%vertex_list(l)%coupl=[gw*fact,0d0]
    enddo

    ! lepton-photon to lepton vertices
    do i=1,3
       l=l+1
       this%vertex_list(l)%type=10
       this%vertex_list(l)%particles(1)=(11+(2*i-2))
       this%vertex_list(l)%particles(2)=22
       this%vertex_list(l)%particles(3)=(11+(2*i-2))
       this%vertex_list(l)%coupl=[ -1d0, -1d0]
       !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo
    ! antilepton-photon to antilepton vertices
    do i=1,3
       l=l+1
       this%vertex_list(l)%type=11
       this%vertex_list(l)%particles(1)=-(11+(2*i-2))
       this%vertex_list(l)%particles(2)=22
       this%vertex_list(l)%particles(3)=-(11+(2*i-2))
       this%vertex_list(l)%coupl=[ -1d0, -1d0]
       !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo

    ! photon-lepton to lepton vertices
    do i=1,3
       l=l+1
       this%vertex_list(l)%type=23
       this%vertex_list(l)%particles(1)=22
       this%vertex_list(l)%particles(2)=(11+(2*i-2))
       this%vertex_list(l)%particles(3)=(11+(2*i-2))
       this%vertex_list(l)%coupl=[ -1d0, -1d0]
       !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo
    ! photon-antilepton to antilepton vertices
    do i=1,3
       l=l+1
       this%vertex_list(l)%type=24
       this%vertex_list(l)%particles(1)=22
       this%vertex_list(l)%particles(2)=-(11+(2*i-2))
       this%vertex_list(l)%particles(3)=-(11+(2*i-2))
       this%vertex_list(l)%coupl=[ -1d0, -1d0]
       !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo

    ! lepton-Zboson to lepton vertices
    do i=1,3
       l=l+1
       this%vertex_list(l)%type=10
       this%vertex_list(l)%particles(1)=(11+(2*i-2))
       this%vertex_list(l)%particles(2)=23
       this%vertex_list(l)%particles(3)=(11+(2*i-2))
       gw=1d0/sw
       fact=1d0/(2d0*sqrt(1d0-sw**2))
       Vf=-0.5d0+2d0*sw**2
       Af=-0.5d0
       this%vertex_list(l)%coupl=[Vf+Af,Vf-Af]*gw*fact
       !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo
    ! antilepton-Zboson to antilepton vertices
    do i=1,3
       l=l+1
       this%vertex_list(l)%type=11
       this%vertex_list(l)%particles(1)=-(11+(2*i-2))
       this%vertex_list(l)%particles(2)=23
       this%vertex_list(l)%particles(3)=-(11+(2*i-2))
       gw=1d0/sw
       fact=1d0/(2d0*sqrt(1d0-sw**2))
       Vf=-0.5d0+2d0*sw**2
       Af=-0.5d0
       this%vertex_list(l)%coupl=[Vf+Af,Vf-Af]*gw*fact
       !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo
   ! Zboson-lepton to lepton vertices
    do i=1,3
       l=l+1
       this%vertex_list(l)%type=23
       this%vertex_list(l)%particles(1)=23
       this%vertex_list(l)%particles(2)=(11+(2*i-2))
       this%vertex_list(l)%particles(3)=(11+(2*i-2))
       gw=1d0/sw
       fact=1d0/(2d0*sqrt(1d0-sw**2))
       Vf=-0.5d0+2d0*sw**2
       Af=-0.5d0
       this%vertex_list(l)%coupl=[Vf+Af,Vf-Af]*gw*fact
       !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo
    ! Zboson-antilepton to antilepton vertices
    do i=1,3
       l=l+1
       this%vertex_list(l)%type=24
       this%vertex_list(l)%particles(1)=23
       this%vertex_list(l)%particles(2)=-(11+(2*i-2))
       this%vertex_list(l)%particles(3)=-(11+(2*i-2))
       gw=1d0/sw
       fact=1d0/(2d0*sqrt(1d0-sw**2))
       Vf=-0.5d0+2d0*sw**2
       Af=-0.5d0
       this%vertex_list(l)%coupl=[Vf+Af,Vf-Af]*gw*fact
       !this%vertex_list(l)%coupl=[0d0,0d0]
    enddo



    write (99,*) l,'interactions loaded'
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
  logical function is_quark(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.ge.1 .and. iPDG.le.6) then
       is_quark=.true.
    else
       is_quark=.false.
    endif
  end function is_quark
  logical function is_antiquark(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.le.-1 .and. iPDG.ge.-6) then
       is_antiquark=.true.
    else
       is_antiquark=.false.
    endif
  end function is_antiquark
  logical function is_lepton(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.ge.11 .and. iPDG.le.16) then
       is_lepton=.true.
    else
       is_lepton=.false.
    endif
  end function is_lepton
  logical function is_antilepton(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.le.-11 .and. iPDG.ge.-16) then
       is_antilepton=.true.
    else
       is_antilepton=.false.
    endif
  end function is_antilepton
  logical function is_lepton_any(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (abs(iPDG).ge.11 .and. abs(iPDG).le.16) then
       is_lepton_any=.true.
    else
       is_lepton_any=.false.
    endif
  end function is_lepton_any
  logical function is_gluon(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.eq.21 .or. iPDG.eq.99) then
       is_gluon=.true.
    else
       is_gluon=.false.
    endif
  end function is_gluon

  logical function is_scalar(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (abs(iPDG).eq.25.or.abs(iPDG).eq.125.or.abs(iPDG).eq.126.or.abs(iPDG).eq.127) then
       is_scalar=.true.
    else
       is_scalar=.false.
    endif
  end function is_scalar
  logical function is_tensor_g(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.eq.-21) then
       is_tensor_g=.true.
    else
       is_tensor_g=.false.
    endif
  end function is_tensor_g
  logical function is_tensor_z(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.eq.-23) then
       is_tensor_z=.true.
    else
       is_tensor_z=.false.
    endif
  end function is_tensor_z
  logical function is_tensor_w(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (abs(iPDG).eq.26) then
       is_tensor_w=.true.
    else
       is_tensor_w=.false.
    endif
  end function is_tensor_w
  logical function is_tensor6(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_tensor6=this%is_tensor_g(iPDG) .or. this%is_tensor_z(iPDG) .or. this%is_tensor_w(iPDG)
  end function is_tensor6
  logical function is_tensor(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_tensor=this%is_tensor_g(iPDG) .or. this%is_tensor_z(iPDG) .or. this%is_tensor_w(iPDG)
  end function is_tensor
  logical function is_singlet(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (abs(iPDG).le.6 .or. iPDG.eq.21) then
       is_singlet=.false.
    else
       is_singlet=.true.
    endif
  end function is_singlet
  logical function is_photon(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.eq.22) then
       is_photon=.true.
    else
       is_photon=.false.
    endif
  end function is_photon
  logical function is_massiveboson(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.eq.23 .or. abs(iPDG).eq.24) then
       is_massiveboson=.true.
    else
       is_massiveboson=.false.
    endif
  end function is_massiveboson
  logical function is_higgs(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.eq.25) then
       is_higgs=.true.
    else
       is_higgs=.false.
    endif
  end function is_higgs
  logical function is_higgsor(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if (iPDG.eq.125.or.iPDG.eq.126.or.iPDG.eq.127) then
       is_higgsor=.true.
    else
       is_higgsor=.false.
    endif
  end function is_higgsor
  logical function is_jet(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    if ((this%is_quark(iPDG).or.this%is_antiquark(iPDG).or.this%is_gluon(iPDG)) .and. &
         this%get_mass(iPDG).eq.0d0) then
       is_jet=.true.
    else
       is_jet=.false.
    endif
  end function is_jet
end module particles
