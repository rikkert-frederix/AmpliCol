module particles
  implicit none
  real(kind=8),parameter :: sw = 0.47143025548407230d0 
  integer,parameter :: model_particle_capacity = 24
  integer,parameter :: model_vertex_capacity = 222
  private :: append_particle, append_vertex, find_particle_index, particle_property_sign&
       &, weak_coupling, weak_cosine, neutral_gauge_coupling&
       &, charged_current_coupling, particle_species_index&
       &, photon_fermion_coupling, z_fermion_coupling&
       &, weak_coupling_squared, weak_cosine_squared&
       &, weak_coupling_over_cosine, higgs_self_coupling&
       &, model_particle_capacity, model_vertex_capacity
  type particle
     ! weak_isospin and weak_hypercharge use index 1=left, 2=right with Q = T3 + Y/2.
     integer :: type,anti_type,spin,dim,color_rep
     real(kind=8) :: mass,width,charge
     real(kind=8),dimension(2) :: weak_isospin,weak_hypercharge
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
          &,get_charge,get_isospin_l,get_isospin_r,get_hypercharge_l&
          &,get_hypercharge_r,get_color_rep,get_color_dim&
          &,is_antiquark,is_lepton,is_antilepton,is_lepton_any,&
          &is_gluon,is_tensor_g,is_tensor_z,is_tensor_w&
          &,is_tensor6,is_tensor,is_singlet,is_photon,is_massiveboson&
          &,is_higgs,is_jet,is_higgsor
     procedure,public :: get_colour_rep => get_color_rep
     procedure,public :: get_colour_dim => get_color_dim
  end type physics_model
contains
  subroutine append_particle(this,l,ipdg,mass,width,spin,anti_type,dim,charge,&
       &weak_isospin,weak_hypercharge,color_rep)
    implicit none
    class(physics_model),intent(inout) :: this
    integer,intent(inout) :: l
    integer,intent(in) :: ipdg,spin,anti_type,dim,color_rep
    real(kind=8),intent(in) :: mass,width,charge
    real(kind=8),dimension(2),intent(in) :: weak_isospin,weak_hypercharge
    l=l+1
    if (l.gt.this%npart) then
       write (*,*) 'ERROR: more particles than allocated',l,this%npart
       stop 1
    endif
    this%particle_list(l)%type=ipdg
    this%particle_list(l)%mass=mass
    this%particle_list(l)%width=width
    this%particle_list(l)%spin=spin
    this%particle_list(l)%anti_type=anti_type
    this%particle_list(l)%dim=dim
    this%particle_list(l)%charge=charge
    this%particle_list(l)%weak_isospin=weak_isospin
    this%particle_list(l)%weak_hypercharge=weak_hypercharge
    this%particle_list(l)%color_rep=color_rep
  end subroutine append_particle

  subroutine append_vertex(this,l,itype,ipdgs,coupl)
    implicit none
    class(physics_model),intent(inout) :: this
    integer,intent(inout) :: l
    integer,intent(in) :: itype
    integer,dimension(3),intent(in) :: ipdgs
    real(kind=8),dimension(2),intent(in) :: coupl
    l=l+1
    if (l.gt.this%nint) then
       write (*,*) 'ERROR: more vertices than allocated',l,this%nint
       stop 1
    endif
    this%vertex_list(l)%type=itype
    this%vertex_list(l)%particles=ipdgs
    this%vertex_list(l)%coupl=coupl
  end subroutine append_vertex

  integer function find_particle_index(this,ipdg)
    implicit none
    class(physics_model),intent(in) :: this
    integer,intent(in) :: ipdg
    integer :: i
    do i=1,this%npart
       if (this%particle_list(i)%type.eq.ipdg .or. this%particle_list(i)%anti_type.eq.ipdg) then
          find_particle_index=i
          return
       endif
    enddo
    find_particle_index=0
  end function find_particle_index

  integer function particle_property_sign(this,ipdg)
    implicit none
    class(physics_model),intent(in) :: this
    integer,intent(in) :: ipdg
    integer :: i
    do i=1,this%npart
       if (this%particle_list(i)%type.eq.ipdg) then
          particle_property_sign=1
          return
       elseif (this%particle_list(i)%anti_type.eq.ipdg) then
          particle_property_sign=-1
          return
       endif
    enddo
    particle_property_sign=0
  end function particle_property_sign

  real(kind=8) function weak_coupling()
    implicit none
    weak_coupling=1d0/sw
  end function weak_coupling

  real(kind=8) function weak_cosine()
    implicit none
    weak_cosine=sqrt(1d0-sw**2)
  end function weak_cosine

  real(kind=8) function neutral_gauge_coupling()
    implicit none
    neutral_gauge_coupling=weak_coupling()*weak_cosine()
  end function neutral_gauge_coupling

  real(kind=8) function weak_coupling_squared()
    implicit none
    weak_coupling_squared=weak_coupling()**2
  end function weak_coupling_squared

  real(kind=8) function weak_cosine_squared()
    implicit none
    weak_cosine_squared=weak_cosine()**2
  end function weak_cosine_squared

  real(kind=8) function weak_coupling_over_cosine()
    implicit none
    weak_coupling_over_cosine=weak_coupling()/weak_cosine()
  end function weak_coupling_over_cosine

  real(kind=8) function charged_current_coupling()
    implicit none
    charged_current_coupling=weak_coupling()*(1d0/(sqrt(2d0)))
  end function charged_current_coupling

  function higgs_self_coupling(this) result(coupl)
    implicit none
    class(physics_model),intent(in) :: this
    real(kind=8),dimension(2) :: coupl
    coupl=[(-3d0/2d0)*weak_coupling()*(this%get_mass(25)**2/this%get_mass(24)),0d0]
  end function higgs_self_coupling

  integer function particle_species_index(this,ipdg,property)
    implicit none
    class(physics_model),intent(in) :: this
    integer,intent(in) :: ipdg
    character(len=*),intent(in) :: property
    particle_species_index=find_particle_index(this,ipdg)
    if (particle_species_index.gt.0) return
    write (*,*) 'Particle not in model ('//trim(property)//')',ipdg
    stop 1
  end function particle_species_index

  function photon_fermion_coupling(this,ipdg) result(coupl)
    implicit none
    class(physics_model),intent(in) :: this
    integer,intent(in) :: ipdg
    real(kind=8),dimension(2) :: coupl
    real(kind=8) :: charge
    integer :: i
    i=particle_species_index(this,ipdg,'photon coupling')
    charge=this%particle_list(i)%charge
    coupl=[charge,charge]
  end function photon_fermion_coupling

  function z_fermion_coupling(this,ipdg) result(coupl)
    implicit none
    class(physics_model),intent(in) :: this
    integer,intent(in) :: ipdg
    real(kind=8),dimension(2) :: coupl
    integer :: i
    real(kind=8) :: charge,left_isospin,right_isospin
    i=particle_species_index(this,ipdg,'Z coupling')
    charge=this%particle_list(i)%charge
    left_isospin=this%particle_list(i)%weak_isospin(1)
    right_isospin=this%particle_list(i)%weak_isospin(2)
    coupl=weak_coupling_over_cosine()*[left_isospin-charge*sw**2,&
         &right_isospin-charge*sw**2]
  end function z_fermion_coupling

  subroutine init_part(this,tmass,twidth,zmass,zwidth,wmass,wwidth,hmass,hwidth)
    implicit none
    class(physics_model) :: this
    integer :: i,l
    real(kind=8) :: tmass,twidth
    real(kind=8) :: zmass,zwidth
    real(kind=8) :: wmass,wwidth
    real(kind=8) :: hmass,hwidth
    l=0
    this%npart=model_particle_capacity ! gluon, quarks, tensors, bosons, Higgs and leptons
    allocate(this%particle_list(this%npart))

    ! 5 massless quarks
    do i=1,5
       if (mod(i,2).eq.0) then
          call append_particle(this,l,i,0d0,0d0,2,-i,4,2d0/3d0,[0.5d0,0d0],[1d0/3d0,4d0/3d0],3)
       else
          call append_particle(this,l,i,0d0,0d0,2,-i,4,-1d0/3d0,[-0.5d0,0d0],[1d0/3d0,-2d0/3d0],3)
       endif
    enddo

    ! top quark
    call append_particle(this,l,6,tmass,twidth,2,-6,4,2d0/3d0,[0.5d0,0d0],[1d0/3d0,4d0/3d0],3)

    ! gluon
    call append_particle(this,l,21,0d0,0d0,2,21,4,0d0,[0d0,0d0],[0d0,0d0],8)

    ! gluon-U1
    call append_particle(this,l,99,0d0,0d0,-1,99,4,0d0,[0d0,0d0],[0d0,0d0],1)

    ! tensor (non-propagator auxiliary particle to decompose 4-gluon interaction)
    call append_particle(this,l,-21,0d0,0d0,-1,-21,6,0d0,[0d0,0d0],[0d0,0d0],8)

    ! photon
    call append_particle(this,l,22,0d0,0d0,2,22,4,0d0,[0d0,0d0],[0d0,0d0],1)

    ! Z-boson
    call append_particle(this,l,23,zmass,zwidth,3,23,4,0d0,[0d0,0d0],[0d0,0d0],1)

    ! Ztensor (non-propagator auxiliary particle to decompose 4-Wboson interaction)
    call append_particle(this,l,-23,0d0,0d0,-1,-23,6,0d0,[0d0,0d0],[0d0,0d0],1)

    ! W-boson
    call append_particle(this,l,24,wmass,wwidth,3,-24,4,1d0,[1d0,1d0],[0d0,0d0],1)

    ! Higgs-boson
    call append_particle(this,l,25,hmass,hwidth,1,25,1,0d0,[-0.5d0,-0.5d0],[1d0,1d0],1)

    ! Higgs"or"A (non-propagator scalar auxiliary particle to decompose 4-boson interactions)
    call append_particle(this,l,125,0d0,0d0,-1,125,1,0d0,[0d0,0d0],[0d0,0d0],1)
    ! Higgs"or"B (non-propagator scalar auxiliary particle to decompose 4-boson interactions)
    call append_particle(this,l,126,0d0,0d0,-1,126,1,0d0,[0d0,0d0],[0d0,0d0],1)
    ! Higgs"or"C (non-propagator scalar auxiliary particle to decompose 4-boson interactions)
    call append_particle(this,l,127,0d0,0d0,-1,127,1,0d0,[0d0,0d0],[0d0,0d0],1)

    ! Wtensor (non-propagator auxiliary particle to decompose 4-boson interactions)
    call append_particle(this,l,26,0d0,0d0,-1,-26,6,1d0,[1d0,1d0],[0d0,0d0],1)

    ! charged leptons
    do i=1,3
       call append_particle(this,l,11+(2*i-2),0d0,0d0,2,-(11+(2*i-2)),4,-1d0,[-0.5d0,0d0],[-1d0,-2d0],1)
    enddo

    ! neutrinos
    do i=1,3
       call append_particle(this,l,12+(2*i-2),0d0,0d0,2,-(12+(2*i-2)),4,0d0,[0.5d0,0d0],[-1d0,0d0],1)
    enddo

    this%npart=l
    write (99,*) l,'particles loaded'
  end subroutine init_part

  subroutine init_vert(this)
    implicit none
    class(physics_model) :: this
    integer :: i,l
    l=0
    this%nint = model_vertex_capacity ! number of vertices
    allocate(this%vertex_list(this%nint))
    ! gluon-gluon to gluon vertex
    call append_vertex(this,l,0,[21,21,21],[1d0,0d0])
    ! gluon-gluon to tensor vertex
    call append_vertex(this,l,1,[21,21,-21],[1d0,0d0])
    ! tensor-gluon to gluon vertex
    call append_vertex(this,l,2,[-21,21,21],[1d0,0d0])
    ! gluon-tensor to gluon vertex
    call append_vertex(this,l,3,[21,-21,21],[1d0,0d0])
    ! gluon-quark to quark vertices
    do i=1,6
       call append_vertex(this,l,4,[21,i,i],[1d0,0d0])
    enddo
    ! gluon-antiquark to antiquark vertices
    do i=1,6
       call append_vertex(this,l,5,[21,-i,-i],[1d0,0d0])
    enddo
    ! quark-gluon to quark vertices
    do i=1,6
       call append_vertex(this,l,6,[i,21,i],[1d0,0d0])
    enddo
    ! antiquark-gluon to quark vertices
    do i=1,6
       call append_vertex(this,l,7,[-i,21,-i],[1d0,0d0])
    enddo
    ! antiquark-quark to gluon vertices
    do i=1,6
       call append_vertex(this,l,9,[-i,i,21],[1d0,0d0])
    enddo
    ! quark-antiquark to gluonU1 vertices
    do i=1,6
       call append_vertex(this,l,8,[i,-i,99],[1d0/3d0,0d0])  ! 1/N_C coupling
    enddo
    ! gluonU1-quark to quark vertices
    do i=1,6
       call append_vertex(this,l,4,[99,i,i],[1d0,0d0])
    enddo
    ! gluonU1-antiquark to antiquark vertices
    do i=1,6
       call append_vertex(this,l,5,[99,-i,-i],[1d0,0d0])
    enddo
    ! quark-gluonU1 to quark vertices
    do i=1,6
       call append_vertex(this,l,6,[i,99,i],[1d0,0d0])
    enddo
    ! antiquark-gluonU1 to quark vertices
    do i=1,6
       call append_vertex(this,l,7,[-i,99,-i],[1d0,0d0])
    enddo
    ! quark-photon to quark vertices
    do i=1,6
       call append_vertex(this,l,10,[i,22,i],photon_fermion_coupling(this,i))
    enddo
    ! antiquark-photon to antiquark vertices
    do i=1,6
       call append_vertex(this,l,11,[-i,22,-i],photon_fermion_coupling(this,-i))
    enddo
    ! quark-Zboson to quark vertices
    do i=1,6
       call append_vertex(this,l,10,[i,23,i],z_fermion_coupling(this,i))
    enddo
    ! antiquark-Zboson to antiquark vertices
    do i=1,6
       call append_vertex(this,l,11,[-i,23,-i],z_fermion_coupling(this,-i))
    enddo
    ! quark-Wboson to quark vertices
    do i=1,6
       if (mod(i,2).eq.0) then
          call append_vertex(this,l,10,[i,-24,i-1],[charged_current_coupling(),0d0])
       else
          call append_vertex(this,l,10,[i,24,i+1],[charged_current_coupling(),0d0])
       endif
    enddo
    ! antiquark-Wboson to antiquark vertices
    do i=1,6
       if (mod(i,2).eq.0) then
          call append_vertex(this,l,11,[-i,+24,-i+1],[charged_current_coupling(),0d0])
       else
          call append_vertex(this,l,11,[-i,-24,-i-1],[charged_current_coupling(),0d0])
       endif
    enddo
    ! Wboson-Wboson to Z-boson
    call append_vertex(this,l,12,[24,-24,23],[-neutral_gauge_coupling(),0d0])
    call append_vertex(this,l,12,[-24,24,23],[neutral_gauge_coupling(),0d0])
    ! Wboson-Wboson to photon
    call append_vertex(this,l,12,[24,-24,22],[-this%get_charge(24),0d0])
    call append_vertex(this,l,12,[-24,24,22],[this%get_charge(24),0d0])
    ! Wboson-photon to Wboson
    call append_vertex(this,l,12,[24,22,24],[this%get_charge(24),0d0])
    call append_vertex(this,l,12,[-24,22,-24],[this%get_charge(-24),0d0])
    ! photon-Wboson to Wboson
    call append_vertex(this,l,12,[22,24,24],[-this%get_charge(24),0d0])
    call append_vertex(this,l,12,[22,-24,-24],[-this%get_charge(-24),0d0])
    ! Wboson-Zboson to Wboson
    call append_vertex(this,l,12,[24,23,24],[neutral_gauge_coupling(),0d0])
    call append_vertex(this,l,12,[-24,23,-24],[-neutral_gauge_coupling(),0d0])
    ! Zboson-Wboson to Wboson
    call append_vertex(this,l,12,[23,24,24],[-neutral_gauge_coupling(),0d0])
    call append_vertex(this,l,12,[23,-24,-24],[neutral_gauge_coupling(),0d0])


    ! Wboson-Wboson to Ztensor
    call append_vertex(this,l,13,[24,-24,-23],[weak_coupling(),0d0]) !!
    call append_vertex(this,l,13,[-24,24,-23],[-weak_coupling(),0d0]) !!
    ! Ztensor-Wboson to Wboson
    call append_vertex(this,l,14,[-23,24,24],[weak_coupling(),0d0]) !!
    call append_vertex(this,l,14,[-23,-24,-24],[-weak_coupling(),0d0]) !!
    ! Wboson-Ztensor to Wboson
    call append_vertex(this,l,15,[24,-23,24],[-weak_coupling(),0d0]) !!
    call append_vertex(this,l,15,[-24,-23,-24],[weak_coupling(),0d0]) !!


    ! Wboson-photon to Wtensor
    call append_vertex(this,l,13,[24,22,26],[this%get_charge(24),0d0])
    call append_vertex(this,l,13,[-24,22,-26],[-this%get_charge(-24),0d0])
    ! photon-Wboson to Wtensor
    call append_vertex(this,l,13,[22,24,26],[-this%get_charge(24),0d0])
    call append_vertex(this,l,13,[22,-24,-26],[this%get_charge(-24),0d0])
    ! Wboson-Zboson to Wtensor
    call append_vertex(this,l,13,[24,23,26],[neutral_gauge_coupling(),0d0]) !!
    call append_vertex(this,l,13,[-24,23,-26],[neutral_gauge_coupling(),0d0]) !!
    ! Zboson-Wboson to Wtensor
    call append_vertex(this,l,13,[23,24,26],[-neutral_gauge_coupling(),0d0]) !!
    call append_vertex(this,l,13,[23,-24,-26],[-neutral_gauge_coupling(),0d0]) !!


    ! Wtensor-photon to Wboson
    call append_vertex(this,l,14,[26,22,24],[this%get_charge(26),0d0])
    call append_vertex(this,l,14,[-26,22,-24],[-this%get_charge(-26),0d0])
    ! Wtensor-Wboson to photon
    call append_vertex(this,l,14,[26,-24,22],[-this%get_charge(26),0d0])
    call append_vertex(this,l,14,[-26,24,22],[this%get_charge(-26),0d0])
    ! Wtensor-Zboson to Wboson
    call append_vertex(this,l,14,[26,23,24],[neutral_gauge_coupling(),0d0])
    call append_vertex(this,l,14,[-26,23,-24],[neutral_gauge_coupling(),0d0])
    ! Wtensor-Wboson to Zboson
    call append_vertex(this,l,14,[26,-24,23],[-neutral_gauge_coupling(),0d0]) !!
    call append_vertex(this,l,14,[-26,24,23],[-neutral_gauge_coupling(),0d0]) !!
    ! photon-Wtensor to Wboson
    call append_vertex(this,l,15,[22,26,24],[-this%get_charge(26),0d0])
    call append_vertex(this,l,15,[22,-26,-24],[this%get_charge(-26),0d0])
    ! Wboson-Wtensor to photon
    call append_vertex(this,l,15,[24,-26,22],[this%get_charge(24),0d0])
    call append_vertex(this,l,15,[-24,26,22],[-this%get_charge(-24),0d0])
    ! Zboson-Wtensor to Wboson
    call append_vertex(this,l,15,[23,26,24],[-neutral_gauge_coupling(),0d0])
    call append_vertex(this,l,15,[23,-26,-24],[-neutral_gauge_coupling(),0d0])
    ! Wboson-Wtensor to Zboson
    call append_vertex(this,l,15,[24,-26,23],[neutral_gauge_coupling(),0d0]) !!
    call append_vertex(this,l,15,[-24,26,23],[neutral_gauge_coupling(),0d0]) !!

    ! quark-Higgs to quark vertices
    do i=1,6
       if (this%get_mass(i).eq.0d0) cycle
       call append_vertex(this,l,16,[i,25,i],[this%get_mass(i)*weak_coupling()/(this%get_mass(24)*2d0),0d0])
    enddo
    ! antiquark-Higgs to antiquark vertices
    do i=1,6
       if (this%get_mass(i).eq.0d0) cycle
       call append_vertex(this,l,16,[-i,25,-i],[this%get_mass(i)*weak_coupling()/(this%get_mass(24)*2d0),0d0])
    enddo
    ! Wboson-Wboson to Higgs
    call append_vertex(this,l,17,[24,-24,25],[this%get_mass(24)*weak_coupling(),0d0]) !!
    call append_vertex(this,l,17,[-24,24,25],[this%get_mass(24)*weak_coupling(),0d0]) !!
    ! Higgs-Wboson to Wboson
    call append_vertex(this,l,18,[25,24,24],[this%get_mass(24)*weak_coupling(),0d0]) !!
    call append_vertex(this,l,18,[25,-24,-24],[this%get_mass(24)*weak_coupling(),0d0]) !!
    ! Wboson-Higgs to Wboson
    call append_vertex(this,l,19,[24,25,24],[this%get_mass(24)*weak_coupling(),0d0]) !!
    call append_vertex(this,l,19,[-24,25,-24],[this%get_mass(24)*weak_coupling(),0d0]) !!
    ! Zboson-Zboson to Higgs
    call append_vertex(this,l,17,[23,23,25],[this%get_mass(23)*weak_coupling_over_cosine(),0d0]) !!
    ! Higgs-Zboson to Zboson
    call append_vertex(this,l,18,[25,23,23],[this%get_mass(23)*weak_coupling_over_cosine(),0d0]) !!
    ! Zboson-Higgs to Zboson
    call append_vertex(this,l,19,[23,25,23],[this%get_mass(23)*weak_coupling_over_cosine(),0d0]) !!
    ! Higgs-Higgs to Higgs
    call append_vertex(this,l,20,[25,25,25],higgs_self_coupling(this)) !!

    ! Wboson-Wboson to HiggsorC
    call append_vertex(this,l,17,[24,-24,127],[1d0/2d0*weak_coupling_squared(),0d0]) !!
    call append_vertex(this,l,17,[-24,24,127],[1d0/2d0*weak_coupling_squared(),0d0]) !!
    !! Zboson-Zboson to HiggsorC
    call append_vertex(this,l,17,[23,23,127],[1d0/2d0*weak_coupling_squared()/weak_cosine_squared(),0d0]) !!
    ! HiggsorA-Higgs to Higgs
    call append_vertex(this,l,20,[125,25,25],[1d0,0d0]) !!
    ! Higgs-HiggsorA to Higgs
    call append_vertex(this,l,20,[25,125,25],[1d0,0d0]) !!
    ! Higgs-Higgs to HiggsorA
    call append_vertex(this,l,20,[25,25,125],[(-3d0/4d0)*weak_coupling_squared()*this%get_mass(25)**2/this%get_mass(24)**2,0d0]) !!
    ! HiggsorC-Higgs to Higgs
    call append_vertex(this,l,20,[127,25,25],[1d0,0d0]) !!
    ! Higgs-HiggsorC to Higgs
    call append_vertex(this,l,20,[25,127,25],[1d0,0d0]) !!
    ! Higgs-Higgs to HiggsorB
    call append_vertex(this,l,20,[25,25,126],[1d0,-10d0]) !!

    ! HiggsorB-Zboson to Zboson
    call append_vertex(this,l,18,[126,23,23],[-1d0/2d0*weak_coupling_squared()/weak_cosine_squared(),0d0]) !!
    ! Zboson-HiggsorB to Zboson
    call append_vertex(this,l,19,[23,126,23],[-1d0/2d0*weak_coupling_squared()/weak_cosine_squared(),0d0]) !!

    ! HiggsorB-Wboson to Wboson
    call append_vertex(this,l,18,[126,24,24],[-1d0/2d0*weak_coupling_squared(),0d0]) !!
    call append_vertex(this,l,18,[126,-24,-24],[-1d0/2d0*weak_coupling_squared(),0d0]) !!
    ! Wboson-HiggsorB to Wboson
    call append_vertex(this,l,19,[24,126,24],[-1d0/2d0*weak_coupling_squared(),0d0]) !!
    call append_vertex(this,l,19,[-24,126,-24],[-1d0/2d0*weak_coupling_squared(),0d0]) !!

    ! lepton-alepton to photon
    do i=1,3
    call append_vertex(this,l,21,[11+(2*i-2),-(11+(2*i-2)),22],photon_fermion_coupling(this,11+(2*i-2)))
    enddo
    ! alepton-lepton to photon
    do i=1,3
    call append_vertex(this,l,22,[-(11+(2*i-2)),11+(2*i-2),22],photon_fermion_coupling(this,-(11+(2*i-2))))
    enddo
    ! lepton-alepton to Zboson
    do i=1,3
    call append_vertex(this,l,21,[(11+(2*i-2)),-(11+(2*i-2)),23],z_fermion_coupling(this,11+(2*i-2)))
    enddo
    ! alepton-lepton to Zboson
    do i=1,3
    call append_vertex(this,l,22,[-(11+(2*i-2)),(11+(2*i-2)),23],z_fermion_coupling(this,-(11+(2*i-2))))
    enddo
    ! charged lepton-lepton to Wboson
    do i=1,3
    call append_vertex(this,l,21,[11+(2*i-2),-(12+(2*i-2)),-24],[charged_current_coupling(),0d0])
    enddo
    ! charged lepton-lepton to Wboson
    do i=1,3
    call append_vertex(this,l,22,[-(11+(2*i-2)),(12+(2*i-2)),24],[charged_current_coupling(),0d0])
    enddo
    ! lepton-photon to lepton vertices
    do i=1,3
       call append_vertex(this,l,10,[(11+(2*i-2)),22,(11+(2*i-2))],photon_fermion_coupling(this,11+(2*i-2)))
    enddo
    ! antilepton-photon to antilepton vertices
    do i=1,3
       call append_vertex(this,l,11,[-(11+(2*i-2)),22,-(11+(2*i-2))],photon_fermion_coupling(this,-(11+(2*i-2))))
    enddo
    ! photon-lepton to lepton vertices
    do i=1,3
       call append_vertex(this,l,23,[22,(11+(2*i-2)),(11+(2*i-2))],photon_fermion_coupling(this,11+(2*i-2)))
    enddo
    ! photon-antilepton to antilepton vertices
    do i=1,3
       call append_vertex(this,l,24,[22,-(11+(2*i-2)),-(11+(2*i-2))],photon_fermion_coupling(this,-(11+(2*i-2))))
    enddo
    ! lepton-Zboson to lepton vertices
    do i=1,3
       call append_vertex(this,l,10,[(11+(2*i-2)),23,(11+(2*i-2))],z_fermion_coupling(this,11+(2*i-2)))
    enddo
    ! antilepton-Zboson to antilepton vertices
    do i=1,3
       call append_vertex(this,l,11,[-(11+(2*i-2)),23,-(11+(2*i-2))],z_fermion_coupling(this,-(11+(2*i-2))))
    enddo
   ! Zboson-lepton to lepton vertices
    do i=1,3
       call append_vertex(this,l,23,[23,(11+(2*i-2)),(11+(2*i-2))],z_fermion_coupling(this,11+(2*i-2)))
    enddo
    ! Zboson-antilepton to antilepton vertices
    do i=1,3
       call append_vertex(this,l,24,[23,-(11+(2*i-2)),-(11+(2*i-2))],z_fermion_coupling(this,-(11+(2*i-2))))
    enddo

    this%nint=l
    write (99,*) l,'interactions loaded'
  end subroutine init_vert

  integer function get_antipart(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    i=find_particle_index(this,ipdg)
    if (i.gt.0 .and. this%particle_list(i)%type.eq.ipdg) then
       get_antipart=this%particle_list(i)%anti_type
       return
    elseif (i.gt.0) then
       get_antipart=this%particle_list(i)%type
       return
    endif
    write (*,*) 'Particle not in model (mass)',ipdg
    stop 1
  end function get_antipart
  real(kind=8) function get_mass(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       get_mass=this%particle_list(i)%mass
       return
    endif
    write (*,*) 'Particle not in model (mass)',ipdg
    stop 1
  end function get_mass
  real(kind=8) function get_width(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       get_width=this%particle_list(i)%width
       return
    endif
    write (*,*) 'Particle not in model (width)',ipdg
    stop 1
  end function get_width
  integer function get_spin(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       get_spin=this%particle_list(i)%spin
       if (get_spin.lt.0) then
          write (*,*) 'Spin ill-defined for particle',ipdg
          stop 1
       endif
       return
    endif
    write (*,*) 'Particle not in model (spin)',ipdg
    stop 1
  end function get_spin
  integer function get_dim(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       get_dim=this%particle_list(i)%dim
       return
    endif
    write (*,*) 'Particle not in model (dim)',ipdg
    stop 1
  end function get_dim
  real(kind=8) function get_charge(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg,sgn
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       sgn=particle_property_sign(this,ipdg)
       get_charge=dble(sgn)*this%particle_list(i)%charge
       return
    endif
    write (*,*) 'Particle not in model (charge)',ipdg
    stop 1
  end function get_charge
  real(kind=8) function get_isospin_l(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg,sgn
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       sgn=particle_property_sign(this,ipdg)
       get_isospin_l=dble(sgn)*this%particle_list(i)%weak_isospin(1)
       return
    endif
    write (*,*) 'Particle not in model (left isospin)',ipdg
    stop 1
  end function get_isospin_l
  real(kind=8) function get_isospin_r(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg,sgn
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       sgn=particle_property_sign(this,ipdg)
       get_isospin_r=dble(sgn)*this%particle_list(i)%weak_isospin(2)
       return
    endif
    write (*,*) 'Particle not in model (right isospin)',ipdg
    stop 1
  end function get_isospin_r
  real(kind=8) function get_hypercharge_l(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg,sgn
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       sgn=particle_property_sign(this,ipdg)
       get_hypercharge_l=dble(sgn)*this%particle_list(i)%weak_hypercharge(1)
       return
    endif
    write (*,*) 'Particle not in model (left hypercharge)',ipdg
    stop 1
  end function get_hypercharge_l
  real(kind=8) function get_hypercharge_r(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg,sgn
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       sgn=particle_property_sign(this,ipdg)
       get_hypercharge_r=dble(sgn)*this%particle_list(i)%weak_hypercharge(2)
       return
    endif
    write (*,*) 'Particle not in model (right hypercharge)',ipdg
    stop 1
  end function get_hypercharge_r
  integer function get_color_rep(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: i,ipdg,sgn
    i=find_particle_index(this,ipdg)
    if (i.gt.0) then
       sgn=particle_property_sign(this,ipdg)
       get_color_rep=this%particle_list(i)%color_rep
       if (sgn.lt.0 .and. abs(get_color_rep).eq.3) get_color_rep=-get_color_rep
       return
    endif
    write (*,*) 'Particle not in model (color representation)',ipdg
    stop 1
  end function get_color_rep
  integer function get_color_dim(this,ipdg)
    implicit none
    class(physics_model) :: this
    integer :: ipdg
    get_color_dim=abs(this%get_color_rep(ipdg))
  end function get_color_dim
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
    is_quark=iPDG.ge.1 .and. iPDG.le.6
  end function is_quark
  logical function is_antiquark(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_antiquark=iPDG.le.-1 .and. iPDG.ge.-6
  end function is_antiquark
  logical function is_lepton(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_lepton=iPDG.ge.11 .and. iPDG.le.16
  end function is_lepton
  logical function is_antilepton(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_antilepton=iPDG.le.-11 .and. iPDG.ge.-16
  end function is_antilepton
  logical function is_lepton_any(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_lepton_any=abs(iPDG).ge.11 .and. abs(iPDG).le.16
  end function is_lepton_any
  logical function is_gluon(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_gluon=iPDG.eq.21 .or. iPDG.eq.99
  end function is_gluon
  logical function is_scalar(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_scalar=abs(iPDG).eq.25.or.abs(iPDG).eq.125.or.abs(iPDG).eq.126.or.abs(iPDG).eq.127
  end function is_scalar
  logical function is_tensor_g(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_tensor_g=iPDG.eq.-21
  end function is_tensor_g
  logical function is_tensor_z(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_tensor_z=iPDG.eq.-23
  end function is_tensor_z
  logical function is_tensor_w(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_tensor_w=abs(iPDG).eq.26
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
    is_singlet=.not.(abs(iPDG).le.6 .or. iPDG.eq.21)
  end function is_singlet
  logical function is_photon(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_photon=iPDG.eq.22
  end function is_photon
  logical function is_massiveboson(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_massiveboson=iPDG.eq.23 .or. abs(iPDG).eq.24
  end function is_massiveboson
  logical function is_higgs(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_higgs=iPDG.eq.25
  end function is_higgs
  logical function is_higgsor(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_higgsor=iPDG.eq.125.or.iPDG.eq.126.or.iPDG.eq.127
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
