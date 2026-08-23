module particles
  use run_parameters, only: sw,configured_flavour_scheme=>flavour_scheme,&
       configured_up_mass=>up_mass,configured_strange_mass=>strange_mass,&
       configured_charm_mass=>charm_mass,configured_bottom_mass=>bottom_mass,&
       configured_top_mass=>top_mass,&
       configured_top_width=>top_width,configured_z_mass=>z_mass,&
       configured_z_width=>z_width,configured_w_mass=>w_mass,&
       configured_w_width=>w_width,configured_higgs_mass=>higgs_mass,&
       configured_higgs_width=>higgs_width,ignore_final_state_width_fix
  implicit none
  integer,parameter :: model_particle_capacity = 24
  integer,parameter :: model_vertex_capacity = 256
  integer,parameter,public :: model_signature_size = 15
  private :: append_particle, append_vertex, find_particle_index, particle_property_sign&
       &, weak_coupling, weak_cosine, neutral_gauge_coupling&
       &, charged_current_coupling, particle_species_index&
       &, photon_fermion_coupling, z_fermion_coupling&
       &, weak_coupling_squared, weak_cosine_squared&
       &, weak_coupling_over_cosine, higgs_self_coupling&
       &, model_particle_capacity, model_vertex_capacity,sw&
       &, configured_flavour_scheme,configured_up_mass,configured_strange_mass&
       &, configured_charm_mass,configured_bottom_mass&
       &, configured_top_mass,configured_top_width,configured_z_mass&
       &, configured_z_width,configured_w_mass,configured_w_width&
       &, configured_higgs_mass,configured_higgs_width&
       &, ignore_final_state_width_fix
  type particle
     ! weak_isospin and weak_hypercharge use index 1=left, 2=right with Q = T3 + Y/2.
     integer :: type,anti_type,spin,dim,color_rep
     real(kind=8) :: mass,width,charge
     real(kind=8),dimension(2) :: weak_isospin,weak_hypercharge
  end type particle
  type vertex
     integer :: type,n_inputs
     integer,dimension(4) :: particles
     integer :: qcd_power,ew_power,heft_power
     real(kind=8),dimension(2) :: coupl
  end type vertex
  type physics_model
     type(particle),dimension(:),allocatable :: particle_list
     type(vertex),dimension(:),allocatable :: vertex_list
     integer :: npart,nint,flavour_scheme=5
     logical :: heft_enabled=.false.
   contains
     procedure,public :: init_part,get_mass,get_width,get_spin&
          &,get_antipart,init_vert,get_dim,get_inter_dim,is_quark&
          &,get_charge,get_isospin_l,get_isospin_r,get_hypercharge_l&
          &,get_hypercharge_r,get_color_rep,get_color_dim&
          &,is_antiquark,is_lepton,is_antilepton,is_lepton_any,&
          &is_gluon,is_colour_flow_vector,is_gluon_aux_tensor,is_z_aux_tensor&
          &,is_w_aux_tensor,is_auxiliary_tensor,is_singlet,is_photon&
          &,is_massive_vector,is_higgs,is_jet,is_auxiliary_scalar,is_fermion,is_massless_fermion&
          &,is_chiral_eligible,set_width,apply_final_state_widths&
          &,model_signature,set_heft_enabled
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

  subroutine append_vertex(this,l,itype,ipdgs,coupl,qcd_power,ew_power,heft_power)
    implicit none
    class(physics_model),intent(inout) :: this
    integer,intent(inout) :: l
    integer,intent(in) :: itype
    integer,dimension(:),intent(in) :: ipdgs
    real(kind=8),dimension(2),intent(in) :: coupl
    integer,intent(in),optional :: qcd_power,ew_power,heft_power
    l=l+1
    if (l.gt.this%nint) then
       write (*,*) 'ERROR: more vertices than allocated',l,this%nint
       stop 1
    endif
    this%vertex_list(l)%type=itype
    this%vertex_list(l)%n_inputs=size(ipdgs)-1
    if (this%vertex_list(l)%n_inputs.lt.2 .or. this%vertex_list(l)%n_inputs.gt.3) then
       write (*,*) 'ERROR: unsupported vertex arity',this%vertex_list(l)%n_inputs
       stop 1
    endif
    this%vertex_list(l)%particles=0
    this%vertex_list(l)%particles(1:size(ipdgs))=ipdgs
    this%vertex_list(l)%coupl=coupl
    this%vertex_list(l)%qcd_power=0
    this%vertex_list(l)%ew_power=0
    this%vertex_list(l)%heft_power=0
    if (itype.ge.0 .and. itype.le.9) then
       this%vertex_list(l)%qcd_power=1
    elseif (itype.ge.10 .and. itype.le.24) then
       this%vertex_list(l)%ew_power=1
    endif
    if (present(qcd_power)) this%vertex_list(l)%qcd_power=qcd_power
    if (present(ew_power)) this%vertex_list(l)%ew_power=ew_power
    if (present(heft_power)) this%vertex_list(l)%heft_power=heft_power
  end subroutine append_vertex

  subroutine set_heft_enabled(this,enabled)
    implicit none
    class(physics_model),intent(inout) :: this
    logical,intent(in) :: enabled
    this%heft_enabled=enabled
  end subroutine set_heft_enabled

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
    class(physics_model),intent(inout) :: this
    integer :: i,l
    real(kind=8),intent(in),optional :: tmass,twidth,zmass,zwidth,wmass,wwidth,hmass,hwidth
    real(kind=8) :: tmass_local,twidth_local,zmass_local,zwidth_local
    real(kind=8) :: wmass_local,wwidth_local,hmass_local,hwidth_local
    real(kind=8),dimension(5) :: configured_quark_masses
    tmass_local=configured_top_mass
    twidth_local=configured_top_width
    zmass_local=configured_z_mass
    zwidth_local=configured_z_width
    wmass_local=configured_w_mass
    wwidth_local=configured_w_width
    hmass_local=configured_higgs_mass
    hwidth_local=configured_higgs_width
    if (present(tmass)) tmass_local=tmass
    if (present(twidth)) twidth_local=twidth
    if (present(zmass)) zmass_local=zmass
    if (present(zwidth)) zwidth_local=zwidth
    if (present(wmass)) wmass_local=wmass
    if (present(wwidth)) wwidth_local=wwidth
    if (present(hmass)) hmass_local=hmass
    if (present(hwidth)) hwidth_local=hwidth
    this%flavour_scheme=configured_flavour_scheme
    configured_quark_masses=[0d0,configured_up_mass,configured_strange_mass,&
         configured_charm_mass,configured_bottom_mass]
    l=0
    if (allocated(this%particle_list)) deallocate(this%particle_list)
    if (allocated(this%vertex_list)) deallocate(this%vertex_list)
    this%npart=model_particle_capacity ! gluon, quarks, tensors, bosons, Higgs and leptons
    allocate(this%particle_list(this%npart))

    ! The active flavour-scheme quarks are exactly massless.  Every heavier
    ! flavour uses its configured kinematic mass, which also switches on its
    ! Higgs Yukawa vertex in init_vert().
    do i=1,5
       if (mod(i,2).eq.0) then
          call append_particle(this,l,i,merge(0d0,configured_quark_masses(i),&
               i.le.this%flavour_scheme),0d0,2,-i,4,2d0/3d0,&
               [0.5d0,0d0],[1d0/3d0,4d0/3d0],3)
       else
          call append_particle(this,l,i,merge(0d0,configured_quark_masses(i),&
               i.le.this%flavour_scheme),0d0,2,-i,4,-1d0/3d0,&
               [-0.5d0,0d0],[1d0/3d0,-2d0/3d0],3)
       endif
    enddo

    ! top quark
    call append_particle(this,l,6,tmass_local,twidth_local,2,-6,4,2d0/3d0,&
         [0.5d0,0d0],[1d0/3d0,4d0/3d0],3)

    ! gluon
    call append_particle(this,l,21,0d0,0d0,2,21,4,0d0,[0d0,0d0],[0d0,0d0],8)

    ! Auxiliary colour-flow U(1) vector
    call append_particle(this,l,99,0d0,0d0,-1,99,4,0d0,[0d0,0d0],[0d0,0d0],1)

    ! Gluon auxiliary tensor (non-propagating field for the four-gluon interaction)
    call append_particle(this,l,-21,0d0,0d0,-1,-21,6,0d0,[0d0,0d0],[0d0,0d0],8)

    ! photon
    call append_particle(this,l,22,0d0,0d0,2,22,4,0d0,[0d0,0d0],[0d0,0d0],1)

    ! Z-boson
    call append_particle(this,l,23,zmass_local,zwidth_local,3,23,4,0d0,&
         [0d0,0d0],[0d0,0d0],1)

    ! Neutral auxiliary tensor (non-propagating field for four-vector interactions)
    call append_particle(this,l,-23,0d0,0d0,-1,-23,6,0d0,[0d0,0d0],[0d0,0d0],1)

    ! W-boson
    call append_particle(this,l,24,wmass_local,wwidth_local,3,-24,4,1d0,&
         [1d0,1d0],[0d0,0d0],1)

    ! Higgs-boson
    call append_particle(this,l,25,hmass_local,hwidth_local,1,25,1,0d0,&
         [-0.5d0,-0.5d0],[1d0,1d0],1)

    ! Auxiliary scalar A (non-propagating field for four-boson interactions)
    call append_particle(this,l,125,0d0,0d0,-1,125,1,0d0,[0d0,0d0],[0d0,0d0],1)
    ! Auxiliary scalar B (non-propagating field for four-boson interactions)
    call append_particle(this,l,126,0d0,0d0,-1,126,1,0d0,[0d0,0d0],[0d0,0d0],1)
    ! Auxiliary scalar C (non-propagating field for four-boson interactions)
    call append_particle(this,l,127,0d0,0d0,-1,127,1,0d0,[0d0,0d0],[0d0,0d0],1)

    ! Charged auxiliary tensor (non-propagating field for four-vector interactions)
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
    if (allocated(this%vertex_list)) deallocate(this%vertex_list)
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
    ! Antiquark-gluon to antiquark vertices
    do i=1,6
       call append_vertex(this,l,7,[-i,21,-i],[1d0,0d0])
    enddo
    ! antiquark-quark to gluon vertices
    do i=1,6
       call append_vertex(this,l,9,[-i,i,21],[1d0,0d0])
    enddo
    ! Quark-antiquark to auxiliary colour-flow U(1) vector vertices
    do i=1,6
       call append_vertex(this,l,8,[i,-i,99],[1d0/3d0,0d0])  ! 1/N_C coupling
    enddo
    ! Auxiliary colour-flow U(1) vector-quark to quark vertices
    do i=1,6
       call append_vertex(this,l,4,[99,i,i],[1d0,0d0])
    enddo
    ! Auxiliary colour-flow U(1) vector-antiquark to antiquark vertices
    do i=1,6
       call append_vertex(this,l,5,[99,-i,-i],[1d0,0d0])
    enddo
    ! Quark-auxiliary colour-flow U(1) vector to quark vertices
    do i=1,6
       call append_vertex(this,l,6,[i,99,i],[1d0,0d0])
    enddo
    ! Antiquark-auxiliary colour-flow U(1) vector to antiquark vertices
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


    ! W-boson-W-boson to neutral auxiliary tensor
    call append_vertex(this,l,13,[24,-24,-23],[weak_coupling(),0d0]) !!
    call append_vertex(this,l,13,[-24,24,-23],[-weak_coupling(),0d0]) !!
    ! Neutral auxiliary tensor-W-boson to W-boson
    call append_vertex(this,l,14,[-23,24,24],[weak_coupling(),0d0]) !!
    call append_vertex(this,l,14,[-23,-24,-24],[-weak_coupling(),0d0]) !!
    ! W-boson-neutral auxiliary tensor to W-boson
    call append_vertex(this,l,15,[24,-23,24],[-weak_coupling(),0d0]) !!
    call append_vertex(this,l,15,[-24,-23,-24],[weak_coupling(),0d0]) !!


    ! W-boson-photon to charged auxiliary tensor
    call append_vertex(this,l,13,[24,22,26],[this%get_charge(24),0d0])
    call append_vertex(this,l,13,[-24,22,-26],[-this%get_charge(-24),0d0])
    ! Photon-W-boson to charged auxiliary tensor
    call append_vertex(this,l,13,[22,24,26],[-this%get_charge(24),0d0])
    call append_vertex(this,l,13,[22,-24,-26],[this%get_charge(-24),0d0])
    ! W-boson-Z-boson to charged auxiliary tensor
    call append_vertex(this,l,13,[24,23,26],[neutral_gauge_coupling(),0d0]) !!
    call append_vertex(this,l,13,[-24,23,-26],[neutral_gauge_coupling(),0d0]) !!
    ! Z-boson-W-boson to charged auxiliary tensor
    call append_vertex(this,l,13,[23,24,26],[-neutral_gauge_coupling(),0d0]) !!
    call append_vertex(this,l,13,[23,-24,-26],[-neutral_gauge_coupling(),0d0]) !!


    ! Charged auxiliary tensor-photon to W-boson
    call append_vertex(this,l,14,[26,22,24],[this%get_charge(26),0d0])
    call append_vertex(this,l,14,[-26,22,-24],[-this%get_charge(-26),0d0])
    ! Charged auxiliary tensor-W-boson to photon
    call append_vertex(this,l,14,[26,-24,22],[-this%get_charge(26),0d0])
    call append_vertex(this,l,14,[-26,24,22],[this%get_charge(-26),0d0])
    ! Charged auxiliary tensor-Z-boson to W-boson
    call append_vertex(this,l,14,[26,23,24],[neutral_gauge_coupling(),0d0])
    call append_vertex(this,l,14,[-26,23,-24],[neutral_gauge_coupling(),0d0])
    ! Charged auxiliary tensor-W-boson to Z-boson
    call append_vertex(this,l,14,[26,-24,23],[-neutral_gauge_coupling(),0d0]) !!
    call append_vertex(this,l,14,[-26,24,23],[-neutral_gauge_coupling(),0d0]) !!
    ! Photon-charged auxiliary tensor to W-boson
    call append_vertex(this,l,15,[22,26,24],[-this%get_charge(26),0d0])
    call append_vertex(this,l,15,[22,-26,-24],[this%get_charge(-26),0d0])
    ! W-boson-charged auxiliary tensor to photon
    call append_vertex(this,l,15,[24,-26,22],[this%get_charge(24),0d0])
    call append_vertex(this,l,15,[-24,26,22],[-this%get_charge(-24),0d0])
    ! Z-boson-charged auxiliary tensor to W-boson
    call append_vertex(this,l,15,[23,26,24],[-neutral_gauge_coupling(),0d0])
    call append_vertex(this,l,15,[23,-26,-24],[-neutral_gauge_coupling(),0d0])
    ! W-boson-charged auxiliary tensor to Z-boson
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

    ! W-boson-W-boson to auxiliary scalar C
    call append_vertex(this,l,17,[24,-24,127],[1d0/2d0*weak_coupling_squared(),0d0],ew_power=2) !!
    call append_vertex(this,l,17,[-24,24,127],[1d0/2d0*weak_coupling_squared(),0d0],ew_power=2) !!
    ! Z-boson-Z-boson to auxiliary scalar C
    call append_vertex(this,l,17,[23,23,127],[1d0/2d0*weak_coupling_squared()/weak_cosine_squared(),0d0],&
         ew_power=2) !!
    ! Auxiliary scalar A-Higgs to Higgs
    call append_vertex(this,l,20,[125,25,25],[1d0,0d0],ew_power=0) !!
    ! Higgs-auxiliary scalar A to Higgs
    call append_vertex(this,l,20,[25,125,25],[1d0,0d0],ew_power=0) !!
    ! Higgs-Higgs to auxiliary scalar A
    call append_vertex(this,l,20,[25,25,125],[(-3d0/4d0)*weak_coupling_squared()*&
         this%get_mass(25)**2/this%get_mass(24)**2,0d0],ew_power=2) !!
    ! Auxiliary scalar C-Higgs to Higgs
    call append_vertex(this,l,20,[127,25,25],[1d0,0d0],ew_power=0) !!
    ! Higgs-auxiliary scalar C to Higgs
    call append_vertex(this,l,20,[25,127,25],[1d0,0d0],ew_power=0) !!
    ! Higgs-Higgs to auxiliary scalar B
    call append_vertex(this,l,20,[25,25,126],[1d0,-10d0],ew_power=0) !!

    ! Auxiliary scalar B-Z-boson to Z-boson
    call append_vertex(this,l,18,[126,23,23],[-1d0/2d0*weak_coupling_squared()/&
         weak_cosine_squared(),0d0],ew_power=2) !!
    ! Z-boson-auxiliary scalar B to Z-boson
    call append_vertex(this,l,19,[23,126,23],[-1d0/2d0*weak_coupling_squared()/&
         weak_cosine_squared(),0d0],ew_power=2) !!

    ! Auxiliary scalar B-W-boson to W-boson
    call append_vertex(this,l,18,[126,24,24],[-1d0/2d0*weak_coupling_squared(),0d0],ew_power=2) !!
    call append_vertex(this,l,18,[126,-24,-24],[-1d0/2d0*weak_coupling_squared(),0d0],ew_power=2) !!
    ! W-boson-auxiliary scalar B to W-boson
    call append_vertex(this,l,19,[24,126,24],[-1d0/2d0*weak_coupling_squared(),0d0],ew_power=2) !!
    call append_vertex(this,l,19,[-24,126,-24],[-1d0/2d0*weak_coupling_squared(),0d0],ew_power=2) !!

    ! Lepton-antilepton to photon
    do i=1,3
    call append_vertex(this,l,21,[11+(2*i-2),-(11+(2*i-2)),22],photon_fermion_coupling(this,11+(2*i-2)))
    enddo
    ! Antilepton-lepton to photon
    do i=1,3
    call append_vertex(this,l,22,[-(11+(2*i-2)),11+(2*i-2),22],photon_fermion_coupling(this,-(11+(2*i-2))))
    enddo
    ! Lepton-antilepton to Z-boson
    do i=1,3
    call append_vertex(this,l,21,[(11+(2*i-2)),-(11+(2*i-2)),23],z_fermion_coupling(this,11+(2*i-2)))
    enddo
    ! Antilepton-lepton to Z-boson
    do i=1,3
    call append_vertex(this,l,22,[-(11+(2*i-2)),(11+(2*i-2)),23],z_fermion_coupling(this,-(11+(2*i-2))))
    enddo
    ! Neutrino-antineutrino and antineutrino-neutrino to Z-boson
    do i=1,3
       call append_vertex(this,l,21,[(12+(2*i-2)),-(12+(2*i-2)),23],&
            z_fermion_coupling(this,12+(2*i-2)))
       call append_vertex(this,l,22,[-(12+(2*i-2)),(12+(2*i-2)),23],&
            z_fermion_coupling(this,-(12+(2*i-2))))
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
    ! Neutrino-Z and Z-neutrino currents (and their antiparticles)
    do i=1,3
       call append_vertex(this,l,10,[(12+(2*i-2)),23,(12+(2*i-2))],&
            z_fermion_coupling(this,12+(2*i-2)))
       call append_vertex(this,l,11,[-(12+(2*i-2)),23,-(12+(2*i-2))],&
            z_fermion_coupling(this,-(12+(2*i-2))))
       call append_vertex(this,l,23,[23,(12+(2*i-2)),(12+(2*i-2))],&
            z_fermion_coupling(this,12+(2*i-2)))
       call append_vertex(this,l,24,[23,-(12+(2*i-2)),-(12+(2*i-2))],&
            z_fermion_coupling(this,-(12+(2*i-2))))
    enddo

    if (this%heft_enabled) then
       ! CP-even HEFT: -g_H/4 h G^a_{mu nu} G^{a,mu nu}.
       ! Coupling powers are applied to each recursive interaction at runtime.
       call append_vertex(this,l,25,[21,21,25],[1d0,0d0],&
            qcd_power=0,ew_power=0,heft_power=1)
       call append_vertex(this,l,26,[25,21,21],[1d0,0d0],&
            qcd_power=0,ew_power=0,heft_power=1)
       call append_vertex(this,l,27,[21,25,21],[1d0,0d0],&
            qcd_power=0,ew_power=0,heft_power=1)

       ! Ordered h-g-g-g contact interaction.  Singlets are canonicalised
       ! after the coloured children, so only the [g,g,h]->g crossing is
       ! needed for an off-shell gluon; adding the other scalar placements
       ! would count the same contact diagram three times.
       call append_vertex(this,l,28,[21,21,21,25],[1d0,0d0],&
            qcd_power=1,ew_power=0,heft_power=1)
       call append_vertex(this,l,31,[21,21,25,21],[1d0,0d0],&
            qcd_power=1,ew_power=0,heft_power=1)

       ! The h-g-g-g-g contact is represented with the existing compact
       ! antisymmetric auxiliary gluon tensor.
       call append_vertex(this,l,32,[-21,-21,25],[1d0,0d0],&
            qcd_power=0,ew_power=0,heft_power=1)
       call append_vertex(this,l,33,[25,-21,-21],[1d0,0d0],&
            qcd_power=0,ew_power=0,heft_power=1)
       call append_vertex(this,l,34,[-21,25,-21],[1d0,0d0],&
            qcd_power=0,ew_power=0,heft_power=1)
    endif

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

  subroutine set_width(this,ipdg,new_width)
    implicit none
    class(physics_model),intent(inout) :: this
    integer,intent(in) :: ipdg
    real(kind=8),intent(in) :: new_width
    integer :: i
    i=find_particle_index(this,ipdg)
    if (i.le.0) then
       write (*,*) 'Particle not in model (set width)',ipdg
       stop 1
    endif
    if (new_width.lt.0d0) then
       write (*,*) 'Cannot assign a negative particle width',ipdg,new_width
       stop 1
    endif
    this%particle_list(i)%width=new_width
  end subroutine set_width

  subroutine apply_final_state_widths(this,n,nprocs,processes)
    implicit none
    class(physics_model),intent(inout) :: this
    integer,intent(in) :: n,nprocs
    integer,dimension(n,nprocs),intent(in) :: processes
    integer :: i,iproc,ipdg
    real(kind=8) :: old_width
    if (ignore_final_state_width_fix) return
    do iproc=1,nprocs
       do i=3,n
          ipdg=processes(i,iproc)
          if (this%get_mass(ipdg).le.0d0) cycle
          old_width=this%get_width(ipdg)
          if (old_width.eq.0d0) cycle
          write (*,'(a,i8,a,es14.6)') 'Setting final-state particle width to zero for PDG ',&
               ipdg,'; configured width = ',old_width
          call this%set_width(ipdg,0d0)
       enddo
    enddo
  end subroutine apply_final_state_widths

  function model_signature(this) result(signature)
    implicit none
    class(physics_model),intent(in) :: this
    integer :: i
    real(kind=8),dimension(model_signature_size) :: signature
    signature=[dble(this%flavour_scheme),sw,(this%get_mass(i),i=1,6),&
         this%get_width(6),this%get_mass(23),&
         this%get_width(23),this%get_mass(24),this%get_width(24),&
         this%get_mass(25),this%get_width(25)]
  end function model_signature
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
          get_inter_dim=this%get_dim(this%vertex_list(i)%particles(&
               this%vertex_list(i)%n_inputs+1))
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
  logical function is_fermion(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_fermion=this%is_quark(iPDG).or.this%is_antiquark(iPDG).or.&
         this%is_lepton(iPDG).or.this%is_antilepton(iPDG)
  end function is_fermion
  logical function is_massless_fermion(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_massless_fermion=this%is_fermion(iPDG).and.this%get_mass(iPDG).eq.0d0
  end function is_massless_fermion
  logical function is_chiral_eligible(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_chiral_eligible=this%is_massless_fermion(iPDG)
  end function is_chiral_eligible
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
    is_gluon=iPDG.eq.21
  end function is_gluon
  logical function is_colour_flow_vector(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_colour_flow_vector=iPDG.eq.21 .or. iPDG.eq.99
  end function is_colour_flow_vector
  logical function is_scalar(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_scalar=abs(iPDG).eq.25.or.abs(iPDG).eq.125.or.abs(iPDG).eq.126.or.abs(iPDG).eq.127
  end function is_scalar
  logical function is_gluon_aux_tensor(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_gluon_aux_tensor=iPDG.eq.-21
  end function is_gluon_aux_tensor
  logical function is_z_aux_tensor(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_z_aux_tensor=iPDG.eq.-23
  end function is_z_aux_tensor
  logical function is_w_aux_tensor(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_w_aux_tensor=abs(iPDG).eq.26
  end function is_w_aux_tensor
  logical function is_auxiliary_tensor(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_auxiliary_tensor=this%is_gluon_aux_tensor(iPDG) .or. &
         this%is_z_aux_tensor(iPDG) .or. this%is_w_aux_tensor(iPDG)
  end function is_auxiliary_tensor
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
  logical function is_massive_vector(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_massive_vector=iPDG.eq.23 .or. abs(iPDG).eq.24
  end function is_massive_vector
  logical function is_higgs(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_higgs=iPDG.eq.25
  end function is_higgs
  logical function is_auxiliary_scalar(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    is_auxiliary_scalar=iPDG.eq.125.or.iPDG.eq.126.or.iPDG.eq.127
  end function is_auxiliary_scalar
  logical function is_jet(this,iPDG)
    implicit none
    class(physics_model) :: this
    integer :: iPDG
    ! A resolved heavy-flavour quark is still a jet for cuts and dynamical
    ! scales.  Beam/inclusive-jet membership is controlled separately by the
    ! process-list flavour scheme; top quarks are deliberately excluded here.
    is_jet=(abs(iPDG).ge.1 .and. abs(iPDG).le.5) .or. this%is_gluon(iPDG)
  end function is_jet
end module particles
