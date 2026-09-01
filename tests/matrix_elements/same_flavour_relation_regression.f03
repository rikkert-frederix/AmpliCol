program same_flavour_relation_regression
  use handling_processes
  implicit none
  integer,parameter :: dp=kind(1d0)
  integer,dimension(4,3) :: two_line_qcd_processes
  integer,dimension(5,3) :: two_line_charged_current_processes
  integer,dimension(6,4) :: three_line_qcd_processes
  integer,dimension(6,7) :: six_term_qcd_processes
  integer,dimension(7,3) :: mixed_three_line_charged_current_processes
  integer,dimension(7,4) :: three_line_charged_current_processes
  integer,dimension(8,3) :: leptonic_three_line_charged_current_processes

  open(unit=99,status='scratch',action='write')
  call phys_model%init_part(173d0,1.4915d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0.006382339d0)
  call phys_model%init_vert()
  keep_processes_separate=.false.

  two_line_qcd_processes(:,1)=[1,3,-1,-3]
  two_line_qcd_processes(:,2)=[1,3,-3,-1]
  two_line_qcd_processes(:,3)=[1,1,-1,-1]
  call check_relation(two_line_qcd_processes,2,'two-line QCD')

  ! Each different-flavour process selects one of the two possible fermion
  ! flows.  A charged current preserves the generation along a line while it
  ! may change weak isospin, so exact quark-PDG matching would miss this sum.
  two_line_charged_current_processes(:,1)=[1,3,-2,-3,24]
  two_line_charged_current_processes(:,2)=[3,1,-2,-3,24]
  two_line_charged_current_processes(:,3)=[1,1,-2,-1,24]
  call check_relation(two_line_charged_current_processes,2,&
       'two-line charged current')

  three_line_qcd_processes(:,1)=[1,1,3,-1,-1,-3]
  three_line_qcd_processes(:,2)=[1,1,3,-1,-3,-1]
  three_line_qcd_processes(:,3)=[1,1,3,-3,-1,-1]
  three_line_qcd_processes(:,4)=[1,1,1,-1,-1,-1]
  call check_relation(three_line_qcd_processes,3,'three-line QCD')

  ! With three distinct flavours each process isolates one of the 3! flows,
  ! exercising the maximum-size cover rather than only the common 2+2+2 one.
  six_term_qcd_processes(:,1)=[1,2,3,-1,-2,-3]
  six_term_qcd_processes(:,2)=[1,2,3,-1,-3,-2]
  six_term_qcd_processes(:,3)=[1,2,3,-2,-1,-3]
  six_term_qcd_processes(:,4)=[1,2,3,-2,-3,-1]
  six_term_qcd_processes(:,5)=[1,2,3,-3,-1,-2]
  six_term_qcd_processes(:,6)=[1,2,3,-3,-2,-1]
  six_term_qcd_processes(:,7)=[1,1,1,-1,-1,-1]
  call check_relation(six_term_qcd_processes,6,'six-term three-line QCD')

  ! Only four of the six generation-compatible flows contain the single
  ! isospin flip permitted by one external W.  Treating all six as possible
  ! hides this otherwise exact two-term relation.
  mixed_three_line_charged_current_processes(:,1)=[1,3,2,-2,-2,-3,24]
  mixed_three_line_charged_current_processes(:,2)=[3,1,2,-2,-2,-3,24]
  mixed_three_line_charged_current_processes(:,3)=[1,1,2,-2,-2,-1,24]
  call check_relation(mixed_three_line_charged_current_processes,2,&
       'mixed three-line charged current')

  three_line_charged_current_processes(:,1)=[1,1,3,-1,-2,-3,24]
  three_line_charged_current_processes(:,2)=[1,3,1,-1,-2,-3,24]
  three_line_charged_current_processes(:,3)=[3,1,1,-1,-2,-3,24]
  three_line_charged_current_processes(:,4)=[1,1,1,-1,-2,-1,24]
  call check_relation(three_line_charged_current_processes,3,&
       'three-line charged current')

  leptonic_three_line_charged_current_processes(:,1)=&
       [1,3,2,-2,-2,-3,-11,12]
  leptonic_three_line_charged_current_processes(:,2)=&
       [3,1,2,-2,-2,-3,-11,12]
  leptonic_three_line_charged_current_processes(:,3)=&
       [1,1,2,-2,-2,-1,-11,12]
  call check_relation(leptonic_three_line_charged_current_processes,2,&
       'resolved-lepton three-line charged current')

  write (*,*) 'Same-flavour relation regression passed'
contains
  subroutine check_relation(outgoing_processes,expected_terms,label)
    implicit none
    integer,dimension(:,:),intent(in) :: outgoing_processes
    integer,intent(in) :: expected_terms
    character(len=*),intent(in) :: label
    type(phase_space_order_group) :: group
    integer :: n,nprocesses,event,iproc,iamp,nhel
    integer,dimension(:),allocatable :: hel,include_hel
    real(kind=dp),dimension(:,:),allocatable :: momenta
    real(kind=dp),dimension(:),allocatable :: amp2
    complex(kind=dp),dimension(:),allocatable :: direct_amplitudes
    real(kind=dp) :: error,scale

    n=size(outgoing_processes,1)
    nprocesses=size(outgoing_processes,2)
    next=n
    group%next=n
    group%nproc=nprocesses
    allocate(group%processes(n,nprocesses))
    allocate(group%color_orders(n,nprocesses))
    allocate(group%amps(1),group%passed(1),hel(n),amp2(nprocesses))
    group%processes=outgoing_processes
    call setup_spin(group)
    call setup_color_order(group)
    do iproc=1,nprocesses
       group%processes(1,iproc)=phys_model%get_antipart(group%processes(1,iproc))
       group%processes(2,iproc)=phys_model%get_antipart(group%processes(2,iproc))
    enddo
    call group%amps(1)%init(1,n,nprocesses,group%processes,&
         group%spin,group%color_orders,phys_model)
    allocate(momenta(0:3,n))
    hel=0
    group%passed=0
    do event=1,10
       call fill_momenta(group%processes(:,1),event,momenta)
       call group%amps(1)%evaluate(n,momenta,hel,.true.,phys_model)
       call squared_process_amplitudes(group,amp2)
       group%passed(1)=event
       if (event.eq.10) then
          allocate(direct_amplitudes(group%amps(1)%n_amps))
          direct_amplitudes=group%amps(1)%amps
       endif
       call find_same_flavour(group,10,amp2)
    enddo

    if (.not.group%amps(1)%same_flav(nprocesses)) then
       write (*,*) trim(label)//': relation was not found'
       stop 1
    endif
    if (count(group%same_flavour(:,nprocesses).gt.0).ne.expected_terms) then
       write (*,*) trim(label)//': unexpected relation size',&
            group%same_flavour(:,nprocesses)
       stop 1
    endif

    call group%amps(1)%evaluate(n,momenta,hel,.true.,phys_model)
    error=0d0
    scale=0d0
    do iamp=group%amps(1)%iproc_start(nprocesses),&
         group%amps(1)%iproc_start(nprocesses+1)-1
       error=max(error,abs(group%amps(1)%amps(iamp)-direct_amplitudes(iamp)))
       scale=max(scale,abs(group%amps(1)%amps(iamp))+abs(direct_amplitudes(iamp)))
    enddo
    if (scale.eq.0d0 .or. error.gt.1d-10*scale) then
       write (*,*) trim(label)//': reconstructed amplitude differs',error,scale
       stop 1
    endif

    ! Exercise the production helicity-filter remapping as well as direct
    ! reconstruction.  Keeping every helicity makes the expected map the
    ! identity while still traversing every daughter slot, including all six.
    nhel=group%amps(1)%n_amps
    allocate(include_hel(nhel))
    include_hel=1
    call group%amps(1)%filter_helicity(n,nhel,include_hel)
    if (nhel.ne.size(direct_amplitudes)) then
       write (*,*) trim(label)//': helicity filter changed the all-included set'
       stop 1
    endif
    call group%amps(1)%evaluate(n,momenta,hel,.true.,phys_model)
    error=maxval(abs(group%amps(1)%amps-direct_amplitudes))
    scale=maxval(abs(group%amps(1)%amps)+abs(direct_amplitudes))
    if (scale.eq.0d0 .or. error.gt.1d-10*scale) then
       write (*,*) trim(label)//': filtered reconstruction differs',error,scale
       stop 1
    endif
  end subroutine check_relation

  subroutine squared_process_amplitudes(group,amp2)
    implicit none
    type(phase_space_order_group),intent(in) :: group
    real(kind=dp),dimension(group%nproc),intent(out) :: amp2
    integer :: iproc,iamp
    amp2=0d0
    do iproc=1,group%nproc
       do iamp=group%amps(1)%iproc_start(iproc),&
            group%amps(1)%iproc_start(iproc+1)-1
          amp2(iproc)=amp2(iproc)+abs(group%amps(1)%amps(iamp))**2
       enddo
    enddo
  end subroutine squared_process_amplitudes

  subroutine fill_momenta(process,event,p)
    implicit none
    integer,dimension(:),intent(in) :: process
    integer,intent(in) :: event
    real(kind=dp),dimension(0:,:),intent(out) :: p
    integer :: n,i,j,pair_start,ipair
    real(kind=dp) :: q,total_energy,mass_i,mass_j,angle
    real(kind=dp),dimension(3) :: direction
    n=size(process)
    p=0d0
    total_energy=0d0
    if (mod(n-2,2).eq.1) then
       q=120d0+3d0*event
       p(1:3,3)=q*[1d0,0d0,0d0]
       p(1:3,4)=q*[-0.5d0,sqrt(0.75d0),0d0]
       p(1:3,5)=q*[-0.5d0,-sqrt(0.75d0),0d0]
       do i=3,5
          p(0,i)=sqrt(q**2+phys_model%get_mass(process(i))**2)
          total_energy=total_energy+p(0,i)
       enddo
       pair_start=6
    else
       pair_start=3
    endif
    ipair=0
    do i=pair_start,n,2
       j=i+1
       ipair=ipair+1
       angle=0.17d0*event+0.41d0*ipair
       direction=[0.6d0*cos(angle),0.6d0*sin(angle),0.8d0]
       q=80d0+7d0*event+19d0*ipair
       p(1:3,i)=q*direction
       p(1:3,j)=-q*direction
       mass_i=phys_model%get_mass(process(i))
       mass_j=phys_model%get_mass(process(j))
       p(0,i)=sqrt(q**2+mass_i**2)
       p(0,j)=sqrt(q**2+mass_j**2)
       total_energy=total_energy+p(0,i)+p(0,j)
    enddo
    p(:,1)=[0.5d0*total_energy,0d0,0d0,0.5d0*total_energy]
    p(:,2)=[0.5d0*total_energy,0d0,0d0,-0.5d0*total_energy]
  end subroutine fill_momenta
end program same_flavour_relation_regression
