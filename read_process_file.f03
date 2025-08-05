module read_process_file
  use mint_module
  use handling_processes
  integer :: sf_nprocs
  integer,dimension(:,:),allocatable :: unique_procs,processes,color_orders,multi_chans &
       ,proc_sf_map,same_flavour_process,same_flavour_order,same_flavour_ichans,same_flavour_process_map
  type(phase_space_order_group),allocatable :: pgl_unique
  real(kind=8),dimension(:),allocatable :: unique_map_value
  integer,dimension(:),allocatable :: unique_map,iden_iproc,same_flavour_idenCOfactor
  integer,dimension(:,:,:),allocatable :: iden_processes
  real(kind=8),dimension(:,:),allocatable :: idenCOandMAPfactor
contains
  subroutine read_processes_from_file(filename)
    implicit none
    character(len=80) :: filename
    integer :: iproc,igroup,icheck,nproc_in_group,max_channels,idenCOfactor
    integer,dimension(:),allocatable :: process,order,ichans,phase_space_orders
    character(len=1024) :: buff
    open(unit=10,file=filename,status='old')
    read (10,*) next,nproc_unique
    ndim=3*(next-2)-4
    if (include_pdf) ndim=ndim+2
    allocate(unique_procs(1:next,1:nproc_unique))
    do iproc=1,nproc_unique
       read(10,*) unique_procs(1:next,iproc)
    enddo
    allocate(pgl_unique)
    pgl_unique%next=next
    call check_unique_processes()
    read(10,*)
    read(10,*)
    read (10,*) ngroups
    allocate(pgl(ngroups))

    allocate(process(1:next))
    allocate(order(1:next))

    read (10,*) 
    do igroup=1,ngroups
       nprocs=0
       sf_nprocs=0
       allocate(phase_space_orders(1:next))
       read(10,*) icheck,nproc_in_group,max_channels,phase_space_orders(1:next)
       if (icheck.ne.igroup) then
          write (*,*) 'ERROR in processes file',icheck,igroup
          stop 1
       endif
       allocate(iden_iproc(nproc_in_group))
       allocate(processes(1:next,nproc_in_group))
       allocate(color_orders(1:next,nproc_in_group))
       allocate(iden_processes(1:next,nproc_in_group,nproc_in_group))
       allocate(idenCOandMAPfactor(nproc_in_group,nproc_in_group))
       allocate(multi_chans(0:max_channels,nproc_in_group))
       allocate(ichans(0:max_channels))
       allocate(same_flavour_process(1:next,nproc_in_group))
       allocate(same_flavour_order(1:next,nproc_in_group))
       allocate(same_flavour_ichans(0:max_channels,nproc_in_group))
       allocate(same_flavour_idenCOfactor(nproc_in_group))
       allocate(same_flavour_process_map(2,nproc_in_group))
       do iproc=1,nproc_in_group
          read(10,'(a)') buff
          read(buff,*) ichans(0)
          read(buff,*) ichans(0),ichans(1:ichans(0)),process(1:next),order(1:next),idenCOfactor
          call add_to_process_list(process,order,idenCOfactor,max_channels,ichans)
       enddo
       call add_same_flavour_processes_to_list(max_channels)
       pgl(igroup)%next=next
       pgl(igroup)%nproc=nprocs
       pgl(igroup)%multichan%max_channels=max_channels
       allocate(pgl(igroup)%processes(1:next,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%color_orders(1:next,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%phase_space_orders(1:next))
       allocate(pgl(igroup)%idenCOandMAPfactor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%iden_iproc(1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%val_procs(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%multichan%channels(1:max_channels,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%multichan%number_of_channels(1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%amps%same_flavour_process_map(2,1:pgl(igroup)%nproc))
       pgl(igroup)%processes(1:next,1:pgl(igroup)%nproc)=processes(1:next,1:pgl(igroup)%nproc)
       pgl(igroup)%color_orders(1:next,1:pgl(igroup)%nproc)=color_orders(1:next,1:pgl(igroup)%nproc)
       pgl(igroup)%phase_space_orders(1:next)=phase_space_orders(1:next)
       pgl(igroup)%idenCOandMAPfactor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)=&
            idenCOandMAPfactor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)
       pgl(igroup)%iden_iproc(1:pgl(igroup)%nproc)=iden_iproc(1:pgl(igroup)%nproc)
       pgl(igroup)%iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)=&
            iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)
       pgl(igroup)%multichan%channels(1:max_channels,1:pgl(igroup)%nproc)=multi_chans(1:max_channels,1:pgl(igroup)%nproc)
       pgl(igroup)%multichan%number_of_channels(1:pgl(igroup)%nproc)=multi_chans(0,1:pgl(igroup)%nproc)
       pgl(igroup)%amps%same_flavour_process_map(1:2,1:pgl(igroup)%nproc)=same_flavour_process_map(1:2,1:pgl(igroup)%nproc)
       deallocate(iden_iproc)
       deallocate(processes)
       deallocate(color_orders)
       deallocate(phase_space_orders)
       deallocate(iden_processes)
       deallocate(idenCOandMAPfactor)
       deallocate(multi_chans)
       deallocate(ichans)
       deallocate(same_flavour_process)
       deallocate(same_flavour_order)
       deallocate(same_flavour_ichans)
       deallocate(same_flavour_idenCOfactor)
       deallocate(same_flavour_process_map)
       write (*,*) '****************************************************'
       do iproc=1,pgl(igroup)%nproc
          write(*,*) iproc,':',pgl(igroup)%processes(1:next,iproc),' ; ',&
               pgl(igroup)%color_orders(1:next,iproc),' ; ',pgl(igroup)%iden_iproc(iproc),' ; ',&
               pgl(igroup)%multichan%channels(1:pgl(igroup)%multichan%number_of_channels(iproc),iproc)
       enddo
       write (*,*) '****************************************************'
       read(10,*)
       read(10,*)
       read(10,*)
    enddo
  end subroutine read_processes_from_file


  
  subroutine check_unique_processes()
    use phase_space_gen23_mod
    use cuts
    implicit none
    integer :: i,iproc,ievent,ih,nqq
    integer,parameter :: nevent=10
    real(kind=8),dimension(:,:),allocatable :: amp2
    complex(kind=8),dimension(:,:),allocatable :: amp
    real(kind=8),dimension(:),allocatable :: mass,width
    real(kind=8),dimension(ndim) :: x
    real(kind=8),external :: ran2
    allocate(phase_space_gen23 :: pgl_unique%phase_space)
    allocate(pgl_unique%processes(next,nproc_unique))
    allocate(pgl_unique%color_orders(next,nproc_unique))
    allocate(pgl_unique%phase_space_orders(next))
    allocate(mass(next))
    allocate(width(next))
    pgl_unique%nproc=nproc_unique
    pgl_unique%processes=unique_procs
    do i=1,pgl_unique%next
        mass(i)=phys_model%get_mass(pgl_unique%processes(i,1))
        width(i)=phys_model%get_width(pgl_unique%processes(i,1))
        do iproc=2,pgl_unique%nproc
           if ( mass(i).ne.phys_model%get_mass(pgl_unique%processes(i,iproc)) .or. &
                width(i).ne.phys_model%get_width(pgl_unique%processes(i,iproc))) then
              write (*,*) 'masses and widths not compatible among processes'
              stop 1
           endif
        enddo
     enddo
     call setup_spin(pgl_unique)
     call setup_color_order(pgl_unique)

     do iproc=1,pgl_unique%nproc
        do i=1,2
           pgl_unique%processes(i,iproc)=phys_model%get_antipart(pgl_unique%processes(i,iproc))
        enddo
     enddo

     ! No multi-channel needed to check: simply use the color_orders for the phase-space order
     pgl_unique%phase_space_orders(1:next)=pgl_unique%color_orders(1:next,1)

     call setup_cuts_for_each_particle(pgl_unique)
     call pgl_unique%phase_space%init(sqrts,next,mass,pgl_unique%phase_space_orders,&
          pgl_unique%pt_min,pgl_unique%eta_max,pgl_unique%DR_min,pgl_unique%sqrt_s_min,.false.,include_pdf)

     call pgl_unique%amps%init(1,next,pgl_unique%nproc,pgl_unique%processes,&
             pgl_unique%spin,pgl_unique%color_orders,phys_model,read_proc_from_file)
     
     allocate(amp2(nevent,pgl_unique%nproc))
     allocate(amp(nevent,pgl_unique%amps%n_amps))

     ievent=0
     do while (ievent.lt.nevent)
        do i=1,ndim
           x(i)=ran2()
        enddo
        call pgl_unique%phase_space%generate_momenta(x)
        if (pgl_unique%phase_space%jac.lt.0d0) cycle
        ievent=ievent+1
        call pgl_unique%amps%evaluate(next,pgl_unique%phase_space%p,pgl_unique%hel,read_proc_from_file,phys_model)
        iproc=0
        amp2(ievent,:)=0d0
        if (use_real_gluons .and. all(pgl_unique%amps%n_qqbar(1:pgl_unique%amps%nprocs).eq.0)) then
           do ih=1,pgl_unique%amps%n_amps
              do while (pgl_unique%amps%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
              amp2(ievent,iproc)=amp2(ievent,iproc)+pgl_unique%amps%amps_r(ih)*pgl_unique%amps%amps_r(ih)
           enddo
        else
           do ih=1,pgl_unique%amps%n_amps
              do while (pgl_unique%amps%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
              amp2(ievent,iproc)=amp2(ievent,iproc)+dble(pgl_unique%amps%amps(ih)*dconjg(pgl_unique%amps%amps(ih)))
           enddo
           amp(ievent,1:pgl_unique%amps%n_amps)=pgl_unique%amps%amps(1:pgl_unique%amps%n_amps)
        endif
     enddo
     allocate(unique_map(1:pgl_unique%nproc))
     allocate(unique_map_value(1:pgl_unique%nproc))
     call find_unique(pgl_unique,nevent,amp2,unique_map,unique_map_value)

     call find_same_flavour(pgl_unique,nevent,amp)
     
     do iproc=1,pgl_unique%nproc
        write (*,*) unique_map(iproc),unique_map_value(iproc),':',pgl_unique%processes(1:pgl_unique%next,iproc),&
                ':',pgl_unique%color_orders(1:pgl_unique%next,iproc)
     enddo
     deallocate(pgl_unique%spin)
     deallocate(pgl_unique%phase_space)
     deallocate(amp2)
     deallocate(amp)
   end subroutine check_unique_processes

   subroutine find_same_flavour(pgl,nevent,amp)
     implicit none
     type(phase_space_order_group),intent(in) :: pgl
     integer,intent(in) :: nevent
     complex(kind=8),dimension(nevent,pgl%amps%n_amps),intent(in) :: amp
     integer :: i,j,k,ii,jj,kk,ievent
     real(kind=8),parameter :: tiny=1d-6
     allocate(proc_sf_map(pgl%nproc,2))
     proc_sf_map=0
     if (.not.decompose_same_flavour_into_two_diff_flavour) return
     do i=1,pgl%nproc
        if (unique_map(i).ne.-1) cycle
        if (pgl_unique%amps%n_qqbar(i).ne.2) cycle
        do j=1,pgl%nproc
           if (i.eq.j) cycle
           if (unique_map(j).ne.-1) cycle
           if (pgl_unique%amps%n_qqbar(j).ne.2) cycle
           do k=1,j-1
              if (k.eq.i) cycle
              if (unique_map(k).ne.-1) cycle
              if (pgl_unique%amps%n_qqbar(k).ne.2) cycle
              do ii=pgl_unique%amps%iproc_start(i),pgl_unique%amps%iproc_start(i+1)-1
                 if (any(amp(1:nevent,ii).eq.(0d0,0d0))) cycle
                 do jj=pgl_unique%amps%iproc_start(j),pgl_unique%amps%iproc_start(j+1)-1
                    if (all(pgl_unique%amps%spins(:,1,ii).eq.pgl_unique%amps%spins(:,1,jj))) exit
                 enddo
                 do kk=pgl_unique%amps%iproc_start(k),pgl_unique%amps%iproc_start(k+1)-1
                    if (all(pgl_unique%amps%spins(:,1,ii).eq.pgl_unique%amps%spins(:,1,kk))) exit
                 enddo
                 if (any((abs(amp(1:nevent,ii)-(amp(1:nevent,jj)+amp(1:nevent,kk)))/&
                          abs(amp(1:nevent,ii)+(amp(1:nevent,jj)+amp(1:nevent,kk)))).gt.tiny)) then
                    exit
                 endif
              enddo
              if (ii.eq.pgl_unique%amps%iproc_start(i+1)) then
                 write (*,'(a,x,i4,x,a,i4,x,a,i4)') &
                      "Found SF amps equal to a sum of DF amps:",i,'=',j,'+',k
                 proc_sf_map(i,1)=j
                 proc_sf_map(i,2)=k
              endif
           enddo
        enddo
     enddo
   end subroutine find_same_flavour

   subroutine find_unique(pgl,nevent,amp2,unique_map,unique_map_value)
     implicit none
     type(phase_space_order_group),intent(in) :: pgl
     integer,intent(in) :: nevent
     real(kind=8),dimension(nevent,pgl%nproc),intent(in) :: amp2
     real(kind=8),dimension(pgl%nproc),intent(out) :: unique_map_value
     integer,dimension(pgl%nproc),intent(out) :: unique_map
     integer :: i,j,n
     real(kind=8),dimension(nevent) :: ratio
     real(kind=8) :: ave
     real(kind=8),parameter :: tiny=1d-6
     unique_map=-1d0
     do i=1,pgl%nproc
        if (all(amp2(1:nevent,i).eq.0d0)) then
           unique_map(i)=0
           unique_map_value(i)=0d0
           cycle
        endif
        do j=1,i-1
           if (all(amp2(1:nevent,j).eq.0d0)) cycle
           ratio(1:nevent)=amp2(1:nevent,i)/amp2(1:nevent,j)
           ave=sum(ratio(1:nevent))/nevent
           if (all(abs(ratio(1:nevent)/ave-1d0).lt.tiny)) then
              unique_map_value(i)=ave
              unique_map(i)=j
              exit
           endif
        enddo
        if (j.eq.i) then
           unique_map(i)=-1
           unique_map_value(i)=1d0
        endif
     enddo
   end subroutine find_unique

   subroutine add_same_flavour_processes_to_list(max_channels)
     implicit none
     integer :: iproc,max_channels
     integer,dimension(next) :: process_mapped,mapping,process_unique
     integer,dimension(4,2) :: quarks_to_map
     real(kind=8) :: idenMAPfactor,idenCOMAPfactor
     logical :: check_sf
     do iproc=1,sf_nprocs
        call map_to_canonical_form(same_flavour_process(1,iproc),process_mapped,mapping)
        call get_unique_process(same_flavour_process(1,iproc),process_mapped,mapping,&
             check_sf,process_unique,idenMAPfactor,quarks_to_map)
        if (.not. check_sf) then
           write (*,*) 'ERROR: current process must be a same-flavour process'
           stop 1
        endif
        idenCOMAPfactor=idenMAPfactor*same_flavour_idenCOfactor(iproc)
        if (idenCOMAPfactor.eq.0d0) cycle
        call add_to_unique_process_list(same_flavour_process(1,iproc),process_unique, &
             same_flavour_order(1,iproc),idenCOMAPfactor,max_channels,same_flavour_ichans(0,iproc),quarks_to_map)
     enddo
   end subroutine add_same_flavour_processes_to_list
  
  subroutine add_to_process_list(process,order,idenCOfactor,max_channels,ichans)
    implicit none
    integer :: max_channels
    integer,dimension(0:max_channels) :: ichans
    integer,dimension(next) :: process,order,process_mapped,process_unique,mapping
    integer,dimension(4,2) :: quarks_to_map
    integer :: idenCOfactor
    real(kind=8) :: idenMAPfactor,idenCOMAPfactor
    logical :: check_sf
    call map_to_canonical_form(process,process_mapped,mapping)
    call get_unique_process(process,process_mapped,mapping,check_sf,process_unique,idenMAPfactor,quarks_to_map)
    if (check_sf) then
       sf_nprocs=sf_nprocs+1
       same_flavour_process(1:next,sf_nprocs)=process(1:next)
       same_flavour_order(1:next,sf_nprocs)=order(1:next)
       same_flavour_ichans(0:ichans(0),sf_nprocs)=ichans(0:ichans(0))
       same_flavour_idenCOfactor(sf_nprocs)=idenCOfactor
       return ! do not add it to the list yet; first finish all the different-flavour processes
    endif
    idenCOMAPfactor=idenMAPfactor*idenCOfactor
    if (idenCOMAPfactor.eq.0d0) return
    call add_to_unique_process_list(process,process_unique,order,idenCOMAPfactor,max_channels,ichans,quarks_to_map)
  end subroutine add_to_process_list

  subroutine map_to_canonical_form(process,part,mapping)
    ! cross the two initial state particle PDGs, order according to
    ! the PDG value, (and reflip the two initial states again)
    implicit none
    integer,dimension(next) :: process,part,mapping
    part(1:next)=process(1:next)
    ! cross the initial state
    part(1)=phys_model%get_antipart(part(1))
    part(2)=phys_model%get_antipart(part(2))
    call sort_with_mapping(next,part,mapping)
    ! cross the initial state
    part(1)=phys_model%get_antipart(part(1))
    part(2)=phys_model%get_antipart(part(2))
  end subroutine map_to_canonical_form
  
  subroutine sort_with_mapping(n,array,mapping)
    !
    ! EXAMPLE:
    !
    ! input:
    ! n=5
    ! array=[4, 1, 8, 2, 3]
    !
    ! output:
    ! array=[1, 2, 3, 4, 8]
    ! mapping=[2, 4, 5, 1, 3]
    !
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(inout) :: array
    integer,dimension(n),intent(out) :: mapping
    integer :: i, j, temp
    ! 'array_map' maps the PDG codes to the 'sort_particle' codes
    ! (see Utilities/process_list.py)
    integer,dimension(-24:25),parameter :: array_map=[83,0,0,0,0,0,0&
         &,0,95,88,93,86,91,84,0,0,0,0,12,11,10,9,8,7,0,1,2,3,4,5,6,0,0&
         &,0,0,85,90,87,92,89,94,0,0,0,0,13,80,81,82,96]
    ! Initialize mapping
    mapping = [(i,i=1,n)]
    ! Sort the array and mapping using a simple bubble sort
    do i=1,n-1
       do j=1,n-i
          if (array_map(array(j)) .gt. array_map(array(j+1))) then
             ! Swap array elements
             temp = array(j)
             array(j) = array(j+1)
             array(j+1) = temp
             ! Swap mapping
             temp = mapping(j)
             mapping(j) = mapping(j+1)
             mapping(j+1) = temp
          endif
       enddo
    enddo
  end subroutine sort_with_mapping

  subroutine get_unique_process(process,process_mapped,mapping,check_sf,process_unique,idenMAPfactor,quarks_to_map)
    implicit none
    integer,dimension(next) :: process,process_mapped,process_unique,mapping
    integer :: iproc,map_from,map_to,i,j
    integer,dimension(4,2) :: quarks_to_map
    real(kind=8) :: idenMAPfactor
    logical :: check_sf
    do iproc=1,pgl_unique%nproc
       if (all(process_mapped(1:next).eq.pgl_unique%processes(1:next,iproc))) exit
    enddo
    if (iproc.gt.pgl_unique%nproc) then
       write (*,*) 'Process not found'
       write (*,*) process(1:next)
       write (*,*) process_mapped(1:next)
       stop 1
    endif
    idenMAPfactor=unique_map_value(iproc)
    if (unique_map(iproc).gt.0) then
       do i=1,next
          if ((i.le.2 .and. mapping(i).gt.2) .or. (i.gt.2 .and. mapping(i).le.2)) then
             process_unique(mapping(i))=phys_model%get_antipart(pgl_unique%processes(i,unique_map(iproc)))
          else
             process_unique(mapping(i))=pgl_unique%processes(i,unique_map(iproc))
          endif
       enddo
    else
       do i=1,next
          if ((i.le.2 .and. mapping(i).gt.2) .or. (i.gt.2 .and. mapping(i).le.2)) then
             process_unique(mapping(i))=phys_model%get_antipart(pgl_unique%processes(i,iproc))
          else
             process_unique(mapping(i))=pgl_unique%processes(i,iproc)
          endif
       enddo
    endif
    ! Check if this is a same-flavour process whose amplitudes are
    ! equal to the sum of the amplitudes of two other processes.
    if (unique_map(iproc).lt.0) then
       call find_to_map(iproc,proc_sf_map(iproc,1:2),quarks_to_map)
    elseif (unique_map(iproc).gt.0) then
       call find_to_map(unique_map(iproc),proc_sf_map(unique_map(iproc),1:2),quarks_to_map)
    endif
    if (any(quarks_to_map(1:4,1:2).ne.0)) then
       check_sf=.true.
    else
       check_sf=.false.
    endif
  end subroutine get_unique_process

  subroutine find_to_map(iproc,sf_map,quarks_to_map)
    implicit none
    integer :: iproc
    integer,dimension(2) :: sf_map
    integer,dimension(4,2),intent(out) :: quarks_to_map
    integer :: i,iq,ia
    quarks_to_map(1:4,1:2)=0
    if (all(sf_map.ne.0)) then
       iq=0
       ia=2
       do i=1,next
          if ((is_quark(pgl_unique%processes(i,iproc)) .and. i.gt.2) .or. &
               (is_antiquark(pgl_unique%processes(i,iproc)) .and. i.le.2)) then
             iq=iq+1
             quarks_to_map(iq,1)=abs(pgl_unique%processes(i,sf_map(1)))
             quarks_to_map(iq,2)=abs(pgl_unique%processes(i,sf_map(2)))
          elseif ((is_antiquark(pgl_unique%processes(i,iproc)) .and. i.gt.2) .or. &
                  (is_quark(pgl_unique%processes(i,iproc)) .and. i.le.2)) then
             ia=ia+1
             quarks_to_map(ia,1)=-abs(pgl_unique%processes(i,sf_map(1)))
             quarks_to_map(ia,2)=-abs(pgl_unique%processes(i,sf_map(2)))
          elseif (pgl_unique%processes(i,iproc).ne.pgl_unique%processes(i,sf_map(1)) .or. &
               pgl_unique%processes(i,iproc).ne.pgl_unique%processes(i,sf_map(2))) then
             write (*,*) 'ERROR in sf_map: inconsistent process definitions',sf_map
             write (*,*) pgl_unique%processes(i,iproc)
             write (*,*) pgl_unique%processes(i,sf_map(1))
             write (*,*) pgl_unique%processes(i,sf_map(2))
             stop 1
          endif
       enddo
    elseif(any(sf_map.ne.0)) then
       write (*,*) 'ERROR in sf_map',sf_map
    endif
  end subroutine find_to_map
  
  subroutine add_to_unique_process_list(process,process_unique,order,idenCOMAPfactor,max_channels,ichans,quarks_to_map)
    implicit none
    integer,intent(in) :: max_channels
    integer,dimension(0:max_channels),intent(in) :: ichans
    integer,dimension(next) :: process,process_unique,order
    integer,dimension(4,2) :: quarks_to_map
    real(kind=8) :: idenCOMAPfactor
    integer :: iproc
    call move_colour_singlet_in_order(process,order)
    do iproc=1,nprocs
       if (all(process_unique(1:next).eq.processes(1:next,iproc)) &
            .and. all(order(1:next).eq.color_orders(1:next,iproc))) exit
    enddo
    if (iproc.gt.nprocs) then
       ! new matrix element to generate
       nprocs=nprocs+1
       processes(1:next,nprocs)=process_unique(1:next)
       color_orders(1:next,nprocs)=order(1:next)
       iden_iproc(nprocs)=1
       iden_processes(1:next,iden_iproc(nprocs),nprocs)=process(1:next)
       idenCOandMAPfactor(iden_iproc(nprocs),nprocs)=idenCOMAPfactor
       multi_chans(0:ichans(0),nprocs)=ichans(0:ichans(0))
       if (any(quarks_to_map.ne.0)) then
          call fill_same_flavour_process_map(nprocs,quarks_to_map)
       else
          same_flavour_process_map(1:2,nprocs)=0
       endif
    else
       ! identical to another matrix element
       iden_iproc(iproc)=iden_iproc(iproc)+1
       iden_processes(1:next,iden_iproc(iproc),iproc)=process(1:next)
       idenCOandMAPfactor(iden_iproc(iproc),iproc)=idenCOMAPfactor
       if (ichans(0).ne.multi_chans(0,iproc)) then
          write (*,*) 'Number of multi-channels not the same among identical processes',&
               ichans(0),multi_chans(0,iproc)
          stop 1
       endif
       if (any(ichans(1:ichans(0)).ne.multi_chans(1:ichans(0),iproc))) then
          write (*,*) 'Multi-channels not the same among identical processes'
          write (*,*) ichans(1:ichans(0))
          write (*,*) multi_chans(1:ichans(0),iproc)
          stop 1
       endif
    endif
  end subroutine add_to_unique_process_list

  subroutine fill_same_flavour_process_map(nprocs,quarks_to_map)
    implicit none
    integer :: nprocs
    integer,dimension(4,2) :: quarks_to_map
    integer,dimension(next) :: proc
    integer :: i,iproc,iq,ia,isum
    do isum=1,2
       iq=0
       ia=2
       do i=1,next
          if ((is_quark(processes(i,nprocs)) .and. i.gt.2) .or. (is_antiquark(processes(i,nprocs)) .and. i.le.2)) then
             iq=iq+1
             if (i.gt.2) then
                proc(i)=quarks_to_map(iq,isum)
             else
                proc(i)=-quarks_to_map(iq,isum)
             endif
          elseif ((is_antiquark(processes(i,nprocs)) .and. i.gt.2) .or. (is_quark(processes(i,nprocs)) .and. i.le.2)) then
             ia=ia+1
             if (i.gt.2) then
                proc(i)=quarks_to_map(ia,isum)
             else
                proc(i)=-quarks_to_map(ia,isum)
             endif
          else
             proc(i)=processes(i,nprocs)
          endif
       enddo
       do iproc=1,nprocs
          if (all(proc(1:next).eq.processes(1:next,iproc)) .and. &
               all(color_orders(1:next,nprocs).eq.color_orders(1:next,iproc))) then
             same_flavour_process_map(isum,nprocs)=iproc
          endif
       enddo
    enddo
    if (same_flavour_process_map(1,nprocs).eq.same_flavour_process_map(2,nprocs)) then
       write (*,*) 'ERROR: found twice the same process',nprocs
       stop 1
    elseif (any(same_flavour_process_map(1:2,nprocs).eq.0)) then
       write (*,*) 'Could not find (one of the) processes',nprocs,same_flavour_process_map(1:2,nprocs)
       stop 1
    endif
  end subroutine fill_same_flavour_process_map
  
  subroutine move_colour_singlet_in_order(process,order)
    ! Move the colour singlet(s) to go *before* the final anti-quark in the order
    implicit none
    integer,dimension(next),intent(in) :: process
    integer,dimension(next),intent(inout) :: order
    integer :: i,iord,aq,iaq,itmp,ipart
    ! find the final anti-quark
    do i=next,1,-1
       iord=order(i)
       ipart=process(iord)
       if (((iord.le.2 .and. is_quark(ipart)).or.(iord.gt.2 .and. is_antiquark(ipart)))) then
          aq=i
          iaq=iord
          exit
       endif
    enddo
    ! move the colour_singlet to after the anti-quark
    i=1
    do while (i.le.next)
       iord=order(i)
       ipart=process(iord)
       if(is_singlet(ipart)) then
          if (i.gt.aq) then
             order(aq:i)=[order(i),order(aq)]
             aq=aq+1
             i=i+1
          else
             order(i:aq)=[order(i+1:aq),order(i)]
             aq=aq-1
          endif
       else
          i=i+1
       endif
    enddo
  end subroutine move_colour_singlet_in_order
          
  
end module read_process_file
