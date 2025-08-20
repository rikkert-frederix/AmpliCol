module read_process_file
  use handling_processes
  integer :: sf_nprocs
  integer,dimension(:,:),allocatable :: unique_procs,processes,color_orders,multi_chans
  type(phase_space_order_group),allocatable :: pgl_unique
  real(kind=8),dimension(:),allocatable :: unique_map_value
  integer,dimension(:),allocatable :: unique_map,iden_iproc,same_flavour_idenCOfactor
  integer,dimension(:,:,:),allocatable :: iden_processes
  real(kind=8),dimension(:,:),allocatable :: idenCOandMAPfactor
contains
  subroutine read_processes_from_file(filename)
    implicit none
    character(len=80) :: filename
    integer :: iproc,igroup,icheck,nproc_in_group,max_channels,idenCOfactor,iflav,ndim
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
    pgl_unique%ndim=ndim
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
       do iflav=1,2
          do iproc=1,nproc_in_group
             read(10,'(a)') buff
             read(buff,*) ichans(0)
             read(buff,*) ichans(0),ichans(1:ichans(0)),process(1:next),order(1:next),idenCOfactor
             if (iflav.eq.1) then
                call add_to_process_list(process,order,idenCOfactor,max_channels,ichans,.true.)
             else
                call add_to_process_list(process,order,idenCOfactor,max_channels,ichans,.false.)
             endif
          enddo
          if (iflav.eq.1) then
             if (.not.decompose_same_flavour_into_two_diff_flavour) exit
             do iproc=1,nproc_in_group
                backspace(10)
             enddo
          endif
       enddo
       pgl(igroup)%next=next
       pgl(igroup)%nproc=nprocs
       pgl(igroup)%ndim=ndim
       pgl(igroup)%multichan%max_channels=max_channels
       if (keep_processes_separate) then
          allocate(pgl(igroup)%amps(1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%nhel(1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%passed(1:pgl(igroup)%nproc))
       else
          allocate(pgl(igroup)%amps(1))
          allocate(pgl(igroup)%nhel(1))
          allocate(pgl(igroup)%passed(1))
       endif
       allocate(pgl(igroup)%processes(1:next,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%color_orders(1:next,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%phase_space_orders(1:next))
       allocate(pgl(igroup)%idenCOandMAPfactor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%iden_iproc(1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%val_procs(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%multichan%channels(1:max_channels,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%multichan%number_of_channels(1:pgl(igroup)%nproc))
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
       pgl(igroup)%passed=0
       deallocate(iden_iproc)
       deallocate(processes)
       deallocate(color_orders)
       deallocate(phase_space_orders)
       deallocate(iden_processes)
       deallocate(idenCOandMAPfactor)
       deallocate(multi_chans)
       deallocate(ichans)
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
    integer :: i,iproc,ih
    integer,parameter :: nevent=10
    real(kind=8),dimension(:,:),allocatable :: amp2
    complex(kind=8),dimension(:,:),allocatable :: amp
    real(kind=8),dimension(:),allocatable :: mass,width
    real(kind=8),dimension(pgl_unique%ndim) :: x
    real(kind=8),external :: ran2
    allocate(phase_space_gen23 :: pgl_unique%phase_space)
    allocate(pgl_unique%processes(next,nproc_unique))
    allocate(pgl_unique%color_orders(next,nproc_unique))
    allocate(pgl_unique%phase_space_orders(next))
    allocate(pgl_unique%amps(1))
    allocate(mass(next))
    allocate(width(next))
    pgl_unique%nproc=nproc_unique
    pgl_unique%processes(1:next,1:nproc_unique)=unique_procs(1:next,1:nproc_unique)
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

    call pgl_unique%amps(1)%init(1,next,pgl_unique%nproc,pgl_unique%processes,&
         pgl_unique%spin,pgl_unique%color_orders,phys_model)

    allocate(amp2(nevent,pgl_unique%nproc))
    allocate(amp(nevent,pgl_unique%amps(1)%n_amps))
    allocate(pgl_unique%passed(1))

    pgl_unique%passed=0
    do while (pgl_unique%passed(1).lt.nevent)
       do i=1,pgl_unique%ndim
          x(i)=ran2()
       enddo
       call pgl_unique%phase_space%generate_momenta(x)
       if (pgl_unique%phase_space%jac.lt.0d0) cycle
       pgl_unique%passed(1)=pgl_unique%passed(1)+1
       call pgl_unique%amps(1)%evaluate(next,pgl_unique%phase_space%p,pgl_unique%hel,read_proc_from_file,phys_model,.false.)
       iproc=0
       amp2(pgl_unique%passed(1),:)=0d0
       if (use_real_gluons .and. all(pgl_unique%amps(1)%n_qqbar(1:pgl_unique%amps(1)%nprocs).eq.0)) then
          do ih=1,pgl_unique%amps(1)%n_amps
             do while (pgl_unique%amps(1)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             amp2(pgl_unique%passed(1),iproc)=amp2(pgl_unique%passed(1),iproc)+&
                  pgl_unique%amps(1)%amps_r(ih)*pgl_unique%amps(1)%amps_r(ih)
          enddo
       else
          do ih=1,pgl_unique%amps(1)%n_amps
             do while (pgl_unique%amps(1)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             amp2(pgl_unique%passed(1),iproc)=amp2(pgl_unique%passed(1),iproc)+&
                  dble(pgl_unique%amps(1)%amps(ih)*dconjg(pgl_unique%amps(1)%amps(ih)))
          enddo
          amp(pgl_unique%passed(1),1:pgl_unique%amps(1)%n_amps)=pgl_unique%amps(1)%amps(1:pgl_unique%amps(1)%n_amps)
       endif
       call find_same_flavour(pgl_unique,nevent)
    enddo
    allocate(unique_map(1:pgl_unique%nproc))
    allocate(unique_map_value(1:pgl_unique%nproc))
    call find_unique(pgl_unique,nevent,amp2,unique_map,unique_map_value)

    do iproc=1,pgl_unique%nproc
       write (*,*) unique_map(iproc),unique_map_value(iproc),':',pgl_unique%processes(1:pgl_unique%next,iproc),&
            ':',pgl_unique%color_orders(1:pgl_unique%next,iproc)
    enddo
    deallocate(pgl_unique%spin)
    deallocate(pgl_unique%phase_space)
    deallocate(amp2)
    deallocate(amp)
  end subroutine check_unique_processes

  subroutine find_unique(pgl,nevent,amp2,unique_map,unique_map_value)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: nevent
    real(kind=8),dimension(nevent,pgl%nproc),intent(in) :: amp2
    real(kind=8),dimension(pgl%nproc),intent(out) :: unique_map_value
    integer,dimension(pgl%nproc),intent(out) :: unique_map
    integer :: i,j
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

  subroutine add_to_process_list(process,order,idenCOfactor,max_channels,ichans,skip_same_flavour)
    implicit none
    integer :: max_channels
    integer,dimension(0:max_channels) :: ichans
    integer,dimension(next) :: process,order,process_unique
    integer :: idenCOfactor
    real(kind=8) :: idenCOMAPfactor,idenMAPfactor
    integer,dimension(0:4) :: quarks
    logical :: skip_same_flavour,is_same_flavour
    call find_quarks(process,order,quarks)
    call get_unique_process_from_quarks(quarks,process,order,process_unique,idenMAPfactor,is_same_flavour)
    if (decompose_same_flavour_into_two_diff_flavour) then
       if (skip_same_flavour .and. is_same_flavour) return
       if ((.not.skip_same_flavour) .and. (.not.is_same_flavour)) return
    endif
    if (reduce_to_unique_matrix_elements) then
       idenCOMAPfactor=idenMAPfactor*idenCOfactor
    else
       idenCOMAPfactor=idenCOfactor
    endif
    if (idenCOMAPfactor.eq.0d0) return
    call add_to_unique_process_list(process,process_unique,order,idenCOMAPfactor,max_channels,ichans)
  end subroutine add_to_process_list


  subroutine get_unique_process_from_quarks(quarks,process,order,process_unique,idenMAPfactor,is_same_flavour)
    implicit none
    integer,dimension(0:4) :: quarks
    integer,dimension(next) :: process,order,process_unique,process_mapped,mapping
    integer :: i,iproc,iq
    real(kind=8) :: idenMAPfactor
    logical :: is_same_flavour
    call map_to_canonical_form(process,process_mapped,mapping)
    do iproc=1,pgl_unique%nproc
       if (pgl_unique%amps(1)%n_qqbar(iproc)*2.ne.quarks(0)) cycle
       if (quarks(0).eq.4) then
          if (all(pgl_unique%processes(1:4,iproc).eq.-abs(quarks(1:4)))) then ! quarks are consistent 
             if (all(process_mapped(5:next).eq.pgl_unique%processes(5:next,iproc))) exit ! and the rest as well
          endif
       else
          if (all(process_mapped(1:next).eq.pgl_unique%processes(1:next,iproc))) exit
       endif
    enddo
    if (iproc.gt.pgl_unique%nproc) then
       write (*,*) 'Process not found',quarks
       stop 1
    endif
    process_unique(1:next)=process(1:next)
    idenMAPfactor=unique_map_value(iproc)
    if (unique_map(iproc).gt.0) then
       iq=0
       do i=1,next
          if (is_quark(process(order(i))) .or. is_antiquark(process(order(i)))) then
             iq=iq+1
             if (quarks(0).eq.4) then
                if (iq.eq.1 .or. iq.eq.4) then
                   process_unique(order(i))=sign(pgl_unique%processes(iq,unique_map(iproc)),quarks(iq))
                elseif (iq.eq.2) then
                   process_unique(order(i))=sign(pgl_unique%processes(3,unique_map(iproc)),quarks(3))
                elseif (iq.eq.3) then
                   process_unique(order(i))=sign(pgl_unique%processes(2,unique_map(iproc)),quarks(2))
                endif
             elseif (quarks(0).eq.2) then
                process_unique(order(i))=sign(pgl_unique%processes(iq,unique_map(iproc)),quarks(iq))
             endif
          endif
       enddo
    endif
    if (pgl_unique%amps(1)%same_flav(iproc)) then
       is_same_flavour=.true.
    else
       is_same_flavour=.false.
    endif
  end subroutine get_unique_process_from_quarks

  subroutine add_to_unique_process_list(process,process_unique,order,idenCOMAPfactor,max_channels,ichans)
    implicit none
    integer,intent(in) :: max_channels
    integer,dimension(0:max_channels),intent(in) :: ichans
    integer,dimension(next) :: process,process_unique,order
    real(kind=8) :: idenCOMAPfactor
    integer :: iproc
    call move_colour_singlet_in_order(process,order)
    if (.not.reduce_to_unique_matrix_elements) then
       ! always add the (unaltered) process
       nprocs=nprocs+1
       processes(1:next,nprocs)=process(1:next)
       color_orders(1:next,nprocs)=order(1:next)
       iden_iproc(nprocs)=1
       iden_processes(1:next,iden_iproc(nprocs),nprocs)=process(1:next)
       idenCOandMAPfactor(iden_iproc(nprocs),nprocs)=idenCOMAPfactor
       multi_chans(0:ichans(0),nprocs)=ichans(0:ichans(0))
       return
    endif
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


  subroutine find_quarks(process,order,quarks)
    implicit none
    integer,dimension(next) :: process,order
    integer,dimension(0:4) :: quarks
    integer :: i
    quarks(0)=0
    do i=1,next
       if (is_quark(process(order(i))) .or. is_antiquark(process(order(i)))) then
          quarks(0)=quarks(0)+1
          quarks(quarks(0))=process(order(i))
       endif
    enddo
    if (quarks(0).eq.4) then
       quarks(1:4)=[quarks(1),quarks(3),quarks(2),quarks(4)]
    endif
  end subroutine find_quarks

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


  subroutine move_colour_singlet_in_order(process,order)
    ! Move the colour singlet(s) to go *before* the final anti-quark in the order
    implicit none
    integer,dimension(next),intent(in) :: process
    integer,dimension(next),intent(inout) :: order
    integer :: i,iord,aq,iaq,ipart
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
