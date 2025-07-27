
program matrix_integrate_QCD
  use common
  use mint_module
  use phase_space_base
  use phase_space_gen23_mod
  use phase_space_genpt_mod
  use phase_space_haag_mod
  use math_functions
  use particles
  use amplitude_QCD_mod
  use cuts
  use pdf_wrap
  use handling_events
  use read_process_file
  use handling_processes
  use multichannel
  implicit none
  integer :: iproc
  real(kind=8) :: weight
  integer :: j,i
  real(kind=8),dimension(:),allocatable :: mass,width
  character(len=80) :: filename
  integer(kind=4) :: PS_choice
  integer,parameter :: nevent_hel_filter=10
  integer :: igroup,integration_step
  logical :: read_amps_from_file=.false.,write_amps_to_file=.false.

  call cpu_time(tTot_B)

  ! relevant input parameters for integration
  ! Number of events to generate. (If negative, start
  ! from a small number of points and double it each
  ! iteration. If positive, this is the number of
  ! points per iteration as well).
!!$  if (integration_step.eq.0 .or. integration_step.eq.2) then
     ncalls0=-500000
!!$  else
!!$     ncalls0=640000
!!$  endif

  itmax=16         ! Number of iterations. (If ncalls0 < 0, the
                   ! integration is aborted if accuracy (next line)
                   ! has been reached.

  ! setting energy
  sqrts=14000.d0

  if (include_pdf) call PDF_initialise

  call phys_model%init_part(173d0,1.491500d0,91.188d0,2.441404d0,80.419002445756163d0,2.0476d0)
!!$  call phys_model%init_part(173d0,0d0,91.188d0,2.441404d0,80.419002445756163d0,2.0476d0)
  call phys_model%init_vert()

  call get_run_arguments()

  if (integration_step.eq.0) then
     accuracy=0.001d0 ! Accuracy of the integration. (Ignored if ncalls0 > 0).
     write_amps_to_file=.true.
  else
     accuracy=max(1d0/sqrt(dble(abs(ncalls0))),0.0005d0)
     read_amps_from_file=.true.
  endif

  
  ! Not so relevant mint-module parameters: only used in special cases.
  call set_mint_module_special_parameters()
  nchans=ngroups ! overwrite the number of mint-channels to the number of
                 ! needed integration channels.

  call create_run_tag(integration_step)
  
  if (read_amps_from_file .or. write_amps_to_file) then
       open(file='Outputs'//trim(adjustl(add_arg))//'/Res_files/amplitudes.bin',&
            unit=32,access='stream',form='unformatted',status='UNKNOWN')
  endif

  do igroup=1,ngroups
     if (pgl(igroup)%nproc.eq.0) cycle
     ! allocate the amplitudes and the phase-space for each of the integration channels
     if (PS_choice.eq.1) then
        allocate(phase_space_gen23 :: pgl(igroup)%phase_space)
     elseif  (PS_choice.eq.2) then
        allocate(phase_space_haag :: pgl(igroup)%phase_space)
     elseif (PS_choice.eq.3) then
        allocate(phase_space_genpt :: pgl(igroup)%phase_space)
     elseif (PS_choice.eq.4) then
        allocate(phase_space_gen23 :: pgl(igroup)%phase_space)
     endif

     allocate(mass(next))
     allocate(width(next))
     do i=1,next
        mass(i)=phys_model%get_mass(pgl(igroup)%processes(i,1))
        width(i)=phys_model%get_width(pgl(igroup)%processes(i,1))
        do iproc=2,pgl(igroup)%nproc
           if (mass(i).ne.phys_model%get_mass(pgl(igroup)%processes(i,iproc)) .or. &
                width(i).ne.phys_model%get_width(pgl(igroup)%processes(i,iproc))) then
              write (*,*) 'masses and widths not compatible among processes'
              stop 1
           endif
        enddo
     enddo
     call setup_spin(pgl(igroup))

     ! Initialise the phase-space parametrisation
     call cpu_time(tBefore)
     call setup_cuts_for_each_particle(pgl(igroup))
     if (PS_choice.ge.1 .and. PS_choice.le.3) then
        call pgl(igroup)%phase_space%init(sqrts,next,mass,pgl(igroup)%phase_space_orders,&
             pgl(igroup)%pt_min,pgl(igroup)%eta_max,pgl(igroup)%DR_min,pgl(igroup)%sqrt_s_min,.false.,include_pdf)
     elseif (PS_choice.eq.4) then
        call pgl(igroup)%phase_space%init(sqrts,next,mass,pgl(igroup)%phase_space_orders,&
             pgl(igroup)%pt_min,pgl(igroup)%eta_max,pgl(igroup)%DR_min,pgl(igroup)%sqrt_s_min,.false.,include_pdf)
     endif
     call cpu_time(tAfter)
     t_PS_init=t_PS_init+tAfter-tBefore
     deallocate(mass)
     deallocate(width)

     allocate(pgl(igroup)%iden(pgl(igroup)%nproc))
     pgl(igroup)%iden(1:pgl(igroup)%nproc)=1
     call set_final_state_identical_particle_factor(pgl(igroup)) ! updates 'iden()'
     call set_initial_state_average_factor(pgl(igroup))          ! updates 'iden()'

     if (include_pdf) then
        call set_ipdgs_for_PDF(pgl(igroup))
     endif

     ! initialize the amplitudes. This creates the whole tree-structure from
     ! which the amps%evaluation() can compute the amplitudes for given
     ! phase-space points.
     call cpu_time(tBefore)
     if (read_amps_from_file) then
        call pgl(igroup)%amps%read_init_amps_from_file(next,32)
     else
        call pgl(igroup)%amps%init(1,next,pgl(igroup)%nproc,pgl(igroup)%processes,&
                pgl(igroup)%spin,pgl(igroup)%color_orders,phys_model,read_proc_from_file)
     endif

     !if (.not.read_proc_from_file .and. pgl(ichan)%amps%same_flav(3)) pgl(igroup)%nproc=3
     
     if (write_amps_to_file) then
        call pgl(igroup)%amps%write_init_amps_to_file(next,32)
     endif

     call cpu_time(tAfter)
     t_amp_init=t_amp_init+tAfter-tBefore

     ! Total number of amplitudes is stored in 'nhel'
     pgl(igroup)%nhel=pgl(igroup)%amps%n_amps

!     if (.not.read_proc_from_file.and.pgl(igroup)%amps%same_flav(3)) then
!        allocate(pgl(igroup)%col_fac(pgl(igroup)%amps%nprocs))
!     else
        allocate(pgl(igroup)%col_fac(pgl(igroup)%nproc))
!     endif

     call compute_LC_colour_factor(pgl(igroup))  ! updates 'col_fac()'

     allocate(pgl(igroup)%amp2(pgl(igroup)%nproc))

     ! number of helicities to sum over
     allocate(pgl(igroup)%amp2_hel(1:pgl(igroup)%nhel))
     allocate(pgl(igroup)%hel(1:next))
     allocate(pgl(igroup)%hel_fac(1:pgl(igroup)%nhel))
     pgl(igroup)%hel_fac(1:pgl(igroup)%nhel)=1

  enddo ! loop over phase-space-order groups
  
  if (read_amps_from_file .or. write_amps_to_file) then
     close(32)
  endif

  if (integration_step.le.1) then
     ! grid setup, or computation of upper bounding envelope
     call mint(integrand)
  else
     ! actual (unweighted) event generation
     call read_grids_from_file
     call gen(integrand,0,-1) ! initialise counters
     filename='Outputs'//trim(adjustl(add_arg))//'/events'//trim(adjustl(tag))//'.lhe'
     open(unit=11,file=filename,status='unknown')
     if (COMMAND_ARGUMENT_COUNT().le.10) &
          call write_unique_in_file_and_deallocate(pgl_unique,unique_map,unique_map_value)
     do j=1,abs(ncalls0)
        call gen(integrand,1,2) ! generate an unweighted event
        call unwgt_process(pgl(ichan))      ! pick a random process
        call unwgt_helicity(pgl(ichan))     ! pick a random helicity for the process picked
        call write_event(11,pgl(ichan),ans(1,0))
     enddo
     close(11)
     call gen(integrand,3,-1) ! print counters
  endif
     
  call cpu_time(tTot_a)
  t_all=tTot_a-tTot_b
  write(*,*) 'Time spent in phase-space initialisation:',t_PS_init 
  write(*,*) 'Time spent in amplitude initialisation',t_Amp_init
  write(*,*) 'Time spent in phase-space generation:',t_PS
  write(*,*) 'Time spent in amplitude evaluation',t_Amp
  write(*,*) 'Time spent in squaring amplitudes',t_mat
  write(*,*) 'Total time:',t_all
  write(*,*) 'Number of events:',pgl(1:ngroups)%all_evt
  write(*,*) 'Number passing cuts:',pgl(1:ngroups)%passed
  write(*,*) 'Fraction passing:',float(pgl(1:ngroups)%passed)/float(pgl(1:ngroups)%all_evt)
 
contains

  real(kind=8) function integrand(x,vol,ifirst,f1)
    implicit none
    integer :: ifirst
    real(kind=8), dimension(ndim) :: x
    real(kind=8), dimension(nintegrals) :: f1
    real(kind=8), dimension(:),allocatable,save :: val,val_abs,vol_ichan
    real(kind=8),dimension(pgl(ichan)%nproc) :: colour_singlet_multichannel_weight
    integer :: ih,iproc,i
    real(kind=8) :: vol
    real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real(kind=4) :: tBefore,tAfter
    if (.not.allocated(val)) then
       allocate(val(1:maxval(pgl(1:ngroups)%nproc)))
       allocate(val_abs(1:maxval(pgl(1:ngroups)%nproc)))
       allocate(vol_ichan(1:ngroups))
    endif
    ! some point-by-point initialisation
    f1(1:nintegrals)=0d0
    if (pgl(ichan)%nproc.eq.0) return
    if (ifirst.eq.2) then
       ! use previously computed integrand
       f1(1)=sum(val_abs(1:pgl(ichan)%nproc))
       f1(2)=sum(val(1:pgl(ichan)%nproc))
       f1(3:pgl(ichan)%nproc+2)=val(1:pgl(ichan)%nproc)
       return
    endif
    new_point=.true.
    pass_cuts_check=.true.
    val_abs(1:pgl(ichan)%nproc)=0d0

    ! Generate phase-space point based on the random numbers 'x(1:ndim)'
    call cpu_time(tBefore)
    call pgl(ichan)%phase_space%generate_momenta(x)
    call cpu_time(tAfter)
    t_PS= t_PS +tAfter-tBefore
    
    if (debug ) then
       write (*,*) pgl(ichan)%phase_space%jac
       stop 1
    endif
    
    pgl(ichan)%all_evt=pgl(ichan)%all_evt+1
    
    if (pgl(ichan)%phase_space%jac.lt.0d0) then
       pass_cuts_check=.false.
       val(1:pgl(ichan)%nproc)=0d0
       return
    endif

    if (.not.pass_cuts(pgl(ichan))) then
       pass_cuts_check=.false.
       val(1:pgl(ichan)%nproc)=0d0
       return
    endif

    pgl(ichan)%passed = pgl(ichan)%passed + 1

    ! compute amplitudes
    call cpu_time(tBefore)

    call pgl(ichan)%amps%evaluate(next,pgl(ichan)%phase_space%p,pgl(ichan)%hel,read_proc_from_file,phys_model)
    call cpu_time(tAfter)
    t_amp=t_amp+tAfter-tBefore
    
    call compute_multichannel_weight(ichan,pgl(ichan)%phase_space%x,pgl(ichan)%phase_space%p, &
                                     pgl(ichan)%phase_space%jac,colour_singlet_multichannel_weight)
    call cpu_time(tBefore)
    iproc=0
    pgl(ichan)%amp2(:)=0d0
    if (use_real_gluons .and. all(pgl(ichan)%amps%n_qqbar(1:pgl(ichan)%amps%nprocs).eq.0)) then
       do ih=1,pgl(ichan)%amps%n_amps
          do while (pgl(ichan)%amps%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
          pgl(ichan)%amp2_hel(ih)=pgl(ichan)%amps%amps_r(ih)*pgl(ichan)%col_fac(iproc)*pgl(ichan)%amps%amps_r(ih)*&
                  pgl(ichan)%hel_fac(ih)
          pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
       enddo
    else
       do ih=1,pgl(ichan)%amps%n_amps
          do while (pgl(ichan)%amps%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
          pgl(ichan)%amp2_hel(ih)=dble(pgl(ichan)%amps%amps(ih)*pgl(ichan)%col_fac(iproc)*dconjg(pgl(ichan)%amps%amps(ih)))*&
               pgl(ichan)%hel_fac(ih)
!          if (.not.read_proc_from_file .and. pgl(ichan)%amps%nprocs.eq.3 .and. &
!                  .not. pgl(ichan)%amps%same_flav(iproc)) then
!               cycle
!          elseif (.not.read_proc_from_file .and. pgl(ichan)%amps%nprocs.eq.3 .and. &
!                  pgl(ichan)%amps%same_flav(iproc)) then 
!               iproc=1
!          endif
          pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
       enddo
    endif
    if (pgl(ichan)%passed.le.nevent_hel_filter) then
       call setup_helicity_filter(pgl(ichan))
       if (integration_step.eq.2 .and. pgl(ichan)%passed.eq.nevent_hel_filter) then
          ! since we update the helicities we need to compute when
          ! passed==nevent_hel_filter, the unweighting of the helicities goes
          ! wrong for this phase-space point. Hence, we need to skip it.
          pgl(ichan)%amp2(1:pgl(ichan)%nproc)=0d0
       endif
    endif

    ! MINT weight, phase-space jacobian and GeV -> pb conversion factor
    weight=vol*pgl(ichan)%phase_space%jac*conv

    ! multiply by the strong coupling
    if (pgl(ichan)%amps%n_sing(1).lt.next-2) then
       weight=weight*(4*pi*alphas)**(next-2-pgl(ichan)%amps%n_sing(1))
    endif
    
    ! multiply by the EW coupling
    if (pgl(ichan)%amps%n_sing(1).ge.1) then
       weight=weight*(2d0*4d0*pi*alphaEW)**pgl(ichan)%amps%n_sing(1)
    endif

    val(1:pgl(ichan)%nproc)=pgl(ichan)%amp2(1:pgl(ichan)%nproc)*weight/dble(pgl(ichan)%iden(1:pgl(ichan)%nproc))
    val(1:pgl(ichan)%nproc)=val(1:pgl(ichan)%nproc)*colour_singlet_multichannel_weight(1:pgl(ichan)%nproc)

    call include_PDF_and_identical_procs(val,val_abs,pgl(ichan))

    ! pass the result to the mint module
    f1(1)=sum(val_abs(1:pgl(ichan)%nproc))
    f1(2)=sum(val(1:pgl(ichan)%nproc))
    f1(3:pgl(ichan)%nproc+2)=val(1:pgl(ichan)%nproc)
    
    integrand=f1(1)
    call cpu_time(tAfter)
    t_mat=t_mat+tAfter-tBefore
  end function integrand

  subroutine setup_helicity_filter(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    real(kind=8) :: max_value
    integer :: ih1,ih2,iproc1,iproc2
    if (.not.allocated(pgl%include_hel)) then
       allocate(pgl%include_hel(pgl%nhel))
       pgl%include_hel(1:pgl%nhel)=0
    endif
    ! filter zero helicities and helicities that are identical
    max_value=maxval(pgl%amp2_hel(1:pgl%nhel))
    do ih1=1,pgl%nhel
       if (pgl%include_hel(ih1).ne.0) cycle
       if (pgl%amp2_hel(ih1)/max_value.gt.1d-10) then
          ! non-zero
          pgl%include_hel(ih1)=1
       else
          cycle
       endif
       do ih2=ih1+1,pgl%nhel
          if (abs(pgl%amp2_hel(ih1)-pgl%amp2_hel(ih2))/abs(pgl%amp2_hel(ih1)+pgl%amp2_hel(ih2)).lt.1d-10) then
             ! identical value. Now check that they belong to the same process
             iproc1=1; do while (iproc1.lt.pgl%nproc .and. (pgl%amps%iproc_start(iproc1+1)-ih1).le.0) ; iproc1=iproc1+1 ; enddo
             iproc2=1; do while (iproc2.lt.pgl%nproc .and. (pgl%amps%iproc_start(iproc2+1)-ih2).le.0) ; iproc2=iproc2+1 ; enddo
             if (iproc1.ne.iproc2) cycle
             ! identical process
             pgl%include_hel(ih2)=-ih1
             pgl%include_hel(ih1)=pgl%include_hel(ih1)+1
          endif
       enddo
    enddo

    if (pgl%passed.lt.nevent_hel_filter) return

    ih2=0
    do ih1=1,pgl%nhel
       if (pgl%include_hel(ih1).gt.0) ih2=ih2+1
    enddo

    call pgl%amps%filter_helicity(next,pgl%nhel,pgl%include_hel) ! this updates 'nhel' and 'include_hel'
    deallocate(pgl%hel_fac)
    allocate(pgl%hel_fac(pgl%nhel))
    pgl%hel_fac(1:pgl%nhel)=pgl%include_hel(1:pgl%nhel)
    deallocate(pgl%include_hel)
  end subroutine setup_helicity_filter


  subroutine get_run_arguments()
    implicit none
    integer :: argc,n_ps
    integer :: i,k,iproc
    character(len=256) :: argv
    logical :: found_1
    integer(kind=8) :: sym_fac
    integer(kind=8) iseed
    integer,dimension(:,:),allocatable :: ps_o
    logical :: same_flavour
    common /to_seed/iseed
    iseed=0
    ! integration steps:
    ! integration_step=0  (Setting up grids)
    ! integration_step=-1 (same as integration_step=0, but starting from existing grids)
    ! integration_step=1  (computing bounding envelope)
    ! integration_step=2  (event generation)
    argc = COMMAND_ARGUMENT_COUNT()
    if (argc.ge.3 .and. argc .le. 4) then
       read_proc_from_file=.true.
       do i=1,argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) read(argv,'(a)') filename
          if (i.eq.2) read(argv,*) PS_choice
          if (i.eq.3) read(argv,*) integration_step
          if (i.eq.4) then
             read(argv,*) iseed
             write(add_arg(1:),*) iseed
          endif
       enddo
       call read_processes_from_file(filename)
       if (integration_step.le.1) then
          ! make sure pgl_unique is deallocated consistently:
          call finalize_phase_space_order_group(pgl_unique)
          deallocate(pgl_unique)
       endif
999    continue
       close(10)
    elseif (argc.le.10) then
       write(*,*) 'Inconsistent arguments:'
       write(*,*) '--------- Should be: --------'
       write(*,*) 'PS_choice, mode, next, *process*, *order*'
       stop 2
    else
       read_proc_from_file=.false.
       do i = 1, argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) read(argv,*) PS_choice
          if (i.eq.2) read(argv,*) integration_step
          if (i.eq.3) then
             read(argv,*) next
             if (next.le.3) then
                write (*,*) 'Need at least 4 particles (2->2 scattering)',next
                stop 1
             endif
             ndim=3*(next-2)-4
             if (include_pdf) ndim=ndim+2
             allocate(part(1:next))
             allocate(o(1:next))
          endif
          do k=0,next-1
             if (i.eq.4+k) then
                read(argv,*) part(k+1)
             endif
          enddo
          do k=0,next-1
             if (i.eq.4+next+k) then
                read(argv,*) o(k+1)
             endif
          enddo
          if (argc.eq.3+2*next +1 .and. i.eq.argc) then
             ! Special case: we have an additional argument. Use it as a special tag
             read(argv,*) add_arg
             read (add_arg((index(add_arg,'S'))+1:(index(add_arg,'I'))-1),*,err=99) iseed
99           continue             
          endif
       enddo
       write (*,*) '******************************************'
       write (*,*) 'Process is     ',part
       write (*,*) 'Colour order is',o
       write (*,*) '******************************************'
       nquarks = 0
       do i=1,next
          if ((abs(part(i)).ge.1).and.(abs(part(i)).le.6)) then
             nquarks = nquarks + 1
          endif
       enddo
       if (nquarks.eq.0) then
          ! count the number of gluons between '1' and '2' in the colour order
          c_o=0
          i=0
          found_1=.false.
          do
             i=i+1
             if (i.gt.next) i=1
             if (found_1) c_o=c_o+1
             if (o(i).eq.1) then
                found_1=.true.
             endif
             if (o(i).eq.2 .and. found_1) then
                c_o=c_o-1
                exit
             endif
          enddo
       endif

       do i=2,next-1
          if ((o(i).gt.2 .and. part(o(i)).le.-1 .and. part(o(i)).ge.-6) .or. &
              (o(i).le.2 .and. part(o(i)).ge. 1 .and. part(o(i)).le. 6) ) then
             if (abs(part(o(i))).eq.abs(part(o(1))) .and. abs(part(o(i))).eq.abs(part(o(next)))) then
                same_flavour=.true.
             else
                same_flavour=.false.
             endif
             exit
          endif
       enddo

       if ((nquarks.ne.0 .and. nquarks.ne.2 .and. nquarks.ne.4) .or. (nquarks.gt.next)) then
          write (*,*) 'Not consistent number of external quarks (up to 2)',nquarks
          stop
       endif
       call compute_multichannel_symmetry_factor(sym_fac)

       call determine_multi_channel_size(part,n_ps)
       
       allocate(ps_o(1:next,1:n_ps))
       call determine_phase_space_orders(part,o,n_ps,ps_o)
       
       ngroups=n_ps
       allocate(pgl(ngroups))
       do i=1,ngroups
          pgl(i)%nproc=1
          allocate(pgl(i)%multichan%channels(1:n_ps,1))
          allocate(pgl(i)%multichan%number_of_channels(1))
          pgl(i)%multichan%number_of_channels(1)=n_ps
          do k=1,n_ps
             pgl(i)%multichan%channels(k,1)=k
          enddo
          allocate(pgl(i)%processes(1:next,pgl(i)%nproc))
          pgl(i)%processes(1:next,1)=part(1:next)
          allocate(pgl(i)%color_orders(1:next,pgl(i)%nproc))
          pgl(i)%color_orders(1:next,1)=o(1:next)
          allocate(pgl(i)%phase_space_orders(1:next))
          pgl(i)%phase_space_orders(1:next)=ps_o(1:next,i)
          allocate(pgl(i)%idenCOandMAPfactor(1,1))
          pgl(i)%idenCOandMAPfactor(1,1)=sym_fac
          allocate(pgl(i)%iden_iproc(1))
          pgl(i)%iden_iproc(1)=1
          allocate(pgl(i)%val_procs(1,1))
          allocate(pgl(i)%iden_processes(1:next,1,1))
          pgl(i)%iden_processes(1:next,1,1)=pgl(i)%processes(1:next,1)
       enddo
    endif

    do i=1,ngroups
       call setup_optimised_multichannel_weight_computation(pgl(i))
    enddo
    
    
    ! basic checks:
    if (next.lt.4) then
       write (*,*) 'Not enough external particles',next
       stop 1
    endif
    if (integration_step.ne.0 .and. integration_step.ne.1 .and. integration_step.ne.2) then
       write (*,*) 'Incorrect integration_step',integration_step
       stop
    else
       imode=integration_step ! imode is used in MINT
    endif
    if (PS_choice.ne.1 .and. PS_choice.ne.2 .and. PS_choice.ne.3 .and. PS_choice.ne.4) then
       write (*,*) 'PS_Choice modes only 1, 2, 3 or 4',PS_choice
       stop
    endif
  end subroutine get_run_arguments

  
end program matrix_integrate_QCD
