
program matrix_integrate_QCD
  use common
!  use mint_module
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
  use amplitude_library
  implicit none
  integer :: ivec
  integer,parameter :: vector_size=128
  integer :: iproc
  integer :: i
  real(kind=8),dimension(:),allocatable :: mass,width
  character(len=80) :: filename,logfile,tag
  integer(kind=4) :: PS_choice
  integer,parameter :: nevent_hel_filter=10
  integer :: igroup
  logical,dimension(vector_size) :: to_write
  integer,dimension(:),allocatable :: nintegrals
  integer :: ichan,iint,itmax,ncalls0,iamp
  real(kind=8),dimension(vector_size) :: f,f_abs
  logical :: done
  real(kind=8),dimension(:,:),allocatable :: wgts
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  character(len=19) :: formatted
  logical :: create_amplitude_library,use_amplitude_library
  call cpu_time(tTot_B)

  call get_run_arguments()
  
  ! setting energy
  sqrts=14000.d0

  if (include_pdf) call PDF_initialise
  
  call phys_model%init_part(173d0,1.491500d0,91.188d0,2.441404d0,&
                           80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call phys_model%init_vert()

  call cpu_time(tBefore)
  if (use_amplitude_library) then
     call read_amplitude_lib(vector_size)
  else
     call read_processes_from_file(filename,vector_size)
     do i=1,ngroups
       call setup_optimised_multichannel_weight_computation(pgl(i))
    enddo
 endif
  call cpu_time(tAfter)
  t_Proc_init=t_Proc_init+tAfter-tBefore
  call date_and_time(date, time, zone)
  write(formatted, '(A4,"-",A2,"-",A2," ",A2,":",A2,":",A2)') &
       date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
  write (*,*)  'Initialise phase-space groups and amplitudes '//trim(formatted)
  write (99,*) 'Initialise phase-space groups and amplitudes '//trim(formatted)
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

     allocate(mass(pgl(igroup)%next))
     allocate(width(pgl(igroup)%next))
     do i=1,pgl(igroup)%next
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
     ! Initialise the phase-space parametrisation
     call cpu_time(tBefore)
     call setup_cuts_for_each_particle(pgl(igroup),igroup)
     if (PS_choice.ge.1 .and. PS_choice.le.3) then
        call pgl(igroup)%phase_space%init(sqrts,pgl(igroup)%next,mass,pgl(igroup)%phase_space_orders,&
             pgl(igroup)%pt_min,pgl(igroup)%eta_max,pgl(igroup)%DR_min,pgl(igroup)%sqrt_s_min,.false.,include_pdf)
     elseif (PS_choice.eq.4) then
        call pgl(igroup)%phase_space%init(sqrts,pgl(igroup)%next,mass,pgl(igroup)%phase_space_orders,&
             pgl(igroup)%pt_min,pgl(igroup)%eta_max,pgl(igroup)%DR_min,pgl(igroup)%sqrt_s_min,.false.,include_pdf)
     endif
     pgl(igroup)%ndim_extra=pgl(igroup)%phase_space%ndim_extra
     allocate(pgl(igroup)%ps(vector_size))
     do ivec=1,vector_size
        allocate(pgl(igroup)%ps(ivec)%x(1:pgl(igroup)%ndim+pgl(igroup)%ndim_extra))
        allocate(pgl(igroup)%ps(ivec)%p(0:3,1:pgl(igroup)%next))
     enddo
     call cpu_time(tAfter)
     t_PS_init=t_PS_init+tAfter-tBefore
     deallocate(mass)
     deallocate(width)

     if (use_amplitude_library) cycle
     
     call setup_spin(pgl(igroup))

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
     if (keep_processes_separate) then
        do iamp=1,pgl(igroup)%nproc
           call pgl(igroup)%amps(iamp)%init(1,pgl(igroup)%next,1,pgl(igroup)%processes(1,iamp),&
                pgl(igroup)%spin,pgl(igroup)%color_orders(1,iamp),phys_model)
        enddo
     else
        call pgl(igroup)%amps(1)%init(1,pgl(igroup)%next,pgl(igroup)%nproc,pgl(igroup)%processes,&
             pgl(igroup)%spin,pgl(igroup)%color_orders,phys_model)
     endif

     call cpu_time(tAfter)
     t_amp_init=t_amp_init+tAfter-tBefore

     ! Total number of amplitudes is stored in 'nhel'
     if (keep_processes_separate) then
        do iamp=1,pgl(igroup)%nproc
           pgl(igroup)%nhel(iamp)=pgl(igroup)%amps(iamp)%n_amps
        enddo
     else
        pgl(igroup)%nhel=pgl(igroup)%amps(1)%n_amps
     endif
     allocate(pgl(igroup)%col_fac(pgl(igroup)%nproc))

     call compute_LC_colour_factor(pgl(igroup))  ! updates 'col_fac()'

     if (keep_processes_separate) then
        allocate(pgl(igroup)%amp2(1,vector_size))
     else
        allocate(pgl(igroup)%amp2(1:pgl(igroup)%nproc,vector_size))
     endif

     ! number of helicities to sum over
     allocate(pgl(igroup)%amp2_hel(1:maxval(pgl(igroup)%nhel),vector_size))
     allocate(pgl(igroup)%hel(1:pgl(igroup)%next))
     if (keep_processes_separate) then
        allocate(pgl(igroup)%hel_fac(1:maxval(pgl(igroup)%nhel),pgl(igroup)%nproc))
     else
        allocate(pgl(igroup)%hel_fac(1:maxval(pgl(igroup)%nhel),1))
     endif
     pgl(igroup)%hel_fac=1

  enddo ! loop over phase-space-order groups

  if (use_amplitude_library) call test_lib
  
  filename='Outputs/'//trim(adjustl(tag))//'events_tmp.lhe'
  open(unit=11,file=filename,action='readwrite',status='unknown')
  call write_unique_in_file(pgl_unique,unique_map,unique_map_value)
  
  allocate(nintegrals(ngroups))
  if (keep_processes_separate) then
     nintegrals(1:ngroups)=pgl(1:ngroups)%nproc
  else
     nintegrals(1:ngroups)=1
  endif
  call simple_integrator%init(ngroups,pgl(1:ngroups)%ndim,pgl(1:ngroups)%ndim_extra,nintegrals,abs(ncalls0),abs(itmax))
  call date_and_time(date, time, zone)
  write(formatted, '(A4,"-",A2,"-",A2," ",A2,":",A2,":",A2)') &
       date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
  write (*,*) 'Start phase-space integration '//trim(formatted)
  write (99,*) 'Start phase-space integration '//trim(formatted)
  call flush(99)
  do
     call simple_integrator%get_points(vector_size,ichan,iint)
     call integrand(ichan,iint,simple_integrator%x,simple_integrator%wgt,f,f_abs)
     call simple_integrator%fill_points(vector_size,f_abs,f,to_write,done)
     if (create_amplitude_library) then
        done=.true.
        do igroup=1,ngroups
           done=done.and.all(pgl(igroup)%amps%lib_created)
        enddo
        if (done) call create_amplitude_lib()
     endif
     do ivec=1,vector_size
        if (to_write(ivec)) then
           call unwgt_process(pgl(ichan),iint,ivec) ! pick a random process
           call unwgt_helicity(pgl(ichan),ivec)     ! pick a random helicity for the process picked
           call write_event(11,pgl(ichan),ivec,1d0)
        endif
     enddo
     if (done) exit
  enddo
  call flush(11)
  call simple_integrator%assign_evnt_wgts(wgts)
  rewind(11)
  filename='Outputs/'//trim(adjustl(tag))//'events.lhe'
  open(unit=12,file=filename,action='write',status='unknown')
  write (*,*) 'Updating event weights...'
  write (99,*) 'Updating event weights...'
  do i=1,size(wgts,dim=2)
     call event_update_wgt(11,12,wgts(1,i))
  enddo
  close(11,status='DELETE')
  close(12)
     
  call cpu_time(tTot_a)
  t_all=tTot_a-tTot_b
  write(*,*) 'Time spent in phase-space initialisation:',t_PS_init 
  write(*,*) 'Time spent in amplitude initialisation',t_Amp_init
  write(*,*) 'Time spent in process initialisation',t_Proc_init
  write(*,*) 'Time spent in phase-space generation:',t_PS
  write(*,*) 'Time spent in amplitude evaluation',t_Amp
  write(*,*) 'Time spent in squaring amplitudes',t_mat
  write(*,*) 'Total time:',t_all
  write(99,*) 'Time spent in phase-space initialisation:',t_PS_init 
  write(99,*) 'Time spent in amplitude initialisation',t_Amp_init
  write(99,*) 'Time spent in process initialisation',t_Proc_init
  write(99,*) 'Time spent in phase-space generation:',t_PS
  write(99,*) 'Time spent in amplitude evaluation',t_Amp
  write(99,*) 'Time spent in squaring amplitudes',t_mat
  write(99,*) 'Total time:',t_all
  close(99)
  
contains

  subroutine integrand(ichan,iint,x,vol,f,f_abs)
    use amp_lib
    implicit none
    integer,intent(in) :: ichan,iint
    real(kind=8),dimension(pgl(ichan)%ndim+pgl(ichan)%ndim_extra,vector_size),intent(in) :: x
    real(kind=8),dimension(vector_size),intent(in) :: vol
    real(kind=8),dimension(vector_size),intent(out) :: f,f_abs
    real(kind=8),dimension(vector_size) :: weight
    real(kind=8),dimension(:,:),allocatable,save :: val,val_abs,vol_ichan
    real(kind=8),dimension(pgl(ichan)%nproc,vector_size) :: colour_singlet_multichannel_weight
    integer :: ih,iproc
    real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real(kind=4) :: tBefore,tAfter
    complex(kind=8),dimension(pgl(ichan)%amps(iint)%n_amps,vector_size) :: amps
    logical :: done
    if (create_amplitude_library) then
       if (pgl(ichan)%amps(iint)%lib_created) return
    endif
    if (.not.allocated(val)) then
       if (keep_processes_separate) then
          allocate(val(1,vector_size))
          allocate(val_abs(1,vector_size))
       else
          allocate(val(1:maxval(pgl(1:ngroups)%nproc),vector_size))
          allocate(val_abs(1:maxval(pgl(1:ngroups)%nproc),vector_size))
       endif
       allocate(vol_ichan(1:ngroups,vector_size))
    endif
    ! some point-by-point initialisation
    f=0d0
    f_abs=0d0
    val=0d0
    val_abs=0d0

    ! Generate phase-space point based on the random numbers 'x(1:ndim)'
    call cpu_time(tBefore)
    !$omp parallel do
    do ivec=1,vector_size
       pgl(ichan)%ps(ivec)%x(:)=x(:,ivec)
       call pgl(ichan)%phase_space%generate_momenta(pgl(ichan)%ps(ivec))
!!$    enddo
!!$    if (debug) then
!!$       write (*,*) pgl(ichan)%ps(1:vector_size)%jac
!!$       stop 1
!!$    endif
!!$    do ivec=1,vector_size
       if (pgl(ichan)%ps(ivec)%jac.lt.0d0) then
          val(:,ivec)=-1d0
          cycle
       endif
       if (.not.pass_cuts(pgl(ichan),pgl(ichan)%ps(ivec))) then
          val(:,ivec)=-1d0
       endif
!!$    enddo
!!$    !$omp parallel do
!!$    do ivec=1,vector_size
       if (val(1,ivec).eq.-1d0) cycle
       call compute_multichannel_weight(ichan,pgl(ichan)%ps(ivec),colour_singlet_multichannel_weight(1,ivec))
!!$    enddo
!!$    call cpu_time(tAfter)
!!$    t_PS= t_PS +tAfter-tBefore
!!$    tBefore=tAfter
!!$    !$omp parallel do
!!$    do ivec=1,vector_size
       if (val(1,ivec).eq.-1d0) cycle
       call compute_the_amps(iint,ichan,ivec,amps(:,ivec))
    enddo
    call cpu_time(tAfter)
    t_amp=t_amp+tAfter-tBefore
    tBefore=tAfter
    do ivec=1,vector_size
       if (val(1,ivec).eq.-1d0) cycle
       call square_the_amps(iint,ichan,ivec,amps(:,ivec))
       if (.not. use_amplitude_library) pgl(ichan)%passed(iint) = pgl(ichan)%passed(iint) + 1
       if ((.not. use_amplitude_library) &
            .and. pgl(ichan)%passed(iint).le.nevent_hel_filter) then
          call optimise_the_amplitudes(iint,ichan,ivec,amps(1,ivec),done)
          if (done) return
       endif
    enddo
    
    ! MINT weight, phase-space jacobian and GeV -> pb conversion factor
    !$omp parallel do
    do ivec=1,vector_size
       if (val(1,ivec).eq.-1d0) cycle
       weight(ivec)=vol(ivec)*pgl(ichan)%ps(ivec)%jac*conv

       ! multiply by the strong coupling
       if (pgl(ichan)%amps(iint)%n_sing(1).lt.pgl(ichan)%next-2) then
          weight(ivec)=weight(ivec)*(4*pi*alphas)**(pgl(ichan)%next-2-pgl(ichan)%amps(iint)%n_sing(1))
       endif
    
       ! multiply by the EW coupling
       if (pgl(ichan)%amps(iint)%n_sing(1).ge.1) then
          weight(ivec)=weight(ivec)*(2d0*4d0*pi*alphaEW)**pgl(ichan)%amps(iint)%n_sing(1)
       endif
       if (keep_processes_separate) then
          val(1,ivec)=pgl(ichan)%amp2(1,ivec)*weight(ivec)/dble(pgl(ichan)%iden(iint))
          val(1,ivec)=val(1,ivec)*colour_singlet_multichannel_weight(iint,ivec)
          call include_PDF_and_identical_procs(val(:,ivec),val_abs(:,ivec),pgl(ichan),pgl(ichan)%ps(ivec),iint,ivec)
          f_abs(ivec)=sum(val_abs(1:1,ivec),dim=1)
          f(ivec)=sum(val(1:1,ivec),dim=1)
       else
          val(1:pgl(ichan)%nproc,ivec)=pgl(ichan)%amp2(1:pgl(ichan)%nproc,ivec)*weight(ivec)/dble(pgl(ichan)%iden(1:pgl(ichan)%nproc))
          val(1:pgl(ichan)%nproc,ivec)=val(1:pgl(ichan)%nproc,ivec)*colour_singlet_multichannel_weight(1:pgl(ichan)%nproc,ivec)
          call include_PDF_and_identical_procs(val(:,ivec),val_abs(:,ivec),pgl(ichan),pgl(ichan)%ps(ivec),-1,ivec)
          f_abs(ivec)=sum(val_abs(1:pgl(ichan)%nproc,ivec),dim=1)
          f(ivec)=sum(val(1:pgl(ichan)%nproc,ivec),dim=1)
       endif
    enddo
    call cpu_time(tAfter)
    t_mat=t_mat+tAfter-tBefore
  end subroutine integrand

  subroutine optimise_the_amplitudes(iint,ichan,ivec,amps,done)
    implicit none
    integer,intent(in) :: iint,ichan,ivec
    logical,intent(out) :: done
    real(kind=8), dimension(:),allocatable :: amp2_save
    complex(kind=8),dimension(pgl(ichan)%amps(iint)%n_amps),intent(inout) :: amps
    done=.false.
    call find_same_flavour(pgl(ichan),nevent_hel_filter,amps)
    call setup_helicity_filter(pgl(ichan),iint,ivec)
    if (pgl(ichan)%passed(iint).eq.nevent_hel_filter) then
       ! recompute the amplitudes (to make sure that helicities are
       ! all filled correctly). We can also check that they are
       ! consistent with the non-optimised ones.
       allocate(amp2_save(1:pgl(ichan)%nproc))
       amp2_save=pgl(ichan)%amp2(:,ivec)
       call compute_the_amps(iint,ichan,ivec,amps)
       call square_the_amps(iint,ichan,ivec,amps)
       if (use_cross_process_optimisation_of_currents) then
          call pgl(ichan)%amps(iint)%optimise_evaluation(pgl(ichan)%next)
          call compute_the_amps(iint,ichan,ivec,amps)
          call square_the_amps(iint,ichan,ivec,amps)
       endif
       if (any(abs(amp2_save-pgl(ichan)%amp2(:,ivec))/(amp2_save+pgl(ichan)%amp2(:,ivec)).gt.1d-8)) then
          write (*,*) 'Find same flavour and helicity filter give different matrix elements'
          write (*,*) amp2_save
          write (*,*) pgl(ichan)%amp2(:,ivec)
          write (*,*) ''
          write (*,*) abs(amp2_save-pgl(ichan)%amp2(:,ivec))/(amp2_save+pgl(ichan)%amp2(:,ivec))
          stop 1
       endif
       deallocate(amp2_save)
       if (create_amplitude_library) then
          call pgl(ichan)%amps(iint)%create_library(pgl(ichan)%next,pgl(ichan)%hel,&
               ichan,iint,phys_model,pgl(ichan)%ps(ivec)%p)
          pgl(ichan)%amps(iint)%lib_created=.true.
          done=.true.
       endif
    endif
  end subroutine optimise_the_amplitudes

  subroutine compute_the_amps(iint,ichan,ivec,amps)
    use amp_lib 
    implicit none
    integer,intent(in) :: iint,ichan,ivec
    complex(kind=8),dimension(pgl(ichan)%amps(iint)%n_amps),intent(out) :: amps
    if (.not. use_amplitude_library) then
       call pgl(ichan)%amps(iint)%evaluate(pgl(ichan)%next,pgl(ichan)%ps(ivec)%p,&
            pgl(ichan)%hel,read_proc_from_file,phys_model)
       amps=pgl(ichan)%amps(iint)%amps
    else
       call evaluate_amp(ichan,iint,pgl(ichan)%ps(ivec)%p,amps)
    endif
  end subroutine compute_the_amps
    
  subroutine square_the_amps(iint,ichan,ivec,amps)
    implicit none
    integer,intent(in) :: iint,ichan,ivec
    complex(kind=8),dimension(pgl(ichan)%amps(iint)%n_amps),intent(in) :: amps
    integer :: iproc,ih
    iproc=0
    pgl(ichan)%amp2(:,ivec)=0d0
    if (keep_processes_separate) then
       if (use_real_gluons) then
          write (*,*) 'Code no longer compatible with using real gluons'
          stop 1
       else
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih,ivec)=dble(amps(ih)*&
                  pgl(ichan)%col_fac(iint)*dconjg(amps(ih)))*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc,ivec)=pgl(ichan)%amp2(iproc,ivec)+pgl(ichan)%amp2_hel(ih,ivec)
          enddo
       endif
    else
       if (use_real_gluons) then
          write (*,*) 'Code no longer compatible with using real gluons'
          stop 1
       else
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih,ivec)=dble(amps(ih)*&
                  pgl(ichan)%col_fac(iproc)*dconjg(amps(ih)))*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc,ivec)=pgl(ichan)%amp2(iproc,ivec)+pgl(ichan)%amp2_hel(ih,ivec)
          enddo
       endif
    endif
  end subroutine square_the_amps
  
  subroutine setup_helicity_filter(pgl,iint,ivec)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer,intent(in) :: iint,ivec
    real(kind=8) :: max_value
    integer :: ih1,ih2,iproc1,iproc2
    if (.not.allocated(pgl%include_hel)) then
       allocate(pgl%include_hel(maxval(pgl%nhel),pgl%nproc))
       pgl%include_hel=0
    endif
    ! filter zero helicities and helicities that are identical
    max_value=maxval(pgl%amp2_hel(1:pgl%nhel(iint),ivec))
    do ih1=1,pgl%nhel(iint)
       if (pgl%include_hel(ih1,iint).ne.0) cycle
       if (pgl%amp2_hel(ih1,ivec)/max_value.gt.1d-12) then
          ! non-zero
          pgl%include_hel(ih1,iint)=1
       else
          cycle
       endif
       do ih2=ih1+1,pgl%nhel(iint)
          if (abs(pgl%amp2_hel(ih1,ivec)-pgl%amp2_hel(ih2,ivec))/abs(pgl%amp2_hel(ih1,ivec)+pgl%amp2_hel(ih2,ivec)).lt.1d-12) then
             ! identical value. Now check that they belong to the same process
             iproc1=1
             do while (iproc1.lt.pgl%nproc .and. (pgl%amps(iint)%iproc_start(iproc1+1)-ih1).le.0)
                iproc1=iproc1+1
             enddo
             iproc2=1
             do while (iproc2.lt.pgl%nproc .and. (pgl%amps(iint)%iproc_start(iproc2+1)-ih2).le.0)
                iproc2=iproc2+1
             enddo
             if (iproc1.ne.iproc2) cycle
             ! identical process
             pgl%include_hel(ih2,iint)=-ih1
             pgl%include_hel(ih1,iint)=pgl%include_hel(ih1,iint)+1
          endif
       enddo
    enddo

    if (pgl%passed(iint).lt.nevent_hel_filter) return

    do i=1,2
       do ih1=1,pgl%nhel(iint)
          if (pgl%amps(iint)%same_flavour_sum(ih1,i).le.0) cycle
          pgl%amps(iint)%same_flavour_sum_operation(ih1,i)=0
          do while (pgl%include_hel(pgl%amps(iint)%same_flavour_sum(ih1,i),iint).lt.0)
             pgl%amps(iint)%same_flavour_sum_operation(ih1,i)= &
                  ieor(pgl%amps(iint)%same_flavour_sum_operation(ih1,i),find_operation(pgl,iint,ih1,i))
             pgl%amps(iint)%same_flavour_sum(ih1,i)=-pgl%include_hel(pgl%amps(iint)%same_flavour_sum(ih1,i),iint)
          enddo
       enddo
    enddo

    ih2=0
    do ih1=1,pgl%nhel(iint)
       if (pgl%include_hel(ih1,iint).gt.0) ih2=ih2+1
    enddo

    call pgl%amps(iint)%filter_helicity(pgl%next,pgl%nhel(iint),pgl%include_hel(1,iint)) ! this updates 'nhel' and 'include_hel'
    pgl%hel_fac(1:pgl%nhel(iint),iint)=pgl%include_hel(1:pgl%nhel(iint),iint)
  end subroutine setup_helicity_filter

  integer function find_operation(pgl,iint,iamp,idau)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: iamp,idau,iint
    complex(kind=8) :: amp1,amp2
    if (use_real_gluons) then
       write (*,*) 'Find operation for same flavour sum only for complex amplitudes'
       stop 1
    endif
    amp1=pgl%amps(iint)%amps(pgl%amps(iint)%same_flavour_sum(iamp,idau))
    amp2=pgl%amps(iint)%amps(-pgl%include_hel(pgl%amps(iint)%same_flavour_sum(iamp,idau),iint))
    find_operation=0
    if (abs((abs(dble(amp1))-abs(dble(amp2)))).lt.abs((abs(dble(amp1))-abs(aimag(amp2))))) then
       ! real==real and iamag==iamag
       if (sign(1d0,dble(amp1)).ne.sign(1d0,dble(amp2))) find_operation=find_operation+1
       if (sign(1d0,aimag(amp1)).ne.sign(1d0,aimag(amp2))) find_operation=find_operation+2
    else
       ! real==iamag and iamag==real
       find_operation=find_operation+4
       if (sign(1d0,dble(amp1)).ne.sign(1d0,aimag(amp2))) find_operation=find_operation+1
       if (sign(1d0,aimag(amp1)).ne.sign(1d0,dble(amp2))) find_operation=find_operation+2
    endif
  end function find_operation
    
  subroutine get_run_arguments()
    use argument_parser
    implicit none
    integer :: argc
    integer :: i
    character(len=256) :: argv
    character(len=80) :: library
    integer(kind=8) iseed
    common /to_seed/iseed
    call parse_argument(filename,ncalls0,itmax,PS_choice,iseed,library,tag)

    logfile="Outputs/"//trim(adjustl(tag))//"log_file.txt"
    open(unit=99,file=logfile,status='unknown')

    if (library.eq.'none') then
       create_amplitude_library=.false.
       use_amplitude_library=.false.
    elseif (library.eq.'create') then
       create_amplitude_library=.true.
       use_amplitude_library=.false.
    elseif (library.eq.'use') then
       create_amplitude_library=.false.
       use_amplitude_library=.true.
       return
    endif

    if (PS_choice.ne.1 .and. PS_choice.ne.2 .and. PS_choice.ne.3 .and. PS_choice.ne.4) then
       write (*,*) 'PS_Choice modes only 1, 2, 3 or 4',PS_choice
       stop 1
    endif
  end subroutine get_run_arguments
  
end program matrix_integrate_QCD
