
program amplicol_generate
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
  use mg_checks
  use subtraction
  use cs_dipole_mappings, only: cs_dipole_topology
  use cs_lc_spin_dipoles, only: dipole_status_is_numerical
  use integrated_dipoles
  use integration_histograms, only: histogram_initialize,histogram_begin_point,&
       histogram_commit_point,histogram_finalize_iteration,histogram_write
  use integration_analysis, only: analysis_begin,analysis_fill,analysis_distinguishes_massless_qcd_flavours
  use real_subtraction_strata
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: n_tail_records=8,n_contribution_components=9
  type :: tail_record
     logical :: valid=.false.,real_pass=.false.
     integer :: ichan=0,iint=0,physical_iint=0,stratum=0,iteration=0
     integer(kind=8) :: point=0_8
     real(kind=8) :: score=0d0,residual_score=0d0,component_score=0d0
     real(kind=8) :: f=0d0,f_abs=0d0,grid_volume=0d0,phase_space_jacobian=0d0
     real(kind=8) :: matrix_real=0d0,dipole_sum=0d0,matrix_residual=0d0
     real(kind=8) :: common_weight=0d0,physical_factor=0d0,counterevent_scale=0d0
     real(kind=8) :: renormalization_scale=0d0,alpha_s_value=0d0
     real(kind=8) :: threshold_distance=-1d0
     real(kind=8),allocatable :: x(:),momenta(:,:),dipole_values(:),alpha_variables(:)
     integer,allocatable :: process(:),dipole_ijk(:,:),dipole_topology(:)
     logical,allocatable :: alpha_active(:),active(:),passes_cuts(:)
  end type tail_record
  integer :: iproc,iident,target_label
  real(kind=8) :: weight
  integer :: i
  real(kind=8),dimension(:),allocatable :: mass,width
  character(len=80) :: filename,real_filename,logfile,limit_logfile,tag
  character(len=256) :: input_filename,tail_logfile,tail_replay_output,&
       tail_residual_replay_output,tail_replay_file
  integer(kind=4) :: PS_choice
  integer,parameter :: nevent_hel_filter=10
  integer :: igroup
  logical,dimension(1) :: to_write
  integer,dimension(:),allocatable :: nintegrals
  integer,dimension(:),allocatable :: integration_ndim_extra
  integer,dimension(:,:),allocatable :: integration_adaptation_classes
  integer :: ichan,iint,itmax,ncalls0,iamp,nborn_groups,nreal_groups,born_flavour_scheme,real_flavour_scheme
  integer :: limit_point
  integer,parameter :: n_limit_points=100,n_limit_failures=5
  integer,dimension(:,:,:),allocatable :: soft_fail,soft_tested
  integer,dimension(:,:,:,:),allocatable :: collinear_fail,collinear_tested
  real(kind=8),dimension(:,:),allocatable :: limit_base
  real(kind=8),dimension(1) :: f,f_abs
  real(kind=8),dimension(n_contribution_components,1) :: f_aux
  logical :: done,iteration_finished,histogram_active,tail_convergence_ok
  real(kind=8),dimension(:,:),allocatable :: wgts
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  character(len=19) :: formatted
  logical :: create_amplitude_library,use_amplitude_library,read_momenta,limit_test,has_real_process
  logical :: limits_ok,replay_tail,tail_tracking_enabled
  logical :: timing_enabled,time_detail_point,time_point_sample
  integer(kind=8) :: timing_point
  real(kind=8) :: tLoopBefore,tLoopAfter,tSampleBefore,tSampleAfter,tFinalBefore,tFinalAfter
  real(kind=8) :: accuracy
  real(kind=8) :: migration_tail_fraction_limit
  character(len=80) :: dim_reg_scheme
  type(tail_record),allocatable :: tail_residual_records(:,:,:),tail_component_records(:,:,:)
  integer(kind=8),allocatable :: tail_npoints(:,:),tail_iteration_npoints(:,:)
  real(kind=8),allocatable :: tail_residual_sum2(:,:),tail_component_sum2(:,:)
  real(kind=8),allocatable :: tail_iteration_residual_sum2(:,:),tail_iteration_max_residual2(:,:)
  integer :: tail_iteration
  real(kind=8) :: tail_global_score,tail_global_residual_score,tail_last_migration_fraction
  real(kind=8) :: tail_last_migration_variance_proxy,tail_last_migration_max_proxy
  logical :: tail_last_migration_converged

  call get_run_arguments()
  call read_run_parameters(input_filename)
  call write_run_parameters(99)
  timing_enabled=timing_mode.ne.timing_none
  if (timing_enabled) call cpu_time(tTot_B)

  call phys_model%init_part()
  if (.not.use_amplitude_library) then
     call apply_final_state_widths_from_process_file(filename)
  endif
  call phys_model%init_vert()

  ! LHAPDF also supplies alpha_s when PDFs themselves are disabled.
  if (use_lhapdf) then
     call InitPDFsetbyname(trim(adjustl(lhapdfset)))
     call initPDF(lhapdf_member)
     call setlhaparm('SILENT')
  elseif (include_pdf) then
     call PDF_initialise(trim(adjustl(internal_pdf_grid)),lhapdf_member)
  endif

  if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
  if (use_amplitude_library) then
     call read_amplitude_lib()
  else
     if (has_real_process) then
        call count_process_groups(filename,nborn_groups)
        call count_process_groups(real_filename,nreal_groups)
        call allocate_process_groups(nborn_groups+nreal_groups)
        call read_processes_from_file(filename,0,born_flavour_scheme)
        call clear_process_file_metadata()
        call read_processes_from_file(real_filename,nborn_groups,real_flavour_scheme)
        if (born_flavour_scheme.eq.0 .or. real_flavour_scheme.eq.0) then
           write(*,*) 'ERROR: --real-process requires v2 process files with an nf= flavour-scheme header'
           write(*,*) ' Regenerate both files with process_list.py.'
           stop 1
        endif
        if (born_flavour_scheme.ne.real_flavour_scheme) then
           write(*,*) 'ERROR: Born and real process files use different flavour schemes:',&
                born_flavour_scheme,real_flavour_scheme
           stop 1
        endif
        if (.not.use_lhapdf .and. born_flavour_scheme.ne.5) then
           write(*,*) 'ERROR: non-FS5 --real-process runs require LHAPDF.'
           write(*,*) ' The internal alpha_s evolution is only supported for FS=5.'
           stop 1
        endif
        pgl(nborn_groups+1:ngroups)%is_subtracted_real=.true.
        call clear_process_file_metadata()
     else
        call read_processes_from_file(filename)
     endif
     do i=1,ngroups
       call setup_optimised_multichannel_weight_computation(pgl(i))
    enddo
 endif
  if (timing_mode.eq.timing_detailed) then
     call cpu_time(tAfter)
     t_Proc_init=t_Proc_init+tAfter-tBefore
  endif
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

     ! Use the cut-aware bounds for Born phase spaces.  The phase-space
     ! generator copies this module setting into each instance during init,
     ! allowing Born and real channels to retain different sampling domains
     ! in a single B+R-D integration.
     if (PS_choice.eq.1 .or. PS_choice.eq.4) then
        call set_use_soft_bounds_as_actual_limits(.not.pgl(igroup)%is_subtracted_real)
     endif

     allocate(mass(pgl(igroup)%next))
     allocate(width(pgl(igroup)%next))
     do iproc=1,pgl(igroup)%nproc
        do iident=1,pgl(igroup)%iden_iproc(iproc)
           do i=1,pgl(igroup)%next
              if (phys_model%get_mass(pgl(igroup)%processes(i,iproc)).ne.&
                   phys_model%get_mass(pgl(igroup)%iden_processes(i,iident,iproc))) then
                 write (*,*) 'Reduced amplitude changes the physical external-mass layout',&
                      igroup,iproc,iident,i
                 stop 1
              endif
           enddo
        enddo
     enddo
     do i=1,pgl(igroup)%next
        ! The phase-space generator uses the canonical density labels, while
        ! iden_processes stores physical subprocesses in the fixed matrix-
        ! element labelling.  Translate the former into the latter before
        ! assigning a mass.  The reduced amplitude representative itself can
        ! have a different crossing and must not be used for this purpose.
        target_label=pgl(igroup)%phase_space_permutations(i,1)
        mass(i)=phys_model%get_mass(pgl(igroup)%iden_processes(target_label,1,1))
        width(i)=phys_model%get_width(pgl(igroup)%iden_processes(target_label,1,1))
        do iproc=1,pgl(igroup)%nproc
           target_label=pgl(igroup)%phase_space_permutations(i,iproc)
           do iident=1,pgl(igroup)%iden_iproc(iproc)
              if (mass(i).ne.phys_model%get_mass(&
                   pgl(igroup)%iden_processes(target_label,iident,iproc)) .or. &
                   width(i).ne.phys_model%get_width(&
                   pgl(igroup)%iden_processes(target_label,iident,iproc))) then
                 write (*,*) 'masses and widths not compatible among physical subprocesses'
                 stop 1
              endif
           enddo
        enddo
     enddo
     ! Initialise the phase-space parametrisation
     if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
     call setup_cuts_for_each_particle(pgl(igroup),igroup)
     if (PS_choice.ge.1 .and. PS_choice.le.3) then
        call pgl(igroup)%phase_space%init(sqrts,pgl(igroup)%next,mass,pgl(igroup)%phase_space_orders,&
             pgl(igroup)%pt_min,pgl(igroup)%eta_max,pgl(igroup)%DR_min,pgl(igroup)%sqrt_s_min,.false.,include_pdf)
     elseif (PS_choice.eq.4) then
        call pgl(igroup)%phase_space%init(sqrts,pgl(igroup)%next,mass,pgl(igroup)%phase_space_orders,&
             pgl(igroup)%pt_min,pgl(igroup)%eta_max,pgl(igroup)%DR_min,pgl(igroup)%sqrt_s_min,.true.,include_pdf)
     endif
     pgl(igroup)%ndim_extra=pgl(igroup)%phase_space%ndim_extra
     allocate(pgl(igroup)%ps(1))
     allocate(pgl(igroup)%ps(1)%x(1:pgl(igroup)%ndim+pgl(igroup)%ndim_extra))
     allocate(pgl(igroup)%ps(1)%p(0:3,1:pgl(igroup)%next))
     if (timing_mode.eq.timing_detailed) then
        call cpu_time(tAfter)
        t_PS_init=t_PS_init+tAfter-tBefore
     endif
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
     if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
     if (keep_processes_separate) then
        do iamp=1,pgl(igroup)%nproc
           if (read_momenta) call run_madgraph_check(pgl(igroup)%next,igroup,iamp,pgl(igroup)%processes(1,iamp))
           call pgl(igroup)%amps(iamp)%init(1,pgl(igroup)%next,1,pgl(igroup)%processes(1,iamp),&
                pgl(igroup)%spin,pgl(igroup)%color_orders(1,iamp),phys_model)
           if (pgl(igroup)%is_subtracted_real .or. limit_test) call initialise_subtraction(igroup,iamp)
           if (read_momenta) then
              if (.not.allocated(p_read)) allocate(p_read(pgl(igroup)%next,0:3))
              call read_in_momenta(pgl(igroup)%next,igroup,iamp,p_read)
              do i=1,pgl(igroup)%next
                 pgl(igroup)%ps(1)%p(:,i)=p_read(i,:)
              enddo
           endif
        enddo
     else
        call pgl(igroup)%amps(1)%init(1,pgl(igroup)%next,pgl(igroup)%nproc,pgl(igroup)%processes,&
             pgl(igroup)%spin,pgl(igroup)%color_orders,phys_model)
     endif

     if (timing_mode.eq.timing_detailed) then
        call cpu_time(tAfter)
        t_amp_init=t_amp_init+tAfter-tBefore
     endif

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
        allocate(pgl(igroup)%amp2(1))
     else
        allocate(pgl(igroup)%amp2(1:pgl(igroup)%nproc))
     endif

     ! number of helicities to sum over
     allocate(pgl(igroup)%amp2_hel(1:maxval(pgl(igroup)%nhel)))
     pgl(igroup)%hel = pgl(igroup)%spin(1,1:pgl(igroup)%next)
     if (keep_processes_separate) then
        allocate(pgl(igroup)%hel_fac(1:maxval(pgl(igroup)%nhel),pgl(igroup)%nproc))
     else
        allocate(pgl(igroup)%hel_fac(1:maxval(pgl(igroup)%nhel),1))
     endif
     pgl(igroup)%hel_fac=1

  enddo ! loop over phase-space-order groups

  if (has_real_process) call initialise_integrated_dipoles(nborn_groups,dim_reg_scheme,born_flavour_scheme)

  if (use_amplitude_library) then
     if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
     call test_lib
     if (timing_mode.eq.timing_detailed) then
        call cpu_time(tAfter)
        t_lib_check=t_lib_check+tAfter-tBefore
     endif
  endif
  
  if (.not.has_real_process) then
     filename='Outputs/'//trim(adjustl(tag))//'events_tmp.lhe'
     open(unit=11,file=filename,action='readwrite',status='unknown')
     if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
     call write_unique_in_file(pgl_unique,unique_map,unique_map_value,abs(ncalls0))
  endif
  
  allocate(nintegrals(ngroups),integration_ndim_extra(ngroups))
  integration_ndim_extra=pgl(1:ngroups)%ndim_extra
  if (keep_processes_separate) then
     nintegrals(1:ngroups)=pgl(1:ngroups)%nproc
  else
     nintegrals(1:ngroups)=1
  endif
  if (has_real_process) then
     do igroup=1,nborn_groups
        nintegrals(igroup)=3*nintegrals(igroup)
        integration_ndim_extra(igroup)=integration_ndim_extra(igroup)+1
     enddo
     do igroup=nborn_groups+1,ngroups
        nintegrals(igroup)=n_real_strata*nintegrals(igroup)
     enddo
  endif
  allocate(integration_adaptation_classes(maxval(nintegrals),ngroups))
  integration_adaptation_classes=1
  if (has_real_process) then
     do igroup=nborn_groups+1,ngroups
        do iint=1,nintegrals(igroup)
           integration_adaptation_classes(iint,igroup)=mod(iint-1,n_real_strata)+1
        enddo
     enddo
  endif
  if (accuracy.gt.0d0) then
     call simple_integrator%init(ngroups,pgl(1:ngroups)%ndim,&
          integration_ndim_extra,nintegrals,abs(ncalls0),abs(itmax),accuracy,naux=n_contribution_components,&
          adaptation_classes=integration_adaptation_classes)
  else
     call simple_integrator%init(ngroups,pgl(1:ngroups)%ndim,&
          integration_ndim_extra,nintegrals,abs(ncalls0),abs(itmax),naux=n_contribution_components,&
          adaptation_classes=integration_adaptation_classes)
  endif
  tail_tracking_enabled=has_real_process .and. .not.replay_tail
  if (tail_tracking_enabled) call initialize_tail_diagnostics()
  histogram_active=accuracy.gt.0d0 .and. .not.limit_test .and. .not.read_momenta
  main_run: block
  if (replay_tail) then
     histogram_active=.false.
     call replay_saved_tail_point()
     exit main_run
  endif
  if (histogram_active) then
     filename='Outputs/'//trim(adjustl(tag))//'histograms.HwU'
     call histogram_initialize(nintegrals,has_real_process,filename)
     call analysis_begin()
  endif
  if (timing_mode.eq.timing_detailed) then
     call cpu_time(tAfter)
     t_Int_init=t_Int_init+tAfter-tBefore
  endif
  call date_and_time(date, time, zone)
  write(formatted, '(A4,"-",A2,"-",A2," ",A2,":",A2,":",A2)') &
       date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
  write (*,*) 'Start phase-space integration '//trim(formatted)
  write (99,*) 'Start phase-space integration '//trim(formatted)
  call flush(99)
  timing_point=0_8
  if (timing_enabled) then
     call cpu_time(tLoopBefore)
     t_Initialise=tLoopBefore-tTot_B
  endif
  do
     timing_point=timing_point+1_8
     time_point_sample=.false.
     time_detail_point=.false.
     if (timing_mode.eq.timing_detailed) then
        time_point_sample=mod(timing_point-1_8,int(timing_sample,kind=8)).eq.0_8
        time_detail_point=time_point_sample
     endif
     if (.not.limit_test) then
        if (time_detail_point) call cpu_time(tSampleBefore)
        call simple_integrator%get_points(1,ichan,iint)
        if (time_detail_point) then
           call cpu_time(tSampleAfter)
           t_Int_get=t_Int_get+(tSampleAfter-tSampleBefore)*dble(timing_sample)
        endif
     endif

     if (limit_test) then
        allocate(soft_fail(ngroups,maxval(pgl(1:ngroups)%nproc),pgl(1)%next))
        allocate(soft_tested(ngroups,maxval(pgl(1:ngroups)%nproc),pgl(1)%next))
        allocate(collinear_fail(ngroups,maxval(pgl(1:ngroups)%nproc),pgl(1)%next,pgl(1)%next))
        allocate(collinear_tested(ngroups,maxval(pgl(1:ngroups)%nproc),pgl(1)%next,pgl(1)%next))
        allocate(limit_base(0:3,pgl(1)%next))
        soft_fail=0
        soft_tested=0
        collinear_fail=0
        collinear_tested=0
        do limit_point=1,n_limit_points
           do igroup=1,ngroups
              call generate_limit_phase_space_point(igroup)
              limit_base=pgl(igroup)%ps(1)%p
              do iint=1,pgl(igroup)%nproc
                 ! Every integral must start from the same accepted hard point;
                 ! the preceding iint scan leaves the phase space deformed.
                 pgl(igroup)%ps(1)%p=limit_base
                 call test_limits_integrand(igroup,iint,limit_point,&
                      soft_fail(igroup,iint,:),soft_tested(igroup,iint,:),&
                      collinear_fail(igroup,iint,:,:),collinear_tested(igroup,iint,:,:),&
                      use_amplitude_library)
              enddo
           enddo
        enddo
        call print_limit_failure_fractions(soft_fail,soft_tested,collinear_fail,collinear_tested,&
             n_limit_failures,limits_ok)
        if (limits_ok) then
           write (*,'(a)') 'Limit diagnostic: PASS'
           write (99,'(a)') 'Limit diagnostic: PASS'
           write (*,'(a)') 'Failed-limit details: '//trim(limit_logfile)
           write (99,'(a)') 'Failed-limit details: '//trim(limit_logfile)
           call flush(100)
           stop 0
        else
           write (*,'(a)') 'Limit diagnostic: FAIL'
           write (99,'(a)') 'Limit diagnostic: FAIL'
           write (*,'(a)') 'Failed-limit details: '//trim(limit_logfile)
           write (99,'(a)') 'Failed-limit details: '//trim(limit_logfile)
           call flush(100)
           stop 1
        endif
     endif
     
     if (histogram_active) call histogram_begin_point(ichan,iint)
     call integrand(ichan,iint,simple_integrator%x(:,1),simple_integrator%wgt(1),&
          f(1),f_abs(1),f_aux(:,1))
     if (histogram_active) call histogram_commit_point()
     tail_convergence_ok=.true.
     if (tail_tracking_enabled) tail_convergence_ok=migration_tail_convergence_ok()
     if (time_detail_point) call cpu_time(tSampleBefore)
     call simple_integrator%fill_points(1,f_abs,f,to_write,done,f_aux=f_aux,&
          iteration_finished=iteration_finished,external_converged=tail_convergence_ok)
     if (histogram_active .and. iteration_finished) call histogram_finalize_iteration()
     if (tail_tracking_enabled .and. iteration_finished) then
        call finalize_tail_iteration()
        tail_iteration=tail_iteration+1
        call write_tail_diagnostics()
        call report_migration_tail_convergence()
        tail_iteration_npoints=0_8
        tail_iteration_residual_sum2=0d0
        tail_iteration_max_residual2=0d0
     endif
     if (time_detail_point) then
        call cpu_time(tSampleAfter)
        t_Int_fill=t_Int_fill+(tSampleAfter-tSampleBefore)*dble(timing_sample)
     endif
     if (create_amplitude_library) then
        done=.true.
        do igroup=1,ngroups
           done=done.and.all(pgl(igroup)%amps%lib_created)
        enddo
        if (done) then
           if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
           call create_amplitude_lib()
           if (timing_mode.eq.timing_detailed) then
              call cpu_time(tAfter)
              t_Amp_opt=t_Amp_opt+tAfter-tBefore
           endif
        endif
     endif
     if (to_write(1) .and. .not.has_real_process) then
        if (timing_mode.eq.timing_detailed) call cpu_time(tSampleBefore)
        call unwgt_process(pgl(ichan),iint) ! pick a random process
        call unwgt_helicity(pgl(ichan))     ! pick a random helicity for the process picked
        call write_event(11,pgl(ichan),1d0)
        if (timing_mode.eq.timing_detailed) then
           call cpu_time(tSampleAfter)
           t_Evt_write=t_Evt_write+tSampleAfter-tSampleBefore
        endif
     endif
     if (done) exit
  enddo
  if (histogram_active) call histogram_write()
  if (tail_tracking_enabled) call write_tail_diagnostics()
  if (timing_enabled) then
     call cpu_time(tLoopAfter)
     t_Int_loop=t_Int_loop+tLoopAfter-tLoopBefore
     call cpu_time(tFinalBefore)
  endif
  if (.not.has_real_process) then
     if (timing_mode.eq.timing_detailed) call cpu_time(tSampleBefore)
     call flush(11)
     call simple_integrator%assign_evnt_wgts(wgts)
     if (timing_mode.eq.timing_detailed) then
        call cpu_time(tSampleAfter)
        t_Evt_wgt_assign=t_Evt_wgt_assign+tSampleAfter-tSampleBefore
     endif
     if (timing_mode.eq.timing_detailed) call cpu_time(tSampleBefore)
     rewind(11)
     filename='Outputs/'//trim(adjustl(tag))//'events.lhe'
     open(unit=12,file=filename,action='write',status='unknown')
     write (*,*) 'Updating event weights...'
     write (99,*) 'Updating event weights...'
     do i=1,size(wgts,dim=2)
        call event_update_wgt(11,12,wgts(1,i))
     enddo
     close(11,status='DELETE')
     write(12,'(a)') '</LesHouchesEvents>'
     close(12)
     if (timing_mode.eq.timing_detailed) then
        call cpu_time(tSampleAfter)
        t_Evt_wgt_update=t_Evt_wgt_update+tSampleAfter-tSampleBefore
     endif
  else
     call print_contribution_results()
  endif
  if (timing_enabled) then
     call cpu_time(tFinalAfter)
     t_Finalise=t_Finalise+tFinalAfter-tFinalBefore
  endif
     
  if (timing_enabled) then
     call cpu_time(tTot_a)
     t_all=tTot_a-tTot_b
     if (timing_mode.eq.timing_detailed) then
        t_other=t_all-(t_PS_init+t_Amp_init+t_Proc_init+t_PS+t_Amp+t_mat+t_Amp_opt+t_weight+&
             t_lib_check+t_Int_init+t_Int_get+t_Int_fill+t_Evt_write+t_Evt_wgt_assign+t_Evt_wgt_update)
     endif
     call print_timing(6)
     call print_timing(99)
  endif
  end block main_run
  close(99)
  
contains

  subroutine print_contribution_results()
    implicit none
    real(kind=8),allocatable :: channel_res(:,:),channel_unc(:,:),aux_res(:,:),aux_unc(:,:)
    real(kind=8) :: component(n_contribution_components),component_unc(n_contribution_components),finite,finite_unc
    character(len=32),parameter :: labels(7)=[character(len=32) ::&
         'Born','Real - local dipoles','Integrated I coefficient -2',&
         'Integrated I coefficient -1','Integrated I finite',&
         'Integrated P','Integrated K']
    integer :: i,j
    call simple_integrator%get_channel_results(channel_res,channel_unc)
    call simple_integrator%get_channel_aux_results(aux_res,aux_unc)
    component=0d0
    component_unc=0d0
    do i=1,ngroups
       component=component+aux_res(:,i)
       component_unc=component_unc+aux_unc(:,i)**2
    enddo
    component_unc=sqrt(component_unc)
    do j=1,7
       write (*,'(a,2x,e14.7,1x,a,1x,e12.5)') trim(labels(j))//':',component(j),'+/-',component_unc(j)
       write (99,'(a,2x,e14.7,1x,a,1x,e12.5)') trim(labels(j))//':',component(j),'+/-',component_unc(j)
    enddo
    write (*,'(a,2x,e14.7,1x,a,1x,e12.5)') 'Real - local dipoles, regular:',component(8),'+/-',component_unc(8)
    write (99,'(a,2x,e14.7,1x,a,1x,e12.5)') 'Real - local dipoles, regular:',component(8),'+/-',component_unc(8)
    write (*,'(a,2x,e14.7,1x,a,1x,e12.5)') 'Real - local dipoles, migration:',component(9),'+/-',component_unc(9)
    write (99,'(a,2x,e14.7,1x,a,1x,e12.5)') 'Real - local dipoles, migration:',component(9),'+/-',component_unc(9)
    if (abs(component(2)-component(8)-component(9)).gt.&
         5d-12*max(1d0,abs(component(2)),abs(component(8))+abs(component(9)))) then
       write(*,*) 'ERROR: regular and migration strata do not close:',component(2),component(8),component(9)
       stop 1
    endif
    finite=sum(component((/1,2,5,6,7/)))
    finite_unc=sqrt(sum(component_unc((/1,2,5,6,7/))**2))
    write (*,'(a,2x,e14.7,1x,a,1x,e12.5)') 'Finite B+(R-D)+I0+P+K:',finite,'+/-',finite_unc
    write (99,'(a,2x,e14.7,1x,a,1x,e12.5)') 'Finite B+(R-D)+I0+P+K:',finite,'+/-',finite_unc
    deallocate(channel_res,channel_unc,aux_res,aux_unc)
  end subroutine print_contribution_results

  subroutine print_timing(iunit)
    implicit none
    integer,intent(in) :: iunit
    logical :: sampled
    character(len=32) :: label_fmt
    sampled=timing_mode.eq.timing_detailed .and. timing_sample.gt.1
    label_fmt=adjustl('Timing summary')
    write(iunit,'(a)') repeat('-',78)
    write(iunit,'(a32,2x,a14,2x,a9,2x,a)') label_fmt,'seconds','percent','note'
    write(iunit,'(a)') repeat('-',78)
    if (timing_mode.eq.timing_basic) then
       call print_timing_row(iunit,'initialisation',t_Initialise,'')
       call print_timing_row(iunit,'main integration loop',t_Int_loop,'')
       call print_timing_row(iunit,'finalisation',t_Finalise,'')
    else
       call print_timing_row(iunit,'phase-space initialisation',t_PS_init,'')
       call print_timing_row(iunit,'amplitude initialisation',t_Amp_init,'')
       call print_timing_row(iunit,'process initialisation',t_Proc_init,'')
       call print_timing_row(iunit,'integrator setup',t_Int_init,'')
       call print_timing_row(iunit,'phase-space generation',t_PS,estimate_note(sampled))
       call print_timing_row(iunit,'amplitude evaluation',t_Amp,estimate_note(sampled))
       call print_timing_row(iunit,'squaring amplitudes',t_mat,estimate_note(sampled))
       call print_timing_row(iunit,'amp optimisation/library',t_Amp_opt,estimate_note(sampled))
       call print_timing_row(iunit,'matrix weight/PDF',t_weight,estimate_note(sampled))
       call print_timing_row(iunit,'amplitude library checks',t_lib_check,'')
       call print_timing_row(iunit,'integrator point generation',t_Int_get,estimate_note(sampled))
       call print_timing_row(iunit,'integrator fill/grid update',t_Int_fill,estimate_note(sampled))
       call print_timing_row(iunit,'event writing',t_Evt_write,'')
       call print_timing_row(iunit,'event weight assignment',t_Evt_wgt_assign,'')
       call print_timing_row(iunit,'event weight update',t_Evt_wgt_update,'')
       call print_timing_row(iunit,'other/overhead',t_other,residual_note(sampled))
    endif
    write(iunit,'(a)') repeat('-',78)
    call print_timing_row(iunit,'total',t_all,'')
    if (sampled) then
       write(iunit,'(a,i0,a)') 'Note: rows marked "estimate" are extrapolated from every ',&
            timing_sample,'th phase-space point.'
       write(iunit,'(a)') '      other/overhead is the residual after subtracting those estimates.'
    endif
    write(iunit,'(a)') repeat('-',78)
  end subroutine print_timing

  subroutine print_timing_row(iunit,label,time,note)
    implicit none
    integer,intent(in) :: iunit
    character(len=*),intent(in) :: label
    real(kind=8),intent(in) :: time
    character(len=*),intent(in) :: note
    real(kind=8) :: pct
    character(len=32) :: label_fmt
    if (t_all.gt.0d0) then
       pct=100d0*time/t_all
    else
       pct=0d0
    endif
    label_fmt=adjustl(label)
    write(iunit,'(a32,2x,f14.6,2x,f8.2,a,2x,a)') label_fmt,time,pct,'%',trim(note)
  end subroutine print_timing_row

  character(len=8) function estimate_note(sampled)
    implicit none
    logical,intent(in) :: sampled
    estimate_note=''
    if (sampled) estimate_note='estimate'
  end function estimate_note

  character(len=8) function residual_note(sampled)
    implicit none
    logical,intent(in) :: sampled
    residual_note=''
    if (sampled) residual_note='residual'
  end function residual_note

  subroutine integrand(ichan,iint,x,vol,f,f_abs,f_components,replay_physical_factor,replay_counter_scale,&
       replay_unsplit)
    use scales
    use amp_lib
    implicit none
    integer,intent(in) :: ichan,iint
    real(kind=8), dimension(:),intent(in) :: x
    real(kind=8),intent(in) :: vol
    real(kind=8),intent(out) :: f,f_abs
    real(kind=8),intent(out) :: f_components(n_contribution_components)
    real(kind=8),intent(in),optional :: replay_physical_factor,replay_counter_scale
    logical,intent(in),optional :: replay_unsplit
    real(kind=8), dimension(:),allocatable,save :: val,val_abs,vol_ichan
    real(kind=8),dimension(pgl(ichan)%nproc) :: colour_singlet_multichannel_weight
    integer :: ih,iproc,a,target_label,eval_iint,integration_role,icopy,idip,requested_stratum,point_stratum
    integer :: resolution_info,resolution_dipole,kernel_info,kernel_dipole
    real(kind=8),dimension(0:3,pgl(ichan)%next) :: p_generated
    real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real(kind=8) :: amp2_integrand,amp2_dip,icoeff(-2:0),pterm,kterm,z,real_hist_factor,real_counter_scale
    real(kind=8) :: full_f,full_f_abs,regular_f,migration_f,real_margin,threshold_distance
    real(kind=8),allocatable :: dipole_values(:),real_hist_weights(:)
    real(kind=8),allocatable :: mapped_margins(:)
    real(kind=8),allocatable :: icoeff_copy(:,:),pterm_copy(:),kterm_copy(:)
    integer :: mapped_process(pgl(ichan)%next-1)
    real(kind=8) :: unresolved_invariant,resolution_tolerance
    real(kind=8),allocatable :: hard_copy(:)
    logical :: done,time_physics,real_pass,unsplit_real,stratum_selected
    logical,allocatable :: alpha_active_flags(:),mapped_pass_flags(:)
    real(kind=8),external :: alphaspdf
    time_physics=(timing_mode.eq.timing_detailed) .and. time_point_sample
    unsplit_real=.false.
    if (present(replay_unsplit)) unsplit_real=replay_unsplit
    integration_role=1
    requested_stratum=0
    eval_iint=iint
    if (has_real_process .and. pgl(ichan)%is_subtracted_real) then
       if (.not.unsplit_real) then
          requested_stratum=mod(iint-1,n_real_strata)+1
          eval_iint=(iint-1)/n_real_strata+1
       endif
    elseif (has_real_process) then
       integration_role=mod(iint-1,3)+1
       if (keep_processes_separate) then
          eval_iint=(iint-1)/3+1
       else
          eval_iint=1
       endif
    endif
    if (create_amplitude_library) then
       if (pgl(ichan)%amps(eval_iint)%lib_created) return
    endif
    if (.not.allocated(val)) then
       if (keep_processes_separate) then
          allocate(val(1))
          allocate(val_abs(1))
       else
          allocate(val(1:maxval(pgl(1:ngroups)%nproc)))
          allocate(val_abs(1:maxval(pgl(1:ngroups)%nproc)))
       endif
       allocate(vol_ichan(1:ngroups))
    endif
    ! some point-by-point initialisation
    f=0d0
    f_abs=0d0
    f_components=0d0
    val_abs=0d0
    if (keep_processes_separate) then
       iproc=eval_iint
    else
       iproc=1
    endif
    point_stratum=real_stratum_regular
    stratum_selected=.true.
    threshold_distance=-1d0
    if (tail_tracking_enabled .and. pgl(ichan)%is_subtracted_real) then
       tail_npoints(ichan,iint)=tail_npoints(ichan,iint)+1_8
       tail_iteration_npoints(ichan,iint)=tail_iteration_npoints(ichan,iint)+1_8
    endif

    ! Generate phase-space point based on the random numbers 'x(1:ndim)'
    if (time_physics) call cpu_time(tBefore)
    if (.not.read_momenta) then
       pgl(ichan)%ps(1)%x=x(1:size(pgl(ichan)%ps(1)%x))
       call pgl(ichan)%phase_space%generate_momenta(pgl(ichan)%ps(1))
    endif
    if (debug ) then
       write (*,*) pgl(ichan)%ps(1)%jac
       stop 1
    endif
    if (pgl(ichan)%ps(1)%jac.lt.0d0) then
       val=0d0
       if (time_physics) then
          call cpu_time(tAfter)
          t_PS= t_PS + (tAfter-tBefore)*dble(timing_sample)
       endif
       return
    endif
    if (.not.read_momenta) then
       ! The adaptive grid uses a canonical labelling of interchangeable
       ! massless-QCD final legs.  Restore the fixed coefficient labelling
       ! before cuts, amplitudes, scales, and event output.
       p_generated=pgl(ichan)%ps(1)%p
       do a=1,pgl(ichan)%next
          target_label=pgl(ichan)%phase_space_permutations(a,iproc)
          pgl(ichan)%ps(1)%p(:,target_label)=p_generated(:,a)
       enddo
    endif
    if (pgl(ichan)%is_subtracted_real) then
       call check_real_subtraction_resolution(eval_iint,ichan,resolution_info,&
            resolution_dipole,unresolved_invariant,resolution_tolerance)
       if (resolution_info.ne.0) then
          call report_subtraction_rejection('invariant',ichan,eval_iint,&
               resolution_dipole,resolution_info,unresolved_invariant,resolution_tolerance)
          return
       endif
    endif
    if (.not.read_momenta .and. .not.pgl(ichan)%is_subtracted_real) then
       if (.not.pass_cuts(pgl(ichan))) then
          val=0d0
          if (time_physics) then
             call cpu_time(tAfter)
             t_PS= t_PS + (tAfter-tBefore)*dble(timing_sample)
          endif
          return
       endif
    endif
    pgl(ichan)%passed(eval_iint) = pgl(ichan)%passed(eval_iint) + 1
    call compute_multichannel_weight(ichan,eval_iint,pgl(ichan)%ps(1),colour_singlet_multichannel_weight)
    if (time_physics) then
       call cpu_time(tAfter)
       t_PS= t_PS + (tAfter-tBefore)*dble(timing_sample)
       tBefore=tAfter
    endif
    call compute_the_amps(eval_iint,ichan,use_amplitude_library)
    if (time_physics) then
       call cpu_time(tAfter)
       t_amp=t_amp+(tAfter-tBefore)*dble(timing_sample)
       tBefore=tAfter
    endif
    call square_the_amps(eval_iint,ichan)
    if (time_physics) then
       call cpu_time(tAfter)
       t_mat=t_mat+(tAfter-tBefore)*dble(timing_sample)
    endif
    if ((.not. use_amplitude_library) &
         .and. pgl(ichan)%passed(eval_iint).le.nevent_hel_filter &
         .and. optimise_amplitudes) then
       if (time_physics) tBefore=tAfter
       call optimise_the_amplitudes(eval_iint,ichan,done)
       if (time_physics) then
          call cpu_time(tAfter)
          t_Amp_opt=t_Amp_opt+(tAfter-tBefore)*dble(timing_sample)
       endif
       if (done) return
    endif

    if (read_momenta) then
        call perform_check(eval_iint,ichan)
        if (pgl(ichan)%passed(eval_iint).gt.me_points) read_momenta=.false.
    endif

    real_pass=.true.
    if (pgl(ichan)%is_subtracted_real) real_pass=pass_real_subtracted_cuts(pgl(ichan),eval_iint)
    amp2_integrand=pgl(ichan)%amp2(1)
    if (pgl(ichan)%is_subtracted_real) then
       allocate(dipole_values(pgl(ichan)%dpl(eval_iint)%ndip))
       call evaluate_real_dipoles(eval_iint,ichan,amp2_dip,kernel_info,kernel_dipole,dipole_values)
       if (kernel_info.ne.0) then
          if (dipole_status_is_numerical(kernel_info)) then
             call report_subtraction_rejection('kernel',ichan,eval_iint,&
                  kernel_dipole,kernel_info,0d0,0d0)
             return
          endif
          write(*,*) 'ERROR: active local dipole kernel failed:',ichan,eval_iint,&
               kernel_dipole,kernel_info
          write(99,*) 'ERROR: active local dipole kernel failed:',ichan,eval_iint,&
               kernel_dipole,kernel_info
          stop 1
       endif
       allocate(alpha_active_flags(size(dipole_values)),mapped_pass_flags(size(dipole_values)))
       allocate(mapped_margins(size(dipole_values)))
       do idip=1,size(dipole_values)
          alpha_active_flags(idip)=pgl(ichan)%dpl(eval_iint)%dl(idip)%alpha_active
          mapped_pass_flags(idip)=pgl(ichan)%dpl(eval_iint)%dl(idip)%passes_cuts
          mapped_margins(idip)=huge(1d0)
          if (alpha_active_flags(idip)) then
             mapped_margins(idip)=mapped_dipole_jet_pt_margin(&
                  pgl(ichan)%dpl(eval_iint)%dl(idip)%p_mapped,&
                  pgl(ichan)%dpl(eval_iint)%dl(idip)%process_r)
          endif
       enddo
       point_stratum=classify_real_subtraction_stratum(real_pass,alpha_active_flags,mapped_pass_flags)
       if (point_stratum.eq.real_stratum_migration) then
          real_margin=real_subtracted_jet_pt_margin(pgl(ichan),eval_iint)
          threshold_distance=migration_pt_distance(real_margin,mapped_margins,real_pass,&
               alpha_active_flags,mapped_pass_flags)
       endif
       stratum_selected=unsplit_real .or. requested_stratum.eq.point_stratum
       if (real_pass) then
          amp2_integrand=amp2_integrand-amp2_dip
       else
          amp2_integrand=-amp2_dip
       endif
    endif

    ! set scales and update alphaS
    if (time_physics) call cpu_time(tBefore)
    call set_scale(scale_choice,pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
         pgl(ichan)%iden_processes(:,1,iproc),scale_ren)
    scale_fac=scale_ren
    scale_shower=scale_ren
    if (use_lhapdf) then
       alphas=alphaspdf(scale_ren)
    else
       alphas=alphas_Q(scale_ren,2,alphas_MZ)
    endif
    
    ! MINT weight, phase-space jacobian and GeV -> pb conversion factor
    weight=vol*pgl(ichan)%ps(1)%jac*conv

    ! multiply by the strong coupling
    if (pgl(ichan)%amps(eval_iint)%n_sing(1).lt.pgl(ichan)%next-2) then
       weight=weight*(4*pi*alphas)**(pgl(ichan)%next-2-pgl(ichan)%amps(eval_iint)%n_sing(1))
    endif
    
    ! multiply by the EW coupling
    if (pgl(ichan)%amps(eval_iint)%n_sing(1).ge.1) then
       weight=weight*(2d0*4d0*pi*alphaEW)**pgl(ichan)%amps(eval_iint)%n_sing(1)
    endif

    if (keep_processes_separate) then
       if (pgl(ichan)%is_subtracted_real) then
          ! R and every local counterevent have the same phase-space,
          ! coupling, PDF, identity, and multichannel prefactor.  Evaluate
          ! that factor once instead of evolving the PDFs for every dipole.
          ! This path is also used by deterministic tail replay.
          allocate(real_hist_weights(pgl(ichan)%iden_iproc(eval_iint)))
          call physical_matrix_weight(ichan,eval_iint,1d0,weight,&
               colour_singlet_multichannel_weight(eval_iint),real_hist_weights)
          real_hist_factor=sum(real_hist_weights)
          real_counter_scale=sum(abs(real_hist_weights))
          if (present(replay_physical_factor) .or. present(replay_counter_scale)) then
             if (.not.(present(replay_physical_factor) .and. present(replay_counter_scale))) then
                write(*,*) 'ERROR: tail replay requires both saved real-emission sampling factors'
                stop 1
             endif
             real_hist_factor=replay_physical_factor
             real_counter_scale=replay_counter_scale
          endif
          f=amp2_integrand*real_hist_factor
          f_abs=abs(amp2_integrand)*real_counter_scale
       else
          val(1)=amp2_integrand*weight/dble(pgl(ichan)%iden(eval_iint))
          val(1)=val(1)*colour_singlet_multichannel_weight(eval_iint)
          allocate(hard_copy(pgl(ichan)%iden_iproc(eval_iint)))
          hard_copy=val(1)*pgl(ichan)%idenCOandMAPfactor(&
               1:pgl(ichan)%iden_iproc(eval_iint),eval_iint)
          call include_PDF_and_identical_procs(val,val_abs,pgl(ichan),eval_iint)
          f_abs=sum(val_abs(1:1))
          f=sum(val(1:1))
       endif
    else
       val(1:pgl(ichan)%nproc)=pgl(ichan)%amp2(1:pgl(ichan)%nproc)*weight/dble(pgl(ichan)%iden(1:pgl(ichan)%nproc))
       val(1:pgl(ichan)%nproc)=val(1:pgl(ichan)%nproc)*colour_singlet_multichannel_weight(1:pgl(ichan)%nproc)
       call include_PDF_and_identical_procs(val,val_abs,pgl(ichan),-1)
       f_abs=sum(val_abs(1:pgl(ichan)%nproc))
       f=sum(val(1:pgl(ichan)%nproc))
    endif
    if (pgl(ichan)%is_subtracted_real) then
       full_f=f
       full_f_abs=f_abs
       call split_real_subtraction_weight(full_f,point_stratum,regular_f,migration_f)
       if (regular_f+migration_f.ne.full_f) then
          write(*,*) 'ERROR: real-subtraction strata do not close pointwise:',full_f,regular_f,migration_f
          stop 1
       endif
       if (.not.unsplit_real) then
          if (requested_stratum.eq.real_stratum_regular) then
             f=regular_f
          elseif (requested_stratum.eq.real_stratum_migration) then
             f=migration_f
          else
             write(*,*) 'ERROR: invalid real-subtraction stratum:',requested_stratum
             stop 1
          endif
          if (.not.stratum_selected) f_abs=0d0
          if (stratum_selected) f_abs=full_f_abs
       endif
       f_components(2)=f
       if (point_stratum.eq.real_stratum_regular .and. stratum_selected) f_components(8)=f
       if (point_stratum.eq.real_stratum_migration .and. stratum_selected) f_components(9)=f
       if (histogram_active .and. stratum_selected) then
          if (real_pass) then
             if (analysis_distinguishes_massless_qcd_flavours) then
                do icopy=1,pgl(ichan)%iden_iproc(eval_iint)
                   call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
                        pgl(ichan)%iden_processes(:,icopy,eval_iint),&
                        pgl(ichan)%amp2(1)*real_hist_weights(icopy),0d0)
                enddo
             else
                call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
                     pgl(ichan)%processes(:,eval_iint),pgl(ichan)%amp2(1)*real_hist_factor,0d0)
             endif
          endif
          do idip=1,size(dipole_values)
             if (.not.pgl(ichan)%dpl(eval_iint)%dl(idip)%active) cycle
             if (analysis_distinguishes_massless_qcd_flavours) then
                do icopy=1,pgl(ichan)%iden_iproc(eval_iint)
                   call make_identical_reduced_process(pgl(ichan)%processes(:,eval_iint),&
                        pgl(ichan)%iden_processes(:,icopy,eval_iint),&
                        pgl(ichan)%dpl(eval_iint)%dl(idip)%process_r,mapped_process)
                   call analysis_fill(pgl(ichan)%next-1,pgl(ichan)%dpl(eval_iint)%dl(idip)%p_mapped,&
                        mapped_process,-dipole_values(idip)*real_hist_weights(icopy),0d0)
                enddo
             else
                call analysis_fill(pgl(ichan)%next-1,pgl(ichan)%dpl(eval_iint)%dl(idip)%p_mapped,&
                     pgl(ichan)%dpl(eval_iint)%dl(idip)%process_r,&
                     -dipole_values(idip)*real_hist_factor,0d0)
             endif
          enddo
       endif
       if (tail_tracking_enabled) then
          call consider_tail_point(ichan,iint,eval_iint,point_stratum,x,vol,f,f_abs,pgl(ichan)%amp2(1),&
               amp2_dip,amp2_integrand,weight,real_hist_factor,real_counter_scale,real_pass,dipole_values,&
               threshold_distance,stratum_selected)
       endif
    elseif (has_real_process) then
       if (integration_role.eq.1) then
          f_components(1)=f
          if (keep_processes_separate) then
             if (histogram_active .and. analysis_distinguishes_massless_qcd_flavours) then
                allocate(icoeff_copy(-2:0,pgl(ichan)%iden_iproc(eval_iint)))
                call integrated_endpoint(ichan,eval_iint,&
                     pgl(ichan)%val_procs(1:pgl(ichan)%iden_iproc(eval_iint),eval_iint),&
                     pgl(ichan)%ps(1)%p,scale_ren,alphas,icoeff,icoeff_copy)
             else
                call integrated_endpoint(ichan,eval_iint,&
                     pgl(ichan)%val_procs(1:pgl(ichan)%iden_iproc(eval_iint),eval_iint),&
                     pgl(ichan)%ps(1)%p,scale_ren,alphas,icoeff)
             endif
          else
             icoeff=0d0
             do iproc=1,pgl(ichan)%nproc
                call add_endpoint_for_process(ichan,iproc,scale_ren,alphas,icoeff)
             enddo
          endif
          f_components(3:5)=icoeff(-2:0)
          f=f+icoeff(0)
          f_abs=abs(f_components(1))+abs(icoeff(0))
          if (histogram_active) then
             if (analysis_distinguishes_massless_qcd_flavours) then
                do icopy=1,pgl(ichan)%iden_iproc(eval_iint)
                   call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
                        pgl(ichan)%iden_processes(:,icopy,eval_iint),&
                        pgl(ichan)%val_procs(icopy,eval_iint)+icoeff_copy(0,icopy),&
                        pgl(ichan)%val_procs(icopy,eval_iint))
                enddo
             else
                call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
                     pgl(ichan)%processes(:,eval_iint),f,f_components(1))
             endif
          endif
       else
          z=x(size(x))
          pterm=0d0
          kterm=0d0
          if (keep_processes_separate) then
             if (histogram_active .and. analysis_distinguishes_massless_qcd_flavours) then
                allocate(pterm_copy(size(hard_copy)),kterm_copy(size(hard_copy)))
                call integrated_beam(ichan,eval_iint,integration_role-1,z,hard_copy,&
                     pgl(ichan)%ps(1)%xbjrk,scale_ren,scale_fac,alphas,&
                     pterm,kterm,pterm_copy,kterm_copy)
             else
                call integrated_beam(ichan,eval_iint,integration_role-1,z,hard_copy,&
                     pgl(ichan)%ps(1)%xbjrk,scale_ren,scale_fac,alphas,pterm,kterm)
             endif
          endif
          f_components(6)=pterm
          f_components(7)=kterm
          f=pterm+kterm
          f_abs=abs(pterm)+abs(kterm)
          if (histogram_active) then
             if (analysis_distinguishes_massless_qcd_flavours) then
                do icopy=1,pgl(ichan)%iden_iproc(eval_iint)
                   call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
                        pgl(ichan)%iden_processes(:,icopy,eval_iint),&
                        pterm_copy(icopy)+kterm_copy(icopy),0d0)
                enddo
             else
                call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
                     pgl(ichan)%processes(:,eval_iint),f,0d0)
             endif
          endif
       endif
    elseif (histogram_active) then
       if (keep_processes_separate) then
          if (analysis_distinguishes_massless_qcd_flavours) then
             do icopy=1,pgl(ichan)%iden_iproc(eval_iint)
                call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
                     pgl(ichan)%iden_processes(:,icopy,eval_iint),&
                     pgl(ichan)%val_procs(icopy,eval_iint),pgl(ichan)%val_procs(icopy,eval_iint))
             enddo
          else
             call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
                  pgl(ichan)%processes(:,eval_iint),f,f)
          endif
       else
          call analysis_fill(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
               pgl(ichan)%processes(:,eval_iint),f,f)
       endif
    endif
    if (allocated(hard_copy)) deallocate(hard_copy)
    if (time_physics) then
       call cpu_time(tAfter)
       t_weight=t_weight+(tAfter-tBefore)*dble(timing_sample)
    endif
  end subroutine integrand

  subroutine physical_matrix_weight(ichan,iint,matrix_element,common_weight,channel_weight,result)
    integer,intent(in) :: ichan,iint
    real(kind=8),intent(in) :: matrix_element,common_weight,channel_weight
    real(kind=8),intent(out) :: result(:)
    real(kind=8) :: value(1),value_abs(1)
    value(1)=matrix_element*common_weight/dble(pgl(ichan)%iden(iint))
    value(1)=value(1)*channel_weight
    call include_PDF_and_identical_procs(value,value_abs,pgl(ichan),iint)
    if (size(result).ne.pgl(ichan)%iden_iproc(iint)) then
       write(*,*) 'ERROR: physical matrix histogram weight array has incompatible size'
       stop 1
    endif
    result=pgl(ichan)%val_procs(1:pgl(ichan)%iden_iproc(iint),iint)
  end subroutine physical_matrix_weight

  subroutine initialize_tail_diagnostics()
    implicit none
    integer :: max_integrals,unit
    max_integrals=maxval(nintegrals)
    allocate(tail_residual_records(n_tail_records,ngroups,max_integrals))
    allocate(tail_component_records(n_tail_records,ngroups,max_integrals))
    allocate(tail_npoints(ngroups,max_integrals))
    allocate(tail_iteration_npoints(ngroups,max_integrals))
    allocate(tail_residual_sum2(ngroups,max_integrals),tail_component_sum2(ngroups,max_integrals))
    allocate(tail_iteration_residual_sum2(ngroups,max_integrals))
    allocate(tail_iteration_max_residual2(ngroups,max_integrals))
    tail_residual_records%valid=.false.
    tail_component_records%valid=.false.
    tail_npoints=0_8
    tail_iteration_npoints=0_8
    tail_residual_sum2=0d0
    tail_component_sum2=0d0
    tail_iteration_residual_sum2=0d0
    tail_iteration_max_residual2=0d0
    tail_iteration=0
    tail_global_score=-1d0
    tail_global_residual_score=-1d0
    tail_last_migration_fraction=-1d0
    tail_last_migration_variance_proxy=0d0
    tail_last_migration_max_proxy=0d0
    tail_last_migration_converged=.true.
    timing_point=0_8
    tail_logfile='Outputs/'//trim(adjustl(tag))//'tail_diagnostics.log'
    tail_replay_output='Outputs/'//trim(adjustl(tag))//'tail_replay.dat'
    tail_residual_replay_output='Outputs/'//trim(adjustl(tag))//'tail_residual_replay.dat'
    ! A run that terminates before accepting a real point must not leave a
    ! replay fixture from an older run with the same tag.
    if (.not.replay_tail) then
       open(newunit=unit,file=trim(tail_replay_output),status='replace',action='write')
       close(unit,status='delete')
       open(newunit=unit,file=trim(tail_residual_replay_output),status='replace',action='write')
       close(unit,status='delete')
    endif
    call write_tail_diagnostics()
  end subroutine initialize_tail_diagnostics

  subroutine consider_tail_point(ichan,iint,physical_iint,stratum,x,vol,f_point,f_abs_point,matrix_real,dipole_sum,&
       matrix_residual,common_weight,physical_factor,counterevent_scale,real_pass,dipole_values,&
       threshold_distance,stratum_selected)
    implicit none
    integer,intent(in) :: ichan,iint,physical_iint,stratum
    real(kind=8),intent(in) :: x(:),vol,f_point,f_abs_point,matrix_real,dipole_sum,matrix_residual
    real(kind=8),intent(in) :: common_weight,physical_factor,counterevent_scale,dipole_values(:)
    real(kind=8),intent(in) :: threshold_distance
    logical,intent(in) :: real_pass,stratum_selected
    integer :: residual_slot,component_slot
    real(kind=8) :: score,residual_score,component_score

    residual_score=abs(f_point)
    component_score=0d0
    if (stratum_selected) then
       if (real_pass) component_score=abs(matrix_real)*counterevent_scale
       if (size(dipole_values).gt.0) then
          component_score=max(component_score,maxval(abs(dipole_values))*counterevent_scale)
       endif
    endif
    score=max(residual_score,component_score)
    if (.not.ieee_is_finite(score)) then
       write(*,*) 'ERROR: non-finite subtracted-real tail score:',ichan,iint,score
       write(99,*) 'ERROR: non-finite subtracted-real tail score:',ichan,iint,score
       stop 1
    endif

    tail_residual_sum2(ichan,iint)=tail_residual_sum2(ichan,iint)+f_point*f_point
    tail_component_sum2(ichan,iint)=tail_component_sum2(ichan,iint)+component_score*component_score
    tail_iteration_residual_sum2(ichan,iint)=tail_iteration_residual_sum2(ichan,iint)+f_point*f_point
    tail_iteration_max_residual2(ichan,iint)=max(tail_iteration_max_residual2(ichan,iint),f_point*f_point)
    if (score.le.0d0) return

    residual_slot=find_tail_record_slot(tail_residual_records(:,ichan,iint),residual_score,.true.)
    component_slot=find_tail_record_slot(tail_component_records(:,ichan,iint),component_score,.false.)
    if (residual_slot.ne.0) then
       call save_tail_record(tail_residual_records(residual_slot,ichan,iint),ichan,iint,physical_iint,stratum,&
            x,vol,f_point,f_abs_point,&
            matrix_real,dipole_sum,matrix_residual,common_weight,physical_factor,counterevent_scale,&
            real_pass,dipole_values,residual_score,component_score,score,threshold_distance)
    endif
    if (component_slot.ne.0) then
       call save_tail_record(tail_component_records(component_slot,ichan,iint),ichan,iint,physical_iint,stratum,&
            x,vol,f_point,f_abs_point,&
            matrix_real,dipole_sum,matrix_residual,common_weight,physical_factor,counterevent_scale,&
            real_pass,dipole_values,residual_score,component_score,score,threshold_distance)
    endif
    if (score.gt.tail_global_score) then
       tail_global_score=score
       if (.not.replay_tail) then
          call write_tail_replay_point(tail_replay_output,ichan,iint,x,vol,f_point,f_abs_point,score,&
               physical_factor,counterevent_scale)
       endif
    endif
    if (residual_score.gt.tail_global_residual_score) then
       tail_global_residual_score=residual_score
       if (.not.replay_tail) then
          call write_tail_replay_point(tail_residual_replay_output,ichan,iint,x,vol,f_point,f_abs_point,&
               residual_score,physical_factor,counterevent_scale)
       endif
    endif
  end subroutine consider_tail_point

  logical function migration_tail_convergence_ok() result(converged)
    implicit none
    real(kind=8) :: fraction,total_proxy,max_proxy

    converged=.true.
    if (migration_tail_fraction_limit.le.0d0) return
    call current_migration_tail_fraction(fraction,total_proxy,max_proxy)
    if (total_proxy.gt.0d0) converged=fraction.le.migration_tail_fraction_limit
  end function migration_tail_convergence_ok

  subroutine current_migration_tail_fraction(fraction,total_proxy,max_proxy)
    implicit none
    real(kind=8),intent(out) :: fraction,total_proxy,max_proxy
    integer :: ichan_local,iint_local,stratum_local
    integer(kind=8) :: npoints_local
    real(kind=8) :: normalization

    total_proxy=0d0
    max_proxy=0d0
    do ichan_local=1,ngroups
       if (.not.pgl(ichan_local)%is_subtracted_real) cycle
       do iint_local=1,nintegrals(ichan_local)
          stratum_local=mod(iint_local-1,n_real_strata)+1
          if (stratum_local.ne.real_stratum_migration) cycle
          npoints_local=tail_iteration_npoints(ichan_local,iint_local)
          if (npoints_local.le.0_8) cycle
          normalization=dble(npoints_local)**2
          total_proxy=total_proxy+tail_iteration_residual_sum2(ichan_local,iint_local)/normalization
          max_proxy=max(max_proxy,tail_iteration_max_residual2(ichan_local,iint_local)/normalization)
       enddo
    enddo
    if (total_proxy.gt.0d0) then
       fraction=max_proxy/total_proxy
    else
       fraction=0d0
    endif
  end subroutine current_migration_tail_fraction

  subroutine finalize_tail_iteration()
    implicit none

    call current_migration_tail_fraction(tail_last_migration_fraction,&
         tail_last_migration_variance_proxy,tail_last_migration_max_proxy)
    tail_last_migration_converged=.true.
    if (migration_tail_fraction_limit.gt.0d0 .and. tail_last_migration_variance_proxy.gt.0d0) then
       tail_last_migration_converged=tail_last_migration_fraction.le.migration_tail_fraction_limit
    endif
  end subroutine finalize_tail_iteration

  subroutine report_migration_tail_convergence()
    implicit none
    character(len=16) :: status

    if (migration_tail_fraction_limit.le.0d0) then
       write(*,'(a)') 'Migration-tail convergence gate: disabled'
       write(99,'(a)') 'Migration-tail convergence gate: disabled'
       return
    endif
    if (tail_last_migration_converged) then
       status='PASS'
    else
       status='NOT SATISFIED'
    endif
    write(*,'(a,es12.4,a,es12.4,2a)') 'Migration largest-point variance fraction: ',&
         tail_last_migration_fraction,' limit ',migration_tail_fraction_limit,' -- ',trim(status)
    write(99,'(a,es12.4,a,es12.4,2a)') 'Migration largest-point variance fraction: ',&
         tail_last_migration_fraction,' limit ',migration_tail_fraction_limit,' -- ',trim(status)
    if (done .and. .not.tail_last_migration_converged) then
       write(*,'(a)') 'WARNING: iteration cap reached before migration-tail convergence.'
       write(99,'(a)') 'WARNING: iteration cap reached before migration-tail convergence.'
    endif
  end subroutine report_migration_tail_convergence

  integer function find_tail_record_slot(records,key,use_residual) result(slot)
    implicit none
    type(tail_record),intent(in) :: records(:)
    real(kind=8),intent(in) :: key
    logical,intent(in) :: use_residual
    integer :: irecord
    real(kind=8) :: record_key,min_key

    slot=0
    if (key.le.0d0) return
    min_key=huge(1d0)
    do irecord=1,size(records)
       if (.not.records(irecord)%valid) then
          slot=irecord
          return
       endif
       if (use_residual) then
          record_key=records(irecord)%residual_score
       else
          record_key=records(irecord)%component_score
       endif
       if (record_key.lt.min_key) then
          min_key=record_key
          slot=irecord
       endif
    enddo
    if (key.le.min_key) slot=0
  end function find_tail_record_slot

  subroutine save_tail_record(record,ichan,iint,physical_iint,stratum,x,vol,f_point,f_abs_point,matrix_real,dipole_sum,&
       matrix_residual,common_weight,physical_factor,counterevent_scale,real_pass,dipole_values,&
       residual_score,component_score,score,threshold_distance)
    implicit none
    type(tail_record),intent(inout) :: record
    integer,intent(in) :: ichan,iint,physical_iint,stratum
    real(kind=8),intent(in) :: x(:),vol,f_point,f_abs_point,matrix_real,dipole_sum,matrix_residual
    real(kind=8),intent(in) :: common_weight,physical_factor,counterevent_scale,dipole_values(:)
    real(kind=8),intent(in) :: residual_score,component_score,score
    real(kind=8),intent(in) :: threshold_distance
    logical,intent(in) :: real_pass
    integer :: idip,ndip

    record%valid=.true.
    record%real_pass=real_pass
    record%ichan=ichan
    record%iint=iint
    record%physical_iint=physical_iint
    record%stratum=stratum
    record%iteration=tail_iteration+1
    record%point=timing_point
    record%score=score
    record%residual_score=residual_score
    record%component_score=component_score
    record%f=f_point
    record%f_abs=f_abs_point
    record%grid_volume=vol
    record%phase_space_jacobian=pgl(ichan)%ps(1)%jac
    record%matrix_real=matrix_real
    record%dipole_sum=dipole_sum
    record%matrix_residual=matrix_residual
    record%common_weight=common_weight
    record%physical_factor=physical_factor
    record%counterevent_scale=counterevent_scale
    record%renormalization_scale=scale_ren
    record%alpha_s_value=alphas
    record%threshold_distance=threshold_distance
    record%x=x
    record%momenta=pgl(ichan)%ps(1)%p
    record%process=pgl(ichan)%processes(:,physical_iint)

    ndip=size(dipole_values)
    if (allocated(record%dipole_values)) deallocate(record%dipole_values)
    if (allocated(record%alpha_variables)) deallocate(record%alpha_variables)
    if (allocated(record%dipole_ijk)) deallocate(record%dipole_ijk)
    if (allocated(record%dipole_topology)) deallocate(record%dipole_topology)
    if (allocated(record%alpha_active)) deallocate(record%alpha_active)
    if (allocated(record%active)) deallocate(record%active)
    if (allocated(record%passes_cuts)) deallocate(record%passes_cuts)
    allocate(record%dipole_values(ndip),record%alpha_variables(ndip),record%dipole_ijk(3,ndip))
    allocate(record%dipole_topology(ndip),record%alpha_active(ndip),record%active(ndip),record%passes_cuts(ndip))
    record%dipole_values=dipole_values
    do idip=1,ndip
       record%alpha_variables(idip)=pgl(ichan)%dpl(physical_iint)%dl(idip)%alpha_variable
       record%dipole_ijk(:,idip)=pgl(ichan)%dpl(physical_iint)%dl(idip)%dip_ijk
       record%dipole_topology(idip)=cs_dipole_topology(pgl(ichan)%dpl(physical_iint)%dl(idip)%dip_ijk)
       record%alpha_active(idip)=pgl(ichan)%dpl(physical_iint)%dl(idip)%alpha_active
       record%active(idip)=pgl(ichan)%dpl(physical_iint)%dl(idip)%active
       record%passes_cuts(idip)=pgl(ichan)%dpl(physical_iint)%dl(idip)%passes_cuts
    enddo
  end subroutine save_tail_record

  subroutine write_tail_replay_point(output_file,ichan,iint,x,vol,f_point,f_abs_point,score,&
       physical_factor,counterevent_scale)
    implicit none
    character(len=*),intent(in) :: output_file
    integer,intent(in) :: ichan,iint
    real(kind=8),intent(in) :: x(:),vol,f_point,f_abs_point,score,physical_factor,counterevent_scale
    integer :: unit
    open(newunit=unit,file=trim(output_file),status='replace',action='write')
    write(unit,'(a)') '# AmpliCol tail replay v3'
    write(unit,*) ichan,iint,size(x),PS_choice,tail_iteration+1,timing_point
    write(unit,'(*(es25.16,1x))') vol,f_point,f_abs_point,score,physical_factor,counterevent_scale
    write(unit,'(*(es25.16,1x))') alpha_dipole
    write(unit,'(*(es25.16,1x))') x
    close(unit)
  end subroutine write_tail_replay_point

  subroutine replay_saved_tail_point()
    implicit none
    character(len=256) :: header
    integer :: unit,replay_ichan,replay_iint,nrandom,replay_ps,replay_iteration
    integer :: replay_physical_iint,replay_stratum
    integer(kind=8) :: replay_point
    real(kind=8) :: replay_vol,expected_f,expected_f_abs,expected_score,replay_alpha(4),tolerance
    real(kind=8) :: replay_physical_factor,replay_counter_scale
    real(kind=8),allocatable :: replay_x(:)
    logical :: replay_v2
    character(len=9),parameter :: stratum_name(2)=[character(len=9) :: 'regular','migration']

    open(newunit=unit,file=trim(tail_replay_file),status='old',action='read')
    read(unit,'(a)') header
    replay_v2=trim(header).eq.'# AmpliCol tail replay v2'
    if (.not.replay_v2 .and. trim(header).ne.'# AmpliCol tail replay v3') then
       write(*,*) 'ERROR: unsupported tail replay file: ',trim(header)
       stop 1
    endif
    read(unit,*) replay_ichan,replay_iint,nrandom,replay_ps,replay_iteration,replay_point
    read(unit,*) replay_vol,expected_f,expected_f_abs,expected_score,replay_physical_factor,replay_counter_scale
    read(unit,*) replay_alpha
    allocate(replay_x(nrandom))
    read(unit,*) replay_x
    close(unit)

    if (replay_ps.ne.PS_choice) then
       write(*,*) 'ERROR: tail replay phase-space choice differs:',replay_ps,PS_choice
       stop 1
    endif
    if (maxval(abs(replay_alpha-alpha_dipole)).gt.1d-14) then
       write(*,*) 'ERROR: tail replay alpha values differ:'
       write(*,*) ' saved:',replay_alpha
       write(*,*) ' active:',alpha_dipole
       stop 1
    endif
    if (replay_ichan.lt.1 .or. replay_ichan.gt.ngroups) then
       write(*,*) 'ERROR: tail replay channel is out of range:',replay_ichan
       stop 1
    endif
    if (.not.pgl(replay_ichan)%is_subtracted_real) then
       write(*,*) 'ERROR: tail replay channel is not subtracted real:',replay_ichan
       stop 1
    endif
    if (replay_v2) then
       if (replay_iint.lt.1 .or. replay_iint.gt.pgl(replay_ichan)%nproc) then
          write(*,*) 'ERROR: v2 tail replay physical integral is out of range:',replay_iint
          stop 1
       endif
    else
       if (replay_iint.lt.1 .or. replay_iint.gt.nintegrals(replay_ichan)) then
          write(*,*) 'ERROR: tail replay integration stratum is out of range:',replay_iint
          stop 1
       endif
    endif
    if (nrandom.ne.size(pgl(replay_ichan)%ps(1)%x)) then
       write(*,*) 'ERROR: tail replay random-coordinate count differs:',nrandom,size(pgl(replay_ichan)%ps(1)%x)
       stop 1
    endif

    if (replay_v2) then
       call integrand(replay_ichan,replay_iint,replay_x,replay_vol,f(1),f_abs(1),f_aux(:,1),&
            replay_physical_factor,replay_counter_scale,.true.)
    else
       call integrand(replay_ichan,replay_iint,replay_x,replay_vol,f(1),f_abs(1),f_aux(:,1),&
            replay_physical_factor,replay_counter_scale)
    endif
    tolerance=5d-11*max(1d0,abs(expected_f),abs(f(1)))
    write(*,'(a,2(i0,1x),a,i0,a,i0)') 'Tail replay channel/integral ',replay_ichan,replay_iint,&
         ' original iteration ',replay_iteration,' point ',replay_point
    if (.not.replay_v2) then
       replay_physical_iint=(replay_iint-1)/n_real_strata+1
       replay_stratum=mod(replay_iint-1,n_real_strata)+1
       write(*,'(a,i0,a,a)') 'Tail replay physical integral ',replay_physical_iint,&
            ' stratum ',trim(stratum_name(replay_stratum))
    endif
    write(*,'(a,es24.16)') 'Saved signed weight:   ',expected_f
    write(*,'(a,es24.16)') 'Replayed signed weight:',f(1)
    write(*,'(a,es24.16)') 'Saved absolute weight: ',expected_f_abs
    write(*,'(a,es24.16)') 'Saved tail score:      ',expected_score
    if (abs(f(1)-expected_f).gt.tolerance) then
       write(*,'(a,es12.4)') 'Tail replay: FAIL, tolerance ',tolerance
       stop 1
    endif
    write(*,'(a)') 'Tail replay: PASS'
    deallocate(replay_x)
  end subroutine replay_saved_tail_point

  subroutine write_tail_diagnostics()
    implicit none
    integer :: unit,ichan_local,iint_local,irecord,nvalid_residual,nvalid_component
    integer :: physical_iint_local,stratum_local
    real(kind=8) :: residual_top2,component_top2
    character(len=9),parameter :: stratum_name(2)=[character(len=9) :: 'regular','migration']

    open(newunit=unit,file=trim(tail_logfile),status='replace',action='write')
    write(unit,'(a)') '# AmpliCol subtracted-real tail diagnostics v2'
    write(unit,'(a,i0)') '# retained records per channel/integral: ',n_tail_records
    write(unit,'(a,4(1x,es16.8))') '# alpha FF FI IF II:',alpha_dipole
    write(unit,'(a)') '# score=max(abs(signed residual weight), largest measured R/D counterevent weight)'
    write(unit,'(a,1x,es16.8)') '# migration_tail_fraction_limit=',migration_tail_fraction_limit
    if (tail_iteration.gt.0) then
       write(unit,'(a,i0)') '# last_completed_iteration=',tail_iteration
       write(unit,'(a,3(1x,es16.8),a,l1)') '# migration variance_proxy max_point_proxy fraction',&
            tail_last_migration_variance_proxy,tail_last_migration_max_proxy,tail_last_migration_fraction,&
            ' converged ',tail_last_migration_converged
    endif
    write(unit,'(a)') '# component_replay_file='//trim(tail_replay_output)
    write(unit,'(a)') '# residual_replay_file='//trim(tail_residual_replay_output)
    do ichan_local=1,ngroups
       if (.not.pgl(ichan_local)%is_subtracted_real) cycle
       do iint_local=1,nintegrals(ichan_local)
          nvalid_residual=count(tail_residual_records(:,ichan_local,iint_local)%valid)
          nvalid_component=count(tail_component_records(:,ichan_local,iint_local)%valid)
          if (nvalid_residual.eq.0 .and. nvalid_component.eq.0) cycle
          residual_top2=0d0
          component_top2=0d0
          do irecord=1,n_tail_records
             if (tail_residual_records(irecord,ichan_local,iint_local)%valid) then
                residual_top2=residual_top2+tail_residual_records(irecord,ichan_local,iint_local)%f**2
             endif
             if (tail_component_records(irecord,ichan_local,iint_local)%valid) then
                component_top2=component_top2+&
                     tail_component_records(irecord,ichan_local,iint_local)%component_score**2
             endif
          enddo
          physical_iint_local=(iint_local-1)/n_real_strata+1
          stratum_local=mod(iint_local-1,n_real_strata)+1
          write(unit,'(a,2(1x,i0),a,i0,a,a,a,i0)') 'leaf channel/integral',ichan_local,iint_local,&
               ' physical ',physical_iint_local,' stratum ',trim(stratum_name(stratum_local)),&
               ' points ',tail_npoints(ichan_local,iint_local)
          write(unit,'(a,es24.16,a,es14.6)') ' residual_sum_w2 ',tail_residual_sum2(ichan_local,iint_local),&
               ' retained_fraction ',safe_fraction(residual_top2,tail_residual_sum2(ichan_local,iint_local))
          write(unit,'(a,es24.16,a,es14.6)') ' component_sum_w2 ',tail_component_sum2(ichan_local,iint_local),&
               ' retained_fraction ',safe_fraction(component_top2,tail_component_sum2(ichan_local,iint_local))
          write(unit,'(a)') ' residual_records'
          call write_tail_record_set(unit,tail_residual_records(:,ichan_local,iint_local),.true.)
          write(unit,'(a)') ' component_records'
          call write_tail_record_set(unit,tail_component_records(:,ichan_local,iint_local),.false.)
       enddo
    enddo
    close(unit)
  end subroutine write_tail_diagnostics

  subroutine write_tail_record_set(unit,records,use_residual)
    implicit none
    integer,intent(in) :: unit
    type(tail_record),intent(in) :: records(:)
    logical,intent(in) :: use_residual
    integer :: irecord,rank,nvalid,ibest
    logical :: used(size(records))
    real(kind=8) :: key,best_key

    nvalid=count(records%valid)
    used=.false.
    do rank=1,nvalid
       ibest=0
       best_key=-1d0
       do irecord=1,size(records)
          if (used(irecord) .or. .not.records(irecord)%valid) cycle
          if (use_residual) then
             key=records(irecord)%residual_score
          else
             key=records(irecord)%component_score
          endif
          if (key.gt.best_key) then
             ibest=irecord
             best_key=key
          endif
       enddo
       used(ibest)=.true.
       call write_tail_record(unit,rank,records(ibest))
    enddo
  end subroutine write_tail_record_set

  subroutine write_tail_record(unit,rank,record)
    implicit none
    integer,intent(in) :: unit,rank
    type(tail_record),intent(in) :: record
    integer :: ileg,idip
    character(len=2),parameter :: topology_name(4)=[character(len=2) :: 'FF','FI','IF','II']
    character(len=9),parameter :: stratum_name(2)=[character(len=9) :: 'regular','migration']

    write(unit,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') ' record rank ',rank,' iteration ',record%iteration,&
         ' point ',record%point,' channel ',record%ichan,' integral ',record%iint,&
         ' physical_integral ',record%physical_iint,' stratum ',trim(stratum_name(record%stratum))
    write(unit,'(a,3(1x,es24.16))') ' scores total residual component',record%score,&
         record%residual_score,record%component_score
    write(unit,'(a,2(1x,es24.16),a,l1)') ' weights signed absolute',record%f,record%f_abs,&
         ' real_pass ',record%real_pass
    write(unit,'(a,3(1x,es24.16))') ' matrix real dipole_sum residual',record%matrix_real,&
         record%dipole_sum,record%matrix_residual
    write(unit,'(a,3(1x,es24.16))') ' factors common physical counterevent_scale',record%common_weight,&
         record%physical_factor,record%counterevent_scale
    write(unit,'(a,4(1x,es24.16))') ' kinematics grid_volume ps_jac muR alphaS',record%grid_volume,&
         record%phase_space_jacobian,record%renormalization_scale,record%alpha_s_value
    write(unit,'(a,1x,es24.16)') ' threshold_distance_pt',record%threshold_distance
    write(unit,'(a,*(1x,es24.16))') ' x',record%x
    write(unit,'(a,*(1x,i0))') ' process',record%process
    do ileg=1,size(record%momenta,2)
       write(unit,'(a,1x,i0,4(1x,es24.16))') ' momentum',ileg,record%momenta(:,ileg)
    enddo
    do idip=1,size(record%dipole_values)
       write(unit,'(a,1x,i0,3(1x,i0),1x,a,2(1x,es24.16),3(1x,l1),1x,es24.16)')&
            ' dipole',idip,record%dipole_ijk(:,idip),topology_name(record%dipole_topology(idip)),&
            record%alpha_variables(idip),record%dipole_values(idip),record%alpha_active(idip),&
            record%active(idip),record%passes_cuts(idip),-record%dipole_values(idip)*record%physical_factor
    enddo
    write(unit,'(a)') ' end_record'
  end subroutine write_tail_record

  real(kind=8) function safe_fraction(numerator,denominator)
    implicit none
    real(kind=8),intent(in) :: numerator,denominator
    if (denominator.gt.0d0) then
       safe_fraction=numerator/denominator
    else
       safe_fraction=0d0
    endif
  end function safe_fraction

  subroutine make_identical_reduced_process(base_real,copy_real,base_reduced,copy_reduced)
    integer,intent(in) :: base_real(:),copy_real(:),base_reduced(:)
    integer,intent(out) :: copy_reduced(:)
    integer :: generic,ileg,target_flavour
    copy_reduced=base_reduced
    do generic=1,2
       target_flavour=0
       do ileg=1,size(base_real)
          if (abs(base_real(ileg)).eq.generic) then
             target_flavour=abs(copy_real(ileg))
             exit
          endif
       enddo
       if (target_flavour.eq.0) cycle
       do ileg=1,size(base_reduced)
          if (abs(base_reduced(ileg)).eq.generic) then
             copy_reduced(ileg)=sign(target_flavour,base_reduced(ileg))
          endif
       enddo
    enddo
  end subroutine make_identical_reduced_process

  subroutine report_subtraction_rejection(stage,ichan,iint,idip,status,invariant,tolerance)
    character(len=*), intent(in) :: stage
    integer, intent(in) :: ichan,iint,idip,status
    real(kind=8), intent(in) :: invariant,tolerance
    integer(kind=8), save :: invariant_rejections=0_8,kernel_rejections=0_8
    integer(kind=8) :: count
    logical :: report

    if (trim(stage).eq.'invariant') then
       invariant_rejections=invariant_rejections+1_8
       count=invariant_rejections
    else
       kernel_rejections=kernel_rejections+1_8
       count=kernel_rejections
    endif
    report=(count.le.10_8 .or. mod(count,1000_8).eq.0_8)
    if (report) then
       write(*,'(a,a,a,i0,a,4(i0,1x),a,2(es12.4,1x))')&
            'INFO: rejected numerically unresolved real point (',trim(stage),&
            '), count=',count,' channel/process/dipole/status=',&
            ichan,iint,idip,status,' invariant/tolerance=',invariant,tolerance
       write(99,'(a,a,a,i0,a,4(i0,1x),a,2(es12.4,1x))')&
            'INFO: rejected numerically unresolved real point (',trim(stage),&
            '), count=',count,' channel/process/dipole/status=',&
            ichan,iint,idip,status,' invariant/tolerance=',invariant,tolerance
    elseif (count.eq.11_8) then
       write(*,'(a,a,a)') 'INFO: further ',trim(stage),&
            ' numerical-rejection messages suppressed (counts reported every 1000)'
       write(99,'(a,a,a)') 'INFO: further ',trim(stage),&
            ' numerical-rejection messages suppressed (counts reported every 1000)'
    endif
  end subroutine report_subtraction_rejection

  subroutine add_endpoint_for_process(igroup,iproc,mu_ren,alpha_s,total)
    integer,intent(in) :: igroup,iproc
    real(kind=8),intent(in) :: mu_ren,alpha_s
    real(kind=8),intent(inout) :: total(-2:0)
    real(kind=8) :: one(-2:0)
    call integrated_endpoint(igroup,iproc,&
         pgl(igroup)%val_procs(1:pgl(igroup)%iden_iproc(iproc),iproc),&
         pgl(igroup)%ps(1)%p,mu_ren,alpha_s,one)
    total=total+one
  end subroutine add_endpoint_for_process

  subroutine optimise_the_amplitudes(iint,ichan,done)
    implicit none
    integer,intent(in) :: iint,ichan
    logical,intent(out) :: done
    real(kind=8), dimension(:),allocatable :: amp2_save
    done=.false.
    call find_same_flavour(pgl(ichan),nevent_hel_filter,pgl(ichan)%amp2)
    call setup_helicity_filter(pgl(ichan),iint)
    if (pgl(ichan)%passed(iint).eq.nevent_hel_filter) then
       ! recompute the amplitudes (to make sure that helicities are
       ! all filled correctly). We can also check that they are
       ! consistent with the non-optimised ones.
       allocate(amp2_save(1:pgl(ichan)%nproc))
       amp2_save=pgl(ichan)%amp2
       call compute_the_amps(iint,ichan,use_amplitude_library)
       call square_the_amps(iint,ichan)
       if (use_cross_process_optimisation_of_currents) then
          call pgl(ichan)%amps(iint)%optimise_evaluation(pgl(ichan)%next)
          call compute_the_amps(iint,ichan,use_amplitude_library)
          call square_the_amps(iint,ichan)
       endif
       if (any(abs(amp2_save-pgl(ichan)%amp2)/(amp2_save+pgl(ichan)%amp2).gt.1d-8)) then
          write (*,*) 'Find same flavour and helicity filter give different matrix elements',ichan,iint
          write (*,*) amp2_save
          write (*,*) pgl(ichan)%amp2
          write (*,*) ''
          write (*,*) abs(amp2_save-pgl(ichan)%amp2)/(amp2_save+pgl(ichan)%amp2)
          stop 1
       endif
       deallocate(amp2_save)
       if (create_amplitude_library) then
          call pgl(ichan)%amps(iint)%create_library(pgl(ichan)%next,pgl(ichan)%hel,&
               ichan,iint,phys_model,pgl(ichan)%ps(1)%p)
          pgl(ichan)%amps(iint)%lib_created=.true.
          done=.true.
       endif
       call flush(99)
    endif
  end subroutine optimise_the_amplitudes

  
  subroutine setup_helicity_filter(pgl,iint)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    real(kind=8) :: max_value
    integer :: ih1,ih2,iproc1,iproc2,iint
    if (.not.allocated(pgl%include_hel)) then
       allocate(pgl%include_hel(maxval(pgl%nhel),pgl%nproc))
       pgl%include_hel=0
    endif
    ! filter zero helicities and helicities that are identical
    max_value=maxval(pgl%amp2_hel(1:pgl%nhel(iint)))
    do ih1=1,pgl%nhel(iint)
       if (pgl%include_hel(ih1,iint).ne.0) cycle
       if (pgl%amp2_hel(ih1)/max_value.gt.1d-12) then
          ! non-zero
          pgl%include_hel(ih1,iint)=1
       else
          cycle
       endif
       do ih2=ih1+1,pgl%nhel(iint)
          if (abs(pgl%amp2_hel(ih1)-pgl%amp2_hel(ih2))/abs(pgl%amp2_hel(ih1)+pgl%amp2_hel(ih2)).lt.1d-12) then
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
!!$    deallocate(pgl%hel_fac)
!!$    allocate(pgl%hel_fac(pgl%nhel))
    pgl%hel_fac(1:pgl%nhel(iint),iint)=pgl%include_hel(1:pgl%nhel(iint),iint)
!!$    deallocate(pgl%include_hel)
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
    character(len=80) :: library,timing_arg
    integer :: timing_sample_arg
    integer(kind=8) iseed
    common /to_seed/iseed
    call parse_argument(filename,real_filename,input_filename,ncalls0,itmax,PS_choice,iseed,library,tag,&
         read_momenta,me_points,&
         limit_test,timing_arg,timing_sample_arg,accuracy,dim_reg_scheme,has_real_process,&
         tail_replay_file,replay_tail,migration_tail_fraction_limit)

    logfile="Outputs/"//trim(adjustl(tag))//"log_file.txt"
    open(unit=99,file=logfile,status='unknown')
    if (limit_test) then
       limit_logfile="Outputs/"//trim(adjustl(tag))//"limit_test_failures.log"
       open(unit=100,file=limit_logfile,action='write',status='replace')
       write (100,'(a)') '# Failed Catani-Seymour limit diagnostic records'
       write (100,'(a)') '# For boundedness checks, dipole and ratio are zero and residual is -1.'
       write (100,'(a,4(es12.4,1x))') '# alpha FF FI IF II: ',alpha_dipole
       write (100,'(a)') '# mapping_status -100: all matching real dipoles excluded by alpha'
       write (100,'(a)') '# mapping_status -101: all alpha-active matching dipoles fail mapped cuts'
    endif

    timing_arg=trim(adjustl(timing_arg))
    if (timing_arg.eq.'none') then
       timing_mode=timing_none
    elseif (timing_arg.eq.'basic') then
       timing_mode=timing_basic
    elseif (timing_arg.eq.'detailed') then
       timing_mode=timing_detailed
    else
       write (*,*) 'Timing mode must be none, basic, or detailed: ',trim(timing_arg)
       stop 1
    endif
    if (timing_sample_arg.lt.1) then
       write (*,*) 'Timing sample must be at least 1: ',timing_sample_arg
       stop 1
    endif
    timing_sample=timing_sample_arg

    dim_reg_scheme=trim(adjustl(dim_reg_scheme))
    if (dim_reg_scheme.ne.'hv' .and. dim_reg_scheme.ne.'fdh') then
       write (*,*) 'Dimensional regularization scheme must be hv or fdh: ',trim(dim_reg_scheme)
       stop 1
    endif

    if (has_real_process) then
       if (len_trim(real_filename).eq.0) then
          write (*,*) '--real-process requires a non-empty file name'
          stop 1
       endif
       if (.not.keep_processes_separate) then
          write (*,*) '--real-process requires keep_processes_separate=true'
          stop 1
       endif
       if (accuracy.le.0d0) then
          write (*,*) '--real-process is integration-only; provide --accuracy=X'
          stop 1
       endif
       if (library.ne.'none') then
          write (*,*) '--real-process is not available with --library=create or --library=use'
          stop 1
       endif
       if (limit_test .or. read_momenta) then
          write (*,*) '--real-process is not available with --limit_test or --me_test'
          stop 1
       endif
    elseif (replay_tail) then
       write (*,*) '--tail-replay requires --real-process=FILE'
       stop 1
    endif
    if (replay_tail .and. len_trim(tail_replay_file).eq.0) then
       write (*,*) '--tail-replay requires a non-empty file name'
       stop 1
    endif

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
    else
       write (*,*) 'library must be none, create or use: ',trim(library)
       stop 1
    endif

    if (PS_choice.ne.1 .and. PS_choice.ne.2 .and. PS_choice.ne.3 .and. PS_choice.ne.4) then
       write (*,*) 'PS_Choice modes only 1, 2, 3 or 4',PS_choice
       stop 1
    endif
  end subroutine get_run_arguments
  
end program amplicol_generate
