
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
  real(kind=8),parameter :: integration_value_limit=0.25d0*huge(1d0)**0.25d0
  real(kind=8),parameter :: integration_value_floor=sqrt(tiny(1d0))
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
  character(len=1024) :: filename,real_filename,logfile,limit_logfile
  character(len=80) :: tag
  character(len=256) :: input_filename,tail_logfile,tail_replay_output,&
       tail_residual_replay_output,tail_replay_file
  integer(kind=4) :: PS_choice
  integer,parameter :: nevent_hel_filter=10
  ! Infer reuse from nine points and validate it on the unseen tenth point.
  integer,parameter :: n_amplitude_optimisation_samples=nevent_hel_filter-1
  integer :: igroup
  logical,dimension(1) :: to_write
  integer,dimension(:),allocatable :: nintegrals
  integer,dimension(:),allocatable :: integration_ndim_extra
  integer,dimension(:,:),allocatable :: integration_adaptation_classes
  integer :: ichan,iint,itmax,ncalls0,iamp,nborn_groups,nreal_groups,born_flavour_scheme,real_flavour_scheme
  integer :: event_tmp_unit,event_output_unit,io_status
  character(len=256) :: io_message
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
  logical :: limits_ok,replay_tail,tail_tracking_enabled,real_requires_jets,scale_process_has_jet
  logical :: matrix_element_test_mode,matrix_element_checks_done
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
  event_tmp_unit=0
  matrix_element_test_mode=read_momenta
  matrix_element_checks_done=.false.
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
        real_requires_jets=.false.
        do igroup=nborn_groups+1,ngroups
           do iproc=1,pgl(igroup)%nproc
              if (real_subtracted_jet_requirement(pgl(igroup)%processes(:,iproc)).gt.0) &
                   real_requires_jets=.true.
           enddo
        enddo
        if (real_requires_jets .and. (pTj_min.le.0d0 .or. DRjj_min.le.0d0)) then
           write(*,*) 'ERROR: an NLO observable with resolved Born jets requires positive pTj_min and DRjj_min.'
           write(*,*) ' configured pTj_min/DRjj_min:',pTj_min,DRjj_min
           stop 1
        endif
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
     if (scale_choice.eq.4) then
        do iproc=1,pgl(igroup)%nproc
           scale_process_has_jet=.false.
           do i=3,pgl(igroup)%next
              if (phys_model%is_jet(pgl(igroup)%iden_processes(i,1,iproc))) &
                   scale_process_has_jet=.true.
           enddo
           if (.not.scale_process_has_jet) then
              write(*,*) 'ERROR: minimum-jet-pT scale requested for a subprocess without a jet:',&
                   igroup,iproc
              stop 1
           endif
        enddo
     endif
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

     if (use_amplitude_library) then
        if (matrix_element_test_mode) then
           do iamp=1,pgl(igroup)%nproc
              call run_madgraph_check(pgl(igroup)%next,igroup,iamp,&
                   pgl(igroup)%processes(1:pgl(igroup)%next,iamp))
           enddo
        endif
        cycle
     endif
     
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

  ! Accuracy-driven runs are integration-only: no unweighted-event candidates
  ! are produced, so neither create nor finalise a temporary LHE file.
  if (.not.has_real_process .and. .not.create_amplitude_library .and. &
       .not.matrix_element_test_mode .and. accuracy.le.0d0) then
     filename='Outputs/'//trim(adjustl(tag))//'events_tmp.lhe'
     io_message=''
     open(newunit=event_tmp_unit,file=trim(filename),action='readwrite',status='replace',&
          iostat=io_status,iomsg=io_message)
     if (io_status.ne.0) then
        write (*,*) 'Could not create temporary event file ',trim(filename),': ',&
             trim(io_message)
        stop 1
     endif
     if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
     call write_unique_in_file(event_tmp_unit,pgl_unique,unique_map,unique_map_value,abs(ncalls0))
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
     timing_point=checked_add8(timing_point,1_8,'integration point counter')
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
     if (matrix_element_test_mode .and. matrix_element_checks_done) then
        write (*,*) ''
        write (*,*) 'Passed all requested MadGraph matrix-element checks.'
        write (99,*) 'Passed all requested MadGraph matrix-element checks.'
        exit main_run
     endif
     if (histogram_active) call histogram_commit_point()
     tail_convergence_ok=.true.
     if (tail_tracking_enabled) tail_convergence_ok=migration_tail_convergence_ok()
     if (time_detail_point) call cpu_time(tSampleBefore)
     call simple_integrator%fill_points(1,f_abs,f,to_write,done,f_aux=f_aux,&
          iteration_finished=iteration_finished,external_converged=tail_convergence_ok)
     if (histogram_active .and. iteration_finished) call histogram_finalize_iteration()
     if (tail_tracking_enabled .and. iteration_finished) then
        call finalize_tail_iteration()
        tail_iteration=checked_add(tail_iteration,1,'tail-diagnostic iteration counter')
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
        call write_event(event_tmp_unit,pgl(ichan),1d0)
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
  if (.not.has_real_process .and. .not.create_amplitude_library .and. &
       .not.matrix_element_test_mode .and. accuracy.le.0d0) then
     if (timing_mode.eq.timing_detailed) call cpu_time(tSampleBefore)
     flush(event_tmp_unit,iostat=io_status,iomsg=io_message)
     if (io_status.ne.0) then
        write (*,*) 'Could not flush temporary event file: ',trim(io_message)
        stop 1
     endif
     call simple_integrator%assign_evnt_wgts(wgts)
     if (timing_mode.eq.timing_detailed) then
        call cpu_time(tSampleAfter)
        t_Evt_wgt_assign=t_Evt_wgt_assign+tSampleAfter-tSampleBefore
     endif
     if (timing_mode.eq.timing_detailed) call cpu_time(tSampleBefore)
     rewind(event_tmp_unit,iostat=io_status,iomsg=io_message)
     if (io_status.ne.0) then
        write (*,*) 'Could not rewind temporary event file: ',trim(io_message)
        stop 1
     endif
     filename='Outputs/'//trim(adjustl(tag))//'events.lhe'
     open(newunit=event_output_unit,file=trim(filename),action='write',status='replace',&
          iostat=io_status,iomsg=io_message)
     if (io_status.ne.0) then
        write (*,*) 'Could not create final event file ',trim(filename),': ',trim(io_message)
        stop 1
     endif
     write (*,*) 'Updating event weights...'
     write (99,*) 'Updating event weights...'
     do i=1,size(wgts,dim=2)
        call event_update_wgt(event_tmp_unit,event_output_unit,wgts(1,i),i.eq.1)
     enddo
     write(event_output_unit,'(a)',iostat=io_status,iomsg=io_message) '</LesHouchesEvents>'
     if (io_status.ne.0) then
        write (*,*) 'Could not finish final event file: ',trim(io_message)
        stop 1
     endif
     close(event_output_unit,iostat=io_status,iomsg=io_message)
     if (io_status.ne.0) then
        write (*,*) 'Could not close final event file: ',trim(io_message)
        stop 1
     endif
     ! Only discard the recoverable candidate file after the final file has
     ! been closed successfully.
     close(event_tmp_unit,status='delete',iostat=io_status,iomsg=io_message)
     if (io_status.ne.0) then
        write (*,*) 'Could not remove temporary event file: ',trim(io_message)
        stop 1
     endif
     if (timing_mode.eq.timing_detailed) then
        call cpu_time(tSampleAfter)
        t_Evt_wgt_update=t_Evt_wgt_update+tSampleAfter-tSampleBefore
     endif
  elseif (has_real_process) then
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
  close(99,iostat=io_status,iomsg=io_message)
  if (io_status.ne.0) then
     write (*,*) 'Could not close run log: ',trim(io_message)
     stop 1
  endif
  
contains

  logical function matrix_element_checks_complete()
    implicit none
    integer :: group_index

    matrix_element_checks_complete=.false.
    if (me_points.lt.1) then
       write (*,*) 'Invalid requested matrix-element check count:',me_points
       stop 1
    endif
    do group_index=1,ngroups
       if (pgl(group_index)%nproc.lt.1) then
          write (*,*) 'Invalid empty process group during matrix-element checks:',&
               group_index
          stop 1
       endif
       if (.not.allocated(pgl(group_index)%passed)) then
          write (*,*) 'Missing process counters during matrix-element checks:',&
               group_index
          stop 1
       endif
       if (size(pgl(group_index)%passed).lt.pgl(group_index)%nproc) then
          write (*,*) 'Process counter array is too short during matrix-element checks:',&
               group_index,size(pgl(group_index)%passed),pgl(group_index)%nproc
          stop 1
       endif
       if (any(pgl(group_index)%passed(1:pgl(group_index)%nproc).lt.me_points)) return
    enddo
    matrix_element_checks_complete=.true.
  end function matrix_element_checks_complete

  subroutine print_contribution_results()
    implicit none
    real(kind=8),allocatable :: aux_res(:,:),aux_unc(:,:)
    real(kind=8) :: component(n_contribution_components),component_unc(n_contribution_components),finite,finite_unc
    real(kind=8) :: closure_scale,closure_residual,closure_tolerance
    character(len=32),parameter :: labels(7)=[character(len=32) ::&
         'Born','Real - local dipoles','Integrated I coefficient -2',&
         'Integrated I coefficient -1','Integrated I finite',&
         'Integrated P','Integrated K']
    integer :: i,j
    call simple_integrator%get_channel_aux_results(aux_res,aux_unc)
    if (ngroups.lt.1 .or. size(aux_res,1).ne.n_contribution_components .or. &
         size(aux_unc,1).ne.n_contribution_components .or. &
         size(aux_res,2).ne.ngroups .or. size(aux_unc,2).ne.ngroups) then
       write(*,*) 'ERROR: contribution results have incompatible dimensions',&
            shape(aux_res),shape(aux_unc),ngroups
       stop 1
    endif
    if (.not.all(ieee_is_finite(aux_res)) .or. &
         .not.all(ieee_is_finite(aux_unc))) then
       write(*,*) 'ERROR: contribution results contain non-finite values'
       stop 1
    endif
    if (any(abs(aux_res).gt.integration_value_limit) .or. &
         any(aux_unc.lt.0d0) .or. any(aux_unc.gt.integration_value_limit) .or. &
         any((aux_unc.ne.0d0) .and. aux_unc.lt.integration_value_floor)) then
       write(*,*) 'ERROR: contribution results exceed the supported numerical range'
       stop 1
    endif
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
    closure_scale=max(1d0,abs(component(2)),abs(component(8)),abs(component(9)))
    closure_residual=abs(component(2)/closure_scale-component(8)/closure_scale-&
         component(9)/closure_scale)
    closure_tolerance=5d-12*max(1d0/closure_scale,abs(component(2))/closure_scale,&
         abs(component(8))/closure_scale+abs(component(9))/closure_scale)
    if (closure_residual.gt.closure_tolerance) then
       write(*,*) 'ERROR: regular and migration strata do not close:',component(2),component(8),component(9)
       stop 1
    endif
    finite=sum(component((/1,2,5,6,7/)))
    finite_unc=sqrt(sum(component_unc((/1,2,5,6,7/))**2))
    write (*,'(a,2x,e14.7,1x,a,1x,e12.5)') 'Finite B+(R-D)+I0+P+K:',finite,'+/-',finite_unc
    write (99,'(a,2x,e14.7,1x,a,1x,e12.5)') 'Finite B+(R-D)+I0+P+K:',finite,'+/-',finite_unc
    deallocate(aux_res,aux_unc)
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
    integer :: iproc,a,target_label,eval_iint,integration_role,icopy,idip,requested_stratum,point_stratum
    integer :: resolution_info,resolution_dipole,kernel_info,kernel_dipole,integrated_info,scale_info,pdf_info,alpha_info
    integer :: amplitude_info
    real(kind=8),dimension(0:3,pgl(ichan)%next) :: p_generated
    real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real(kind=8) :: amp2_integrand,amp2_dip,icoeff(-2:0),pterm,kterm,z,real_hist_factor,real_counter_scale
    real(kind=8) :: coupling_factor,coupling_base
    integer :: coupling_power
    real(kind=8) :: full_f,full_f_abs,regular_f,migration_f,real_margin,threshold_distance
    real(kind=8),allocatable :: dipole_values(:),real_hist_weights(:)
    real(kind=8),allocatable :: mapped_margins(:)
    real(kind=8),allocatable :: icoeff_copy(:,:),pterm_copy(:),kterm_copy(:)
    integer :: mapped_process(pgl(ichan)%next-1)
    real(kind=8) :: unresolved_invariant,resolution_tolerance
    real(kind=8),allocatable :: hard_copy(:)
    logical :: done,time_physics,real_pass,unsplit_real,stratum_selected,analysis_weights_valid
    logical,allocatable :: alpha_active_flags(:),mapped_pass_flags(:)
    f=0d0
    f_abs=0d0
    f_components=0d0
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
    val_abs=0d0
    if (.not.all(numerical_value_is_safe(x)) .or. &
         .not.numerical_value_is_safe(vol)) then
       call report_numerical_rejection('integration input',ichan,eval_iint)
       return
    endif
    if (any(x.lt.0d0) .or. any(x.gt.1d0)) then
       call report_numerical_rejection('integration input',ichan,eval_iint)
       return
    endif
    if (keep_processes_separate) then
       iproc=eval_iint
    else
       iproc=1
    endif
    point_stratum=real_stratum_regular
    stratum_selected=.true.
    threshold_distance=-1d0
    if (tail_tracking_enabled .and. pgl(ichan)%is_subtracted_real) then
       tail_npoints(ichan,iint)=checked_add8(tail_npoints(ichan,iint),1_8,&
            'tail-diagnostic point counter')
       tail_iteration_npoints(ichan,iint)=checked_add8(&
            tail_iteration_npoints(ichan,iint),1_8,&
            'tail-diagnostic iteration point counter')
    endif

    ! Generate phase-space point based on the random numbers 'x(1:ndim)'
    if (time_physics) call cpu_time(tBefore)
    if (read_momenta) then
       if (allocated(p_read)) then
          if (size(p_read,1).ne.pgl(ichan)%next .or. size(p_read,2).ne.4) deallocate(p_read)
       endif
       if (.not.allocated(p_read)) allocate(p_read(pgl(ichan)%next,0:3))
       call read_in_momenta(pgl(ichan)%next,ichan,eval_iint,p_read)
       do a=1,pgl(ichan)%next
          pgl(ichan)%ps(1)%p(:,a)=p_read(a,:)
       enddo
       pgl(ichan)%ps(1)%jac=1d0
       pgl(ichan)%ps(1)%xbjrk(1)=&
            (pgl(ichan)%ps(1)%p(0,1)+pgl(ichan)%ps(1)%p(3,1))/sqrts
       pgl(ichan)%ps(1)%xbjrk(2)=&
            (pgl(ichan)%ps(1)%p(0,2)-pgl(ichan)%ps(1)%p(3,2))/sqrts
    else
       pgl(ichan)%ps(1)%x=x(1:size(pgl(ichan)%ps(1)%x))
       call pgl(ichan)%phase_space%generate_momenta(pgl(ichan)%ps(1))
    endif
    if (debug ) then
       write (*,*) pgl(ichan)%ps(1)%jac
       stop 1
    endif
    if (.not.numerical_value_is_safe(pgl(ichan)%ps(1)%jac) .or. &
         .not.all(numerical_value_is_safe(pgl(ichan)%ps(1)%p)) .or. &
         (include_pdf .and. .not.all(numerical_value_is_safe(pgl(ichan)%ps(1)%xbjrk)))) then
       call report_numerical_rejection('phase space',ichan,eval_iint)
       return
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
    if (.not.generated_momenta_are_valid(pgl(ichan)%ps(1)%p,&
         pgl(ichan)%phase_space%masses,pgl(ichan)%ps(1)%xbjrk,include_pdf)) then
       call report_numerical_rejection('physical momenta',ichan,eval_iint)
       return
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
    if (pgl(ichan)%passed(eval_iint).lt.huge(pgl(ichan)%passed(eval_iint))) &
         pgl(ichan)%passed(eval_iint)=pgl(ichan)%passed(eval_iint)+1
    call compute_multichannel_weight(ichan,eval_iint,pgl(ichan)%ps(1),&
         colour_singlet_multichannel_weight,max(1,requested_stratum))
    if (.not.all(numerical_value_is_safe(colour_singlet_multichannel_weight))) then
       call report_numerical_rejection('multichannel weight',ichan,eval_iint)
       return
    endif
    if (time_physics) then
       call cpu_time(tAfter)
       t_PS= t_PS + (tAfter-tBefore)*dble(timing_sample)
       tBefore=tAfter
    endif
    call compute_the_amps(eval_iint,ichan,use_amplitude_library,amplitude_info)
    if (time_physics) then
       call cpu_time(tAfter)
       t_amp=t_amp+(tAfter-tBefore)*dble(timing_sample)
       tBefore=tAfter
    endif
    if (amplitude_info.ne.0) then
       call report_numerical_rejection('matrix-element current',ichan,eval_iint)
       return
    endif
    call square_the_amps(eval_iint,ichan)
    if (time_physics) then
       call cpu_time(tAfter)
       t_mat=t_mat+(tAfter-tBefore)*dble(timing_sample)
    endif
    if (.not.all(numerical_value_is_safe(pgl(ichan)%amp2)) .or. &
         .not.all(numerical_value_is_safe(pgl(ichan)%amp2_hel(1:pgl(ichan)%nhel(eval_iint))))) then
       call report_numerical_rejection('matrix element',ichan,eval_iint)
       return
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
       if (.not.all(numerical_value_is_safe(pgl(ichan)%amp2)) .or. &
            .not.all(numerical_value_is_safe(pgl(ichan)%amp2_hel(1:pgl(ichan)%nhel(eval_iint))))) then
          call report_numerical_rejection('optimised matrix element',ichan,eval_iint)
          return
       endif
    endif

    if (read_momenta) then
        if (pgl(ichan)%passed(eval_iint).le.me_points) &
             call perform_check(eval_iint,ichan)
        if (matrix_element_checks_complete()) then
           read_momenta=.false.
           matrix_element_checks_done=.true.
        endif
    endif

    real_pass=.true.
    if (pgl(ichan)%is_subtracted_real) real_pass=pass_real_subtracted_cuts(pgl(ichan),eval_iint)
    amp2_integrand=pgl(ichan)%amp2(1)
    if (pgl(ichan)%is_subtracted_real) then
       allocate(dipole_values(pgl(ichan)%dpl(eval_iint)%ndip))
       call evaluate_real_dipoles(eval_iint,ichan,amp2_dip,kernel_info,kernel_dipole,dipole_values)
       if (kernel_info.ne.0) then
          if (dipole_status_is_numerical(kernel_info)) then
             call report_subtraction_rejection('local dipole',ichan,eval_iint,&
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
       if (.not.numerical_value_is_safe(amp2_dip) .or. &
            .not.numerical_value_is_safe(amp2_integrand) .or. &
            .not.all(numerical_value_is_safe(dipole_values)) .or. &
            .not.all(ieee_is_finite(mapped_margins)) .or. &
            .not.ieee_is_finite(threshold_distance)) then
          call report_numerical_rejection('local subtraction',ichan,eval_iint)
          return
       endif
    endif

    ! set scales and update alphaS
    if (time_physics) call cpu_time(tBefore)
    scale_info=0
    call set_scale(scale_choice,pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
         pgl(ichan)%iden_processes(:,1,iproc),scale_ren,scale_info)
    if (scale_info.ne.0) then
       call report_numerical_rejection('dynamic scale',ichan,eval_iint)
       return
    endif
    if (.not.numerical_value_is_safe(scale_ren)) then
       call report_numerical_rejection('dynamic scale',ichan,eval_iint)
       return
    endif
    scale_fac=scale_ren
    scale_shower=scale_ren
    alpha_info=0
    if (use_lhapdf) then
       alphas=alphaspdf(scale_ren)
    else
       alphas=alphas_Q(scale_ren,2,alphas_MZ,alpha_info)
    endif
    if (alpha_info.ne.0 .or. .not.numerical_value_is_safe(alphas)) then
       call report_numerical_rejection('running coupling',ichan,eval_iint)
       return
    endif
    if (alphas.le.0d0) then
       call report_numerical_rejection('running coupling',ichan,eval_iint)
       return
    endif
    
    ! MINT weight, phase-space jacobian and GeV -> pb conversion factor
    if (.not.product_is_safe(vol,pgl(ichan)%ps(1)%jac)) then
       call report_numerical_rejection('common weight',ichan,eval_iint)
       return
    endif
    weight=vol*pgl(ichan)%ps(1)%jac
    if (.not.product_is_safe(weight,conv)) then
       call report_numerical_rejection('common weight',ichan,eval_iint)
       return
    endif
    weight=weight*conv

    ! multiply by the strong coupling
    coupling_power=pgl(ichan)%next-2-pgl(ichan)%amps(eval_iint)%n_sing(1)
    if (coupling_power.gt.0) then
       if (.not.product_is_safe(4d0*pi,alphas)) then
          call report_numerical_rejection('strong-coupling factor',ichan,eval_iint)
          return
       endif
       coupling_base=4d0*pi*alphas
       coupling_factor=1d0
       do a=1,coupling_power
          if (.not.product_is_safe(coupling_factor,coupling_base)) then
             call report_numerical_rejection('strong-coupling factor',ichan,eval_iint)
             return
          endif
          coupling_factor=coupling_factor*coupling_base
       enddo
       if (.not.product_is_safe(weight,coupling_factor)) then
          call report_numerical_rejection('strong-coupling weight',ichan,eval_iint)
          return
       endif
       weight=weight*coupling_factor
    endif
    
    ! multiply by the EW coupling
    if (pgl(ichan)%amps(eval_iint)%n_sing(1).ge.1) then
       coupling_base=2d0*4d0*pi*alphaEW
       coupling_factor=1d0
       do a=1,pgl(ichan)%amps(eval_iint)%n_sing(1)
          if (.not.product_is_safe(coupling_factor,coupling_base)) then
             call report_numerical_rejection('electroweak-coupling factor',ichan,eval_iint)
             return
          endif
          coupling_factor=coupling_factor*coupling_base
       enddo
       if (.not.product_is_safe(weight,coupling_factor)) then
          call report_numerical_rejection('electroweak-coupling weight',ichan,eval_iint)
          return
       endif
       weight=weight*coupling_factor
    endif
    if (.not.numerical_value_is_safe(weight)) then
       call report_numerical_rejection('scale or common weight',ichan,eval_iint)
       return
    endif

    if (keep_processes_separate) then
       if (pgl(ichan)%is_subtracted_real) then
          ! R and every local counterevent have the same phase-space,
          ! coupling, PDF, identity, and multichannel prefactor.  Evaluate
          ! that factor once instead of evolving the PDFs for every dipole.
          ! This path is also used by deterministic tail replay.
          allocate(real_hist_weights(pgl(ichan)%iden_iproc(eval_iint)))
          pdf_info=0
          call physical_matrix_weight(ichan,eval_iint,1d0,weight,&
               colour_singlet_multichannel_weight(eval_iint),real_hist_weights,pdf_info)
          if (pdf_info.ne.0) then
             call report_numerical_rejection('real-emission PDF',ichan,eval_iint)
             return
          endif
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
          if (.not.product_is_safe(amp2_integrand,weight)) then
             call report_numerical_rejection('real matrix weight',ichan,eval_iint)
             return
          endif
          val(1)=amp2_integrand*weight/dble(pgl(ichan)%iden(eval_iint))
          if (.not.numerical_value_is_safe(val(1)) .or. &
               .not.product_is_safe(val(1),colour_singlet_multichannel_weight(eval_iint))) then
             call report_numerical_rejection('real multichannel weight',ichan,eval_iint)
             return
          endif
          val(1)=val(1)*colour_singlet_multichannel_weight(eval_iint)
          if (has_real_process .and. integration_role.ne.1) then
             ! P/K needs the PDF-free reduced Born copies.  integrated_beam
             ! performs the three required table evolutions itself, so do
             ! not evolve the two Born PDFs here only to discard them.
             allocate(hard_copy(pgl(ichan)%iden_iproc(eval_iint)))
             hard_copy=val(1)*pgl(ichan)%idenCOandMAPfactor(&
                  1:pgl(ichan)%iden_iproc(eval_iint),eval_iint)
             f=0d0
             f_abs=0d0
          else
             pdf_info=0
             call include_PDF_and_identical_procs(val,val_abs,pgl(ichan),eval_iint,pdf_info)
             if (pdf_info.ne.0) then
                call report_numerical_rejection('Born PDF',ichan,eval_iint)
                return
             endif
             f_abs=sum(val_abs(1:1))
             f=sum(val(1:1))
          endif
       endif
    else
       if (.not.all(product_is_safe(pgl(ichan)%amp2(1:pgl(ichan)%nproc),weight))) then
          call report_numerical_rejection('Born matrix weight',ichan,eval_iint)
          return
       endif
       val(1:pgl(ichan)%nproc)=pgl(ichan)%amp2(1:pgl(ichan)%nproc)*weight/&
            dble(pgl(ichan)%iden(1:pgl(ichan)%nproc))
       if (.not.all(numerical_value_is_safe(val(1:pgl(ichan)%nproc))) .or. &
            .not.all(product_is_safe(val(1:pgl(ichan)%nproc),&
            colour_singlet_multichannel_weight(1:pgl(ichan)%nproc)))) then
          call report_numerical_rejection('Born multichannel weight',ichan,eval_iint)
          return
       endif
       val(1:pgl(ichan)%nproc)=val(1:pgl(ichan)%nproc)*colour_singlet_multichannel_weight(1:pgl(ichan)%nproc)
       pdf_info=0
       call include_PDF_and_identical_procs(val,val_abs,pgl(ichan),-1,pdf_info)
       if (pdf_info.ne.0) then
          call report_numerical_rejection('combined PDF',ichan,eval_iint)
          return
       endif
       f_abs=sum(val_abs(1:pgl(ichan)%nproc))
       f=sum(val(1:pgl(ichan)%nproc))
    endif
    if (.not.numerical_value_is_safe(f) .or. .not.numerical_value_is_safe(f_abs)) then
       call report_numerical_rejection('PDF-weighted integrand',ichan,eval_iint)
       f=0d0
       f_abs=0d0
       return
    endif
    if (f_abs.lt.0d0) then
       call report_numerical_rejection('PDF-weighted integrand',ichan,eval_iint)
       f=0d0
       f_abs=0d0
       return
    endif
    if (allocated(hard_copy)) then
       if (.not.all(numerical_value_is_safe(hard_copy))) then
          call report_numerical_rejection('beam hard weight',ichan,eval_iint)
          f=0d0
          f_abs=0d0
          return
       endif
    endif
    if (allocated(real_hist_weights)) then
       if (.not.all(numerical_value_is_safe(real_hist_weights)) .or. &
            .not.numerical_value_is_safe(real_hist_factor) .or. &
            .not.numerical_value_is_safe(real_counter_scale)) then
          call report_numerical_rejection('real-emission physical weight',ichan,eval_iint)
          f=0d0
          f_abs=0d0
          return
       endif
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
       analysis_weights_valid=all(numerical_value_is_safe(f_components)) .and. &
            numerical_value_is_safe(f) .and. numerical_value_is_safe(f_abs)
       if (histogram_active .and. stratum_selected) then
          if (real_pass) then
             if (analysis_distinguishes_massless_qcd_flavours) then
                analysis_weights_valid=analysis_weights_valid .and. &
                     all(product_is_safe(pgl(ichan)%amp2(1),real_hist_weights))
             else
                analysis_weights_valid=analysis_weights_valid .and. &
                     product_is_safe(pgl(ichan)%amp2(1),real_hist_factor)
             endif
          endif
          do idip=1,size(dipole_values)
             if (.not.pgl(ichan)%dpl(eval_iint)%dl(idip)%active) cycle
             if (analysis_distinguishes_massless_qcd_flavours) then
                analysis_weights_valid=analysis_weights_valid .and. &
                     all(product_is_safe(dipole_values(idip),real_hist_weights))
             else
                analysis_weights_valid=analysis_weights_valid .and. &
                     product_is_safe(dipole_values(idip),real_hist_factor)
             endif
          enddo
       endif
       if (tail_tracking_enabled .and. stratum_selected) then
          if (real_pass) analysis_weights_valid=analysis_weights_valid .and. &
               product_is_safe(pgl(ichan)%amp2(1),real_counter_scale)
          analysis_weights_valid=analysis_weights_valid .and. &
               all(product_is_safe(dipole_values,real_counter_scale))
       endif
       if (.not.analysis_weights_valid) then
          call report_numerical_rejection('real-emission analysis weight',ichan,eval_iint)
          f=0d0
          f_abs=0d0
          f_components=0d0
          return
       endif
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
          integrated_info=0
          if (keep_processes_separate) then
             if (histogram_active .and. analysis_distinguishes_massless_qcd_flavours) then
                allocate(icoeff_copy(-2:0,pgl(ichan)%iden_iproc(eval_iint)))
                call integrated_endpoint(ichan,eval_iint,&
                     pgl(ichan)%val_procs(1:pgl(ichan)%iden_iproc(eval_iint),eval_iint),&
                     pgl(ichan)%ps(1)%p,scale_ren,alphas,icoeff,icoeff_copy,integrated_info)
             else
                call integrated_endpoint(ichan,eval_iint,&
                     pgl(ichan)%val_procs(1:pgl(ichan)%iden_iproc(eval_iint),eval_iint),&
                     pgl(ichan)%ps(1)%p,scale_ren,alphas,icoeff,status=integrated_info)
             endif
          else
             icoeff=0d0
             do iproc=1,pgl(ichan)%nproc
                call add_endpoint_for_process(ichan,iproc,scale_ren,alphas,icoeff,integrated_info)
                if (integrated_info.ne.0) exit
             enddo
          endif
          if (integrated_info.ne.0) then
             call report_integrated_rejection('endpoint',ichan,eval_iint,integrated_info)
             f=0d0
             f_abs=0d0
             f_components=0d0
             return
          endif
          f_components(3:5)=icoeff(-2:0)
          f=f+icoeff(0)
          f_abs=abs(f_components(1))+abs(icoeff(0))
          analysis_weights_valid=numerical_value_is_safe(f) .and. &
               numerical_value_is_safe(f_abs) .and. all(numerical_value_is_safe(f_components))
          if (allocated(icoeff_copy)) then
             analysis_weights_valid=analysis_weights_valid .and. &
                  all(numerical_value_is_safe(icoeff_copy)) .and. &
                  all(numerical_value_is_safe(&
                  pgl(ichan)%val_procs(1:size(icoeff_copy,2),eval_iint)+icoeff_copy(0,:)))
          endif
          if (.not.analysis_weights_valid) then
             call report_numerical_rejection('integrated endpoint total',ichan,eval_iint)
             f=0d0
             f_abs=0d0
             f_components=0d0
             return
          endif
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
          integrated_info=0
          if (keep_processes_separate) then
             if (histogram_active .and. analysis_distinguishes_massless_qcd_flavours) then
                allocate(pterm_copy(size(hard_copy)),kterm_copy(size(hard_copy)))
                call integrated_beam(ichan,eval_iint,integration_role-1,z,hard_copy,&
                     pgl(ichan)%ps(1)%xbjrk,scale_ren,scale_fac,alphas,&
                     pterm,kterm,pterm_copy,kterm_copy,integrated_info)
             else
                call integrated_beam(ichan,eval_iint,integration_role-1,z,hard_copy,&
                     pgl(ichan)%ps(1)%xbjrk,scale_ren,scale_fac,alphas,pterm,kterm,&
                     status=integrated_info)
             endif
          endif
          if (integrated_info.ne.0) then
             call report_integrated_rejection('beam',ichan,eval_iint,integrated_info)
             f=0d0
             f_abs=0d0
             f_components=0d0
             return
          endif
          f_components(6)=pterm
          f_components(7)=kterm
          f=pterm+kterm
          f_abs=abs(pterm)+abs(kterm)
          analysis_weights_valid=numerical_value_is_safe(f) .and. &
               numerical_value_is_safe(f_abs) .and. all(numerical_value_is_safe(f_components))
          if (allocated(pterm_copy)) then
             analysis_weights_valid=analysis_weights_valid .and. &
                  all(numerical_value_is_safe(pterm_copy)) .and. &
                  all(numerical_value_is_safe(kterm_copy)) .and. &
                  all(numerical_value_is_safe(pterm_copy+kterm_copy))
          endif
          if (.not.analysis_weights_valid) then
             call report_numerical_rejection('integrated beam total',ichan,eval_iint)
             f=0d0
             f_abs=0d0
             f_components=0d0
             return
          endif
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

  subroutine physical_matrix_weight(ichan,iint,matrix_element,common_weight,channel_weight,result,status)
    integer,intent(in) :: ichan,iint
    real(kind=8),intent(in) :: matrix_element,common_weight,channel_weight
    real(kind=8),intent(out) :: result(:)
    integer,intent(out) :: status
    real(kind=8) :: value(1),value_abs(1)
    result=0d0
    status=0
    if (size(result).ne.pgl(ichan)%iden_iproc(iint)) then
       write(*,*) 'ERROR: physical matrix histogram weight array has incompatible size'
       stop 1
    endif
    if (.not.product_is_safe(matrix_element,common_weight)) then
       status=-20
       return
    endif
    value(1)=matrix_element*common_weight/dble(pgl(ichan)%iden(iint))
    if (.not.numerical_value_is_safe(value(1)) .or. &
         .not.product_is_safe(value(1),channel_weight)) then
       value=0d0
       status=-20
       return
    endif
    value(1)=value(1)*channel_weight
    call include_PDF_and_identical_procs(value,value_abs,pgl(ichan),iint,status)
    if (status.ne.0) return
    result=pgl(ichan)%val_procs(1:pgl(ichan)%iden_iproc(iint),iint)
  end subroutine physical_matrix_weight

  subroutine initialize_tail_diagnostics()
    implicit none
    integer :: max_integrals,unit,diagnostic_io_status,allocation_status
    integer(kind=8) :: leaf_count,bytes_per_leaf,workspace_bytes
    character(len=256) :: diagnostic_io_message,allocation_message
    if (ngroups.lt.1 .or. .not.allocated(nintegrals)) then
       write(*,*) 'ERROR: cannot initialise tail diagnostics before integration metadata'
       stop 1
    endif
    if (lbound(nintegrals,1).ne.1 .or. size(nintegrals).ne.ngroups .or. &
         any(nintegrals.lt.1)) then
       write(*,*) 'ERROR: invalid integral counts for tail diagnostics',nintegrals
       stop 1
    endif
    max_integrals=maxval(nintegrals)
    leaf_count=checked_multiply8(int(ngroups,kind=8),int(max_integrals,kind=8),&
         'tail-diagnostic leaf count')
    bytes_per_leaf=(2_8*int(n_tail_records,kind=8)*&
         int(storage_size(tail_record()),kind=8)+&
         2_8*int(storage_size(0_8),kind=8)+&
         4_8*int(storage_size(0d0),kind=8)+7_8)/8_8
    workspace_bytes=checked_multiply8(leaf_count,bytes_per_leaf,&
         'tail-diagnostic workspace')
    if (workspace_bytes.gt.max_process_workspace_bytes) then
       write(*,*) 'ERROR: tail-diagnostic workspace exceeds the supported limit',&
            workspace_bytes,max_process_workspace_bytes
       stop 1
    endif
    if (allocated(tail_residual_records)) deallocate(tail_residual_records)
    if (allocated(tail_component_records)) deallocate(tail_component_records)
    if (allocated(tail_npoints)) deallocate(tail_npoints)
    if (allocated(tail_iteration_npoints)) deallocate(tail_iteration_npoints)
    if (allocated(tail_residual_sum2)) deallocate(tail_residual_sum2)
    if (allocated(tail_component_sum2)) deallocate(tail_component_sum2)
    if (allocated(tail_iteration_residual_sum2)) deallocate(tail_iteration_residual_sum2)
    if (allocated(tail_iteration_max_residual2)) deallocate(tail_iteration_max_residual2)
    allocation_message=''
    allocate(tail_residual_records(n_tail_records,ngroups,max_integrals),&
         tail_component_records(n_tail_records,ngroups,max_integrals),&
         tail_npoints(ngroups,max_integrals),&
         tail_iteration_npoints(ngroups,max_integrals),&
         tail_residual_sum2(ngroups,max_integrals),&
         tail_component_sum2(ngroups,max_integrals),&
         tail_iteration_residual_sum2(ngroups,max_integrals),&
         tail_iteration_max_residual2(ngroups,max_integrals),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write(*,*) 'ERROR: could not allocate tail diagnostics: ',trim(allocation_message)
       stop 1
    endif
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
       diagnostic_io_message=''
       open(newunit=unit,file=trim(tail_replay_output),status='replace',action='write',&
            iostat=diagnostic_io_status,iomsg=diagnostic_io_message)
       call require_diagnostic_io(diagnostic_io_status,'creating stale-replay placeholder',&
            tail_replay_output,diagnostic_io_message)
       close(unit,status='delete',iostat=diagnostic_io_status,iomsg=diagnostic_io_message)
       call require_diagnostic_io(diagnostic_io_status,'deleting stale replay',tail_replay_output,&
            diagnostic_io_message)
       diagnostic_io_message=''
       open(newunit=unit,file=trim(tail_residual_replay_output),status='replace',action='write',&
            iostat=diagnostic_io_status,iomsg=diagnostic_io_message)
       call require_diagnostic_io(diagnostic_io_status,'creating stale-replay placeholder',&
            tail_residual_replay_output,diagnostic_io_message)
       close(unit,status='delete',iostat=diagnostic_io_status,iomsg=diagnostic_io_message)
       call require_diagnostic_io(diagnostic_io_status,'deleting stale replay',&
            tail_residual_replay_output,diagnostic_io_message)
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
    integer :: residual_slot,component_slot,idip
    real(kind=8) :: score,residual_score,component_score,component_candidate
    real(kind=8) :: residual_square,component_square

    if (.not.all(numerical_value_is_safe(x)) .or. &
         .not.numerical_value_is_safe(vol) .or. &
         .not.numerical_value_is_safe(f_point) .or. &
         .not.numerical_value_is_safe(f_abs_point) .or. &
         .not.numerical_value_is_safe(matrix_real) .or. &
         .not.numerical_value_is_safe(dipole_sum) .or. &
         .not.numerical_value_is_safe(matrix_residual) .or. &
         .not.numerical_value_is_safe(common_weight) .or. &
         .not.numerical_value_is_safe(physical_factor) .or. &
         .not.numerical_value_is_safe(counterevent_scale) .or. &
         .not.all(numerical_value_is_safe(dipole_values)) .or. &
         .not.numerical_value_is_safe(threshold_distance)) then
       write(*,*) 'ERROR: invalid input to subtracted-real tail diagnostics',ichan,iint
       stop 1
    endif
    if (f_abs_point.lt.0d0 .or. counterevent_scale.lt.0d0) then
       write(*,*) 'ERROR: negative magnitude in subtracted-real tail diagnostics',ichan,iint
       stop 1
    endif
    residual_score=abs(f_point)
    component_score=0d0
    if (stratum_selected) then
       if (real_pass) then
          if (.not.product_is_safe(matrix_real,counterevent_scale)) then
             write(*,*) 'ERROR: unsafe real component in tail diagnostics',ichan,iint
             stop 1
          endif
          component_score=abs(matrix_real)*counterevent_scale
       endif
       do idip=1,size(dipole_values)
          if (.not.product_is_safe(dipole_values(idip),counterevent_scale)) then
             write(*,*) 'ERROR: unsafe dipole component in tail diagnostics',ichan,iint,idip
             stop 1
          endif
          component_candidate=abs(dipole_values(idip))*counterevent_scale
          component_score=max(component_score,component_candidate)
       enddo
    endif
    score=max(residual_score,component_score)
    if (.not.ieee_is_finite(score)) then
       write(*,*) 'ERROR: non-finite subtracted-real tail score:',ichan,iint,score
       write(99,*) 'ERROR: non-finite subtracted-real tail score:',ichan,iint,score
       stop 1
    endif

    call checked_diagnostic_square(f_point,residual_square,'tail residual square')
    call checked_diagnostic_square(component_score,component_square,'tail component square')
    call checked_diagnostic_accumulate(tail_residual_sum2(ichan,iint),&
         residual_square,'tail residual sum')
    call checked_diagnostic_accumulate(tail_component_sum2(ichan,iint),&
         component_square,'tail component sum')
    call checked_diagnostic_accumulate(tail_iteration_residual_sum2(ichan,iint),&
         residual_square,'tail iteration residual sum')
    if (.not.ieee_is_finite(tail_iteration_max_residual2(ichan,iint))) then
       write(*,*) 'ERROR: invalid tail maximum accumulator',ichan,iint
       stop 1
    endif
    if (tail_iteration_max_residual2(ichan,iint).lt.0d0) then
       write(*,*) 'ERROR: invalid tail maximum accumulator',ichan,iint
       stop 1
    endif
    tail_iteration_max_residual2(ichan,iint)=max(&
         tail_iteration_max_residual2(ichan,iint),residual_square)
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
    real(kind=8) :: normalization,proxy_contribution

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
          if (.not.ieee_is_finite(tail_iteration_residual_sum2(ichan_local,iint_local)) .or. &
               .not.ieee_is_finite(tail_iteration_max_residual2(ichan_local,iint_local))) then
             write(*,*) 'ERROR: invalid migration-tail accumulator',ichan_local,iint_local
             stop 1
          endif
          if (tail_iteration_residual_sum2(ichan_local,iint_local).lt.0d0 .or. &
               tail_iteration_max_residual2(ichan_local,iint_local).lt.0d0) then
             write(*,*) 'ERROR: negative migration-tail accumulator',ichan_local,iint_local
             stop 1
          endif
          proxy_contribution=tail_iteration_residual_sum2(ichan_local,iint_local)/normalization
          call checked_diagnostic_accumulate(total_proxy,proxy_contribution,&
               'migration-tail variance proxy')
          proxy_contribution=tail_iteration_max_residual2(ichan_local,iint_local)/normalization
          max_proxy=max(max_proxy,proxy_contribution)
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
    integer :: idip,ndip,allocation_status
    character(len=256) :: allocation_message

    record%valid=.true.
    record%real_pass=real_pass
    record%ichan=ichan
    record%iint=iint
    record%physical_iint=physical_iint
    record%stratum=stratum
    record%iteration=checked_add(tail_iteration,1,'tail-record iteration')
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
    ndip=size(dipole_values)
    if (allocated(record%x)) deallocate(record%x)
    if (allocated(record%momenta)) deallocate(record%momenta)
    if (allocated(record%process)) deallocate(record%process)
    if (allocated(record%dipole_values)) deallocate(record%dipole_values)
    if (allocated(record%alpha_variables)) deallocate(record%alpha_variables)
    if (allocated(record%dipole_ijk)) deallocate(record%dipole_ijk)
    if (allocated(record%dipole_topology)) deallocate(record%dipole_topology)
    if (allocated(record%alpha_active)) deallocate(record%alpha_active)
    if (allocated(record%active)) deallocate(record%active)
    if (allocated(record%passes_cuts)) deallocate(record%passes_cuts)
    allocation_message=''
    allocate(record%x(lbound(x,1):ubound(x,1)),&
         record%momenta(lbound(pgl(ichan)%ps(1)%p,1):ubound(pgl(ichan)%ps(1)%p,1),&
         lbound(pgl(ichan)%ps(1)%p,2):ubound(pgl(ichan)%ps(1)%p,2)),&
         record%process(lbound(pgl(ichan)%processes,1):ubound(pgl(ichan)%processes,1)),&
         record%dipole_values(ndip),record%alpha_variables(ndip),&
         record%dipole_ijk(3,ndip),record%dipole_topology(ndip),&
         record%alpha_active(ndip),record%active(ndip),record%passes_cuts(ndip),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write(*,*) 'ERROR: could not allocate a tail-diagnostic record: ',&
            trim(allocation_message)
       stop 1
    endif
    record%x=x
    record%momenta=pgl(ichan)%ps(1)%p
    record%process=pgl(ichan)%processes(:,physical_iint)
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
    integer :: unit,replay_io_status,replay_iteration_out
    character(len=256) :: replay_io_message
    replay_iteration_out=checked_add(tail_iteration,1,'tail-replay iteration')
    replay_io_message=''
    open(newunit=unit,file=trim(output_file),status='replace',action='write',&
         iostat=replay_io_status,iomsg=replay_io_message)
    call require_diagnostic_io(replay_io_status,'opening replay output',output_file,replay_io_message)
    write(unit,'(a)',iostat=replay_io_status,iomsg=replay_io_message) '# AmpliCol tail replay v3'
    call require_diagnostic_io(replay_io_status,'writing replay header',output_file,replay_io_message)
    write(unit,*,iostat=replay_io_status,iomsg=replay_io_message) &
         ichan,iint,size(x),PS_choice,replay_iteration_out,timing_point
    call require_diagnostic_io(replay_io_status,'writing replay metadata',output_file,replay_io_message)
    write(unit,'(*(es25.16,1x))',iostat=replay_io_status,iomsg=replay_io_message) &
         vol,f_point,f_abs_point,score,physical_factor,counterevent_scale
    call require_diagnostic_io(replay_io_status,'writing replay weights',output_file,replay_io_message)
    write(unit,'(*(es25.16,1x))',iostat=replay_io_status,iomsg=replay_io_message) alpha_dipole
    call require_diagnostic_io(replay_io_status,'writing replay alpha values',output_file,replay_io_message)
    write(unit,'(*(es25.16,1x))',iostat=replay_io_status,iomsg=replay_io_message) x
    call require_diagnostic_io(replay_io_status,'writing replay coordinates',output_file,replay_io_message)
    close(unit,iostat=replay_io_status,iomsg=replay_io_message)
    call require_diagnostic_io(replay_io_status,'closing replay output',output_file,replay_io_message)
  end subroutine write_tail_replay_point

  subroutine replay_saved_tail_point()
    implicit none
    character(len=256) :: header
    character(len=256) :: io_message
    integer :: unit,replay_ichan,replay_iint,nrandom,replay_ps,replay_iteration,io_status
    integer :: replay_physical_iint,replay_stratum
    integer(kind=8) :: replay_point
    real(kind=8) :: replay_vol,expected_f,expected_f_abs,expected_score,replay_alpha(4),tolerance
    real(kind=8) :: comparison_scale,comparison_residual
    real(kind=8) :: replay_physical_factor,replay_counter_scale
    real(kind=8),allocatable :: replay_x(:)
    logical :: replay_v2
    character(len=9),parameter :: stratum_name(2)=[character(len=9) :: 'regular','migration']

    io_message=''
    open(newunit=unit,file=trim(tail_replay_file),status='old',action='read',&
         iostat=io_status,iomsg=io_message)
    call require_tail_replay_io(io_status,'opening',io_message)
    read(unit,'(a)',iostat=io_status,iomsg=io_message) header
    call require_tail_replay_io(io_status,'reading the header from',io_message)
    replay_v2=trim(header).eq.'# AmpliCol tail replay v2'
    if (.not.replay_v2 .and. trim(header).ne.'# AmpliCol tail replay v3') then
       write(*,*) 'ERROR: unsupported tail replay file: ',trim(header)
       stop 1
    endif
    read(unit,*,iostat=io_status,iomsg=io_message) &
         replay_ichan,replay_iint,nrandom,replay_ps,replay_iteration,replay_point
    call require_tail_replay_io(io_status,'reading metadata from',io_message)
    read(unit,*,iostat=io_status,iomsg=io_message) &
         replay_vol,expected_f,expected_f_abs,expected_score,replay_physical_factor,replay_counter_scale
    call require_tail_replay_io(io_status,'reading weights from',io_message)
    read(unit,*,iostat=io_status,iomsg=io_message) replay_alpha
    call require_tail_replay_io(io_status,'reading alpha values from',io_message)

    if (replay_ps.ne.PS_choice) then
       write(*,*) 'ERROR: tail replay phase-space choice differs:',replay_ps,PS_choice
       stop 1
    endif
    if (.not.all(numerical_value_is_safe(replay_alpha))) then
       write(*,*) 'ERROR: tail replay contains invalid alpha values:',replay_alpha
       stop 1
    endif
    comparison_scale=max(1d0,maxval(abs(replay_alpha)),maxval(abs(alpha_dipole)))
    comparison_residual=maxval(abs(replay_alpha/comparison_scale-&
         alpha_dipole/comparison_scale))
    if (comparison_residual.gt.1d-14/comparison_scale) then
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
    if (nrandom.lt.1) then
       write(*,*) 'ERROR: tail replay random-coordinate count is invalid:',nrandom
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
    if (replay_iteration.lt.1 .or. replay_point.lt.1_8) then
       write(*,*) 'ERROR: tail replay iteration/point metadata is invalid:',&
            replay_iteration,replay_point
       stop 1
    endif
    if (.not.numerical_value_is_safe(replay_vol) .or. &
         .not.numerical_value_is_safe(expected_f) .or. &
         .not.numerical_value_is_safe(expected_f_abs) .or. &
         .not.numerical_value_is_safe(expected_score) .or. &
         .not.numerical_value_is_safe(replay_physical_factor) .or. &
         .not.numerical_value_is_safe(replay_counter_scale)) then
       write(*,*) 'ERROR: tail replay contains invalid saved weights'
       stop 1
    endif
    if (replay_vol.le.0d0 .or. expected_f_abs.lt.0d0 .or. expected_score.lt.0d0 .or. &
         replay_counter_scale.lt.0d0) then
       write(*,*) 'ERROR: tail replay contains invalid saved weights'
       stop 1
    endif
    allocate(replay_x(nrandom))
    read(unit,*,iostat=io_status,iomsg=io_message) replay_x
    call require_tail_replay_io(io_status,'reading coordinates from',io_message)
    close(unit,iostat=io_status,iomsg=io_message)
    call require_tail_replay_io(io_status,'closing',io_message)
    if (.not.all(numerical_value_is_safe(replay_x))) then
       write(*,*) 'ERROR: tail replay contains invalid random coordinates'
       stop 1
    endif
    if (any(replay_x.lt.0d0) .or. any(replay_x.gt.1d0)) then
       write(*,*) 'ERROR: tail replay contains invalid random coordinates'
       stop 1
    endif

    if (replay_v2) then
       call integrand(replay_ichan,replay_iint,replay_x,replay_vol,f(1),f_abs(1),f_aux(:,1),&
            replay_physical_factor,replay_counter_scale,.true.)
    else
       call integrand(replay_ichan,replay_iint,replay_x,replay_vol,f(1),f_abs(1),f_aux(:,1),&
            replay_physical_factor,replay_counter_scale)
    endif
    comparison_scale=max(1d0,abs(expected_f),abs(f(1)))
    tolerance=5d-11*comparison_scale
    comparison_residual=abs(f(1)/comparison_scale-expected_f/comparison_scale)
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
    if (comparison_residual.gt.5d-11) then
       write(*,'(a,es12.4)') 'Tail replay: FAIL, tolerance ',tolerance
       stop 1
    endif
    write(*,'(a)') 'Tail replay: PASS'
    deallocate(replay_x)
  end subroutine replay_saved_tail_point

  subroutine require_tail_replay_io(io_status,operation,io_message)
    implicit none
    integer,intent(in) :: io_status
    character(len=*),intent(in) :: operation,io_message
    if (io_status.ne.0) then
       write(*,*) 'ERROR: failed while ',trim(operation),' tail replay file ',&
            trim(tail_replay_file),': ',trim(io_message)
       stop 1
    endif
  end subroutine require_tail_replay_io

  subroutine require_diagnostic_io(io_status,operation,output_file,io_message)
    implicit none
    integer,intent(in) :: io_status
    character(len=*),intent(in) :: operation,output_file,io_message
    if (io_status.ne.0) then
       write(*,*) 'ERROR: failed while ',trim(operation),' ',trim(output_file),': ',trim(io_message)
       stop 1
    endif
  end subroutine require_diagnostic_io

  subroutine write_tail_diagnostics()
    implicit none
    integer :: unit,ichan_local,iint_local,irecord,nvalid_residual,nvalid_component
    integer :: diagnostic_io_status
    integer :: physical_iint_local,stratum_local
    real(kind=8) :: residual_top2,component_top2,record_square
    character(len=256) :: diagnostic_io_message
    character(len=9),parameter :: stratum_name(2)=[character(len=9) :: 'regular','migration']

    diagnostic_io_message=''
    open(newunit=unit,file=trim(tail_logfile),status='replace',action='write',&
         iostat=diagnostic_io_status,iomsg=diagnostic_io_message)
    call require_diagnostic_io(diagnostic_io_status,'opening tail diagnostics',tail_logfile,&
         diagnostic_io_message)
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
                call checked_diagnostic_square(&
                     tail_residual_records(irecord,ichan_local,iint_local)%f,&
                     record_square,'retained tail residual square')
                call checked_diagnostic_accumulate(residual_top2,record_square,&
                     'retained tail residual sum')
             endif
             if (tail_component_records(irecord,ichan_local,iint_local)%valid) then
                call checked_diagnostic_square(&
                     tail_component_records(irecord,ichan_local,iint_local)%component_score,&
                     record_square,'retained tail component square')
                call checked_diagnostic_accumulate(component_top2,record_square,&
                     'retained tail component sum')
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
    close(unit,iostat=diagnostic_io_status,iomsg=diagnostic_io_message)
    call require_diagnostic_io(diagnostic_io_status,'closing tail diagnostics',tail_logfile,&
         diagnostic_io_message)
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
    real(kind=8) :: scale
    safe_fraction=huge(1d0)
    if (.not.ieee_is_finite(numerator) .or. &
         .not.ieee_is_finite(denominator)) return
    if (numerator.lt.0d0 .or. denominator.lt.0d0) return
    if (denominator.eq.0d0) then
       if (numerator.eq.0d0) safe_fraction=0d0
       return
    endif
    scale=max(numerator,denominator)
    if (numerator/scale.gt.denominator/scale+1d-12) return
    safe_fraction=numerator/denominator
  end function safe_fraction

  subroutine make_identical_reduced_process(base_real,copy_real,base_reduced,copy_reduced)
    integer,intent(in) :: base_real(:),copy_real(:),base_reduced(:)
    integer,intent(out) :: copy_reduced(:)
    integer :: generic,ileg,target_flavour
    if (size(copy_real).ne.size(base_real) .or. &
         size(copy_reduced).ne.size(base_reduced)) then
       write (*,*) 'Incompatible identical-process remapping dimensions'
       stop 1
    endif
    copy_reduced=base_reduced
    do generic=1,2
       target_flavour=0
       do ileg=1,size(base_real)
          if (base_real(ileg).eq.generic .or. base_real(ileg).eq.-generic) then
             if (copy_real(ileg).ge.1 .and. copy_real(ileg).le.6) then
                target_flavour=copy_real(ileg)
             elseif (copy_real(ileg).le.-1 .and. copy_real(ileg).ge.-6) then
                target_flavour=-copy_real(ileg)
             else
                write (*,*) 'Invalid identical-process target flavour:',&
                     ileg,copy_real(ileg)
                stop 1
             endif
             exit
          endif
       enddo
       if (target_flavour.eq.0) cycle
       do ileg=1,size(base_reduced)
          if (base_reduced(ileg).eq.generic .or. &
               base_reduced(ileg).eq.-generic) then
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
       invariant_rejections=checked_add8(invariant_rejections,1_8,&
            'subtraction invariant rejection counter')
       count=invariant_rejections
    else
       kernel_rejections=checked_add8(kernel_rejections,1_8,&
            'subtraction kernel rejection counter')
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

  elemental logical function numerical_value_is_safe(value)
    real(kind=8),intent(in) :: value
    numerical_value_is_safe=.false.
    if (.not.ieee_is_finite(value)) return
    if (abs(value).gt.integration_value_limit) return
    if (value.ne.0d0 .and. abs(value).lt.integration_value_floor) return
    numerical_value_is_safe=.true.
  end function numerical_value_is_safe

  elemental logical function product_is_safe(first,second)
    real(kind=8),intent(in) :: first,second
    product_is_safe=.false.
    if (.not.numerical_value_is_safe(first) .or. .not.numerical_value_is_safe(second)) return
    if (first.eq.0d0 .or. second.eq.0d0) then
       product_is_safe=.true.
    else
       if (abs(first).gt.integration_value_limit/abs(second)) return
       if (abs(second).lt.1d0) then
          if (abs(first).lt.integration_value_floor/abs(second)) return
       elseif (abs(first).lt.1d0) then
          if (abs(second).lt.integration_value_floor/abs(first)) return
       endif
       product_is_safe=.true.
    endif
  end function product_is_safe

  subroutine checked_diagnostic_square(value,square,context)
    real(kind=8),intent(in) :: value
    real(kind=8),intent(out) :: square
    character(len=*),intent(in) :: context
    square=0d0
    if (.not.ieee_is_finite(value)) then
       write(*,*) 'ERROR: non-finite value in ',trim(context),value
       stop 1
    endif
    if (abs(value).gt.sqrt(huge(1d0))) then
       write(*,*) 'ERROR: overflow in ',trim(context),value
       stop 1
    endif
    square=value*value
    if (.not.ieee_is_finite(square)) then
       write(*,*) 'ERROR: non-finite result in ',trim(context),value
       stop 1
    endif
  end subroutine checked_diagnostic_square

  subroutine checked_diagnostic_accumulate(total,increment,context)
    real(kind=8),intent(inout) :: total
    real(kind=8),intent(in) :: increment
    character(len=*),intent(in) :: context
    if (.not.ieee_is_finite(total) .or. .not.ieee_is_finite(increment)) then
       write(*,*) 'ERROR: non-finite accumulator in ',trim(context),total,increment
       stop 1
    endif
    if (total.lt.0d0 .or. increment.lt.0d0) then
       write(*,*) 'ERROR: negative accumulator in ',trim(context),total,increment
       stop 1
    endif
    if (increment.gt.huge(total)-total) then
       write(*,*) 'ERROR: accumulator overflow in ',trim(context),total,increment
       stop 1
    endif
    total=total+increment
  end subroutine checked_diagnostic_accumulate

  subroutine report_numerical_rejection(stage,ichan,iint)
    character(len=*),intent(in) :: stage
    integer,intent(in) :: ichan,iint
    integer(kind=8),save :: rejection_count=0_8

    rejection_count=checked_add8(rejection_count,1_8,&
         'numerical rejection counter')
    if (rejection_count.le.10_8 .or. mod(rejection_count,1000_8).eq.0_8) then
       write(*,'(a,a,a,i0,a,2(i0,1x))') 'INFO: rejected invalid numerical point (',&
            trim(stage),'), count=',rejection_count,' channel/process=',ichan,iint
       write(99,'(a,a,a,i0,a,2(i0,1x))') 'INFO: rejected invalid numerical point (',&
            trim(stage),'), count=',rejection_count,' channel/process=',ichan,iint
    elseif (rejection_count.eq.11_8) then
       write(*,'(a)') 'INFO: further invalid-numerical-point messages suppressed'
       write(99,'(a)') 'INFO: further invalid-numerical-point messages suppressed'
    endif
  end subroutine report_numerical_rejection

  subroutine report_integrated_rejection(stage,ichan,iint,status)
    character(len=*),intent(in) :: stage
    integer,intent(in) :: ichan,iint,status
    integer(kind=8),save :: rejection_count=0_8

    rejection_count=checked_add8(rejection_count,1_8,&
         'integrated rejection counter')
    if (rejection_count.le.10_8 .or. mod(rejection_count,1000_8).eq.0_8) then
       write(*,'(a,a,a,i0,a,3(i0,1x))') 'INFO: rejected numerical integrated ',&
            trim(stage),' point, count=',rejection_count,&
            ' channel/process/status=',ichan,iint,status
       write(99,'(a,a,a,i0,a,3(i0,1x))') 'INFO: rejected numerical integrated ',&
            trim(stage),' point, count=',rejection_count,&
            ' channel/process/status=',ichan,iint,status
    elseif (rejection_count.eq.11_8) then
       write(*,'(a)') 'INFO: further integrated numerical-rejection messages suppressed'
       write(99,'(a)') 'INFO: further integrated numerical-rejection messages suppressed'
    endif
  end subroutine report_integrated_rejection

  subroutine add_endpoint_for_process(igroup,iproc,mu_ren,alpha_s,total,status)
    integer,intent(in) :: igroup,iproc
    real(kind=8),intent(in) :: mu_ren,alpha_s
    real(kind=8),intent(inout) :: total(-2:0)
    integer,intent(out) :: status
    real(kind=8) :: one(-2:0)
    call integrated_endpoint(igroup,iproc,&
         pgl(igroup)%val_procs(1:pgl(igroup)%iden_iproc(iproc),iproc),&
         pgl(igroup)%ps(1)%p,mu_ren,alpha_s,one,status=status)
    if (status.ne.0) return
    total=total+one
  end subroutine add_endpoint_for_process

  subroutine optimise_the_amplitudes(iint,ichan,done)
    implicit none
    integer,intent(in) :: iint,ichan
    logical,intent(out) :: done
    integer,dimension(:),allocatable :: helicity_filter
    real(kind=8),dimension(:),allocatable :: amp2_save,amp2_hel_save,amps_r_save
    complex(kind=8),dimension(:),allocatable :: amps_save
    done=.false.
    if (pgl(ichan)%passed(iint).le.n_amplitude_optimisation_samples) then
       call record_helicity_optimisation_sample(pgl(ichan),iint)
       call pgl(ichan)%amps(iint)%record_optimisation_sample(&
            pgl(ichan)%passed(iint),n_amplitude_optimisation_samples)
    endif
    call find_same_flavour(pgl(ichan),nevent_hel_filter,pgl(ichan)%amp2)
    if (pgl(ichan)%passed(iint).eq.nevent_hel_filter) then
       ! Validate current sharing and helicity filtering independently at
       ! the last warm-up point.  This keeps failures attributable and also
       ! avoids the non-conforming nproc-sized copy used by the old
       ! keep_processes_separate path.
       allocate(amp2_save(size(pgl(ichan)%amp2)))
       amp2_save=pgl(ichan)%amp2
       allocate(amp2_hel_save(pgl(ichan)%nhel(iint)))
       amp2_hel_save=pgl(ichan)%amp2_hel(1:pgl(ichan)%nhel(iint))
       if (use_real_gluons) then
          allocate(amps_r_save(size(pgl(ichan)%amps(iint)%amps_r)))
          amps_r_save=pgl(ichan)%amps(iint)%amps_r
       else
          allocate(amps_save(size(pgl(ichan)%amps(iint)%amps)))
          amps_save=pgl(ichan)%amps(iint)%amps
       endif
       if (use_cross_process_optimisation_of_currents) then
          call pgl(ichan)%amps(iint)%optimise_evaluation(pgl(ichan)%next)
          call compute_the_amps(iint,ichan,use_amplitude_library)
          call square_the_amps(iint,ichan)
          call check_optimised_matrix_element('current sharing',amp2_save,&
               pgl(ichan)%amp2,ichan,iint)
          call check_optimised_matrix_element('current sharing by helicity',&
               amp2_hel_save,pgl(ichan)%amp2_hel(1:size(amp2_hel_save)),ichan,iint)
          if (use_real_gluons) then
             call check_optimised_real_amplitudes('current sharing',amps_r_save,&
                  pgl(ichan)%amps(iint)%amps_r,ichan,iint)
          else
             call check_optimised_complex_amplitudes('current sharing',amps_save,&
                  pgl(ichan)%amps(iint)%amps,ichan,iint)
          endif
       endif
       call setup_helicity_filter(pgl(ichan),iint,helicity_filter)
       call pgl(ichan)%amps(iint)%clear_optimisation_samples()
       call compute_the_amps(iint,ichan,use_amplitude_library)
       call square_the_amps(iint,ichan)
       call check_optimised_matrix_element('helicity filtering',amp2_save,&
            pgl(ichan)%amp2,ichan,iint)
       call check_filtered_helicity_weights(helicity_filter,amp2_hel_save,&
            pgl(ichan)%amp2_hel(1:pgl(ichan)%nhel(iint)),ichan,iint)
       if (use_real_gluons) then
          call check_filtered_real_amplitudes(helicity_filter,amps_r_save,&
               pgl(ichan)%amps(iint)%amps_r,ichan,iint)
       else
          call check_filtered_complex_amplitudes(helicity_filter,amps_save,&
               pgl(ichan)%amps(iint)%amps,ichan,iint)
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

  subroutine record_helicity_optimisation_sample(pgl,iint)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer,intent(in) :: iint
    integer :: isample
    real(kind=8) :: sample_scale
    if (.not.allocated(pgl%amp2_hel_samples)) then
       allocate(pgl%amp2_hel_samples(maxval(pgl%nhel),size(pgl%nhel),&
            n_amplitude_optimisation_samples))
       pgl%amp2_hel_samples=0d0
    endif
    isample=pgl%passed(iint)
    if (.not.all(ieee_is_finite(pgl%amp2_hel(1:pgl%nhel(iint))))) then
       write (*,*) 'Non-finite helicity weight encountered during amplitude optimisation',&
            iint,isample
       stop 1
    endif
    sample_scale=maxval(abs(pgl%amp2_hel(1:pgl%nhel(iint))))
    pgl%amp2_hel_samples(:,iint,isample)=0d0
    if (sample_scale.gt.tiny(1d0)) then
       pgl%amp2_hel_samples(1:pgl%nhel(iint),iint,isample)=&
            pgl%amp2_hel(1:pgl%nhel(iint))/sample_scale
    endif
  end subroutine record_helicity_optimisation_sample

  subroutine check_optimised_matrix_element(stage,reference,value,ichan,iint)
    implicit none
    character(len=*),intent(in) :: stage
    real(kind=8),intent(in) :: reference(:),value(:)
    integer,intent(in) :: ichan,iint
    real(kind=8) :: scale(size(reference)),difference(size(reference)),&
         baseline(size(reference))
    if (size(reference).ne.size(value)) then
       write (*,*) 'Amplitude optimisation changed matrix-element shape during ',trim(stage)
       stop 1
    endif
    if (.not.all(ieee_is_finite(reference)) .or. &
         .not.all(ieee_is_finite(value))) then
       write (*,*) 'Non-finite matrix element encountered during ',trim(stage),ichan,iint
       stop 1
    endif
    scale=max(abs(reference),abs(value),tiny(1d0))
    difference=abs(reference/scale-value/scale)
    baseline=max(abs(reference/scale),abs(value/scale),tiny(1d0)/scale)
    if (any(difference.gt.1d-9*baseline)) then
       write (*,*) 'Amplitude optimisation changed matrix element during ',trim(stage),ichan,iint
       write (*,*) 'reference:',reference
       write (*,*) 'optimised:',value
       write (*,*) 'relative difference:',difference/baseline
       stop 1
    endif
  end subroutine check_optimised_matrix_element

  subroutine check_optimised_complex_amplitudes(stage,reference,value,ichan,iint)
    implicit none
    character(len=*),intent(in) :: stage
    complex(kind=8),intent(in) :: reference(:),value(:)
    integer,intent(in) :: ichan,iint
    real(kind=8) :: scale(size(reference)),difference(size(reference)),&
         baseline(size(reference))
    if (size(reference).ne.size(value)) then
       write (*,*) 'Amplitude optimisation changed amplitude shape during ',trim(stage)
       stop 1
    endif
    if (.not.complex_array_is_finite(reference) .or. &
         .not.complex_array_is_finite(value)) then
       write (*,*) 'Non-finite complex amplitude encountered during ',trim(stage),ichan,iint
       stop 1
    endif
    scale=max(abs(real(reference,kind=8)),abs(aimag(reference)),&
         abs(real(value,kind=8)),abs(aimag(value)),tiny(1d0))
    difference=abs(reference/scale-value/scale)
    baseline=max(abs(reference/scale),abs(value/scale),tiny(1d0)/scale)
    if (any(difference.gt.1d-9*baseline)) then
       write (*,*) 'Amplitude optimisation changed a complex amplitude during ',&
            trim(stage),ichan,iint,maxval(difference/baseline)
       stop 1
    endif
  end subroutine check_optimised_complex_amplitudes

  subroutine check_optimised_real_amplitudes(stage,reference,value,ichan,iint)
    implicit none
    character(len=*),intent(in) :: stage
    real(kind=8),intent(in) :: reference(:),value(:)
    integer,intent(in) :: ichan,iint
    real(kind=8) :: scale(size(reference)),difference(size(reference)),&
         baseline(size(reference))
    if (size(reference).ne.size(value)) then
       write (*,*) 'Amplitude optimisation changed real-amplitude shape during ',trim(stage)
       stop 1
    endif
    if (.not.all(ieee_is_finite(reference)) .or. &
         .not.all(ieee_is_finite(value))) then
       write (*,*) 'Non-finite real amplitude encountered during ',trim(stage),ichan,iint
       stop 1
    endif
    scale=max(abs(reference),abs(value),tiny(1d0))
    difference=abs(reference/scale-value/scale)
    baseline=max(abs(reference/scale),abs(value/scale),tiny(1d0)/scale)
    if (any(difference.gt.1d-9*baseline)) then
       write (*,*) 'Amplitude optimisation changed a real amplitude during ',&
            trim(stage),ichan,iint,maxval(difference/baseline)
       stop 1
    endif
  end subroutine check_optimised_real_amplitudes

  subroutine check_filtered_helicity_weights(filter,reference,value,ichan,iint)
    implicit none
    integer,intent(in) :: filter(:),ichan,iint
    real(kind=8),intent(in) :: reference(:),value(:)
    integer :: old_hel,new_hel,member,nmembers
    real(kind=8) :: expected,reference_scale,scale
    if (size(filter).ne.size(reference) .or. count(filter.gt.0).ne.size(value)) then
       write (*,*) 'Helicity-filter validation has inconsistent dimensions',ichan,iint
       stop 1
    endif
    if (.not.all(ieee_is_finite(reference)) .or. &
         .not.all(ieee_is_finite(value))) then
       write (*,*) 'Non-finite helicity-filter weight encountered',ichan,iint
       stop 1
    endif
    reference_scale=max(maxval(abs(reference)),tiny(1d0))
    do old_hel=1,size(filter)
       if (filter(old_hel).eq.0 .and. &
            abs(reference(old_hel)).gt.100d0*helicity_zero_tolerance*reference_scale) then
          write (*,*) 'A discarded helicity is non-zero at the holdout point',&
               ichan,iint,old_hel,reference(old_hel)/reference_scale
          stop 1
       endif
    enddo
    new_hel=0
    do old_hel=1,size(filter)
       if (filter(old_hel).le.0) cycle
       new_hel=new_hel+1
       nmembers=1
       do member=old_hel+1,size(filter)
          if (filter(member).eq.-old_hel) then
             nmembers=nmembers+1
          endif
       enddo
       if (nmembers.ne.filter(old_hel)) then
          write (*,*) 'Helicity-filter multiplicity is inconsistent',&
               ichan,iint,old_hel,nmembers,filter(old_hel)
          stop 1
       endif
       scale=max(maxval(abs(reference)),abs(value(new_hel)),tiny(1d0))
       expected=reference(old_hel)/scale
       do member=old_hel+1,size(filter)
          if (filter(member).eq.-old_hel) &
               expected=expected+reference(member)/scale
       enddo
       if (abs(expected-value(new_hel)/scale).gt.&
            1d-9*max(abs(expected),abs(value(new_hel)/scale),tiny(1d0)/scale)) then
          write (*,*) 'Helicity filtering changed a grouped helicity weight',&
               ichan,iint,old_hel,expected,value(new_hel)/scale
          stop 1
       endif
    enddo
  end subroutine check_filtered_helicity_weights

  subroutine check_filtered_complex_amplitudes(filter,reference,value,ichan,iint)
    implicit none
    integer,intent(in) :: filter(:),ichan,iint
    complex(kind=8),intent(in) :: reference(:),value(:)
    integer :: old_hel,new_hel
    real(kind=8) :: scale,baseline,difference
    if (size(filter).ne.size(reference) .or. count(filter.gt.0).ne.size(value)) then
       write (*,*) 'Complex helicity-filter validation has inconsistent dimensions',ichan,iint
       stop 1
    endif
    if (.not.complex_array_is_finite(reference) .or. &
         .not.complex_array_is_finite(value)) then
       write (*,*) 'Non-finite filtered complex amplitude encountered',ichan,iint
       stop 1
    endif
    new_hel=0
    do old_hel=1,size(filter)
       if (filter(old_hel).le.0) cycle
       new_hel=new_hel+1
       scale=max(abs(real(reference(old_hel),kind=8)),&
            abs(aimag(reference(old_hel))),abs(real(value(new_hel),kind=8)),&
            abs(aimag(value(new_hel))),tiny(1d0))
       difference=abs(reference(old_hel)/scale-value(new_hel)/scale)
       baseline=max(abs(reference(old_hel)/scale),abs(value(new_hel)/scale),&
            tiny(1d0)/scale)
       if (difference.gt.1d-9*baseline) then
          write (*,*) 'Helicity filtering changed a retained complex amplitude',&
               ichan,iint,old_hel,new_hel
          stop 1
       endif
    enddo
  end subroutine check_filtered_complex_amplitudes

  subroutine check_filtered_real_amplitudes(filter,reference,value,ichan,iint)
    implicit none
    integer,intent(in) :: filter(:),ichan,iint
    real(kind=8),intent(in) :: reference(:),value(:)
    integer :: old_hel,new_hel
    real(kind=8) :: scale,baseline,difference
    if (size(filter).ne.size(reference) .or. count(filter.gt.0).ne.size(value)) then
       write (*,*) 'Real helicity-filter validation has inconsistent dimensions',ichan,iint
       stop 1
    endif
    if (.not.all(ieee_is_finite(reference)) .or. &
         .not.all(ieee_is_finite(value))) then
       write (*,*) 'Non-finite filtered real amplitude encountered',ichan,iint
       stop 1
    endif
    new_hel=0
    do old_hel=1,size(filter)
       if (filter(old_hel).le.0) cycle
       new_hel=new_hel+1
       scale=max(abs(reference(old_hel)),abs(value(new_hel)),tiny(1d0))
       difference=abs(reference(old_hel)/scale-value(new_hel)/scale)
       baseline=max(abs(reference(old_hel)/scale),abs(value(new_hel)/scale),&
            tiny(1d0)/scale)
       if (difference.gt.1d-9*baseline) then
          write (*,*) 'Helicity filtering changed a retained real amplitude',&
               ichan,iint,old_hel,new_hel
          stop 1
       endif
    enddo
  end subroutine check_filtered_real_amplitudes

  logical function complex_array_is_finite(values)
    implicit none
    complex(kind=8),intent(in) :: values(:)
    complex_array_is_finite=all(ieee_is_finite(real(values,kind=8))) .and. &
         all(ieee_is_finite(aimag(values)))
  end function complex_array_is_finite


  subroutine setup_helicity_filter(pgl,iint,original_filter)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer,dimension(:),allocatable,intent(out),optional :: original_filter
    integer :: i,ih1,iint,current_index,target_index,link_value,chain_steps,&
         allocation_status
    character(len=256) :: allocation_message
    if (iint.lt.1 .or. iint.gt.pgl%nproc .or. pgl%nproc.lt.1 .or. &
         .not.allocated(pgl%nhel) .or. .not.allocated(pgl%amps)) then
       write (*,*) 'Invalid process state for helicity filtering',iint,pgl%nproc
       stop 1
    endif
    if (size(pgl%nhel).lt.pgl%nproc .or. size(pgl%amps).lt.iint) then
       write (*,*) 'Incomplete amplitude arrays for helicity filtering',&
            iint,size(pgl%nhel),pgl%nproc,size(pgl%amps)
       stop 1
    endif
    if (any(pgl%nhel(1:pgl%nproc).lt.1)) then
       write (*,*) 'Invalid stored helicity counts for filtering',&
            pgl%nhel(1:pgl%nproc)
       stop 1
    endif
    if (pgl%nhel(iint).ne.pgl%amps(iint)%n_amps) then
       write (*,*) 'Inconsistent amplitude dimensions for helicity filtering',&
            iint,pgl%nhel(iint),pgl%amps(iint)%n_amps
       stop 1
    endif
    if (.not.allocated(pgl%include_hel)) then
       allocate(pgl%include_hel(maxval(pgl%nhel(1:pgl%nproc)),pgl%nproc),&
            stat=allocation_status,errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate helicity filter: ',trim(allocation_message)
          stop 1
       endif
       pgl%include_hel=0
    elseif (size(pgl%include_hel,1).lt.maxval(pgl%nhel(1:pgl%nproc)) .or. &
         size(pgl%include_hel,2).lt.pgl%nproc) then
       write (*,*) 'Stored helicity-filter workspace has incompatible dimensions'
       stop 1
    endif
    if (.not.allocated(pgl%amp2_hel_samples)) then
       write (*,*) 'Missing helicity warm-up samples',iint
       stop 1
    endif
    if (.not.allocated(pgl%hel_fac) .or. &
         .not.allocated(pgl%amps(iint)%same_flavour_sum) .or. &
         .not.allocated(pgl%amps(iint)%same_flavour_sum_operation)) then
       write (*,*) 'Incomplete helicity-filter warm-up state',iint
       stop 1
    endif
    if (size(pgl%amp2_hel_samples,1).lt.pgl%nhel(iint) .or. &
         size(pgl%amp2_hel_samples,2).lt.iint .or. &
         size(pgl%amp2_hel_samples,3).lt.n_amplitude_optimisation_samples .or. &
         size(pgl%hel_fac,1).lt.pgl%nhel(iint) .or. &
         size(pgl%hel_fac,2).lt.iint .or. &
         size(pgl%amps(iint)%same_flavour_sum,1).lt.pgl%nhel(iint) .or. &
         size(pgl%amps(iint)%same_flavour_sum,2).ne.2 .or. &
         size(pgl%amps(iint)%same_flavour_sum_operation,1).lt.pgl%nhel(iint) .or. &
         size(pgl%amps(iint)%same_flavour_sum_operation,2).ne.2) then
       write (*,*) 'Incomplete helicity-filter warm-up state',iint
       stop 1
    endif
    call build_helicity_filter(pgl%amps(iint),&
         pgl%amp2_hel_samples(1:pgl%nhel(iint),iint,1:n_amplitude_optimisation_samples),&
         pgl%include_hel(1:pgl%nhel(iint),iint),keep_processes_separate)
    if (present(original_filter)) then
       allocate(original_filter(pgl%nhel(iint)),stat=allocation_status,&
            errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not preserve the original helicity filter: ',&
               trim(allocation_message)
          stop 1
       endif
       original_filter=pgl%include_hel(1:pgl%nhel(iint),iint)
    endif

    do i=1,2
       do ih1=1,pgl%nhel(iint)
          current_index=pgl%amps(iint)%same_flavour_sum(ih1,i)
          if (current_index.le.0) cycle
          pgl%amps(iint)%same_flavour_sum_operation(ih1,i)=0
          chain_steps=0
          do
             if (current_index.lt.1 .or. current_index.gt.pgl%nhel(iint)) then
                write (*,*) 'Invalid same-flavour helicity index',ih1,i,current_index
                stop 1
             endif
             link_value=pgl%include_hel(current_index,iint)
             if (link_value.ge.0) exit
             if (link_value.lt.-pgl%nhel(iint)) then
                write (*,*) 'Invalid helicity-filter representative',current_index,link_value
                stop 1
             endif
             target_index=-link_value
             chain_steps=chain_steps+1
             if (chain_steps.gt.pgl%nhel(iint)) then
                write (*,*) 'Cyclic helicity-filter representative chain',ih1,i
                stop 1
             endif
             pgl%amps(iint)%same_flavour_sum_operation(ih1,i)= &
                  ieor(pgl%amps(iint)%same_flavour_sum_operation(ih1,i),&
                  find_operation(pgl,iint,current_index,target_index))
             current_index=target_index
          enddo
          pgl%amps(iint)%same_flavour_sum(ih1,i)=current_index
       enddo
    enddo

    call pgl%amps(iint)%filter_helicity(pgl%next,pgl%nhel(iint),&
         pgl%include_hel(1:pgl%nhel(iint),iint)) ! updates nhel and include_hel
    pgl%hel_fac(1:pgl%nhel(iint),iint)=pgl%include_hel(1:pgl%nhel(iint),iint)
  end subroutine setup_helicity_filter

  integer function find_operation(pgl,iint,first_index,second_index)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: first_index,second_index,iint
    complex(kind=8) :: amp1,amp2
    real(kind=8) :: scale,direct_distance,swapped_distance
    if (use_real_gluons) then
       write (*,*) 'Find operation for same flavour sum only for complex amplitudes'
       stop 1
    endif
    if (iint.lt.1 .or. .not.allocated(pgl%amps)) then
       write (*,*) 'Invalid amplitude state while determining a helicity operation'
       stop 1
    endif
    if (iint.gt.size(pgl%amps)) then
       write (*,*) 'Invalid amplitude state while determining a helicity operation'
       stop 1
    endif
    if (.not.allocated(pgl%amps(iint)%amps)) then
       write (*,*) 'Missing amplitudes while determining a helicity operation'
       stop 1
    endif
    if (first_index.lt.1 .or. second_index.lt.1 .or. &
         first_index.gt.size(pgl%amps(iint)%amps) .or. &
         second_index.gt.size(pgl%amps(iint)%amps)) then
       write (*,*) 'Invalid amplitude index while determining a helicity operation',&
            first_index,second_index
       stop 1
    endif
    amp1=pgl%amps(iint)%amps(first_index)
    amp2=pgl%amps(iint)%amps(second_index)
    if (.not.ieee_is_finite(real(amp1,kind=8)) .or. &
         .not.ieee_is_finite(aimag(amp1)) .or. &
         .not.ieee_is_finite(real(amp2,kind=8)) .or. &
         .not.ieee_is_finite(aimag(amp2))) then
       write(*,*) 'Cannot determine a same-flavour operation from non-finite amplitudes'
       stop 1
    endif
    scale=max(1d0,abs(real(amp1,kind=8)),abs(aimag(amp1)),&
         abs(real(amp2,kind=8)),abs(aimag(amp2)))
    direct_distance=abs(abs(real(amp1,kind=8))/scale-&
         abs(real(amp2,kind=8))/scale)
    swapped_distance=abs(abs(real(amp1,kind=8))/scale-abs(aimag(amp2))/scale)
    find_operation=0
    if (direct_distance.lt.swapped_distance) then
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
    character(len=80) :: library,timing_arg
    integer :: timing_sample_arg
    integer(kind=8) iseed
    common /to_seed/iseed
    call parse_argument(filename,real_filename,input_filename,ncalls0,itmax,PS_choice,iseed,library,tag,&
         read_momenta,me_points,&
         limit_test,timing_arg,timing_sample_arg,accuracy,dim_reg_scheme,has_real_process,&
         tail_replay_file,replay_tail,migration_tail_fraction_limit)

    logfile="Outputs/"//trim(adjustl(tag))//"log_file.txt"
    io_message=''
    open(unit=99,file=trim(logfile),status='replace',action='write',&
         iostat=io_status,iomsg=io_message)
    if (io_status.ne.0) then
       write (*,*) 'Could not create run log ',trim(logfile),': ',trim(io_message)
       stop 1
    endif
    if (limit_test) then
       limit_logfile="Outputs/"//trim(adjustl(tag))//"limit_test_failures.log"
       open(unit=100,file=trim(limit_logfile),action='write',status='replace',&
            iostat=io_status,iomsg=io_message)
       if (io_status.ne.0) then
          write (*,*) 'Could not create limit-test log ',trim(limit_logfile),': ',&
               trim(io_message)
          stop 1
       endif
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
    else
       write (*,*) 'library must be none, create or use: ',trim(library)
       stop 1
    endif

    if (limit_test .and. read_momenta) then
       write (*,*) '--limit_test and --me_test are mutually exclusive diagnostics'
       stop 1
    endif
    if (read_momenta .and. accuracy.gt.0d0) then
       write (*,*) '--me_test cannot be combined with --accuracy'
       stop 1
    endif
    if (create_amplitude_library .and. (limit_test .or. read_momenta)) then
       write (*,*) '--library=create cannot be combined with --limit_test or --me_test'
       stop 1
    endif
    if (use_amplitude_library .and. limit_test) then
       write (*,*) '--limit_test requires direct amplitudes; use --library=none'
       stop 1
    endif
    if (read_momenta .and. .not.keep_processes_separate) then
       write (*,*) '--me_test requires keep_processes_separate=true'
       stop 1
    endif

    if (PS_choice.ne.1 .and. PS_choice.ne.2 .and. PS_choice.ne.3 .and. PS_choice.ne.4) then
       write (*,*) 'PS_Choice modes only 1, 2, 3 or 4',PS_choice
       stop 1
    endif
  end subroutine get_run_arguments
  
end program amplicol_generate
