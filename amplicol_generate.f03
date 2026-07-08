
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
  implicit none
  integer :: iproc
  real(kind=8) :: weight
  integer :: i
  real(kind=8),dimension(:),allocatable :: mass,width
  character(len=80) :: filename,logfile,tag
  integer(kind=4) :: PS_choice
  integer,parameter :: nevent_hel_filter=10
  integer :: igroup
  logical,dimension(1) :: to_write
  integer,dimension(:),allocatable :: nintegrals
  integer :: ichan,iint,itmax,ncalls0,iamp
  real(kind=8),dimension(1) :: f,f_abs
  logical :: done
  real(kind=8),dimension(:,:),allocatable :: wgts
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  character(len=19) :: formatted
  logical :: create_amplitude_library,use_amplitude_library,read_momenta
  logical :: timing_enabled,time_detail_point,time_point_sample
  integer(kind=8) :: timing_point
  real(kind=8) :: tLoopBefore,tLoopAfter,tSampleBefore,tSampleAfter,tFinalBefore,tFinalAfter
  real(kind=8) :: accuracy

  call get_run_arguments()
  timing_enabled=timing_mode.ne.timing_none
  if (timing_enabled) call cpu_time(tTot_B)
  
  if (include_pdf) then
     if (use_lhapdf) then
        call InitPDFsetbyname(trim(adjustl(lhapdfset)))
        call initPDF(0)
        call setlhaparm('SILENT')
     else
        call PDF_initialise
     endif
  endif
  
  call phys_model%init_part(173d0,1.491500d0,91.188d0,2.441404d0,&
                           80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call phys_model%init_vert()

  if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
  if (use_amplitude_library) then
     call read_amplitude_lib()
  else
     call read_processes_from_file(filename)
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
                pgl(igroup)%spin,pgl(igroup)%color_orders(1,iamp),phys_model,&
                pgl(igroup)%lepton_list(1),pgl(igroup)%lepton_list)
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
             pgl(igroup)%spin,pgl(igroup)%color_orders,phys_model,&
             pgl(igroup)%lepton_list(1),pgl(igroup)%lepton_list)
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
     allocate(pgl(igroup)%hel(1:pgl(igroup)%next))
     if (keep_processes_separate) then
        allocate(pgl(igroup)%hel_fac(1:maxval(pgl(igroup)%nhel),pgl(igroup)%nproc))
     else
        allocate(pgl(igroup)%hel_fac(1:maxval(pgl(igroup)%nhel),1))
     endif
     pgl(igroup)%hel_fac=1

  enddo ! loop over phase-space-order groups

  if (use_amplitude_library) then
     if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
     call test_lib
     if (timing_mode.eq.timing_detailed) then
        call cpu_time(tAfter)
        t_lib_check=t_lib_check+tAfter-tBefore
     endif
  endif
  
  filename='Outputs/'//trim(adjustl(tag))//'events_tmp.lhe'
  open(unit=11,file=filename,action='readwrite',status='unknown')
  if (timing_mode.eq.timing_detailed) call cpu_time(tBefore)
  call write_unique_in_file(pgl_unique,unique_map,unique_map_value,abs(ncalls0))
  
  allocate(nintegrals(ngroups))
  if (keep_processes_separate) then
     nintegrals(1:ngroups)=pgl(1:ngroups)%nproc
  else
     nintegrals(1:ngroups)=1
  endif
  if (accuracy.gt.0d0) then
     call simple_integrator%init(ngroups,pgl(1:ngroups)%ndim,&
          pgl(1:ngroups)%ndim_extra,nintegrals,abs(ncalls0),abs(itmax),accuracy)
  else
     call simple_integrator%init(ngroups,pgl(1:ngroups)%ndim,&
          pgl(1:ngroups)%ndim_extra,nintegrals,abs(ncalls0),abs(itmax))
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
     if (time_detail_point) call cpu_time(tSampleBefore)
     call simple_integrator%get_points(1,ichan,iint)
     if (time_detail_point) then
        call cpu_time(tSampleAfter)
        t_Int_get=t_Int_get+(tSampleAfter-tSampleBefore)*dble(timing_sample)
     endif
     call integrand(ichan,iint,simple_integrator%x(1,1),simple_integrator%wgt(1),f(1),f_abs(1))
     if (time_detail_point) call cpu_time(tSampleBefore)
     call simple_integrator%fill_points(1,f_abs,f,to_write,done)
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
     if (to_write(1)) then
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
  if (timing_enabled) then
     call cpu_time(tLoopAfter)
     t_Int_loop=t_Int_loop+tLoopAfter-tLoopBefore
     call cpu_time(tFinalBefore)
  endif
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
  close(99)
  
contains

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

  subroutine integrand(ichan,iint,x,vol,f,f_abs)
    use scales
    use amp_lib
    implicit none
    integer,intent(in) :: ichan,iint
    real(kind=8), dimension(pgl(ichan)%ndim+pgl(ichan)%ndim_extra),intent(in) :: x
    real(kind=8),intent(in) :: vol
    real(kind=8),intent(out) :: f,f_abs
    real(kind=8), dimension(:),allocatable,save :: val,val_abs,vol_ichan
    real(kind=8),dimension(pgl(ichan)%nproc) :: colour_singlet_multichannel_weight
    integer :: ih,iproc
    real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    logical :: done,time_physics
    real(kind=8),external :: alphaspdf
    time_physics=(timing_mode.eq.timing_detailed) .and. time_point_sample
    if (create_amplitude_library) then
       if (pgl(ichan)%amps(iint)%lib_created) return
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
    val_abs=0d0

    ! Generate phase-space point based on the random numbers 'x(1:ndim)'
    if (time_physics) call cpu_time(tBefore)
    if (.not.read_momenta) then
       pgl(ichan)%ps(1)%x=x
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
    if (.not.pass_cuts(pgl(ichan))) then
       val=0d0
       if (time_physics) then
          call cpu_time(tAfter)
          t_PS= t_PS + (tAfter-tBefore)*dble(timing_sample)
       endif
       return
    endif
    endif
    pgl(ichan)%passed(iint) = pgl(ichan)%passed(iint) + 1
    call compute_multichannel_weight(ichan,pgl(ichan)%ps(1),colour_singlet_multichannel_weight)
    if (time_physics) then
       call cpu_time(tAfter)
       t_PS= t_PS + (tAfter-tBefore)*dble(timing_sample)
       tBefore=tAfter
    endif
    call compute_the_amps(iint,ichan)
    if (time_physics) then
       call cpu_time(tAfter)
       t_amp=t_amp+(tAfter-tBefore)*dble(timing_sample)
       tBefore=tAfter
    endif
    call square_the_amps(iint,ichan)
    if (time_physics) then
       call cpu_time(tAfter)
       t_mat=t_mat+(tAfter-tBefore)*dble(timing_sample)
    endif
    if ((.not. use_amplitude_library) &
         .and. pgl(ichan)%passed(iint).le.nevent_hel_filter) then
       if (time_physics) tBefore=tAfter
       call optimise_the_amplitudes(iint,ichan,done)
       if (time_physics) then
          call cpu_time(tAfter)
          t_Amp_opt=t_Amp_opt+(tAfter-tBefore)*dble(timing_sample)
       endif
       if (done) return
    endif

    if (read_momenta) then
        call perform_check(iint,ichan)
        if (pgl(ichan)%passed(iint).gt.me_points) read_momenta=.false.
    endif

    ! set scales and update alphaS
    if (time_physics) call cpu_time(tBefore)
    call set_scale(scale_choice,pgl(ichan)%next,pgl(ichan)%ps(1)%p,pgl(ichan)%processes(:,1),scale_ren)
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
    if (pgl(ichan)%amps(iint)%n_sing(1).lt.pgl(ichan)%next-2) then
       weight=weight*(4*pi*alphas)**(pgl(ichan)%next-2-pgl(ichan)%amps(iint)%n_sing(1))
    endif
    
    ! multiply by the EW coupling
    if (pgl(ichan)%amps(iint)%n_sing(1).ge.1) then
       weight=weight*(2d0*4d0*pi*alphaEW)**pgl(ichan)%amps(iint)%n_sing(1)
    endif

    if (keep_processes_separate) then
       val(1)=pgl(ichan)%amp2(1)*weight/dble(pgl(ichan)%iden(iint))
       val(1)=val(1)*colour_singlet_multichannel_weight(iint)
       call include_PDF_and_identical_procs(val,val_abs,pgl(ichan),iint)
       f_abs=sum(val_abs(1:1))
       f=sum(val(1:1))
    else
       val(1:pgl(ichan)%nproc)=pgl(ichan)%amp2(1:pgl(ichan)%nproc)*weight/dble(pgl(ichan)%iden(1:pgl(ichan)%nproc))
       val(1:pgl(ichan)%nproc)=val(1:pgl(ichan)%nproc)*colour_singlet_multichannel_weight(1:pgl(ichan)%nproc)
       call include_PDF_and_identical_procs(val,val_abs,pgl(ichan),-1)
       f_abs=sum(val_abs(1:pgl(ichan)%nproc))
       f=sum(val(1:pgl(ichan)%nproc))
    endif
    if (time_physics) then
       call cpu_time(tAfter)
       t_weight=t_weight+(tAfter-tBefore)*dble(timing_sample)
    endif
  end subroutine integrand

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
       call compute_the_amps(iint,ichan)
       call square_the_amps(iint,ichan)
       if (use_cross_process_optimisation_of_currents) then
          call pgl(ichan)%amps(iint)%optimise_evaluation(pgl(ichan)%next)
          call compute_the_amps(iint,ichan)
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

  subroutine compute_the_amps(iint,ichan)
    use amp_lib 
    implicit none
    integer,intent(in) :: iint,ichan
    if (.not. use_amplitude_library) then
       call pgl(ichan)%amps(iint)%evaluate(pgl(ichan)%next,pgl(ichan)%ps(1)%p,&
            pgl(ichan)%hel,read_proc_from_file,phys_model)
    else
       call evaluate_amp(ichan,iint,pgl(ichan)%ps(1)%p,pgl(ichan)%amps(iint)%amps)
    endif
  end subroutine compute_the_amps
    
  subroutine square_the_amps(iint,ichan)
    implicit none
    integer,intent(in) :: iint,ichan
    integer :: iproc,ih
    iproc=0
    pgl(ichan)%amp2=0d0
    if (keep_processes_separate) then
       if (use_real_gluons .and. all(pgl(ichan)%amps(iint)%n_qqbar(1:1).eq.0)) then
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih)=pgl(ichan)%amps(iint)%amps_r(ih)*&
                  pgl(ichan)%col_fac(iint)*pgl(ichan)%amps(iint)%amps_r(ih)*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
          enddo
       else
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih)=dble(pgl(ichan)%amps(iint)%amps(ih)*&
                  pgl(ichan)%col_fac(iint)*dconjg(pgl(ichan)%amps(iint)%amps(ih)))*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
          enddo
       endif
    else
       if (use_real_gluons .and. all(pgl(ichan)%amps(iint)%n_qqbar(1:1).eq.0)) then
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih)=pgl(ichan)%amps(iint)%amps_r(ih)*&
                  pgl(ichan)%col_fac(iproc)*pgl(ichan)%amps(iint)%amps_r(ih)*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
          enddo
       else
          do ih=1,pgl(ichan)%amps(iint)%n_amps
             do while (pgl(ichan)%amps(iint)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             pgl(ichan)%amp2_hel(ih)=dble(pgl(ichan)%amps(iint)%amps(ih)*&
                  pgl(ichan)%col_fac(iproc)*dconjg(pgl(ichan)%amps(iint)%amps(ih)))*pgl(ichan)%hel_fac(ih,iint)
             pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
          enddo
       endif
    endif
  end subroutine square_the_amps
  
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
    call parse_argument(filename,ncalls0,itmax,PS_choice,iseed,library,tag,read_momenta,me_points,&
         timing_arg,timing_sample_arg,accuracy)

    logfile="Outputs/"//trim(adjustl(tag))//"log_file.txt"
    open(unit=99,file=logfile,status='unknown')

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
