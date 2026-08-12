program coupling_order_pruning_regression
  use amplitude_QCD_mod
  use particles
  implicit none

  integer,parameter :: dp=kind(1d0)
  integer,parameter :: select_pure_qcd=1
  integer,parameter :: select_interference=2
  integer,parameter :: select_pure_ew=3
  integer,parameter :: select_as_range=4
  real(kind=dp),parameter :: pi=3.14159265358979323846d0
  real(kind=dp),parameter :: alpha_s=0.118d0
  real(kind=dp),parameter :: alpha_ew=1d0/132.507d0
  real(kind=dp),parameter :: relative_tolerance=2d-11
  integer,dimension(6),parameter :: flow_with_all_sectors=[3,4,2,1,5,6]
  integer,dimension(6),parameter :: flow_without_pure_ew=[5,4,3,1,2,6]
  type(physics_model) :: model
  character(len=64) :: option

  open(unit=99,file='/dev/null',status='unknown',action='write')
  call model%init_part(173d0,0d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call model%init_vert()
  if (command_argument_count().gt.0) then
     call get_command_argument(1,option)
     if (trim(option).eq.'--emit-empty-library') then
        call check_pruning_case('locally empty pure EW flow',select_pure_ew,&
             flow_without_pure_ew,.true.,.true.)
        write (*,'(a)') 'Generated empty coupling-order pruning library'
        stop
     endif
     call fail('command line','unknown option')
  endif

  call check_pruning_case('pure QCD',select_pure_qcd,&
       flow_with_all_sectors,.false.)
  call check_pruning_case('QCD/EW interference only',select_interference,&
       flow_with_all_sectors,.false.)
  call check_pruning_case('pure EW',select_pure_ew,&
       flow_with_all_sectors,.false.)
  call check_pruning_case('aS squared-order range',select_as_range,&
       flow_with_all_sectors,.false.)
  call check_pruning_case('locally empty pure EW flow',select_pure_ew,&
       flow_without_pure_ew,.true.)
  call check_full_colour_pruning()

  write (*,'(a)') 'Coupling-order pruning regression passed'

contains

  subroutine check_full_colour_pruning()
    implicit none
    integer,parameter :: n=6
    type(amplitude_QCD) :: reference,pruned
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    integer,dimension(:),allocatable :: saved_iproc_start,saved_n_col_vals,&
         saved_col_index
    integer,dimension(:,:),allocatable :: saved_perm,saved_powers,saved_i_col_i
    integer,dimension(:,:,:),allocatable :: saved_spins,saved_row_index
    logical,dimension(:,:),allocatable :: allowed,expected_present
    real(kind=dp),dimension(:,:),allocatable :: saved_diff_col_vals
    real(kind=dp),dimension(0:3,n) :: p_before,p_after
    real(kind=dp),dimension(3) :: reference_square,pruned_square
    real(kind=dp) :: scale
    logical :: had_spins
    integer :: old_n_amps,old_n_sectors,old_n_cur,old_n_vert,old_nterms
    integer :: iamp,isector,old_terms,new_terms,expected_nterms,iqcd,iew,iacc

    part(:,1)=[1,-1,2,-2,3,-3]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    hel=[-1,1,-1,1,-1,1]
    call fill_three_line_momenta(p_before)
    p_after=1.137d0*p_before

    call reference%init(2,n,1,part,spin,orders,model)
    call reference%init_col(n,20)
    call reference%evaluate(n,p_after,hel,.false.,model)
    call pruned%init(2,n,1,part,spin,orders,model)
    call pruned%init_col(n,20)
    call pruned%evaluate(n,p_before,hel,.false.,model)

    allocate(allowed(pruned%n_sectors,pruned%n_sectors))
    allowed=.false.
    iqcd=pruned%sector_index(4,0)
    iew=pruned%sector_index(0,4)
    if (min(iqcd,iew).le.0) call fail('full-colour cross-flow interference',&
         'fixture is missing a required endpoint sector')
    ! The EW endpoint is absent from several colour rows.  This pair therefore
    ! forces imode=2 pruning to use the global endpoint union, including Gram
    ! terms between a QCD-only row and an EW-containing row.
    allowed(iqcd,iew)=.true.
    allocate(expected_present(pruned%n_amps,pruned%n_sectors))
    expected_present=.false.
    do isector=1,pruned%n_sectors
       if (isector.ne.iqcd .and. isector.ne.iew) cycle
       expected_present(:,isector)=pruned%sector_present(:,isector)
    enddo
    if (.not.any(expected_present) .or.&
         all(expected_present.eqv.pruned%sector_present)) &
         call fail('full-colour cross-flow interference','fixture cannot exercise pruning')

    do iacc=1,3
       reference_square(iacc)=selected_full_colour_square(reference,allowed,iacc)
    enddo
    if (maxval(abs(reference_square)).le.1d-30) call fail(&
         'full-colour cross-flow interference','selected reference square vanishes')
    old_n_amps=pruned%n_amps
    old_n_sectors=pruned%n_sectors
    old_n_cur=pruned%n_cur
    old_n_vert=pruned%n_vert
    old_nterms=size(pruned%sector_term_sign)
    saved_perm=pruned%perm
    had_spins=allocated(pruned%spins)
    if (had_spins) saved_spins=pruned%spins
    saved_iproc_start=pruned%iproc_start
    saved_powers=pruned%sector_powers
    saved_n_col_vals=pruned%n_col_vals
    saved_col_index=pruned%col_index
    saved_i_col_i=pruned%i_col_i
    saved_row_index=pruned%row_index
    saved_diff_col_vals=pruned%diff_col_vals

    call pruned%prune_coupling_sectors(allowed)

    if (pruned%n_amps.ne.old_n_amps .or. pruned%n_sectors.ne.old_n_sectors .or.&
         pruned%nColOrd.ne.reference%nColOrd) &
         call fail('full-colour cross-flow interference',&
              'changed an amplitude/colour axis')
    if (allocated(pruned%spins).neqv.had_spins) call fail(&
         'full-colour cross-flow interference','changed spin metadata allocation')
    if (any(pruned%perm.ne.saved_perm) .or.&
         any(pruned%iproc_start.ne.saved_iproc_start) .or.&
         any(pruned%sector_powers.ne.saved_powers) .or.&
         any(pruned%n_col_vals.ne.saved_n_col_vals) .or.&
         any(pruned%col_index.ne.saved_col_index) .or.&
         any(pruned%i_col_i.ne.saved_i_col_i) .or.&
         any(pruned%row_index.ne.saved_row_index) .or.&
         any(transfer(pruned%diff_col_vals,[0_8]).ne.&
         transfer(saved_diff_col_vals,[0_8]))) &
         call fail('full-colour cross-flow interference',&
              'changed stable colour metadata')
    if (had_spins) then
       if (any(pruned%spins.ne.saved_spins)) call fail(&
            'full-colour cross-flow interference','changed spin metadata')
    endif
    if (any(pruned%sector_present.neqv.expected_present)) &
         call fail('full-colour cross-flow interference',&
              'kept the wrong global endpoint union')

    expected_nterms=0
    do isector=1,pruned%n_sectors
       do iamp=1,pruned%n_amps
          old_terms=reference%sector_term_start(iamp,isector)-&
               reference%sector_term_start(iamp-1,isector)
          new_terms=pruned%sector_term_start(iamp,isector)-&
               pruned%sector_term_start(iamp-1,isector)
          if (expected_present(iamp,isector)) then
             if (new_terms.ne.old_terms) call fail(&
                  'full-colour cross-flow interference',&
                  'changed the terminal multiplicity of a kept sector')
             expected_nterms=expected_nterms+old_terms
          elseif (new_terms.ne.0) then
             call fail('full-colour cross-flow interference',&
                  'left an unused terminal root')
          endif
       enddo
    enddo
    if (size(pruned%sector_term_sign).ne.expected_nterms .or.&
         size(pruned%sector_term_sign).ge.old_nterms) &
         call fail('full-colour cross-flow interference',&
              'did not compact terminal roots')
    if (pruned%n_cur.ge.old_n_cur .or. pruned%n_vert.ge.old_n_vert) &
         call fail('full-colour cross-flow interference',&
              'did not reduce dead recursion trees')

    call pruned%evaluate(n,p_after,hel,.false.,model)
    if (any(pruned%sector_present.neqv.expected_present)) &
         call fail('full-colour cross-flow interference',&
              'evaluation repopulated a pruned sector')
    call compare_kept_coefficients(reference,pruned,expected_present,&
         'full-colour cross-flow interference')
    do iacc=1,3
       pruned_square(iacc)=selected_full_colour_square(pruned,allowed,iacc)
       scale=max(1d-30,abs(reference_square(iacc)),abs(pruned_square(iacc)))
       if (abs(pruned_square(iacc)-reference_square(iacc)).gt.&
            relative_tolerance*scale) then
          write (*,'(a,i2,2es24.16)') 'Full-colour accuracy changed:',iacc,&
               reference_square(iacc),pruned_square(iacc)
          call fail('full-colour cross-flow interference',&
               'selected LC/NLC/FC square changed')
       endif
    enddo

    write (*,'(a,2(a,i0,a,i0),a,i0,a,i0)') &
         'PRUNING_CASE full-colour cross-flow interference',&
         ' currents=',old_n_cur,'->',pruned%n_cur,&
         ' vertices=',old_n_vert,'->',pruned%n_vert,&
         ' roots=',old_nterms,'->',size(pruned%sector_term_sign)
  end subroutine check_full_colour_pruning

  subroutine check_pruning_case(label,selector,colour_order,expect_empty,emit_library)
    implicit none
    character(len=*),intent(in) :: label
    integer,intent(in) :: selector
    integer,dimension(6),intent(in) :: colour_order
    logical,intent(in) :: expect_empty
    logical,intent(in),optional :: emit_library
    integer,parameter :: n=6
    type(amplitude_QCD) :: reference,pruned
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    integer,dimension(:,:),allocatable :: saved_perm,saved_powers
    integer,dimension(:,:,:),allocatable :: saved_spins
    integer,dimension(:),allocatable :: saved_iproc_start
    logical,dimension(:,:),allocatable :: allowed,original_present,expected_present
    logical,dimension(:),allocatable :: expected_global,requested_global
    real(kind=dp),dimension(0:3,n) :: p_before,p_after
    real(kind=dp) :: reference_square,pruned_square,scale
    integer :: old_n_amps,old_n_sectors,old_n_cur,old_n_vert,old_nterms
    integer :: expected_nterms,iamp,isector,jsector,old_terms,new_terms
    logical :: write_library

    write_library=.false.
    if (present(emit_library)) write_library=emit_library

    part(:,1)=[1,-1,2,-2,3,-3]
    orders(:,1)=colour_order
    spin=0
    spin(0,:)=1
    spin(1,:)=[-1,1,-1,1,-1,1]
    hel=spin(1,:)
    call fill_three_line_momenta(p_before)
    ! Fill every candidate work array at a different point before pruning.
    ! The locally empty case then detects stale coefficients or currents.
    p_after=1.137d0*p_before

    call reference%init(1,n,1,part,spin,orders,model)
    call reference%evaluate(n,p_after,hel,.false.,model)
    call pruned%init(1,n,1,part,spin,orders,model)
    call pruned%evaluate(n,p_before,hel,.false.,model)

    if (reference%n_amps.ne.1 .or. pruned%n_amps.ne.1 .or.&
         reference%n_sectors.ne.pruned%n_sectors) &
         call fail(label,'unexpected fixed-flow amplitude layout')
    if (any(reference%sector_powers.ne.pruned%sector_powers) .or.&
         any(reference%sector_present.neqv.pruned%sector_present)) &
         call fail(label,'fresh reference and pruning candidate disagree')

    allocate(allowed(pruned%n_sectors,pruned%n_sectors))
    call build_selector_mask(pruned,selector,allowed)
    allocate(original_present(pruned%n_amps,pruned%n_sectors))
    allocate(expected_present(pruned%n_amps,pruned%n_sectors))
    original_present=pruned%sector_present
    expected_present=.false.
    do iamp=1,pruned%n_amps
       do isector=1,pruned%n_sectors
          if (.not.original_present(iamp,isector)) cycle
          do jsector=1,pruned%n_sectors
             if (.not.original_present(iamp,jsector)) cycle
             if (allowed(isector,jsector) .or. allowed(jsector,isector)) then
                expected_present(iamp,isector)=.true.
                exit
             endif
          enddo
       enddo
    enddo
    allocate(expected_global(pruned%n_sectors),requested_global(pruned%n_sectors))
    expected_global=any(expected_present,dim=1)
    call requested_sector_support(pruned,selector,requested_global)
    if (any(expected_global.neqv.requested_global)) &
         call fail(label,'did not expose exactly the requested sector support')
    if (expect_empty.neqv.(.not.any(expected_present))) &
         call fail(label,'does not exercise the expected empty/nonempty path')

    reference_square=selected_fixed_square(reference,allowed)
    if (.not.expect_empty .and. abs(reference_square).le.1d-30) &
         call fail(label,'selected reference square unexpectedly vanishes')

    old_n_amps=pruned%n_amps
    old_n_sectors=pruned%n_sectors
    old_n_cur=pruned%n_cur
    old_n_vert=pruned%n_vert
    old_nterms=size(pruned%sector_term_sign)
    if (allocated(pruned%perm)) saved_perm=pruned%perm
    saved_spins=pruned%spins
    saved_iproc_start=pruned%iproc_start
    saved_powers=pruned%sector_powers

    call pruned%prune_coupling_sectors(allowed)

    if (pruned%n_amps.ne.old_n_amps .or. pruned%n_sectors.ne.old_n_sectors) &
         call fail(label,'pruning changed a public amplitude/sector axis')
    if (allocated(saved_perm)) then
       if (.not.allocated(pruned%perm)) &
            call fail(label,'pruning deallocated the fixed colour order')
       if (any(pruned%perm.ne.saved_perm)) &
            call fail(label,'pruning changed the fixed colour order')
    endif
    if (any(pruned%spins.ne.saved_spins) .or.&
         any(pruned%iproc_start.ne.saved_iproc_start) .or.&
         any(pruned%sector_powers.ne.saved_powers)) &
         call fail(label,'pruning changed stable amplitude metadata')
    if (any(pruned%sector_present.neqv.expected_present)) &
         call fail(label,'kept the wrong per-amplitude coupling sectors')
    if (any(any(pruned%sector_present,dim=1).neqv.expected_global)) &
         call fail(label,'kept the wrong global coupling sectors')

    expected_nterms=0
    do isector=1,pruned%n_sectors
       do iamp=1,pruned%n_amps
          old_terms=reference%sector_term_start(iamp,isector)-&
               reference%sector_term_start(iamp-1,isector)
          new_terms=pruned%sector_term_start(iamp,isector)-&
               pruned%sector_term_start(iamp-1,isector)
          if (expected_present(iamp,isector)) then
             if (new_terms.ne.old_terms) &
                  call fail(label,'changed the terminal multiplicity of a kept sector')
             expected_nterms=expected_nterms+old_terms
          elseif (new_terms.ne.0) then
             call fail(label,'left terminal roots in an unused sector')
          endif
       enddo
    enddo
    if (size(pruned%sector_term_sign).ne.expected_nterms) &
         call fail(label,'sparse terminal arrays have the wrong compacted size')
    if (size(pruned%sector_term_sign).ge.old_nterms .and. .not.expect_empty) &
         call fail(label,'did not remove any unused terminal roots')
    if (pruned%n_cur.ge.old_n_cur) &
         call fail(label,'did not reduce the dead-tree current count')
    if (pruned%n_vert.ge.old_n_vert) &
         call fail(label,'did not reduce the dead-tree interaction count')

    call pruned%evaluate(n,p_after,hel,.false.,model)
    if (any(pruned%sector_present.neqv.expected_present)) &
         call fail(label,'evaluation repopulated a pruned sector')
    call compare_kept_coefficients(reference,pruned,expected_present,label)
    call check_empty_amplitudes(pruned,expected_present,label)
    pruned_square=selected_fixed_square(pruned,allowed)
    scale=max(1d-30,abs(reference_square),abs(pruned_square))
    if (abs(pruned_square-reference_square).gt.relative_tolerance*scale) then
       write (*,'(a,1x,a,2es24.16)') trim(label),&
            'selected fixed-flow square changed after pruning:',&
            reference_square,pruned_square
       stop 1
    endif

    if (.not.expect_empty .and. selector.eq.select_pure_qcd .and.&
         all(colour_order.eq.flow_with_all_sectors)) &
         call check_cache_round_trip(pruned,n,p_after,hel,label,.false.)

    if (expect_empty) then
       if (.not.pruned%sectors_pruned_empty) &
            call fail(label,'did not mark the locally empty recursion')
       if (pruned%n_vert.ne.0 .or. size(pruned%sector_term_sign).ne.0) &
            call fail(label,'empty recursion retained interactions or roots')
       if (any(pruned%amps.ne.cmplx(0d0,0d0,kind=dp)) .or.&
            any(pruned%amps_by_order.ne.cmplx(0d0,0d0,kind=dp))) &
            call fail(label,'empty recursion did not evaluate to exact zero')
       call pruned%optimise_evaluation(n)
       call check_cache_round_trip(pruned,n,p_after,hel,label,.true.)
       if (write_library) call pruned%create_library(n,hel,1,1,model,p_after)
    endif

    write (*,'(a,1x,a,2(a,i0,a,i0),a,i0,a,i0)') 'PRUNING_CASE',trim(label),&
         ' currents=',old_n_cur,'->',pruned%n_cur,&
         ' vertices=',old_n_vert,'->',pruned%n_vert,&
         ' roots=',old_nterms,'->',size(pruned%sector_term_sign)
  end subroutine check_pruning_case

  subroutine check_cache_round_trip(source,n,p,hel,label,expect_empty)
    implicit none
    type(amplitude_QCD),intent(inout) :: source
    integer,intent(in) :: n
    real(kind=dp),dimension(0:3,n),intent(in) :: p
    integer,dimension(n),intent(in) :: hel
    character(len=*),intent(in) :: label
    logical,intent(in) :: expect_empty
    type(amplitude_QCD) :: loaded
    complex(kind=dp),dimension(:),allocatable :: expected_amps
    complex(kind=dp),dimension(:,:),allocatable :: expected_by_order
    logical,dimension(:,:),allocatable :: expected_present,expected_retained
    integer :: cache_unit,iamp,isector
    real(kind=dp) :: scale

    expected_amps=source%amps
    expected_by_order=source%amps_by_order
    expected_present=source%sector_present
    if (.not.source%sectors_pruned .or. .not.allocated(source%sector_retained)) &
         call fail(label,'source cache lacks the retained-sector gate')
    expected_retained=source%sector_retained

    open(newunit=cache_unit,status='scratch',form='unformatted',&
         access='stream',action='readwrite')
    call source%write_init_amps_to_file(n,cache_unit)
    rewind(cache_unit)
    call loaded%read_init_amps_from_file(n,cache_unit)
    close(cache_unit)

    if (.not.loaded%sectors_pruned .or.&
         loaded%sectors_pruned_empty.neqv.expect_empty .or.&
         .not.allocated(loaded%sector_retained)) &
         call fail(label,'cache did not restore the pruning lifecycle flags')
    if (loaded%n_amps.ne.source%n_amps .or.&
         loaded%n_sectors.ne.source%n_sectors .or.&
         any(loaded%iproc_start.ne.source%iproc_start) .or.&
         any(loaded%sector_powers.ne.source%sector_powers) .or.&
         any(loaded%sector_present.neqv.expected_present) .or.&
         any(loaded%sector_retained.neqv.expected_retained)) &
         call fail(label,'cache changed pruned amplitude metadata')

    call loaded%evaluate(n,p,hel,.false.,model)
    if (any(loaded%sector_present.neqv.expected_present)) &
         call fail(label,'cached evaluation resurrected a pruned sector')
    do isector=1,loaded%n_sectors
       do iamp=1,loaded%n_amps
          scale=max(1d-30,abs(expected_by_order(iamp,isector)),&
               abs(loaded%amps_by_order(iamp,isector)))
          if (abs(loaded%amps_by_order(iamp,isector)-&
               expected_by_order(iamp,isector)).gt.relative_tolerance*scale) &
               call fail(label,'cached sector coefficient changed')
       enddo
    enddo
    scale=max(1d-30,maxval(abs(expected_amps)),maxval(abs(loaded%amps)))
    if (maxval(abs(loaded%amps-expected_amps)).gt.relative_tolerance*scale) &
         call fail(label,'cached coherent coefficient changed')
    if (expect_empty) then
       if (any(loaded%amps.ne.cmplx(0d0,0d0,kind=dp)) .or.&
            any(loaded%amps_by_order.ne.cmplx(0d0,0d0,kind=dp))) &
            call fail(label,'cached empty recursion is not exact zero')
       call loaded%optimise_evaluation(n)
    endif
  end subroutine check_cache_round_trip

  subroutine build_selector_mask(amp,selector,allowed)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer,intent(in) :: selector
    logical,dimension(:,:),intent(out) :: allowed
    integer :: isector,jsector,as2,aew2

    allowed=.false.
    do isector=1,amp%n_sectors
       do jsector=isector,amp%n_sectors
          as2=amp%sector_powers(1,isector)+amp%sector_powers(1,jsector)
          aew2=amp%sector_powers(2,isector)+amp%sector_powers(2,jsector)
          select case(selector)
          case(select_pure_qcd)
             allowed(isector,jsector)=as2.eq.8 .and. aew2.eq.0
          case(select_interference)
             ! Use only the reverse orientation to cover unordered-mask semantics.
             if (as2.eq.6 .and. aew2.eq.2) allowed(jsector,isector)=.true.
          case(select_pure_ew)
             allowed(isector,jsector)=as2.eq.0 .and. aew2.eq.8
          case(select_as_range)
             ! Squared powers are doubled integers: this is aS >= 3.
             allowed(isector,jsector)=as2.ge.6
          case default
             call fail('selector mask','unknown test selector')
          end select
       enddo
    enddo
  end subroutine build_selector_mask

  subroutine requested_sector_support(amp,selector,support)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer,intent(in) :: selector
    logical,dimension(:),intent(out) :: support
    integer :: isector

    support=.false.
    do isector=1,amp%n_sectors
       select case(selector)
       case(select_pure_qcd)
          support(isector)=all(amp%sector_powers(:,isector).eq.[4,0])
       case(select_interference,select_as_range)
          support(isector)=all(amp%sector_powers(:,isector).eq.[4,0]) .or.&
               all(amp%sector_powers(:,isector).eq.[2,2])
       case(select_pure_ew)
          support(isector)=all(amp%sector_powers(:,isector).eq.[0,4])
       end select
       support(isector)=support(isector) .and. any(amp%sector_present(:,isector))
    enddo
  end subroutine requested_sector_support

  subroutine compare_kept_coefficients(reference,pruned,expected_present,label)
    implicit none
    type(amplitude_QCD),intent(in) :: reference,pruned
    logical,dimension(:,:),intent(in) :: expected_present
    character(len=*),intent(in) :: label
    integer :: iamp,isector
    real(kind=dp) :: scale

    do isector=1,pruned%n_sectors
       do iamp=1,pruned%n_amps
          if (.not.expected_present(iamp,isector)) cycle
          scale=max(1d-30,abs(reference%amps_by_order(iamp,isector)),&
               abs(pruned%amps_by_order(iamp,isector)))
          if (abs(reference%amps_by_order(iamp,isector)-&
               pruned%amps_by_order(iamp,isector)).gt.relative_tolerance*scale) &
               call fail(label,'changed a retained fixed-flow coefficient')
       enddo
    enddo
  end subroutine compare_kept_coefficients

  subroutine check_empty_amplitudes(amp,expected_present,label)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    logical,dimension(:,:),intent(in) :: expected_present
    character(len=*),intent(in) :: label
    integer :: iamp

    do iamp=1,amp%n_amps
       if (any(expected_present(iamp,:))) cycle
       if (any(amp%curr2amp(:,iamp).ne.0)) &
            call fail(label,'empty amplitude retained a legacy terminal root')
       if (amp%amps(iamp).ne.cmplx(0d0,0d0,kind=dp) .or.&
            any(amp%amps_by_order(iamp,:).ne.cmplx(0d0,0d0,kind=dp))) &
            call fail(label,'empty amplitude did not evaluate to exact zero')
    enddo
  end subroutine check_empty_amplitudes

  real(kind=dp) function selected_fixed_square(amp,allowed)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    logical,dimension(:,:),intent(in) :: allowed
    integer :: left_sector,right_sector
    complex(kind=dp),dimension(amp%n_amps) :: left,right

    selected_fixed_square=0d0
    do left_sector=1,amp%n_sectors
       left=scaled_sector(amp,left_sector)
       do right_sector=left_sector,amp%n_sectors
          if (.not.allowed(left_sector,right_sector) .and.&
               .not.allowed(right_sector,left_sector)) cycle
          right=scaled_sector(amp,right_sector)
          if (right_sector.eq.left_sector) then
             selected_fixed_square=selected_fixed_square+&
                  sum(dble(left*conjg(right)))
          else
             selected_fixed_square=selected_fixed_square+&
                  2d0*sum(dble(left*conjg(right)))
          endif
       enddo
    enddo
  end function selected_fixed_square

  real(kind=dp) function selected_full_colour_square(amp,allowed,accuracy)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    logical,dimension(:,:),intent(in) :: allowed
    integer,intent(in) :: accuracy
    integer :: left_sector,right_sector
    complex(kind=dp),dimension(amp%n_amps) :: left,right,combined
    real(kind=dp) :: left_square,right_square

    selected_full_colour_square=0d0
    do left_sector=1,amp%n_sectors
       left=scaled_sector(amp,left_sector)
       left_square=colour_quadratic(amp,left,accuracy)
       do right_sector=left_sector,amp%n_sectors
          if (.not.allowed(left_sector,right_sector) .and.&
               .not.allowed(right_sector,left_sector)) cycle
          right=scaled_sector(amp,right_sector)
          if (right_sector.eq.left_sector) then
             selected_full_colour_square=selected_full_colour_square+left_square
          else
             right_square=colour_quadratic(amp,right,accuracy)
             combined=left+right
             selected_full_colour_square=selected_full_colour_square+&
                  colour_quadratic(amp,combined,accuracy)-left_square-right_square
          endif
       enddo
    enddo
  end function selected_full_colour_square

  real(kind=dp) function colour_quadratic(amp,values,accuracy)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    complex(kind=dp),dimension(amp%n_amps),intent(in) :: values
    integer,intent(in) :: accuracy
    integer :: row,value,entry,colour,offset
    complex(kind=dp) :: weighted,sum_for_factor

    colour_quadratic=0d0
    offset=amp%iproc_start(amp%nprocs)-1
    do row=1,amp%nColOrd
       weighted=(0d0,0d0)
       do value=1,amp%n_col_vals(accuracy)
          sum_for_factor=(0d0,0d0)
          do entry=amp%row_index(row-1,value,accuracy)+1,&
               amp%row_index(row,value,accuracy)
             colour=amp%col_index(amp%i_col_i(value,accuracy)+entry)
             sum_for_factor=sum_for_factor+values(offset+colour)
          enddo
          weighted=weighted+sum_for_factor*amp%diff_col_vals(value,accuracy)
       enddo
       colour_quadratic=colour_quadratic+&
            dble(weighted*conjg(values(offset+row)))
    enddo
  end function colour_quadratic

  function scaled_sector(amp,sector) result(values)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer,intent(in) :: sector
    complex(kind=dp),dimension(amp%n_amps) :: values
    real(kind=dp) :: factor

    factor=sqrt(4d0*pi*alpha_s)**amp%sector_powers(1,sector)*&
         sqrt(8d0*pi*alpha_ew)**amp%sector_powers(2,sector)
    values=amp%amps_by_order(:,sector)*factor
  end function scaled_sector

  subroutine fill_three_line_momenta(p)
    implicit none
    real(kind=dp),dimension(0:3,6),intent(out) :: p

    p(:,1)=[4671.200996478833d0,0d0,0d0,4671.200996478833d0]
    p(:,2)=[5452.624496459750d0,0d0,0d0,-5452.624496459750d0]
    p(:,3)=[3848.769685279069d0,-1105.1428951212934d0,&
         -1737.2006769281086d0,-3251.7412381317486d0]
    p(:,4)=[3660.1488063565753d0,2228.296688374595d0,&
         1630.5712325370819d0,2402.6278548445193d0]
    p(:,5)=[1275.9167404758375d0,-470.2433927087442d0,&
         -436.62760397349416d0,1102.8105076070963d0]
    p(:,6)=[1338.9902608271016d0,-652.9104005445573d0,&
         543.2570483645209d0,-1035.1206243007837d0]
  end subroutine fill_three_line_momenta

  subroutine fail(label,message)
    implicit none
    character(len=*),intent(in) :: label,message

    write (*,'(a,1x,a,2a)') 'Coupling-order pruning failure in',trim(label),': ',&
         trim(message)
    stop 1
  end subroutine fail

end program coupling_order_pruning_regression
