program coupling_order_pruning_regression
  use amplitude_QCD_mod
  use particles
  implicit none

  integer,parameter :: dp=kind(1d0)
  integer,parameter :: select_pure_qcd=1
  integer,parameter :: select_interference=2
  integer,parameter :: select_pure_ew=3
  integer,parameter :: select_as_range=4
  ! The general mixed-order recursion leaves 1124 currents and 1240 vertices
  ! after max-aS pruning for each production flow below.  Keep a conservative
  ! margin above the direct pure-QCD representation while still detecting a
  ! return of route histories to the numerical-current identity.
  integer,parameter :: max_compacted_q3_currents=400
  integer,parameter :: max_compacted_q3_vertices=700
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
     if (trim(option).eq.'--maxas-q3') then
        call check_max_as_fixed_flow_compaction()
        write (*,'(a)') 'Max-aS three-quark-line compaction regression passed'
        stop
     endif
     if (trim(option).eq.'--maxas-q3-mixed') then
        call check_max_as_mixed_three_line_compaction()
        write (*,'(a)') 'Mixed max-aS three-quark-line compaction regression passed'
        stop
     endif
     if (trim(option).eq.'--maxas-lower-lines') then
        call check_max_as_lower_line_compaction()
        write (*,'(a)') 'Lower-line max-aS compaction regression passed'
        stop
     endif
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
  call check_max_as_fixed_flow_compaction()
  call check_max_as_mixed_three_line_compaction()
  call check_max_as_lower_line_compaction()
  call check_full_colour_pruning()

  write (*,'(a)') 'Coupling-order pruning regression passed'

contains

  subroutine check_max_as_fixed_flow_compaction()
    implicit none
    integer,parameter :: n=6,nflows=3,nhelicities=2,npoints=2
    integer,dimension(n),parameter :: process=[2,-2,2,2,-2,-2]
    integer,dimension(n,nflows),parameter :: colour_orders=reshape([&
         4,1,2,5,3,6,&
         4,1,3,5,2,6,&
         2,1,3,5,4,6],[n,nflows])
    integer,dimension(n,nhelicities),parameter :: helicities=reshape([&
         -1,1,-1,-1,1,1,&
         1,-1,1,1,-1,-1],[n,nhelicities])
    real(kind=dp),dimension(npoints),parameter :: point_scale=[1d0,1.137d0]
    type(amplitude_QCD),allocatable :: reference,pruned
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    logical,dimension(:,:),allocatable :: allowed,reference_allowed
    integer,dimension(:),allocatable :: amp_map
    real(kind=dp),dimension(0:3,n) :: base_p,p
    real(kind=dp) :: coefficient_scale,square_scale,reference_square,pruned_square
    integer :: flow,helicity,point,isector,iqcd,reference_iqcd,iamp,jamp
    integer :: old_ncur,old_nvert,old_nroots,max_ncur,max_nvert
    logical :: found_nonzero

    ! These are the three fixed-flow recursions generated for the identical-
    ! flavour benchmark u u~ > u u u~ u~.  It is the case that exposed the
    ! large max-aS DAG after coupling-sector pruning.
    part(:,1)=process
    spin=0
    spin(0,:)=2
    spin(1,:)=-1
    spin(2,:)=1
    call fill_three_line_momenta(base_p)
    max_ncur=0
    max_nvert=0

    do flow=1,nflows
       orders(:,1)=colour_orders(:,flow)
       allocate(reference,pruned)
       call reference%init(1,n,1,part,spin,orders,model)
       call pruned%init(1,n,1,part,spin,orders,model,max_as_ew_order=0)
       if (reference%n_amps.ne.pruned%n_amps) &
            call fail('max-aS identical q3','pure-QCD build changed the helicity axis')
       allocate(amp_map(reference%n_amps))
       amp_map=0
       do iamp=1,reference%n_amps
          do jamp=1,pruned%n_amps
             if (any(reference%spins(:,1,iamp).ne.pruned%spins(:,1,jamp))) cycle
             amp_map(iamp)=jamp
             exit
          enddo
          if (amp_map(iamp).eq.0) call fail('max-aS identical q3',&
               'pure-QCD build lost a helicity coefficient')
       enddo
       do iamp=1,reference%n_amps
          if (count(amp_map.eq.amp_map(iamp)).ne.1) call fail(&
               'max-aS identical q3','pure-QCD helicity map is not bijective')
       enddo

       reference_iqcd=reference%sector_index(4,0)
       iqcd=pruned%sector_index(4,0)
       if (iqcd.le.0 .or. .not.any(pruned%sector_present(:,iqcd))) &
            call fail('max-aS identical q3','pure-QCD endpoint sector is missing')
       allocate(allowed(pruned%n_sectors,pruned%n_sectors))
       allocate(reference_allowed(reference%n_sectors,reference%n_sectors))
       allowed=.false.
       allowed(iqcd,iqcd)=.true.
       reference_allowed=.false.
       reference_allowed(reference_iqcd,reference_iqcd)=.true.
       old_ncur=reference%n_cur
       old_nvert=reference%n_vert
       old_nroots=size(reference%sector_term_sign)
       call pruned%prune_coupling_sectors(allowed)

       if (.not.pruned%sectors_pruned .or. pruned%sectors_pruned_empty) &
            call fail('max-aS identical q3','pure-QCD recursion was not retained')
       do isector=1,pruned%n_sectors
          if (isector.eq.iqcd) then
             if (any(pruned%sector_present(amp_map,isector).neqv.&
                  reference%sector_present(:,reference_iqcd))) &
                  call fail('max-aS identical q3','changed pure-QCD sector support')
          elseif (any(pruned%sector_present(:,isector))) then
             call fail('max-aS identical q3','retained a non-max-aS sector')
          endif
       enddo

       found_nonzero=.false.
       do point=1,npoints
          p=point_scale(point)*base_p
          do helicity=1,nhelicities
             hel=helicities(:,helicity)
             call reference%evaluate(n,p,hel,.false.,model)
             call pruned%evaluate(n,p,hel,.false.,model)
             coefficient_scale=max(1d-30,&
                  maxval(abs(reference%amps_by_order(:,reference_iqcd))),&
                  maxval(abs(pruned%amps_by_order(amp_map,iqcd))))
             if (maxval(abs(reference%amps_by_order(:,reference_iqcd))).gt.1d-30) &
                  found_nonzero=.true.
             if (maxval(abs(reference%amps_by_order(:,reference_iqcd)-&
                  pruned%amps_by_order(amp_map,iqcd))).gt.&
                  relative_tolerance*coefficient_scale) then
                  call fail('max-aS identical q3',&
                  'changed a pure-QCD coefficient after compaction')
             endif
             do isector=1,pruned%n_sectors
                if (isector.eq.iqcd) cycle
                if (any(pruned%amps_by_order(:,isector).ne.&
                     cmplx(0d0,0d0,kind=dp))) &
                     call fail('max-aS identical q3',&
                     'evaluation repopulated a discarded sector')
             enddo
             reference_square=selected_fixed_square(reference,reference_allowed)
             pruned_square=selected_fixed_square(pruned,allowed)
             square_scale=max(1d-30,abs(reference_square),abs(pruned_square))
             if (abs(reference_square-pruned_square).gt.&
                  relative_tolerance*square_scale) then
                  call fail('max-aS identical q3',&
                  'changed the max-aS fixed-flow square')
             endif
          enddo
       enddo
       if (.not.found_nonzero) call fail('max-aS identical q3',&
            'test helicities do not exercise a nonzero QCD coefficient')

       max_ncur=max(max_ncur,pruned%n_cur)
       max_nvert=max(max_nvert,pruned%n_vert)
       write (*,'(a,i0,2(a,i0,a,i0),a,i0,a,i0)') 'MAXAS_Q3_FLOW ',flow,&
            ' currents=',old_ncur,'->',pruned%n_cur,&
            ' vertices=',old_nvert,'->',pruned%n_vert,&
            ' roots=',old_nroots,'->',size(pruned%sector_term_sign)
       deallocate(allowed,reference_allowed,amp_map)
       deallocate(reference,pruned)
    enddo

    if (max_ncur.gt.max_compacted_q3_currents) call fail(&
         'max-aS identical q3','compacted recursion exceeds the current-count budget')
    if (max_nvert.gt.max_compacted_q3_vertices) call fail(&
         'max-aS identical q3','compacted recursion exceeds the interaction-count budget')
    call check_max_as_three_line_plus_gluon()
  end subroutine check_max_as_fixed_flow_compaction

  subroutine check_max_as_three_line_plus_gluon()
    implicit none
    integer,parameter :: n=7
    type(amplitude_QCD) :: reference,pure_qcd
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp) :: scale
    integer :: reference_sector,pure_sector

    part(:,1)=[1,-1,2,-2,3,-3,21]
    ! Three open colour strings, with the emitted gluon on the first string.
    orders(:,1)=[2,7,1,3,4,5,6]
    hel=[-1,1,-1,1,-1,1,-1]
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    p(:,1)=[500d0,0d0,0d0,500d0]
    p(:,2)=[500d0,0d0,0d0,-500d0]
    p(:,3)=[80.89693031749577d0,-35.363962752382264d0,&
         25.38917257558109d0,-68.18426056773842d0]
    p(:,4)=[251.36895499042657d0,-179.6909561218757d0,&
         131.2696963946096d0,-116.90072125291724d0]
    p(:,5)=[344.0807935954407d0,219.66991720940467d0,&
         -259.57359242327067d0,52.519235628094506d0]
    p(:,6)=[131.41282860492473d0,125.97195147403914d0,&
         -1.2352831219401352d0,37.401511191104284d0]
    p(:,7)=[192.24049249171213d0,-130.58694980918583d0,&
         104.15000657501999d0,95.16423500145686d0]

    call reference%init(1,n,1,part,spin,orders,model)
    call pure_qcd%init(1,n,1,part,spin,orders,model,max_as_ew_order=0)
    reference_sector=reference%sector_index(5,0)
    pure_sector=pure_qcd%sector_index(5,0)
    if (reference_sector.le.0 .or. pure_sector.le.0) call fail(&
         'max-aS q3 plus gluon','pure-QCD sector is missing')
    call reference%evaluate(n,p,hel,.false.,model)
    call pure_qcd%evaluate(n,p,hel,.false.,model)
    scale=max(1d-30,abs(reference%amps_by_order(1,reference_sector)),&
         abs(pure_qcd%amps_by_order(1,pure_sector)))
    if (abs(reference%amps_by_order(1,reference_sector)-&
         pure_qcd%amps_by_order(1,pure_sector)).gt.relative_tolerance*scale) &
         call fail('max-aS q3 plus gluon','pure-QCD coefficient changed')
  end subroutine check_max_as_three_line_plus_gluon

  subroutine check_max_as_mixed_three_line_compaction()
    implicit none

    ! Neutral radiation is the common mixed max-aS case, while the charged
    ! and Yukawa fixtures exercise the two ways in which the mandatory EW
    ! vertex can alter a fermion line.  Include an extra gluon explicitly so
    ! the compact construction is not accidentally limited to seven legs.
    call check_compact_max_as_case(&
         [2,1,2,1,3,-3,23,0],&
         [3,1,4,2,5,7,6,0],7,1,600,2300,'u d > u d s s~ z')
    call check_compact_max_as_case(&
         [2,-2,1,-1,3,-3,22,0],&
         [2,1,5,6,3,7,4,0],7,1,450,1700,'u u~ > d d~ s s~ a')
    call check_compact_max_as_case(&
         [2,-2,1,-1,3,-3,23,21],&
         [2,8,1,5,6,3,7,4],8,1,1200,5500,'u u~ > d d~ s s~ z g')
    call check_compact_max_as_case(&
         [-1,1,1,3,-3,-2,24,0],&
         [1,2,3,5,4,7,6,0],7,1,350,900,'d d~ > d s s~ u~ w+')
    call check_compact_max_as_case(&
         [5,-5,5,-5,6,-6,25,0],&
         [2,1,3,4,5,7,6,0],7,1,300,700,'b b~ > b b~ t t~ h')
  end subroutine check_max_as_mixed_three_line_compaction

  subroutine check_max_as_lower_line_compaction()
    implicit none

    ! No tree-level process in this recursion with only external gluons can
    ! radiate a colour-singlet SM boson, so the zero-line case has only a QCD
    ! endpoint.  Purely colourless amplitudes use a different interface.
    call check_compact_max_as_case(&
         [21,21,21,21,21,0,0,0],&
         [1,2,3,4,5,0,0,0],5,0,100,250,'g g > g g g',.false.)
    ! One open colour line: pure QCD, neutral/charged/Yukawa radiation, and
    ! Drell--Yan with a resolved gluon.  The last case requires the compact
    ! construction to retain the q-qbar -> neutral-vector closure which is
    ! forbidden only for a multi-line QCD-backbone representation.
    call check_compact_max_as_case(&
         [2,-2,21,21,21,0,0,0],&
         [2,3,4,5,1,0,0,0],5,0,80,120,'u u~ > g g g')
    call check_compact_max_as_case(&
         [2,-2,21,21,23,0,0,0],&
         [2,3,4,5,1,0,0,0],5,1,90,160,'u u~ > g g z')
    call check_compact_max_as_case(&
         [2,-1,21,21,24,0,0,0],&
         [2,3,4,5,1,0,0,0],5,1,60,90,'u d~ > g g w+')
    call check_compact_max_as_case(&
         [6,-6,21,21,25,0,0,0],&
         [2,3,4,5,1,0,0,0],5,1,60,80,'t t~ > g g h')
    call check_compact_max_as_case(&
         [2,-2,-11,11,21,0,0,0],&
         [2,5,3,4,1,0,0,0],5,2,60,100,'u u~ > e+ e- g',.false.)

    ! Two colour lines cover the dominant jet-production shapes.  Include
    ! single neutral, charged-current, and top-Yukawa emissions as well as a
    ! genuine diboson VBS-like external state.
    call check_compact_max_as_case(&
         [2,1,2,1,21,0,0,0],&
         [3,1,4,5,2,0,0,0],5,0,60,80,'u d > u d g')
    call check_compact_max_as_case(&
         [2,2,2,2,21,0,0,0],&
         [3,1,4,5,2,0,0,0],5,0,90,140,'u u > u u g')
    call check_compact_max_as_case(&
         [2,1,2,1,21,22,0,0],&
         [3,1,4,5,6,2,0,0],6,1,160,350,'u d > u d g a')
    call check_compact_max_as_case(&
         [2,-2,2,-2,23,0,0,0],&
         [2,1,3,5,4,0,0,0],5,1,140,290,'u u~ > u u~ z')
    call check_compact_max_as_case(&
         [2,1,1,1,21,24,0,0],&
         [3,1,4,5,6,2,0,0],6,1,190,330,'u d > d d g w+')
    call check_compact_max_as_case(&
         [6,2,6,2,21,25,0,0],&
         [3,1,4,5,6,2,0,0],6,1,70,90,'t u > t u g h')
    call check_compact_max_as_case(&
         [2,1,2,1,24,-24,0,0],&
         [3,1,4,5,6,2,0,0],6,2,280,700,'u d > u d w+ w-')
    call check_compact_max_as_case(&
         [2,1,2,1,-11,11,0,0],&
         [3,1,4,5,6,2,0,0],6,2,220,500,'u d > u d e+ e-')
  end subroutine check_max_as_lower_line_compaction

  subroutine check_compact_max_as_case(process_storage,order_storage,n,target_ew,&
       max_currents,max_vertices,label,require_reduction)
    implicit none
    integer,dimension(8),intent(in) :: process_storage,order_storage
    integer,intent(in) :: n,target_ew,max_currents,max_vertices
    character(len=*),intent(in) :: label
    logical,optional,intent(in) :: require_reduction
    integer,parameter :: npoints=2,nhelicities=2
    type(amplitude_QCD),allocatable :: reference,compact
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    integer,dimension(n,nhelicities) :: helicities
    integer,dimension(:),allocatable :: amp_map
    logical,dimension(:,:),allocatable :: reference_allowed,compact_allowed
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(npoints),parameter :: point_scale=[1d0,1.071d0]
    real(kind=dp) :: coefficient_scale,square_scale,reference_square,compact_square
    integer :: i,iamp,jamp,point,helicity,reference_sector,compact_sector
    integer :: expected_gs,reference_currents,reference_vertices
    logical :: found_nonzero,must_reduce

    must_reduce=.true.
    if (present(require_reduction)) must_reduce=require_reduction
    part(:,1)=process_storage(1:n)
    orders(:,1)=order_storage(1:n)
    helicities=0
    spin=0
    do i=1,n
       select case(model%get_spin(part(i,1)))
       case(1)
          spin(0,i)=1
       case(2)
          spin(0,i)=2
          spin(1:2,i)=[-1,1]
          helicities(i,1)=-1
          if (part(i,1).lt.0) helicities(i,1)=1
          helicities(i,2)=-helicities(i,1)
       case(3)
          spin(0,i)=3
          spin(1:3,i)=[-1,0,1]
          helicities(i,:)=[-1,0]
       case default
          call fail(label,'fixture contains an unsupported spin')
       end select
    enddo
    allocate(reference,compact)
    call reference%init(1,n,1,part,spin,orders,model)
    call compact%init(1,n,1,part,spin,orders,model,&
         max_as_ew_order=target_ew)
    expected_gs=n-2-target_ew
    reference_sector=reference%sector_index(expected_gs,target_ew)
    compact_sector=compact%sector_index(expected_gs,target_ew)
    if (reference_sector.le.0 .or. compact_sector.le.0) &
         call fail(label,'maximum-aS endpoint sector is missing')
    if (.not.any(reference%sector_present(:,reference_sector)) .or.&
         .not.any(compact%sector_present(:,compact_sector))) &
         call fail(label,'maximum-aS endpoint sector has no coefficients')

    allocate(amp_map(reference%n_amps))
    amp_map=0
    do iamp=1,reference%n_amps
       if (.not.reference%sector_present(iamp,reference_sector)) cycle
       do jamp=1,compact%n_amps
          if (.not.compact%sector_present(jamp,compact_sector)) cycle
          if (any(reference%spins(:,1,iamp).ne.compact%spins(:,1,jamp))) cycle
          amp_map(iamp)=jamp
          exit
       enddo
       if (amp_map(iamp).eq.0) call fail(label,&
            'compact build lost a helicity coefficient')
    enddo
    do iamp=1,reference%n_amps
       if (amp_map(iamp).eq.0) cycle
       if (count(amp_map.eq.amp_map(iamp)).ne.1) call fail(label,&
            'compact helicity map is not bijective')
    enddo
    if (count(amp_map.gt.0).ne.count(compact%sector_present(:,compact_sector))) &
         call fail(label,'compact build added an unmatched helicity coefficient')

    allocate(reference_allowed(reference%n_sectors,reference%n_sectors))
    allocate(compact_allowed(compact%n_sectors,compact%n_sectors))
    reference_allowed=.false.
    reference_allowed(reference_sector,reference_sector)=.true.
    compact_allowed=.false.
    compact_allowed(compact_sector,compact_sector)=.true.
    call compact%prune_coupling_sectors(compact_allowed)
    if (.not.compact%sectors_pruned .or. compact%sectors_pruned_empty) &
         call fail(label,'compact maximum-aS recursion was not retained')
    do iamp=1,reference%n_amps
       if (.not.reference%sector_present(iamp,reference_sector)) cycle
       if (amp_map(iamp).le.0 .or.&
            .not.compact%sector_present(amp_map(iamp),compact_sector)) &
            call fail(label,'compact build changed maximum-aS helicity support')
    enddo

    found_nonzero=.false.
    do point=1,npoints
       call fill_generic_mixed_momenta(n,part(:,1),point_scale(point),p)
       do helicity=1,nhelicities
          hel=helicities(:,helicity)
          call reference%evaluate(n,p,hel,.false.,model)
          call compact%evaluate(n,p,hel,.false.,model)
          coefficient_scale=max(1d-30,&
               maxval(abs(reference%amps_by_order(:,reference_sector))))
          do iamp=1,reference%n_amps
             if (amp_map(iamp).le.0) cycle
             coefficient_scale=max(coefficient_scale,&
                  abs(compact%amps_by_order(amp_map(iamp),compact_sector)))
          enddo
          if (coefficient_scale.gt.1d-30) found_nonzero=.true.
          do iamp=1,reference%n_amps
             if (amp_map(iamp).le.0) cycle
             if (abs(reference%amps_by_order(iamp,reference_sector)-&
                  compact%amps_by_order(amp_map(iamp),compact_sector)).gt.&
                  relative_tolerance*coefficient_scale) &
                  call fail(label,'compact build changed a maximum-aS coefficient')
          enddo
          reference_square=selected_fixed_square(reference,reference_allowed)
          compact_square=selected_fixed_square(compact,compact_allowed)
          square_scale=max(1d-30,abs(reference_square),abs(compact_square))
          if (abs(reference_square-compact_square).gt.&
               relative_tolerance*square_scale) &
               call fail(label,'compact build changed the maximum-aS square')
          ! The production generator performs this numerical DAG
          ! optimisation after its first helicity samples.  Continue the
          ! multi-point comparison on the optimised graph so accidental
          ! one-point mergers of distinct particle/momentum currents are
          ! caught here as well.
          if (point.eq.1 .and. helicity.eq.1) &
               call compact%optimise_evaluation(n)
       enddo
    enddo
    if (.not.found_nonzero) call fail(label,&
         'fixture helicity does not exercise a nonzero maximum-aS coefficient')

    reference_currents=reference%n_cur
    reference_vertices=reference%n_vert
    if (compact%n_cur.gt.reference_currents .or.&
         compact%n_vert.gt.reference_vertices) &
         call fail(label,'compact build enlarged the recursion DAG')
    if (must_reduce .and. (compact%n_cur.eq.reference_currents .or.&
         compact%n_vert.eq.reference_vertices)) &
         call fail(label,'compact build did not reduce the recursion DAG')
    if (compact%n_cur.gt.max_currents) &
         call fail(label,'compact build exceeds the current-count budget')
    if (compact%n_vert.gt.max_vertices) &
         call fail(label,'compact build exceeds the interaction-count budget')
    write (*,'(a,1x,a,2(a,i0,a,i0))') 'MAXAS_COMPACT',trim(label),&
         ' currents=',reference_currents,'->',compact%n_cur,&
         ' vertices=',reference_vertices,'->',compact%n_vert
    deallocate(reference_allowed,compact_allowed,amp_map)
    deallocate(reference,compact)
  end subroutine check_compact_max_as_case

  subroutine fill_generic_mixed_momenta(n,process,momentum_scale,p)
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(in) :: process
    real(kind=dp),intent(in) :: momentum_scale
    real(kind=dp),dimension(0:3,n),intent(out) :: p
    real(kind=dp),dimension(3,3),parameter :: directions=reshape([&
         0.3d0,0.4d0,sqrt(0.75d0),&
         0.6d0,-0.2d0,sqrt(0.60d0),&
         -0.1d0,0.7d0,sqrt(0.50d0)],[3,3])
    real(kind=dp) :: q,total_energy,mass1,mass2,sqrts,s,lambda
    integer :: i,j,pair

    p=0d0
    total_energy=0d0
    if (mod(n-2,2).eq.1) then
       q=momentum_scale*197d0
       p(1:3,3)=q*[1d0,0d0,0d0]
       p(1:3,4)=q*[-0.5d0,sqrt(0.75d0),0d0]
       p(1:3,5)=q*[-0.5d0,-sqrt(0.75d0),0d0]
       do i=3,5
          p(0,i)=sqrt(q**2+model%get_mass(process(i))**2)
          total_energy=total_energy+p(0,i)
       enddo
       i=6
    else
       i=3
    endif
    pair=0
    do while(i.le.n)
       pair=pair+1
       j=i+1
       if (j.gt.n .or. pair.gt.size(directions,2)) &
            call fail('mixed max-aS momenta','invalid generic point layout')
       q=momentum_scale*(113d0+41d0*dble(pair)+7d0*dble(n))
       p(1:3,i)=q*directions(:,pair)
       p(1:3,j)=-p(1:3,i)
       p(0,i)=sqrt(q**2+model%get_mass(process(i))**2)
       p(0,j)=sqrt(q**2+model%get_mass(process(j))**2)
       total_energy=total_energy+p(0,i)+p(0,j)
       i=i+2
    enddo
    ! Keep massive incoming fixtures (the top-Yukawa controls) on shell as
    ! well.  The generated final state fixes sqrt(s); solve the two-body
    ! incoming centre-of-mass kinematics for the corresponding beam momentum.
    mass1=model%get_mass(process(1))
    mass2=model%get_mass(process(2))
    sqrts=total_energy
    s=sqrts**2
    lambda=(s-(mass1+mass2)**2)*(s-(mass1-mass2)**2)
    if (lambda.lt.-1d-10*s**2) call fail(&
         'mixed max-aS momenta','final-state energy is below incoming threshold')
    q=sqrt(max(0d0,lambda))/(2d0*sqrts)
    p(:,1)=[(s+mass1**2-mass2**2)/(2d0*sqrts),0d0,0d0,q]
    p(:,2)=[(s+mass2**2-mass1**2)/(2d0*sqrts),0d0,0d0,-q]
  end subroutine fill_generic_mixed_momenta

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
