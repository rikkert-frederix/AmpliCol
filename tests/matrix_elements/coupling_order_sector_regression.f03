program coupling_order_sector_regression
  use amplitude_QCD_mod
  use particles
  implicit none

  integer,parameter :: dp=kind(1d0)
  real(kind=dp),parameter :: pi=3.14159265358979323846d0
  real(kind=dp),parameter :: alpha_s=0.118d0
  real(kind=dp),parameter :: alpha_ew=1d0/132.507d0
  real(kind=dp),parameter :: relative_tolerance=2d-11
  type(physics_model) :: model
  character(len=64) :: option
  logical :: generator_lc_mismatch,fixed_order_map_mismatch

  open(unit=99,file='/dev/null',status='unknown',action='write')
  call model%init_part(173d0,0d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call model%init_vert()
  generator_lc_mismatch=.false.
  fixed_order_map_mismatch=.false.
  if (command_argument_count().gt.0) then
     call get_command_argument(1,option)
     if (trim(option).eq.'--emit-library') then
        call emit_three_line_library()
        stop
     endif
     write (*,*) 'Unknown option:',trim(option)
     stop 1
  endif

  call check_vertex_orders()
  call check_charged_closure_crossing()
  call check_higgs_vector_boson_fusion()
  call check_same_sign_vbs_reachability()
  call check_three_quark_lines()
  call check_ud_scattering()
  call check_uu_scattering()
  if (generator_lc_mismatch .or. fixed_order_map_mismatch) then
     if (fixed_order_map_mismatch) write (*,'(a)') &
          'Fixed-order coefficients disagree with physical imode=2 flows'
     if (generator_lc_mismatch) write (*,'(a)') &
          'Generator fixed-order leading-colour slices disagree with imode=2'
     stop 1
  endif
  write (*,'(a)') 'Coupling-order sector regression passed'

contains

  subroutine emit_three_line_library()
    implicit none
    integer,parameter :: n=6
    type(amplitude_QCD) :: amp
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    real(kind=dp),dimension(0:3,n) :: p

    part(:,1)=[1,-1,2,-2,3,-3]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    hel=[-1,1,-1,1,-1,1]
    call fill_three_line_momenta(p)
    call amp%init(2,n,1,part,spin,orders,model)
    call amp%evaluate(n,p,hel,.false.,model)
    call amp%create_library(n,hel,1,1,model,p)
    write (*,'(a)') 'Generated coupling-order multi-root regression library'
  end subroutine emit_three_line_library

  subroutine check_vertex_orders()
    implicit none

    call require_vertex_order([-2,2,21],1,0,'QCD quark current')
    call require_vertex_order([2,-2,23],0,1,'crossed EW quark current')
    call require_vertex_order([24,-24,127],0,2,'EW-squared auxiliary-scalar vertex')
    call require_vertex_order([127,25,25],0,0,'order-zero auxiliary-scalar vertex')
  end subroutine check_vertex_orders

  subroutine require_vertex_order(particles,n_gs,n_ew,label)
    implicit none
    integer,dimension(3),intent(in) :: particles
    integer,intent(in) :: n_gs,n_ew
    character(len=*),intent(in) :: label
    integer :: vertex

    do vertex=1,model%nint
       if (.not.all(model%vertex_list(vertex)%particles.eq.particles)) cycle
       if (model%vertex_list(vertex)%n_gs.ne.n_gs .or.&
            model%vertex_list(vertex)%n_ew.ne.n_ew) then
          write (*,*) trim(label),' has the wrong coupling order:',&
               model%vertex_list(vertex)%n_gs,model%vertex_list(vertex)%n_ew
          stop 1
       endif
       return
    enddo
    write (*,*) 'Required model vertex is missing:',trim(label),particles
    stop 1
  end subroutine require_vertex_order

  subroutine check_uu_scattering()
    implicit none
    integer,parameter :: n=4
    type(amplitude_QCD) :: amp
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    integer,dimension(2,2) :: expected_powers
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(3) :: slices
    ! MadGraph 3.2.0 standalone SMATRIX_SPLITORDERS at the point below.
    real(kind=dp),dimension(3),parameter :: madgraph_slices=&
         [1.650449737914214d3,1.129825144296226d2,4.317013764726654d1]
    real(kind=dp) :: coherent_squared,scale
    integer :: slice

    part(:,1)=[2,2,2,2]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    hel=[-1,-1,-1,-1]
    p(:,1)=[500d0,0d0,0d0,500d0]
    p(:,2)=[500d0,0d0,0d0,-500d0]
    p(:,3)=[500d0,300d0,0d0,400d0]
    p(:,4)=[500d0,-300d0,0d0,-400d0]

    call amp%init(2,n,1,part,spin,orders,model)
    call amp%init_col(n,20)
    call amp%evaluate(n,p,hel,.false.,model)

    expected_powers(:,1)=[2,0]
    expected_powers(:,2)=[0,2]
    call check_sector_layout(amp,expected_powers,'u u > u u')
    call check_fixed_order_sector_map(amp,n,part,p,hel,'u u > u u')
    slices(1)=exact_squared_slice(amp,4,0)
    slices(2)=exact_squared_slice(amp,2,2)
    slices(3)=exact_squared_slice(amp,0,4)
    coherent_squared=coherent_sector_square(amp)
    scale=max(1d-30,abs(coherent_squared),sum(abs(slices)))
    if (abs(sum(slices)-coherent_squared).gt.relative_tolerance*scale) then
       write (*,*) 'u u sector slices do not reconstruct the coherent square:',&
            slices,coherent_squared
       stop 1
    endif
    if (slices(1).le.0d0 .or. slices(3).le.0d0) then
       write (*,*) 'u u diagonal coupling-order slices are not positive:',slices
       stop 1
    endif
    if (abs(slices(2)).le.1d-12*max(slices(1),slices(3))) then
       write (*,*) 'u u QCD-EW interference slice unexpectedly vanished:',slices
       stop 1
    endif
    do slice=1,size(slices)
       if (abs(slices(slice)/madgraph_slices(slice)-1d0).gt.relative_tolerance) then
          write (*,*) 'u u exact-order slice disagrees with MadGraph:',&
               slice,slices(slice),madgraph_slices(slice)
          stop 1
       endif
    enddo
    write (*,'(a,3es24.16)') 'UU_EXACT_AS2_AEW2_SLICES=',slices
  end subroutine check_uu_scattering

  subroutine check_ud_scattering()
    implicit none
    integer,parameter :: n=4
    type(amplitude_QCD) :: amp
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    integer,dimension(2,2) :: expected_powers
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(3) :: slices
    ! MadGraph 3.2.0 standalone SMATRIX_SPLITORDERS at the point below.
    real(kind=dp),dimension(3),parameter :: madgraph_slices=&
         [1.759031957513834d3,1.116852647290016d2,3.037408816836559d1]
    integer :: slice

    part(:,1)=[2,1,2,1]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    hel=[-1,-1,-1,-1]
    p(:,1)=[500d0,0d0,0d0,500d0]
    p(:,2)=[500d0,0d0,0d0,-500d0]
    p(:,3)=[500d0,300d0,0d0,400d0]
    p(:,4)=[500d0,-300d0,0d0,-400d0]

    call amp%init(2,n,1,part,spin,orders,model)
    call amp%init_col(n,20)
    call amp%evaluate(n,p,hel,.false.,model)
    expected_powers(:,1)=[2,0]
    expected_powers(:,2)=[0,2]
    call check_sector_layout(amp,expected_powers,'u d > u d')
    call check_fixed_order_sector_map(amp,n,part,p,hel,'u d > u d')
    slices(1)=exact_squared_slice(amp,4,0)
    slices(2)=exact_squared_slice(amp,2,2)
    slices(3)=exact_squared_slice(amp,0,4)
    do slice=1,size(slices)
       if (abs(slices(slice)/madgraph_slices(slice)-1d0).gt.relative_tolerance) then
          write (*,*) 'u d exact-order slice disagrees with MadGraph:',&
               slice,slices(slice),madgraph_slices(slice)
          stop 1
       endif
    enddo
    write (*,'(a,3es24.16)') 'UD_EXACT_AS2_AEW2_SLICES=',slices
  end subroutine check_ud_scattering

  subroutine check_higgs_vector_boson_fusion()
    implicit none
    integer,parameter :: n=5
    ! MadGraph 3.7.2 standalone pure-EW result at the point below.
    real(kind=dp),parameter :: madgraph_pure_ew=2.246144714146901d-6
    type(amplitude_QCD) :: amp
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    integer,dimension(2,1) :: expected_powers
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp) :: pure_ew

    part(:,1)=[2,2,2,2,25]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    hel=[-1,-1,-1,-1,0]
    p(:,1)=[500d0,0d0,0d0,500d0]
    p(:,2)=[500d0,0d0,0d0,-500d0]
    p(:,3)=[437.5d0,300d0,0d0,sqrt(101406.25d0)]
    p(:,4)=[437.5d0,-300d0,0d0,-sqrt(101406.25d0)]
    p(:,5)=[125d0,0d0,0d0,0d0]

    call amp%init(2,n,1,part,spin,orders,model)
    call amp%init_col(n,20)
    call amp%evaluate(n,p,hel,.false.,model)
    expected_powers(:,1)=[0,3]
    call check_sector_layout(amp,expected_powers,'u u > u u h')
    pure_ew=exact_squared_slice(amp,0,6)
    if (pure_ew.le.0d0 .or.&
         abs(pure_ew/madgraph_pure_ew-1d0).gt.relative_tolerance) then
       write (*,*) 'Pure-EW Higgs VBF fixed point disagrees with MadGraph:',&
            pure_ew,madgraph_pure_ew
       stop 1
    endif
    write (*,'(a,es24.16)') 'UUH_EXACT_AEW3_SQUARE=',pure_ew
  end subroutine check_higgs_vector_boson_fusion

  subroutine check_same_sign_vbs_reachability()
    implicit none
    integer,parameter :: n=6
    ! MadGraph 3.7.2 standalone split orders. Its parameter card rounds MW,
    ! so use a looser comparison than for identical model parameters.
    real(kind=dp),dimension(3),parameter :: madgraph_slices=&
         [3.451386752919821d-14,1.219336701511415d-15,3.359791832459152d-16]
    type(amplitude_QCD) :: amp
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    integer,dimension(2,2) :: expected_powers
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(3) :: slices
    real(kind=dp) :: coherent_squared,scale
    integer :: sector,slice
    logical :: found_qcd_induced,found_pure_ew

    part(:,1)=[2,2,1,1,24,24]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    hel=[-1,-1,-1,-1,0,0]
    p(:,1)=[1000d0,0d0,0d0,1000d0]
    p(:,2)=[1000d0,0d0,0d0,-1000d0]
    p(:,3)=[400d0,400d0,0d0,0d0]
    p(:,4)=[400d0,-400d0,0d0,0d0]
    p(:,5)=[600d0,0d0,sqrt(600d0**2-model%get_mass(24)**2),0d0]
    p(:,6)=[600d0,0d0,-sqrt(600d0**2-model%get_mass(24)**2),0d0]
    call amp%init(2,n,1,part,spin,orders,model)
    found_qcd_induced=.false.
    found_pure_ew=.false.
    do sector=1,amp%n_sectors
       if (sum(amp%sector_powers(:,sector)).ne.n-2) then
          write (*,*) 'VBS sector violates the tree-order invariant:',&
               amp%sector_powers(:,sector),n-2
          stop 1
       endif
       if (all(amp%sector_powers(:,sector).eq.[2,2])) found_qcd_induced=.true.
       if (all(amp%sector_powers(:,sector).eq.[0,4])) found_pure_ew=.true.
    enddo
    if (.not.found_qcd_induced .or. .not.found_pure_ew) then
       write (*,*) 'Missing QCD-induced or pure-EW same-sign VBS sector:'
       do sector=1,amp%n_sectors
          write (*,*) amp%sector_powers(:,sector)
       enddo
       stop 1
    endif
    call amp%init_col(n,20)
    call amp%evaluate(n,p,hel,.false.,model)
    expected_powers(:,1)=[2,2]
    expected_powers(:,2)=[0,4]
    call check_sector_layout(amp,expected_powers,'u u > d d w+ w+')
    call check_fixed_order_sector_map(amp,n,part,p,hel,&
         'u u > d d w+ w+')
    slices(1)=exact_squared_slice(amp,4,4)
    slices(2)=exact_squared_slice(amp,2,6)
    slices(3)=exact_squared_slice(amp,0,8)
    coherent_squared=coherent_sector_square(amp)
    scale=max(1d-30,abs(coherent_squared),sum(abs(slices)))
    if (abs(sum(slices)-coherent_squared).gt.relative_tolerance*scale .or.&
         slices(1).le.0d0 .or. slices(3).le.0d0) then
       write (*,*) 'Same-sign VBS sectors do not form a valid coherent square:',&
            slices,coherent_squared
       stop 1
    endif
    do slice=1,size(slices)
       if (abs(slices(slice)/madgraph_slices(slice)-1d0).gt.1d-6) then
          write (*,*) 'Same-sign VBS exact-order slice disagrees with MadGraph:',&
               slice,slices(slice),madgraph_slices(slice)
          stop 1
       endif
    enddo
    write (*,'(a,3es24.16)') 'VBS_EXACT_AS2_AEW2_SLICES=',slices
  end subroutine check_same_sign_vbs_reachability

  subroutine check_three_quark_lines()
    implicit none
    integer,parameter :: n=6
    ! MadGraph 3.2.0 standalone SMATRIX_SPLITORDERS at the point below.
    real(kind=dp),dimension(5),parameter :: madgraph_slices=&
         [2.105741599135726d-13,-1.812427662010086d-13,&
         5.779090350078999d-14,-1.248928669293392d-15,&
         3.588064623776670d-16]
    type(amplitude_QCD) :: amp
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    integer,dimension(2,3) :: expected_powers
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(5) :: slices
    real(kind=dp) :: coherent_squared,scale
    integer :: slice

    part(:,1)=[1,-1,2,-2,3,-3]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    hel=[-1,1,-1,1,-1,1]
    call fill_three_line_momenta(p)

    call amp%init(2,n,1,part,spin,orders,model)
    call amp%init_col(n,20)
    call amp%evaluate(n,p,hel,.false.,model)

    expected_powers(:,1)=[4,0]
    expected_powers(:,2)=[2,2]
    expected_powers(:,3)=[0,4]
    call check_sector_layout(amp,expected_powers,'d d~ > u u~ s s~')
    call check_fixed_order_sector_map(amp,n,part,p,hel,&
         'd d~ > u u~ s s~')
    slices(1)=exact_squared_slice(amp,8,0)
    slices(2)=exact_squared_slice(amp,6,2)
    slices(3)=exact_squared_slice(amp,4,4)
    slices(4)=exact_squared_slice(amp,2,6)
    slices(5)=exact_squared_slice(amp,0,8)
    coherent_squared=coherent_sector_square(amp)
    scale=max(1d-30,abs(coherent_squared),sum(abs(slices)))
    if (abs(sum(slices)-coherent_squared).gt.relative_tolerance*scale) then
       write (*,*) 'Three-line sector slices do not reconstruct the coherent square:',&
            slices,coherent_squared
       stop 1
    endif
    if (slices(1).le.0d0 .or. slices(3).le.0d0 .or. slices(5).le.0d0) then
       write (*,*) 'Three-line diagonal-order slices are not positive:',slices
       stop 1
    endif
    do slice=1,size(slices)
       if (abs(slices(slice)/madgraph_slices(slice)-1d0).gt.relative_tolerance) then
          write (*,*) 'Three-line exact-order slice disagrees with MadGraph:',&
               slice,slices(slice),madgraph_slices(slice)
          stop 1
       endif
    enddo
    write (*,'(a,5es24.16)') 'THREE_LINE_EXACT_AS2_AEW2_SLICES=',slices
  end subroutine check_three_quark_lines

  subroutine check_charged_closure_crossing()
    implicit none
    integer,parameter :: n=3
    real(kind=dp),parameter :: energy1=100d0
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp) :: energy2,mw

    mw=model%get_mass(24)
    energy2=mw**2/(4d0*energy1)
    p(:,1)=[energy1,0d0,0d0,energy1]
    p(:,2)=[energy2,0d0,0d0,-energy2]
    p(:,3)=[energy1+energy2,0d0,0d0,energy1-energy2]

    ! In imode=2 the W closes the recursion through a type-21/22
    ! fermion-antifermion current.  The fixed colour order closes the same
    ! amplitude on a fermion through type 10/11.  Cover both charges and both
    ! spinor-chain orientations.
    call compare_charged_roots([2,-1,24],[-1,1,-1],[2,3,1],p,&
         'u d~ > w+')
    call compare_charged_roots([-1,2,24],[1,-1,1],[1,3,2],p,&
         'd~ u > w+')
    call compare_charged_roots([1,-2,-24],[-1,1,-1],[2,3,1],p,&
         'd u~ > w-')
    call compare_charged_roots([-2,1,-24],[1,-1,1],[1,3,2],p,&
         'u~ d > w-')
  end subroutine check_charged_closure_crossing

  subroutine compare_charged_roots(particles_in,helicities,colour_order,p,label)
    implicit none
    integer,parameter :: n=3
    integer,dimension(n),intent(in) :: particles_in,helicities,colour_order
    real(kind=dp),dimension(0:3,n),intent(in) :: p
    character(len=*),intent(in) :: label
    type(amplitude_QCD) :: all_orders,fixed_order
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    real(kind=dp) :: scale

    part(:,1)=particles_in
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    call all_orders%init(2,n,1,part,spin,orders,model)
    call all_orders%evaluate(n,p,helicities,.false.,model)

    orders(:,1)=colour_order
    spin(1,:)=helicities
    call fixed_order%init(1,n,1,part,spin,orders,model)
    call fixed_order%evaluate(n,p,helicities,.false.,model)
    if (all_orders%n_amps.ne.1 .or. all_orders%n_sectors.ne.1 .or.&
         fixed_order%n_amps.ne.1 .or. fixed_order%n_sectors.ne.1 .or.&
         any(all_orders%sector_powers(:,1).ne.fixed_order%sector_powers(:,1))) then
       write (*,*) trim(label),' charged-current root layouts disagree'
       stop 1
    endif
    scale=max(1d-30,abs(all_orders%amps_by_order(1,1)),&
         abs(fixed_order%amps_by_order(1,1)))
    if (abs(all_orders%amps_by_order(1,1)-fixed_order%amps_by_order(1,1)).gt.&
         relative_tolerance*scale) then
       write (*,*) trim(label),' type-21/22 closure disagrees with type 10/11:',&
            all_orders%amps_by_order(1,1),fixed_order%amps_by_order(1,1)
       stop 1
    endif
  end subroutine compare_charged_roots

  subroutine check_fixed_order_sector_map(all_orders,n,part,p,helicities,label)
    implicit none
    type(amplitude_QCD),intent(in) :: all_orders
    integer,intent(in) :: n
    integer,dimension(n,1),intent(in) :: part
    real(kind=dp),dimension(0:3,n),intent(in) :: p
    integer,dimension(n),intent(in) :: helicities
    character(len=*),intent(in) :: label
    type(amplitude_QCD),allocatable :: fixed_order
    integer,dimension(n,1) :: orders
    integer,dimension(0:3,n) :: spin
    complex(kind=dp),dimension(all_orders%n_amps,all_orders%n_sectors) :: fixed
    complex(kind=dp),allocatable :: fixed_components(:,:)
    logical,allocatable :: fixed_present(:,:)
    real(kind=dp),allocatable :: generator_slices(:,:)
    integer :: flow,sector,fixed_sector,ncoloured,position,singlet
    real(kind=dp) :: scale

    if (all_orders%n_amps.ne.all_orders%nColOrd) then
       write (*,*) trim(label),' cannot compare fixed orders to physical flows:',&
            all_orders%n_amps,all_orders%nColOrd
       stop 1
    endif
    ncoloured=size(all_orders%perm,1)
    fixed=(0d0,0d0)
    allocate(generator_slices(0:2*maxval(all_orders%sector_powers(1,:)),&
         0:2*maxval(all_orders%sector_powers(2,:))))
    generator_slices=0d0
    spin=0
    spin(0,:)=1
    spin(1,:)=helicities
    do flow=1,all_orders%nColOrd
       orders=0
       orders(1:ncoloured-1,1)=all_orders%perm(1:ncoloured-1,flow)
       position=ncoloured
       do singlet=1,n
          if (.not.model%is_singlet(part(singlet,1))) cycle
          orders(position,1)=singlet
          position=position+1
       enddo
       orders(n,1)=all_orders%perm(ncoloured,flow)
       if (position.ne.n) then
          write (*,*) trim(label),' failed to construct a fixed colour order:',&
               orders(:,1)
          stop 1
       endif
       allocate(fixed_order)
       call fixed_order%init(1,n,1,part,spin,orders,model)
       call fixed_order%evaluate(n,p,helicities,.false.,model)
       if (fixed_order%n_amps.ne.1) then
          write (*,*) trim(label),' fixed order did not compact to one physical coefficient:',&
               flow,fixed_order%n_amps
          fixed_order_map_mismatch=.true.
       endif
       allocate(fixed_components(fixed_order%n_amps,all_orders%n_sectors))
       allocate(fixed_present(fixed_order%n_amps,all_orders%n_sectors))
       fixed_components=(0d0,0d0)
       fixed_present=.false.
       do fixed_sector=1,fixed_order%n_sectors
          do sector=1,all_orders%n_sectors
             if (all(fixed_order%sector_powers(:,fixed_sector).eq.&
                  all_orders%sector_powers(:,sector))) exit
          enddo
          if (sector.gt.all_orders%n_sectors) then
             write (*,*) trim(label),' fixed order contains an unknown sector:',&
                  fixed_order%sector_powers(:,fixed_sector)
             stop 1
          endif
          fixed_components(:,sector)=fixed_order%amps_by_order(:,fixed_sector)
          fixed_present(:,sector)=fixed_order%sector_present(:,fixed_sector)
       enddo
       fixed(flow,:)=sum(fixed_components,dim=1)
       call accumulate_generator_lc_slices(all_orders,fixed_components,&
            fixed_present,part(:,1),generator_slices)
       deallocate(fixed_components,fixed_present)
       deallocate(fixed_order)
    enddo

    do sector=1,all_orders%n_sectors
       do flow=1,all_orders%n_amps
          scale=max(1d-30,abs(all_orders%amps_by_order(flow,sector)),&
               abs(fixed(flow,sector)))
          if (abs(all_orders%amps_by_order(flow,sector)-fixed(flow,sector)).gt.&
               relative_tolerance*scale) then
             write (*,*) trim(label),' imode=1 coefficient disagrees with imode=2:',&
                  flow,all_orders%sector_powers(:,sector),&
                  all_orders%amps_by_order(flow,sector),fixed(flow,sector)
             fixed_order_map_mismatch=.true.
          endif
       enddo
    enddo
    call check_generator_lc_slices(all_orders,generator_slices,label)
    deallocate(generator_slices)
  end subroutine check_fixed_order_sector_map

  subroutine check_generator_lc_slices(all_orders,generator_slices,label)
    implicit none
    type(amplitude_QCD),intent(in) :: all_orders
    real(kind=dp),dimension(0:,0:),intent(in) :: generator_slices
    character(len=*),intent(in) :: label
    integer :: as2,aew2,left_sector,right_sector
    real(kind=dp) :: fixed_slice,physical_slice,scale
    logical :: order_exists

    do as2=0,2*maxval(all_orders%sector_powers(1,:))
       do aew2=0,2*maxval(all_orders%sector_powers(2,:))
          order_exists=.false.
          do left_sector=1,all_orders%n_sectors
             do right_sector=left_sector,all_orders%n_sectors
                if (all_orders%sector_powers(1,left_sector)+&
                     all_orders%sector_powers(1,right_sector).ne.as2) cycle
                if (all_orders%sector_powers(2,left_sector)+&
                     all_orders%sector_powers(2,right_sector).ne.aew2) cycle
                order_exists=.true.
             enddo
          enddo
          if (.not.order_exists) cycle
          fixed_slice=generator_slices(as2,aew2)
          physical_slice=exact_squared_slice(all_orders,as2,aew2,1)
          scale=max(1d-30,abs(fixed_slice),abs(physical_slice))
          if (abs(fixed_slice-physical_slice).gt.relative_tolerance*scale) then
             write (*,'(a,1x,a,2i4,2es24.16)') trim(label),&
                  'generator/imode2 LC slice mismatch:',as2,aew2,&
                  fixed_slice,physical_slice
             generator_lc_mismatch=.true.
          endif
       enddo
    enddo
  end subroutine check_generator_lc_slices

  subroutine accumulate_generator_lc_slices(all_orders,fixed,fixed_present,&
       particles,generator_slices)
    implicit none
    type(amplitude_QCD),intent(in) :: all_orders
    complex(kind=dp),dimension(:,:),intent(in) :: fixed
    logical,dimension(:,:),intent(in) :: fixed_present
    integer,dimension(:),intent(in) :: particles
    real(kind=dp),dimension(0:,0:),intent(inout) :: generator_slices
    integer :: amplitude,left_sector,right_sector,colour_factor,as2,aew2
    complex(kind=dp) :: left_amp,right_amp

    colour_factor=generator_leading_colour_factor(particles)
    do amplitude=1,size(fixed,1)
       do left_sector=1,all_orders%n_sectors
          if (.not.fixed_present(amplitude,left_sector)) cycle
          left_amp=fixed(amplitude,left_sector)*sector_coupling(all_orders,left_sector)
          do right_sector=left_sector,all_orders%n_sectors
             if (.not.fixed_present(amplitude,right_sector)) cycle
             as2=all_orders%sector_powers(1,left_sector)+&
                  all_orders%sector_powers(1,right_sector)
             aew2=all_orders%sector_powers(2,left_sector)+&
                  all_orders%sector_powers(2,right_sector)
             right_amp=fixed(amplitude,right_sector)*&
                  sector_coupling(all_orders,right_sector)
             if (right_sector.eq.left_sector) then
                generator_slices(as2,aew2)=generator_slices(as2,aew2)+colour_factor*&
                     dble(left_amp*dconjg(right_amp))
             else
                generator_slices(as2,aew2)=generator_slices(as2,aew2)+2d0*colour_factor*&
                     dble(left_amp*dconjg(right_amp))
             endif
          enddo
       enddo
    enddo
  end subroutine accumulate_generator_lc_slices

  integer function generator_leading_colour_factor(particles)
    implicit none
    integer,dimension(:),intent(in) :: particles
    integer :: particle
    real(kind=dp) :: exponent

    exponent=0d0
    do particle=1,size(particles)
       if (particles(particle).eq.21) then
          exponent=exponent+1d0
       elseif (abs(particles(particle)).ge.1 .and.&
            abs(particles(particle)).le.6) then
          exponent=exponent+0.5d0
       endif
    enddo
    if (abs(exponent-dble(nint(exponent))).gt.1d-12) then
       write (*,*) 'Non-integer generator leading-colour exponent:',exponent
       stop 1
    endif
    generator_leading_colour_factor=3**nint(exponent)
  end function generator_leading_colour_factor

  subroutine check_sector_layout(amp,expected_powers,label)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer,dimension(:,:),intent(in) :: expected_powers
    character(len=*),intent(in) :: label
    integer :: sector,physical_amp
    real(kind=dp) :: value_scale

    if (.not.allocated(amp%sector_powers) .or.&
         .not.allocated(amp%sector_present) .or.&
         .not.allocated(amp%amps_by_order)) then
       write (*,*) trim(label),' did not allocate coupling-sector data'
       stop 1
    endif
    if (amp%n_sectors.ne.size(expected_powers,2)) then
       write (*,*) trim(label),' has the wrong sector count:',amp%n_sectors,&
            size(expected_powers,2)
       stop 1
    endif
    if (size(amp%sector_powers,1).ne.2 .or.&
         size(amp%sector_powers,2).lt.amp%n_sectors .or.&
         size(amp%sector_present,1).ne.amp%n_amps .or.&
         size(amp%sector_present,2).ne.amp%n_sectors .or.&
         size(amp%amps_by_order,1).ne.amp%n_amps .or.&
         size(amp%amps_by_order,2).ne.amp%n_sectors) then
       write (*,*) trim(label),' has inconsistent coupling-sector dimensions'
       write (*,*) 'n_amps, n_sectors:',amp%n_amps,amp%n_sectors
       write (*,*) 'sector_powers:',shape(amp%sector_powers)
       write (*,*) 'sector_present:',shape(amp%sector_present)
       write (*,*) 'amps_by_order:',shape(amp%amps_by_order)
       stop 1
    endif
    if (any(amp%sector_powers(:,1:amp%n_sectors).ne.expected_powers)) then
       write (*,*) trim(label),' has unexpected amplitude-sector powers:'
       do sector=1,amp%n_sectors
          write (*,*) sector,amp%sector_powers(:,sector)
       enddo
       stop 1
    endif
    do sector=1,amp%n_sectors
       if (.not.any(amp%sector_present(:,sector))) then
          write (*,*) trim(label),' has an empty advertised sector:',sector
          stop 1
       endif
       value_scale=max(1d-30,maxval(abs(amp%amps_by_order(:,sector))))
       do physical_amp=1,amp%n_amps
          if (amp%sector_present(physical_amp,sector)) cycle
          if (abs(amp%amps_by_order(physical_amp,sector)).gt.1d-14*value_scale) then
             write (*,*) trim(label),' populated an absent sector component:',&
                  physical_amp,sector,amp%amps_by_order(physical_amp,sector)
             stop 1
          endif
       enddo
    enddo
  end subroutine check_sector_layout

  real(kind=dp) function exact_squared_slice(amp,as2,aew2,accuracy)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer,intent(in) :: as2,aew2
    integer,intent(in),optional :: accuracy
    integer :: left_sector,right_sector,iacc
    complex(kind=dp),dimension(amp%n_amps) :: left,right,combined
    real(kind=dp) :: left_square,right_square

    exact_squared_slice=0d0
    iacc=3
    if (present(accuracy)) iacc=accuracy
    do left_sector=1,amp%n_sectors
       left=scaled_sector(amp,left_sector)
       left_square=colour_quadratic(amp,left,iacc)
       do right_sector=left_sector,amp%n_sectors
          if (amp%sector_powers(1,left_sector)+&
               amp%sector_powers(1,right_sector).ne.as2) cycle
          if (amp%sector_powers(2,left_sector)+&
               amp%sector_powers(2,right_sector).ne.aew2) cycle
          right=scaled_sector(amp,right_sector)
          if (right_sector.eq.left_sector) then
             exact_squared_slice=exact_squared_slice+left_square
          else
             right_square=colour_quadratic(amp,right,iacc)
             combined=left+right
             exact_squared_slice=exact_squared_slice+colour_quadratic(amp,combined,iacc)-&
                  left_square-right_square
          endif
       enddo
    enddo
  end function exact_squared_slice

  real(kind=dp) function sector_coupling(amp,sector)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer,intent(in) :: sector

    sector_coupling=sqrt(4d0*pi*alpha_s)**amp%sector_powers(1,sector)*&
         sqrt(8d0*pi*alpha_ew)**amp%sector_powers(2,sector)
  end function sector_coupling

  real(kind=dp) function coherent_sector_square(amp)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer :: sector
    complex(kind=dp),dimension(amp%n_amps) :: coherent

    coherent=(0d0,0d0)
    do sector=1,amp%n_sectors
       coherent=coherent+scaled_sector(amp,sector)
    enddo
    coherent_sector_square=colour_quadratic(amp,coherent)
  end function coherent_sector_square

  function scaled_sector(amp,sector) result(values)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer,intent(in) :: sector
    complex(kind=dp),dimension(amp%n_amps) :: values
    real(kind=dp) :: factor

    factor=sector_coupling(amp,sector)
    values=amp%amps_by_order(:,sector)*factor
  end function scaled_sector

  real(kind=dp) function colour_quadratic(amp,values,accuracy)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    complex(kind=dp),dimension(amp%n_amps),intent(in) :: values
    integer,intent(in),optional :: accuracy
    integer :: row,value,entry,colour,offset,iacc
    complex(kind=dp) :: weighted,sum_for_factor

    colour_quadratic=0d0
    iacc=3
    if (present(accuracy)) iacc=accuracy
    offset=amp%iproc_start(amp%nprocs)-1
    do row=1,amp%nColOrd
       weighted=(0d0,0d0)
       do value=1,amp%n_col_vals(iacc)
          sum_for_factor=(0d0,0d0)
          do entry=amp%row_index(row-1,value,iacc)+1,&
               amp%row_index(row,value,iacc)
             colour=amp%col_index(amp%i_col_i(value,iacc)+entry)
             sum_for_factor=sum_for_factor+values(offset+colour)
          enddo
          weighted=weighted+sum_for_factor*amp%diff_col_vals(value,iacc)
       enddo
       colour_quadratic=colour_quadratic+&
            dble(weighted*conjg(values(offset+row)))
    enddo
  end function colour_quadratic

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

end program coupling_order_sector_regression
