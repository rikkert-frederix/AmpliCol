program process_color_order_regression
  use common, only: phys_model
  use handling_processes, only: next,phase_space_order_group,pgl
  use read_process_file, only: move_colour_singlet_in_order,find_unique
  use integrated_dipoles, only: integrated_history,history_matches_born,&
       integrated_endpoint,integrated_beam
  use pdf_wrap, only: pdf_table_flavour,pdf_flavour_is_supported
  implicit none
  integer :: process(5),order(5)
  integer :: minimum_integer,i
  type(phase_space_order_group) :: group
  real(kind=8) :: samples(4,5),map_value(5)
  integer :: map(5)

  open(unit=99,status='scratch',action='readwrite')
  call phys_model%init_part()
  next=5

  process=[1,-1,21,22,23]
  order=[1,3,2,4,5]
  call move_colour_singlet_in_order(process,order)
  call assert_equal(order,[1,3,4,5,2],&
       'quark order does not close on its final coloured leg')

  process=[21,22,21,23,21]
  order=[2,1,4,3,5]
  call move_colour_singlet_in_order(process,order)
  call assert_equal(order,[1,3,2,4,5],&
       'gluon-only colour string mishandles singlets')

  process=[1,-1,21,2,-2]
  order=[5,2,3,1,4]
  call move_colour_singlet_in_order(process,order)
  call assert_equal(order,[5,2,3,1,4],'all-coloured order changed')

  process=[22,23,24,-24,25]
  order=[5,4,3,2,1]
  call move_colour_singlet_in_order(process,order)
  call assert_equal(order,[5,4,3,2,1],'all-singlet order changed')

  group%next=4
  group%nproc=5
  allocate(group%processes(4,5))
  group%processes=spread([2,-2,21,21],2,5)
  samples(:,1)=[1d0,2d0,4d0,8d0]
  samples(:,2)=7d0*samples(:,1)
  samples(:,3)=[1d0,2d0,4d0,8.00001d0]
  samples(:,4)=0d0
  samples(:,5)=2d0*samples(:,3)
  call find_unique(group,4,samples,map,map_value)
  call assert_equal(map,[-1,1,-1,0,3],&
       'hashed process reduction changed representative mapping')
  if (abs(map_value(2)-7d0).gt.1d-13 .or. &
       abs(map_value(5)-2d0).gt.1d-13 .or. map_value(4).ne.0d0) then
     write(*,*) 'FAIL: hashed process reduction changed mapping factors',map_value
     stop 1
  endif

  call check_integrated_history_colour_order()
  minimum_integer=-huge(0)-1
  if (pdf_flavour_is_supported(minimum_integer) .or. &
       pdf_table_flavour([(0d0,i=-6,7)],minimum_integer).ne.0d0) then
     write(*,*) 'FAIL: minimum integer was accepted as a PDF flavour'
     stop 1
  endif

  close(99)
  write(*,'(a)') 'process colour-order regression: PASS'

contains

  subroutine check_integrated_history_colour_order()
    type(integrated_history) :: history
    real(kind=8) :: born_copy(1),momentum(0:3,4),coeff(-2:0),pterm,kterm
    integer :: status

    allocate(pgl(1))
    pgl(1)%next=4
    pgl(1)%nproc=1
    allocate(pgl(1)%iden_iproc(1),pgl(1)%iden_processes(4,1,1),&
         pgl(1)%color_orders(4,1))
    pgl(1)%iden_iproc=1
    pgl(1)%iden_processes(:,1,1)=[2,-2,21,22]
    pgl(1)%color_orders(:,1)=[1,3,4,2]
    history%born_flavours=[2,-2,21,22]
    history%born_colour_order=[1,3,4,2]
    if (.not.history_matches_born(history,1,1)) then
       write(*,*) 'FAIL: valid integrated-history colour order was rejected'
       stop 1
    endif
    history%born_colour_order=[1,3,5,2]
    if (history_matches_born(history,1,1)) then
       write(*,*) 'FAIL: out-of-range integrated-history colour label was accepted'
       stop 1
    endif
    history%born_colour_order=[1,3,3,2]
    if (history_matches_born(history,1,1)) then
       write(*,*) 'FAIL: duplicate integrated-history colour label was accepted'
       stop 1
    endif
    if (history_matches_born(history,0,1) .or. &
         history_matches_born(history,1,2)) then
       write(*,*) 'FAIL: invalid integrated-history group or process was accepted'
       stop 1
    endif
    born_copy=1d0
    momentum=0d0
    call integrated_endpoint(0,1,born_copy,momentum,1d0,0.1d0,coeff,status=status)
    if (status.ne.-4 .or. any(coeff.ne.0d0)) then
       write(*,*) 'FAIL: invalid integrated endpoint group was not rejected safely',status
       stop 1
    endif
    call integrated_beam(0,1,1,0.5d0,born_copy,[0.1d0,0.1d0],1d0,1d0,0.1d0,&
         pterm,kterm,status=status)
    if (status.ne.-4 .or. pterm.ne.0d0 .or. kterm.ne.0d0) then
       write(*,*) 'FAIL: invalid integrated beam group was not rejected safely',status
       stop 1
    endif
    deallocate(pgl)
  end subroutine check_integrated_history_colour_order

  subroutine assert_equal(actual,expected,label)
    integer,intent(in) :: actual(:),expected(:)
    character(len=*),intent(in) :: label
    if (size(actual).ne.size(expected)) then
       write(*,*) 'FAIL: ',trim(label),' (shape mismatch)'
       stop 1
    endif
    if (any(actual.ne.expected)) then
       write(*,*) 'FAIL: ',trim(label)
       write(*,*) 'actual:  ',actual
       write(*,*) 'expected:',expected
       stop 1
    endif
  end subroutine assert_equal

end program process_color_order_regression
