program find_unique_coupling_sector_regression
  use read_process_file, only: find_unique
  use handling_processes, only: phase_space_order_group
  use common, only: phys_model
  implicit none

  integer,parameter :: dp=kind(1d0),nevent=2,nslots=4,nsectors=2
  type(phase_space_order_group) :: group
  complex(kind=dp),dimension(nevent,nslots,nsectors) :: amplitudes
  integer,dimension(2) :: unique_map
  real(kind=dp),dimension(2) :: unique_value

  call phys_model%init_part()
  call initialise_group(group)

  call fill_proportional(amplitudes,(2d0,0d0))
  call find_unique(group,nevent,amplitudes,unique_map,unique_value)
  call require_result(unique_map,unique_value,1,4d0,&
       'real amplitude proportionality')

  call fill_proportional(amplitudes,(0d0,1d0))
  call find_unique(group,nevent,amplitudes,unique_map,unique_value)
  call require_result(unique_map,unique_value,1,1d0,&
       'complex-phase amplitude proportionality')

  call fill_base(amplitudes)
  amplitudes(:,3:4,1)=amplitudes(:,1:2,2)
  amplitudes(:,3:4,2)=amplitudes(:,1:2,1)
  call find_unique(group,nevent,amplitudes,unique_map,unique_value)
  call require_result(unique_map,unique_value,-1,1d0,&
       'coupling-sector permutation')

  call fill_base(amplitudes)
  amplitudes(:,3:4,:)=amplitudes(:,1:2,:)
  group%amps(1)%sector_present=.true.
  group%amps(1)%sector_present(2,2)=.false.
  amplitudes(:,2,2)=(0d0,0d0)
  amplitudes(:,4,2)=(0d0,0d0)
  call find_unique(group,nevent,amplitudes,unique_map,unique_value)
  call require_result(unique_map,unique_value,-1,1d0,&
       'absent-sector mismatch')
  group%amps(1)%sector_present=.true.

  call fill_base(amplitudes)
  amplitudes(:,1:2,2)=1d-8*amplitudes(:,1:2,1)
  amplitudes(:,3:4,1)=amplitudes(:,1:2,1)
  amplitudes(:,3:4,2)=2d-8*amplitudes(:,1:2,1)
  call find_unique(group,nevent,amplitudes,unique_map,unique_value)
  call require_result(unique_map,unique_value,-1,1d0,&
       'subleading-sector scale hierarchy')

  write (*,'(a)') 'Coupling-sector unique-process regression passed'

contains

  subroutine initialise_group(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl

    pgl%next=4
    pgl%nproc=2
    allocate(pgl%processes(pgl%next,pgl%nproc))
    pgl%processes(:,1)=[2,1,2,1]
    pgl%processes(:,2)=pgl%processes(:,1)
    allocate(pgl%amps(1))
    pgl%amps(1)%n_amps=nslots
    pgl%amps(1)%n_sectors=nsectors
    allocate(pgl%amps(1)%iproc_start(3))
    pgl%amps(1)%iproc_start=[1,3,5]
    allocate(pgl%amps(1)%sector_present(nslots,nsectors))
    pgl%amps(1)%sector_present=.true.
    allocate(pgl%amps(1)%sector_powers(2,nsectors))
    pgl%amps(1)%sector_powers(:,1)=[2,0]
    pgl%amps(1)%sector_powers(:,2)=[0,2]
  end subroutine initialise_group

  subroutine fill_base(values)
    implicit none
    complex(kind=dp),dimension(nevent,nslots,nsectors),intent(out) :: values

    values=(0d0,0d0)
    values(1,1,1)=(1d0,0.2d0)
    values(1,2,1)=(-0.4d0,0.5d0)
    values(2,1,1)=(0.6d0,-0.8d0)
    values(2,2,1)=(1.1d0,0.3d0)
    values(1,1,2)=(0.3d0,-0.1d0)
    values(1,2,2)=(0.7d0,0.2d0)
    values(2,1,2)=(-0.2d0,0.9d0)
    values(2,2,2)=(0.5d0,-0.6d0)
  end subroutine fill_base

  subroutine fill_proportional(values,factor)
    implicit none
    complex(kind=dp),dimension(nevent,nslots,nsectors),intent(out) :: values
    complex(kind=dp),intent(in) :: factor

    call fill_base(values)
    values(:,3:4,:)=factor*values(:,1:2,:)
    group%amps(1)%sector_present=.true.
  end subroutine fill_proportional

  subroutine require_result(mapping,factors,expected_map,expected_factor,label)
    implicit none
    integer,dimension(2),intent(in) :: mapping
    real(kind=dp),dimension(2),intent(in) :: factors
    integer,intent(in) :: expected_map
    real(kind=dp),intent(in) :: expected_factor
    character(len=*),intent(in) :: label

    if (mapping(1).ne.-1 .or. abs(factors(1)-1d0).gt.1d-12 .or.&
         mapping(2).ne.expected_map .or.&
         abs(factors(2)-expected_factor).gt.1d-12) then
       write (*,*) trim(label),' produced the wrong unique map:',&
            mapping,factors
       stop 1
    endif
  end subroutine require_result

end program find_unique_coupling_sector_regression
