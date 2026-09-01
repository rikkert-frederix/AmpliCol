module mg_checks
  use common
  use amplitude_QCD_mod
  use handling_processes
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: max_madgraph_reference_values=1000000
  real(kind=8),dimension(:,:),allocatable :: p_read
  real(kind=8),parameter  :: alpha_check=0.118d0
  real(kind=8), parameter :: pi=3.14159265358979323846d0
  integer :: me_points=0

contains

  logical function checked_finite_product(first,second,value) result(valid)
    implicit none
    real(kind=8),intent(in) :: first,second
    real(kind=8),intent(out) :: value
    valid=.false.
    value=0d0
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    if (first.eq.0d0 .or. second.eq.0d0) then
       valid=.true.
       return
    endif
    if (abs(second).gt.1d0) then
       if (abs(first).gt.huge(1d0)/abs(second)) return
    elseif (abs(first).gt.1d0) then
       if (abs(second).gt.huge(1d0)/abs(first)) return
    endif
    value=first*second
    valid=ieee_is_finite(value)
    if (.not.valid) value=0d0
  end function checked_finite_product

  logical function checked_positive_integer_power(base,power,value) result(valid)
    implicit none
    real(kind=8),intent(in) :: base
    integer,intent(in) :: power
    real(kind=8),intent(out) :: value
    real(kind=8) :: updated_value
    integer :: i
    valid=.false.
    value=0d0
    if (.not.ieee_is_finite(base) .or. base.le.0d0 .or. power.lt.0) return
    value=1d0
    do i=1,power
       if (.not.checked_finite_product(value,base,updated_value)) then
          value=0d0
          return
       endif
       value=updated_value
       if (value.le.0d0) then
          value=0d0
          return
       endif
    enddo
    valid=.true.
  end function checked_positive_integer_power

  subroutine get_madgraph_results(n,ichan,iint,me,nlines)
    implicit none
    integer,intent(in) :: n,ichan,iint
    real(kind=8),dimension(:),allocatable,intent(out) :: me
    integer,intent(out) :: nlines
    real(kind=8),dimension(0:3) :: dum
    character(len=256) :: filename,io_message,allocation_message
    character(len=1024) :: trailing_line
    integer :: i,iunit,io_status,allocation_status

    if (n.lt.1 .or. n.gt.max_amplitude_external_particles .or. &
         ichan.lt.1 .or. iint.lt.1) then
       write (*,*) 'Invalid external multiplicity in MadGraph result reader:',n
       stop 1
    endif
    write(filename, '("Utilities/ME_checks/momenta_", I0, "_",I0,".txt")') ichan,iint
    open(newunit=iunit,file=trim(filename),status='old',action='read',&
         iostat=io_status,iomsg=io_message)
    if (io_status.ne.0) then
       write (*,*) 'Could not open MadGraph reference file ',trim(filename),': ',&
            trim(io_message)
       stop 1
    endif
    do i=1,n
       read(iunit,*,iostat=io_status,iomsg=io_message) dum
       if (io_status.ne.0) then
          write (*,*) 'Could not read MadGraph reference momentum ',i,': ',&
               trim(io_message)
          stop 1
       endif
       if (.not.all(ieee_is_finite(dum))) then
          write (*,*) 'MadGraph reference file contains a non-finite momentum:',i
          stop 1
       endif
    enddo
    nlines=0
    read(iunit,*,iostat=io_status,iomsg=io_message) nlines
    if (io_status.ne.0) then
       write (*,*) 'Could not read MadGraph reference count: ',trim(io_message)
       stop 1
    endif
    if (nlines.lt.1 .or. nlines.gt.max_madgraph_reference_values) then
       write (*,*) 'Invalid MadGraph reference count:',nlines
       stop 1
    endif
    allocate(me(nlines),stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate MadGraph reference values: ',&
            trim(allocation_message)
       stop 1
    endif
    do i=1,nlines
       read(iunit,*,iostat=io_status,iomsg=io_message) me(i)
       if (io_status.ne.0) then
          write (*,*) 'Could not read MadGraph reference value ',i,': ',&
               trim(io_message)
          stop 1
       endif
    enddo
    do
       trailing_line=''
       read(iunit,'(A)',iostat=io_status,iomsg=io_message) trailing_line
       if (io_status.lt.0) exit
       if (io_status.gt.0) then
          write (*,*) 'Could not validate the end of the MadGraph reference file: ',&
               trim(io_message)
          stop 1
       endif
       if (len_trim(trailing_line).gt.0) then
          write (*,*) 'MadGraph reference file contains trailing data'
          stop 1
       endif
    enddo
    close(iunit,iostat=io_status,iomsg=io_message)
    if (io_status.ne.0) then
       write (*,*) 'Could not close MadGraph reference file: ',trim(io_message)
       stop 1
    endif
    if (.not.all(ieee_is_finite(me))) then
       write (*,*) 'MadGraph reference file contains invalid matrix elements'
       stop 1
    endif
    if (any(me.lt.0d0)) then
       write (*,*) 'MadGraph reference file contains invalid matrix elements'
       stop 1
    endif
  end subroutine get_madgraph_results

  subroutine read_in_momenta(n,igroup,iamp,p_in)
    implicit none
    integer,intent(in) :: n,iamp,igroup
    real(kind=8),dimension(n,0:3),intent(out) :: p_in
    character(len=256) :: filename,io_message
    integer :: i,iunit,io_status

    if (n.lt.1 .or. n.gt.max_amplitude_external_particles .or. &
         igroup.lt.1 .or. iamp.lt.1) then
       write (*,*) 'Invalid external multiplicity in MadGraph momentum reader:',n
       stop 1
    endif
    write(filename, '("Utilities/ME_checks/momenta_", I0, "_",I0,".txt")') igroup,iamp
    open(newunit=iunit,file=trim(filename),status='old',action='read',&
         iostat=io_status,iomsg=io_message)
    if (io_status.ne.0) then
       write (*,*) 'Could not open MadGraph momentum file ',trim(filename),': ',&
            trim(io_message)
       stop 1
    endif
    do i=1,n
       read(iunit,*,iostat=io_status,iomsg=io_message) p_in(i,:)
       if (io_status.ne.0) then
          write (*,*) 'Could not read MadGraph momentum ',i,': ',trim(io_message)
          stop 1
       endif
    enddo
    close(iunit,iostat=io_status,iomsg=io_message)
    if (io_status.ne.0) then
       write (*,*) 'Could not close MadGraph momentum file: ',trim(io_message)
       stop 1
    endif
    if (.not.all(ieee_is_finite(p_in))) then
       write (*,*) 'MadGraph momentum file contains non-finite values'
       stop 1
    endif
  end subroutine read_in_momenta

  subroutine run_madgraph_check(n,igroup,iamp,list)
    implicit none
    integer,intent(in) :: n,igroup,iamp
    integer,dimension(n),intent(in) :: list
    integer :: exit_status,command_status,io_status,i
    character(len=1024) :: command
    character(len=256) :: command_message

    if (n.lt.1 .or. n.gt.max_amplitude_external_particles .or. &
         igroup.lt.1 .or. iamp.lt.1) then
       write (*,*) 'Invalid external multiplicity for MadGraph command:',n
       stop 1
    endif
    command=''
    write(command,'(a,i0,1x,i0)',iostat=io_status) &
         './Utilities/ME_checks/mg_wrapper.sh ',igroup,iamp
    if (io_status.ne.0) then
       write (*,*) 'Could not construct MadGraph command'
       stop 1
    endif
    do i=1,n
       if (len_trim(command)+32.gt.len(command)) then
          write (*,*) 'MadGraph command exceeds supported length'
          stop 1
       endif
       write(command(len_trim(command)+1:),'(1x,i0)',iostat=io_status) list(i)
       if (io_status.ne.0) then
          write (*,*) 'Could not append particle to MadGraph command:',i,list(i)
          stop 1
       endif
    enddo
    command_message=''
    exit_status=-1
    command_status=-1
    call execute_command_line(trim(command),wait=.true.,exitstat=exit_status,&
         cmdstat=command_status,cmdmsg=command_message)
    if (command_status.ne.0 .or. exit_status.ne.0) then
       write (*,*) 'MadGraph command failed:',command_status,exit_status,&
            trim(command_message)
       stop 1
    endif
    write (*,*) 'MadGraph script finished successfully for group/amplitude:',&
         igroup,iamp
  end subroutine run_madgraph_check

  subroutine perform_check(iint,ichan)
    implicit none
    integer,intent(in) :: iint,ichan
    integer :: i,nreference,amp_index,singlet_count
    real(kind=8),dimension(:),allocatable :: mg_check
    real(kind=8) :: me_code,normalization,comparison_scale,qcd_base,ew_base,&
         qcd_factor,ew_factor,normalization_numerator
    logical :: match

    if (ichan.lt.1 .or. ichan.gt.ngroups) then
       write (*,*) 'Invalid process state for MadGraph comparison:',ichan,iint
       stop 1
    endif
    if (.not.allocated(pgl)) then
       write (*,*) 'Missing process groups for MadGraph comparison'
       stop 1
    endif
    if (ngroups.lt.1 .or. size(pgl).lt.ngroups .or. ichan.gt.size(pgl)) then
       write (*,*) 'Inconsistent process-group state for MadGraph comparison:',&
            ngroups,size(pgl),ichan
       stop 1
    endif
    if (pgl(ichan)%next.lt.3 .or. &
         pgl(ichan)%next.gt.max_amplitude_external_particles .or. &
         pgl(ichan)%nproc.lt.1) then
       write (*,*) 'Invalid process dimensions for MadGraph comparison:',&
            pgl(ichan)%next,pgl(ichan)%nproc
       stop 1
    endif
    if (iint.lt.1 .or. iint.gt.pgl(ichan)%nproc) then
       write (*,*) 'Invalid process index for MadGraph comparison:',ichan,iint
       stop 1
    endif
    if (.not.allocated(pgl(ichan)%amp2) .or. &
         .not.allocated(pgl(ichan)%amps) .or. .not.allocated(pgl(ichan)%iden)) then
       write (*,*) 'Incomplete process state for MadGraph comparison:',ichan,iint
       stop 1
    endif
    if (size(pgl(ichan)%amp2).eq.1) then
       amp_index=1
    elseif (size(pgl(ichan)%amp2).ge.iint) then
       amp_index=iint
    else
       write (*,*) 'Matrix-element result array is too short for MadGraph comparison'
       stop 1
    endif
    if (size(pgl(ichan)%amps).eq.1) then
       if (.not.allocated(pgl(ichan)%amps(1)%n_sing)) then
          write (*,*) 'Amplitude singlet metadata is missing for MadGraph comparison'
          stop 1
       endif
       if (size(pgl(ichan)%amps(1)%n_sing).lt.iint) then
          write (*,*) 'Amplitude singlet metadata is incomplete for MadGraph comparison'
          stop 1
       endif
       singlet_count=pgl(ichan)%amps(1)%n_sing(iint)
    elseif (size(pgl(ichan)%amps).ge.iint) then
       if (.not.allocated(pgl(ichan)%amps(iint)%n_sing)) then
          write (*,*) 'Amplitude singlet metadata is missing for MadGraph comparison'
          stop 1
       endif
       if (size(pgl(ichan)%amps(iint)%n_sing).lt.1) then
          write (*,*) 'Amplitude singlet metadata is incomplete for MadGraph comparison'
          stop 1
       endif
       singlet_count=pgl(ichan)%amps(iint)%n_sing(1)
    else
       write (*,*) 'Amplitude array is too short for MadGraph comparison'
       stop 1
    endif
    if (size(pgl(ichan)%iden).lt.iint) then
       write (*,*) 'Process normalization array is too short for MadGraph comparison'
       stop 1
    endif
    if (singlet_count.lt.0 .or. singlet_count.gt.pgl(ichan)%next-2 .or. &
         pgl(ichan)%iden(iint).le.0_8 .or. &
         .not.ieee_is_finite(pgl(ichan)%amp2(amp_index))) then
       write (*,*) 'Invalid normalization metadata for MadGraph comparison:',&
            singlet_count,pgl(ichan)%iden(iint)
       stop 1
    endif
    if (pgl(ichan)%amp2(amp_index).lt.0d0) then
       write (*,*) 'Invalid normalization metadata for MadGraph comparison:',&
            singlet_count,pgl(ichan)%iden(iint)
       stop 1
    endif
    if (.not.ieee_is_finite(alphaEW)) then
       write (*,*) 'Invalid electromagnetic coupling for MadGraph comparison:',alphaEW
       stop 1
    endif
    if (alphaEW.le.0d0 .or. alphaEW.gt.huge(alphaEW)/(8d0*pi)) then
       write (*,*) 'Invalid electromagnetic coupling for MadGraph comparison:',alphaEW
       stop 1
    endif
    qcd_base=4d0*pi*alpha_check
    ew_base=8d0*pi*alphaEW
    if (.not.checked_positive_integer_power(qcd_base,&
         pgl(ichan)%next-2-singlet_count,qcd_factor) .or. &
         .not.checked_positive_integer_power(ew_base,singlet_count,ew_factor)) then
       write (*,*) 'MadGraph comparison normalization cannot be represented'
       stop 1
    endif
    if (.not.checked_finite_product(qcd_factor,ew_factor,&
         normalization_numerator)) then
       write (*,*) 'MadGraph comparison normalization cannot be represented'
       stop 1
    endif
    normalization=normalization_numerator/dble(pgl(ichan)%iden(iint))
    if (.not.ieee_is_finite(normalization) .or. normalization.le.0d0) then
       write (*,*) 'MadGraph comparison normalization cannot be represented'
       stop 1
    endif
    if (.not.checked_finite_product(pgl(ichan)%amp2(amp_index),&
         normalization,me_code)) then
       write (*,*) 'Matrix element overflows MadGraph comparison normalization:',&
            pgl(ichan)%amp2(amp_index),normalization
       stop 1
    endif
    call get_madgraph_results(pgl(ichan)%next,ichan,iint,mg_check,nreference)
    match=.false.
    do i=1,nreference
       comparison_scale=max(abs(mg_check(i)),abs(me_code))
       if (comparison_scale.eq.0d0) then
          match=.true.
       elseif (abs(mg_check(i)/comparison_scale-&
            me_code/comparison_scale).le.1d-4) then
          match=.true.
       endif
       if (match) exit
    enddo
    if (.not.match) then
       write (*,*) 'ERROR: disagreement with MadGraph in ME-level check:',&
            ichan,iint,me_code
       write (*,*) 'MadGraph reference values:',mg_check
       stop 4
    endif
  end subroutine perform_check


end module mg_checks
