program matrix_element_regression
  use amplitude_QCD_mod
  use particles
  implicit none

  integer,parameter :: dp=kind(1d0)
  integer,parameter :: max_name=64,max_failures_to_print=20,max_errors=100
  real(kind=dp),parameter :: abs_tol=1d-12,rel_tol=5d-10

  type(physics_model) :: model
  character(len=256) :: mode,cases_file,golden_file

  call parse_arguments(mode,cases_file,golden_file)

  open(unit=99,file='/dev/null',status='unknown',action='write')
  call model%init_part(173d0,1.491500d0,91.188d0,2.441404d0,&
                       80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call model%init_vert()

  if (trim(mode).eq.'--write') then
     call write_goldens(trim(cases_file),trim(golden_file))
  elseif (trim(mode).eq.'--check') then
     call check_goldens(trim(cases_file),trim(golden_file))
  else
     call usage()
     stop 2
  endif

contains

  subroutine parse_arguments(mode,cases_file,golden_file)
    implicit none
    character(len=*),intent(out) :: mode,cases_file,golden_file
    if (command_argument_count().ne.3) then
       call usage()
       stop 2
    endif
    call get_command_argument(1,mode)
    call get_command_argument(2,cases_file)
    call get_command_argument(3,golden_file)
  end subroutine parse_arguments

  subroutine usage()
    implicit none
    write (*,'(a)') 'Usage: matrix_element_regression --check cases.dat golden.dat'
    write (*,'(a)') '   or: matrix_element_regression --write cases.dat golden.dat'
  end subroutine usage

  subroutine read_header(iunit,version,ncases)
    implicit none
    integer,intent(in) :: iunit
    integer,intent(out) :: version,ncases
    read(iunit,*) version,ncases
    if (version.ne.1) then
       write (*,*) 'Unsupported matrix-element fixture version:',version
       stop 2
    endif
  end subroutine read_header

  subroutine read_case(iunit,case_id,family,n,group_id,row_id,point,process,order)
    implicit none
    integer,intent(in) :: iunit
    integer,intent(out) :: case_id,n,group_id,row_id
    character(len=*),intent(out) :: family,point
    integer,dimension(:),allocatable,intent(out) :: process,order
    read(iunit,*) case_id,family,n,group_id,row_id,point
    allocate(process(n),order(n))
    read(iunit,*) process(1:n)
    read(iunit,*) order(1:n)
  end subroutine read_case

  subroutine write_goldens(cases_file,golden_file)
    implicit none
    character(len=*),intent(in) :: cases_file,golden_file
    integer :: cunit,gunit,version,ncases,icase
    integer :: case_id,n,group_id,row_id,n_amps
    integer,dimension(:),allocatable :: process,order
    integer,dimension(:,:),allocatable :: spins
    complex(kind=dp),dimension(:),allocatable :: amps
    character(len=max_name) :: family,point

    open(newunit=cunit,file=cases_file,status='old',action='read')
    open(newunit=gunit,file=golden_file,status='replace',action='write')
    call read_header(cunit,version,ncases)
    write(gunit,'(i0,1x,i0)') version,ncases

    do icase=1,ncases
       call read_case(cunit,case_id,family,n,group_id,row_id,point,process,order)
       call evaluate_case(n,process,order,point,n_amps,spins,amps)
       write(gunit,'(i0,1x,a,1x,i0,1x,i0,1x,i0,1x,i0)') &
            case_id,trim(family),n,group_id,row_id,n_amps
       call write_int_array(gunit,process)
       call write_int_array(gunit,order)
       call write_amplitude_rows(gunit,spins,amps)
       if (mod(icase,100).eq.0) write (*,*) 'Wrote matrix-element goldens:',icase,'/',ncases
       deallocate(process,order,spins,amps)
    enddo

    close(cunit)
    close(gunit)
    write (*,*) 'Wrote matrix-element golden file:',trim(golden_file)
  end subroutine write_goldens

  subroutine check_goldens(cases_file,golden_file)
    implicit none
    character(len=*),intent(in) :: cases_file,golden_file
    integer :: cunit,gunit,version,ncases,gversion,gncases,icase
    integer :: case_id,n,group_id,row_id,n_amps
    integer :: g_case_id,g_n,g_group_id,g_row_id,g_n_amps
    integer :: total_failures,total_checks,total_cases
    integer,dimension(:),allocatable :: process,order,g_process,g_order
    integer,dimension(:,:),allocatable :: spins,g_spins
    complex(kind=dp),dimension(:),allocatable :: amps,g_amps
    character(len=max_name) :: family,point,g_family

    open(newunit=cunit,file=cases_file,status='old',action='read')
    open(newunit=gunit,file=golden_file,status='old',action='read')
    call read_header(cunit,version,ncases)
    call read_header(gunit,gversion,gncases)
    if (gncases.ne.ncases) then
       write (*,*) 'Golden case count mismatch:',gncases,ncases
       stop 1
    endif

    total_failures=0
    total_checks=0
    total_cases=0
    do icase=1,ncases
       call read_case(cunit,case_id,family,n,group_id,row_id,point,process,order)
       call read_golden_case(gunit,g_case_id,g_family,g_n,g_group_id,g_row_id,&
            g_n_amps,g_process,g_order,g_spins,g_amps)

       if (.not.same_case_metadata(case_id,family,n,group_id,row_id,process,order,&
            g_case_id,g_family,g_n,g_group_id,g_row_id,g_process,g_order)) then
          write (*,*) 'Golden metadata mismatch at case:',case_id
          stop 1
       endif

       call evaluate_case(n,process,order,point,n_amps,spins,amps)
       call compare_case(case_id,family,group_id,row_id,n,n_amps,spins,amps,&
            g_n_amps,g_spins,g_amps,total_failures,total_checks)
       total_cases=total_cases+1
       if (mod(icase,100).eq.0) write (*,*) 'Checked matrix-element goldens:',icase,'/',ncases

       deallocate(process,order,g_process,g_order,spins,g_spins,amps,g_amps)
    enddo

    close(cunit)
    close(gunit)

    if (total_failures.gt.0) then
       write (*,*) 'Matrix-element regression failures:',total_failures,'of',total_checks,'amplitudes'
       stop 1
    endif
    write (*,*) 'Matrix-element regression passed:',total_cases,'cases and',total_checks,'helicity amplitudes'
  end subroutine check_goldens

  subroutine read_golden_case(iunit,case_id,family,n,group_id,row_id,n_amps,process,order,spins,amps)
    implicit none
    integer,intent(in) :: iunit
    integer,intent(out) :: case_id,n,group_id,row_id,n_amps
    character(len=*),intent(out) :: family
    integer,dimension(:),allocatable,intent(out) :: process,order
    integer,dimension(:,:),allocatable,intent(out) :: spins
    complex(kind=dp),dimension(:),allocatable,intent(out) :: amps
    integer :: i
    real(kind=dp) :: re,im

    read(iunit,*) case_id,family,n,group_id,row_id,n_amps
    allocate(process(n),order(n),spins(n,n_amps),amps(n_amps))
    read(iunit,*) process(1:n)
    read(iunit,*) order(1:n)
    do i=1,n_amps
       read(iunit,*) spins(1:n,i),re,im
       amps(i)=cmplx(re,im,kind=dp)
    enddo
  end subroutine read_golden_case

  logical function same_case_metadata(case_id,family,n,group_id,row_id,process,order,&
       g_case_id,g_family,g_n,g_group_id,g_row_id,g_process,g_order)
    implicit none
    integer,intent(in) :: case_id,n,group_id,row_id,g_case_id,g_n,g_group_id,g_row_id
    character(len=*),intent(in) :: family,g_family
    integer,dimension(:),intent(in) :: process,order,g_process,g_order
    same_case_metadata=.false.
    if (case_id.ne.g_case_id) return
    if (trim(family).ne.trim(g_family)) return
    if (n.ne.g_n) return
    if (group_id.ne.g_group_id .or. row_id.ne.g_row_id) return
    if (size(process).ne.size(g_process) .or. size(order).ne.size(g_order)) return
    if (any(process.ne.g_process) .or. any(order.ne.g_order)) return
    same_case_metadata=.true.
  end function same_case_metadata

  subroutine evaluate_case(n,process,order,point,n_amps,spins,amps)
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(in) :: process,order
    character(len=*),intent(in) :: point
    integer,intent(out) :: n_amps
    integer,dimension(:,:),allocatable,intent(out) :: spins
    complex(kind=dp),dimension(:),allocatable,intent(out) :: amps

    type(amplitude_QCD) :: amp
    integer,dimension(:,:),allocatable :: part,orders
    integer,dimension(:,:),allocatable :: spin
    integer,dimension(:),allocatable :: hel
    real(kind=dp),dimension(:,:),allocatable :: p
    integer :: i

    allocate(part(n,1),orders(n,1),hel(n))
    part(1:n,1)=process(1:n)
    orders(1:n,1)=order(1:n)
    hel=0

    call setup_spin(n,process,spin)
    call fill_reference_momenta(n,process,point,p)

    call amp%init(1,n,1,part,spin,orders,model)
    call amp%evaluate(n,p,hel,.false.,model)

    n_amps=amp%n_amps
    allocate(spins(n,n_amps),amps(n_amps))
    do i=1,n_amps
       spins(1:n,i)=amp%spins(1:n,1,i)
    enddo
    amps(1:n_amps)=amp%amps(1:n_amps)
  end subroutine evaluate_case

  subroutine setup_spin(n,process,spin)
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(in) :: process
    integer,dimension(:,:),allocatable,intent(out) :: spin
    integer :: i,nspin
    allocate(spin(0:3,n))
    spin=0
    do i=1,n
       nspin=model%get_spin(process(i))
       spin(0,i)=nspin
       if (nspin.eq.1) then
          spin(1,i)=0
       elseif (nspin.eq.2) then
          spin(1,i)=-1
          spin(2,i)=1
       elseif (nspin.eq.3) then
          spin(1,i)=-1
          spin(2,i)=0
          spin(3,i)=1
       else
          write (*,*) 'Unsupported spin state:',process(i),nspin
          stop 1
       endif
    enddo
  end subroutine setup_spin

  subroutine fill_reference_momenta(n,process,point,p)
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(in) :: process
    character(len=*),intent(in) :: point
    real(kind=dp),dimension(:,:),allocatable,intent(out) :: p
    integer :: i,nmassive,imassive

    allocate(p(0:3,n))
    if (trim(point).eq.'generic') then
       call fill_generic_momenta(n,process,p)
       return
    endif

    nmassive=0
    imassive=0
    do i=3,n
       if (model%get_mass(process(i)).gt.0d0) then
          nmassive=nmassive+1
          imassive=i
       endif
    enddo

    if (nmassive.eq.1) then
       if (trim(point).ne.'massive5' .or. n.ne.5) call unsupported_point(n,process,point)
       call fill_massive5(process,imassive,p)
    elseif (nmassive.eq.0 .and. n.eq.4) then
       if (trim(point).ne.'massless4') call unsupported_point(n,process,point)
       call fill_massless4(p)
    elseif (nmassive.eq.0 .and. n.eq.5) then
       if (trim(point).ne.'massless5') call unsupported_point(n,process,point)
       call fill_massless5(p)
    elseif (nmassive.eq.0 .and. n.eq.6) then
       if (trim(point).ne.'massless6') call unsupported_point(n,process,point)
       call fill_massless6(p)
    else
       call unsupported_point(n,process,point)
    endif
  end subroutine fill_reference_momenta

  subroutine fill_generic_momenta(n,process,p)
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(in) :: process
    real(kind=dp),dimension(0:3,n),intent(out) :: p
    integer :: nfinal,pair_start,i,j,ipair
    real(kind=dp) :: total_energy,q,mass_i,mass_j
    real(kind=dp),dimension(3,5) :: directions

    nfinal=n-2
    if (nfinal.lt.2) then
       write (*,*) 'Generic matrix-element point needs at least two final particles'
       stop 1
    endif

    directions(:,1)=[0.3d0,0.4d0,sqrt(0.75d0)]
    directions(:,2)=[0.6d0,-0.2d0,sqrt(0.60d0)]
    directions(:,3)=[-0.1d0,0.7d0,sqrt(0.50d0)]
    directions(:,4)=[sqrt(0.30d0),sqrt(0.20d0),-sqrt(0.50d0)]
    directions(:,5)=[-sqrt(0.20d0),sqrt(0.55d0),-0.5d0]

    p=0d0
    total_energy=0d0

    if (mod(nfinal,2).eq.1) then
       q=120d0+11d0*dble(n)
       p(1:3,3)=q*[1d0,0d0,0d0]
       p(1:3,4)=q*[-0.5d0,sqrt(0.75d0),0d0]
       p(1:3,5)=q*[-0.5d0,-sqrt(0.75d0),0d0]
       do i=3,5
          p(0,i)=sqrt(q**2+model%get_mass(process(i))**2)
          total_energy=total_energy+p(0,i)
       enddo
       pair_start=6
    else
       pair_start=3
    endif

    ipair=0
    do i=pair_start,n,2
       j=i+1
       ipair=ipair+1
       if (j.gt.n) then
          write (*,*) 'Internal error in generic matrix-element point pairing'
          stop 1
       endif
       if (ipair.gt.size(directions,2)) then
          write (*,*) 'Not enough generic directions for matrix-element point:',n
          stop 1
       endif
       q=70d0+23d0*dble(ipair)+5d0*dble(n)
       p(1:3,i)= q*directions(1:3,ipair)
       p(1:3,j)=-q*directions(1:3,ipair)
       mass_i=model%get_mass(process(i))
       mass_j=model%get_mass(process(j))
       p(0,i)=sqrt(q**2+mass_i**2)
       p(0,j)=sqrt(q**2+mass_j**2)
       total_energy=total_energy+p(0,i)+p(0,j)
    enddo

    p(0:3,1)=[0.5d0*total_energy,0d0,0d0, 0.5d0*total_energy]
    p(0:3,2)=[0.5d0*total_energy,0d0,0d0,-0.5d0*total_energy]

    if (maxval(abs(sum(p(1:3,3:n),dim=2))).gt.1d-9) then
       write (*,*) 'Internal error: generic final-state momenta do not sum to zero'
       stop 1
    endif

  end subroutine fill_generic_momenta

  subroutine fill_massless4(p)
    implicit none
    real(kind=dp),dimension(0:3,4),intent(out) :: p
    p(0:3,1)=[500.0000d0,  0.00000d0,  0.0000000d0,  500.0000d0]
    p(0:3,2)=[500.0000d0,  0.00000d0,  0.0000000d0, -500.0000d0]
    p(0:3,3)=[500.0000d0,  110.9243d0,  444.8308d0, -199.5529d0]
    p(0:3,4)=[500.0000d0, -110.9243d0, -444.8308d0,  199.5529d0]
  end subroutine fill_massless4

  subroutine fill_massless5(p)
    implicit none
    real(kind=dp),dimension(0:3,5),intent(out) :: p
    p(0:3,1)=[500.0000d0,  0.0000000d0,  0.0000000d0,  500.0000d0]
    p(0:3,2)=[500.0000d0,  0.0000000d0,  0.0000000d0, -500.0000d0]
    p(0:3,3)=[495.9179d0, -19.99048d0,  79.81352d0, -489.0448d0]
    p(0:3,4)=[132.8543d0, -26.48250d0, -44.16981d0,  122.4662d0]
    p(0:3,5)=[371.2277d0,  46.47298d0, -35.64372d0,  366.5785d0]
  end subroutine fill_massless5

  subroutine fill_massless6(p)
    implicit none
    real(kind=dp),dimension(0:3,6),intent(out) :: p
    p(0:3,1)=[500.0000000d0,   0.00000000d0,  0.0000000000d0,   500.0000000d0]
    p(0:3,2)=[500.0000000d0,   0.00000000d0,  0.0000000000d0,  -500.0000000d0]
    p(0:3,3)=[88.55133305d0,  -22.1006902d0,  40.080353191d0,  -75.80543095d0]
    p(0:3,4)=[328.3294192d0,  -103.849611d0, -301.93375538d0,   76.49492138d0]
    p(0:3,5)=[152.3581094d0,  -105.880959d0, -97.709638326d0,   49.54838522d0]
    p(0:3,6)=[430.7611382d0,   231.831261d0,  359.56304052d0,  -50.23787565d0]
  end subroutine fill_massless6

  subroutine fill_massive5(process,imassive,p)
    implicit none
    integer,dimension(5),intent(in) :: process
    integer,intent(in) :: imassive
    real(kind=dp),dimension(0:3,5),intent(out) :: p
    real(kind=dp) :: mass,ejet,zdir,sgn
    integer :: i,njets
    mass=model%get_mass(process(imassive))
    ejet=(1000d0-mass)/2d0
    if (ejet.le.0d0) then
       write (*,*) 'Massive reference point has non-positive jet energy:',mass
       stop 1
    endif
    zdir=sqrt(0.75d0)
    p=0d0
    p(0:3,1)=[500d0,0d0,0d0,500d0]
    p(0:3,2)=[500d0,0d0,0d0,-500d0]
    njets=0
    do i=3,5
       if (i.eq.imassive) then
          p(0:3,i)=[mass,0d0,0d0,0d0]
       else
          njets=njets+1
          if (njets.eq.1) then
             sgn=1d0
          else
             sgn=-1d0
          endif
          p(0:3,i)=[ejet,sgn*0.3d0*ejet,sgn*0.4d0*ejet,sgn*zdir*ejet]
       endif
    enddo
    if (njets.ne.2) then
       write (*,*) 'Expected two massless final particles for massive5 point'
       stop 1
    endif
  end subroutine fill_massive5

  subroutine unsupported_point(n,process,point)
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(in) :: process
    character(len=*),intent(in) :: point
    write (*,*) 'Unsupported matrix-element reference point:',trim(point),'n=',n
    write (*,*) process(1:n)
    stop 1
  end subroutine unsupported_point

  subroutine compare_case(case_id,family,group_id,row_id,n,n_amps,spins,amps,&
       g_n_amps,g_spins,g_amps,total_failures,total_checks)
    implicit none
    integer,intent(in) :: case_id,group_id,row_id,n,n_amps,g_n_amps
    character(len=*),intent(in) :: family
    integer,dimension(n,n_amps),intent(in) :: spins
    integer,dimension(n,g_n_amps),intent(in) :: g_spins
    complex(kind=dp),dimension(n_amps),intent(in) :: amps
    complex(kind=dp),dimension(g_n_amps),intent(in) :: g_amps
    integer,intent(inout) :: total_failures,total_checks
    logical,dimension(:),allocatable :: used
    integer :: i,j,found
    real(kind=dp) :: diff,scale,rel
    allocate(used(n_amps))
    used=.false.
    do i=1,g_n_amps
       total_checks=total_checks+1
       found=0
       do j=1,n_amps
          if ((.not.used(j)).and.all(spins(1:n,j).eq.g_spins(1:n,i))) then
             found=j
             exit
          endif
       enddo
       if (found.eq.0) then
          if (abs(g_amps(i)).gt.abs_tol) then
             if (should_print_failure(total_failures)) then
                call print_failure_header(case_id,family,group_id,row_id)
                write (*,*) 'Missing helicity:',g_spins(1:n,i),g_amps(i)
             endif
             total_failures=total_failures+1
             cycle
          else
             cycle
          endif
       endif
       used(found)=.true.
       diff=abs(amps(found)-g_amps(i))
       scale=max(abs(amps(found)),abs(g_amps(i)))
       if (scale.gt.0d0) then
          rel=diff/scale
       else
          rel=0d0
       endif
       if (diff.gt.abs_tol+rel_tol*scale) then
          if (should_print_failure(total_failures)) then
             call print_failure_header(case_id,family,group_id,row_id)
             write (*,*) 'Helicity:',g_spins(1:n,i)
             write (*,*) 'Expected:',dble(g_amps(i)),aimag(g_amps(i))
             write (*,*) 'Actual  :',dble(amps(found)),aimag(amps(found))
             write (*,*) 'Diff/rel:',diff,rel
          endif
          total_failures=total_failures+1
       endif
    enddo
    do i=1,n_amps
       if (.not.used(i)) then
          if (abs(amps(i)).gt.abs_tol) then
             if (should_print_failure(total_failures)) then
                call print_failure_header(case_id,family,group_id,row_id)
                write (*,*) 'Missing helicity:',spins(1:n,i),amps(i)
             endif
             total_failures=total_failures+1
          endif
       endif
    enddo
  end subroutine compare_case

  logical function should_print_failure(total_failures)
    implicit none
    integer,intent(in) :: total_failures
    should_print_failure=total_failures.lt.max_failures_to_print
    if (total_failures.gt.max_errors) then
       write (*,*) 'Found ',total_failures,' errors. Quitting.'
       stop 1
    endif
  end function should_print_failure

  subroutine print_failure_header(case_id,family,group_id,row_id)
    implicit none
    integer,intent(in) :: case_id,group_id,row_id
    character(len=*),intent(in) :: family
    write (*,*) '--- Matrix-element mismatch ---'
    write (*,*) 'case/family/group/row:',case_id,trim(family),group_id,row_id
  end subroutine print_failure_header

  subroutine write_int_array(iunit,values)
    implicit none
    integer,intent(in) :: iunit
    integer,dimension(:),intent(in) :: values
    integer :: i
    do i=1,size(values)
       if (i.lt.size(values)) then
          write(iunit,'(i0,1x)',advance='no') values(i)
       else
          write(iunit,'(i0)') values(i)
       endif
    enddo
  end subroutine write_int_array

  subroutine write_amplitude_rows(iunit,spins,amps)
    implicit none
    integer,intent(in) :: iunit
    integer,dimension(:,:),intent(in) :: spins
    complex(kind=dp),dimension(:),intent(in) :: amps
    integer :: i,j
    do i=1,size(amps)
       do j=1,size(spins,1)
          write(iunit,'(i0,1x)',advance='no') spins(j,i)
       enddo
       write(iunit,'(es26.17e3,1x,es26.17e3)') dble(amps(i)),aimag(amps(i))
    enddo
  end subroutine write_amplitude_rows

end program matrix_element_regression
