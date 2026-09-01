module handling_processes
  use common
  use amplitude_QCD_mod
  use phase_space_base
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: max_process_records=20000000
  integer(kind=8),parameter :: max_process_workspace_bytes=2147483648_8
  integer(kind=8),parameter :: max_process_colour_orders=1000000_8
  integer(kind=8),parameter :: max_process_order_comparisons=100000000_8
  private :: handling_light_quark_code,handling_same_light_flavour
  type :: multichan_info
     ! if adding variables here, also update the finalize_multichan_info subroutine
     integer,dimension(:,:),allocatable :: channels,unique_channelgroup_list
     integer,dimension(:,:,:),allocatable :: channel_permutations
     integer,dimension(:),allocatable :: unique_channel_list,map_proc_to_channelgroup,number_of_channels
     integer :: max_channels=0,n_unique_channels=0,n_unique_channelgroups=0
   contains
     final :: finalize_multichan_info
  end type multichan_info
  type dipole
     integer,dimension(:),allocatable :: process_r,dip_map,reduced_color_order
     integer,dimension(3) :: dip_ijk=0,dip_ijk_f=0
     integer,dimension(2) :: dip_r_ijk=0,dip_r_ijk_f=0
     integer :: dipole_type=0 ! bit convention: 0:II, 1:FI, 2:IF, 3:FF
     integer :: col_fac=1
     real(kind=8) :: lc_weight=1d0
     real(kind=8) :: alpha_variable=huge(1d0)
     type(amplitude_QCD) :: amp
     real(kind=8),dimension(0:3) :: p_mapped_ij=0d0
     real(kind=8),dimension(:,:),allocatable :: p_mapped
     logical :: active=.true.,alpha_active=.true.,passes_cuts=.true.
     integer,dimension(:),allocatable :: rho_lookup_ih1,rho_lookup_ih2
     logical :: rho_lookup_upper=.false.,rho_hermitian_checked=.false.
   contains
     final :: finalize_dipole
  end type dipole
  type dipole_set
     type(dipole),dimension(:),allocatable :: dl
     integer :: ndip=0
   contains
     final :: finalize_dipole_set
  end type dipole_set
  type phase_space_order_group
     ! if adding variables here, also update the finalize_phase_space_order_group subroutine
     type(amplitude_QCD),dimension(:),allocatable :: amps
     class(phase_space_type),allocatable :: phase_space
     type(multichan_info) :: multichan
     type(dipole_set),dimension(:),allocatable :: dpl
     type(psv),dimension(:),allocatable :: ps
     integer,dimension(:,:),allocatable :: processes,color_orders
     integer,dimension(:,:),allocatable :: phase_space_permutations
     integer,dimension(:),allocatable :: iden_iproc,phase_space_orders,nhel
     integer :: nproc=0
     real(kind=8),dimension(:,:),allocatable :: val_procs,idenCOandMAPfactor
     integer,dimension(:,:,:),allocatable :: iden_processes,same_flavour
     integer(kind=4),dimension(:,:),allocatable :: spin,hel_fac
     integer(kind=8),dimension(:),allocatable :: iden
     logical,dimension(-6:7,2) :: ipdgs=.false.
     integer(kind=4) :: next=0,ndim=0,ndim_extra=0
     logical :: is_subtracted_real=.false.
     integer(kind=8) :: event_selection_epoch=0_8
     integer,dimension(:),allocatable :: col_fac
     real(kind=8),dimension(:),allocatable :: amp2,amp2_hel
     real(kind=8),dimension(:,:,:),allocatable :: amp2_hel_samples
     integer(kind=4),dimension(:),allocatable :: hel,passed
     integer,dimension(:,:),allocatable :: include_hel
     ! cuts
     double precision,dimension(:),allocatable :: pT_min,eta_max
     double precision,dimension(:,:),allocatable :: DR_min,sqrt_s_min
   contains
     final :: finalize_phase_space_order_group
  end type phase_space_order_group
  integer :: next=0,nproc_unique=0,ngroups=0,nprocs=0,c_o=0,nquarks=0
  type(phase_space_order_group),dimension(:),allocatable :: pgl
  logical :: read_proc_from_file=.false.
  integer(kind=4),dimension(:),allocatable :: o,part
contains
  pure logical function handling_light_quark_code(flavour)
    integer,intent(in) :: flavour
    handling_light_quark_code=(flavour.ge.1 .and. flavour.le.6) .or. &
         (flavour.le.-1 .and. flavour.ge.-6)
  end function handling_light_quark_code

  pure logical function handling_same_light_flavour(first,second)
    integer,intent(in) :: first,second
    handling_same_light_flavour=.false.
    if (.not.handling_light_quark_code(first) .or. &
         .not.handling_light_quark_code(second)) return
    handling_same_light_flavour=first.eq.second .or. first.eq.-second
  end function handling_same_light_flavour

  subroutine determine_phase_space_orders(part,col_o,n_ps,PS_o)
    use math_functions
    implicit none
    integer,dimension(1:next),intent(in) :: part,col_o
    integer,intent(in) :: n_ps
    integer,dimension(1:next,1:n_ps),intent(inout) :: ps_o
    integer :: i,n_sing,i_sing,j
    integer,dimension(:),allocatable :: jmap,ips,ips_out
    if (n_ps.eq.1) then
       PS_o(1:next,1)=col_o(1:next)
       return
    endif
    n_sing=0
    do i=1,next
       if (phys_model%is_singlet(part(i))) n_sing=n_sing+1
    enddo
    allocate(jmap(n_sing))
    j=0
    do i=1,next
       if (phys_model%is_singlet(part(o(i)))) then
          j=j+1
          jmap(j)=i
       endif
    enddo
    allocate(ips(n_sing))
    allocate(ips_out(n_sing))
    do i=1,n_sing
       ips(i)=i
    enddo
    do j=1,n_ps
       i_sing=0
       do i=1,next
          if (phys_model%is_singlet(part(o(i)))) then
             i_sing=i_sing+1
             ps_o(i,j)=col_o(jmap(ips(i_sing)))
          else
             ps_o(i,j)=col_o(i)
          endif
       enddo
       if (j.eq.n_ps) exit
       call get_next_iperm(n_sing,ips,ips_out,n_sing)
       ips(1:n_sing)=ips_out(1:n_sing)
    enddo
  end subroutine determine_phase_space_orders

  subroutine find_same_flavour(pgl,nevent,amp2)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    real(kind=8),dimension(pgl%nproc),intent(in) :: amp2
    integer,intent(in) :: nevent
    integer :: i,j,k,ii,jj,kk,allocation_status
    integer(kind=8) :: workspace_bytes
    character(len=256) :: allocation_message
    real(kind=8),parameter :: tiny=1d-8
    real(kind=8) :: amplitude_scale,normalized_norm,normalized_residual
    complex(kind=8) :: amp_i,amp_j,amp_k
    if (keep_processes_separate) return
    if (.not.decompose_same_flavour_into_two_diff_flavour) return
    if (nevent.lt.2 .or. pgl%nproc.lt.1 .or. &
         pgl%nproc.gt.max_process_records) then
       write (*,*) 'Invalid same-flavour reduction dimensions:',nevent,pgl%nproc
       stop 1
    endif
    if (.not.all(ieee_is_finite(amp2))) then
       write (*,*) 'Invalid matrix elements in same-flavour reduction'
       stop 1
    endif
    if (any(amp2.lt.0d0)) then
       write (*,*) 'Invalid matrix elements in same-flavour reduction'
       stop 1
    endif
    if (.not.allocated(pgl%amps) .or. .not.allocated(pgl%passed)) then
       write (*,*) 'Incomplete amplitude state in same-flavour reduction'
       stop 1
    endif
    if (size(pgl%amps).lt.1) then
       write (*,*) 'Empty amplitude state in same-flavour reduction'
       stop 1
    endif
    if (size(pgl%passed).lt.1) then
       write (*,*) 'Empty sample state in same-flavour reduction'
       stop 1
    endif
    if (pgl%passed(1).lt.1 .or. pgl%passed(1).gt.nevent) then
       write (*,*) 'Invalid sample index in same-flavour reduction:',&
            pgl%passed(1),nevent
       stop 1
    endif
    if (.not.allocated(pgl%amps(1)%n_qqbar)) then
       write (*,*) 'Missing quark-line metadata in same-flavour reduction'
       stop 1
    endif
    if (size(pgl%amps(1)%n_qqbar).ne.pgl%nproc) then
       write (*,*) 'Incompatible quark-line metadata in same-flavour reduction'
       stop 1
    endif
    if (all(pgl%amps(1)%n_qqbar.lt.2)) return
    if (.not.allocated(pgl%amps(1)%iproc_start) .or. &
         .not.allocated(pgl%amps(1)%amps) .or. &
         .not.allocated(pgl%amps(1)%spins) .or. &
         .not.allocated(pgl%amps(1)%same_flav) .or. &
         .not.allocated(pgl%amps(1)%same_flavour_sum)) then
       write (*,*) 'Incomplete amplitude metadata in same-flavour reduction'
       stop 1
    endif
    if (size(pgl%amps(1)%iproc_start).ne.pgl%nproc+1 .or. &
         size(pgl%amps(1)%same_flav).ne.pgl%nproc .or. &
         size(pgl%amps(1)%amps).lt.pgl%amps(1)%n_amps .or. &
         size(pgl%amps(1)%spins,3).lt.pgl%amps(1)%n_amps .or. &
         size(pgl%amps(1)%same_flavour_sum,1).lt.pgl%amps(1)%n_amps) then
       write (*,*) 'Incompatible amplitude metadata in same-flavour reduction'
       stop 1
    endif
    workspace_bytes=8_8*int(nevent,kind=8)*int(pgl%nproc,kind=8)
    if (workspace_bytes.gt.max_process_workspace_bytes) then
       write (*,*) 'Same-flavour reduction exceeds the supported workspace:',&
            workspace_bytes,max_process_workspace_bytes
       stop 1
    endif
    if (.not.allocated(pgl%same_flavour)) then
       allocate(pgl%same_flavour(nevent,pgl%nproc,2),stat=allocation_status,&
            errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate same-flavour reduction samples: ',&
               trim(allocation_message)
          stop 1
       endif
       pgl%same_flavour=0
    elseif (size(pgl%same_flavour,1).ne.nevent .or. &
         size(pgl%same_flavour,2).ne.pgl%nproc .or. &
         size(pgl%same_flavour,3).ne.2) then
       write (*,*) 'Incompatible same-flavour reduction sample storage'
       stop 1
    endif
    do i=1,pgl%nproc
       if (pgl%amps(1)%n_qqbar(i).lt.2) cycle
       do j=1,pgl%nproc
          if (i.eq.j) cycle
          if (pgl%amps(1)%n_qqbar(j).ne.pgl%amps(1)%n_qqbar(i)) cycle
          do k=1,j-1
             if (k.eq.i) cycle
             if (pgl%amps(1)%n_qqbar(k).ne.pgl%amps(1)%n_qqbar(i)) cycle
             do ii=pgl%amps(1)%iproc_start(i),pgl%amps(1)%iproc_start(i+1)-1
                if (pgl%amps(1)%amps(ii).eq.(0d0,0d0)) cycle
                do jj=pgl%amps(1)%iproc_start(j),pgl%amps(1)%iproc_start(j+1)-1
                   if (all(pgl%amps(1)%spins(:,1,ii).eq.pgl%amps(1)%spins(:,1,jj))) exit
                enddo
                if (jj.ge.pgl%amps(1)%iproc_start(j+1)) exit
                do kk=pgl%amps(1)%iproc_start(k),pgl%amps(1)%iproc_start(k+1)-1
                   if (all(pgl%amps(1)%spins(:,1,ii).eq.pgl%amps(1)%spins(:,1,kk))) exit
                enddo
                if (kk.ge.pgl%amps(1)%iproc_start(k+1)) exit
                amp_i=pgl%amps(1)%amps(ii)
                amp_j=pgl%amps(1)%amps(jj)
                amp_k=pgl%amps(1)%amps(kk)
                if (.not.ieee_is_finite(real(amp_i,kind=8)) .or. &
                     .not.ieee_is_finite(aimag(amp_i)) .or. &
                     .not.ieee_is_finite(real(amp_j,kind=8)) .or. &
                     .not.ieee_is_finite(aimag(amp_j)) .or. &
                     .not.ieee_is_finite(real(amp_k,kind=8)) .or. &
                     .not.ieee_is_finite(aimag(amp_k))) then
                   write (*,*) 'Non-finite amplitude in same-flavour reduction:',i,j,k
                   stop 1
                endif
                amplitude_scale=max(abs(real(amp_i,kind=8)),abs(aimag(amp_i)),&
                     abs(real(amp_j,kind=8)),abs(aimag(amp_j)),&
                     abs(real(amp_k,kind=8)),abs(aimag(amp_k)))
                if (amplitude_scale.eq.0d0) cycle
                normalized_norm=abs(amp_i/amplitude_scale)+abs(amp_j/amplitude_scale)+&
                     abs(amp_k/amplitude_scale)
                normalized_residual=abs(amp_i/amplitude_scale-&
                     (amp_j/amplitude_scale+amp_k/amplitude_scale))
                if (normalized_norm.le.0d0 .or. &
                     .not.ieee_is_finite(normalized_norm) .or. &
                     .not.ieee_is_finite(normalized_residual)) then
                   write (*,*) 'Invalid normalized amplitude in same-flavour reduction:',i,j,k
                   stop 1
                endif
                if (normalized_residual/normalized_norm.gt.tiny) then
                   exit
                endif
             enddo
             if (ii.eq.pgl%amps(1)%iproc_start(i+1)) then
                pgl%same_flavour(pgl%passed(1),i,1)=j
                pgl%same_flavour(pgl%passed(1),i,2)=k
             endif
          enddo
       enddo
    enddo
    if (pgl%passed(1).lt.nevent) return
    do i=1,pgl%nproc
       if ( any(pgl%same_flavour(1,i,1).ne.pgl%same_flavour(2:nevent,i,1)) .or. &
            any(pgl%same_flavour(1,i,2).ne.pgl%same_flavour(2:nevent,i,2)) ) then
          write (*,*) 'Inconsistent same flavour decomposition'
          write (*,*) i
          write (*,*) pgl%same_flavour(1:nevent,i,1)
          write (*,*) pgl%same_flavour(1:nevent,i,2)
          stop 1
       endif
       if (pgl%same_flavour(1,i,1).ne.0 .or. pgl%same_flavour(1,i,2).ne.0) then
          j=pgl%same_flavour(1,i,1)
          k=pgl%same_flavour(1,i,2)
          write (99,'(a,x,i4,x,a,i4,x,a,i4)') &
               "Found SF amps equal to a sum of DF amps:",i,'=',j,'+',k
          pgl%amps(1)%same_flav(i)=.true.
          do ii=pgl%amps(1)%iproc_start(i),pgl%amps(1)%iproc_start(i+1)-1
             do jj=pgl%amps(1)%iproc_start(j),pgl%amps(1)%iproc_start(j+1)-1
                if (all(pgl%amps(1)%spins(:,1,ii).eq.pgl%amps(1)%spins(:,1,jj))) exit
             enddo
             if (jj.ge.pgl%amps(1)%iproc_start(j+1)) then
                write (*,*) 'Missing first helicity partner in same-flavour decomposition',i,j,ii
                stop 1
             endif
             do kk=pgl%amps(1)%iproc_start(k),pgl%amps(1)%iproc_start(k+1)-1
                if (all(pgl%amps(1)%spins(:,1,ii).eq.pgl%amps(1)%spins(:,1,kk))) exit
             enddo
             if (kk.ge.pgl%amps(1)%iproc_start(k+1)) then
                write (*,*) 'Missing second helicity partner in same-flavour decomposition',i,k,ii
                stop 1
             endif
             pgl%amps(1)%same_flavour_sum(ii,1)=jj
             pgl%amps(1)%same_flavour_sum(ii,2)=kk
          enddo
       endif
    enddo
  end subroutine find_same_flavour

  subroutine setup_spin(pgl)
    ! Use the first process in the processes() array to setup all the possible
    ! spin states. Note that this assumes that all the processes() have the
    ! same number of spin states
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,iproc,allocation_status
    character(len=256) :: allocation_message
    if (pgl%next.lt.3 .or. pgl%next.gt.max_amplitude_external_particles .or. &
         pgl%nproc.lt.1 .or. .not.allocated(pgl%processes)) then
       write (*,*) 'Invalid process dimensions while setting up spins:',&
            pgl%next,pgl%nproc
       stop 1
    endif
    if (size(pgl%processes,1).ne.pgl%next .or. &
         size(pgl%processes,2).ne.pgl%nproc) then
       write (*,*) 'Incompatible process array while setting up spins'
       stop 1
    endif
    if (.not.allocated(pgl%spin)) then
       allocate(pgl%spin(0:3,1:pgl%next),stat=allocation_status,&
            errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate process spin table: ',&
               trim(allocation_message)
          stop 1
       endif
    elseif (size(pgl%spin,1).ne.4 .or. size(pgl%spin,2).ne.pgl%next) then
       write (*,*) 'Existing process spin table has incompatible dimensions'
       stop 1
    endif
    pgl%spin=0
    do i=1,pgl%next
       pgl%spin(0,i)=phys_model%get_spin(pgl%processes(i,1))
       if (pgl%spin(0,i).eq.2) then
          pgl%spin(1,i)=-1
          pgl%spin(2,i)=1
       elseif (pgl%spin(0,i).eq.3) then
          pgl%spin(1,i)=-1
          pgl%spin(2,i)=0
          pgl%spin(3,i)=1
       elseif (pgl%spin(0,i).eq.1) then
          pgl%spin(1,i)=0
       else
          write (*,*) 'spin state not known',i,pgl%processes(i,1),pgl%spin(0,i)
          stop 1
       endif
    enddo
    do iproc=2,pgl%nproc
       do i=1,pgl%next
          if (pgl%spin(0,i).ne.phys_model%get_spin(pgl%processes(i,iproc))) then
             write (*,*) 'Spin states of particles in different processes not compatible',iproc
             stop 1
          endif
       enddo
    enddo
  end subroutine setup_spin

  subroutine setup_color_order(pgl_unique)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl_unique
    integer :: i,iproc,nq,ng,nsing,iq,iaq,is,ig
    do iproc=1,pgl_unique%nproc
       nq=0
       ng=0
       nsing=0
       do i=1,pgl_unique%next
          if (handling_light_quark_code(pgl_unique%processes(i,iproc))) then
             nq=nq+1
          elseif(phys_model%is_gluon(pgl_unique%processes(i,iproc))) then
             ng=ng+1
          elseif(phys_model%is_singlet(pgl_unique%processes(i,iproc))) then
             nsing=nsing+1
          else
             write (*,*) 'unknown particle type:',pgl_unique%processes(i,iproc)
             stop 1
          endif
       enddo

       if (nq.eq.0 .and. nsing.ne.0) then
          ig=nsing+1
          is=1
          do i=1,pgl_unique%next
             if (phys_model%is_singlet(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(is,iproc)=i
                is=is+1
             else
                pgl_unique%color_orders(ig,iproc)=i
                ig=ig+1
             endif
          enddo
!!$          write (*,*) 'when there are colour singlets, there should be quarks'
!!$          stop 1
       elseif (nq.eq.0) then
          do i=1,pgl_unique%next
             pgl_unique%color_orders(i,iproc)=i
          enddo
       elseif (nq.eq.2) then
          ig=2
          is=ng+2
          do i=1,pgl_unique%next
             if (phys_model%is_quark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(1,iproc)=i
             elseif (phys_model%is_antiquark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(pgl_unique%next,iproc)=i
             elseif (phys_model%is_gluon(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(ig,iproc)=i
                ig=ig+1
             elseif (phys_model%is_singlet(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(is,iproc)=i
                is=is+1
             endif
          enddo
       elseif (nq.eq.4 .or. nq.eq.6) then
          iq=1
          iaq=2
          ! The first nq-1 slots hold alternating q/qbar endpoints.  Put
          ! gluons next, then singlets, while the last antiquark closes the
          ! order in slot next.  Using a fixed offset of four here overwrote
          ! the third quark for six-quark processes with a singlet.
          ig=nq
          is=ng+nq
          do i=1,pgl_unique%next
             if (phys_model%is_quark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(iq,iproc)=i
                iq=iq+2
             elseif (phys_model%is_antiquark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(iaq,iproc)=i
                if (iaq.eq.nq-2) then
                   iaq=pgl_unique%next
                else
                   iaq=iaq+2
                endif
             elseif (phys_model%is_gluon(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(ig,iproc)=i
                ig=ig+1
             elseif (phys_model%is_singlet(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(is,iproc)=i
                is=is+1
             endif
          enddo
       else
          write (*,*) 'Unknown number of quarks and anti-quarks'
          write (*,*) iproc,':',pgl_unique%processes(:,iproc)
          stop 1
       endif
    enddo
  end subroutine setup_color_order

  subroutine set_initial_state_average_factor(pgl)
    use math_functions, only: checked_multiply8
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,iproc
    do iproc=1,pgl%nproc
       do i=1,2
          if (pgl%processes(i,iproc).eq.21) then
             ! gluon: two polarisations and 8 colours
             pgl%iden(iproc)=checked_multiply8(pgl%iden(iproc),16_8,&
                  'initial-state average factor')
          elseif (handling_light_quark_code(pgl%processes(i,iproc))) then
             ! (anti-)quark: two helicities and 3 colours
             pgl%iden(iproc)=checked_multiply8(pgl%iden(iproc),6_8,&
                  'initial-state average factor')
          else
             ! assume two spin states and no colour:
             pgl%iden(iproc)=checked_multiply8(pgl%iden(iproc),2_8,&
                  'initial-state average factor')
          endif
       enddo
    enddo
  end subroutine set_initial_state_average_factor

  subroutine set_final_state_identical_particle_factor(pgl)
    use math_functions
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,j,ni,iproc,allocation_status
    integer,dimension(:,:),allocatable :: iden_part
    character(len=256) :: allocation_message
    if (pgl%next.lt.3 .or. pgl%next.gt.max_amplitude_external_particles .or. &
         pgl%nproc.lt.1 .or. .not.allocated(pgl%processes) .or. &
         .not.allocated(pgl%iden)) then
       write (*,*) 'Invalid process state for identical-particle factors'
       stop 1
    endif
    if (size(pgl%processes,1).ne.pgl%next .or. &
         size(pgl%processes,2).ne.pgl%nproc .or. &
         size(pgl%iden).ne.pgl%nproc) then
       write (*,*) 'Incompatible process state for identical-particle factors'
       stop 1
    endif
    allocate(iden_part(1:pgl%next,2),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate identical-particle workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    iden_part=0
    do iproc=1,pgl%nproc
       ni=0
       do i=3,pgl%next
          do j=1,ni
             if (iden_part(j,1).eq.pgl%processes(i,iproc)) then
                iden_part(j,2)=iden_part(j,2)+1
                exit
             endif
          enddo
          if (j.eq.ni+1) then
             ni=ni+1
             iden_part(j,1)=pgl%processes(i,iproc)
             iden_part(j,2)=1
          endif
       enddo
       do i=1,ni
          pgl%iden(iproc)=checked_multiply8(pgl%iden(iproc),factorial8(iden_part(i,2)),&
               'identical-particle factor')
       enddo
    enddo
    deallocate(iden_part)
  end subroutine set_final_state_identical_particle_factor

  subroutine compute_LC_colour_factor(pgl)
    use math_functions, only: checked_integer_power
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,ifac,iproc
    real(kind=8) :: fac
    do iproc=1,pgl%nproc
       fac=0d0
       do i=1,pgl%next
          if (pgl%processes(i,iproc).eq.21) then
             fac=fac+1d0
          elseif (handling_light_quark_code(pgl%processes(i,iproc))) then
             fac=fac+0.5d0
          endif
       enddo
       ifac=nint(fac)
       if (dble(ifac).ne.fac) then
          write (*,*) 'There is some issue with the LC colour factor computation: '// &
               'colour factor is not an integer',ifac,fac
          stop 1
       endif
       pgl%col_fac(iproc)=checked_integer_power(3,ifac,'leading-colour factor')
    enddo
  end subroutine compute_LC_colour_factor

  subroutine define_identical_procs(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: iproc,ip,n,max_identical,allocation_status
    integer(kind=8) :: workspace_bytes
    character(len=256) :: allocation_message
    if (pgl%next.lt.3 .or. pgl%next.gt.max_amplitude_external_particles .or. &
         pgl%nproc.lt.1 .or. pgl%nproc.gt.max_process_records .or. &
         .not.allocated(pgl%processes)) then
       write (*,*) 'Invalid dimensions while defining identical processes:',&
            pgl%next,pgl%nproc
       stop 1
    endif
    if (size(pgl%processes,1).ne.pgl%next .or. &
         size(pgl%processes,2).ne.pgl%nproc) then
       write (*,*) 'Incompatible process array while defining identical processes'
       stop 1
    endif
    if (allocated(pgl%iden_iproc)) deallocate(pgl%iden_iproc)
    if (allocated(pgl%val_procs)) deallocate(pgl%val_procs)
    if (allocated(pgl%iden_processes)) deallocate(pgl%iden_processes)
    ! first fill the number of identical processes per iproc (so that we can
    ! allocate the array with the right size)
    allocate(pgl%iden_iproc(1:pgl%nproc),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate identical-process counts: ',&
            trim(allocation_message)
       stop 1
    endif
    do iproc=1,pgl%nproc
       pgl%iden_iproc(iproc)=1
       if (any(pgl%processes(1:pgl%next,iproc).eq.1 .or. &
            pgl%processes(1:pgl%next,iproc).eq.-1)) then
          pgl%iden_iproc(iproc)=pgl%iden_iproc(iproc)*5
       endif
       if (any(pgl%processes(1:pgl%next,iproc).eq.2 .or. &
            pgl%processes(1:pgl%next,iproc).eq.-2)) then
          pgl%iden_iproc(iproc)=pgl%iden_iproc(iproc)*4
       endif
    enddo
    max_identical=maxval(pgl%iden_iproc)
    workspace_bytes=4_8*int(pgl%nproc,kind=8)+&
         8_8*int(max_identical,kind=8)*int(pgl%nproc,kind=8)+&
         4_8*int(pgl%next,kind=8)*int(max_identical,kind=8)*&
         int(pgl%nproc,kind=8)
    if (workspace_bytes.lt.0_8 .or. &
         workspace_bytes.gt.max_process_workspace_bytes) then
       write (*,*) 'Identical-process expansion exceeds the supported workspace:',&
            workspace_bytes,max_process_workspace_bytes
       stop 1
    endif
    allocate(pgl%val_procs(1:max_identical,1:pgl%nproc),&
         pgl%iden_processes(1:pgl%next,1:max_identical,1:pgl%nproc),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate identical-process expansion: ',&
            trim(allocation_message)
       stop 1
    endif
    pgl%val_procs=0d0
    pgl%iden_processes=0
    ! Loop again and actually fill the iden_processes()
    do iproc=1,pgl%nproc
       do ip=0,pgl%iden_iproc(iproc)-1
          do n=1,pgl%next
             if (pgl%processes(n,iproc).eq.1 .or. &
                  pgl%processes(n,iproc).eq.-1) then
                pgl%iden_processes(n,ip+1,iproc)=sign(mod(ip,5)+1,pgl%processes(n,iproc))
             elseif (pgl%processes(n,iproc).eq.2 .or. &
                  pgl%processes(n,iproc).eq.-2) then
                if (mod(ip,5)+1.eq.ip/5+2) then
                   pgl%iden_processes(n,ip+1,iproc)=sign(1,pgl%processes(n,iproc))
                else
                   pgl%iden_processes(n,ip+1,iproc)=sign(ip/5+2,pgl%processes(n,iproc))
                endif
             else
                pgl%iden_processes(n,ip+1,iproc)=pgl%processes(n,iproc)
             endif
          enddo
       enddo
    enddo
  end subroutine define_identical_procs

  subroutine compute_multichannel_symmetry_factor(sym_fac)
    use math_functions
    implicit none
    integer(kind=8),intent(out) :: sym_fac
    integer :: ngl,ngl_tot,n_sing
    integer,dimension(6) :: nq,naq
    integer :: i,j
    integer(kind=8) :: tot_ord
    integer :: ic,i_qq,iaq,i_ini,i_inv,i_swap,k,l,ios
    integer(kind=8) :: ip,workspace_bytes,comparison_count
    integer,dimension(next) :: fgluons,ips,ips_out
    logical :: same_flavour
    logical,dimension(next) :: fgluon
    integer,dimension(:,:),allocatable :: io_list,io
    character(len=256) :: allocation_message

    if (next.lt.4 .or. next.gt.max_amplitude_external_particles .or. &
         nquarks.lt.0 .or. nquarks.gt.next) then
       write (*,*) 'Invalid multichannel process dimensions:',next,nquarks
       stop 1
    endif
    if (.not.allocated(part) .or. .not.allocated(o)) then
       write (*,*) 'Missing process data for multichannel symmetry factor'
       stop 1
    endif
    if (size(part).lt.next .or. size(o).lt.next) then
       write (*,*) 'Short process data for multichannel symmetry factor'
       stop 1
    endif
    if (any(o(1:next).lt.1) .or. any(o(1:next).gt.next)) then
       write (*,*) 'Out-of-range colour order for multichannel symmetry factor'
       stop 1
    endif

    ngl=0
    ngl_tot=0
    n_sing=0
    tot_ord=0_8
    sym_fac=0_8
    nq=0
    naq=0
    ! count the number of final state gluons and quarks
    do i=3,next
       if (part(i).eq.21) then
          ngl=ngl+1
       endif
       do j=1,6
          if (part(i).eq.j) nq(j)=nq(j)+1
          if (part(i).eq.-j) naq(j)=naq(j)+1
       enddo
    enddo
    do i=1,next
       if (part(i).eq.21) then
          ngl_tot=ngl_tot+1
       endif
    enddo

    ! Since we only need to include a subset of all the colour-orderings, we
    ! need to compensate with a symmetry factor
    if (nquarks.eq.0) then
       tot_ord=factorial8(max(ngl_tot-1,0))
       ! All gluon process. This assumes that the only channels we are
       ! including are strictly different. We distinguish them by considering
       ! how many (final state) gluons are attached to the two colour lines
       ! that link the two incoming gluons. Hence, we only include
       ! floor(next/2) channels, e.g., for next=6 we only consider:
       ! i   --> 1,2,3,4,5,6   (0 and 4 gluons on the two lines)
       ! ii  --> 1,3,2,4,5,6   (1 and 3 gluons on the two lines)
       ! iii --> 1,3,4,2,5,6   (2 and 2 gluons on the two lines)
       ! And, e.g., for next=9, we only consider:
       ! i   --> 1,2,3,4,5,6,7,8,9   (0 and 7 gluons on the two lines)
       ! ii  --> 1,3,2,4,5,6,7,8,9   (1 and 6 gluons on the two lines)
       ! iii --> 1,3,4,2,5,6,7,8,9   (2 and 5 gluons on the two lines)
       ! iv  --> 1,3,4,5,2,6,7,8,9   (3 and 4 gluons on the two lines)
       ! This means that the sym_fac should be equal to the number of final
       ! state gluon permutations, multiplied by 2 (except if we have an equal
       ! number of gluons on both colour lines that attached the two incoming
       ! gluons).
       if (c_o*2.eq.(ngl)) then
          sym_fac=factorial8(ngl)
       else
          sym_fac=checked_multiply8(2_8,factorial8(ngl),'all-gluon multichannel factor')
       endif
    elseif (nquarks.eq.2) then
       tot_ord=factorial8(ngl_tot)
       if (handling_light_quark_code(part(1)) .and. &
            handling_light_quark_code(part(2))) then
          ! quark and anti-quark are incoming. Only 1 channel needed,
          ! which would result in the following symmetry factor:
          sym_fac=factorial8(ngl)
       elseif (handling_light_quark_code(part(1)) .or. &
            handling_light_quark_code(part(2))) then
          ! one incoming quark (or anti-quark). There are ngluons
          ! channels needed: they correspond to having the incoming
          ! gluon at all possible positions between the quark and
          ! anti-quark in the colour order. Hence, each channel comes
          ! with an (ngluons-1)! symmetry factor:
          sym_fac=factorial8(ngl)
       else
          ! both quark and anti-quark are final state. This is similar
          ! to the all-gluon case above, treating the q-qbar pair as
          ! another gluon. This special gluon is identifiable! So, for
          ! next=6 (and assuming that the qqbar pair are particles 5
          ! and 6) one has the following possibilities:
          !
          ! ia   --> 1,2,3,4,(5,6) = 5,4,3,2,1,6 ---- : both gluons on the same
          ! ib   --> 1,2,3,(5,6),4 = 5,3,2,1,4,6 --/         line as the qqbar pair
          ! ic   --> 1,2,(5,6),3,4 = 5,2,1,4,3,6 -/
          ! iia  --> 1,3,2,4,(5,6) = 5,4,2,3,1,6 ---- : one gluon on the same 
          ! iib  --> 1,3,2,(5,6),4 = 5,2,3,1,4,6 -/          line as the qqbar pair
          ! iii  --> 1,3,4,2,(5,6) = 5,2,4,3,1,6 ---- : both gluons on the other quark line
          !
          ! Note, however, that in the above situation ia and ic results in
          ! the same cross section. Also iia and iib result in the same
          ! rate. Furthermore, there is an additional factor 2, since one can
          ! interchange the two incoming particles. Hence, there are only
          ! trully 4 independent colour orders to consider:
          ! ia, with symmetry factor (ngluon-2)!*2*2
          ! ib, with symmetry factor (ngluon-2)!*2
          ! iaa, with symmetry factor (ngluon-2)!*2*2
          ! iii, with symmetry factor (ngluon-2)!*2
          ! Hence
          sym_fac=checked_multiply8(factorial8(ngl),2_8,'two-quark multichannel factor')
          if (ngl.gt.0) then
             sym_fac=checked_multiply8(factorial8(ngl),2_8,'two-quark multichannel factor')
          endif
          if (ifindloc(o,next,1).ne.next-ifindloc(o,next,2)+1) then
             sym_fac=checked_multiply8(sym_fac,2_8,'two-quark multichannel factor')
          endif
       endif
    elseif (nquarks.eq.4) then
       ! total number of potentially different orders: two ways of connecting
       ! quarks, (n-4)! orderings for the gluons, n-3 ways for an order to
       ! distribute the gluons among the two quark lines
       tot_ord=checked_multiply8(2_8,factorial8(ngl_tot),'four-quark colour-order count')
       tot_ord=checked_multiply8(tot_ord,int(ngl_tot,kind=8)+1_8,&
            'four-quark colour-order count')
       n_sing=next-4-ngl_tot

       same_flavour=.false.
       do i=2,next-1
          if ((o(i).gt.2 .and. part(o(i)).le.-1 .and. part(o(i)).ge.-6) .or. &
               (o(i).le.2 .and. part(o(i)).ge. 1 .and. part(o(i)).le. 6) ) then
             if (handling_same_light_flavour(part(o(i)),part(o(1))) .and. &
                  handling_same_light_flavour(part(o(i)),part(o(next)))) then
                same_flavour=.true.
             else
                same_flavour=.false.
             endif
             exit
          endif
       enddo

       if (tot_ord.gt.int(huge(ic),kind=8)) then
          write (*,*) 'ERROR: multichannel colour-order list exceeds supported integer size:',tot_ord
          stop 1
       endif
       if (tot_ord.lt.1_8 .or. tot_ord.gt.max_process_colour_orders) then
          write (*,*) 'ERROR: multichannel colour-order count exceeds the supported limit:',&
               tot_ord,max_process_colour_orders
          stop 1
       endif
       workspace_bytes=4_8*int(next,kind=8)*(tot_ord+5_8)
       if (workspace_bytes.gt.max_process_workspace_bytes) then
          write (*,*) 'ERROR: multichannel colour orders exceed the supported workspace:',&
               workspace_bytes,max_process_workspace_bytes
          stop 1
       endif
       allocate(io_list(1:next,int(tot_ord)),stat=ios,&
            errmsg=allocation_message)
       if (ios.ne.0) then
          write (*,*) 'ERROR: could not allocate multichannel colour-order list:',&
               next,tot_ord,trim(allocation_message)
          stop 1
       endif
       allocate(io(1:next,5),stat=ios,errmsg=allocation_message)
       if (ios.ne.0) then
          write (*,*) 'ERROR: could not allocate multichannel order workspace:',&
               next,trim(allocation_message)
          stop 1
       endif
       ic=0
       comparison_count=0_8
       ! 1. two ways of connecting quarks with anti-quarks
       do_i_qq: do i_qq=1,2
          io(1:next,1)=o(1:next)
          if (i_qq.eq.2) then
             if (.not. same_flavour) cycle
             do i=2,next-1
                if ((o(i).gt.2 .and. part(o(i)).le.-1 .and. part(o(i)).ge.-6) .or. &
                     (o(i).le.2 .and. part(o(i)).ge. 1 .and. part(o(i)).le. 6) ) then
                   ! check if it would give the same result:
                   if ( ((o(1   ).le.2 .and. o(i+1).gt.2) .or. (o(1   ).gt.2 .and. o(i+1).le.2)) .and. &
                        ((o(next).le.2 .and. o(i  ).gt.2) .or. (o(next).gt.2 .and. o(i  ).le.2)) ) then
                      ! not both initial or both final state
                      cycle do_i_qq
                   endif
                   iaq=o(i)
                   io(i,1)=o(next)
                   io(next,1)=iaq
                   exit
                endif
             enddo
          endif

          ! 2. invert order of two initial states
          do i_ini=1,2
             io(1:next,2)=io(1:next,1)
             if (part(1).ne.part(2) .and. i_ini.eq.2) cycle
             if (i_ini.eq.2) then
                do i=1,next
                   if (io(i,1).eq.1) then
                      io(i,2)=2
                   elseif (io(i,1).eq.2) then
                      io(i,2)=1
                   endif
                enddo
             endif
             ! 3. invert order of both quark lines
             do_i_inv: do i_inv=1,2
                io(1:next,3)=io(1:next,2)
                if (i_inv.eq.2) then
                   ! invert order
                   do i=2,next-n_sing-1
                      if ((io(i,2).gt.2 .and. part(io(i,2)).le.-1 .and. part(io(i,2)).ge.-6) .or. &
                           (io(i,2).le.2 .and. part(io(i,2)).ge. 1 .and. part(io(i,2)).le. 6) ) then
                         if ( ((io(1   ,2).le.2 .and. io(i  ,2).gt.2) .or. (io(1   ,2).gt.2 .and. io(i  ,2).le.2)) .or. &
                              ((io(next,2).le.2 .and. io(i+1,2).gt.2) .or. (io(next,2).gt.2 .and. io(i+1,2).le.2)) ) then
                            ! not both initial or both final state. Cannot invert order.
                            cycle do_i_inv
                         endif
                         io(2:i-1,3)=io(i-1:2:-1,2)
                         io(i+2:next-n_sing-1,3)=io(next-n_sing-1:i+2:-1,2)
                         exit
                      endif
                   enddo
                endif
                ! 4. If both quark lines are similar (identical quarks and
                ! FF+FF or IF+FI or FI+IF or FI+IF or IF+FI), we can swap the gluons from one line to the
                ! other
                do i_swap=1,2
                   io(1:next,4)=io(1:next,3)
                   if (i_swap.eq.2) then
                      if (.not.same_flavour) cycle
                      do i=2,next-1
                         if ((io(i,3).gt.2 .and. part(io(i,3)).le.-1 .and. part(io(i,3)).ge.-6) .or. &
                              (io(i,3).le.2 .and. part(io(i,3)).ge. 1 .and. part(io(i,3)).le. 6) ) then
                            if ((io(1,3).gt.2 .and. io(i,3).gt.2 .and. io(i+1,3).gt.2 .and. io(next,3).gt.2) .or. &
                                 (io(1,3).gt.2 .and. io(i,3).le.2 .and. io(i+1,3).le.2 .and. io(next,3).gt.2) .or. &
                                 (io(1,3).le.2 .and. io(i,3).gt.2 .and. io(i+1,3).gt.2 .and. io(next,3).le.2) .or. &
                                 (io(1,3).gt.2 .and. io(i,3).le.2 .and. io(i+1,3).gt.2 .and. io(next,3).le.2) .or. &
                                 (io(1,3).le.2 .and. io(i,3).gt.2 .and. io(i+1,3).le.2 .and. io(next,3).gt.2) ) then
                               io(2:next-i-1,4)=io(i+2:next-1,3) ! gluons
                               io(next-i:next-i+1,4)=io(i:i+1,3) ! qbarq
                               io(2+next-i:next-1,4)=io(2:i-1,3) ! gluons
                               exit
                            endif
                         endif
                      enddo
                      if (i.eq.next) cycle
                   endif
                   ! 5. permute all final state gluons
                   k=0
                   l=0
                   do i=1,next
                      if (part(io(i,4)).eq.21 .and. io(i,4).gt.2) then
                         k=k+1
                         fgluons(k)=io(i,4)
                         fgluon(i)=.true.
                      else
                         fgluon(i)=.false.
                      endif
                   enddo
                   io(1:next,5)=io(1:next,4)
                   do ip=1_8,factorial8(k)
                      if (ip.eq.1) then
                         do i=1,k
                            ips(i)=i
                         enddo
                      else
                         call get_next_iperm(k,ips,ips_out,k)
                         ips(1:k)=ips_out(1:k)
                      endif
                      i=0
                      l=1
                      do
                         i=i+1
                         if (i.eq.next) exit
                         if (fgluon(i)) then
                            j=0
                            do
                               if (fgluon(i+j+1)) then
                                  j=j+1
                               else
                                  exit
                               endif
                            enddo
                            io(i:i+j,5)=fgluons(ips(l:l+j))
                            l=l+j+1
                            i=i+j
                         endif
                      enddo
                      ! ---> if not yet in list of identical contributions, add it!
                      do i=1,ic
                         if (comparison_count.ge.max_process_order_comparisons) then
                            write (*,*) 'ERROR: multichannel order de-duplication exceeds ',&
                                 'the supported comparison budget:',&
                                 max_process_order_comparisons
                            stop 1
                         endif
                         comparison_count=comparison_count+1_8
                         if (all(io(1:next,5).eq.io_list(1:next,i))) exit
                      enddo
                      if (i.eq.ic+1) then
                         ! new identical contribution
                         io_list(1:next,i)=io(1:next,5)
                         ic=i
                      endif
                   enddo
                enddo
             enddo do_i_inv
          enddo
       enddo do_i_qq
       sym_fac=int(ic,kind=8)
       deallocate(io_list,io)
    else        
       write (*,*) 'ERROR: multichannel symmetry factor is unavailable for quark count',nquarks
       stop 1
    endif
    write (99,*) 'total number of orders is',tot_ord,' and multi-channel symmetry factor is',sym_fac,&
         '. They should be the same when including all channels.'
  end subroutine compute_multichannel_symmetry_factor

  subroutine finalize_multichan_info(mi)
    type(multichan_info),intent(inout) :: mi
    if (allocated(mi%channels)) deallocate(mi%channels)
    if (allocated(mi%channel_permutations)) deallocate(mi%channel_permutations)
    if (allocated(mi%unique_channelgroup_list)) deallocate(mi%unique_channelgroup_list)
    if (allocated(mi%unique_channel_list)) deallocate(mi%unique_channel_list)
    if (allocated(mi%map_proc_to_channelgroup)) deallocate(mi%map_proc_to_channelgroup)
    if (allocated(mi%number_of_channels)) deallocate(mi%number_of_channels)
    mi%max_channels=0
    mi%n_unique_channels=0
    mi%n_unique_channelgroups=0
  end subroutine finalize_multichan_info

  subroutine finalize_dipole(di)
    type(dipole),intent(inout) :: di
    if (allocated(di%process_r)) deallocate(di%process_r)
    if (allocated(di%dip_map)) deallocate(di%dip_map)
    if (allocated(di%reduced_color_order)) deallocate(di%reduced_color_order)
    if (allocated(di%rho_lookup_ih1)) deallocate(di%rho_lookup_ih1)
    if (allocated(di%rho_lookup_ih2)) deallocate(di%rho_lookup_ih2)
    if (allocated(di%p_mapped)) deallocate(di%p_mapped)
    call finalize_amplitude_QCD(di%amp)
    di%dip_ijk=0
    di%dip_ijk_f=0
    di%dip_r_ijk=0
    di%dip_r_ijk_f=0
    di%dipole_type=0
    di%col_fac=1
    di%lc_weight=1d0
    di%alpha_variable=huge(1d0)
    di%p_mapped_ij=0d0
    di%active=.true.
    di%alpha_active=.true.
    di%passes_cuts=.true.
    di%rho_lookup_upper=.false.
    di%rho_hermitian_checked=.false.
  end subroutine finalize_dipole

  subroutine finalize_dipole_set(ds)
    type(dipole_set),intent(inout) :: ds
    integer :: i
    if (allocated(ds%dl)) then
       do i=1,size(ds%dl)
          call finalize_dipole(ds%dl(i))
       enddo
       deallocate(ds%dl)
    endif
    ds%ndip=0
  end subroutine finalize_dipole_set

  subroutine finalize_phase_space_order_group(pgl)
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i
    if (allocated(pgl%amps)) then
       do i=1,size(pgl%amps)
          call finalize_amplitude_QCD(pgl%amps(i))
       enddo
       deallocate(pgl%amps)
    endif
    if(allocated(pgl%phase_space)) then
       call pgl%phase_space%cleanup()
       deallocate(pgl%phase_space)
    endif
    if (allocated(pgl%ps)) then
       do i=1,size(pgl%ps)
          if (allocated(pgl%ps(i)%p)) deallocate(pgl%ps(i)%p)
          if (allocated(pgl%ps(i)%x)) deallocate(pgl%ps(i)%x)
       enddo
       deallocate(pgl%ps)
    endif
    call finalize_multichan_info(pgl%multichan)
    if (allocated(pgl%processes)) deallocate(pgl%processes)
    if (allocated(pgl%color_orders)) deallocate(pgl%color_orders)
    if (allocated(pgl%phase_space_permutations)) deallocate(pgl%phase_space_permutations)
    if (allocated(pgl%iden_iproc)) deallocate(pgl%iden_iproc)
    if (allocated(pgl%phase_space_orders)) deallocate(pgl%phase_space_orders)
    if (allocated(pgl%nhel)) deallocate(pgl%nhel)
    if (allocated(pgl%val_procs)) deallocate(pgl%val_procs)
    if (allocated(pgl%idenCOandMAPfactor)) deallocate(pgl%idenCOandMAPfactor)
    if (allocated(pgl%iden_processes)) deallocate(pgl%iden_processes)
    if (allocated(pgl%spin)) deallocate(pgl%spin)
    if (allocated(pgl%iden)) deallocate(pgl%iden)
    if (allocated(pgl%col_fac)) deallocate(pgl%col_fac)
    if (allocated(pgl%amp2)) deallocate(pgl%amp2)
    if (allocated(pgl%amp2_hel)) deallocate(pgl%amp2_hel)
    if (allocated(pgl%amp2_hel_samples)) deallocate(pgl%amp2_hel_samples)
    if (allocated(pgl%hel)) deallocate(pgl%hel)
    if (allocated(pgl%passed)) deallocate(pgl%passed)
    if (allocated(pgl%hel_fac)) deallocate(pgl%hel_fac)
    if (allocated(pgl%include_hel)) deallocate(pgl%include_hel)
    if (allocated(pgl%pT_min)) deallocate(pgl%pT_min)
    if (allocated(pgl%eta_max)) deallocate(pgl%eta_max)
    if (allocated(pgl%DR_min)) deallocate(pgl%DR_min)
    if (allocated(pgl%sqrt_s_min)) deallocate(pgl%sqrt_s_min)
    if (allocated(pgl%same_flavour)) deallocate(pgl%same_flavour)
    if (allocated(pgl%dpl)) then
       do i=1,size(pgl%dpl)
          call finalize_dipole_set(pgl%dpl(i))
       enddo
       deallocate(pgl%dpl)
    endif
    pgl%nproc=0
    pgl%next=0
    pgl%ndim=0
    pgl%ndim_extra=0
    pgl%ipdgs=.false.
    pgl%is_subtracted_real=.false.
    pgl%event_selection_epoch=0_8
  end subroutine finalize_phase_space_order_group
end module handling_processes
