module handling_processes
  use common
  use amplitude_QCD_mod
  use phase_space_base
  type :: multichan_info
     ! if adding variables here, also update the finalize_multichan_info subroutine
     integer,dimension(:,:),allocatable :: channels,unique_channelgroup_list
     integer,dimension(:,:,:),allocatable :: channel_permutations
     integer,dimension(:),allocatable :: unique_channel_list,map_proc_to_channelgroup,number_of_channels
     integer :: max_channels,n_unique_channels,n_unique_channelgroups
   contains
     final :: finalize_multichan_info
  end type multichan_info
  type phase_space_order_group
     ! if adding variables here, also update the finalize_phase_space_order_group subroutine
     type(amplitude_QCD),dimension(:),allocatable :: amps
     class(phase_space_type),allocatable :: phase_space
     type(multichan_info) :: multichan
     type(psv),dimension(:),allocatable :: ps
     integer,dimension(:,:),allocatable :: processes,color_orders
     integer,dimension(:,:),allocatable :: phase_space_permutations
     integer,dimension(:),allocatable :: iden_iproc,phase_space_orders,nhel
     integer :: nproc
     real(kind=8),dimension(:,:),allocatable :: val_procs,idenCOandMAPfactor
     integer,dimension(:,:,:),allocatable :: iden_processes
     integer,dimension(:,:),allocatable :: same_flavour
     integer(kind=4),dimension(:,:),allocatable :: spin,hel_fac
     integer(kind=8),dimension(:),allocatable :: iden
     logical,dimension(-6:7,2) :: ipdgs
     integer(kind=4) :: next,ndim,ndim_extra
     integer,dimension(:),allocatable :: col_fac
     real(kind=8),dimension(:),allocatable :: amp2,amp2_hel
     integer(kind=4),dimension(:),allocatable :: hel,passed
     integer,dimension(:,:),allocatable :: include_hel
     ! cuts
     double precision,dimension(:),allocatable :: pT_min,eta_max
     double precision,dimension(:,:),allocatable :: DR_min,sqrt_s_min
   contains
     final :: finalize_phase_space_order_group
  end type phase_space_order_group
  integer :: next,nproc_unique,ngroups,nprocs,c_o,nquarks
  type(phase_space_order_group),dimension(:),allocatable :: pgl
  logical :: read_proc_from_file
  integer(kind=4),dimension(:),allocatable :: o,part
contains
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
    ! Identical external fermions admit several Wick contractions.  A
    ! different-flavour process selects a subset of those contractions.  Find
    ! an exact, disjoint cover of the target contractions and validate the
    ! resulting amplitude identity at every sampled phase-space point.
    !
    ! For two quark lines the cover has two terms.  With three lines it can
    ! contain as many as 3!=6 terms (three two-flow terms are common).  The
    ! candidates are restricted to earlier processes, making the dependency
    ! graph acyclic and allowing decomposed amplitudes to depend on one another.
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer,intent(in) :: nevent
    real(kind=8),dimension(pgl%nproc),intent(in) :: amp2
    integer :: i,ii,idau,jj,nterms

    if (keep_processes_separate) return
    if (.not.decompose_same_flavour_amplitudes) return
    if (.not.allocated(pgl%same_flavour)) then
       allocate(pgl%same_flavour(max_same_flavour_terms,pgl%nproc))
       pgl%same_flavour=0
       call build_same_flavour_relations(pgl)
    endif

    do i=1,pgl%nproc
       if (pgl%same_flavour(1,i).le.0) cycle
       if (amp2(i).le.0d0 .or. .not.relation_holds(pgl,i,pgl%same_flavour(:,i))) then
          pgl%same_flavour(:,i)=0
       endif
    enddo
    if (pgl%passed(1).lt.nevent) return

    do i=1,pgl%nproc
       if (pgl%same_flavour(1,i).le.0) cycle
       nterms=count(pgl%same_flavour(:,i).gt.0)
       write (99,'(a,i0,a)',advance='no') &
            'Found same-flavour amplitude relation: ',i,' ='
       do idau=1,nterms
          write (99,'(1x,i0)',advance='no') pgl%same_flavour(idau,i)
          if (idau.lt.nterms) write (99,'(a)',advance='no') ' +'
       enddo
       write (99,*)
       pgl%amps(1)%same_flav(i)=.true.
       do ii=pgl%amps(1)%iproc_start(i),pgl%amps(1)%iproc_start(i+1)-1
          do idau=1,nterms
             jj=matching_helicity(pgl%amps(1),pgl%same_flavour(idau,i),ii)
             pgl%amps(1)%same_flavour_sum(ii,idau)=jj
          enddo
       enddo
    enddo
  end subroutine find_same_flavour

  subroutine build_same_flavour_relations(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: target,candidate,nlines,target_mask,candidate_mask
    integer :: ncandidates,wanted_terms
    integer,dimension(pgl%nproc) :: candidate_processes,candidate_masks
    integer,dimension(max_same_flavour_terms) :: terms
    logical :: found

    do target=1,pgl%nproc
       nlines=pgl%amps(1)%n_qqbar(target)
       if (nlines.lt.2 .or. nlines.gt.3) cycle
       if (.not.has_repeated_external_quark(pgl,target)) cycle
       target_mask=fermion_flow_mask(pgl,target,target,nlines)
       if (popcnt(target_mask).lt.2) cycle

       ncandidates=0
       do candidate=1,target-1
          if (pgl%amps(1)%n_qqbar(candidate).ne.nlines) cycle
          if (.not.compatible_flavour_refinement(pgl,target,candidate)) cycle
          candidate_mask=fermion_flow_mask(pgl,target,candidate,nlines)
          if (candidate_mask.eq.0 .or. candidate_mask.eq.target_mask) cycle
          if (iand(candidate_mask,not(target_mask)).ne.0) cycle
          ncandidates=ncandidates+1
          candidate_processes(ncandidates)=candidate
          candidate_masks(ncandidates)=candidate_mask
       enddo
       if (ncandidates.lt.2) cycle

       ! Several flavour assignments can select the same flow subset.  Their
       ! electroweak amplitudes need not be interchangeable, particularly for
       ! charged currents.  Try the small set of exact covers and retain only
       ! one whose complex helicity amplitudes satisfy the relation.
       found=.false.
       terms=0
       do wanted_terms=2,max_same_flavour_terms
          call search_cover(1,0)
          if (found) exit
       enddo
       if (found) pgl%same_flavour(:,target)=terms
    enddo
  contains
    recursive subroutine search_cover(depth,coverage)
      implicit none
      integer,intent(in) :: depth,coverage
      integer :: first_bit,icandidate,mask,new_coverage
      if (found) return
      do first_bit=0,max_same_flavour_terms-1
         if (btest(target_mask,first_bit).and..not.btest(coverage,first_bit)) exit
      enddo
      do icandidate=1,ncandidates
         mask=candidate_masks(icandidate)
         if (.not.btest(mask,first_bit)) cycle
         if (iand(mask,coverage).ne.0) cycle
         new_coverage=ior(coverage,mask)
         terms(depth)=candidate_processes(icandidate)
         if (depth.eq.wanted_terms) then
            if (new_coverage.ne.target_mask) cycle
            terms(depth+1:max_same_flavour_terms)=0
            if (relation_holds(pgl,target,terms)) then
               found=.true.
               return
            endif
         elseif (new_coverage.ne.target_mask) then
            call search_cover(depth+1,new_coverage)
            if (found) return
         endif
      enddo
      terms(depth)=0
    end subroutine search_cover
  end subroutine build_same_flavour_relations

  logical function has_repeated_external_quark(pgl,iproc)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: iproc
    integer :: i,j,pi,pj
    has_repeated_external_quark=.false.
    do i=1,pgl%next-1
       pi=outgoing_particle(pgl,i,iproc)
       if (.not.(phys_model%is_quark(pi).or.phys_model%is_antiquark(pi))) cycle
       do j=i+1,pgl%next
          pj=outgoing_particle(pgl,j,iproc)
          if (pi.eq.pj) then
             has_repeated_external_quark=.true.
             return
          endif
       enddo
    enddo
  end function has_repeated_external_quark

  logical function compatible_flavour_refinement(pgl,target,candidate)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: target,candidate
    integer :: i,pt,pc
    logical :: pure_qcd
    pure_qcd=process_is_pure_qcd(pgl,target)
    compatible_flavour_refinement=.false.
    do i=1,pgl%next
       pt=outgoing_particle(pgl,i,target)
       pc=outgoing_particle(pgl,i,candidate)
       if (phys_model%is_quark(pt)) then
          if (.not.phys_model%is_quark(pc)) return
          if (.not.pure_qcd .and. mod(abs(pt),2).ne.mod(abs(pc),2)) return
       elseif (phys_model%is_antiquark(pt)) then
          if (.not.phys_model%is_antiquark(pc)) return
          if (.not.pure_qcd .and. mod(abs(pt),2).ne.mod(abs(pc),2)) return
       elseif (pt.ne.pc) then
          return
       endif
       if (phys_model%get_mass(pt).ne.phys_model%get_mass(pc)) return
       if (phys_model%get_spin(pt).ne.phys_model%get_spin(pc)) return
    enddo
    compatible_flavour_refinement=.true.
  end function compatible_flavour_refinement

  logical function process_is_pure_qcd(pgl,iproc)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: iproc
    integer :: i,p
    process_is_pure_qcd=.true.
    do i=1,pgl%next
       p=outgoing_particle(pgl,i,iproc)
       if (phys_model%is_quark(p).or.phys_model%is_antiquark(p).or.&
            phys_model%is_gluon(p)) cycle
       process_is_pure_qcd=.false.
       return
    enddo
  end function process_is_pure_qcd

  integer function number_of_charged_current_insertions(pgl,iproc)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: iproc
    integer :: i,p,charged_lepton
    integer :: positive_charged,negative_charged,positive_neutrino,negative_neutrino
    number_of_charged_current_insertions=0
    do i=1,pgl%next
       p=abs(outgoing_particle(pgl,i,iproc))
       if (p.eq.24) number_of_charged_current_insertions=&
            number_of_charged_current_insertions+1
    enddo
    ! A resolved W decay is represented by an oppositely signed charged
    ! lepton/neutrino pair from the same generation.  Count those pairs rather
    ! than treating every neutrino as charged-current evidence: nu/anti-nu is
    ! also a perfectly ordinary neutral-current final state.
    do charged_lepton=11,15,2
       positive_charged=0
       negative_charged=0
       positive_neutrino=0
       negative_neutrino=0
       do i=1,pgl%next
          p=outgoing_particle(pgl,i,iproc)
          if (p.eq.charged_lepton) positive_charged=positive_charged+1
          if (p.eq.-charged_lepton) negative_charged=negative_charged+1
          if (p.eq.charged_lepton+1) positive_neutrino=positive_neutrino+1
          if (p.eq.-(charged_lepton+1)) negative_neutrino=negative_neutrino+1
       enddo
       number_of_charged_current_insertions=&
            number_of_charged_current_insertions+&
            min(positive_charged,negative_neutrino)+&
            min(negative_charged,positive_neutrino)
    enddo
  end function number_of_charged_current_insertions

  integer function outgoing_particle(pgl,ileg,iproc)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: ileg,iproc
    outgoing_particle=pgl%processes(ileg,iproc)
    if (ileg.le.2) outgoing_particle=phys_model%get_antipart(outgoing_particle)
  end function outgoing_particle

  integer function fermion_flow_mask(pgl,target,candidate,nlines)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: target,candidate,nlines
    integer :: i,iperm,iq,iaq,nflips,ncharged_currents
    integer,dimension(3) :: qlegs,aqlegs
    integer,dimension(3,6),parameter :: permutations=reshape([&
         1,2,3, 2,1,3, 1,3,2, 2,3,1, 3,1,2, 3,2,1],[3,6])
    logical :: charged_current,allowed

    qlegs=0
    aqlegs=0
    iq=0
    iaq=0
    do i=1,pgl%next
       if (phys_model%is_quark(outgoing_particle(pgl,i,target))) then
          iq=iq+1
          qlegs(iq)=i
       elseif (phys_model%is_antiquark(outgoing_particle(pgl,i,target))) then
          iaq=iaq+1
          aqlegs(iaq)=i
       endif
    enddo
    if (iq.ne.nlines .or. iaq.ne.nlines) then
       fermion_flow_mask=0
       return
    endif
    ncharged_currents=number_of_charged_current_insertions(pgl,target)
    charged_current=ncharged_currents.gt.0
    fermion_flow_mask=0
    do iperm=1,merge(2,6,nlines.eq.2)
       allowed=.true.
       nflips=0
       do i=1,nlines
          iq=abs(outgoing_particle(pgl,qlegs(i),candidate))
          iaq=abs(outgoing_particle(pgl,aqlegs(permutations(i,iperm)),candidate))
          if (charged_current) then
             if ((iq+1)/2.ne.(iaq+1)/2) allowed=.false.
             if (mod(iq,2).ne.mod(iaq,2)) nflips=nflips+1
          elseif (iq.ne.iaq) then
             allowed=.false.
          endif
       enddo
       ! A charged-current insertion changes weak isospin each time it is
       ! attached to a fermion line.  With one W (explicit or decayed to
       ! leptons), for example, flows containing three simultaneous flavour
       ! flips are not diagrams of this amplitude even when every endpoint
       ! belongs to the same generation.
       if (ncharged_currents.gt.0) then
          if (nflips.gt.ncharged_currents .or. &
               mod(nflips,2).ne.mod(ncharged_currents,2)) allowed=.false.
       endif
       if (allowed) fermion_flow_mask=ibset(fermion_flow_mask,iperm-1)
    enddo
  end function fermion_flow_mask

  logical function relation_holds(pgl,target,terms)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: target
    integer,dimension(:),intent(in) :: terms
    integer :: ii,idau,jj
    complex(kind=8) :: rhs
    real(kind=8) :: scale,error
    real(kind=8),parameter :: relative_tolerance=1d-9
    logical :: has_signal
    relation_holds=.false.
    has_signal=.false.
    do ii=pgl%amps(1)%iproc_start(target),pgl%amps(1)%iproc_start(target+1)-1
       rhs=(0d0,0d0)
       scale=abs(pgl%amps(1)%amps(ii))
       do idau=1,size(terms)
          if (terms(idau).le.0) exit
          jj=matching_helicity(pgl%amps(1),terms(idau),ii)
          if (jj.le.0) cycle
          rhs=rhs+pgl%amps(1)%amps(jj)
          scale=scale+abs(pgl%amps(1)%amps(jj))
       enddo
       error=abs(pgl%amps(1)%amps(ii)-rhs)
       if (scale.gt.0d0 .and. error.gt.relative_tolerance*scale) then
          relation_holds=.false.
          return
       endif
       if (scale.gt.0d0) has_signal=.true.
    enddo
    relation_holds=has_signal
  end function relation_holds

  integer function matching_helicity(amp,iproc,reference)
    implicit none
    type(amplitude_QCD),intent(in) :: amp
    integer,intent(in) :: iproc,reference
    integer :: candidate
    matching_helicity=0
    do candidate=amp%iproc_start(iproc),amp%iproc_start(iproc+1)-1
       if (all(amp%spins(:,1,reference).eq.amp%spins(:,1,candidate))) then
          matching_helicity=candidate
          return
       endif
    enddo
  end function matching_helicity

  subroutine setup_spin(pgl)
    ! Use the first process in the processes() array to setup all the possible
    ! spin states. Note that this assumes that all the processes() have the
    ! same number of spin states
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,iproc
    if (.not. allocated(pgl%spin)) allocate(pgl%spin(0:3,1:next))
    do i=1,next
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
       do i=1,next
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
    integer :: i,iproc,nq,ng,nsing,iq,iaq,is,ig,naq
    do iproc=1,pgl_unique%nproc
       nq=0
       ng=0
       nsing=0
       do i=1,next
          if (phys_model%is_quark(abs(pgl_unique%processes(i,iproc)))) then
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
          do i=1,next
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
          do i=1,next
             pgl_unique%color_orders(i,iproc)=i
          enddo
       elseif (nq.eq.2) then
          ig=2
          is=ng+2
          do i=1,next
             if (phys_model%is_quark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(1,iproc)=i
             elseif (phys_model%is_antiquark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(next,iproc)=i
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
          do i=1,next
             if (phys_model%is_quark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(iq,iproc)=i
                iq=iq+2
             elseif (phys_model%is_antiquark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(iaq,iproc)=i
                if (iaq.eq.nq-2) then
                   iaq=next
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
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,iproc
    do iproc=1,pgl%nproc
       do i=1,2
          if (pgl%processes(i,iproc).eq.21) then
             ! gluon: two polarisations and 8 colours
             pgl%iden(iproc)=pgl%iden(iproc)*2*8
          elseif (abs(pgl%processes(i,iproc)).ge.1 .and. abs(pgl%processes(i,iproc)).le.6) then
             ! (anti-)quark: two helicities and 3 colours
             pgl%iden(iproc)=pgl%iden(iproc)*2*3
          else
             ! assume two spin states and no colour:
             pgl%iden(iproc)=pgl%iden(iproc)*2
          endif
       enddo
    enddo
  end subroutine set_initial_state_average_factor

  subroutine set_final_state_identical_particle_factor(pgl)
    use math_functions
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,j,ni,iproc
    integer,dimension(:,:),allocatable :: iden_part
    allocate(iden_part(1:next,2))
    do iproc=1,pgl%nproc
       ni=0
       do i=3,next
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
          pgl%iden(iproc)=pgl%iden(iproc)*factorial8(iden_part(i,2))
       enddo
    enddo
    deallocate(iden_part)
  end subroutine set_final_state_identical_particle_factor

  subroutine compute_LC_colour_factor(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,ifac,iproc
    real(kind=8) :: fac
    do iproc=1,pgl%nproc
       fac=0d0
       do i=1,next
          if (pgl%processes(i,iproc).eq.21) then
             fac=fac+1d0
          elseif (abs(pgl%processes(i,iproc)).ge.1 .and. abs(pgl%processes(i,iproc)).le.6) then
             fac=fac+0.5d0
          endif
       enddo
       ifac=nint(fac)
       if (dble(ifac).ne.fac) then
          write (*,*) 'There is some issue with the LC colour factor computation: '// &
               'colour factor is not an integer',ifac,fac
          stop 1
       endif
       pgl%col_fac(iproc)=3**ifac
    enddo
  end subroutine compute_LC_colour_factor

  subroutine define_identical_procs(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: iproc,ip,n
    ! first fill the number of identical processes per iproc (so that we can
    ! allocate the array with the right size)
    allocate(pgl%iden_iproc(1:pgl%nproc))
    do iproc=1,pgl%nproc
       pgl%iden_iproc(iproc)=1
       if (any(abs(pgl%processes(1:next,iproc)).eq.1)) then
          pgl%iden_iproc(iproc)=pgl%iden_iproc(iproc)*5
       endif
       if (any(abs(pgl%processes(1:next,iproc)).eq.2)) then
          pgl%iden_iproc(iproc)=pgl%iden_iproc(iproc)*4
       endif
    enddo
    allocate(pgl%val_procs(1:maxval(pgl%iden_iproc(1:pgl%nproc)),1:pgl%nproc))
    allocate(pgl%iden_processes(1:next,1:maxval(pgl%iden_iproc(1:pgl%nproc)),1:pgl%nproc))
    ! Loop again and actually fill the iden_processes()
    do iproc=1,pgl%nproc
       do ip=0,pgl%iden_iproc(iproc)-1
          do n=1,next
             if (abs(pgl%processes(n,iproc)).eq.1) then
                pgl%iden_processes(n,ip+1,iproc)=sign(mod(ip,5)+1,pgl%processes(n,iproc))
             elseif (abs(pgl%processes(n,iproc)).eq.2) then
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
    integer :: ngl=0,ngl_tot=0,n_sing=0
    integer,dimension(6) :: nq,naq
    integer :: i,j
    integer(kind=8) :: tot_ord
    integer :: ic,i_qq,iaq,i_ini,i_inv,i_swap,k,l,ip
    integer,dimension(next) :: fgluons,ips,ips_out
    logical :: same_flavour
    logical,dimension(next) :: fgluon
    integer,dimension(:,:),allocatable :: io_list,io

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
       tot_ord=factorial8(ngl_tot-1)
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
          sym_fac=2*factorial8(ngl)
       endif
    elseif (nquarks.eq.2) then
       tot_ord=factorial8(ngl_tot)
       if ((abs(part(1)).ge.1 .and. abs(part(1)).le.6) .and. &
            (abs(part(2)).ge.1 .and. abs(part(2)).le.6) )then
          ! quark and anti-quark are incoming. Only 1 channel needed,
          ! which would result in the following symmetry factor:
          sym_fac=factorial8(ngl)
       elseif ((abs(part(1)).ge.1 .and. abs(part(1)).le.6) .or. &
            (abs(part(2)).ge.1 .and. abs(part(2)).le.6) )then
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
          sym_fac=factorial8(ngl)*2
          if (ngl.gt.0) then
             sym_fac=factorial8(ngl)*2
          endif
          if (ifindloc(o,next,1).ne.next-ifindloc(o,next,2)+1) then
             sym_fac=sym_fac*2
          endif
       endif
    elseif (nquarks.eq.4) then
       ! total number of potentially different orders: two ways of connecting
       ! quarks, (n-4)! orderings for the gluons, n-3 ways for an order to
       ! distribute the gluons among the two quark lines
       tot_ord=2*factorial8(ngl_tot)*(ngl_tot+1)
       n_sing=next-4-ngl_tot

       do i=2,next-1
          if ((o(i).gt.2 .and. part(o(i)).le.-1 .and. part(o(i)).ge.-6) .or. &
               (o(i).le.2 .and. part(o(i)).ge. 1 .and. part(o(i)).le. 6) ) then
             if (abs(part(o(i))).eq.abs(part(o(1))) .and. abs(part(o(i))).eq.abs(part(o(next)))) then
                same_flavour=.true.
             else
                same_flavour=.false.
             endif
             exit
          endif
       enddo

       allocate(io_list(1:next,tot_ord))
       allocate(io(1:next,5))
       ic=0
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
                   do ip=1,factorial(k)
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
       sym_fac=ic
    else        
       write (*,*) 'WARNING: symmetry factor missing',nquarks
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
  end subroutine finalize_multichan_info

  subroutine finalize_phase_space_order_group(pgl)
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i
    do i=1,size(pgl%amps)
       call finalize_amplitude_QCD(pgl%amps(i))
    enddo
    if (allocated(pgl%amps)) deallocate(pgl%amps)
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
    if (allocated(pgl%val_procs)) deallocate(pgl%val_procs)
    if (allocated(pgl%idenCOandMAPfactor)) deallocate(pgl%idenCOandMAPfactor)
    if (allocated(pgl%iden_processes)) deallocate(pgl%iden_processes)
    if (allocated(pgl%spin)) deallocate(pgl%spin)
    if (allocated(pgl%iden)) deallocate(pgl%iden)
    if (allocated(pgl%col_fac)) deallocate(pgl%col_fac)
    if (allocated(pgl%amp2)) deallocate(pgl%amp2)
    if (allocated(pgl%amp2_hel)) deallocate(pgl%amp2_hel)
    if (allocated(pgl%hel)) deallocate(pgl%hel)
    if (allocated(pgl%hel_fac)) deallocate(pgl%hel_fac)
    if (allocated(pgl%include_hel)) deallocate(pgl%include_hel)
    if (allocated(pgl%pT_min)) deallocate(pgl%pT_min)
    if (allocated(pgl%eta_max)) deallocate(pgl%eta_max)
    if (allocated(pgl%DR_min)) deallocate(pgl%DR_min)
    if (allocated(pgl%sqrt_s_min)) deallocate(pgl%sqrt_s_min)
    if (allocated(pgl%same_flavour)) deallocate(pgl%same_flavour)
  end subroutine finalize_phase_space_order_group
end module handling_processes
