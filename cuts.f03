module cuts
  use common
  use particles
  use handling_processes
contains
  logical function pass_cuts(pgl)
    ! Cuts on the phase-space point. 
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer :: i,j
    ! cuts on single particles
    pass_cuts=.true.
    do i=1,pgl%next
       if (pgl%pT_min(i).gt.0d0) then
          if (pt(pgl%ps(1)%p(0,i)).lt.pgl%pT_min(i)) then
             pass_cuts=.false.
             return
          endif
       endif
       if (pgl%eta_max(i).gt.0d0) then
          if (abs(eta(pgl%ps(1)%p(0,i))).gt.pgl%eta_max(i)) then
             pass_cuts=.false.
             return
          endif
       endif
    enddo
    ! cuts on pairs of particles
    do i=1,pgl%next-1
       do j=i+1,pgl%next
          if (pgl%sqrt_s_min(i,j).gt.0d0) then
             if (abs(sumdot(pgl%ps(1)%p(0,i),pgl%ps(1)%p(0,j))).lt.pgl%sqrt_s_min(i,j)**2) then
                pass_cuts=.false.
                return
             endif
          endif
          if (pgl%DR_min(i,j).gt.0d0) then
             if (abs(deltaR(pgl%ps(1)%p(0,i),pgl%ps(1)%p(0,j))).lt.pgl%DR_min(i,j)) then
                pass_cuts=.false.
                return
             endif
          endif
       enddo
    enddo
!!$    if (abs(sumdot(pgl%ps(1)%p(0,pgl%next-1),pgl%ps(1)%p(0,pgl%next))).lt.50d0**2) then
!!$       pass_cuts=.false.
!!$       return
!!$    endif
  end function pass_cuts

  logical function pass_real_subtracted_cuts(pgl,iint)
    ! The real matrix element is measured with inclusive kT jets.  The
    ! required multiplicity is one below the real final-state jet count.
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: iint
    integer :: njet_required

    njet_required=real_subtracted_jet_requirement(pgl%processes(:,iint))
    pass_real_subtracted_cuts=pass_clustered_jet_cuts(pgl%ps(1)%p,pgl%processes(:,iint),njet_required)
  end function pass_real_subtracted_cuts

  integer function real_subtracted_jet_requirement(process) result(njet_required)
    implicit none
    integer,intent(in) :: process(:)
    integer :: i

    njet_required=0
    do i=3,size(process)
       if (phys_model%is_jet(process(i))) njet_required=njet_required+1
    enddo
    njet_required=max(0,njet_required-1)
  end function real_subtracted_jet_requirement

  logical function pass_mapped_dipole_cuts(p,process)
    ! A mapped dipole has the Born jet multiplicity, so every final-state
    ! jet parton is required to pass the jet cuts directly.
    implicit none
    real(kind=8),intent(in) :: p(0:,:)
    integer,intent(in) :: process(:)
    integer :: i,j,kind_i,kind_j

    pass_mapped_dipole_cuts=.false.
    if (size(p,2).ne.size(process)) return
    do i=3,size(process)
       kind_i=object_kind(process(i))
       if (kind_i.eq.1) then
          if (.not.pass_jet_cuts(p(:,i))) return
       elseif (.not.pass_nonjet_cuts(p(:,i),kind_i)) then
          return
       endif
    enddo
    do i=3,size(process)-1
       kind_i=object_kind(process(i))
       do j=i+1,size(process)
          kind_j=object_kind(process(j))
          if (.not.pass_object_pair_cuts(p(:,i),kind_i,p(:,j),kind_j,.true.)) return
       enddo
    enddo
    pass_mapped_dipole_cuts=.true.
  end function pass_mapped_dipole_cuts

  logical function pass_clustered_jet_cuts(p,process,njet_required)
    implicit none
    real(kind=8),intent(in) :: p(0:,:)
    integer,intent(in) :: process(:),njet_required
    real(kind=8) :: jets(0:3,max(1,size(process)-2))
    logical :: selected(max(1,size(process)-2))
    integer :: i,j,ijet,njets,nselected,kind_i,kind_j

    pass_clustered_jet_cuts=.false.
    if (size(p,2).ne.size(process)) return
    call cluster_kt_jets(p,process,jets,njets)
    selected=.false.
    nselected=0
    do ijet=1,njets
       if (pass_jet_cuts(jets(:,ijet))) then
          selected(ijet)=.true.
          nselected=nselected+1
       endif
    enddo
    if (nselected.lt.njet_required) return

    do i=3,size(process)
       kind_i=object_kind(process(i))
       if (kind_i.eq.1) cycle
       if (.not.pass_nonjet_cuts(p(:,i),kind_i)) return
    enddo
    do i=1,njets-1
       if (.not.selected(i)) cycle
       do j=i+1,njets
          if (.not.selected(j)) cycle
          ! The inclusive kT algorithm has already imposed the jet radius.
          if (.not.pass_object_pair_cuts(jets(:,i),1,jets(:,j),1,.false.)) return
       enddo
    enddo
    do ijet=1,njets
       if (.not.selected(ijet)) cycle
       do i=3,size(process)
          kind_i=object_kind(process(i))
          if (kind_i.eq.1) cycle
          if (.not.pass_object_pair_cuts(jets(:,ijet),1,p(:,i),kind_i,.true.)) return
       enddo
    enddo
    do i=3,size(process)-1
       kind_i=object_kind(process(i))
       if (kind_i.eq.1) cycle
       do j=i+1,size(process)
          kind_j=object_kind(process(j))
          if (kind_j.eq.1) cycle
          if (.not.pass_object_pair_cuts(p(:,i),kind_i,p(:,j),kind_j,.true.)) return
       enddo
    enddo
    pass_clustered_jet_cuts=.true.
  end function pass_clustered_jet_cuts

  subroutine cluster_kt_jets(p,process,jets,njets)
    ! Inclusive longitudinally invariant kT clustering with E-scheme sums.
    implicit none
    real(kind=8),intent(in) :: p(0:,:)
    integer,intent(in) :: process(:)
    real(kind=8),intent(out) :: jets(0:,:)
    integer,intent(out) :: njets
    real(kind=8) :: work(0:3,max(1,size(process)-2))
    real(kind=8) :: dmin,dbeam,dij,pti,ptj
    integer :: i,j,imin,jmin,nwork

    jets=0d0
    work=0d0
    nwork=0
    do i=3,size(process)
       if (.not.phys_model%is_jet(process(i))) cycle
       nwork=nwork+1
       work(:,nwork)=p(:,i)
    enddo
    njets=0
    do while (nwork.gt.0)
       dmin=huge(1d0)
       imin=0
       jmin=0
       do i=1,nwork
          pti=pt(work(:,i))
          if (pti .le. sqrt(tiny(1d0))) then
             dbeam=0d0
          else
             dbeam=pti*pti
          endif
          if (dbeam.lt.dmin) then
             dmin=dbeam
             imin=i
             jmin=0
          endif
       enddo
       do i=1,nwork-1
          pti=pt(work(:,i))
          if (pti .le. sqrt(tiny(1d0))) cycle
          do j=i+1,nwork
             ptj=pt(work(:,j))
             if (ptj .le. sqrt(tiny(1d0))) cycle
             dij=min(pti*pti,ptj*ptj)*deltaR(work(:,i),work(:,j))**2/DRjj_min**2
             if (dij.lt.dmin) then
                dmin=dij
                imin=i
                jmin=j
             endif
          enddo
       enddo
       if (jmin.eq.0) then
          njets=njets+1
          jets(:,njets)=work(:,imin)
          call remove_cluster(work,nwork,imin)
       else
          work(:,imin)=work(:,imin)+work(:,jmin)
          call remove_cluster(work,nwork,jmin)
       endif
    enddo
  end subroutine cluster_kt_jets

  real(kind=8) function real_subtracted_jet_pt_margin(pgl,iint) result(margin)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: iint
    real(kind=8) :: jets(0:3,max(1,pgl%next-2))
    integer :: njets,njet_required

    njet_required=real_subtracted_jet_requirement(pgl%processes(:,iint))
    call cluster_kt_jets(pgl%ps(1)%p,pgl%processes(:,iint),jets,njets)
    margin=jet_pt_cut_margin(jets,njets,njet_required)
  end function real_subtracted_jet_pt_margin

  real(kind=8) function mapped_dipole_jet_pt_margin(p,process) result(margin)
    implicit none
    real(kind=8),intent(in) :: p(0:,:)
    integer,intent(in) :: process(:)
    real(kind=8) :: jets(0:3,max(1,size(process)-2))
    integer :: i,njets

    jets=0d0
    njets=0
    do i=3,size(process)
       if (.not.phys_model%is_jet(process(i))) cycle
       njets=njets+1
       jets(:,njets)=p(:,i)
    enddo
    margin=jet_pt_cut_margin(jets,njets,njets)
  end function mapped_dipole_jet_pt_margin

  real(kind=8) function jet_pt_cut_margin(jets,njets,njet_required) result(margin)
    implicit none
    real(kind=8),intent(in) :: jets(0:,:)
    integer,intent(in) :: njets,njet_required
    real(kind=8) :: transverse_momenta(max(1,njets)),value
    integer :: i,j,neligible

    if (njet_required.le.0 .or. pTj_min.le.0d0) then
       margin=huge(1d0)
       return
    endif
    neligible=0
    do i=1,njets
       value=pt(jets(:,i))
       if (value.le.sqrt(tiny(1d0))) cycle
       if (etaj_max.gt.0d0) then
          if (abs(eta(jets(:,i))).ge.etaj_max) cycle
       endif
       neligible=neligible+1
       transverse_momenta(neligible)=value
    enddo
    if (neligible.lt.njet_required) then
       margin=-pTj_min
       return
    endif
    do i=2,neligible
       value=transverse_momenta(i)
       j=i-1
       do while (j.ge.1)
          if (transverse_momenta(j).ge.value) exit
          transverse_momenta(j+1)=transverse_momenta(j)
          j=j-1
       enddo
       transverse_momenta(j+1)=value
    enddo
    margin=transverse_momenta(njet_required)-pTj_min
  end function jet_pt_cut_margin

  subroutine remove_cluster(work,nwork,idx)
    implicit none
    real(kind=8),intent(inout) :: work(0:,:)
    integer,intent(inout) :: nwork
    integer,intent(in) :: idx
    if (idx.lt.nwork) work(:,idx:nwork-1)=work(:,idx+1:nwork)
    work(:,nwork)=0d0
    nwork=nwork-1
  end subroutine remove_cluster

  integer function object_kind(ipdg)
    implicit none
    integer,intent(in) :: ipdg
    if (phys_model%is_jet(ipdg)) then
       object_kind=1
    elseif (phys_model%is_photon(ipdg)) then
       object_kind=2
    elseif (phys_model%is_lepton_any(ipdg)) then
       object_kind=3
    else
       object_kind=0
    endif
  end function object_kind

  logical function pass_jet_cuts(p)
    implicit none
    real(kind=8),intent(in) :: p(0:3)
    pass_jet_cuts=pt(p).gt.pTj_min
    if (pass_jet_cuts .and. etaj_max.gt.0d0) pass_jet_cuts=abs(eta(p)).lt.etaj_max
  end function pass_jet_cuts

  logical function pass_nonjet_cuts(p,kind)
    implicit none
    real(kind=8),intent(in) :: p(0:3)
    integer,intent(in) :: kind
    pass_nonjet_cuts=.true.
    if (kind.eq.2) then
       if (pTa_min.gt.0d0 .and. pt(p).lt.pTa_min) pass_nonjet_cuts=.false.
       if (pass_nonjet_cuts .and. etaa_max.gt.0d0) pass_nonjet_cuts=abs(eta(p)).lt.etaa_max
    elseif (kind.eq.3) then
       if (pTl_min.gt.0d0 .and. pt(p).lt.pTl_min) pass_nonjet_cuts=.false.
       if (pass_nonjet_cuts .and. etal_max.gt.0d0) pass_nonjet_cuts=abs(eta(p)).lt.etal_max
    endif
  end function pass_nonjet_cuts

  logical function pass_object_pair_cuts(p1,kind1,p2,kind2,apply_jet_dr)
    implicit none
    real(kind=8),intent(in) :: p1(0:3),p2(0:3)
    integer,intent(in) :: kind1,kind2
    logical,intent(in) :: apply_jet_dr
    real(kind=8) :: drmin,smin

    pass_object_pair_cuts=.true.
    drmin=-1d0
    smin=-1d0
    if (kind1.eq.1 .and. kind2.eq.1) then
       smin=sqrt_sjj_min
       if (apply_jet_dr) drmin=DRjj_min
    elseif (kind1.eq.2 .and. kind2.eq.2) then
       smin=sqrt_saa_min
       drmin=DRaa_min
    elseif (kind1.eq.3 .and. kind2.eq.3) then
       smin=sqrt_sll_min
       drmin=DRll_min
    elseif ((kind1.eq.1 .and. kind2.eq.2) .or. (kind1.eq.2 .and. kind2.eq.1)) then
       smin=sqrt_sja_min
       drmin=DRja_min
    elseif ((kind1.eq.1 .and. kind2.eq.3) .or. (kind1.eq.3 .and. kind2.eq.1)) then
       smin=sqrt_sjl_min
       drmin=DRjl_min
    elseif ((kind1.eq.2 .and. kind2.eq.3) .or. (kind1.eq.3 .and. kind2.eq.2)) then
       smin=sqrt_sla_min
       drmin=DRla_min
    else
       return
    endif
    if (smin.gt.0d0 .and. abs(sumdot(p1,p2)).lt.smin**2) pass_object_pair_cuts=.false.
    if (pass_object_pair_cuts .and. drmin.gt.0d0) then
       pass_object_pair_cuts=deltaR(p1,p2).ge.drmin
    endif
  end function pass_object_pair_cuts
  
  real(kind=8) function pt(p)
    ! transverse momentum of 'p'
    implicit none
    real(kind=8), dimension(0:3) :: p
    real(kind=8) :: scale
    scale=max(abs(p(1)),abs(p(2)))
    if (scale .le. sqrt(tiny(1d0))) then
       pt=0d0
    else
       pt=scale*sqrt((p(1)/scale)**2+(p(2)/scale)**2)
    endif
  end function pt
  
  real(kind=8) function dot(p1,p2)
    ! Inner product between two 4-vectors
    implicit none
    real(kind=8),intent(in),dimension(0:3) :: p1,p2
    dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
  end function dot

  real(kind=8) function sumdot(p1,p2)
    ! Inner product between two 4-vectors
    implicit none
    real(kind=8),intent(in),dimension(0:3) :: p1,p2
    real(kind=8),dimension(0:3) :: p
    p=p1+p2
    sumdot=dot(p,p)
  end function sumdot

  real(kind=8) function eta(p)
    ! pseudo-rapidity of 'p'
    implicit none
    real(kind=8), dimension(0:3) :: p
    real(kind=8) :: pt_value
    pt_value=pt(p)
    if (pt_value .le. sqrt(tiny(1d0))) then
       eta=sign(huge(1d0),p(3))
    else
       ! asinh(pz/pT) is equivalent to pseudorapidity and remains stable
       ! for nearly beam-collinear momenta.
       eta=asinh(p(3)/pt_value)
    endif
  end function eta

  real(kind=8) function delta_phi(p1,p2)
    ! azimuthal difference of 'p1' and 'p2'
    implicit none
    real(kind=8), dimension(0:3) :: p1,p2
    real(kind=8) :: pt1,pt2,arg
    real(kind=8),parameter :: angle_tolerance=1d-8
    pt1=pt(p1)
    pt2=pt(p2)
    if (pt1 .le. sqrt(tiny(1d0)) .or. pt2 .le. sqrt(tiny(1d0))) then
       delta_phi=0d0
       return
    endif
    arg=(p1(1)/pt1)*(p2(1)/pt2)+(p1(2)/pt1)*(p2(2)/pt2)
    if (arg.lt.-1d0-angle_tolerance) then
       write (*,*) 'cosine is complex'
       stop 1
    elseif (arg.lt.-1d0) then
       arg=-1d0
    elseif(arg.gt.1d0+angle_tolerance) then
       write (*,*) 'cosine is complex'
       stop 1
    elseif(arg.gt.1d0) then
       arg=1d0
    endif
    delta_phi=acos(arg)
  end function delta_phi

  real(kind=8) function deltaR(p1,p2)
    ! Distance (Delta-R) between 'p1' and 'p2'
    implicit none
    real(kind=8), dimension(0:3) :: p1,p2
    deltaR=sqrt(delta_phi(p1,p2)**2+(eta(p1)-eta(p2))**2)
  end function deltaR

  subroutine setup_cuts_for_each_particle(pgl,ichan)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer,intent(in) :: ichan
    integer :: i,j,iproc,iident
    integer,dimension(:),allocatable :: reference_process
    if (allocated(pgl%pT_min)) then
       write (*,*) 'ERROR: setting-up phase space cuts already'//&
            ' done for this phase-space group'
       stop 1
    endif
    allocate(reference_process(pgl%next))
    if (ichan.gt.0 .and. allocated(pgl%iden_processes)) then
       ! Cuts act on the physical, target-labelled momenta.  Reduced matrix-
       ! element representatives may use a different flavour crossing, so
       ! derive the particle classes from the physical subprocess aliases.
       reference_process=pgl%iden_processes(:,1,1)
       do i=1,pgl%next
          do iproc=1,pgl%nproc
             do iident=1,pgl%iden_iproc(iproc)
                if (cut_particle_class(reference_process(i)).ne.&
                     cut_particle_class(pgl%iden_processes(i,iident,iproc))) then
                   write (*,*) 'inconsistent physical subprocesses and cuts',i,iproc,iident
                   stop 1
                endif
             enddo
          enddo
       enddo
    else
       reference_process=pgl%processes(:,1)
    endif
    ! initialize all:
    allocate(pgl%pT_min(1:pgl%next))
    allocate(pgl%eta_max(1:pgl%next))
    allocate(pgl%DR_min(1:pgl%next,1:pgl%next))
    allocate(pgl%sqrt_s_min(1:pgl%next,1:pgl%next))
    pgl%pT_min(1:pgl%next)=-1d0
    pgl%eta_max(1:pgl%next)=-1d0
    pgl%DR_min(1:pgl%next,1:pgl%next)=-1d0
    pgl%sqrt_s_min(1:pgl%next,1:pgl%next)=-1d0
    ! cuts on single jets
    do i=3,pgl%next
       if (.not. phys_model%is_jet(reference_process(i))) cycle
       pgl%pT_min(i)=ptj_min
       pgl%eta_max(i)=etaj_max
    enddo
    ! cuts on single photons
    do i=3,pgl%next
       if (.not. phys_model%is_photon(reference_process(i))) cycle
       pgl%pT_min(i)=pta_min
       pgl%eta_max(i)=etaa_max
    enddo
    ! cuts on single leptons
    do i=3,pgl%next
       if (.not. phys_model%is_lepton_any(reference_process(i))) cycle
       pgl%pT_min(i)=ptl_min
       pgl%eta_max(i)=etal_max
    enddo
    ! cuts on pair of jets
    do i=1,pgl%next
       if (.not. phys_model%is_jet(reference_process(i))) cycle
          do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not. phys_model%is_jet(reference_process(j))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_sjj_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRjj_min
          endif
       enddo
    enddo
    ! cuts on pair of photons
    do i=1,pgl%next
       if (.not. phys_model%is_photon(reference_process(i))) cycle
          do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not. phys_model%is_photon(reference_process(j))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_saa_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRaa_min
          endif
       enddo
    enddo
    ! cuts on pair of leptons
    do i=1,pgl%next
       if (.not. phys_model%is_lepton_any(reference_process(i))) cycle
          do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not. phys_model%is_lepton_any(reference_process(j))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_sll_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRll_min
          endif
       enddo
    enddo
    ! cuts on jet-photon pair
    do i=1,pgl%next
       do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not.((phys_model%is_jet(reference_process(i)) .and. &
               phys_model%is_photon(reference_process(j))) .or. &
               (phys_model%is_photon(reference_process(i)) .and. &
               phys_model%is_jet(reference_process(j))))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_sja_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRja_min
          endif
       enddo
    enddo
    ! cuts on jet-lepton pair
    do i=1,pgl%next
       do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not.((phys_model%is_jet(reference_process(i)) .and. &
               phys_model%is_lepton_any(reference_process(j))) .or. &
               (phys_model%is_lepton_any(reference_process(i)) .and. &
               phys_model%is_jet(reference_process(j))))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_sjl_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRjl_min
          endif
       enddo
    enddo
    ! cuts on lepton-photon pair
    do i=1,pgl%next
       do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not.((phys_model%is_lepton_any(reference_process(i)) .and. &
               phys_model%is_photon(reference_process(j))) .or. &
               (phys_model%is_photon(reference_process(i)) .and. &
               phys_model%is_lepton_any(reference_process(j))))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_sla_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRla_min
          endif
       enddo
    enddo
    write (99,*) '****************************************************'
    write (99,*) 'CUTS for channel',ichan
    do i=1,pgl%next
       write (99,*) i,'pT_min:',pgl%pT_min(i),'eta_max',pgl%eta_max(i)
    enddo
    do i=1,pgl%next-1
       do j=i+1,pgl%next
          write (99,*) i,j,'sqrt_s_min:',pgl%sqrt_s_min(i,j),'DR_min',pgl%DR_min(i,j)
       enddo
    enddo
    deallocate(reference_process)
    write (99,*) '****************************************************'
  end subroutine setup_cuts_for_each_particle

  integer function cut_particle_class(ipdg)
    implicit none
    integer,intent(in) :: ipdg
    if (phys_model%is_jet(ipdg)) then
       cut_particle_class=1
    elseif (phys_model%is_photon(ipdg)) then
       cut_particle_class=2
    elseif (phys_model%is_lepton_any(ipdg)) then
       cut_particle_class=3
    else
       cut_particle_class=0
    endif
  end function cut_particle_class
  
end module cuts
