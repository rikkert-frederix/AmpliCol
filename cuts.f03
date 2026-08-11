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
  
  real(kind=8) function pt(p)
    ! transverse momentum of 'p'
    implicit none
    real(kind=8), dimension(0:3) :: p
    pt=sqrt(p(1)**2+p(2)**2)
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
    real(kind=8) :: theta
    theta=acos(p(3)/sqrt(p(1)**2+p(2)**2+p(3)**2))
    eta=-log(dtan(theta/2d0))
  end function eta

  real(kind=8) function delta_phi(p1,p2)
    ! azimuthal difference of 'p1' and 'p2'
    implicit none
    real(kind=8), dimension(0:3) :: p1,p2
    real(kind=8) :: denom,arg
    real(kind=8),parameter :: tiny=1d-8
    denom=pt(p1)*pt(p2)
    arg=(p1(1)*p2(1)+p1(2)*p2(2))/denom
    if (arg.lt.-1d0-tiny) then
       write (*,*) 'cosine is complex'
       stop 1
    elseif (arg.lt.-1d0) then
       arg=-1d0
    elseif(arg.gt.1d0+tiny) then
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
