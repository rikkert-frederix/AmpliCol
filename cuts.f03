module cuts
  use common
  use particles
  
  ! setup cuts. Set to '-1d0' means do not apply any cut on this variable
  ! jets:
  real(kind=8),parameter :: pTj_min      = 30d0
  real(kind=8),parameter :: DRjj_min     = 0.4d0  ! max allowed value: Drjj_min=1d0
  real(kind=8),parameter :: etaj_max     = 6d0
  real(kind=8),parameter :: sqrt_sjj_min = -1d0

  ! photons:
  real(kind=8),parameter :: pTa_min      = 30d0
  real(kind=8),parameter :: DRaa_min     = 0.4d0  ! max allowed value: Draa_min=1d0
  real(kind=8),parameter :: etaa_max     = 6d0
  real(kind=8),parameter :: sqrt_saa_min = -1d0

  ! jets+photons:
  real(kind=8),parameter :: DRja_min     = 0.4d0  ! max allowed value: Drja_min=1d0
  real(kind=8),parameter :: sqrt_sja_min = -1d0

contains
  logical function pass_cuts(pgl)
    ! Cuts on the phase-space point. Note that these cuts need to be symmetric
    ! under pz -> -pz.
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer :: i,j
    ! cuts on single particles
    pass_cuts=.true.
    do i=1,pgl%next
       if (pgl%pT_min(i).gt.0d0) then
          if (pt(pgl%phase_space%p(0,i)).lt.pgl%pT_min(i)) then
             pass_cuts=.false.
             return
          endif
       endif
       if (pgl%eta_max(i).gt.0d0) then
          if (abs(eta(pgl%phase_space%p(0,i))).gt.pgl%eta_max(i)) then
             pass_cuts=.false.
             return
          endif
       endif
    enddo
    ! cuts on pairs of particles
    do i=1,pgl%next-1
       do j=i+1,pgl%next
          if (pgl%sqrt_s_min(i,j).gt.0d0) then
             if (abs(sumdot(pgl%phase_space%p(0,i),pgl%phase_space%p(0,j))).lt.pgl%sqrt_s_min(i,j)**2) then
                pass_cuts=.false.
                return
             endif
          endif
          if (pgl%DR_min(i,j).gt.0d0) then
             if (abs(deltaR(pgl%phase_space%p(0,i),pgl%phase_space%p(0,j))).lt.pgl%DR_min(i,j)) then
                pass_cuts=.false.
                return
             endif
          endif
       enddo
    enddo
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


  subroutine setup_cuts_for_each_particle(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,j
    if (allocated(pgl%pT_min)) then
       write (*,*) 'ERROR: setting-up phase space cuts already'//&
            ' done for this phase-space group'
       stop 1
    endif
    ! check consistency among processes
    do i=1,pgl%next
       if (is_jet(pgl%processes(i,1))) then
          do j=2,pgl%nproc
             if (.not.is_jet(pgl%processes(i,j))) then
                write (*,*) 'inconsistent processes and cuts #1'
                stop 1
             endif
          enddo
       elseif(is_photon(pgl%processes(i,1))) then
          do j=2,pgl%nproc
             if (.not.is_photon(pgl%processes(i,j))) then
                write (*,*) 'inconsistent processes and cuts #2'
                stop 1
             endif
          enddo
       else
          do j=2,pgl%nproc
             if (is_jet(pgl%processes(i,j)) .or. is_photon(pgl%processes(i,j))) then
                write (*,*) 'inconsistent processes and cuts #3'
                stop 1
             endif
          enddo
       endif
    enddo
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
       if (.not. is_jet(pgl%processes(i,1))) cycle
       pgl%pT_min(i)=ptj_min
       pgl%eta_max(i)=etaj_max
    enddo
    ! cuts on single photons
    do i=3,pgl%next
       if (.not. is_photon(pgl%processes(i,1))) cycle
       pgl%pT_min(i)=pta_min
       pgl%eta_max(i)=etaa_max
    enddo
    ! cuts on pair of jets
    do i=1,pgl%next
       if (.not. is_jet(pgl%processes(i,1))) cycle
          do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not. is_jet(pgl%processes(j,1))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_sjj_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRjj_min
          endif
       enddo
    enddo
    ! cuts on pair of photons
    do i=1,pgl%next
       if (.not. is_photon(pgl%processes(i,1))) cycle
          do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not. is_photon(pgl%processes(j,1))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_saa_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRaa_min
          endif
       enddo
    enddo
    ! cuts on jet-photon pair
    do i=1,pgl%next
       do j=1,pgl%next
          if (i.eq.j) cycle
          if (.not.((is_jet(pgl%processes(i,1)) .and. is_photon(pgl%processes(j,1))) .or. &
                    (is_photon(pgl%processes(i,1)) .and. is_jet(pgl%processes(j,1))))) cycle
          pgl%sqrt_s_min(i,j)=sqrt_sja_min
          if (i.ge.3 .and. j.ge.3) then
             pgl%DR_min(i,j)=DRja_min
          endif
       enddo
    enddo
  end subroutine setup_cuts_for_each_particle
  
end module cuts
