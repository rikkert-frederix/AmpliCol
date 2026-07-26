module integration_analysis
  ! Process-independent parton-level SM analysis.
  !
  ! Jets use the same inclusive longitudinally invariant kT clustering,
  ! E-scheme recombination, radius DRjj_min, and pT/eta selection as the
  ! subtraction cuts.  Charged leptons and photons are bare final-state
  ! particles.  Missing transverse momentum is the vector sum of final-state
  ! neutrinos.  Ranked objects are ordered by decreasing transverse momentum.
  use integration_histograms, only: histogram_book,histogram_fill
  use common, only: pTj_min,etaj_max
  use cuts, only: cluster_kt_jets
  implicit none
  private

  integer,parameter :: max_ranked_jets=6,max_ranked_leptons=2
  integer,parameter :: h_inclusive=1,h_njets=2,h_jet_pt=10,h_jet_eta=20
  integer,parameter :: h_ht_jets=30,h_mjj=31,h_drjj=32
  integer,parameter :: h_nleptons=40,h_lepton_pt=41,h_lepton_eta=43,h_mll=45
  integer,parameter :: h_nphotons=50,h_photon_pt=51,h_photon_eta=52,h_maa=53
  integer,parameter :: h_met=60,h_ht_visible=61,h_m_final=62
  integer,parameter :: h_nheavy=70,h_heavy_pt=71,h_heavy_y=72
  ! The default observables depend only on momenta and broad object classes,
  ! not on the flavour of a massless QCD parton.  Set this to true in a custom
  ! analysis that distinguishes, for example, u-, d-, or b-initiated terms.
  logical,parameter,public :: analysis_distinguishes_massless_qcd_flavours=.false.

  public :: analysis_begin,analysis_fill

contains

  subroutine analysis_begin()
    integer :: i
    character(len=96) :: title

    call histogram_book(h_inclusive,'inclusive cross section [pb/bin]',1,0d0,1d0)
    call histogram_book(h_njets,'selected jet multiplicity',15,-0.5d0,14.5d0)

    do i=1,max_ranked_jets
       write(title,'(a,i0,a)') 'jet ',i,' pT [GeV]'
       call histogram_book(h_jet_pt+i-1,trim(title),100,0d0,2000d0)
       write(title,'(a,i0,a)') 'jet ',i,' eta'
       call histogram_book(h_jet_eta+i-1,trim(title),60,-6d0,6d0)
    enddo
    call histogram_book(h_ht_jets,'jet HT [GeV]',100,0d0,5000d0)
    call histogram_book(h_mjj,'leading dijet invariant mass [GeV]',100,0d0,5000d0)
    call histogram_book(h_drjj,'leading dijet deltaR',60,0d0,6d0)

    call histogram_book(h_nleptons,'charged-lepton multiplicity',9,-0.5d0,8.5d0)
    do i=1,max_ranked_leptons
       write(title,'(a,i0,a)') 'charged lepton ',i,' pT [GeV]'
       call histogram_book(h_lepton_pt+i-1,trim(title),100,0d0,1000d0)
       write(title,'(a,i0,a)') 'charged lepton ',i,' eta'
       call histogram_book(h_lepton_eta+i-1,trim(title),60,-6d0,6d0)
    enddo
    call histogram_book(h_mll,'leading dilepton invariant mass [GeV]',100,0d0,2000d0)

    call histogram_book(h_nphotons,'photon multiplicity',9,-0.5d0,8.5d0)
    call histogram_book(h_photon_pt,'leading photon pT [GeV]',100,0d0,1000d0)
    call histogram_book(h_photon_eta,'leading photon eta',60,-6d0,6d0)
    call histogram_book(h_maa,'leading diphoton invariant mass [GeV]',100,0d0,2000d0)

    call histogram_book(h_met,'missing transverse momentum [GeV]',100,0d0,1000d0)
    call histogram_book(h_ht_visible,'visible scalar HT [GeV]',100,0d0,5000d0)
    call histogram_book(h_m_final,'final-state invariant mass [GeV]',100,0d0,7000d0)

    call histogram_book(h_nheavy,'top/W/Z/H multiplicity',9,-0.5d0,8.5d0)
    call histogram_book(h_heavy_pt,'leading top/W/Z/H pT [GeV]',100,0d0,2000d0)
    call histogram_book(h_heavy_y,'leading top/W/Z/H rapidity',60,-6d0,6d0)
  end subroutine analysis_begin

  subroutine analysis_fill(nexternal,p,ipdg,wgt_nlo,wgt_born)
    integer,intent(in) :: nexternal
    real(kind=8),intent(in) :: p(0:3,nexternal)
    integer,intent(in) :: ipdg(nexternal)
    real(kind=8),intent(in) :: wgt_nlo,wgt_born
    real(kind=8) :: clustered(0:3,max(1,nexternal-2))
    real(kind=8) :: jets(0:3,max(1,nexternal-2))
    real(kind=8) :: leptons(0:3,max(1,nexternal-2))
    real(kind=8) :: photons(0:3,max(1,nexternal-2))
    real(kind=8) :: heavy(0:3,max(1,nexternal-2))
    real(kind=8) :: p_invisible(0:3),p_final(0:3)
    real(kind=8) :: ht_jets,ht_visible
    integer :: i,nclustered,njets,nleptons,nphotons,nheavy

    call histogram_fill(h_inclusive,0.5d0,wgt_nlo,wgt_born)

    call cluster_kt_jets(p,ipdg,clustered,nclustered)
    jets=0d0
    njets=0
    do i=1,nclustered
       if (transverse_momentum(clustered(:,i)).lt.pTj_min) cycle
       if (etaj_max.gt.0d0 .and. abs(pseudorapidity(clustered(:,i))).gt.etaj_max) cycle
       njets=njets+1
       jets(:,njets)=clustered(:,i)
    enddo
    call sort_by_pt(jets,njets)

    leptons=0d0
    photons=0d0
    heavy=0d0
    p_invisible=0d0
    p_final=0d0
    ht_visible=0d0
    nleptons=0
    nphotons=0
    nheavy=0
    do i=3,nexternal
       p_final=p_final+p(:,i)
       if (is_neutrino(ipdg(i))) then
          p_invisible=p_invisible+p(:,i)
       else
          ht_visible=ht_visible+transverse_momentum(p(:,i))
       endif
       if (is_charged_lepton(ipdg(i))) then
          nleptons=nleptons+1
          leptons(:,nleptons)=p(:,i)
       elseif (ipdg(i).eq.22) then
          nphotons=nphotons+1
          photons(:,nphotons)=p(:,i)
       elseif (is_heavy_sm_object(ipdg(i))) then
          nheavy=nheavy+1
          heavy(:,nheavy)=p(:,i)
       endif
    enddo
    call sort_by_pt(leptons,nleptons)
    call sort_by_pt(photons,nphotons)
    call sort_by_pt(heavy,nheavy)

    call histogram_fill(h_njets,dble(njets),wgt_nlo,wgt_born)
    ht_jets=0d0
    do i=1,njets
       ht_jets=ht_jets+transverse_momentum(jets(:,i))
       if (i.le.max_ranked_jets) then
          call histogram_fill(h_jet_pt+i-1,transverse_momentum(jets(:,i)),wgt_nlo,wgt_born)
          call histogram_fill(h_jet_eta+i-1,pseudorapidity(jets(:,i)),wgt_nlo,wgt_born)
       endif
    enddo
    call histogram_fill(h_ht_jets,ht_jets,wgt_nlo,wgt_born)
    if (njets.ge.2) then
       call histogram_fill(h_mjj,invariant_mass(jets(:,1)+jets(:,2)),wgt_nlo,wgt_born)
       call histogram_fill(h_drjj,delta_r(jets(:,1),jets(:,2)),wgt_nlo,wgt_born)
    endif

    call histogram_fill(h_nleptons,dble(nleptons),wgt_nlo,wgt_born)
    do i=1,min(nleptons,max_ranked_leptons)
       call histogram_fill(h_lepton_pt+i-1,transverse_momentum(leptons(:,i)),wgt_nlo,wgt_born)
       call histogram_fill(h_lepton_eta+i-1,pseudorapidity(leptons(:,i)),wgt_nlo,wgt_born)
    enddo
    if (nleptons.ge.2) then
       call histogram_fill(h_mll,invariant_mass(leptons(:,1)+leptons(:,2)),wgt_nlo,wgt_born)
    endif

    call histogram_fill(h_nphotons,dble(nphotons),wgt_nlo,wgt_born)
    if (nphotons.ge.1) then
       call histogram_fill(h_photon_pt,transverse_momentum(photons(:,1)),wgt_nlo,wgt_born)
       call histogram_fill(h_photon_eta,pseudorapidity(photons(:,1)),wgt_nlo,wgt_born)
    endif
    if (nphotons.ge.2) then
       call histogram_fill(h_maa,invariant_mass(photons(:,1)+photons(:,2)),wgt_nlo,wgt_born)
    endif

    call histogram_fill(h_met,transverse_momentum(p_invisible),wgt_nlo,wgt_born)
    call histogram_fill(h_ht_visible,ht_visible,wgt_nlo,wgt_born)
    call histogram_fill(h_m_final,invariant_mass(p_final),wgt_nlo,wgt_born)

    call histogram_fill(h_nheavy,dble(nheavy),wgt_nlo,wgt_born)
    if (nheavy.ge.1) then
       call histogram_fill(h_heavy_pt,transverse_momentum(heavy(:,1)),wgt_nlo,wgt_born)
       call histogram_fill(h_heavy_y,rapidity(heavy(:,1)),wgt_nlo,wgt_born)
    endif
  end subroutine analysis_fill

  subroutine sort_by_pt(objects,nobjects)
    real(kind=8),intent(inout) :: objects(0:,:)
    integer,intent(in) :: nobjects
    real(kind=8) :: tmp(0:3)
    integer :: i,j,imax
    do i=1,nobjects-1
       imax=i
       do j=i+1,nobjects
          if (transverse_momentum(objects(:,j)).gt.transverse_momentum(objects(:,imax))) imax=j
       enddo
       if (imax.ne.i) then
          tmp=objects(:,i)
          objects(:,i)=objects(:,imax)
          objects(:,imax)=tmp
       endif
    enddo
  end subroutine sort_by_pt

  real(kind=8) function transverse_momentum(momentum)
    real(kind=8),intent(in) :: momentum(0:3)
    transverse_momentum=sqrt(max(0d0,momentum(1)**2+momentum(2)**2))
  end function transverse_momentum

  real(kind=8) function pseudorapidity(momentum)
    real(kind=8),intent(in) :: momentum(0:3)
    real(kind=8) :: pabs,pplus,pminus
    pabs=sqrt(max(0d0,sum(momentum(1:3)**2)))
    pplus=pabs+momentum(3)
    pminus=pabs-momentum(3)
    if (pplus.le.tiny(1d0) .or. pminus.le.tiny(1d0)) then
       pseudorapidity=sign(huge(1d0),momentum(3))
    else
       pseudorapidity=0.5d0*log(pplus/pminus)
    endif
  end function pseudorapidity

  real(kind=8) function rapidity(momentum)
    real(kind=8),intent(in) :: momentum(0:3)
    real(kind=8) :: eplus,eminus
    eplus=momentum(0)+momentum(3)
    eminus=momentum(0)-momentum(3)
    if (eplus.le.tiny(1d0) .or. eminus.le.tiny(1d0)) then
       rapidity=sign(huge(1d0),momentum(3))
    else
       rapidity=0.5d0*log(eplus/eminus)
    endif
  end function rapidity

  real(kind=8) function invariant_mass(momentum)
    real(kind=8),intent(in) :: momentum(0:3)
    real(kind=8) :: mass2
    mass2=momentum(0)**2-sum(momentum(1:3)**2)
    invariant_mass=sqrt(max(0d0,mass2))
  end function invariant_mass

  real(kind=8) function delta_r(momentum1,momentum2)
    real(kind=8),intent(in) :: momentum1(0:3),momentum2(0:3)
    real(kind=8) :: dphi,deta,pi
    pi=acos(-1d0)
    dphi=atan2(momentum1(2),momentum1(1))-atan2(momentum2(2),momentum2(1))
    if (dphi.gt.pi) dphi=dphi-2d0*pi
    if (dphi.lt.-pi) dphi=dphi+2d0*pi
    deta=pseudorapidity(momentum1)-pseudorapidity(momentum2)
    delta_r=sqrt(deta**2+dphi**2)
  end function delta_r

  logical function is_charged_lepton(ipdg)
    integer,intent(in) :: ipdg
    is_charged_lepton=abs(ipdg).eq.11 .or. abs(ipdg).eq.13 .or. abs(ipdg).eq.15
  end function is_charged_lepton

  logical function is_neutrino(ipdg)
    integer,intent(in) :: ipdg
    is_neutrino=abs(ipdg).eq.12 .or. abs(ipdg).eq.14 .or. abs(ipdg).eq.16
  end function is_neutrino

  logical function is_heavy_sm_object(ipdg)
    integer,intent(in) :: ipdg
    is_heavy_sm_object=abs(ipdg).eq.6 .or. ipdg.eq.23 .or. abs(ipdg).eq.24 .or. ipdg.eq.25
  end function is_heavy_sm_object

end module integration_analysis
