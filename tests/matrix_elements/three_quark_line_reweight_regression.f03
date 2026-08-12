program three_quark_line_reweight_regression
  use amplitude_QCD_mod
  use particles
  implicit none
  integer,parameter :: dp=kind(1d0),n=6
  type(physics_model) :: model
  type(amplitude_QCD) :: amp
  integer,dimension(n,1) :: part,orders
  integer,dimension(0:3,n) :: spin
  integer,dimension(n) :: hel
  integer :: flow,leading_sector
  logical :: found_sparse_partner
  real(kind=dp),dimension(0:3,n) :: p
  real(kind=dp),dimension(3) :: matrix2,normalized
  real(kind=dp),parameter :: pi=3.14159265358979323846d0
  real(kind=dp),parameter :: alpha_s=0.118d0
  real(kind=dp),parameter :: alpha_ew=1d0/132.507d0
  ! Stock MadGraph 3.2.0 MATRIX(P,NHEL,IC=+1) for
  ! d d~ > u u~ s s~ QED=0 QCD=4 at this labelled fixed-helicity point.
  ! MATRIX omits the initial-state average; divide both results by 36.
  real(kind=dp),parameter :: madgraph_full=2.1057415991357259d-13
  real(kind=dp),parameter :: expected_lc=1.3039052156690463d-14

  open(unit=99,file='/dev/null',status='unknown',action='write')
  call model%init_part(173d0,0d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call model%init_vert()

  part(:,1)=[1,-1,2,-2,3,-3]
  orders=0
  hel=[-1,1,-1,1,-1,1]
  spin=0
  spin(0,:)=1
  spin(1,:)=-9
  call fill_momenta(p)

  call amp%init(2,n,1,part,spin,orders,model)
  call amp%init_col(n,20)
  if (amp%nColOrd.ne.6 .or. amp%n_amps.ne.6) then
     write (*,*) 'Unexpected three-line colour basis size:',amp%nColOrd,amp%n_amps
     stop 1
  endif
  call check_colour_factor_count(amp,1,27d0,6)
  call check_colour_factor_count(amp,3,27d0,6)
  call check_colour_factor_count(amp,3,-18d0,9)
  call check_colour_factor_count(amp,3,6d0,6)
  if (.not.allocated(amp%sector_term_start) .or.&
       .not.allocated(amp%sector_term_curr2amp) .or.&
       .not.allocated(amp%sector_term_sign)) then
     write (*,*) 'Three-line sparse coupling-sector closures are unavailable'
     stop 1
  endif
  if (size(amp%sector_term_curr2amp,1).ne.2 .or.&
       size(amp%sector_term_curr2amp,2).ne.size(amp%sector_term_sign) .or.&
       amp%sector_term_start(amp%n_amps,amp%n_sectors).ne.&
       size(amp%sector_term_sign)) then
     write (*,*) 'Three-line sparse coupling-sector metadata is inconsistent'
     stop 1
  endif
  leading_sector=maxloc(amp%sector_powers(1,1:amp%n_sectors),dim=1)
  found_sparse_partner=.false.
  do flow=1,amp%n_amps
     if (amp%sector_term_start(flow,leading_sector)-&
          amp%sector_term_start(flow-1,leading_sector).gt.1) &
          found_sparse_partner=.true.
  enddo
  if (.not.found_sparse_partner) then
     write (*,*) 'Three-line duplicate QCD closures were not retained sparsely'
     stop 1
  endif

  call amp%evaluate(n,p,hel,.false.,model)
  call check_lc_amplitudes(amp,p,hel)
  call colour_squared(amp,matrix2)
  normalized=matrix2*(4d0*pi*alpha_s)**4/36d0
  write (*,'(a,3es24.16)') 'THREE_LINE_LC_NLC_FC=',normalized
  if (abs(normalized(1)/expected_lc-1d0).gt.1d-11) then
     write (*,*) 'Three-line leading-colour regression mismatch:',normalized(1),expected_lc
     stop 1
  endif
  if (abs(normalized(3)/(madgraph_full/36d0)-1d0).gt.1d-11) then
     write (*,*) 'Three-line full-colour result disagrees with MadGraph:',&
          normalized(3),madgraph_full/36d0
     stop 1
  endif
  if (abs(matrix2(2)/matrix2(3)-1d0).gt.1d-13) then
     write (*,*) 'Three-line NLC and full-colour results unexpectedly differ:',matrix2(2:3)
     stop 1
  endif
  if (abs(matrix2(3)/matrix2(1)-1d0).lt.0.1d0) then
     write (*,*) 'Three-line colour regression is too close to leading colour'
     stop 1
  endif
  call check_one_gluon()
  call check_photon()
  call check_identical_quarks()
  call check_safe_multichannel_partner()
  write (*,'(a)') 'Three-quark-line reweight regression passed'

contains

  subroutine check_colour_factor_count(amplitude,accuracy,factor,expected_count)
    implicit none
    type(amplitude_QCD),intent(in) :: amplitude
    integer,intent(in) :: accuracy,expected_count
    real(kind=dp),intent(in) :: factor
    integer :: value_index

    do value_index=1,amplitude%n_col_vals(accuracy)
       if (amplitude%diff_col_vals(value_index,accuracy).eq.factor) then
          if (amplitude%row_index(amplitude%nColOrd,value_index,accuracy)&
               .ne.expected_count) then
             write (*,*) 'Unexpected colour-factor multiplicity:',accuracy,&
                  factor,amplitude%row_index(amplitude%nColOrd,value_index,accuracy),&
                  expected_count
             stop 1
          endif
          return
       endif
    enddo
    write (*,*) 'Missing colour factor:',accuracy,factor
    stop 1
  end subroutine check_colour_factor_count

  subroutine check_lc_amplitudes(all_colour_amp,momenta,helicities)
    implicit none
    type(amplitude_QCD),intent(in) :: all_colour_amp
    real(kind=dp),dimension(0:3,n),intent(in) :: momenta
    integer,dimension(n),intent(in) :: helicities
    type(amplitude_QCD),allocatable :: lc_amp
    integer :: flow
    integer,dimension(n,1) :: lc_part,lc_order
    integer,dimension(0:3,n) :: lc_spin
    integer :: all_colour_sector,lc_sector
    real(kind=dp) :: scale
    complex(kind=dp) :: lc_value

    lc_part(:,1)=part(:,1)
    lc_spin=0
    lc_spin(0,:)=1
    lc_spin(1,:)=helicities
    all_colour_sector=leading_qcd_sector(all_colour_amp)
    do flow=1,all_colour_amp%nColOrd
       lc_order(:,1)=all_colour_amp%perm(:,flow)
       allocate(lc_amp)
       call lc_amp%init(1,n,1,lc_part,lc_spin,lc_order,model)
       call lc_amp%evaluate(n,momenta,helicities,.false.,model)
       if (lc_amp%n_amps.lt.1) then
          write (*,*) 'Missing fixed-order amplitude:',flow
          stop 1
       endif
       lc_sector=leading_qcd_sector(lc_amp)
       lc_value=sum(lc_amp%amps_by_order(:,lc_sector))
       scale=max(1d-30,abs(all_colour_amp%amps_by_order(flow,all_colour_sector)))
       if (abs(lc_value-&
            all_colour_amp%amps_by_order(flow,all_colour_sector)).gt.1d-11*scale) then
          write (*,*) 'Grouped imode=2 flow disagrees with imode=1:',flow,&
               lc_value,&
               all_colour_amp%amps_by_order(flow,all_colour_sector)
          stop 1
       endif
       deallocate(lc_amp)
    enddo
  end subroutine check_lc_amplitudes

  subroutine check_one_gluon()
    implicit none
    integer,parameter :: n7=7
    type(amplitude_QCD) :: amp7
    integer,dimension(n7,1) :: part7,orders7
    integer,dimension(0:3,n7) :: spin7
    integer,dimension(n7) :: hel7
    real(kind=dp),dimension(0:3,n7) :: p7
    real(kind=dp),dimension(3) :: raw7,normalized7
    real(kind=dp),parameter :: expected_lc7=1.7878308543514494d-12
    ! Stock MadGraph 3.2.0 MATRIX result for QED=0 QCD=5.
    real(kind=dp),parameter :: madgraph_full7=2.9863921562813620d-11

    part7(:,1)=[1,-1,2,-2,3,-3,21]
    orders7=0
    hel7=[-1,1,-1,1,-1,1,-1]
    spin7=0
    spin7(0,:)=1
    spin7(1,:)=-9
    p7(:,1)=[500d0,0d0,0d0,500d0]
    p7(:,2)=[500d0,0d0,0d0,-500d0]
    p7(:,3)=[80.89693031749577d0,-35.363962752382264d0,&
         25.38917257558109d0,-68.18426056773842d0]
    p7(:,4)=[251.36895499042657d0,-179.6909561218757d0,&
         131.2696963946096d0,-116.90072125291724d0]
    p7(:,5)=[344.0807935954407d0,219.66991720940467d0,&
         -259.57359242327067d0,52.519235628094506d0]
    p7(:,6)=[131.41282860492473d0,125.97195147403914d0,&
         -1.2352831219401352d0,37.401511191104284d0]
    p7(:,7)=[192.24049249171213d0,-130.58694980918583d0,&
         104.15000657501999d0,95.16423500145686d0]

    call amp7%init(2,n7,1,part7,spin7,orders7,model)
    call amp7%init_col(n7,20)
    if (amp7%nColOrd.ne.18 .or. amp7%n_amps.ne.18) then
       write (*,*) 'Unexpected three-line plus gluon basis size:',&
            amp7%nColOrd,amp7%n_amps
       stop 1
    endif
    call check_colour_factor_count(amp7,1,81d0,18)
    call check_colour_factor_count(amp7,3,72d0,18)
    call check_colour_factor_count(amp7,3,-48d0,45)
    call check_colour_factor_count(amp7,3,16d0,54)
    call amp7%evaluate(n7,p7,hel7,.false.,model)
    call colour_squared(amp7,raw7)
    normalized7=raw7*(4d0*pi*alpha_s)**5/36d0
    write (*,'(a,3es24.16)') 'THREE_LINE_G_LC_NLC_FC=',normalized7
    if (abs(normalized7(1)/expected_lc7-1d0).gt.1d-11 .or.&
         abs(normalized7(3)/(madgraph_full7/36d0)-1d0).gt.1d-11) then
       write (*,*) 'Three-line plus gluon regression mismatch:',normalized7
       stop 1
    endif
  end subroutine check_one_gluon

  subroutine check_photon()
    implicit none
    integer,parameter :: n7=7
    type(amplitude_QCD) :: amp7
    integer,dimension(n7,1) :: part7,orders7
    integer,dimension(0:3,n7) :: spin7
    integer,dimension(n7) :: hel7
    real(kind=dp),dimension(0:3,n7) :: p7
    real(kind=dp),dimension(3) :: raw7,normalized7
    ! Stock MadGraph 3.2.0 MATRIX result for QED=1 QCD=4.
    real(kind=dp),parameter :: madgraph_full7=3.3601672982535744d-13

    part7(:,1)=[1,-1,2,-2,3,-3,22]
    orders7=0
    hel7=[-1,1,-1,1,-1,1,-1]
    spin7=0
    spin7(0,:)=1
    spin7(1,:)=-9
    p7(:,1)=[500d0,0d0,0d0,500d0]
    p7(:,2)=[500d0,0d0,0d0,-500d0]
    p7(:,3)=[80.89693031749577d0,-35.363962752382264d0,&
         25.38917257558109d0,-68.18426056773842d0]
    p7(:,4)=[251.36895499042657d0,-179.6909561218757d0,&
         131.2696963946096d0,-116.90072125291724d0]
    p7(:,5)=[344.0807935954407d0,219.66991720940467d0,&
         -259.57359242327067d0,52.519235628094506d0]
    p7(:,6)=[131.41282860492473d0,125.97195147403914d0,&
         -1.2352831219401352d0,37.401511191104284d0]
    p7(:,7)=[192.24049249171213d0,-130.58694980918583d0,&
         104.15000657501999d0,95.16423500145686d0]

    call amp7%init(2,n7,1,part7,spin7,orders7,model)
    call amp7%init_col(n7,20)
    if (amp7%nColOrd.ne.6 .or. amp7%n_amps.ne.6) then
       write (*,*) 'Unexpected three-line plus photon basis size:',&
            amp7%nColOrd,amp7%n_amps
       stop 1
    endif
    call amp7%evaluate(n7,p7,hel7,.false.,model)
    call colour_squared(amp7,raw7)
    normalized7=raw7*(4d0*pi*alpha_s)**4*&
         (2d0*4d0*pi*alpha_ew)/36d0
    write (*,'(a,3es24.16)') 'THREE_LINE_A_LC_NLC_FC=',normalized7
    if (abs(normalized7(3)/(madgraph_full7/36d0)-1d0).gt.1d-11) then
       write (*,*) 'Three-line plus photon result disagrees with MadGraph:',&
            normalized7(3),madgraph_full7/36d0
       stop 1
    endif
  end subroutine check_photon

  subroutine check_identical_quarks()
    implicit none
    ! These MadGraph MATRIX values are unaveraged. Dividing by 144 applies
    ! the 36-state initial average and 2! factors for the identical final
    ! quarks and antiquarks.
    call check_quark_flavour_case([1,-1,2,-2,2,-2],&
         1.9288834655530365d-12,144d0,'two identical final quark lines')
    call check_quark_flavour_case([2,-2,2,-2,2,-2],&
         5.6827348236142209d-11,144d0,'three identical quark lines')
    ! Two initial quarks exercise the other initial-state crossing branch.
    ! Here 72 is the initial average times 2! for the two final u quarks.
    call check_quark_flavour_case([2,2,2,2,1,-1],&
         9.5640936205091913d-14,72d0,'two initial quarks')
  end subroutine check_identical_quarks

  subroutine check_safe_multichannel_partner()
    implicit none
    type(amplitude_QCD) :: first_amp,second_amp
    integer,dimension(n,1) :: partner_part,first_order,second_order
    integer,dimension(0:3,n) :: partner_spin
    integer,dimension(n) :: partner_hel
    real(kind=dp),dimension(0:3,n) :: partner_p
    real(kind=dp) :: first_squared,second_squared,scale

    ! These are the two phase-space representatives paired by
    ! process_list.py for u u > u u d d~.  Their exact coloured-leg string
    ! contents are equal, so the imode=1 coefficients must agree at the same
    ! labelled momentum point without permuting any external legs.
    partner_part(:,1)=[2,2,2,2,1,-1]
    first_order(:,1)=[4,1,3,2,5,6]
    second_order(:,1)=[4,1,5,6,3,2]
    partner_hel=[-1,1,-1,1,-1,1]
    partner_spin=0
    partner_spin(0,:)=1
    partner_spin(1,:)=partner_hel
    call fill_momenta(partner_p)

    call first_amp%init(1,n,1,partner_part,partner_spin,first_order,model)
    call second_amp%init(1,n,1,partner_part,partner_spin,second_order,model)
    call first_amp%evaluate(n,partner_p,partner_hel,.false.,model)
    call second_amp%evaluate(n,partner_p,partner_hel,.false.,model)
    if (first_amp%n_amps.lt.1 .or. second_amp%n_amps.lt.1) then
       write (*,*) 'Missing multichannel-partner amplitude:',&
            first_amp%n_amps,second_amp%n_amps
       stop 1
    endif
    first_squared=abs(sum(first_amp%amps_by_order(:,&
         leading_qcd_sector(first_amp))))**2
    second_squared=abs(sum(second_amp%amps_by_order(:,&
         leading_qcd_sector(second_amp))))**2
    write (*,'(a,2es24.16)') 'THREE_LINE_SAFE_MC_PARTNER=',&
         first_squared,second_squared
    scale=max(first_squared,second_squared)
    if (scale.lt.1d-30 .or. abs(first_squared-second_squared).gt.1d-12*scale) then
       write (*,*) 'Safe three-line multichannel partners disagree:',&
            first_squared,second_squared
       stop 1
    endif
  end subroutine check_safe_multichannel_partner

  subroutine check_quark_flavour_case(process,madgraph_value,normalization,label)
    implicit none
    integer,dimension(n),intent(in) :: process
    real(kind=dp),intent(in) :: madgraph_value,normalization
    character(len=*),intent(in) :: label
    type(amplitude_QCD) :: identical_amp
    integer,dimension(n,1) :: identical_part,identical_orders
    integer,dimension(0:3,n) :: identical_spin
    integer,dimension(n) :: identical_hel
    real(kind=dp),dimension(0:3,n) :: identical_p
    real(kind=dp),dimension(3) :: raw_identical,normalized_identical

    identical_part(:,1)=process
    identical_orders=0
    identical_hel=[-1,1,-1,1,-1,1]
    identical_spin=0
    identical_spin(0,:)=1
    identical_spin(1,:)=-9
    call fill_momenta(identical_p)
    call identical_amp%init(2,n,1,identical_part,identical_spin,&
         identical_orders,model)
    call identical_amp%init_col(n,20)
    call identical_amp%evaluate(n,identical_p,identical_hel,.false.,model)
    call colour_squared(identical_amp,raw_identical)
    normalized_identical=raw_identical*(4d0*pi*alpha_s)**4/normalization
    write (*,'(a,a,a,3es24.16)') 'THREE_LINE_IDENTICAL[',trim(label),&
         ']_LC_NLC_FC=',normalized_identical
    if (abs(normalized_identical(3)/(madgraph_value/normalization)-1d0).gt.1d-11) then
       write (*,*) trim(label),' result disagrees with MadGraph:',&
            normalized_identical(3),madgraph_value/normalization
       stop 1
    endif
  end subroutine check_quark_flavour_case

  subroutine colour_squared(amplitude,result)
    implicit none
    type(amplitude_QCD),intent(in) :: amplitude
    real(kind=dp),dimension(3),intent(out) :: result
    integer :: iacc,irow,ival,ic,icol,ioff
    complex(kind=dp) :: weighted,sum_for_factor
    complex(kind=dp),dimension(:),allocatable :: leading_amps

    result=0d0
    allocate(leading_amps(amplitude%n_amps))
    leading_amps=amplitude%amps_by_order(:,leading_qcd_sector(amplitude))
    ioff=amplitude%iproc_start(amplitude%nprocs)-1
    do iacc=1,3
       do irow=1,amplitude%nColOrd
          weighted=(0d0,0d0)
          do ival=1,amplitude%n_col_vals(iacc)
             sum_for_factor=(0d0,0d0)
             do ic=amplitude%row_index(irow-1,ival,iacc)+1,&
                  amplitude%row_index(irow,ival,iacc)
                icol=amplitude%col_index(amplitude%i_col_i(ival,iacc)+ic)
                sum_for_factor=sum_for_factor+leading_amps(ioff+icol)
             enddo
             weighted=weighted+sum_for_factor*amplitude%diff_col_vals(ival,iacc)
          enddo
          result(iacc)=result(iacc)+dble(weighted*&
               conjg(leading_amps(ioff+irow)))
       enddo
    enddo
    deallocate(leading_amps)
  end subroutine colour_squared

  integer function leading_qcd_sector(amplitude)
    implicit none
    type(amplitude_QCD),intent(in) :: amplitude

    if (.not.allocated(amplitude%amps_by_order) .or. amplitude%n_sectors.le.0) then
       write (*,*) 'Coupling-sector amplitudes are unavailable'
       stop 1
    endif
    leading_qcd_sector=maxloc(amplitude%sector_powers(1,1:amplitude%n_sectors),dim=1)
  end function leading_qcd_sector

  subroutine fill_momenta(momenta)
    implicit none
    real(kind=dp),dimension(0:3,n),intent(out) :: momenta
    momenta(:,1)=[4671.200996478833d0,0d0,0d0,4671.200996478833d0]
    momenta(:,2)=[5452.624496459750d0,0d0,0d0,-5452.624496459750d0]
    momenta(:,3)=[3848.769685279069d0,-1105.1428951212934d0,&
         -1737.2006769281086d0,-3251.7412381317486d0]
    momenta(:,4)=[3660.1488063565753d0,2228.296688374595d0,&
         1630.5712325370819d0,2402.6278548445193d0]
    momenta(:,5)=[1275.9167404758375d0,-470.2433927087442d0,&
         -436.62760397349416d0,1102.8105076070963d0]
    momenta(:,6)=[1338.9902608271016d0,-652.9104005445573d0,&
         543.2570483645209d0,-1035.1206243007837d0]
  end subroutine fill_momenta

end program three_quark_line_reweight_regression
