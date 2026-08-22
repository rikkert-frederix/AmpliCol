program heft_regression
  use amplitude_QCD_mod
  use particles
  implicit none
  integer,parameter :: dp=kind(1d0)
  real(kind=dp),parameter :: gs=1.2177157847767197d0
  real(kind=dp),parameter :: gew=sqrt(8d0*acos(-1d0)/132.507d0)
  real(kind=dp),parameter :: gheft=5.248659993424408d-5
  type(physics_model) :: model
  character(len=64) :: option

  open(unit=99,file='/dev/null',status='unknown',action='write')
  call model%init_part(173d0,1.4915d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0d0)
  call model%set_heft_enabled(.true.)
  call model%init_vert()

  if (command_argument_count().gt.0) then
     call get_command_argument(1,option)
     if (trim(option).eq.'--emit-library') then
        call emit_heft_libraries()
        stop
     endif
     write (*,*) 'Unknown option: ',trim(option)
     stop 1
  endif

  call check_hgg_helicity_rule()
  call check_gg_to_gh_madgraph()
  call check_gg_to_ggh_madgraph()
  call check_gg_to_h3g_madgraph()
  call check_uubar_to_h4g_madgraph()
  call check_ttbar_to_gh_madgraph()
  write (*,'(a)') 'HEFT direct and MadGraph regression passed'

contains

  subroutine check_hgg_helicity_rule()
    implicit none
    integer,parameter :: n=3
    type(amplitude_QCD) :: amplitude
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    real(kind=dp),dimension(0:3,n) :: p

    part(:,1)=[21,21,25]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    p(:,1)=[62.5d0,0d0,0d0,62.5d0]
    p(:,2)=[62.5d0,0d0,0d0,-62.5d0]
    p(:,3)=[125d0,0d0,0d0,0d0]
    call amplitude%init(2,n,1,part,spin,orders,model)
    call amplitude%init_col(n,20)

    hel=[-1,-1,0]
    call amplitude%evaluate(n,p,hel,.false.,model,1d0,1d0,1d0)
    call assert_complex_close('Hgg equal-helicity amplitude',amplitude%amps(1),&
         cmplx(0d0,-7812.5d0,kind=dp),1d-13)

    hel=[-1,1,0]
    call amplitude%evaluate(n,p,hel,.false.,model,1d0,1d0,1d0)
    call assert_complex_close('Hgg mixed-helicity zero',amplitude%amps(1),&
         cmplx(0d0,0d0,kind=dp),1d-13)
  end subroutine check_hgg_helicity_rule

  subroutine check_gg_to_gh_madgraph()
    ! Fixed-helicity, unaveraged, full-colour |M|^2 values from
    ! MadGraph5_aMC@NLO 3.7.2 with the bundled sm__hgg_plugin model.
    implicit none
    integer,parameter :: n=4,nchecks=4
    type(amplitude_QCD) :: amplitude
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n,nchecks) :: helicities
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(nchecks),parameter :: expected=[&
         5.57734416993195814d-2,8.51035182179562960d-7,&
         2.69272850611501071d-3,2.69272850611500941d-3]
    complex(kind=dp),dimension(2) :: nominal,doubled
    integer :: i

    part(:,1)=[21,21,21,25]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    helicities(:,1)=[-1,-1,-1,0]
    helicities(:,2)=[-1,-1,1,0]
    helicities(:,3)=[-1,1,-1,0]
    helicities(:,4)=[-1,1,1,0]
    p(:,1)=[250d0,0d0,0d0,250d0]
    p(:,2)=[250d0,0d0,0d0,-250d0]
    p(:,3)=[234.375d0,234.375d0,0d0,0d0]
    p(:,4)=[265.625d0,-234.375d0,0d0,0d0]
    call amplitude%init(2,n,1,part,spin,orders,model)
    call amplitude%init_col(n,20)
    if (amplitude%n_amps.ne.2) then
       write (*,*) 'gg > gh should have two colour-ordered coefficients'
       stop 1
    endif

    do i=1,nchecks
       call amplitude%evaluate(n,p,helicities(:,i),.false.,model,gs,gew,gheft)
       call assert_real_close('MadGraph gg > gh',full_colour_squared(amplitude),&
            expected(i),2d-12)
       if (i.eq.1) nominal=amplitude%amps
    enddo

    call amplitude%evaluate(n,p,helicities(:,1),.false.,model,gs,gew,2d0*gheft)
    doubled=amplitude%amps
    do i=1,2
       call assert_complex_close('one HEFT insertion in gg > gh',doubled(i),&
            2d0*nominal(i),2d-12)
    enddo
  end subroutine check_gg_to_gh_madgraph

  subroutine check_gg_to_ggh_madgraph()
    ! Same MadGraph setup as check_gg_to_gh_madgraph at a symmetric 2->3 point.
    implicit none
    integer,parameter :: n=5,nchecks=4
    type(amplitude_QCD) :: amplitude
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n,nchecks) :: helicities
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(nchecks),parameter :: expected=[&
         4.83602074181558366d-5,2.14613251468425786d-6,&
         4.36127146077883778d-7,1.24044308952597881d-6]
    real(kind=dp) :: energy
    integer :: i

    part(:,1)=[21,21,21,21,25]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    helicities(:,1)=[-1,-1,-1,-1,0]
    helicities(:,2)=[-1,-1,-1,1,0]
    helicities(:,3)=[-1,1,-1,1,0]
    helicities(:,4)=[-1,1,1,1,0]
    energy=(2000d0-sqrt(1187500d0))/6d0
    p(:,1)=[250d0,0d0,0d0,250d0]
    p(:,2)=[250d0,0d0,0d0,-250d0]
    p(:,3)=[energy,0.5d0*energy,0.5d0*sqrt(3d0)*energy,0d0]
    p(:,4)=[energy,0.5d0*energy,-0.5d0*sqrt(3d0)*energy,0d0]
    p(:,5)=[500d0-2d0*energy,-energy,0d0,0d0]
    call amplitude%init(2,n,1,part,spin,orders,model)
    call amplitude%init_col(n,20)
    if (amplitude%n_amps.ne.6) then
       write (*,*) 'gg > ggh should have six colour-ordered coefficients'
       stop 1
    endif

    do i=1,nchecks
       call amplitude%evaluate(n,p,helicities(:,i),.false.,model,gs,gew,gheft)
       call assert_real_close('MadGraph gg > ggh',full_colour_squared(amplitude),&
            expected(i),2d-12)
    enddo
  end subroutine check_gg_to_ggh_madgraph

  subroutine check_gg_to_h3g_madgraph()
    ! One generated-event point from the cross-section comparison.  The
    ! MadGraph hgg_plugin top-mass expansion was put in the strict heavy-top
    ! limit for this value.
    implicit none
    integer,parameter :: n=6
    real(kind=dp),parameter :: gs_event=1.13730550681465803d0
    real(kind=dp),parameter :: gheft_event=4.43558169905496733d-5
    real(kind=dp),parameter :: expected=1.58479536096905880d-6
    real(kind=dp),parameter :: expected_sum=5.16448908303389631d-6
    type(amplitude_QCD) :: amplitude
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp) :: helicity_sum
    integer :: h1,h2,h4,h5,h6

    part(:,1)=[21,21,25,21,21,21]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    hel=[-1,-1,0,-1,-1,-1]
    p(:,1)=[373.75312326566302d0,0d0,0d0,373.75312326566302d0]
    p(:,2)=[403.47072902643276d0,0d0,0d0,-403.47072902643276d0]
    p(:,3)=[216.54829708118396d0,-42.341877309465147d0,&
         147.54305366710500d0,87.785976723259580d0]
    p(:,4)=[200.03287849382329d0,-34.923666820043678d0,&
         -20.350832369065799d0,-195.90644092590080d0]
    p(:,5)=[220.58980288703404d0,76.031960247160669d0,&
         -203.34116228757387d0,-39.132772432492416d0]
    p(:,6)=[140.05287383005444d0,1.2335838823481566d0,&
         76.148940989534665d0,117.53563087436389d0]
    call amplitude%init(2,n,1,part,spin,orders,model)
    call amplitude%init_col(n,20)
    if (amplitude%n_amps.ne.24) then
       write (*,*) 'gg > hggg should have 24 colour coefficients'
       stop 1
    endif
    call amplitude%evaluate(n,p,hel,.false.,model,gs_event,gew,gheft_event)
    call assert_real_close('MadGraph gg > hggg event point',&
         full_colour_squared(amplitude),expected,3d-11)

    helicity_sum=0d0
    do h1=-1,1,2
       do h2=-1,1,2
          do h4=-1,1,2
             do h5=-1,1,2
                do h6=-1,1,2
                   hel=[h1,h2,0,h4,h5,h6]
                   call amplitude%evaluate(n,p,hel,.false.,model,gs_event,gew,&
                        gheft_event)
                   helicity_sum=helicity_sum+full_colour_squared(amplitude)
                enddo
             enddo
          enddo
       enddo
    enddo
    call assert_real_close('MadGraph gg > hggg helicity sum',helicity_sum,&
         expected_sum,3d-11)
  end subroutine check_gg_to_h3g_madgraph

  subroutine check_uubar_to_h4g_madgraph()
    ! Non-planar 2->5 point for the complete H+4g amplitude on a quark line.
    ! MadGraph's fixed-helicity matrix routine includes every colour
    ! interference term but no initial-state spin or colour average.
    implicit none
    integer,parameter :: n=7,nchecks=8
    type(amplitude_QCD) :: amplitude
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n,nchecks) :: helicities
    integer,dimension(n) :: hel
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(nchecks),parameter :: expected=[&
         1.08933201738183408d-14,1.31856701884300264d-13,&
         2.66148301159155096d-14,1.42066241393605224d-13,&
         1.55692668111654504d-13,2.10973848619359794d-15,&
         1.98400952820559632d-13,1.87816924306901880d-13]
    real(kind=dp) :: ecm,ein,value
    integer :: i

    part(:,1)=[2,-2,25,21,21,21,21]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    helicities(:,1)=[-1,1,0,-1,-1,-1,-1]
    helicities(:,2)=[-1,1,0,-1,-1,-1,1]
    helicities(:,3)=[-1,1,0,-1,-1,1,-1]
    helicities(:,4)=[-1,1,0,-1,1,-1,-1]
    helicities(:,5)=[-1,1,0,1,-1,-1,1]
    helicities(:,6)=[-1,1,0,1,-1,1,-1]
    helicities(:,7)=[-1,1,0,1,1,-1,-1]
    helicities(:,8)=[1,-1,0,-1,-1,1,-1]
    ecm=sqrt(8900d0)+sqrt(13400d0)+sqrt(16500d0)+&
         sqrt(11475d0)+sqrt(19500d0)
    ein=0.5d0*ecm
    p(:,1)=[ein,0d0,0d0,ein]
    p(:,2)=[ein,0d0,0d0,-ein]
    p(:,3)=[sqrt(19500d0),-25d0,-15d0,-55d0]
    p(:,4)=[sqrt(8900d0),80d0,30d0,40d0]
    p(:,5)=[sqrt(13400d0),-20d0,110d0,-30d0]
    p(:,6)=[sqrt(16500d0),-70d0,-40d0,100d0]
    p(:,7)=[sqrt(11475d0),35d0,-85d0,-55d0]
    call amplitude%init(2,n,1,part,spin,orders,model)
    call amplitude%init_col(n,20)
    if (amplitude%n_amps.ne.24) then
       write (*,*) 'u ubar > hgggg should have 24 colour coefficients'
       stop 1
    endif

    do i=1,nchecks
       call amplitude%evaluate(n,p,helicities(:,i),.false.,model,gs,gew,gheft)
       call assert_real_close('MadGraph u ubar > hgggg',&
            full_colour_squared(amplitude),expected(i),3d-11)
    enddo

    ! A massless vector current only connects opposite incoming helicities.
    hel=[-1,-1,0,-1,1,-1,1]
    call amplitude%evaluate(n,p,hel,.false.,model,gs,gew,gheft)
    value=full_colour_squared(amplitude)
    if (abs(value).gt.1d-25) then
       write (*,*) 'u ubar > hgggg equal-helicity value should vanish:',value
       stop 1
    endif
  end subroutine check_uubar_to_h4g_madgraph

  subroutine check_ttbar_to_gh_madgraph()
    ! This point contains both the SM top-Yukawa graphs and the HEFT graphs.
    ! Agreement therefore fixes their relative phase as well as each sector.
    implicit none
    integer,parameter :: n=4,nchecks=8
    type(amplitude_QCD) :: amplitude
    integer,dimension(n,1) :: part,orders
    integer,dimension(0:3,n) :: spin
    integer,dimension(n) :: hel
    real(kind=dp),dimension(0:3,n) :: p
    real(kind=dp),dimension(nchecks),parameter :: expected=[&
         8.69725749136089910d0,1.01017541752548023d0,&
         1.87579268266963517d0,1.87579268266963517d0,&
         1.87579268266963517d0,1.87579268266963517d0,&
         1.01017541752548023d0,8.69725749136089910d0]
    complex(kind=dp),dimension(1) :: sm,total,doubled
    real(kind=dp) :: pz,pg
    integer :: h1,h2,h3,icheck,i

    part(:,1)=[6,-6,21,25]
    orders=0
    spin=0
    spin(0,:)=1
    spin(1,:)=-9
    pz=sqrt(300d0**2-173d0**2)
    pg=(600d0**2-125d0**2)/(2d0*600d0)
    p(:,1)=[300d0,0d0,0d0,pz]
    p(:,2)=[300d0,0d0,0d0,-pz]
    p(:,3)=[pg,pg,0d0,0d0]
    p(:,4)=[600d0-pg,-pg,0d0,0d0]
    call amplitude%init(2,n,1,part,spin,orders,model)
    call amplitude%init_col(n,20)
    if (amplitude%n_amps.ne.1) then
       write (*,*) 'ttbar > gh should have one colour-ordered coefficient'
       stop 1
    endif

    icheck=0
    do h1=-1,1,2
       do h2=-1,1,2
          do h3=-1,1,2
             icheck=icheck+1
             hel=[h1,h2,h3,0]
             call amplitude%evaluate(n,p,hel,.false.,model,gs,gew,gheft)
             call assert_real_close('MadGraph ttbar > gh interference',&
                  full_colour_squared(amplitude),expected(icheck),2d-12)
          enddo
       enddo
    enddo

    hel=[-1,-1,-1,0]
    call amplitude%evaluate(n,p,hel,.false.,model,gs,gew,0d0)
    sm=amplitude%amps
    call amplitude%evaluate(n,p,hel,.false.,model,gs,gew,gheft)
    total=amplitude%amps
    call amplitude%evaluate(n,p,hel,.false.,model,gs,gew,2d0*gheft)
    doubled=amplitude%amps
    if (abs(total(1)-sm(1)).lt.1d-12) then
       write (*,*) 'ttbar > gh has no HEFT contribution to interfere with the SM'
       stop 1
    endif
    do i=1,size(total)
       call assert_complex_close('SM plus one-insertion HEFT sector',doubled(i),&
            2d0*total(i)-sm(i),2d-12)
    enddo
  end subroutine check_ttbar_to_gh_madgraph

  subroutine emit_heft_libraries()
    implicit none
    integer,parameter :: n3=3,n4=4,n5=5
    type(amplitude_QCD) :: amp5,amp4,amp3
    integer,dimension(n5,1) :: part5,order5
    integer,dimension(0:3,n5) :: spin5
    integer,dimension(n5) :: hel5
    real(kind=dp),dimension(0:3,n5) :: p5
    integer,dimension(n4,1) :: part4,order4
    integer,dimension(0:3,n4) :: spin4
    integer,dimension(n4) :: hel4
    real(kind=dp),dimension(0:3,n4) :: p4
    integer,dimension(n3,1) :: part3,order3
    integer,dimension(0:3,n3) :: spin3
    integer,dimension(n3) :: hel3
    real(kind=dp),dimension(0:3,n3) :: p3
    real(kind=dp) :: energy,pz,pg

    part5(:,1)=[21,21,21,21,25]
    order5(:,1)=[5,1,3,4,2]
    spin5=0
    spin5(0,:)=1
    hel5=[-1,-1,-1,-1,0]
    spin5(1,:)=hel5
    energy=(2000d0-sqrt(1187500d0))/6d0
    p5(:,1)=[250d0,0d0,0d0,250d0]
    p5(:,2)=[250d0,0d0,0d0,-250d0]
    p5(:,3)=[energy,0.5d0*energy,0.5d0*sqrt(3d0)*energy,0d0]
    p5(:,4)=[energy,0.5d0*energy,-0.5d0*sqrt(3d0)*energy,0d0]
    p5(:,5)=[500d0-2d0*energy,-energy,0d0,0d0]
    call amp5%init(1,n5,1,part5,spin5,order5,model)
    call amp5%evaluate(n5,p5,hel5,.false.,model,gs,gew,gheft)
    call amp5%create_library(n5,hel5,1,1,model,p5,gs,gew,gheft)

    part4(:,1)=[6,-6,21,25]
    order4(:,1)=[2,3,4,1]
    spin4=0
    spin4(0,:)=1
    hel4=[-1,-1,-1,0]
    spin4(1,:)=hel4
    pz=sqrt(300d0**2-173d0**2)
    pg=(600d0**2-125d0**2)/(2d0*600d0)
    p4(:,1)=[300d0,0d0,0d0,pz]
    p4(:,2)=[300d0,0d0,0d0,-pz]
    p4(:,3)=[pg,pg,0d0,0d0]
    p4(:,4)=[600d0-pg,-pg,0d0,0d0]
    call amp4%init(1,n4,1,part4,spin4,order4,model)
    call amp4%evaluate(n4,p4,hel4,.false.,model,gs,gew,gheft)
    call amp4%create_library(n4,hel4,2,1,model,p4,gs,gew,gheft)

    ! Keep the Higgs as the closing external current.  This specifically
    ! exercises the one-component scalar contraction in emitted code.
    part3(:,1)=[21,21,25]
    order3(:,1)=[1,2,3]
    spin3=0
    spin3(0,:)=1
    hel3=[-1,-1,0]
    spin3(1,:)=hel3
    p3(:,1)=[62.5d0,0d0,0d0,62.5d0]
    p3(:,2)=[62.5d0,0d0,0d0,-62.5d0]
    p3(:,3)=[125d0,0d0,0d0,0d0]
    call amp3%init(1,n3,1,part3,spin3,order3,model)
    call amp3%evaluate(n3,p3,hel3,.false.,model,gs,gew,gheft)
    call amp3%create_library(n3,hel3,3,1,model,p3,gs,gew,gheft)
  end subroutine emit_heft_libraries

  real(kind=dp) function full_colour_squared(amplitude)
    implicit none
    type(amplitude_QCD),intent(in) :: amplitude
    integer :: irow,ival,ic,icol
    complex(kind=dp) :: weighted,sum_for_factor

    full_colour_squared=0d0
    do irow=1,amplitude%nColOrd
       weighted=(0d0,0d0)
       do ival=1,amplitude%n_col_vals(3)
          sum_for_factor=(0d0,0d0)
          do ic=amplitude%row_index(irow-1,ival,3)+1,&
               amplitude%row_index(irow,ival,3)
             icol=amplitude%col_index(amplitude%i_col_i(ival,3)+ic)
             sum_for_factor=sum_for_factor+amplitude%amps(icol)
          enddo
          weighted=weighted+sum_for_factor*amplitude%diff_col_vals(ival,3)
       enddo
       full_colour_squared=full_colour_squared+&
            dble(weighted*conjg(amplitude%amps(irow)))
    enddo
  end function full_colour_squared

  subroutine assert_real_close(label,value,reference,tolerance)
    implicit none
    character(len=*),intent(in) :: label
    real(kind=dp),intent(in) :: value,reference,tolerance
    real(kind=dp) :: relative_difference

    relative_difference=abs(value-reference)/max(1d-30,abs(value)+abs(reference))
    if (relative_difference.gt.tolerance) then
       write (*,*) trim(label),' mismatch:',value,reference
       write (*,*) 'relative difference:',relative_difference
       stop 1
    endif
  end subroutine assert_real_close

  subroutine assert_complex_close(label,value,reference,tolerance)
    implicit none
    character(len=*),intent(in) :: label
    complex(kind=dp),intent(in) :: value,reference
    real(kind=dp),intent(in) :: tolerance
    real(kind=dp) :: difference,scale

    difference=abs(value-reference)
    scale=max(1d-30,abs(value)+abs(reference))
    if (difference.gt.tolerance*scale .and. difference.gt.tolerance) then
       write (*,*) trim(label),' mismatch:',value,reference
       write (*,*) 'absolute/relative difference:',difference,difference/scale
       stop 1
    endif
  end subroutine assert_complex_close

end program heft_regression
