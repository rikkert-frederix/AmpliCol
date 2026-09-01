!***
!*
!* NNPDF Fortran Driver
!*
!* Stefano Carrazza for the NNPDF Collaboration
!* email: stefano.carrazza@mi.infn.it
!*
!* February 2013
!*
!* Usage:
!*
!*  NNPDFDriver("gridname.LHgrid");
!*
!*  NNinitPDF(0); // select replica [0,Mem]
!*
!*  NNevolvePDF(x,Q,pdf); // -> returns double array (-6,7)
!*     
subroutine NNPDFDriver(gridfilename)
  use pdf_internal_interface, only: readPDFSet
  use pdf_internal_state, only: nfl,nx,nq2,mem,rep,alphas,hasphoton,pdf_grid_loaded
  implicit none

  character(len=*),intent(in) :: gridfilename
  
  nfl = 13
  nx = 100
  nq2 = 60
  mem = 1
  rep = 0
  alphas = 0d0
  hasphoton=.false.
  pdf_grid_loaded=.false.
  write(99,*) " ****************************************"
  write(99,*) ""
  write(99,*) "      NNPDFDriver version 1.0.3"
  write(99,*) "  Grid: ", gridfilename
  write(99,*) " ****************************************"

  call readPDFSet(gridfilename)
      
end subroutine NNPDFDriver

subroutine NNinitPDF(irep)
  use pdf_internal_state, only: mem,rep,pdf_grid_loaded
  implicit none
  integer,intent(in) :: irep

  if (.not.pdf_grid_loaded) then
     write(*,*) 'Error: cannot select a PDF replica before loading a grid'
     stop 1
  endif
  
  if (irep.gt.mem.or.irep.lt.0) then
     write(*,*) 'Error: PDF replica out of range [0,',mem,']:',irep
     stop 1
  endif
  rep = irep
  
end subroutine NNinitPDF

subroutine readPDFSet(gridfilename)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use pdf_internal_state
  implicit none

  integer :: i,ix,iq,fl,imem,iu,ios,close_status
  character(len=*),intent(in) :: gridfilename
  character(len=100) :: line
  character(len=1024) :: grid_path
  character(len=512) :: io_message
  logical :: found,grid_open

  pdf_grid_loaded=.false.
  grid_open=.false.
  if (len_trim(gridfilename).eq.0) then
     call grid_error('empty PDF grid filename')
  endif
  io_message=''
  open(newunit=iu,file=trim(gridfilename),status='old',action='read',&
       iostat=ios,iomsg=io_message)
  if (ios.ne.0) then
     if (len_trim(gridfilename).gt.len(grid_path)-4) then
        call grid_error('PDF grid filename is too long')
     endif
     grid_path='PDF/'//trim(gridfilename)
     io_message=''
     open(newunit=iu,file=trim(grid_path),status='old',action='read',&
          iostat=ios,iomsg=io_message)
     if (ios.ne.0) call grid_io_error('could not open PDF grid',io_message)
  endif
  grid_open=.true.

  found=.false.
  do i=1,1000
     io_message=''
     read(iu,*,iostat=ios,iomsg=io_message) line
     if (ios.ne.0) call grid_io_error('reading the PDF parameter header',io_message)
     if (line(1:14).eq.'Parameterlist:') then
        io_message=''
        read(iu,*,iostat=ios,iomsg=io_message) line,mem,line,alphas
        if (ios.ne.0) call grid_io_error('reading PDF member metadata',io_message)
        found=.true.
        exit
     endif
  enddo
  if (.not.found) call grid_error('missing Parameterlist header in PDF grid')
  if (mem.lt.0 .or. mem.gt.max_pdf_members) then
     call grid_error('PDF member count exceeds the supported range')
  endif
  if (.not.ieee_is_finite(alphas)) then
     call grid_error('invalid alpha_s value in PDF grid')
  endif
  if (alphas.le.0d0) then
     call grid_error('invalid alpha_s value in PDF grid')
  endif

  found=.false.
  do i=1,1000
     io_message=''
     read(iu,*,iostat=ios,iomsg=io_message) line
     if (ios.ne.0) call grid_io_error('reading the PDF interpolation header',io_message)
     if (line(1:13).eq.'NNPDF20intqed') then
        hasphoton = .true.
        nfl = 14
        io_message=''
        read(iu,*,iostat=ios,iomsg=io_message) line,line
        if (ios.ne.0) call grid_io_error('reading the QED PDF header',io_message)
        found=.true.
        exit
     endif
     if (line(1:13).eq.'NNPDF20int') then
        hasphoton = .false.
        nfl = 13
        io_message=''
        read(iu,*,iostat=ios,iomsg=io_message) line,line
        if (ios.ne.0) call grid_io_error('reading the PDF header',io_message)
        found=.true.
        exit
     endif
  enddo
  if (.not.found) call grid_error('missing NNPDF interpolation header in PDF grid')

  io_message=''
  read(iu,*,iostat=ios,iomsg=io_message) nx
  if (ios.ne.0) call grid_io_error('reading the PDF x-grid size',io_message)
  if (nx.lt.4 .or. nx.gt.max_pdf_x_points) then
     call grid_error('PDF x-grid size exceeds the supported range')
  endif
  io_message=''
  read(iu,*,iostat=ios,iomsg=io_message) xgrid(1)
  if (ios.ne.0) call grid_io_error('reading the PDF x grid',io_message)
  if (.not.ieee_is_finite(xgrid(1))) then
     call grid_error('PDF x grid contains a value outside (0,1]')
  endif
  if (xgrid(1).le.0d0 .or. xgrid(1).gt.1d0) then
     call grid_error('PDF x grid contains a value outside (0,1]')
  endif
  logxgrid(1)=dlog(xgrid(1))
  do ix=2,nx
     io_message=''
     read(iu,*,iostat=ios,iomsg=io_message) xgrid(ix)
     if (ios.ne.0) call grid_io_error('reading the PDF x grid',io_message)
     if (.not.ieee_is_finite(xgrid(ix))) then
        call grid_error('PDF x grid contains a value outside (0,1]')
     endif
     if (xgrid(ix).le.0d0 .or. xgrid(ix).gt.1d0) then
        call grid_error('PDF x grid contains a value outside (0,1]')
     endif
     if (xgrid(ix).le.xgrid(ix-1)) call grid_error('PDF x grid is not strictly increasing')
     logxgrid(ix) = dlog(xgrid(ix))
  enddo
  io_message=''
  read(iu,*,iostat=ios,iomsg=io_message) nq2
  if (ios.ne.0) call grid_io_error('reading the PDF Q2-grid size',io_message)
  if (nq2.lt.2 .or. nq2.gt.max_pdf_q2_points) then
     call grid_error('PDF Q2-grid size exceeds the supported range')
  endif
  io_message=''
  read(iu,*,iostat=ios,iomsg=io_message) line
  if (ios.ne.0) call grid_io_error('reading the PDF Q2-grid header',io_message)
  io_message=''
  read(iu,*,iostat=ios,iomsg=io_message) q2grid(1)
  if (ios.ne.0) call grid_io_error('reading the PDF Q2 grid',io_message)
  if (.not.ieee_is_finite(q2grid(1))) then
     call grid_error('PDF Q2 grid contains a nonpositive or non-finite value')
  endif
  if (q2grid(1).le.0d0) then
     call grid_error('PDF Q2 grid contains a nonpositive or non-finite value')
  endif
  logq2grid(1)=dlog(q2grid(1))
  do iq=2,nq2
     io_message=''
     read(iu,*,iostat=ios,iomsg=io_message) q2grid(iq)
     if (ios.ne.0) call grid_io_error('reading the PDF Q2 grid',io_message)
     if (.not.ieee_is_finite(q2grid(iq))) then
        call grid_error('PDF Q2 grid contains a nonpositive or non-finite value')
     endif
     if (q2grid(iq).le.0d0) then
        call grid_error('PDF Q2 grid contains a nonpositive or non-finite value')
     endif
     if (q2grid(iq).le.q2grid(iq-1)) call grid_error('PDF Q2 grid is not strictly increasing')
     logq2grid(iq) = dlog(q2grid(iq))
  enddo
  io_message=''
  read(iu,*,iostat=ios,iomsg=io_message) line
  if (ios.ne.0) call grid_io_error('reading the PDF data header',io_message)
  do imem=0,mem
     do ix=1,nx
        do iq=1,nq2
           io_message=''
           read(iu,*,iostat=ios,iomsg=io_message) (pdfgrid(imem,fl,ix,iq),fl=1,nfl)
           if (ios.ne.0) call grid_io_error('reading PDF interpolation data',io_message)
           if (.not.all(ieee_is_finite(pdfgrid(imem,1:nfl,ix,iq)))) then
              call grid_error('PDF interpolation data contain a non-finite value')
           endif
        enddo
     enddo
  enddo

  io_message=''
  close(iu,iostat=close_status,iomsg=io_message)
  grid_open=.false.
  if (close_status.ne.0) call grid_io_error('closing the PDF grid',io_message)
  rep=0
  pdf_grid_loaded=.true.

contains
  subroutine grid_io_error(context,detail)
    character(len=*),intent(in) :: context,detail
    call grid_error(trim(context)//': '//trim(detail))
  end subroutine grid_io_error

  subroutine grid_error(message)
    character(len=*),intent(in) :: message
    logical :: is_open
    integer :: cleanup_status
    character(len=256) :: cleanup_message
    write(*,*) 'Error: invalid internal PDF grid: ',trim(message)
    if (grid_open) then
       inquire(unit=iu,opened=is_open)
       if (is_open) close(iu,iostat=cleanup_status,iomsg=cleanup_message)
       grid_open=.false.
    endif
    stop 1
  end subroutine grid_error
end subroutine readPDFSet

subroutine NNevolvePDF(ipdgs,x,Q,xpdf)
  use pdf_internal_interface, only: lh_polin2
  use pdf_internal_state
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
      
  logical,intent(in) :: ipdgs(-6:7)
  integer i,j,ix,iq2,ipdf,fmax
  integer minx,maxx,midx
  integer minq,maxq,midq
  double precision,intent(in) :: x,Q
  double precision,intent(out) :: xpdf(-6:7)
  double precision Q2,xx
  double precision x2,x1,dy,y
  integer,parameter :: m=4,n=2
  double precision,parameter :: xmingrid=1d-7, xch=1d-1

  integer ix1a(m), ix2a(n)
  double precision x1a(m), x2a(n)
  double precision ya(m,n)
  
  xpdf=0d0
  if (.not.pdf_grid_loaded) then
     write(*,*) 'Error: cannot evaluate an internal PDF before loading its grid'
     stop 1
  endif
  if (rep.lt.0 .or. rep.gt.mem) then
     write(*,*) 'Error: invalid active internal PDF replica:',rep
     stop 1
  endif
  if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(Q)) then
     write(*,*) 'Error: invalid x or scale supplied to the internal PDF:',x,Q
     stop 1
  endif
  if (x.le.0d0 .or. Q.le.0d0) then
     write(*,*) 'Error: invalid x or scale supplied to the internal PDF:',x,Q
     stop 1
  endif
  if (Q.gt.sqrt(huge(1d0))) then
     write(*,*) 'Error: internal PDF scale is too large to square:',Q
     stop 1
  endif
  Q2 = Q*Q
  xx=x
  if (xx.lt.xmingrid.or.xx.lt.xgrid(1).or.xx.gt.xgrid(nx)) then
!         write(6,*) "Parton interpolation: x out of range -- freezed"
     if (xx.lt.xgrid(1)) xx = xgrid(1)
     if (xx.gt.xgrid(nx))xx = xgrid(nx)
  endif
  if (Q2.lt.q2grid(1).or.Q2.gt.q2grid(nq2)) then
!         write(6,*) "Parton interpolation: Q2 out of range -- freezed"
!         write(6,*) "Q2 = ",Q2, " GeV2", q2grid(1)
     if (Q2.lt.q2grid(1)) Q2 = q2grid(1)
     if (Q2.gt.q2grid(nq2)) Q2 = q2grid(nq2)
  endif
  minx = 1
  maxx = NX+1
10 continue
  midx = (minx+maxx)/2
  if (xx.lt.xgrid(midx)) then
     maxx=midx
  else
     minx=midx
  endif
  if ((maxx-minx).gt.1) go to 10
  ix = minx

  minq = 1
  maxq = nq2+1
20 continue
  midq = (minq+maxq)/2
  if (Q2.lt.q2grid(midq)) then
     maxq=midq
  else
     minq=midq
  endif
  if ((maxq-minq).gt.1) go to 20
  iq2 = minq

  do I=1,M
     if(IX.ge.M/2.and.IX.le.(NX-M/2)) IX1A(I) = IX - M/2 + I
     if(IX.lt.M/2) IX1A(I) = I
     if(IX.gt.(NX-M/2)) IX1A(I) = (NX - M) + I
         
     if(IX1A(I).le.0.or.IX1A(I).gt.NX) then
        write(6,*) "Error in grids! "
        write(6,*) "I, IXIA(I) = ",I, IX1A(I)
        stop 10
     endif
  enddo
  
  do J=1,N
     if(IQ2.ge.N/2.and.IQ2.le.(NQ2-N/2)) IX2A(J) = IQ2 - N/2 + J
     if(IQ2.lt.N/2) IX2A(J) = J
     if(IQ2.gt.(NQ2-N/2)) IX2A(J) = (NQ2 - N) + J
     if(IX2A(J).le.0.or.IX2A(J).gt.NQ2) then
        write(6,*) "Error in grids! "
        write(6,*) "J, IXIA(J) = ",J,IX2A(J)
        stop 10
     endif
  enddo
            
  IF(xx.LT.XCH)THEN
     X1=dlog(xx)          
  ELSE
     X1=xx
  ENDIF
  X2=dlog(Q2)
  
  fmax = 6
  if (nfl.eq.14) fmax=7

  DO IPDF = -6,fmax,1
     if (.not.ipdgs(ipdf)) cycle
     DO I=1,M
        IF(xx.LT.XCH)THEN
           X1A(I)= logxgrid(IX1A(I))
        ELSE
           X1A(I)= xgrid(IX1A(I))
        ENDIF
        DO J=1,N
           X2A(J) = logq2grid(IX2A(J))
           YA(I,J) = pdfgrid(REP,IPDF+7,IX1A(I),IX2A(J))               
        enddo
     enddo
     
     !     2D polynomial interpolation
     call lh_polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
     if (.not.ieee_is_finite(y)) then
        write(*,*) 'Error: internal PDF interpolation returned a non-finite value'
        stop 1
     endif
     XPDF(IPDF) = y
  enddo
  
end subroutine NNevolvePDF

subroutine lh_polin2(x1a,x2a,ya,m,n,x1,x2,y,dy) 
  use pdf_internal_interface, only: lh_polint
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none 
!
  integer,intent(in) :: m,n
  integer,parameter :: max_interpolation_points=1024
  integer j,k,allocation_status
  double precision,intent(in):: x1,x2,x1a(m),x2a(n),ya(m,n)
  double precision,intent(out):: dy,y
  double precision,allocatable :: ymtmp(:),yntmp(:)
  character(len=256) :: allocation_message
  y=0d0
  dy=0d0
  if (m.lt.1 .or. n.lt.1 .or. m.gt.max_interpolation_points .or. &
       n.gt.max_interpolation_points .or. &
       int(m,kind=8)*int(n,kind=8).gt.int(max_interpolation_points,kind=8)**2) then
     write(*,*) 'failure in polin2: invalid interpolation dimensions',m,n
     stop 1
  endif
  if (.not.ieee_is_finite(x1) .or. .not.ieee_is_finite(x2) .or. &
       .not.all(ieee_is_finite(x1a)) .or. .not.all(ieee_is_finite(x2a)) .or. &
       .not.all(ieee_is_finite(ya))) then
     write(*,*) 'failure in polin2: non-finite interpolation input'
     stop 1
  endif
  allocate(ymtmp(m),yntmp(n),stat=allocation_status,errmsg=allocation_message)
  if (allocation_status.ne.0) then
     write(*,*) 'failure in polin2: could not allocate interpolation workspace: ',&
          trim(allocation_message)
     stop 1
  endif
  do j=1,m 
     do k=1,n 
        yntmp(k)=ya(j,k) 
     enddo
     call lh_polint(x2a,yntmp,n,x2,ymtmp(j),dy)
  enddo
  call lh_polint(x1a,ymtmp,m,x1,y,dy) 
END subroutine lh_polin2

subroutine lh_polint(xa,ya,n,x,y,dy) 
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none 
  integer,intent(in) :: n
  double precision,intent(in) :: x,xa(n),ya(n)
  double precision,intent(out) :: dy,y
  integer,parameter :: max_interpolation_points=1024
  integer i,m,ns,allocation_status
  double precision den,dif,dift,ho,hp,w,difference,quotient,product
  double precision,allocatable :: c(:),d(:)
  character(len=256) :: allocation_message
  y=0d0
  dy=0d0
  if (n.lt.1 .or. n.gt.max_interpolation_points) then
     write(*,*) 'failure in polint: invalid interpolation dimension',n
     stop 1
  endif
  if (.not.ieee_is_finite(x) .or. .not.all(ieee_is_finite(xa)) .or. &
       .not.all(ieee_is_finite(ya))) then
     write(*,*) 'failure in polint: non-finite interpolation input'
     stop 1
  endif
  allocate(c(n),d(n),stat=allocation_status,errmsg=allocation_message)
  if (allocation_status.ne.0) then
     write(*,*) 'failure in polint: could not allocate interpolation workspace: ',&
          trim(allocation_message)
     stop 1
  endif
  ns=1 
  if (.not.safe_subtract(x,xa(1),difference)) call unsafe_arithmetic()
  dif=abs(difference)
  do i=1,n 
     if (.not.safe_subtract(x,xa(i),difference)) call unsafe_arithmetic()
     dift=abs(difference)
     if(dift.lt.dif) then 
        ns=i 
        dif=dift 
     endif
     c(i)=ya(i) 
     d(i)=ya(i) 
  enddo
  y=ya(ns) 
  ns=ns-1 
  do m=1,n-1 
     do i=1,n-m 
        if (.not.safe_subtract(xa(i),x,ho)) call unsafe_arithmetic()
        if (.not.safe_subtract(xa(i+m),x,hp)) call unsafe_arithmetic()
        if (.not.safe_subtract(c(i+1),d(i),w)) call unsafe_arithmetic()
        if (.not.safe_subtract(ho,hp,den)) call unsafe_arithmetic()
        if(den.eq.0) then 
           write(*,*)'failure in polint' 
           stop 1
        endif
        if (.not.safe_divide(w,den,quotient)) call unsafe_arithmetic()
        if (.not.safe_multiply(hp,quotient,product)) call unsafe_arithmetic()
        d(i)=product
        if (.not.safe_multiply(ho,quotient,product)) call unsafe_arithmetic()
        c(i)=product
     enddo
     if(2*ns.lt.(n-m)) then 
        dy=c(ns+1) 
     else 
        dy=d(ns) 
        ns=ns-1 
     endif
     if (.not.safe_add(y,dy,product)) call unsafe_arithmetic()
     y=product
  enddo
  if (.not.ieee_is_finite(y) .or. .not.ieee_is_finite(dy)) then
     write(*,*) 'failure in polint: non-finite interpolation result'
     stop 1
  endif
contains
  logical function safe_subtract(first,second,value) result(ok)
    double precision,intent(in) :: first,second
    double precision,intent(out) :: value
    value=0d0
    ok=.false.
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    if (second.lt.0d0) then
       if (first.gt.huge(first)+second) return
    elseif (second.gt.0d0) then
       if (first.lt.-huge(first)+second) return
    endif
    value=first-second
    ok=ieee_is_finite(value)
  end function safe_subtract

  logical function safe_add(first,second,value) result(ok)
    double precision,intent(in) :: first,second
    double precision,intent(out) :: value
    ok=safe_subtract(first,-second,value)
  end function safe_add

  logical function safe_multiply(first,second,value) result(ok)
    double precision,intent(in) :: first,second
    double precision,intent(out) :: value
    double precision :: abs_first,abs_second
    value=0d0
    ok=.false.
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    abs_first=abs(first)
    abs_second=abs(second)
    if (abs_first.eq.0d0 .or. abs_second.eq.0d0) then
       ok=.true.
       return
    endif
    if (abs_second.gt.1d0) then
       if (abs_first.gt.huge(first)/abs_second) return
    endif
    value=first*second
    ok=ieee_is_finite(value)
  end function safe_multiply

  logical function safe_divide(numerator,denominator,value) result(ok)
    double precision,intent(in) :: numerator,denominator
    double precision,intent(out) :: value
    double precision :: abs_numerator,abs_denominator
    value=0d0
    ok=.false.
    if (.not.ieee_is_finite(numerator) .or. &
         .not.ieee_is_finite(denominator)) return
    if (denominator.eq.0d0) return
    abs_numerator=abs(numerator)
    abs_denominator=abs(denominator)
    if (abs_denominator.lt.1d0) then
       if (abs_numerator.gt.huge(numerator)*abs_denominator) return
    endif
    value=numerator/denominator
    ok=ieee_is_finite(value)
  end function safe_divide

  subroutine unsafe_arithmetic()
    write(*,*) 'failure in polint: unsafe interpolation arithmetic'
    stop 1
  end subroutine unsafe_arithmetic
END subroutine lh_polint
