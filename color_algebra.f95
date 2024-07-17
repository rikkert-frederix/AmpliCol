module color_algebra
! By calling Tr_full_simplify(result) reduces a sum of products of traces of
! colour matrices (in the fundamental representation). All indices must be
! contracted, so that the results is a single (possibly complex) number.
!
! Tr(0,0,0)         : number of terms in the sum of traces
! Tr(0,0,iterm)     : number traces multiplied in iterm
! Tr(0,iprod,iterm) : number of lambda matrices in the trace
! coef(iterm)       : coefficient multiplying the iterm in the sum of traces
!
! Note: there is no check on the consistency of the input string. If
! the input is not consistent, this package might end up in an
! infinite loop!
!  
! Written by Rikkert Frederix 2018-2019
! Additional functions added by Timea Vitos, 2021
!
! 2019-12-05: Added conversion subroutine for string of matrices in
!             adjoint representation (F matrices) to lambda matrices
!             (fundamental representation).
!
! 2019-12-07: Collected computation of color factor in color flow
!             representation into this package (using the improved
!             version written by Johannes Bellm)
!
! 2019-12-20: Added the computation of the analytic 1/Nc expansion (saved in
!             coef_Nc(:,0)). Also added an alternative 'check_NLC' subroutine
!             that quickly computes the (N)LC contributions in the fundamental
!             basis for all-gluon amplitudes.
!
! 2021       : Added Kronecker delta contraction for 1qq, 2qqDF and 2qqSF cases
!               Added check_NLC for 1qq and 2qqDF 2qqSF cases (Timea)
!
! This package also contains to useful helper functions, 'ipnext' and
! 'get_next_iperm' that can be used to loop over permuation of
! objects.
!
  private
  integer,dimension(:,:,:),allocatable,public :: Tr
  real*16,dimension(:),allocatable,public :: coef
  integer,dimension(:,:),allocatable,public :: coef_Nc
  integer,dimension(:),allocatable,public :: F
  integer,parameter :: Nc=3
  public &
       & Tr_allocate, &         ! allocate memory for lambda matrices (and F-string)
       & Tr_deallocate, &       ! deallocate memory for lambda matrices (and F-string)
       & Tr_full_simplify,&     ! recursively simply the lambda matrices to get colour factor
       & Tr_complex_conjugate,& ! Take the complex conjugate of a string of lambda matrices
       & Tr_print_string,&      ! Print the colour string in human readable format
       & convert_Fs_to_Tr,&     ! Convert a string of F-matrices to lambda matrices
       & color_flow_factor,&    ! Compute color factor in color flow basis
       & color_flow_factor_1qqbar, &
       & color_flow_factor_2qqbar, &
       & color_flow_factor_2qqbar_sf, & 
       & ipnext,&    
       & ipnextgen,&       ! Helper function: get next permutation
       & check_NLC,&            ! gives order in NLC based on permutations for all-gluon
       & check_NLC_1qqbar,&     ! gives order in NLC for 1qq
       & check_NLC_2qqbar, &       ! gives order in NLC based on permutations for 2qq
       & check_NLC_2qqbar_SF
contains

  subroutine Tr_allocate(n)
! if input is in the format Tr(...)*Tr(...) with both traces containing a
! permutation of 'n' different matrices, declares the Tr() and coef() arrays
! to cover the maximum possible size in all intermediate steps.
    implicit none
    integer :: n
    if (allocated(Tr)) deallocate(Tr)
    if (allocated(coef)) deallocate(coef)
    if (allocated(coef_Nc)) deallocate(coef_Nc)
    if (allocated(F)) deallocate(F)
    allocate(coef(2**n))
    allocate(coef_Nc(-n:n,0:2**n))
    allocate(Tr(0:2*(n-1),0:n,0:2**n))
    allocate(F(n-4))
    coef_Nc(:,:)=0
  end subroutine Tr_allocate

  subroutine Tr_deallocate
    implicit none
    if (allocated(Tr)) deallocate(Tr)
    if (allocated(coef)) deallocate(coef)
    if (allocated(coef_Nc)) deallocate(coef_Nc)
    if (allocated(F)) deallocate(F)
  end subroutine Tr_deallocate

  subroutine Tr_full_simplify(res)
! Calls Tr_simplify repeatedly until no traces left. Returns the colour
! factor. Note there is no fail-safe, so might go into an infinite loop if
! input cannot be completely simplified.
    implicit none
    logical :: done
    integer :: iterm
    real*16 :: res
    ! Tr(a)=0
    call simplify_Tr1()
    ! Tr()=Nc
    call simplify_Tr0()
    done =.false.
    do while (.not.done)
       call Tr_simplify()
       done=.true.
       do iterm=1,Tr(0,0,0)
          if (Tr(0,0,iterm).ne.0) then
             done =.false.
             exit
          endif
       enddo
    enddo
    res=0d0
    coef_Nc(:,0)=0
    do iterm=1,Tr(0,0,0)
       res=res+coef(iterm)
       coef_Nc(:,0)=coef_Nc(:,0)+coef_Nc(:,iterm)
    enddo
  end subroutine Tr_full_simplify

  subroutine Tr_complex_conjugate(iprod,iterm)
    ! takes the complex conjugate, i.e., just reverse the order of the
    ! elements in the trace
    implicit none
    integer :: iprod,iterm,i
    integer, allocatable,dimension(:) :: temp
    allocate(temp(1:Tr(0,iprod,iterm)))
    temp(1:Tr(0,iprod,iterm))=Tr(1:Tr(0,iprod,iterm),iprod,iterm)
    do i=1,Tr(0,iprod,iterm)
       Tr(i,iprod,iterm)=temp(Tr(0,iprod,iterm)-i+1)
    enddo
    deallocate(temp)
  end subroutine Tr_complex_conjugate

  subroutine Tr_print_string()
    ! Prints the colour string to the screen in some human-readable format.
    implicit none
    integer :: iterm,iprod,ele
    do iterm=1,Tr(0,0,0)
       write (*,fmt="(a1,f8.4,f8.4,a1)",advance="no") '(',coef(iterm),')'
       do iprod=1,Tr(0,0,iterm)
          write (*,fmt="(a4)",advance="no") "*Tr("
          do ele=1,Tr(0,iprod,iterm)
             if (ele.ne.Tr(0,iprod,iterm)) then
                write (*,fmt="(i2,a1)",advance="no") Tr(ele,iprod,iterm),','
             else
                write (*,fmt="(i2)",advance="no") Tr(ele,iprod,iterm)
             endif
          enddo
          write (*,fmt="(a1)",advance="no") ")"
       enddo
       if (iterm.ne.Tr(0,0,0)) then
          write (*,*) '+'
       else
          write(*,*) ''
       endif
    enddo
  end subroutine Tr_print_string


  subroutine Tr_simplify()
! simplifies a sum of products of traces of T^a_ij colour factors and
! simplifies single elements.
    implicit none
    ! Tr(a,x,b)*Tr(c,x,d)=(Tr(a,d,c,b)-1/Nc*Tr(a,b)*Tr(c,d))
    call Tr_pair_simplify()
    ! Tr(a)=0
    call simplify_Tr1()
    ! Tr()=Nc
    call simplify_Tr0()
!!$    ! Use cyclicity to bring to canonical order [always smallest number first]
!!$    call simplify_cyc()
    ! tr(a,x,b,x,c)=(tr(a,c)*tr(b) - 1/Nc*tr(a,b,c))
    call simplify_Traxbxc()
    ! Tr(a)=0
    call simplify_Tr1()
    ! Tr()=Nc
    call simplify_Tr0()
!!$    ! Use cyclicity to bring to canonical order [always smallest number first]
!!$    call simplify_cyc()
  end subroutine Tr_simplify

  subroutine simplify_Tr1()
    ! The trace of a single lambda matrix is equal to zero
    implicit none
    integer :: iterm,iprod
    logical :: removed
    iterm=1
    removed=.false.
    do while (iterm.le.Tr(0,0,0))
       do iprod=1,Tr(0,0,iterm)
          ! trace of 1 element:
          if (Tr(0,iprod,iterm).eq.1) then
             coef(iterm)=0d0
             coef_Nc(:,iterm)=0
             removed=.true.
             call remove_iterm(iterm)
             exit
          endif
       enddo
       if (removed) then
          ! if a term has been removed, do no increment 'iterm' since all
          ! terms have been shifted by remove_iterm().
          removed=.false.
       else
          iterm=iterm+1
       endif
    enddo
  end subroutine simplify_Tr1

  subroutine update_coef_Nc(op,iterm)
    implicit none
    integer :: op,iterm,il,i
    il=(size(coef_Nc,1)-1)/2
    if (op.eq.1) then
       if (coef_Nc(il,iterm).ne.0) then
          write (*,*) 'out of bounds in updating coef_Nc',il,iterm,op
          stop 1
       endif
       do i=il,-il+op,-1
          coef_Nc(i,iterm)=coef_Nc(i-op,iterm)
       enddo
    elseif (op.eq.-1) then
       if (coef_Nc(-il,iterm).ne.0) then
          write (*,*) 'out of bounds in updating coef_Nc',-il,iterm,op
          stop 1
       endif
       do i=-il-op,il
          coef_Nc(i+op,iterm)=coef_Nc(i,iterm)
       enddo
    else
       write (*,*) 'UNKNOWN OPERATION',op
       stop 1
    endif
  end subroutine update_coef_Nc

  subroutine simplify_Tr0()
    ! The trace of zero lambda matrices (i.e., the trace of the identity
    ! matrix) is equal to Nc
    implicit none
    integer :: iterm,iprod
    do iterm=1,Tr(0,0,0)
       iprod=1
       do while (iprod.le.Tr(0,0,iterm))
          ! trace of zero elements:
          if (Tr(0,iprod,iterm).eq.0) then
             coef(iterm)=coef(iterm)*Nc
             call update_coef_Nc(+1,iterm)
             call remove_iprod(iprod,iterm)
             cycle ! do not increment iprod, since all has been shifted by
                   ! remove_iprod()
          endif
          iprod=iprod+1
       enddo
    enddo
  end subroutine simplify_Tr0

  subroutine simplify_cyc()
    implicit none
    integer :: iterm,iprod
    integer,dimension(1) :: min_loc
    integer,allocatable,dimension(:) :: temp
    do iterm=1,Tr(0,0,0)
       do iprod=1,Tr(0,0,iterm)
          ! Use cyclicity to bring to canonical order [always smallest number
          ! first]
          min_loc=minloc(Tr(1:Tr(0,iprod,iterm),iprod,iterm))
          if (min_loc(1).ne.1) then
             allocate(temp(Tr(0,iprod,iterm)))
             temp(1:Tr(0,iprod,iterm))=Tr(1:Tr(0,iprod,iterm),iprod,iterm)
             Tr(1:Tr(0,iprod,iterm)-min_loc(1)+1,iprod,iterm)= &
                                         temp(min_loc(1):Tr(0,iprod,iterm))
             Tr(Tr(0,iprod,iterm)-min_loc(1)+2:Tr(0,iprod,iterm),iprod,iterm)= &
                                         temp(1:min_loc(1)-1)
             deallocate(temp)
          endif
       enddo
    enddo
  end subroutine simplify_cyc


  subroutine Tr_pair_simplify()
! simplifies products of traces:
! Tr(a,x,b)*Tr(c,x,d)=(Tr(a,d,c,b)-1/Nc*Tr(a,b)*Tr(c,d))
    implicit none
    integer :: iterm,iprod1,iprod2,ele1,ele2,iprod
    iterm=1
    do while (iterm.le.Tr(0,0,0))
       iprod1=1
       do while (iprod1.le.Tr(0,0,iterm))
          iprod2=iprod1+1
          do while (iprod2.le.Tr(0,0,iterm))
             elements : do ele1=1,Tr(0,iprod1,iterm)
                do ele2=1,Tr(0,iprod2,iterm)
                   if (Tr(ele1,iprod1,iterm).eq.Tr(ele2,iprod2,iterm)) then
                      ! Add one extra term corresponding to Tr(a,d,c,b)
                      Tr(0,0,0)=Tr(0,0,0)+1
                      Tr(0,0,Tr(0,0,0))=Tr(0,0,iterm)-1
                      do iprod=1,Tr(0,0,iterm)
                         if (iprod.eq.iprod1) then
                            Tr(0,iprod1,Tr(0,0,0))=Tr(0,iprod1,iterm)-1+Tr(0,iprod2,iterm)-1
                            Tr(1:Tr(0,iprod1,Tr(0,0,0)),iprod1,Tr(0,0,0))= &
                                                [Tr(1:ele1-1,iprod1,iterm), &
                                                 Tr(ele2+1:Tr(0,iprod2,iterm),iprod2,iterm), &
                                                 Tr(1:ele2-1,iprod2,iterm), &
                                                 Tr(ele1+1:Tr(0,iprod1,iterm),iprod1,iterm)]
                         elseif (iprod.eq.iprod2) then
                            cycle
                         elseif (iprod.gt.iprod2) then
                            Tr(0,iprod-1,Tr(0,0,0))=Tr(0,iprod,iterm)
                            Tr(1:Tr(0,iprod-1,Tr(0,0,0)),iprod-1,Tr(0,0,0))= &
                                                 Tr(1:Tr(0,iprod,iterm),iprod,iterm)
                         else
                            Tr(0,iprod,Tr(0,0,0))=Tr(0,iprod,iterm)
                            Tr(1:Tr(0,iprod,Tr(0,0,0)),iprod,Tr(0,0,0))= &
                                                 Tr(1:Tr(0,iprod,iterm),iprod,iterm)
                         endif
                      enddo
                      coef(Tr(0,0,0))=coef(iterm)
                      coef_Nc(:,Tr(0,0,0))=coef_Nc(:,iterm)
                      ! Replace the prod1 and prod2 terms by -1/Nc*Tr(a,b)*Tr(c,d)
                      Tr(0,iprod1,iterm)=Tr(0,iprod1,iterm)-1
                      Tr(1:Tr(0,iprod1,iterm),iprod1,iterm)=[Tr(1:ele1-1,iprod1,iterm), &
                                               Tr(ele1+1:Tr(0,iprod1,iterm)+1,iprod1,iterm)]
                      Tr(0,iprod2,iterm)=Tr(0,iprod2,iterm)-1
                      Tr(1:Tr(0,iprod2,iterm),iprod2,iterm)=[Tr(1:ele2-1,iprod2,iterm), &
                                               Tr(ele2+1:Tr(0,iprod2,iterm)+1,iprod2,iterm)]
                      coef(iterm)=-coef(iterm)/Nc
                      call update_coef_Nc(-1,iterm)
                      coef_Nc(:,iterm)=-coef_Nc(:,iterm)
                      exit elements
                   endif
                enddo
             enddo elements
             iprod2=iprod2+1
          enddo
          iprod1=iprod1+1
       enddo
       iterm=iterm+1
    enddo
  end subroutine Tr_pair_simplify

  subroutine simplify_Traxbxc()
    ! Tr(a,x,b,x,c)=(Tr(a,c)*Tr(b) - 1/Nc*Tr(a,b,c)).  We do this by
    ! replacing the Tr(a,x,b,x,c) by the 1/Nc*Tr(a,b,c) and adding an
    ! additional term equal to Tr(a,c)*Tr(b)
    implicit none
    integer :: iterm,iprod,ele1,ele2,i
    logical :: replace
    replace=.false.
    iterm=1
    do while (iterm.le.Tr(0,0,0)) ! need while loop since Tr(0,0,0) can be updated in the loop
       iprod=1
       do while (iprod.le.Tr(0,0,iterm))
          elements : do ele1=1,Tr(0,iprod,iterm)
             do ele2=ele1+1,Tr(0,iprod,iterm)
                if (Tr(ele1,iprod,iterm).eq.Tr(ele2,iprod,iterm)) then

                   ! First, add an extra term equal to (Tr(a,c)*Tr(b))
                   Tr(0,0,0)=Tr(0,0,0)+1
                   Tr(0,0,Tr(0,0,0))=Tr(0,0,iterm)+1 ! new term has one more product
                   do i=1,Tr(0,0,Tr(0,0,0)) ! loop over all the terms in the product
                      if (i.eq.iprod) then
                         ! put here the Tr(a,c) term
                         Tr(0,i,Tr(0,0,0))=(ele1-1)+(Tr(0,iprod,iterm)-ele2)
                         Tr(1:Tr(0,i,Tr(0,0,0)),i,Tr(0,0,0))= &
                                                   [Tr(1:ele1-1,iprod,iterm), & 
                                                    Tr(ele2+1:Tr(0,iprod,iterm),iprod,iterm)]
                      elseif(i.eq.Tr(0,0,Tr(0,0,0))) then
                         ! add the Tr(b) term
                         Tr(0,i,Tr(0,0,0))=ele2-ele1-1
                         Tr(1:Tr(0,i,Tr(0,0,0)),i,Tr(0,0,0))=Tr(ele1+1:ele2-1,iprod,iterm)
                      else
                         ! just copy the relevant info
                         Tr(0,i,Tr(0,0,0))=Tr(0,i,iterm)
                         Tr(1:Tr(0,i,Tr(0,0,0)),i,Tr(0,0,0))=Tr(1:Tr(0,i,iterm),i,iterm)
                      endif
                   enddo
                   coef(Tr(0,0,0))=coef(iterm)
                   coef_Nc(:,Tr(0,0,0))=coef_Nc(:,iterm)

                   ! Second, replace current term (Tr(a,x,b,x,c) by -1/Nc*Tr(a,b,c))
                   Tr(0,iprod,iterm)=Tr(0,iprod,iterm)-2
                   Tr(1:Tr(0,iprod,iterm),iprod,iterm)=[Tr(1:ele1-1,iprod,iterm), &
                                                        Tr(ele1+1:ele2-1,iprod,iterm), &
                                                        Tr(ele2+1:Tr(0,iprod,iterm)+2,iprod,iterm)]
                   coef(iterm)=-coef(iterm)/Nc
                   call update_coef_Nc(-1,iterm)
                   coef_Nc(:,iterm)=-coef_Nc(:,iterm)
                   replace=.true.
                   exit elements
                endif
             enddo
          enddo elements
          if (replace) then
             replace=.false.
          else
             iprod=iprod+1
          endif
       enddo
       iterm=iterm+1
    enddo
  end subroutine simplify_Traxbxc

  subroutine convert_Fs_to_Tr(nF)
    ! Convert a trace of F-matrices to lambda matrices using:
    ! (F^a)_bc= Tr(b,a,c)-Tr(a,b,c)
    ! 'nF' is the number of F matrices in the trace
    implicit none
    integer :: maxl,i,iterm,nF
    ! number of terms in the string of Fs: nF
    if (Tr(0,0,0).ne.0) then
       write (*,*) 'Can only convert Fs to Tr when there is no Tr() yet'
       stop 1
    endif
    maxl=maxval(F(1:nF))
    do i=1,nF
       if (i.eq.1) then
          ! create the first traces of fundamental matrices
          Tr(0,0,0)=2
          coef(2)=-coef(1)
!          write(*,*) coef(1)
          coef_Nc(:,2)=-coef_Nc(:,1)
          Tr(0,0,1)=1
          Tr(0,0,2)=1
          Tr(0,1,1)=3
          Tr(0,1,2)=3
          Tr(1:3,1,1)=(/ maxl+i , F(i) , maxl+i+1 /) ! 'Tr(b,a,c)'
          Tr(1:3,1,2)=(/ F(i) , maxl+i , maxl+i+1 /) ! 'Tr(a,b,c)'
        else
          ! Multiply the existing traces of fundamental matrices by
          ! '[Tr(b,a,c)-Tr(a,b,c)]':
          call double_trs()
          do iterm=1,Tr(0,0,0)
             Tr(0,0,iterm)=Tr(0,0,iterm)+1
             Tr(0,Tr(0,0,iterm),iterm)=3
             if (iterm.le.Tr(0,0,0)/2) then
                if (i.ne.nF) then
                   ! The new trace that's added is Tr(b,a,c)
                   Tr(1:3,Tr(0,0,iterm),iterm)=(/ maxl+i , F(i) , maxl+i+1 /)
                else
                   ! Final F to close the Tr(F..F): maxl+i+1 must now
                   ! be equal to the first new label, i.e., maxl+1
                   Tr(1:3,Tr(0,0,iterm),iterm)=(/ maxl+i , F(i) , maxl+1 /)
                endif
              else
                if (i.ne.nF) then
                   ! The new trace that's added is Tr(a,b,c)
                   Tr(1:3,Tr(0,0,iterm),iterm)=(/ F(i) , maxl+i , maxl+i+1 /)
                else
                   ! Final F to close the Tr(F..F): maxl+i+1 must now
                   ! be equal to the first new label, i.e.,  maxl+1
                   Tr(1:3,Tr(0,0,iterm),iterm)=(/ F(i) , maxl+i , maxl+1 /)
                endif
                coef(iterm)=-coef(iterm)
                coef_Nc(:,iterm)=-coef_Nc(:,iterm)
             endif
          enddo
       endif
       ! A quick Simplify every time to not let the number of traces
       ! of fundamental matrices run out of control
       call Tr_simplify()
       call check_identical_and_remove()
    enddo
  end subroutine convert_Fs_to_Tr

  subroutine check_identical_and_remove()
    ! Searches for identical terms in the sum and combines them
    implicit none
    integer iterm,jterm
    logical :: identical
    if (Tr(0,0,0).le.1) return
    do iterm=Tr(0,0,0),2,-1
       do jterm=iterm-1,1,-1
          call check_identical(iterm,jterm,identical)
          if (identical) then
             coef(jterm)=coef(jterm)+coef(iterm)
             coef_Nc(:,jterm)=coef_Nc(:,jterm)+coef_Nc(:,iterm)
             call remove_iterm(iterm)
             exit
          endif
       enddo
    enddo
    do iterm=Tr(0,0,0),1,-1
       if (coef(iterm).eq.0d0) then
          if (all(coef_Nc(:,iterm).eq.0)) then
             call remove_iterm(iterm)
          endif
       endif
    enddo
  end subroutine check_identical_and_remove
  
  subroutine check_identical(iterm,jterm,identical)
    ! Check if two terms are identical
    implicit none
    logical :: identical
    integer :: iterm,jterm,iprod,ele1
    identical = .false.
    if (Tr(0,0,iterm).ne.Tr(0,0,jterm)) return
    do iprod=1,Tr(0,0,iterm)
       if (Tr(0,iprod,iterm).ne.Tr(0,iprod,jterm)) return
       do ele1=1,Tr(0,iprod,iterm)
          if (Tr(ele1,iprod,iterm).ne.Tr(ele1,iprod,jterm)) return
       enddo
    enddo
    identical=.True.
  end subroutine check_identical
  
  subroutine double_trs()
    ! Simply dublicate all the terms
    implicit none
    integer :: iterm,iprod,ele1
    do iterm=1,Tr(0,0,0)
       Tr(0,0,Tr(0,0,0)+iterm)=Tr(0,0,iterm)
       do iprod=1,Tr(0,0,iterm)
          Tr(0,iprod,Tr(0,0,0)+iterm)=Tr(0,iprod,iterm)
          do ele1=1,Tr(0,iprod,iterm)
             Tr(ele1,iprod,Tr(0,0,0)+iterm)=Tr(ele1,iprod,iterm)
          enddo
       enddo
       coef(Tr(0,0,0)+iterm)=coef(iterm)
       coef_Nc(:,Tr(0,0,0)+iterm)=coef_Nc(:,iterm)
    enddo
    Tr(0,0,0)=2*Tr(0,0,0)
  end subroutine double_trs
  
  subroutine remove_iterm(iterm)
    implicit none
    integer :: iterm,i,j
    ! remove one term in the sum of (products of) traces
    Tr(0,0,0)=Tr(0,0,0)-1
    ! shift everything else
    do i=iterm,Tr(0,0,0)
       Tr(0,0,i)=Tr(0,0,i+1)
       do j=1,Tr(0,0,i)
          Tr(0,j,i)=Tr(0,j,i+1)
          Tr(1:Tr(0,j,i),j,i)=Tr(1:Tr(0,j,i),j,i+1)
       enddo
       ! also shift the coefficients multiplying the terms
       coef(i)=coef(i+1)
       coef_Nc(:,i)=coef_Nc(:,i+1)
    enddo
  end subroutine remove_iterm

  subroutine remove_iprod(iprod,iterm)
    implicit none
    integer :: iprod,iterm,i
    Tr(0,0,iterm)=Tr(0,0,iterm)-1 ! one less product in iterm
    do i=iprod,Tr(0,0,iterm)      ! shift the rest
       Tr(0,i,iterm)=Tr(0,i+1,iterm)
       Tr(1:Tr(0,i,iterm),i,iterm)=Tr(1:Tr(0,i,iterm),i+1,iterm)
    enddo
  end subroutine remove_iprod


  subroutine color_flow_factor(n,iper,jper,col_fac)
    ! Given two permutations (iper and jper) compute corresponding
    ! colour factor using the color flow decomposition (i.e., string
    ! of kronecker delta's -- eq.3 of hep-ph/0209271 [Maltoni, Paul,
    ! Stelzer & Willenbrock]).
    ! It sets up the string of delta's, contracts all the indices,
    ! counts the closed loops ('colfac'), and set the color factor
    ! equal to Nc**colfac.
    !
    ! Improved version, by Johannes Bellm.
    !
    implicit none
    integer :: n                 ! number of gluons
    integer,dimension(n) :: iper ! order in amplitude
    integer,dimension(n) :: jper ! order in conjugate amplitude
    integer :: col_fac           ! color factor
    integer :: index,colfac,i,flip,m,k
    integer, dimension(n) :: it,jt
    do i=1,n
       it(iper(mod(i,n)+1))=iper(i)
       jt(jper(i))=jper(mod(i,n)+1)
    enddo
    index=1
    colfac=-1
    flip=0
    do
       do
          if (index.gt.n) exit
          if (it(index).ge.1) exit
          index=index+1
       enddo
       if (index .ge. n) exit
       k=index
       do while (it(k) .ge. 1)
          m=k
          k=jt(it(k))
          it(m)=colfac
          flip=flip+1
       enddo
       if (flip.lt.n) colfac=colfac-1
    enddo
!!$    col_fac=Nc**(-colfac)
    col_fac=-colfac
  end subroutine color_flow_factor

  subroutine color_flow_factor_1qqbar(n,iper,jper,ri,rj,col_fac)
    ! Given two permutations (iper and jper) compute corresponding
    ! colour factor using the color flow decomposition (i.e., string
    ! of kronecker delta's -- eq.3 of hep-ph/0209271 [Maltoni, Paul,
    ! Stelzer & Willenbrock]).
    ! It sets up the string of delta's, contracts all the indices,
    ! counts the closed loops ('colfac'), and set the color factor
    ! equal to Nc**colfac.
    !
    ! by Timea Vitos
    !
    implicit none
    integer :: n                 ! number of external particles
    integer,dimension(n-2) :: iper ! order in amplitude
    integer,dimension(n-2) :: jper ! order in conjugate amplitude
    integer,dimension(n-1) :: iper_q_up,iper_q_down,iper_q ! order in amplitude
    integer,dimension(n-1) :: jper_q_up,jper_q_down,jper_q
    integer :: col_fac           ! color factor
    integer :: index,colfac,i,flip,m,k
    integer :: ri,rj
    integer, dimension(n-1) :: itemp

    iper_q_up(1) = n-1
    iper_q_up(2:n-1) = iper(1:n-2)

    jper_q_up(1) = n-1
    jper_q_up(2:n-1) = jper(1:n-2)

    iper_q_down(1:ri) = iper(1:ri)
    iper_q_down(ri+1) = n-1
    iper_q_down(ri+2:n-1) = iper(ri+1:n-2)

    jper_q_down(1:rj) = jper(1:rj)
    jper_q_down(rj+1) = n-1
    jper_q_down(rj+2:n-1) = jper(rj+1:n-2)
    
    itemp = jper_q_down
    jper_q_down = jper_q_up
    jper_q_up = itemp

    do i=1,size(iper_q_down)
       if (iper_q_down(i) .eq. i) then
               iper_q(i) = iper_q_up(i)
       else 
          iper_q(iper_q_down(i)) = iper_q_up(i)
       endif
    enddo

    do i=1,size(jper_q_down)
       if (jper_q_down(i) .eq. i) then
               jper_q(i) = jper_q_up(i)
       else
          jper_q(jper_q_down(i)) = jper_q_up(i)
       endif
    enddo
    
    index=1
    colfac=-1
    flip=0
    do
       do
          if (index.gt.n-1) exit
          if (iper_q(index).ge.1) exit
          index=index+1
       enddo
       if (index .ge. n-1) exit
       k=index
       do while (iper_q(k) .ge. 1)
          m=k
          k=jper_q(iper_q(k))
          iper_q(m)=colfac
          flip=flip+1
       enddo
       if (flip.lt.n-1) colfac=colfac-1
    enddo
!!$    col_fac=Nc**(-colfac)    
    col_fac=-colfac
    col_fac = col_fac - (n-2-ri) - (n-2-rj)
    col_fac = col_fac
    
  end subroutine color_flow_factor_1qqbar



  subroutine color_flow_factor_2qqbar(n,iper_t,jper_t,ni_t,nj_t,ii_t,jj_t,col_fac)
    ! Given two permutations (iper and jper) compute corresponding
    ! colour factor using the color flow decomposition (i.e., string
    ! of kronecker delta's -- eq.3 of hep-ph/0209271 [Maltoni, Paul,
    ! Stelzer & Willenbrock]).
    ! It sets up the string of delta's, contracts all the indices,
    ! counts the closed loops ('colfac'), and set the color factor
    ! equal to Nc**colfac.
    !
    ! by Timea Vitos
    !
    implicit none
    integer :: n                 ! number of external particles
    integer,dimension(n-3) :: iper_t ! order in amplitude
    integer,dimension(n-3) :: jper_t ! order in conjugate amplitude
    integer,dimension(n-2) :: iper_q_up,iper_q_down,iper_q ! order in amplitude
    integer,dimension(n-2) :: jper_q_up,jper_q_down,jper_q
    integer :: col_fac           ! color factor
    integer :: index,colfac,u,flip,m,w
    integer :: ii_t,jj_t,ni_t,nj_t
    integer, dimension(n-2) :: itemp

    if (ii_t .eq. 1) then

    iper_q_up(1) = n-2
    iper_q_up(2:n-2) = iper_t(1:n-3)

    iper_q_down(1:ni_t+1) = iper_t(1:ni_t+1)
    iper_q_down(ni_t+2) = n-2
    iper_q_down(ni_t+3:n-2) = iper_t(ni_t+2:n-3)


    elseif (ii_t .eq. 2) then
       
    iper_q_up(1) = n-2
    iper_q_up(2:n-2) = iper_t(1:n-3)

    do u = 1,ni_t+1
    if (iper_t(u) .eq. n-3) then
            iper_q_down(u) = n-2
    else
            iper_q_down(u) = iper_t(u)
    endif
    enddo

    iper_q_down(ni_t+2) = n-3
    iper_q_down(ni_t+3:n-2) = iper_t(ni_t+2:n-3)

    endif

    if (jj_t .eq. 1) then

    jper_q_up(1) = n-2
    jper_q_up(2:n-2) = jper_t(1:n-3)

    jper_q_down(1:nj_t+1) = jper_t(1:nj_t+1)
    jper_q_down(nj_t+2) = n-2
    jper_q_down(nj_t+3:n-2) = jper_t(nj_t+2:n-3)

    elseif (jj_t .eq. 2) then

    jper_q_up(1) = n-2
    jper_q_up(2:n-2) = jper_t(1:n-3)

    do u = 1,nj_t+1
    if (jper_t(u) .eq. n-3) then
            jper_q_down(u) = n-2
    else
            jper_q_down(u) = jper_t(u)
    endif
    enddo

    jper_q_down(nj_t+2) = n-3
    jper_q_down(nj_t+3:n-2) = jper_t(nj_t+2:n-3)

    endif 

    itemp = jper_q_down
    jper_q_down = jper_q_up
    jper_q_up = itemp

    do u=1,size(iper_q_down)
       if (iper_q_down(u) .eq. u) then
               iper_q(u) = iper_q_up(u)
       else
          iper_q(iper_q_down(u)) = iper_q_up(u)
       endif
    enddo

    do u=1,size(jper_q_down)
       if (jper_q_down(u) .eq. u) then
               jper_q(u) = jper_q_up(u)
       else
          jper_q(jper_q_down(u)) = jper_q_up(u)
       endif
    enddo

    index=1
    colfac=-1
    flip=0
    do
       do
          if (index.gt.n-2) exit
          if (iper_q(index).ge.1) exit
          index=index+1
       enddo
       if (index .ge. n-2) exit
       w=index
       do while (iper_q(w) .ge. 1)
          m=w
          w=jper_q(iper_q(w))
          iper_q(m)=colfac
          flip=flip+1
       enddo
       if (flip.lt.n-2) colfac=colfac-1
    enddo
!!$    col_fac=Nc**(-colfac)    
    col_fac=-colfac
    col_fac = col_fac - (n-4-ni_t) - (n-4-nj_t)+(1-ii_t)+ (1-jj_t)



  end subroutine color_flow_factor_2qqbar




  subroutine color_flow_factor_2qqbar_sf(n,iper,jper,ni,nj,ii,jj,col_fac)
    ! Given two permutations (iper and jper) compute corresponding
    ! colour factor using the color flow decomposition (i.e., string
    ! of kronecker delta's -- eq.3 of hep-ph/0209271 [Maltoni, Paul,
    ! Stelzer & Willenbrock]).
    ! It sets up the string of delta's, contracts all the indices,
    ! counts the closed loops ('colfac'), and set the color factor
    ! equal to Nc**colfac.
    !
    ! by Timea Vitos
    !

    implicit none
    integer :: n                 ! number of external particles
    integer,dimension(n-3) :: iper ! order in amplitude
    integer,dimension(n-3) :: jper ! order in conjugate amplitude
    integer,dimension(n-2) :: iper_q_up,iper_q_down,iper_q ! order in amplitude
    integer,dimension(n-2) :: jper_q_up,jper_q_down,jper_q
    integer :: col_fac           ! color factor
    integer :: index,colfac1,i,flip,m,k
    integer :: ii,jj,ni,nj
    integer, dimension(n-2) :: itemp

!!! It is assumed that the last next-4-ri/rj indices of iper,jper
!!! are the U(1) gluon indices, so be careful with the input!!!! 

    if (ii .eq. 1) then

    iper_q_up(1) = n-2
    iper_q_up(2:n-2) = iper(1:n-3)

    iper_q_down(1:ni+1) = iper(1:ni+1)
    iper_q_down(ni+2) = n-2
    iper_q_down(ni+3:n-2) = iper(ni+2:n-3)


    elseif (ii .eq. 2) then

    iper_q_up(1) = n-2
    iper_q_up(2:n-2) = iper(1:n-3)

    do i = 1,ni+1
    if (iper(i) .eq. n-3) then
            iper_q_down(i) = n-2
    else
            iper_q_down(i) = iper(i)
    endif
    enddo


    iper_q_down(ni+2) = n-3
    iper_q_down(ni+3:n-2) = iper(ni+2:n-3)


    endif

    if (jj .eq. 1) then

    jper_q_up(1) = n-2
    jper_q_up(2:n-2) = jper(1:n-3)

    jper_q_down(1:nj+1) = jper(1:nj+1)
    jper_q_down(nj+2) = n-2
    jper_q_down(nj+3:n-2) = jper(nj+2:n-3)


    elseif (jj .eq. 2) then

    jper_q_up(1) = n-2
    jper_q_up(2:n-2) = jper(1:n-3)

    do i = 1,nj+1
    if (jper(i) .eq. n-3) then
            jper_q_down(i) = n-2
    else
            jper_q_down(i) = jper(i)
    endif
    enddo

    jper_q_down(nj+2) = n-3
    jper_q_down(nj+3:n-2) = jper(nj+2:n-3)

    endif

    itemp = jper_q_down
    jper_q_down = jper_q_up
    jper_q_up = itemp

    do i=1,size(iper_q_down)
       if (iper_q_down(i) .eq. i) then
               iper_q(i) = iper_q_up(i)
       else
          iper_q(iper_q_down(i)) = iper_q_up(i)
       endif
    enddo

    do i=1,size(jper_q_down)
       if (jper_q_down(i) .eq. i) then
               jper_q(i) = jper_q_up(i)
       else
          jper_q(jper_q_down(i)) = jper_q_up(i)
       endif
    enddo

    index=1
    colfac1=-1
    flip=0
    do
       do
          if (index.gt.n-2) exit
          if (iper_q(index).ge.1) exit
          index=index+1
       enddo
       if (index .ge. n-2) exit
       k=index
       do while (iper_q(k) .ge. 1)
          m=k
          k=jper_q(iper_q(k))
          iper_q(m)=colfac1
          flip=flip+1
       enddo
       if (flip.lt.n-2) colfac1=colfac1-1
    enddo
    col_fac=-colfac1

    col_fac = col_fac - (n-4-ni) - (n-4-nj)
  
  end subroutine color_flow_factor_2qqbar_sf




  subroutine check_NLC(n,iper,jper,acc)
    ! alternative to compute only the (N)LC terms for all-gloun amplitude
    ! squared.
    implicit none
    integer :: n                 ! number of gluons
    integer,dimension(n) :: iper ! order in amplitude
    integer,dimension(n) :: jper ! order in conjugate amplitude
    integer :: acc               ! is equal to 0,1,-1 or 99.
                                 ! 99 : LC contributions (NLC coefficient of that term is '-n')
                                 ! 1,-1 ; NLC contribution with positive/negative sign
                                 ! 0 : not a NLC contribution, but NNLC or further suppressed.
    integer :: i,i1,i2,i3,i4,i5,sign
    integer,dimension(n) :: itemp
    itemp=0
    acc=0
    ! find i1, i.e. the location where jper and iper start to differ
    do i=1,n
       if (jper(i).ne.iper(i)) exit
    enddo
    if (i.eq.n+1) then
       acc=99
       return
    endif
    i1=i
    ! find i2, i.e. the location in jper, such that it's equal to the iper(i1)
    ! (i.e., the iper(i) which no longer equal to jper(i))
    do i=i1+1,n
       if (jper(i).eq.iper(i1)) exit
    enddo
    i2=i
    ! find the max i3 and i4 such that jper(i2:i3) == iper(i1:i4)
    do i=1,n-i2
       if (jper(i2+i).ne.iper(i1+i)) exit
    enddo
    i3=i2+i-1
    i4=i1+i-1
    ! start from jper(i1), and find where it's equal to iper(i4+1)
    do i=i1,n
       if (jper(i).eq.iper(i4+1)) exit
    enddo
    i5=i
    if (i5.gt.i3) return
    ! hence jper(i2:i3) should be switched with jper(i1:i5-1) [which are not
    ! necessarily of the same size!]
    ! Do the swtiching:
    sign=1
    itemp(1:i1-1)=jper(1:i1-1)
    itemp(i1:i4)=jper(i2:i3)
    itemp(i4+1:i4+i2-i5)=jper(i5:i2-1)
    if (i1.gt.i5-1) then
       ! check when switching neighbouring (strings of) elements
       if (i4-i1.eq.0 .and. i2-1-i5.eq.0) then
          ! switch of neighbouring single elements. Here we need a minus sign
          sign=-1
          continue
       elseif ((i4-i1.eq.0   .and. i2-1-i5.eq.n-3) .or. &
            &  (i4-i1.eq.n-3 .and. i2-1-i5.eq.0  )) then
          ! switch of all elements: equal to single neighbouring elements +
          ! cyclicity
          sign=-1
          continue
       elseif (i4-i1.eq.0 .or. i2-1-i5.eq.0) then
          ! not a NLC contribution. Don't switch a single string with string
          return
       elseif (i4-i1 + i2-1-i5 .le. n-4) then
          ! a potential NLC contribution. Not many elements switched
          continue
       else
          ! not a NLC contribution
          return
       endif
    endif
    itemp(i4+i2-i5+1:i4+i2-i1)=jper(i1:i5-1)
    itemp(i4+i2-i1+1:n)=jper(i3+1:n)
    ! Did all the switching. If equal to iper, then it's an NLC contribution
    if (all(itemp.eq.iper)) then
       acc=sign
    endif
  end subroutine check_NLC


  subroutine check_NLC_1qqbar(next,iper,jper,acc)
    implicit none
    integer :: next                 ! number of external particles

    integer,dimension(next-2) :: iper ! order in amplitude
    integer,dimension(next-2) :: jper ! order in conjugate amplitude
    integer :: acc               ! is equal to 0,1,-1 or 99.
                                 ! 99 : LC contributions (NLC coefficient of that term is '-n')
                                 ! 1,-1 ; NLC contribution with positive/negative sign
                                 ! 0 : not a NLC contribution, but NNLC or further suppressed.
    integer :: i,i1,i2,i3,i4,i5,sign

    integer,dimension(next-2) :: itemp


    acc = 0
          do i=1,next-2
            if (jper(i).ne.iper(i)) exit
          enddo
          i1=i
          if (i1 .gt. size(iper)) then
                  acc = 99
                  return
          endif
          do i=i1+1,next-2
            if (jper(i).eq.iper(i1)) exit
          enddo
          i2=i
          do i=1,next-2-i2
             if (jper(i2+i).ne.iper(i1+i)) exit
          enddo
          i3=i2+i-1
          i4=i1+i-1
          do i=i1,next-2
            if (jper(i).eq.iper(i4+1)) exit
          enddo
          i5=i
           if (i5.gt.i3) return
            sign=1
            itemp(1:i1-1)=jper(1:i1-1)
            itemp(i1:i4)=jper(i2:i3)
            itemp(i4+1:i4+i2-i5)=jper(i5:i2-1)
             if (i1.gt.i5-1) then
               if (i4-i1.eq.0 .and. i2-1-i5.eq.0) then
                   sign=-1
                   continue
               elseif (i4-i1.eq.0 .or. i2-1-i5.eq.0) then
                   return
               elseif (i4-i1 + i2-1-i5 .le. next-2-4) then
                   continue
               else
                   sign = 1
             endif
           endif
           itemp(i4+i2-i5+1:i4+i2-i1)=jper(i1:i5-1)
           itemp(i4+i2-i1+1:next-2)=jper(i3+1:next-2)
           if ((all(itemp.eq.iper))) then
             acc=sign
          endif


    end subroutine check_NLC_1qqbar  


  subroutine check_NLC_2qqbar(next,iper,jper,rri,rrj,iii,jjj,acc)
    implicit none
    integer :: next                 ! number of external particles
 
    integer,dimension(next-4) :: iper ! order in amplitude
    integer,dimension(next-4) :: jper ! order in conjugate amplitude
    integer :: acc               ! is equal to 0,1,-1 or 99.
                                 ! 99 : LC contributions (NLC coefficient of that term is '-n')
                                 ! 1,-1 ; NLC contribution with positive/negative sign
                                 ! 0 : not a NLC contribution, but NNLC or further suppressed.
    integer :: i,i1,i2,i3,i4,i5,sign,ii,jj,yy
    integer rri,rrj,iii,jjj
    integer,dimension(2*(next-4)) :: temp
    integer,dimension(next-4,2) :: index_i

    integer, dimension(rri+rrj) :: temp2
    integer, dimension(2*(next-4)-rri-rrj) :: temp3

    integer,dimension(rri) :: itemp
    integer,dimension(next-4-rrj) :: itemp2
    
    integer,dimension(1:rri-1) :: itemp4
    integer,dimension(1:next-4-rrj-1) :: itemp5
    integer,dimension(1:rrj-1) :: itemp6
    integer,dimension(1:next-4-rri-1) :: itemp7

    logical disjoint
    integer :: ind_i,ind_j
    integer,dimension(2*(next-4)-2) :: perm
    integer skipped

    integer, dimension(rri) :: Aa
    integer, dimension(next-4-rri) :: Bb
    integer, dimension(rrj) :: Cc
    integer, dimension(next-4-rrj) :: Dd

    itemp=0

    acc = 0

    if ((iii .eq. 2) .and. (jjj .eq. 2)) then
        if (rri .eq. rrj) then
               do ii=1,next-4
                 if (iper(ii) .eq. jper(ii)) then
                    continue  
                 else
                     acc = 0
                     return
                 endif
               enddo
               acc = 1
        else 
           acc = 0
           return
        endif
    endif  


    if (((iii .eq. 2) .and. (jjj .eq. 1)) .or. ((iii .eq. 1) .and. (jjj .eq. 2))) then

        temp(1:rri)= iper(1:rri)
        temp(rri+1:rri+(next-4-rrj)) = jper(next-4:rrj+1:-1)
        temp(rri+(next-4-rrj)+1:2*(next-4)-rrj) = iper(rri+1:next-4)
        temp(2*(next-4)-rrj+1:2*(next-4)) = jper(rrj:1:-1)

       do jj = 1,next-4
         yy = 1
          do ii=1,2*(next-4) 
             if (temp(ii) .eq. jj) then
                  index_i(jj,yy) = ii
                  yy=yy+1
             endif
          enddo 

         if (mod(abs(index_i(jj,1)-index_i(jj,2)),2) .eq. 1) then
             continue
         else
             acc = 0
             return
         endif
       enddo


     do jj = 1,next-4
       do ii = jj+1,next-4
          if ( (((  (  ( ((index_i(ii,1) .lt. index_i(jj,1)) .and. (index_i(ii,2) .lt. index_i(jj,1))) ) .or.  &
                  (  ( ((index_i(ii,1) .gt. index_i(jj,2)) .and. (index_i(ii,2) .gt. index_i(jj,2)))) .or.  &
                  (  ( ((index_i(ii,1) .lt. index_i(jj,2)) .and. (index_i(ii,2) .lt. index_i(jj,2)))) .and. &
                    (  ((index_i(ii,1) .gt. index_i(jj,1)) .and. (index_i(ii,2) .gt. index_i(jj,1)))   ))))))) .or. &
                    (((  (  ( ((index_i(jj,1) .lt. index_i(ii,1)) .and. (index_i(jj,2) .lt. index_i(ii,1))) ) .or.  &
                  (  ( ((index_i(jj,1) .gt. index_i(ii,2)) .and. (index_i(jj,2) .gt. index_i(ii,2)))) .or.  &
                  (  ( ((index_i(jj,1) .lt. index_i(ii,2)) .and. (index_i(jj,2) .lt. index_i(ii,2)))) .and. &
                    (  ((index_i(jj,1) .gt. index_i(ii,1)) .and. (index_i(jj,2) .gt. index_i(ii,1)))   ))))))) ) then
              continue
          else
               acc = 0 
               return
          endif
       enddo
     enddo   
      acc = -1     
    endif


   if ((iii .eq. 1) .and. (jjj .eq. 1)) then

        if ((rri .eq. rrj) .and. ((all(iper .eq. jper)) .or. ((size(iper).eq.0) &
                .and. (size(jper).eq.0)))) then
               acc = 99
               return
        endif


      Aa(1:rri) = iper(1:rri)
      Bb(1:next-4-rri) = iper(rri+1:next-4)
      Cc(1:rrj) = jper(1:rrj)
      Dd(1:next-4-rrj) = jper(rrj+1:next-4)

      temp2(1:rri) = Aa(1:rri)
      temp2(rri+1:rri+rrj) = Cc(rrj:1:-1)
      
      temp3(1:next-4-rri) = Bb(1:next-4-rri)
      temp3(next-4-rri+1:2*(next-4)-rrj-rri) = Dd(next-4-rrj:1:-1)

      disjoint = .true.

    do ii=1,rri+rrj
       if (.not. any(temp3(:) == temp2(ii))) then 
           disjoint = .true.
       else 
           disjoint = .false.
           exit
       endif
    enddo

    if (disjoint) then

          if ((size(Aa) .ne. 0) .and. (all(Bb .eq. Dd) .or. size(Bb) .eq. 0)) then
 
          do i=1,rri
            if (CC(i).ne.Aa(i)) exit
          enddo
          i1=i
          if (i1 .gt. size(AA)) then
                  acc = 99
                  return
          endif
          do i=i1+1,rri
            if (Cc(i).eq.Aa(i1)) exit
          enddo
          i2=i
          do i=1,rri-i2
             if (Cc(i2+i).ne.Aa(i1+i)) exit
          enddo
          i3=i2+i-1
          i4=i1+i-1
          do i=i1,rri
            if (cc(i).eq.Aa(i4+1)) exit
          enddo
          i5=i
           if (i5.gt.i3) return
            sign=1
            itemp(1:i1-1)=Cc(1:i1-1)
            itemp(i1:i4)=Cc(i2:i3)
            itemp(i4+1:i4+i2-i5)=Cc(i5:i2-1)
             if (i1.gt.i5-1) then
               if (i4-i1.eq.0 .and. i2-1-i5.eq.0) then
                   sign=-1
                   continue
               elseif (i4-i1.eq.0 .or. i2-1-i5.eq.0) then
                   return
               elseif (i4-i1 + i2-1-i5 .le. next-4-4) then
                   continue
               else
                   sign = 1
             endif
           endif
           itemp(i4+i2-i5+1:i4+i2-i1)=Cc(i1:i5-1)
           itemp(i4+i2-i1+1:rri)=Cc(i3+1:rri)
           if ((all(itemp.eq.Aa))) then
             acc=sign
          endif 
 
          elseif ((size(Bb) .ne. 0) .and. (all(Aa .eq. Cc) .or. size(Aa) .eq. 0)) then 

          do i=1,next-4-rrj
             if (Dd(i).ne.Bb(i)) exit
          enddo
          i1=i
          if (i1 .gt. size(Bb)) then
                  acc = 99
                  return
          endif
          do i=i1+1,next-4-rrj
            if (Dd(i).eq.Bb(i1)) exit
          enddo
          i2=i
          do i=1,next-4-rrj-i2
             if (Dd(i2+i).ne.Bb(i1+i)) exit
          enddo
          i3=i2+i-1
          i4=i1+i-1
          do i=i1,next-4-rrj
            if (Dd(i).eq.Bb(i4+1)) exit
          enddo
          i5=i
           if (i5.gt.i3) return
            sign=1
            itemp2(1:i1-1)=Dd(1:i1-1)
            itemp2(i1:i4)=Dd(i2:i3)
            itemp2(i4+1:i4+i2-i5)=Dd(i5:i2-1)
             if (i1.gt.i5-1) then
               if (i4-i1.eq.0 .and. i2-1-i5.eq.0) then
                   sign=-1
                   continue
               elseif (i4-i1.eq.0 .or. i2-1-i5.eq.0) then
                   return
               elseif (i4-i1 + i2-1-i5 .le. next-4-rrj-4) then
                   continue
               else
                 sign = 1
               endif
           endif
         itemp2(i4+i2-i5+1:i4+i2-i1)=Dd(i1:i5-1)
         itemp2(i4+i2-i1+1:next-4-rrj)=Dd(i3+1:next-4-rrj)
         if ((all(itemp2.eq.Bb))) then
            acc=sign
         endif
         else
             acc = 0
         endif



    elseif (.not. disjoint) then

      do ii=1,rri+rrj
       if (any(temp3(:) == temp2(ii))) then
           do jj=1,size(temp3)
              if (temp3(jj) .eq. temp2(ii)) then 
                        ind_i = ii
                        ind_j = jj
                        skipped = temp2(ii)
                        goto 1
              endif
           enddo
       endif
      enddo


1   if ((size(temp2) .eq. 1) .or. (size(temp3) .eq. 1)) then
            acc = 0 
            return
    endif

    perm=0

    !!!   If common generator in A-D pair    
     if ((ind_i .le. rri) .and. (ind_j .gt. next-4-rri)) then
    
      perm(1:rri-ind_i) = temp2(ind_i+1:rri)
      perm(rri-ind_i+1:rri-ind_i+rrj) = temp2(rri+1:rri+rrj)
      perm(rri+rrj-ind_i+1:rri+rrj-1) = temp2(1:ind_i-1)

      perm(rri+rrj: -1+ 2*(next-4)-ind_j  ) = temp3(ind_j+1:2*(next-4)-rri-rrj) 
      perm( 2*(next-4)-ind_j   :  -1+ 2*(next-4)-ind_j + next-4-rri   ) =temp3(1:next-4-rri)
      perm(2*(next-4)-ind_j + next-4-rri:2*(next-4)-2) = temp3(next-4-rri+1:ind_j-1)

      itemp4(1:ind_i-1) = temp2(1:ind_i-1)
      itemp4(ind_i:rri-1) = temp2(ind_i+1:rri)

      itemp5(1:2*(next-4)-rri-rrj-ind_j) = temp3(2*(next-4)-rri-rrj:ind_j+1:-1)
      itemp5(2*(next-4)-rri-rrj-ind_j+1:next-4-rrj-1) = temp3(ind_j-1:next-4-rri+1:-1)

      if (rrj .eq. rri-1) then
              if (all(itemp4 .eq. Cc)) then
                      acc = 0
                      return
              endif
      endif

      if (rrj+1 .eq. rri) then
              if (all(itemp5 .eq. Bb)) then
                      acc = 0
                      return
              endif
      endif

   

      !!!    If common generator in B-C pair      
      elseif  ((ind_i .gt. rri) .and. (ind_j .le. next-4-rri)) then   
        perm(1:rri+rrj-ind_i) = temp2(ind_i+1:rri+rrj)
        perm(rri+rrj-ind_i+1:rri+rrj-ind_i+rri) = temp2(1:rri)
        perm(rri+rrj-ind_i+rri+1:rri+rrj-1) = temp2(rri+1:ind_i-1)

        perm(rri+rrj:rrj-1+next-4-ind_j) = temp3(ind_j+1:next-4-rri)       
        perm(rrj-1+next-4-ind_j+1:-1-ind_j+2*(next-4)) = temp3(next-4-rri+1:2*(next-4)-rri-rrj)
        perm(-1-ind_j+2*(next-4)+1: 2*(next-4)-2) = temp3(1:ind_j-1)

        itemp6(1:rri+rrj-ind_i) = temp2(rri+rrj:ind_i+1:-1)
        itemp6(rri+rrj-ind_i+1:rrj-1) = temp2(ind_i-1:rri+1:-1)

        itemp7(1:ind_j-1) = temp3(1:ind_j-1)
        itemp7(ind_j -1 +1 : next-4-rri-1) = temp3(ind_j+1:next-4-rri)


      if (rrj-1 .eq. rri) then
              if (all(itemp6 .eq. Aa)) then
                      acc = 0 
                      return
              endif
      endif

      if (rrj .eq. rri+1) then
              if (all(itemp7 .eq. Dd)) then
                      acc = 0
                      return
              endif
      endif


    endif 


    do jj = 1,next-4
       if (jj .eq. skipped) then
               index_i(jj,1) = 0
               index_i(jj,2) = 0
               cycle
       endif 
         yy = 1
          do ii=1,2*(next-4)-2
             if (perm(ii) .eq. jj) then
                  index_i(jj,yy) = ii
                  yy=yy+1
             endif
          enddo
       
         if (mod(abs(index_i(jj,1)-index_i(jj,2)),2) .eq. 1) then
            continue
         else
             acc = 0
             return
         endif
     enddo
   
     do jj = 1,next-4
       do ii = jj+1,next-4
          if ((ii .eq. skipped) .or. (jj .eq. skipped)) then
               cycle
          endif
          if ((    ((index_i(ii,1) .lt. index_i(jj,1)) .and. (index_i(ii,2) .lt. index_i(jj,1)))  .or.  &
                  ((index_i(ii,1) .gt. index_i(jj,2)) .and. (index_i(ii,2) .gt. index_i(jj,2)))  .or.  &
                ( ((index_i(ii,1) .lt. index_i(jj,2)) .and. (index_i(ii,2) .lt. index_i(jj,2)))  .and. &
                  ((index_i(ii,1) .gt. index_i(jj,1)) .and. (index_i(ii,2) .gt. index_i(jj,1))))   ) .or. & 
               (    ((index_i(jj,1) .lt. index_i(ii,1)) .and. (index_i(jj,2) .lt. index_i(ii,1)))  .or.  &
                  ((index_i(jj,1) .gt. index_i(ii,2)) .and. (index_i(jj,2) .gt. index_i(ii,2)))  .or.  &
                ( ((index_i(jj,1) .lt. index_i(ii,2)) .and. (index_i(jj,2) .lt. index_i(ii,2)))  .and. &
                  ((index_i(jj,1) .gt. index_i(ii,1)) .and. (index_i(jj,2) .gt. index_i(ii,1))))   ) ) then
              continue
          else
               acc = 0
               return
          endif
       enddo
     enddo

     acc = 1

    else 
       acc = 0

    endif

    endif

  end subroutine check_NLC_2qqbar





  subroutine check_NLC_2qqbar_SF(next,iper,jper,ri,rj,ii,jj,acc)
    implicit none
    integer :: next                 ! number of external particles

    integer,dimension(next-4) :: iper ! order in amplitude
    integer,dimension(next-4) :: jper ! order in conjugate amplitude
    integer :: acc               ! is equal to 0,1,-1 or 99.
                                 ! 99 : LC contributions (NLC coefficient of that term is '-n')
                                 ! 1,-1 ; NLC contribution with positive/negative sign
                                 ! 0 : not a NLC contribution, but NNLC or further suppressed.
    integer :: i1,i2,i3,i4,i5,sign,i,j,yy
    integer ri,rj,ii,jj
    integer,dimension(2*(next-4)) :: temp
    integer,dimension(next-4,2) :: index_i

    integer, dimension(ri+rj) :: temp2
    integer, dimension(2*(next-4)-ri-rj) :: temp3
    integer,dimension(ri) :: itemp
    integer,dimension(next-4-rj) :: itemp2

    integer,dimension(1:ri-1) :: itemp4
    integer,dimension(1:next-4-rj-1) :: itemp5
    integer,dimension(1:rj-1) :: itemp6
    integer,dimension(1:next-4-ri-1) :: itemp7

    logical disjoint
    integer :: ind_i,ind_j
    integer,dimension(2*(next-4)-2) :: perm
    integer skipped

    integer, dimension(ri) :: Aa
    integer, dimension(next-4-ri) :: Bb
    integer, dimension(rj) :: Cc
    integer, dimension(next-4-rj) :: Dd

    if (((ii .eq. 1) .and. (jj .eq. 2)) .or. ((ii .eq. 2) .and. (jj .eq. 1)))  then

        temp(1:ri)= iper(1:ri)
        temp(ri+1:ri+(next-4-rj)) = jper(next-4:rj+1:-1)
        temp(ri+(next-4-rj)+1:2*(next-4)-rj) = iper(ri+1:next-4)
        temp(2*(next-4)-rj+1:2*(next-4)) = jper(rj:1:-1)

        do j = 1,next-4
         yy = 1
          do i=1,2*(next-4)
             if (temp(i) .eq. j) then
                  index_i(j,yy) = i
                  yy=yy+1
             endif
          enddo

         if (mod(abs(index_i(j,1)-index_i(j,2)),2) .eq. 1) then
             continue
         else
             acc = 0
             return
         endif
       enddo

       do j = 1,next-4
       do i = j+1,next-4
          if ( (((  (  ( ((index_i(i,1) .lt. index_i(j,1)) .and. (index_i(i,2) .lt. index_i(j,1))) ) .or.  &
                  (  ( ((index_i(i,1) .gt. index_i(j,2)) .and. (index_i(i,2) .gt. index_i(j,2)))) .or.  &
                  (  ( ((index_i(i,1) .lt. index_i(j,2)) .and. (index_i(i,2) .lt. index_i(j,2)))) .and. &
                    (  ((index_i(i,1) .gt. index_i(j,1)) .and. (index_i(i,2) .gt. index_i(j,1)))   ))))))) .or. &
                    (((  (  ( ((index_i(j,1) .lt. index_i(i,1)) .and. (index_i(j,2) .lt. index_i(i,1))) ) .or.  &
                  (  ( ((index_i(j,1) .gt. index_i(i,2)) .and. (index_i(j,2) .gt. index_i(i,2)))) .or.  &
                  (  ( ((index_i(j,1) .lt. index_i(i,2)) .and. (index_i(j,2) .lt. index_i(i,2)))) .and. &
                    (  ((index_i(j,1) .gt. index_i(i,1)) .and. (index_i(j,2) .gt. index_i(i,1)))   ))))))) ) then
              continue
          else
               acc = 0
               return
          endif
       enddo
       enddo

       acc = -1
    
!************************************************************

  elseif (((ii .eq. 1) .and. (jj .eq. 1)) .or. ((ii .eq. 2) .and. (jj .eq. 2))) then

      if ((ri .eq. rj) .and. ((all(iper .eq. jper)) .or. ((size(iper).eq.0) &
                .and. (size(jper).eq.0)))) then
               acc = 99
               return
      endif

      Aa(1:ri) = iper(1:ri)
      Bb(1:next-4-ri) = iper(ri+1:next-4)
      Cc(1:rj) = jper(1:rj)
      Dd(1:next-4-rj) = jper(rj+1:next-4)

      temp2(1:ri) = Aa(1:ri)
      temp2(ri+1:ri+rj) = Cc(rj:1:-1)

      temp3(1:next-4-ri) = Bb(1:next-4-ri)
      temp3(next-4-ri+1:2*(next-4)-rj-ri) = Dd(next-4-rj:1:-1)

      disjoint = .true.

      do i=1,ri+rj
       if (.not. any(temp3(:) == temp2(i))) then
           disjoint = .true.
       else
           disjoint = .false.
           exit
       endif
      enddo

  if (disjoint) then

        if ((size(Aa) .ne. 0) .and. (all(Bb .eq. Dd) .or. size(Bb) .eq. 0)) then

          do i=1,ri
            if (CC(i).ne.Aa(i)) exit
          enddo
          i1=i
          if (i1 .gt. size(AA)) then
                  acc = 99
                  return
          endif
          do i=i1+1,ri
            if (Cc(i).eq.Aa(i1)) exit
          enddo
          i2=i
          do i=1,ri-i2
             if (Cc(i2+i).ne.Aa(i1+i)) exit
          enddo
          i3=i2+i-1
          i4=i1+i-1
          do i=i1,ri
            if (cc(i).eq.Aa(i4+1)) exit
          enddo
          i5=i
        if (i5.gt.i3) then
            acc = 0
            return
        endif
            sign=1
            itemp(1:i1-1)=Cc(1:i1-1)
            itemp(i1:i4)=Cc(i2:i3)
            itemp(i4+1:i4+i2-i5)=Cc(i5:i2-1)
             if (i1.gt.i5-1) then
               if (i4-i1.eq.0 .and. i2-1-i5.eq.0) then
                   sign=-1
                   continue
               elseif (i4-i1.eq.0 .or. i2-1-i5.eq.0) then
                   acc = 0
                   return
               elseif (i4-i1 + i2-1-i5 .le. next-4-4) then
                   continue
               else
                   sign = 1
             endif
           endif
           itemp(i4+i2-i5+1:i4+i2-i1)=Cc(i1:i5-1)
           itemp(i4+i2-i1+1:ri)=Cc(i3+1:ri)
           if ((all(itemp.eq.Aa))) then                  
             acc=sign
             return
           endif
            acc = 0

        elseif ((size(Bb) .ne. 0) .and. (all(Aa .eq. Cc) .or. size(Aa) .eq. 0)) then

          do i=1,next-4-rj
             if (Dd(i).ne.Bb(i)) exit
          enddo

          i1=i
          if (i1 .gt. size(Bb)) then
                  acc = 99
                  return
          endif
          do i=i1+1,next-4-rj
            if (Dd(i).eq.Bb(i1)) exit
          enddo
          i2=i
          do i=1,next-4-rj-i2
             if (Dd(i2+i).ne.Bb(i1+i)) exit
          enddo
          i3=i2+i-1
          i4=i1+i-1
          do i=i1,next-4-rj
            if (Dd(i).eq.Bb(i4+1)) exit
          enddo
          i5=i
          if (i5.gt.i3) then 
            acc = 0
            return
          endif
            sign=1
            itemp2(1:i1-1)=Dd(1:i1-1)
            itemp2(i1:i4)=Dd(i2:i3)
            itemp2(i4+1:i4+i2-i5)=Dd(i5:i2-1)
             if (i1.gt.i5-1) then
               if (i4-i1.eq.0 .and. i2-1-i5.eq.0) then
                   sign=-1
                   continue
               elseif (i4-i1.eq.0 .or. i2-1-i5.eq.0) then
                   acc = 0
                   return
               elseif (i4-i1 + i2-1-i5 .le. next-4-rj-4) then
                   continue
               else
                 sign = 1
               endif
           endif
         itemp2(i4+i2-i5+1:i4+i2-i1)=Dd(i1:i5-1)
         itemp2(i4+i2-i1+1:next-4-rj)=Dd(i3+1:next-4-rj)
         if ((all(itemp2.eq.Bb))) then
            acc=sign
            return
         endif
            acc = 0

         else
             acc = 0
             return
         endif



     elseif (.not. disjoint) then

      do i=1,ri+rj
       if (any(temp3(:) == temp2(i))) then
           do j=1,size(temp3)
              if (temp3(j) .eq. temp2(i)) then
                        ind_i = i
                        ind_j = j
                        skipped = temp2(i)
                        goto 2
              endif
           enddo
       endif
      enddo

2     if ((size(temp2) .eq. 1) .or. (size(temp3) .eq. 1)) then
            acc = 0
            return
      endif

    perm=0
    !!!   If common generator in A-D pair
     if ((ind_i .le. ri) .and. (ind_j .gt. next-4-ri)) then
      perm(1:ri-ind_i) = temp2(ind_i+1:ri)
      perm(ri-ind_i+1:ri-ind_i+rj) = temp2(ri+1:ri+rj)
      perm(ri+rj-ind_i+1:ri+rj-1) = temp2(1:ind_i-1)

      perm(ri+rj: -1+ 2*(next-4)-ind_j  ) = temp3(ind_j+1:2*(next-4)-ri-rj)
      perm( 2*(next-4)-ind_j   :  -1+ 2*(next-4)-ind_j + next-4-ri   ) =temp3(1:next-4-ri)
      perm(2*(next-4)-ind_j + next-4-ri:2*(next-4)-2) = temp3(next-4-ri+1:ind_j-1)

      itemp4(1:ind_i-1) = temp2(1:ind_i-1)
      itemp4(ind_i:ri-1) = temp2(ind_i+1:ri)

      itemp5(1:2*(next-4)-ri-rj-ind_j) = temp3(2*(next-4)-ri-rj:ind_j+1:-1)
      itemp5(2*(next-4)-ri-rj-ind_j+1:next-4-rj-1) = temp3(ind_j-1:next-4-ri+1:-1)
      if (rj .eq. ri-1) then
              if (all(itemp4 .eq. Cc)) then
                      acc = 0
                      return
              endif
      endif

      if (rj+1 .eq. ri) then
              if (all(itemp5 .eq. Bb)) then
                      acc = 0
                      return
              endif
      endif

      !!!    If common generator in B-C pair
      elseif  ((ind_i .gt. ri) .and. (ind_j .le. next-4-ri)) then
        perm(1:ri+rj-ind_i) = temp2(ind_i+1:ri+rj)
        perm(ri+rj-ind_i+1:ri+rj-ind_i+ri) = temp2(1:ri)
        perm(ri+rj-ind_i+ri+1:ri+rj-1) = temp2(ri+1:ind_i-1)

        perm(ri+rj:rj-1+next-4-ind_j) = temp3(ind_j+1:next-4-ri)
        perm(rj-1+next-4-ind_j+1:-1-ind_j+2*(next-4)) = temp3(next-4-ri+1:2*(next-4)-ri-rj)
        perm(-1-ind_j+2*(next-4)+1: 2*(next-4)-2) = temp3(1:ind_j-1)

        itemp6(1:ri+rj-ind_i) = temp2(ri+rj:ind_i+1:-1)
        itemp6(ri+rj-ind_i+1:rj-1) = temp2(ind_i-1:ri+1:-1)

        itemp7(1:ind_j-1) = temp3(1:ind_j-1)
        itemp7(ind_j -1 +1 : next-4-ri-1) = temp3(ind_j+1:next-4-ri)


      if (rj-1 .eq. ri) then
              if (all(itemp6 .eq. Aa)) then
                      acc = 0
                      return
              endif
      endif

      if (rj .eq. ri+1) then
              if (all(itemp7 .eq. Dd)) then
                      acc = 0
                      return
              endif
      endif

   endif

  do j = 1,next-4
       if (j .eq. skipped) then
               index_i(j,1) = 0
               index_i(j,2) = 0
               cycle
       endif
         yy = 1
          do i=1,2*(next-4)-2
             if (perm(i) .eq. j) then
                  index_i(j,yy) = i
                  yy=yy+1
             endif
          enddo


         if (mod(abs(index_i(j,1)-index_i(j,2)),2) .eq. 1) then
            continue
         else
             acc = 0
             return
         endif

  enddo


     do j = 1,next-4
       do i = j+1,next-4
          if ((i .eq. skipped) .or. (j .eq. skipped)) then
               cycle
          endif
          if ((    ((index_i(i,1) .lt. index_i(j,1)) .and. (index_i(i,2) .lt. index_i(j,1)))  .or.  &
                  ((index_i(i,1) .gt. index_i(j,2)) .and. (index_i(i,2) .gt. index_i(j,2)))  .or.  &
                ( ((index_i(i,1) .lt. index_i(j,2)) .and. (index_i(i,2) .lt. index_i(j,2)))  .and. &
                  ((index_i(i,1) .gt. index_i(j,1)) .and. (index_i(i,2) .gt. index_i(j,1))))   ) .or. &
               (    ((index_i(j,1) .lt. index_i(i,1)) .and. (index_i(j,2) .lt. index_i(i,1)))  .or.  &
                  ((index_i(j,1) .gt. index_i(i,2)) .and. (index_i(j,2) .gt. index_i(i,2)))  .or.  &
                ( ((index_i(j,1) .lt. index_i(i,2)) .and. (index_i(j,2) .lt. index_i(i,2)))  .and. &
                  ((index_i(j,1) .gt. index_i(i,1)) .and. (index_i(j,2) .gt. index_i(i,1))))   ) ) then
              continue
          else
               acc = 0
               return
          endif
       enddo
     enddo

     acc = 1

    endif

    else
       write(*,*) 'ERROR: only for 2qqbar processes'
    endif

   end subroutine check_NLC_2qqbar_SF

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! HELPER FUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ipnext(ia, n)
    ! Compute next permutation starting from the list 'ia' (of length
    ! 'n'). Similar to get_next_iperm (see below), but with ip=n there.
    implicit none
    integer ::n,i,j,itemp
    integer, dimension(n) :: ia
    i=n-1
    do
       if (i.lt.1) exit
       if (ia(i).le.ia(i+1)) exit
       i=i-1
    enddo
    if (i.lt.1) then
       ! go back to the beginning
       do i=1,n
          ia(i)=i
       enddo
       return
    endif
    j = n
    do while (ia(i).gt.ia(j))
       j = j-1
    enddo
    itemp = ia(i)
    ia(i) = ia(j)
    ia(j) = itemp
    i = i+1
    j = n
    do while (i.lt.j)
       itemp = ia(i)
       ia(i) = ia(j)
       ia(j) = itemp
       i = i+1
       j = j-1
    enddo
    return
  end subroutine ipnext
 

  subroutine ipnextgen(ia, n)
    ! Compute next permutation starting from the list 'ia' (of length
    ! 'n'). Similar to get_next_iperm (see below), but with ip=n there.
    implicit none
    integer ::n,i,j,itemp
    integer, dimension(n) :: ia,ia_new
    i=n-1


    do
       if (i.lt.1) exit
       if (ia(i).le.ia(i+1)) exit
       i=i-1
    enddo
    if (i.lt.1) then
       ! go back to the beginning
       do i=1,n
          ia_new(i)=ia(n-i+1)
       enddo
       ia=ia_new
       return
    endif
    j = n
    do while (ia(i).gt.ia(j))
       j = j-1
    enddo
    itemp = ia(i)
    ia(i) = ia(j)
    ia(j) = itemp
    i = i+1
    j = n
    do while (ia(i).gt.ia(j))
       j = j-1
    enddo
    itemp = ia(i)
    ia(i) = ia(j)
    ia(j) = itemp
    i = i+1
    j = n
    do while (i.lt.j)
       itemp = ia(i)
       ia(i) = ia(j)
       ia(j) = itemp
       i = i+1
       j = j-1
    enddo

    return


  end subroutine ipnextgen
  
end module color_algebra
