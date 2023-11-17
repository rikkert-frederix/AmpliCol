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
       & ipnext,&               ! Helper function: get next permutation
       & get_next_iperm,&       ! Helper function: get next permutation (advanced)
       & check_NLC,&
       & check_NLC_1qqbar
contains

  subroutine Tr_allocate(n)
! if input is in the format Tr(...)*Tr(...) with both traces containing a
! permutation of 'n' different matrices, declares the Tr() and coef() arrays
! to cover the maximum possible size in all intermediate steps.
    implicit none
    integer :: n
    call Tr_deallocate
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

!!$  subroutine simplify_cyc()
!!$    implicit none
!!$    integer :: iterm,iprod
!!$    integer,dimension(1) :: min_loc
!!$    integer,allocatable,dimension(:) :: temp
!!$    do iterm=1,Tr(0,0,0)
!!$       do iprod=1,Tr(0,0,iterm)
!!$          ! Use cyclicity to bring to canonical order [always smallest number
!!$          ! first]
!!$          min_loc=minloc(Tr(1:Tr(0,iprod,iterm),iprod,iterm))
!!$          if (min_loc(1).ne.1) then
!!$             allocate(temp(Tr(0,iprod,iterm)))
!!$             temp(1:Tr(0,iprod,iterm))=Tr(1:Tr(0,iprod,iterm),iprod,iterm)
!!$             Tr(1:Tr(0,iprod,iterm)-min_loc(1)+1,iprod,iterm)= &
!!$                                         temp(min_loc(1):Tr(0,iprod,iterm))
!!$             Tr(Tr(0,iprod,iterm)-min_loc(1)+2:Tr(0,iprod,iterm),iprod,iterm)= &
!!$                                         temp(1:min_loc(1)-1)
!!$             deallocate(temp)
!!$          endif
!!$       enddo
!!$    enddo
!!$  end subroutine simplify_cyc

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
  
  subroutine get_next_iperm(ip,ips_in,ips,n)
    ! Given a permutation ips_in, find the next one and return it through ips.
    ! For example for ip=3 (length of permutation list), n=4 (elements to be
    ! considered in the permutation) this gives
    !
    !    ips_in        ips
    !-------------------------
    !    1,2,3   -->   1,2,4
    !    1,2,4   -->   1,3,2
    !    1,3,2   -->   1,3,4
    !    1,3,4   -->   1,4,2
    !    1,4,2   -->   1,4,3
    !    1,4,3   -->   2,1,3
    !    2,1,3   -->   2,1,4
    !    2,1,4   -->   2,3,1
    !    2,3,1   -->   2,3,4
    !    2,3,4   -->   2,4,1
    !    2,4,1   -->   2,4,3
    !    2,4,3   -->   3,1,2
    !    3,1,2   -->   3,1,4
    !    3,1,4   -->   3,2,1
    !    3,2,1   -->   3,2,4
    !    3,2,4   -->   4,1,2
    !    4,1,2   -->   4,1,3
    !    4,1,3   -->   4,2,1
    !    4,2,1   -->   4,2,3
    !    4,2,3   -->   4,3,1
    !    4,3,1   -->   4,3,2
    !    4,3,2   -->   XXXXX
    !
    ! Note that when giving non-sensical inputs (e.g., the last one in the
    ! list above), the code either goes into an infinite loop, or returns some
    ! bogus result. There is no check on the consistency of the input.
    implicit none
    integer :: ip,n,i_up,i,j
    integer,dimension(ip) :: ips,ips_in
    logical :: found
    
    found=.false.
    ips(1:ip)=ips_in(1:ip)
    do i_up=ip,1,-1
       do while (ips(i_up).lt.n) 
          ips(i_up)=ips(i_up)+1
          if (any(ips(1:i_up-1).eq.ips(i_up))) cycle
          found=.true.
          exit
       enddo
       if (found) exit
    enddo
    do i=i_up+1,ip
       do j=1,n
          if (any(ips(1:i).eq.j)) then
             continue
          else
             ips(i)=j
             exit
          endif
       enddo
    enddo
  end subroutine get_next_iperm
  
  
end module color_algebra
