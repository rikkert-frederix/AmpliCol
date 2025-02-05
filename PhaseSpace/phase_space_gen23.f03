module phase_space_gen23_mod
!  use common
  use phase_space_base
  implicit none
  type,extends(phase_space_type),public :: phase_space_gen23
   contains
     procedure :: init => gen23_init
     procedure :: generate_momenta => gen23_generate_momenta
     procedure :: compute_x_from_momenta => gen23_compute_x_from_momenta
  end type phase_space_gen23
  private
  logical :: includePDF
  ! TECHNIAL PARAMETERS
  ! vebose:
  logical,parameter :: verbose=.true.,eff_photon_int=.true.
  logical,parameter,public :: debug=.false.
  ! importance sampling (0d0=flat transformation; -1d0=1/x transformation):
  real(kind=8),parameter :: ip=-1d0,ip_shat=-1.2d0
  ! tiny parameter cutoff to prevent/reduce numerical instabilities:
  real(kind=8),parameter :: vtiny=1d-12,tiny=1d-8
  real(kind=8),parameter :: pi=3.1415926535897932d0
  logical,parameter :: use_t_channel_at_start=.true.

contains
  subroutine gen23_init(this,sqrts,n,m,o,s_cut,pt_cut,rap_cut,dr_cut,sqrt_s_min,t_chan,include_pdf,part)
    ! Phase-space initialisation routines.
    implicit none
    class(phase_space_gen23),intent(inout) :: this
    ! INPUT
    ! Sqrt(s-hat), i.e, the collision energy
    real(kind=8),intent(in) :: sqrts
    ! number of particles (initial state + final state)
    integer(kind=4),intent(in) :: n
    ! the colour order:
    integer(kind=4),dimension(n),intent(in) :: o,part
    integer(kind=4),dimension(n) :: ord_temp
    ! cut on the minimum invariant mass of all pairs of particles,
    ! abs((p_i+p_j)^2)>s_cut, (initial and final state). s_cut(1) is between
    ! an initial and a final state particle; s_cut(2) is between two final
    ! state particles.
    real(kind=8),intent(in) :: s_cut(2),pt_cut,dr_cut,rap_cut,sqrt_s_min
    ! masses of all the particles. The two incoming particles must be
    ! massless.
    real(kind=8),dimension(n),intent(in) :: m
    ! ** With t_chan=.false., use the Byckling Kayantie 2->3 phase-space
    ! generation, i.e., eq (8) of E.~Byckling and K.~Kajantie,
    ! ``Reductions of the phase-space integral in terms of simpler
    ! processes,'' Phys. Rev. 187 (1969), 2008-2016,
    ! doi:10.1103/PhysRev.187.2008.
    ! ** With t_chan=.true., use equation (13) instead.
    ! ** Note, for the first splitting, there is an additional logical
    ! parameter in the generate_momenta() subtroutine that allows one
    ! to choose between s-channel and t-channel.
    logical,intent(in) :: t_chan
    ! Should we include a PDF set? Currently, only the NNPDF2.3 NLO QED is available.
    logical,intent(in) :: include_pdf
    integer(kind=4) :: i,j,ns,nns
    real(kind=8) :: ptcut
    this%sqrtshat=sqrts
    this%sqrts=sqrts
    this%t_channel=t_chan
    if (verbose) then
       write (*,*) 'Setting up',n,'particle phase-space'
       write (*,*) 'Total available energy, sqrt(s-hat) =',this%sqrtshat
       write (*,*) 'Cut on invariants used in the phase-space generation: abs((p_i+p_j)^2) >=',s_cut
       write (*,*) 'Use the simple t-channel?',this%t_channel
    endif
    includePDF=include_pdf
    call gen23_deallocate
    this%next=n
    this%ndim=3*(this%next-2)-4
    if (includePDF) this%ndim=this%ndim+2 ! the two Bjorken x's
    allocate(this%order(this%next))
    allocate(this%invm(maskr(this%next)))
    allocate(this%invm_min(maskr(this%next)))
    allocate(this%ETmin(maskr(this%next)))
    allocate(this%invm_max(maskr(this%next)))
    allocate(this%pp(0:3,0:maskr(this%next)))
    this%pp(0:3,0:maskr(this%next))=0d0
    allocate(this%p(0:3,this%next))
    allocate(this%x(this%ndim))
    allocate(this%sets(0:this%next-2,3))
    ! masses of external particles
    do i=1,n
       if ((i.eq.1 .or. i.eq.2) .and. m(i).ne.0d0) then
          write (*,*) 'ERROR in gen23_init() -- ', &
               & 'incoming particles should be massless'
          write (*,*) m
          stop 1
       endif
       this%invm(ibset(0,i-1))=m(i)**2
       this%invm(ibclr(maskr(this%next),i-1))=m(i)**2
    enddo
    if (verbose) write (*,*) 'masses:',m(1:n)
    if (pt_cut.gt.0d0) then
       this%drcut=dr_cut
       ptcut=pt_cut
    else
       this%drcut=0d0
       ptcut=0d0
    endif
    call setup_PS_cuts(s_cut)
    ! Bring the colour order to a canonical order (first in the list
    ! should be particle 1, i.e., the first incoming particle).
    do i=1,this%next
       if (o(i).eq.1) then
          do j=0,this%next-1
             this%order(j+1)=o(1+mod(i+j-1,this%next))
          enddo
          exit
       endif
    enddo

    ns=0
    nns=1
    this%sets=0
    ! move all singlets to the end of the order
    if (eff_photon_int) then
    ord_temp=this%order
    do i=1,this%next
       if (part(this%order(i)).ne.21.and.abs(part(this%order(i))).gt.6.and.&
           this%order(i).gt.2) then
            ord_temp(this%next-ns)=this%order(i)
            ns=ns+1
       else
            ord_temp(nns)=this%order(i)
            nns=nns+1
       endif
    enddo


    !this%sets(1:ns,3)=ord_temp(this%next-ns+1:this%next)
    !do j=this%next-ns+1,this%next
    !    this%sets(0,3)=ibset(this%sets(0,3),ord_temp(j)-1)
    !enddo
    j=1
    do i=1,this%next
       if (part(this%order(i)).eq.21.or.abs(part(this%order(i))).le.6.) then
            ord_temp(j)=this%order(i)
            j=j+1
       endif
    enddo
    this%order=ord_temp
    endif

   
    ns=0 ! needed for the new ordering of the phase-space order 
    i=0
    do i=2,this%next-ns
       if (this%order(i).eq.2) then
          do j=i+1,this%next-ns
             this%sets(0,2)=ibset(this%sets(0,2),this%order(j)-1)
          enddo
          this%sets(1:i-2,1)=this%order(2:i-1)
          this%sets(1:this%next-ns-i,2)=this%order(i+1:this%next-ns)
          exit
       endif
       this%sets(0,1)=ibset(this%sets(0,1),this%order(i)-1)
    enddo
    if (verbose) then
       write (*,*) "set 1:",this%sets(:,1)
       write (*,*) "set 2:",this%sets(:,2)
       write (*,*) "set 3:",this%sets(:,3)
    endif
    if (verbose) then
       write (*,*) "Power in importance sampling:",ip
    endif
  contains
    subroutine setup_PS_cuts(s_cut)
      ! Given s_cut = abs((p_i+p_j)^2), fills the minimum (s-channel)
      ! and/or maximum (t-channel) values the invariants can be in the
      ! phase-space generation. Does not apply these cuts on invariants
      ! not used in the phase-space generation.
      implicit none
      real(kind=8),intent(in) :: s_cut(2)
      real(kind=8) :: mass
      integer(kind=4) :: i,j,npart
      this%invm_min=0d0
      this%invm_max=0d0
      do i=1,maskr(this%next)
         npart=popcnt(i)
         if (btest(i,0).and.btest(i,1)) then
            this%invm_min(i)=0d0
         elseif (btest(i,0).or.btest(i,1)) then
            if (npart.eq.2 .or. npart.eq.this%next-2) then
               this%invm_max(i)=-s_cut(1)
            endif
         else
            mass=0d0
            do j=0,this%next-1
               if (btest(i,j)) then
                  mass=mass+sqrt(this%invm(ibset(0,j)))
               endif
            enddo
            if (npart.eq.this%next-2) then
               this%invm_min(i)=max(s_cut(2)*npart**2,mass**2)
            else
               this%invm_min(i)=max(s_cut(2)*(npart)*(npart-1)/2d0,mass**2)
            endif
         endif
      enddo
      call setup_ETmin
    end subroutine setup_PS_cuts

    subroutine setup_ETmin
      ! Setup the minimum required (transverse) energy for each (combination of)
      ! final state particle(s) in the collision c.o.m. frame. Based on the
      ! masses and the ptcut (i.e., assumes that all pz=0 and pT=pTcut)
      integer :: i,j
      this%ETmin(1:maskr(this%next))=0d0
      do i=1,maskr(this%next)
         if (btest(i,0).or.btest(i,1)) cycle ! skip the ones that include incoming particles
         do j=0,this%next-1
            if (btest(i,j)) this%ETmin(i)=this%ETmin(i)+sqrt(this%invm(ibset(0,j))+ptcut**2)
         enddo
         this%ETmin(i)=max(this%ETmin(i),sqrt(this%invm_min(i)))
      enddo
    end subroutine setup_ETmin

    subroutine gen23_deallocate
      implicit none
      if (allocated(this%order)) deallocate(this%order)
      if (allocated(this%invm)) deallocate(this%invm)
      if (allocated(this%invm_min)) deallocate(this%invm_min)
      if (allocated(this%invm_max)) deallocate(this%invm_max)
      if (allocated(this%ETmin)) deallocate(this%ETmin)
      if (allocated(this%pp)) deallocate(this%pp)
      if (allocated(this%p)) deallocate(this%p)
      if (allocated(this%x)) deallocate(this%x)
      if (allocated(this%sets)) deallocate(this%sets)
    end subroutine gen23_deallocate
  end subroutine gen23_init

  subroutine gen23_generate_momenta(this,xx)
    ! Wrapper for the routine that generates the momenta.
    implicit none
    class(phase_space_gen23),intent(inout) :: this
    real(kind=8),dimension(99),intent(in) :: xx
    integer(kind=4) :: i,ix
    real(kind=8) :: ycm
    this%x(1:this%ndim)=xx(1:this%ndim)
    this%jac=1d0
    ix=0
    if (includePDF) call generate_initial_state
    call generate_momenta
    do i=1,this%next
       if (includePDF) then
          ! Note: 'ycm' is the rapidity needed to go from lab to CM
          ! frame. Hence, here we boost from CM to lab frame with '-ycm'
          call boostz(this%pp(0:3,ibset(0,i-1)),-ycm,this%p(0:3,i))
       else
          this%p(0:3,i)=this%pp(0:3,ibset(0,i-1))
       endif
    enddo
  contains
    subroutine generate_initial_state
      implicit none
      real(kind=8) :: tau
      call generate_tau(tau)
      call generate_y(tau)
      this%sqrtshat=sqrt(tau)*this%sqrts
      this%xbjrk(1)=sqrt(tau)*exp(ycm)
      this%xbjrk(2)=sqrt(tau)*exp(-ycm)
      if (debug) write (*,*) 'sqrtshat :',this%sqrtshat,this%xbjrk(1:2),this%sqrtshat**2
    end subroutine generate_initial_state

    subroutine generate_tau(tau)
      implicit none
      real(kind=8),intent(out) :: tau
      real(kind=8) :: smin,smax,shat
      smin=max(this%invm_min(maskr(this%next)-3),this%ETmin(maskr(this%next)-3)**2)
      smax=this%sqrts**2
      ix=ix+1
      call random_to_var(this%x(ix),ip_shat,smin,smax,shat,this%jac)
      tau=shat/smax
      this%jac=this%jac/smax
    end subroutine generate_tau

    subroutine generate_y(tau)
      implicit none
      real(kind=8),intent(in) :: tau
      real(kind=8) ::  ymin,ymax
      ymin= log(tau)/2d0
      ymax=-log(tau)/2d0
      ix=ix+1
      call random_to_var(this%x(ix),0d0,ymin,ymax,ycm,this%jac)
    end subroutine generate_y

    subroutine generate_momenta
      implicit none
      integer(kind=4) :: i,j,inext,im1
      integer(kind=4),dimension(3) :: set
      ! incoming momenta
      this%pp(0,ibset(0,0))=this%sqrtshat/2d0
      this%pp(0,ibset(0,1))=this%sqrtshat/2d0
      this%pp(1:2,ibset(0,0))=0d0
      this%pp(1:2,ibset(0,1))=0d0
      this%pp(3,ibset(0,0))= this%pp(0,ibset(0,0))
      this%pp(3,ibset(0,1))=-this%pp(0,ibset(0,1))
      ! incoming momenta when labeled from the other side
      this%pp(0:3,maskr(this%next)-ibset(0,0))=this%pp(0:3,1)
      this%pp(0:3,maskr(this%next)-ibset(0,1))=this%pp(0:3,2)
      ! invariant mass of all final state particles combined
      this%invm(ibset(0,0)+ibset(0,1))=this%sqrtshat**2
      this%invm(maskr(this%next)-ibset(0,0)-ibset(0,1))=this%sqrtshat**2
      set(1)=this%sets(0,1)
      set(2)=this%sets(0,2)
      set(3)=this%sets(0,3)
      do i=1,popcnt(set(3))
         inext=ibset(0,this%sets(i,3)-1)
         set(3)=set(3)-inext
         if (popcnt(set(1)+set(2)+set(3)).eq.1) then
            call gent_one_step(inext,set(1)+set(2)+set(3),1)
            if (this%jac.le.0d0) return
            this%pp(0:3,set(1)+set(2)+set(3)+2)=this%pp(0:3,1)-this%pp(0:3,inext)
            this%invm(set(1)+set(2)+set(3)+2)=dot(this%pp(0:3,set(1)+set(2)+set(3)+2),this%pp(0:3,set(1)+set(2)+set(3)+2))
            !do j=1,this%next
            !   write(*,*) j,this%pp(0:3,2**(j-1))
            !   write(*,*) dot(this%pp(0:3,2**(j-1)),this%pp(0:3,2**(j-1)))
            !enddo
         else
            call double_t(inext,set(1)+set(2)+set(3),1,2)
            if (this%jac.le.0d0) return
            this%pp(0:3,set(1)+set(2)+set(3)+1)=this%pp(0:3,2)-this%pp(0:3,inext)
            this%invm(set(1)+set(2)+set(3)+1)=dot(this%pp(0:3,set(1)+set(2)+set(3)+1),this%pp(0:3,set(1)+set(2)+set(3)+1))
            !do j=1,this%next
            !   write(*,*) j,this%pp(0:3,2**(j-1))
            !   write(*,*) dot(this%pp(0:3,2**(j-1)),this%pp(0:3,2**(j-1)))
            !enddo
         endif
      enddo
      ! Generate the central 2->2 process in case both set(1) and set(2) are not empty
      if (popcnt(set(1)).gt.1 .and. popcnt(set(2)).gt.1) then
         if (debug) write (*,*) 'two sets with at least two ',&
              & 'particles',popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         if (use_t_channel_at_start) then
            call gent_one_step(set(2),set(1),1)
         else
            call gens_one_step(set(2),set(1))
         endif
         if (this%jac.le.0d0) return
         this%pp(0:3,set(2)+2)=this%pp(0:3,1)-this%pp(0:3,set(1))
         this%invm(set(2)+2)=dot(this%pp(0:3,set(2)+2),this%pp(0:3,set(2)+2))
      elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).gt.1) then
         if (debug) write (*,*) 'special double t-channel (1)'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call double_t(set(1),set(2),1,2)
         if (this%jac.le.0d0) return
         this%pp(0:3,set(2)+2)=this%pp(0:3,1)-this%pp(0:3,set(1))
         this%invm(set(2)+2)=dot(this%pp(0:3,set(2)+2),this%pp(0:3,set(2)+2))
      elseif (popcnt(set(1)).gt.1 .and. popcnt(set(2)).eq.1) then
         if (debug) write (*,*) 'special double t-channel (2)'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call double_t(set(2),set(1),1,2)
         if (this%jac.le.0d0) return
         this%pp(0:3,set(1)+1)=this%pp(0:3,2)-this%pp(0:3,set(2))
         this%invm(set(1)+1)=dot(this%pp(0:3,set(1)+1),this%pp(0:3,set(1)+1))
      elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).eq.1) then
         if (debug) write (*,*) '2->2 scattering with one particle in each set'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call gent_one_step(set(2),set(1),1)
         if (this%jac.le.0d0) return
         this%pp(0:3,set(2)+2)=this%pp(0:3,1)-this%pp(0:3,set(1))
         this%invm(set(2)+2)=dot(this%pp(0:3,set(2)+2),this%pp(0:3,set(2)+2))
      endif
      do i=1,2
         if (popcnt(set(i)).le.1) cycle ! at least 2 particles in a set
         inext=ibset(0,this%sets(1,i)-1)
         set(i)=set(i)-inext
         if (popcnt(set(i)).ge.2) then
            ! at least 3 particles in a set
            if (debug) write (*,*) 'At least 3 particles in a set',&
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            call gent_one_step(set(i),inext,i)
            if (this%jac.le.0d0) return
            this%pp(0:3,(3-i)+set(i)+inext)=this%pp(0:3,i)-this%pp(0:3,this%sets(0,(3-i)))
            this%invm((3-i)+set(i)+inext)=dot(this%pp(0:3,(3-i)+set(i)+inext),this%pp(0:3,(3-i)+set(i)+inext))
            this%pp(0:3,(3-i)+set(i))=this%pp(0:3,i)-this%pp(0:3,this%sets(0,(3-i)))-this%pp(0:3,inext)
            this%invm((3-i)+set(i))=dot(this%pp(0:3,(3-i)+set(i)),this%pp(0:3,(3-i)+set(i)))
            do j=2,popcnt(set(i))-1
               ! loop over the remaining particles in the set
               inext=ibset(0,this%sets(j,i)-1)
               im1=ibset(0,this%sets(j-1,i)-1)
               set(i)=set(i)-inext
               if (this%t_channel) then
                  call gent_one_step(inext,set(i),3-i)
               else
                  call gen23_one_step(inext,set(i),3-i,im1)
               endif
               if (this%jac.le.0d0) return
            enddo
            inext=ibset(0,this%sets(j,i)-1)
            im1=ibset(0,this%sets(j-1,i)-1)
            set(i)=set(i)-inext
            if (this%t_channel) then
               call gent_one_step(inext,set(i),3-i)
            else
               call gen23_one_step(inext,set(i),3-i,im1)
            endif
            if (this%jac.le.0d0) return
         elseif (popcnt(set(i)).eq.1 .and. popcnt(this%sets(0,3-i)).ne.0) then
            ! Exactly 2 particles in a set (and the other set contains at least one)
            if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
                 & 'the other set contains at least one)', &
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            im1=3-i
            this%pp(0:3,set(i)+inext+im1)=this%pp(0:3,set(i)+inext)-this%pp(0:3,im1)
            this%pp(0:3,set(i)+inext+i+im1)=this%pp(0:3,set(i)+inext+im1)-this%pp(0:3,i)
            this%invm(set(i)+inext+im1)=dot(this%pp(0:3,set(i)+inext+im1),this%pp(0:3,set(i)+inext+im1))
            this%invm(set(i)+inext+i+im1)=dot(this%pp(0:3,set(i)+inext+im1+i),this%pp(0:3,set(i)+inext+im1+i))
            if (this%t_channel) then
               call gent_one_step(set(i),inext,i)
            else
               call gen23_one_step(set(i),inext,i,im1)
            endif
            if (this%jac.le.0d0) return
         elseif (popcnt(set(i)).eq.1 .and. popcnt(this%sets(0,3-i)).eq.0) then
            ! Exactly 2 particles in a set (and the other set contains none)
            if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
                 & 'the other set contains none)', &
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            call gent_one_step(set(i),inext,i)
            if (this%jac.le.0d0) return
         else
            write (*,*) 'Inconsistent sets'
            write (*,*) i,':',this%sets(:,i)
            write (*,*) 3-i,':',this%sets(:,3-i)
            stop 
         endif
         ! We need to get the momentum of the final particle of the set.
         this%pp(0:3,set(i))=this%pp(0:3,set(i)+inext+(3-i))+this%pp(0:3,(3-i))-this%pp(0:3,inext)
      enddo
      if (debug) call test_momenta
      ! Add factors of 2*pi
      this%jac=this%jac/((2d0*pi)**(3*(this%next-2)-4))
      ! Add flux factor
      this%jac=this%jac/(2d0*this%sqrtshat**2)
    end subroutine generate_momenta

    subroutine test_momenta
      ! Writes the momenta to the screen.
      implicit none
      integer(kind=4) :: i,i1,i2,i1b,i2b
      real(kind=8),dimension(0:3) :: ptot
      write (*,*) 'Momenta check:'
      if (ix.ne.this%ndim) then
         write (*,*) 'ERROR: number of random numbers used not consistent',ix,this%ndim
         stop 1
      endif
      ptot(0:3)=0d0
      do i=0,this%next-1
         ptot(0:3)=ptot(0:3)+this%pp(0:3,ibset(0,i))
         write (*,*) i+1,this%pp(0:3,ibset(0,i)),dot(this%pp(0,ibset(0,i)),this%pp(0,ibset(0,i)))
      enddo
      write (*,*) 'ptot',ptot(0:3)
      do i=1,this%next
         i1=this%order(i)
         i2=this%order(mod(i,this%next)+1)
         i1b=ibset(0,i1-1)
         i2b=ibset(0,i2-1)
         write (*,*) '(pi+pj)^2',i1,i2,this%invm(i1b+i2b),this%invm(maskr(this%next)-(i1b+i2b))&
              &,dot(this%pp(0:3,i1b)+this%pp(0:3,i2b),this%pp(0:3,i1b)+this%pp(0:3,i2b))
      enddo
      write (*,*) 'number of dimensions',this%ndim
      write (*,*) ''
      write (*,*) ''
    end subroutine test_momenta


    subroutine double_t(i,ir,ia,ib)
      ! Generates the p_a+p_b->p_i+p_ir process, with (p_a+p_i)^2 and
      ! (p_b+p_i)^2 (and phi) as integration variables. Note that this means
      ! that the mass of the system ir is derived from the integration
      ! variables.
      !
      ! p_a and p_b are assumed to be massless, and p_a and p_b
      ! back-to-back incoming particles. The mass of i is fixed.
      implicit none
      integer(kind=4),intent(in) :: i,ir,ia,ib
      real(kind=8) :: tmin,tmax,phi,pt2,yr,Eimax,pzmax
      if (popcnt(i).ne.1 .or. popcnt(ir).le.1) then
         write (*,*) 'Subroutine only for i is a single particle '&
              //'and ir is more than 1',i,ir,popcnt(i),popcnt(ir)
         stop 1
      endif
      yr=sqrt(lambda(this%invm(ia+ib),this%invm(i),this%invm_min(ir)))
      tmin=(-this%invm(ia+ib)+this%invm(i)+this%invm_min(ir)-yr)/2d0
      tmax=(-this%invm(ia+ib)+this%invm(i)+this%invm_min(ir)+yr)/2d0
      if (this%invm_max(ir+ib).ne.0d0) tmax=min(tmax,this%invm_max(ir+ib))
      if (this%invm_min(ir+ib).ne.0d0) tmin=max(tmin,this%invm_min(ir+ib))
      ! Additional constraints on tmin and tmax due to pp(0,i) and pp(0,ir)
      ! being larger than ETmin(i) and ETmin(ir), respectively:
      pzmax=sqrt(lambda(this%sqrtshat**2,this%ETmin(i)**2,this%ETmin(ir)**2))/(2d0*this%sqrtshat)
      Eimax=this%sqrtshat-sqrt(this%ETmin(ir)**2+pzmax**2)
      tmin=max(tmin,this%invm(i)-this%sqrtshat*(Eimax+pzmax))
      tmax=min(tmax,this%invm(i)-this%sqrtshat*(Eimax-pzmax))
      if (tmin.ge.tmax) then
         this%jac=-1d0
         if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
         return
      endif
      ix=ix+1
      call random_to_var(this%x(ix),ip,tmin,tmax,this%invm(i+ia),this%jac)
      if (debug) then
         write (*,*) 'dt- i+ia',i+ia,this%invm(i+ia),tmin,tmax
      endif
      tmin=-this%invm(ia+ib)-this%invm(i+ia)+this%invm(i)+this%invm_min(ir)
      tmax=this%invm(i)*(this%invm(i)-this%invm(ia+ib)-this%invm(i+ia))/(this%invm(i)-this%invm(i+ia))
      if (this%invm_max(ir+ib).ne.0d0) tmax=min(tmax,this%invm_max(ir+ib))
      if (this%invm_min(ir+ib).ne.0d0) tmin=max(tmin,this%invm_min(ir+ib))
      ! Additional constraints on tmin and tmax due to pp(0,i) and pp(0,ir)
      ! being larger than ETmin(i) and ETmin(ir), respectively:
      tmin=max(tmin,this%invm(i)-this%sqrtshat**2*(1-this%ETmin(ir)**2/(this%sqrtshat**2+this%invm(i+ia)-this%invm(i))))
      tmax=min(tmax,this%invm(i)-this%sqrtshat**2*(this%ETmin(i)**2/(this%invm(i)-this%invm(i+ia))))
      if (tmin.ge.tmax) then
         this%jac=-2d0
         if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
         return
      endif
      ix=ix+1
      call random_to_var(this%x(ix),ip,tmin,tmax,this%invm(i+ib),this%jac)
      if (debug) then
         write (*,*) 'dt- i+ib',i+ib,this%invm(i+ib),tmin,tmax
      endif
      ix=ix+1
      call random_to_var(this%x(ix),0d0,0d0,2d0*pi,phi,this%jac)
      if (debug) then
         write (*,*) 'dt- phi',phi
      endif
      pt2=this%invm(i+ia)*this%invm(i+ib)/this%invm(ia+ib)+ &
           & this%invm(i)**2/this%invm(ia+ib)-(this%invm(i+ia)+this%invm(i+ib))*this%invm(i)/this%invm(ia+ib)-this%invm(i)
      if (pt2/this%invm(ia+ib).lt.-vtiny) then
         write (*,*) "ERROR in double_t: pt2 should be larger than zero", &
              & pt2,this%invm(i+ia),this%invm(i+ib),this%invm(ia+ib),this%invm(i)
         write (*,*) this%invm(i+ia)*this%invm(i+ib)/this%invm(ia+ib), &
              & this%invm(i)**2/this%invm(ia+ib),(this%invm(i+ia)+this%invm(i+ib))*this%invm(i)/this%invm(ia+ib),this%invm(i)
         stop 1
      elseif (pt2.lt.0d0) then
         pt2=0d0
      endif
      this%pp(0,i)=(-this%invm(i+ia)-this%invm(i+ib)+2d0*this%invm(i))/(2d0*this%sqrtshat)
      this%pp(1,i)=sqrt(pt2)*cos(phi)
      this%pp(2,i)=sqrt(pt2)*sin(phi)
      this%pp(3,i)=(this%invm(i+ia)-this%invm(i+ib))/(2d0*this%sqrtshat)
      this%pp(0,ir)=this%sqrtshat-this%pp(0,i)
      this%pp(1:3,ir)=-this%pp(1:3,i)
      this%invm(ir)=dot(this%pp(0,ir),this%pp(0,ir))
      if (this%invm(ir).le.0d0) then
         write (*,*) "ERROR in double_t: invariant mass of system", &
              & " must be larger than zero",ir,this%invm(ir),i
         write (*,*) this%invm(ir),this%invm(ia+ib)+this%invm(i+ia)+this%invm(i+ib)-this%invm(i)&
              &,this%invm(ia+ib),this%invm(i +ia),this%invm(i+ib),this%invm(i)
         stop
      endif
      this%jac = this%jac/(4d0*sqrt(lambda(this%invm(ir+i),0d0,0d0)))
    end subroutine double_t

    subroutine gen23_one_step(i,ir,ib,im1)
      ! Generates one step using the 2->3 setup, using the invariants
      ! shat(i), s(i) and t(i) as defined in E.~Byckling and K.~Kajantie,
      ! ``Reductions of the phase-space integral in terms of simpler
      ! processes,'' Phys. Rev. 187 (1969), 2008-2016,
      ! doi:10.1103/PhysRev.187.2008.  Assumes massless incoming particles.
      implicit none
      integer(kind=4),intent(in) :: im1,i,ir,ib
      real(kind=8) :: tmin,tmax,smin,smax,phi1,phi2,gram4,V,sqrtGG,shatmin,shatmax,y,base,root,phi_rot,&
           etminir,etmini
      real(kind=8),dimension(0:3) :: pi1,pr1,ppibir1,pi2,pr2,ppibir2,piir,pib,pim1,piirr,pim1r
      real(kind=8),external :: ran2
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) this%invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax)
         call generate_mass(i,shatmin,shatmax)
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax)
         call generate_mass(ir,shatmin,shatmax)
      endif
      if (this%jac.le.0d0) return
      if (debug) then
         write (*,*) '23- i    ',i,this%invm(i)
         write (*,*) '23- ir   ',ir,this%invm(ir)
      endif
      call tminmax(this%invm(ir+i),this%invm(ir+i+ib),this%invm(ir),this%invm(i),0d0,tmin,tmax)
      if (this%invm_max(ir+ib).ne.0d0) tmax=min(tmax,this%invm_max(ir+ib))
      if (this%invm_min(ir+ib).ne.0d0) tmin=max(tmin,this%invm_min(ir+ib))
      ! Make sure that the t-range is compatible with the pT cut. Since t is an
      ! invariant we can compute it in any frame. Let's use the frame in which
      ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
      ! ETmin() is boost invariant in the z-direction)
      this%pp(0:3,i+ir)=this%pp(0:3,i+ir+ib)+this%pp(0:3,ib)
      y=log((this%pp(0,i+ir)+this%pp(3,i+ir))/(this%pp(0,i+ir)-this%pp(3,i+ir)))/2d0
      call boostz(this%pp(0,i+ir),y,piir)
      call boostz(this%pp(0,ib),y,pib)
      if ( piir(1)**2+piir(2)**2.lt.this%ETmin(i)**2-this%invm(i) .and. popcnt(i).eq.1 ) then
         etminir=max(this%ETmin(ir),sqrt(this%invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(this%ETmin(i)**2-this%invm(i)))**2) )
      else
         etminir=max(this%ETmin(ir),sqrt(this%invm(ir)))
      endif
      if ( piir(1)**2+piir(2)**2.lt.this%ETmin(ir)**2-this%invm(ir) .and. popcnt(ir).eq.1 ) then
         etmini=max(this%ETmin(i),sqrt(this%invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(this%ETmin(ir)**2-this%invm(ir)))**2) )
      else
         etmini=max(this%ETmin(i),sqrt(this%invm(i)))
      endif
      base=piir(0)**2-ETmini**2+ETminir**2
      ! Note, root=lambda(piir(0)**2,this%ETmin(i)**2,this%ETmin(ir)**2), but the
      ! following is more stable:
      root=(piir(0)-ETmini-ETminir)*(piir(0)+ETmini-ETminir)*&
           (piir(0)-ETmini+ETminir)*(piir(0)+ETmini+ETminir)
      if (root.lt.0d0) then
         this%jac=-33d0
         if (debug) write (*,*) 'root.lt.0d0',root
         return
      endif
      tmin=max(tmin,this%invm(ir)-pib(0)/piir(0)*(base+sqrt(root)))
      tmax=min(tmax,this%invm(ir)-pib(0)/piir(0)*(base-sqrt(root)))
      if (tmin.ge.tmax) then
         this%jac=-3d0
         if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
         return
      endif
      ix=ix+1
      call random_to_var(this%x(ix),ip,tmin,tmax,this%invm(ir+ib),this%jac)
      if (debug) then
         write (*,*) '23- ir+ib',ir+ib,this%invm(ir+ib),tmin,tmax
      endif
      call sminmax(this%invm(ir+i),this%invm(ir),this%invm(ir+i+im1),this%invm(ir+i+ib)&
           &,this%invm(ir+ib),this%invm(ir+ib+i+im1),this%invm(i),this%invm(im1),smin,smax,V,sqrtGG)
      if (this%invm_min(i+im1).ne.0d0) smin=max(smin,this%invm_min(i+im1))
      if (this%invm_max(i+im1).ne.0d0) smax=min(smax,this%invm_max(i+im1))
      if (im1.gt.2) then
         ! Boost and rotate in z-direction such that pp(:,im1) goes in the x-direction.
         y=log((this%pp(0,im1)+this%pp(3,im1))/(this%pp(0,im1)-this%pp(3,im1)))/2d0
         call boostz(this%pp(0,i+ir),y,piirr)
         call boostz(this%pp(0,im1),y,pim1r)
         call boostz(this%pp(0,ib),y,pib)
         phi_rot=atan(this%pp(2,im1)/this%pp(1,im1))
         if(this%pp(1,im1).lt.0d0) phi_rot=phi_rot+pi
         call rotz(piirr,-phi_rot,piir)
         call rotz(pim1r,-phi_rot,pim1)
         ! Eir > Etmin(ir) + constraint coming from t
         etminir=max(pib(0)*this%ETmin(ir)**2/(this%invm(ir)-this%invm(ir+ib))+(this%invm(ir)-this%invm(ir+ib))/(4d0*pib(0)),&
              this%ETmin(ir))
         smax=min(smax,&
              this%invm(i)+this%invm(im1)+2d0*(piir(0)-etminir)*pim1(0)+2d0*sqrt((piir(0)-etminir)**2-this%invm(i))*pim1(1))

         if(this%invm(i).eq.0d0) then
            smin=max(smin,2d0*this%ETmin(i)*(pim1(0)-pim1(1)*cos(this%drcut)))
         endif

      endif
      if (smin.ge.smax) then
         this%jac=-4d0
         if (debug) write (*,*) 'smin.ge.smax',smin,smax
         return
      endif
      ix=ix+1
      call random_to_var(this%x(ix),ip,smin,smax,this%invm(i+im1),this%jac)
      if (debug) then
         write (*,*) '23- i+im1',i+im1,this%invm(i+im1),smin,smax
      endif
      ! Generate the momenta from the integration variables. Since there is an
      ! ambiguity in phi, get both of them and pick the one that passes the cuts
      ! (if it's only one). If both pass, simply pick one of the two at random
      ! with a flat prior.
      phi1=getphifroms(this%invm(i+im1),this%invm(ir+i),this%invm(ir),this%invm(ir+i+im1)&
           &,this%invm(ir+i+ib),V,sqrtGG,1d0)
      call gentcms2(this%pp(0,ib),this%pp(0,ib+ir+i),this%pp(0,ib+ir+i+im1),this%invm(ir+ib),phi1 &
           &,sqrt(this%invm(i)),sqrt(this%invm(ir)),pi1,ppibir1)
      pr1(0:3)=this%pp(0:3,ir+i)-pi1(0:3)
      phi2=getphifroms(this%invm(i+im1),this%invm(ir+i),this%invm(ir),this%invm(ir+i+im1)&
           &,this%invm(ir+i+ib),V,sqrtGG,0d0)
      call gentcms2(this%pp(0,ib),this%pp(0,ib+ir+i),this%pp(0,ib+ir+i+im1),this%invm(ir+ib),phi2 &
           &,sqrt(this%invm(i)),sqrt(this%invm(ir)),pi2,ppibir2)
      pr2(0:3)=this%pp(0:3,ir+i)-pi2(0:3)
      if ( pi1(0)**2-pi1(3)**2.ge.this%ETmin(i)**2 .and. pr1(0)**2-pr1(3)**2.ge.this%ETmin(ir)**2 .and. &
           pi2(0)**2-pi2(3)**2.ge.this%ETmin(i)**2 .and. pr2(0)**2-pr2(3)**2.ge.this%ETmin(ir)**2 ) then
         if(ran2().gt.0.5d0) then
            this%pp(0:3,i)=pi1(0:3)
            this%pp(0:3,ir)=pr1(0:3)
            this%pp(0:3,ib+ir)=ppibir1(0:3)
         else
            this%pp(0:3,i)=pi2(0:3)
            this%pp(0:3,ir)=pr2(0:3)
            this%pp(0:3,ib+ir)=ppibir2(0:3)
         endif
      elseif (pi1(0)**2-pi1(3)**2.ge.this%ETmin(i)**2 .and. pr1(0)**2-pr1(3)**2.ge.this%ETmin(ir)**2) then
         this%pp(0:3,i)=pi1(0:3)
         this%pp(0:3,ir)=pr1(0:3)
         this%pp(0:3,ib+ir)=ppibir1(0:3)
         this%jac=this%jac/2d0
      elseif (pi2(0)**2-pi2(3)**2.ge.this%ETmin(i)**2 .and. pr2(0)**2-pr2(3)**2.ge.this%ETmin(ir)**2) then
         this%pp(0:3,i)=pi2(0:3)
         this%pp(0:3,ir)=pr2(0:3)
         this%pp(0:3,ib+ir)=ppibir2(0:3)
         this%jac=this%jac/2d0
      else
         this%jac=-19d0
         if (debug) then
            write (*,*) 'piir',this%pp(0:3,i+ir)
            write (*,*) 'pim1',this%pp(0:3,im1)
            write (*,*) '1:',phi1,(phi1+phi2)/(2d0*pi)
            write (*,*) 'i',i,this%ETmin(i),':',pi1(0:3)
            write (*,*) 'ir',ir,this%ETmin(ir),':',pr1(0:3)
            write (*,*) '2:',phi2
            write (*,*) 'i',i,this%ETmin(i),':',pi2(0:3)
            write (*,*) 'ir',ir,this%ETmin(ir),':',pr2(0:3)
            write (*,*) ''
         endif
         return
      endif
      ! Compute the This%Jacobian
      gram4=gram_determinant4(this%invm(ir+i+im1),this%invm(ir+ib),this%invm(ir+i+ib)&
           &,this%invm(ir+i),this%invm(i+im1),this%invm(ir+ib+i+im1),this%invm(ir),this%invm(i)&
           &,this%invm(im1))
      if (gram4.ge.0d0) then 
         write (*,*) 'error, gram4 greater than or equal to zero',gram4,i,ir
         this%jac=-5d0
         return
      endif
      this%jac=this%jac/(8d0*sqrt(-gram4))
    end subroutine gen23_one_step





    subroutine gent_one_step(i,ir,ib)
      ! One step in the usual MadGraph t-channel phase-space generation.
      implicit none
      integer(kind=4),intent(in) :: i,ir,ib
      real(kind=8) :: tmin,tmax,phi,Eimax,shatmin,shatmax,base,etminir,root,y,etmini
      real(kind=8),dimension(0:3) :: piir,pib
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) this%invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax)
         if (popcnt(i+ir).eq.this%next-2) then
            ! The energy of i will be
            ! Ei=(this%sqrtshat+(this%invm(i)-this%invm(ir))/this%sqrtshat)/2d0. This gives a
            ! constraint on the allowed value of this%invm(i), since Ei>ETmin(i)
            shatmin=max(shatmin,max(this%invm(ir),this%invm_min(ir))+this%sqrtshat*(2d0*this%ETmin(i)-this%sqrtshat))
            Eimax=this%sqrtshat-this%ETmin(ir) ! maximum energy for i
            shatmax=min(shatmax,max(this%invm(ir),this%invm_min(ir))+this%sqrtshat*(2d0*Eimax-this%sqrtshat))
         endif
         call generate_mass(i,shatmin,shatmax)
         if (debug) then
            write (*,*) 't - i    ',i,this%invm(i),shatmin,shatmax
         endif
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax)
         if (popcnt(i+ir).eq.this%next-2) then
            ! The energy of ir will be
            ! Eir=(this%sqrtshat+(this%invm(ir)-this%invm(i))/this%sqrtshat)/2d0. This gives a
            ! constraint on the allowed value of this%invm(ir), since Eir>ETmin(ir)
            shatmin=max(shatmin,this%invm(i)+this%sqrtshat*(2d0*this%ETmin(ir)-this%sqrtshat))
            shatmax=min(shatmax,this%invm(i)+this%sqrtshat*(this%sqrtshat-2d0*max(sqrt(this%invm(i)),this%ETmin(i))))
         endif
         call generate_mass(ir,shatmin,shatmax)
         if (debug) then
            write (*,*) 't - ir   ',ir,this%invm(ir),shatmin,shatmax
         endif
      endif
      if (this%jac.le.0d0) return
      call tminmax(this%invm(ir+i),this%invm(ir+i+ib),this%invm(ir),this%invm(i),0d0,tmin,tmax)
      if (this%invm_max(ir+ib).ne.0d0) tmax=min(tmax,this%invm_max(ir+ib))
      if (this%invm_min(ir+ib).ne.0d0) tmin=max(tmin,this%invm_min(ir+ib))
      ! Make sure that the t-range is compatible with the pT cut. Since t is an
      ! invariant we can compute it in any frame. Let's use the frame in which
      ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
      ! ETmin() is boost invariant in the z-direction)
      this%pp(0:3,i+ir)=this%pp(0:3,i+ir+ib)+this%pp(0:3,ib)
      y=log((this%pp(0,i+ir)+this%pp(3,i+ir))/(this%pp(0,i+ir)-this%pp(3,i+ir)))/2d0
      call boostz(this%pp(0,i+ir),y,piir)
      call boostz(this%pp(0,ib),y,pib)
      if ( piir(1)**2+piir(2)**2.lt.this%ETmin(i)**2-this%invm(i) .and. popcnt(i).eq.1 ) then
         etminir=max(this%ETmin(ir),sqrt(this%invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(this%ETmin(i)**2-this%invm(i)))**2) )
      else
         etminir=max(this%ETmin(ir),sqrt(this%invm(ir)))
      endif
      if ( piir(1)**2+piir(2)**2.lt.this%ETmin(ir)**2-this%invm(ir) .and. popcnt(ir).eq.1 ) then
         etmini=max(this%ETmin(i),sqrt(this%invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(this%ETmin(ir)**2-this%invm(ir)))**2) )
      else
         etmini=max(this%ETmin(i),sqrt(this%invm(i)))
      endif
      base=piir(0)**2-ETmini**2+ETminir**2
      ! Note, root=lambda(piir(0)**2,this%ETmin(i)**2,this%ETmin(ir)**2), but the
      ! following is more stable:
      root=(piir(0)-ETmini-ETminir)*(piir(0)+ETmini-ETminir)*&
           (piir(0)-ETmini+ETminir)*(piir(0)+ETmini+ETminir)
      if (root.lt.0d0) then
         this%jac=-33d0
         if (debug) write (*,*) 'root.lt.0d0',root
         return
      endif
      tmin=max(tmin,this%invm(ir)-pib(0)/piir(0)*(base+sqrt(root)))
      tmax=min(tmax,this%invm(ir)-pib(0)/piir(0)*(base-sqrt(root)))
      if (tmin.ge.tmax) then
         this%jac=-3d0
         if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
         return
      endif
      ix=ix+1
      call random_to_var(this%x(ix),ip,tmin,tmax,this%invm(ir+ib),this%jac)
      if (debug) then
         write (*,*) 't- ir+ib',ir+ib,this%invm(ir+ib),tmin,tmax
      endif
      ix=ix+1
      call random_to_var(this%x(ix),0d0,0d0,2d0*pi,phi,this%jac)
      if (debug) then
         write (*,*) 't - phi  ',i,phi,0d0,2d0*pi
      endif
      call gentcms(this%pp(0,ib+ir+i),this%pp(0,ib),this%invm(ib+ir),phi,sqrt(this%invm(i)) &
           &,sqrt(this%invm(ir)),this%pp(0,i),this%pp(0,ib+ir))
      this%pp(0:3,ir)=this%pp(0:3,ib+ir+i)+this%pp(0:3,ib)-this%pp(0:3,i)
      this%jac = this%jac/(4d0*sqrt(lambda(this%invm(ir+i),0d0,this%invm(ir+i+ib))))
    end subroutine gent_one_step

    subroutine gens_one_step(i,ir)
      implicit none
      integer(kind=4),intent(in) :: i,ir
      real(kind=8) :: esum,costh,phi,shatmin,shatmax
      real(kind=8),dimension(0:3) :: p_i,p_ir
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) this%invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax)
         call generate_mass(i,shatmin,shatmax)
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax)
         call generate_mass(ir,shatmin,shatmax)
      endif
      if (this%jac.le.0d0) return
      esum=sqrt(this%invm(i+ir))
      ix=ix+1
      call random_to_var(this%x(ix),0d0,-1d0,1d0,costh,this%jac)
      ix=ix+1
      call random_to_var(this%x(ix),0d0,0d0,2d0*pi,phi,this%jac)
      call mom2cx(esum,sqrt(this%invm(i)),sqrt(this%invm(ir)),costh,phi,p_i,p_ir)
      call boostm(p_i,this%pp(0:3,i+ir),esum,this%pp(0:3,i))
      call boostm(p_ir,this%pp(0:3,i+ir),esum,this%pp(0:3,ir))
      this%jac=this%jac*sqrt(lambda(this%invm(i+ir),this%invm(i),this%invm(ir)))/(8d0*this%invm(i+ir))
      ! fill t-channel stuff to be safe.
      this%pp(0:3,i+1)=this%pp(0:3,i)-this%pp(0:3,1)
      this%pp(0:3,i+2)=this%pp(0:3,i)-this%pp(0:3,2)
      this%pp(0:3,ir+1)=this%pp(0:3,ir)-this%pp(0:3,1)
      this%pp(0:3,ir+2)=this%pp(0:3,ir)-this%pp(0:3,2)
      this%invm(i+1)=dot(this%pp(0:3,i+1),this%pp(0:3,i+1))
      this%invm(i+2)=dot(this%pp(0:3,i+2),this%pp(0:3,i+2))
      this%invm(ir+1)=dot(this%pp(0:3,ir+1),this%pp(0:3,ir+1))
      this%invm(ir+2)=dot(this%pp(0:3,ir+2),this%pp(0:3,ir+2))
    end subroutine gens_one_step

    subroutine generate_mass(i,shatmin,shatmax)
      implicit none
      integer :: i
      real(kind=8) :: shatmin,shatmax
      if (this%invm_min(i).ne.0d0) shatmin=max(shatmin,this%invm_min(i))
      if (this%invm_max(i).ne.0d0) shatmax=min(shatmax,this%invm_max(i))
      if (shatmin.ge.shatmax) then
         this%jac=-7d0
         if (debug) write (*,*) 'shatmin.ge.shatmax',i,shatmin,shatmax
         return
      endif
      ix=ix+1
      call random_to_var(this%x(ix),-0.5d0,shatmin,shatmax,this%invm(i),this%jac)
    end subroutine generate_mass

    subroutine mom2cx(esum,mass1,mass2,costh1,phi1,p1,p2)
      ! This subroutine sets up two four-momenta in the two particle rest
      ! frame.
      ! input:
      !       real    esum           : energy sum of particle 1 and 2
      !       real    mass1          : mass            of particle 1
      !       real    mass2          : mass            of particle 2
      !       real    costh1         : cos(theta)      of particle 1
      !       real    phi1           : azimuthal angle of particle 1
      ! output:
      !       real    p1(0:3)        : four-momentum of particle 1
      !       real    p2(0:3)        : four-momentum of particle 2
      implicit none
      real(kind=8),dimension(0:3) :: p1,p2
      real(kind=8) :: esum,mass1,mass2,costh1,phi1,md2,ed,pp,sinth1
      real(kind=8),parameter :: rZero = 0d0,rHalf = 0.5d0,rOne = 1d0,rTwo = 2d0
      md2 = (mass1-mass2)*(mass1+mass2)
      ed = md2/esum
      if ( mass1*mass2.eq.rZero ) then
         pp = (esum-abs(ed))*rHalf
      else
         pp = sqrt((md2/esum)**2-rTwo*(mass1**2+mass2**2)+esum**2)*rHalf
      endif
      sinth1 = sqrt((rOne-costh1)*(rOne+costh1))
      p1(0) = max((esum+ed)*rHalf,rZero)
      p1(1) = pp*sinth1*cos(phi1)
      p1(2) = pp*sinth1*sin(phi1)
      p1(3) = pp*costh1
      p2(0) = max((esum-ed)*rHalf,rZero)
      p2(1:3) = -p1(1:3)
    end subroutine mom2cx

    subroutine gentcms(pa,pb,t,phi,m1,m2,p1,pr)
      ! Generates 4 momentum for particle p1, and remainder pr=pa-p1=p2-pb given the
      ! values t, and phi in the process pa+pb -> p1+p2.  Assuming incoming
      ! particles with momenta pa, pb and outgoing particles with mass m1,m2;
      ! t=(pb-p2)^2 ; phi is the azimuthal angle between p1 and pa in the
      ! pa+pb rest frame, with pa aligned with the positive z-axis.
      implicit none
      real(kind=8),intent(in) :: t,phi,m1,m2
      real(kind=8),intent(in),dimension(0:3) :: pa,pb
      real(kind=8),intent(out),dimension(0:3) :: p1,pr
      real(kind=8) :: E_acms,p_acms,esum,esum2,ed,pp2,md2,ma2,pt,pt2
      real(kind=8),dimension(0:3) :: ptot,pa_cms,ptotm,p1_rot
      real(kind=8),parameter :: tiny=1d-8
      ptot(0:3)=pa(0:3)+pb(0:3)
      ptotm(0)=ptot(0)
      ptotm(1:3)=-ptot(1:3)
      ma2=dot(pa,pa)
      ! determine magnitude of p1 in cms frame (from dhelas routine mom2cx)
      ESUM2 = dot(ptot,ptot)
      if (esum2 .le. 0d0) then
         write (*,*) "error :: must be time-like momentum in gentcms",esum2
         stop 1
      endif
      esum=sqrt(esum2)
      MD2=(M1-M2)*(M1+M2)
      ED=MD2/ESUM
      IF (M1*M2.EQ.0.d0) THEN
         PP2=0.25d0*(ESUM-ABS(ED))**2
      ELSE
         PP2=0.25d0*((MD2/ESUM)**2-2d0*(M1**2+M2**2)+ESUM**2)
         if(pp2.lt.0d0) then
            write(*,*) 'Error #12 in genps_fks.f: magnitude^2 of '/&
                 &/'3-momentum smaller than 0',pp2
            stop 1
         endif
      ENDIF
      call boostm(pa,ptotm,esum,pa_cms)
      E_acms = pa_cms(0)
      p_acms = sqrt(pa_cms(1)**2+pa_cms(2)**2+pa_cms(3)**2)
      ! define p1 in the frame where pa_cms is aligned with the positive z axis.
      p1(0) = MAX((ESUM+ED)*0.5d0,0.d0)
      if (esum+ed.le.0d0) then
         write (*,*) 'Error #14 in genps_fks.f: negative energy',esum,ed
         this%jac=-8d0
         return
         write (*,*) pa(0:3)
         write (*,*) pb(0:3)
         write (*,*) m1,m2,t,phi
         write (*,*) pa_cms(0:3)
         stop 1
      endif
      p1(3) = -(m1**2+ma2-t-2d0*p1(0)*E_acms)/(2d0*p_acms)
      pt2=pp2-p1(3)**2
      if (pt2/esum2.lt.-tiny) then
         write (*,*) 'Error #13 in genps_fks.f: relative pt^2 smaller than 0',pt2,esum2
         stop 1
      elseif (pt2.lt.0d0) then
         pt2=0d0
      endif
      pt = sqrt(pt2)
      p1(1) = pt*cos(phi)
      p1(2) = pt*sin(phi)
      call rotxxx(p1,pa_cms,p1_rot)       !Rotate p1 to the pa_cms frame
      call boostm(p1_rot,ptot,esum,p1)    !boost back to lab frame
      pr(0:3)=pa(0:3)-p1(0:3)         !Return remainder of momentum
    end subroutine gentcms

    subroutine random_to_var(x,power_in,var_min,var_max,var,jac)
      ! Given a random number x, it generates var in the range var_min
      ! <= var <= var_max according to var^(power)
      implicit none
      real(kind=8),intent(in) :: x,power_in,var_min,var_max
      real(kind=8),intent(out) :: var
      real(kind=8),intent(inout) :: jac
      integer(kind=4) :: ip
      real(kind=8) :: varmin,varmax,power
      if (var_min.lt.0d0 .and. var_max.le.0d0) then
         power=power_in
         varmin=-var_max
         varmax=-var_min
      elseif (var_min.lt.0d0 .and. var_max.gt.0d0 .and. (abs(power_in).gt.vtiny)) then
         write (*,*) 'ERROR: in random_to_var one of the two limits '/&
              &/'is negative',var_min,var_max,power_in,jac,x
         write (*,*) 'using flat transformation'
         power=0d0
         varmin=var_min
         varmax=var_max
      else
         power=power_in
         varmin=var_min
         varmax=var_max
      endif
      ip=nint(power)
      if (dble(ip).eq.power) then
         ! integer
         if (ip.eq.-1) then
            var=varmax**x * varmin**(1d0-x)
            jac=jac*var*log(varmax/varmin)
         elseif (ip.eq.-2) then
            var=1d0/( (1d0-x)/varmin + x/varmax )
            jac=jac*var**2*(varmax-varmin)/(varmax*varmin)
         elseif (ip.eq.0) then
            var=varmin+x*(varmax-varmin)
            jac=jac*(varmax-varmin)
         else
            var=(varmin**(1+ip)*(1d0-x)+varmax**(1+ip)*x)**(1d0/(1d0+power))
            jac=jac*(varmax**(1+ip)-varmin**(1+ip))* &
                 (varmin**(1+ip)*(1d0-x)+varmax**(1+ip)*x)**(-power/(1d0+power))/ &
                 (1d0+power)
         endif
      else
         var=(varmin**(1d0+power)*(1d0-x)+varmax**(1d0+power)*x)**(1d0/(1d0+power))
         jac=jac*(varmax**(1d0+power)-varmin**(1d0+power))* &
              (varmin**(1d0+power)*(1d0-x)+varmax**(1d0+power)*x)**(-power/(1d0+power))/&
              (1d0+power)
      endif
      if (var_min.le.0d0 .and. var_max.le.0d0) then
         var=-var
      endif
    end subroutine random_to_var
  end subroutine gen23_generate_momenta

  subroutine gen23_compute_x_from_momenta(this,p)
    implicit none
    class(phase_space_gen23),intent(inout) :: this
    real(kind=8),dimension(0:3,this%next),intent(in) :: p
    real(kind=8) :: ycm
    integer :: ix
    this%jac=1d0
    ix=0
    ! Fill the full momentum array, including all possible
    ! intermediate states:
    call fill_momentum_array
    ! get the two random number corresponding to the initial state
    if (includePDF) call compute_x_initial_state
    ! The final-state momenta configuration gives all the other random numbers
    call compute_x_final_state
  contains
    subroutine compute_x_initial_state
      implicit none
      real(kind=8) :: tau
      call compute_x_from_tau(tau)
      call compute_x_from_y(tau)
    end subroutine compute_x_initial_state
    subroutine compute_x_from_tau(tau)
      implicit none
      real(kind=8),intent(out) :: tau
      real(kind=8) :: smin,smax,shat
      smin=max(this%invm_min(maskr(this%next)-3),this%ETmin(maskr(this%next)-3)**2)
      smax=this%sqrts**2
      shat=dot(this%pp(0:3,3),this%pp(0:3,3))
      this%sqrtshat=sqrt(shat)
      ix=ix+1
      call var_to_random(shat,ip_shat,smin,smax,this%x(ix),this%jac)
      tau=shat/smax
      this%jac=this%jac/smax
    end subroutine compute_x_from_tau
    subroutine compute_x_from_y(tau)
      implicit none
      real(kind=8),intent(in) :: tau
      real(kind=8) ::  ymin,ymax
      ymin= log(tau)/2d0
      ymax=-log(tau)/2d0
      ix=ix+1
      call var_to_random(ycm,0d0,ymin,ymax,this%x(ix),this%jac)
    end subroutine compute_x_from_y
    subroutine compute_x_final_state
      implicit none
      integer(kind=4) :: i,j,inext,im1
      integer(kind=4),dimension(2) :: set
      ! invariant mass of all final state particles combined
      this%invm(ibset(0,0)+ibset(0,1))=this%sqrtshat**2
      this%invm(maskr(this%next)-ibset(0,0)-ibset(0,1))=this%sqrtshat**2
      set(1)=this%sets(0,1)
      set(2)=this%sets(0,2)
      ! Generate the central 2->2 process in case both set(1) and set(2) are not empty
      if (popcnt(set(1)).gt.1 .and. popcnt(set(2)).gt.1) then
         if (debug) write (*,*) 'two sets with at least two ',&
              & 'particles',popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         if (use_t_channel_at_start) then
            call gent_one_step_inverse(set(2),set(1),1)
         else
            call gens_one_step_inverse(set(2),set(1))
         endif
         this%invm(set(2)+2)=dot(this%pp(0:3,set(2)+2),this%pp(0:3,set(2)+2))
      elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).gt.1) then
         if (debug) write (*,*) 'special double t-channel (1)'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call double_t_inverse(set(1),set(2),1,2)
         this%invm(set(2)+2)=dot(this%pp(0:3,set(2)+2),this%pp(0:3,set(2)+2))
      elseif (popcnt(set(1)).gt.1 .and. popcnt(set(2)).eq.1) then
         if (debug) write (*,*) 'special double t-channel (2)'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call double_t_inverse(set(2),set(1),1,2)
         this%invm(set(1)+1)=dot(this%pp(0:3,set(1)+1),this%pp(0:3,set(1)+1))
      elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).eq.1) then
         if (debug) write (*,*) '2->2 scattering with one particle in each set'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call gent_one_step_inverse(set(2),set(1),1)
         this%invm(set(2)+2)=dot(this%pp(0:3,set(2)+2),this%pp(0:3,set(2)+2))
      endif

      do i=1,2
         if (popcnt(set(i)).le.1) cycle ! at least 2 particles in a set
         inext=ibset(0,this%sets(1,i)-1)
         set(i)=set(i)-inext
         if (popcnt(set(i)).ge.2) then
            ! at least 3 particles in a set
            if (debug) write (*,*) 'At least 3 particles in a set',&
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            call gent_one_step_inverse(set(i),inext,i)
            this%invm((3-i)+set(i)+inext)=dot(this%pp(0:3,(3-i)+set(i)+inext),this%pp(0:3,(3-i)+set(i)+inext))
            this%invm((3-i)+set(i))=dot(this%pp(0:3,(3-i)+set(i)),this%pp(0:3,(3-i)+set(i)))
            do j=2,popcnt(set(i))-1
               ! loop over the remaining particles in the set
               inext=ibset(0,this%sets(j,i)-1)
               im1=ibset(0,this%sets(j-1,i)-1)
               set(i)=set(i)-inext
               if (this%t_channel) then
                  call gent_one_step_inverse(inext,set(i),3-i)
               else
                  call gen23_one_step_inverse(inext,set(i),3-i,im1)
               endif
            enddo
            inext=ibset(0,this%sets(j,i)-1)
            im1=ibset(0,this%sets(j-1,i)-1)
            set(i)=set(i)-inext
            if (this%t_channel) then
               call gent_one_step_inverse(inext,set(i),3-i)
            else
               call gen23_one_step_inverse(inext,set(i),3-i,im1)
            endif
         elseif (popcnt(set(i)).eq.1 .and. popcnt(this%sets(0,3-i)).ne.0) then
            ! Exactly 2 particles in a set (and the other set contains at least one)
            if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
                 & 'the other set contains at least one)', &
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            im1=3-i
            this%invm(set(i)+inext+im1)=dot(this%pp(0:3,set(i)+inext+im1),this%pp(0:3,set(i)+inext+im1))
            this%invm(set(i)+inext+i+im1)=dot(this%pp(0:3,set(i)+inext+im1+i),this%pp(0:3,set(i)+inext+im1+i))
            if (this%t_channel) then
               call gent_one_step_inverse(set(i),inext,i)
            else
               call gen23_one_step_inverse(set(i),inext,i,im1)
            endif
         elseif (popcnt(set(i)).eq.1 .and. popcnt(this%sets(0,3-i)).eq.0) then
            ! Exactly 2 particles in a set (and the other set contains none)
            if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
                 & 'the other set contains none)', &
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            call gent_one_step_inverse(set(i),inext,i)
         else
            write (*,*) 'Inconsistent sets'
            write (*,*) i,':',this%sets(:,i)
            write (*,*) 3-i,':',this%sets(:,3-i)
            stop 
         endif
      enddo

      ! Add factors of 2*pi
      this%jac=this%jac/((2d0*pi)**(3*(this%next-2)-4))
      ! Add flux factor
      this%jac=this%jac/(2d0*this%sqrtshat**2)
      
    end subroutine compute_x_final_state
    subroutine gent_one_step_inverse(i,ir,ib)
      implicit none
      integer(kind=4),intent(in) :: i,ir,ib
      real(kind=8) :: tmin,tmax,phi,Eimax,shatmin,shatmax,base,etminir,root,y,etmini,esum
      real(kind=8),dimension(0:3) :: piir,pib,p_boost,pi_rot,pa_cms
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) this%invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax)
         if (popcnt(i+ir).eq.this%next-2) then
            shatmin=max(shatmin,max(this%invm(ir),this%invm_min(ir))+this%sqrtshat*(2d0*this%ETmin(i)-this%sqrtshat))
            Eimax=this%sqrtshat-this%ETmin(ir) ! maximum energy for i
            shatmax=min(shatmax,max(this%invm(ir),this%invm_min(ir))+this%sqrtshat*(2d0*Eimax-this%sqrtshat))
         endif
         call generate_mass_inverse(i,shatmin,shatmax)
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax)
         if (popcnt(i+ir).eq.this%next-2) then
            shatmin=max(shatmin,this%invm(i)+this%sqrtshat*(2d0*this%ETmin(ir)-this%sqrtshat))
            shatmax=min(shatmax,this%invm(i)+this%sqrtshat*(this%sqrtshat-2d0*max(sqrt(this%invm(i)),this%ETmin(i))))
         endif
         call generate_mass_inverse(ir,shatmin,shatmax)
      endif
      call tminmax(this%invm(ir+i),this%invm(ir+i+ib),this%invm(ir),this%invm(i),0d0,tmin,tmax)
      if (this%invm_max(ir+ib).ne.0d0) tmax=min(tmax,this%invm_max(ir+ib))
      if (this%invm_min(ir+ib).ne.0d0) tmin=max(tmin,this%invm_min(ir+ib))
      ! Make sure that the t-range is compatible with the pT cut. Since t is an
      ! invariant we can compute it in any frame. Let's use the frame in which
      ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
      ! ETmin() is boost invariant in the z-direction)
      y=log((this%pp(0,i+ir)+this%pp(3,i+ir))/(this%pp(0,i+ir)-this%pp(3,i+ir)))/2d0
      call boostz(this%pp(0,i+ir),y,piir)
      call boostz(this%pp(0,ib),y,pib)
      if ( piir(1)**2+piir(2)**2.lt.this%ETmin(i)**2-this%invm(i) .and. popcnt(i).eq.1 ) then
         etminir=max(this%ETmin(ir),sqrt(this%invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(this%ETmin(i)**2-this%invm(i)))**2) )
      else
         etminir=max(this%ETmin(ir),sqrt(this%invm(ir)))
      endif
      if ( piir(1)**2+piir(2)**2.lt.this%ETmin(ir)**2-this%invm(ir) .and. popcnt(ir).eq.1 ) then
         etmini=max(this%ETmin(i),sqrt(this%invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(this%ETmin(ir)**2-this%invm(ir)))**2) )
      else
         etmini=max(this%ETmin(i),sqrt(this%invm(i)))
      endif
      base=piir(0)**2-ETmini**2+ETminir**2
      ! Note, root=lambda(piir(0)**2,this%ETmin(i)**2,this%ETmin(ir)**2), but the
      ! following is more stable:
      root=(piir(0)-ETmini-ETminir)*(piir(0)+ETmini-ETminir)*&
           (piir(0)-ETmini+ETminir)*(piir(0)+ETmini+ETminir)
      if (root.lt.0d0) then
          write (*,*) 'root.lt.0d0 in gent_one_step_inverse',root
         stop 1
      endif
      tmin=max(tmin,this%invm(ir)-pib(0)/piir(0)*(base+sqrt(root)))
      tmax=min(tmax,this%invm(ir)-pib(0)/piir(0)*(base-sqrt(root)))
      this%invm(ir+ib)=dot(this%pp(0:3,ir+ib),this%pp(0:3,ir+ib))
      if (tmin.ge.tmax) then
         write (*,*) 'tmin.ge.tmax in gent_one_step_inverse',tmin,tmax,this%invm(ir+ib)
         stop 1
      endif
      ix=ix+1
      call var_to_random(this%invm(ir+ib),ip,tmin,tmax,this%x(ix),this%jac)
      ! inverse of boosts and rotation from gentcms()
      esum=sqrt(this%invm(i+ir))
      p_boost(0)=this%pp(0,i+ir)
      p_boost(1:3)=-this%pp(1:3,i+ir)
      call boostm(this%pp(0:3,i),p_boost,esum,pib)
      call boostm(this%pp(0:3,ib+ir+i),p_boost,esum,pa_cms)
      call rotxxx_inv(pib,pa_cms,pi_rot)
      phi=atan(pi_rot(2)/pi_rot(1))
      if(pi_rot(1).lt.0d0) phi=phi+pi
      if(phi.lt.0d0) phi=phi+2d0*pi
      ix=ix+1
      call var_to_random(phi,0d0,0d0,2d0*pi,this%x(ix),this%jac)
      this%jac = this%jac/(4d0*sqrt(lambda(this%invm(ir+i),0d0,this%invm(ir+i+ib))))
    end subroutine gent_one_step_inverse
    subroutine gens_one_step_inverse(i,ir)
      implicit none
      integer(kind=4),intent(in) :: i,ir
      real(kind=8) :: esum,costh,phi,shatmin,shatmax
      real(kind=8),dimension(0:3) :: p_i,p_boost
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) this%invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax)
         call generate_mass_inverse(i,shatmin,shatmax)
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax)
         call generate_mass_inverse(ir,shatmin,shatmax)
      endif
      ! boost p(i) and p(ir) to the p(i+ir) rest frame
      esum=sqrt(this%invm(i+ir))
      p_boost(0)=-this%pp(0,i+ir)
      p_boost(1:3)=-this%pp(1:3,i+ir)
      call boostm(this%pp(0:3,i),p_boost,esum,p_i)
      ! compute the angles from the momenta and use that to get the random numbers
      costh=p_i(3)/sqrt(p_i(1)**2+p_i(2)**2+p_i(3)**2)
      ix=ix+1
      call var_to_random(costh,0d0,-1d0,1d0,this%x(ix),this%jac)
      phi=atan(p_i(2)/p_i(1))
      if(p_i(1).lt.0d0) phi=phi+pi
      if(phi.lt.0d0) phi=phi+2d0*pi
      ix=ix+1
      call var_to_random(phi,0d0,0d0,2d0*pi,this%x(ix),this%jac)
      ! update the Jacobian
      this%jac=this%jac*sqrt(lambda(this%invm(i+ir),this%invm(i),this%invm(ir)))/(8d0*this%invm(i+ir))
      ! compute some t-channel invariants just to make sure they are filled. 
      this%invm(i+1)=dot(this%pp(0:3,i+1),this%pp(0:3,i+1))
      this%invm(i+2)=dot(this%pp(0:3,i+2),this%pp(0:3,i+2))
      this%invm(ir+1)=dot(this%pp(0:3,ir+1),this%pp(0:3,ir+1))
      this%invm(ir+2)=dot(this%pp(0:3,ir+2),this%pp(0:3,ir+2))
    end subroutine gens_one_step_inverse
    subroutine double_t_inverse(i,ir,ia,ib)
      implicit none
      integer(kind=4),intent(in) :: i,ir,ia,ib
      real(kind=8) :: tmin,tmax,phi,yr,Eimax,pzmax
      if (popcnt(i).ne.1 .or. popcnt(ir).le.1) then
         write (*,*) 'Subroutine only for i is a single particle '&
              //'and ir is more than 1',i,ir,popcnt(i),popcnt(ir)
         stop 1
      endif
      yr=sqrt(lambda(this%invm(ia+ib),this%invm(i),this%invm_min(ir)))
      tmin=(-this%invm(ia+ib)+this%invm(i)+this%invm_min(ir)-yr)/2d0
      tmax=(-this%invm(ia+ib)+this%invm(i)+this%invm_min(ir)+yr)/2d0
      if (this%invm_max(ir+ib).ne.0d0) tmax=min(tmax,this%invm_max(ir+ib))
      if (this%invm_min(ir+ib).ne.0d0) tmin=max(tmin,this%invm_min(ir+ib))
      pzmax=sqrt(lambda(this%sqrtshat**2,this%ETmin(i)**2,this%ETmin(ir)**2))/(2d0*this%sqrtshat)
      Eimax=this%sqrtshat-sqrt(this%ETmin(ir)**2+pzmax**2)
      tmin=max(tmin,this%invm(i)-this%sqrtshat*(Eimax+pzmax))
      tmax=min(tmax,this%invm(i)-this%sqrtshat*(Eimax-pzmax))
      if (tmin.ge.tmax) then
         write (*,*) 'tmin.ge.tmax in double_t_inverse',tmin,tmax
         stop 1
      endif
      this%invm(i+ia)=dot(this%pp(0:3,i+ia),this%pp(0:3,i+ia))
      ix=ix+1
      call var_to_random(this%invm(i+ia),ip,tmin,tmax,this%x(ix),this%jac)
      tmin=-this%invm(ia+ib)-this%invm(i+ia)+this%invm(i)+this%invm_min(ir)
      tmax=this%invm(i)*(this%invm(i)-this%invm(ia+ib)-this%invm(i+ia))/(this%invm(i)-this%invm(i+ia))
      if (this%invm_max(ir+ib).ne.0d0) tmax=min(tmax,this%invm_max(ir+ib))
      if (this%invm_min(ir+ib).ne.0d0) tmin=max(tmin,this%invm_min(ir+ib))
      ! Additional constraints on tmin and tmax due to pp(0,i) and pp(0,ir)
      ! being larger than ETmin(i) and ETmin(ir), respectively:
      tmin=max(tmin,this%invm(i)-this%sqrtshat**2*(1-this%ETmin(ir)**2/(this%sqrtshat**2+this%invm(i+ia)-this%invm(i))))
      tmax=min(tmax,this%invm(i)-this%sqrtshat**2*(this%ETmin(i)**2/(this%invm(i)-this%invm(i+ia))))
      if (tmin.ge.tmax) then
         write (*,*) 'tmin.ge.tmax in double_t_inverse',tmin,tmax
         stop 1
      endif
      this%invm(i+ib)=dot(this%pp(0:3,i+ib),this%pp(0:3,i+ib))
      ix=ix+1
      call var_to_random(this%invm(i+ib),ip,tmin,tmax,this%x(ix),this%jac)
      phi=atan(this%pp(2,i)/this%pp(1,i))
      if(this%pp(1,i).lt.0d0) phi=phi+pi
      if(phi.lt.0d0) phi=phi+2d0*pi
      ix=ix+1
      call var_to_random(phi,0d0,0d0,2d0*pi,this%x(ix),this%jac)
      this%invm(ir)=dot(this%pp(0,ir),this%pp(0,ir))
      if (this%invm(ir).le.0d0) then
         write (*,*) "ERROR in double_t: invariant mass of system", &
              & " must be larger than zero",ir,this%invm(ir),i
         write (*,*) this%invm(ir),this%invm(ia+ib)+this%invm(i+ia)+this%invm(i+ib)-this%invm(i)&
              &,this%invm(ia+ib),this%invm(i +ia),this%invm(i+ib),this%invm(i)
         stop
      endif
      this%jac = this%jac/(4d0*sqrt(lambda(this%invm(ir+i),0d0,0d0)))
    end subroutine double_t_inverse
    subroutine gen23_one_step_inverse(i,ir,ib,im1)
      implicit none
      integer(kind=4),intent(in) :: im1,i,ir,ib
      real(kind=8) :: tmin,tmax,smin,smax,phi1,phi2,gram4,V,sqrtGG,shatmin,shatmax,y,base,root,phi_rot,&
           etminir,etmini
      real(kind=8),dimension(0:3) :: pi1,pr1,ppibir1,pi2,pr2,ppibir2,piir,pib,pim1,piirr,pim1r
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) this%invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax)
         call generate_mass_inverse(i,shatmin,shatmax)
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax)
         call generate_mass_inverse(ir,shatmin,shatmax)
      endif
      call tminmax(this%invm(ir+i),this%invm(ir+i+ib),this%invm(ir),this%invm(i),0d0,tmin,tmax)
      if (this%invm_max(ir+ib).ne.0d0) tmax=min(tmax,this%invm_max(ir+ib))
      if (this%invm_min(ir+ib).ne.0d0) tmin=max(tmin,this%invm_min(ir+ib))
      this%pp(0:3,i+ir)=this%pp(0:3,i+ir+ib)+this%pp(0:3,ib)
      y=log((this%pp(0,i+ir)+this%pp(3,i+ir))/(this%pp(0,i+ir)-this%pp(3,i+ir)))/2d0
      call boostz(this%pp(0,i+ir),y,piir)
      call boostz(this%pp(0,ib),y,pib)
      if ( piir(1)**2+piir(2)**2.lt.this%ETmin(i)**2-this%invm(i) .and. popcnt(i).eq.1 ) then
         etminir=max(this%ETmin(ir),sqrt(this%invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(this%ETmin(i)**2-this%invm(i)))**2) )
      else
         etminir=max(this%ETmin(ir),sqrt(this%invm(ir)))
      endif
      if ( piir(1)**2+piir(2)**2.lt.this%ETmin(ir)**2-this%invm(ir) .and. popcnt(ir).eq.1 ) then
         etmini=max(this%ETmin(i),sqrt(this%invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(this%ETmin(ir)**2-this%invm(ir)))**2) )
      else
         etmini=max(this%ETmin(i),sqrt(this%invm(i)))
      endif
      base=piir(0)**2-ETmini**2+ETminir**2
      ! Note, root=lambda(piir(0)**2,this%ETmin(i)**2,this%ETmin(ir)**2), but the
      ! following is more stable:
      root=(piir(0)-ETmini-ETminir)*(piir(0)+ETmini-ETminir)*&
           (piir(0)-ETmini+ETminir)*(piir(0)+ETmini+ETminir)
      if (root.lt.0d0) then
         write (*,*) 'root.lt.0d0 in gen23_one_step_inverse',root
         stop 1
      endif
      tmin=max(tmin,this%invm(ir)-pib(0)/piir(0)*(base+sqrt(root)))
      tmax=min(tmax,this%invm(ir)-pib(0)/piir(0)*(base-sqrt(root)))
      if (tmin.ge.tmax) then
         if (debug) write (*,*) 'tmin.ge.tmax in gen23_one_step_inverse',tmin,tmax
         stop 1
      endif
      this%invm(ir+ib)=dot(this%pp(0:3,ir+ib),this%pp(0:3,ir+ib))
      ix=ix+1
      call var_to_random(this%invm(ir+ib),ip,tmin,tmax,this%x(ix),this%jac)
      call sminmax(this%invm(ir+i),this%invm(ir),this%invm(ir+i+im1),this%invm(ir+i+ib)&
           &,this%invm(ir+ib),this%invm(ir+ib+i+im1),this%invm(i),this%invm(im1),smin,smax,V,sqrtGG)
      if (this%invm_min(i+im1).ne.0d0) smin=max(smin,this%invm_min(i+im1))
      if (this%invm_max(i+im1).ne.0d0) smax=min(smax,this%invm_max(i+im1))
      if (im1.gt.2) then
         ! Boost and rotate in z-direction such that pp(:,im1) goes in the x-direction.
         y=log((this%pp(0,im1)+this%pp(3,im1))/(this%pp(0,im1)-this%pp(3,im1)))/2d0
         call boostz(this%pp(0,i+ir),y,piirr)
         call boostz(this%pp(0,im1),y,pim1r)
         call boostz(this%pp(0,ib),y,pib)
         phi_rot=atan(this%pp(2,im1)/this%pp(1,im1))
         if(this%pp(1,im1).lt.0d0) phi_rot=phi_rot+pi
         call rotz(piirr,-phi_rot,piir)
         call rotz(pim1r,-phi_rot,pim1)
         ! Eir > Etmin(ir) + constraint coming from t
         etminir=max(pib(0)*this%ETmin(ir)**2/(this%invm(ir)-this%invm(ir+ib))+(this%invm(ir)-this%invm(ir+ib))/(4d0*pib(0)),&
              this%ETmin(ir))
         smax=min(smax,&
              this%invm(i)+this%invm(im1)+2d0*(piir(0)-etminir)*pim1(0)+2d0*sqrt((piir(0)-etminir)**2-this%invm(i))*pim1(1))
         if(this%invm(i).eq.0d0) then
            smin=max(smin,2d0*this%ETmin(i)*(pim1(0)-pim1(1)*cos(this%drcut)))
         endif
      endif
      if (smin.ge.smax) then
         write (*,*) 'smin.ge.smax in gen23_one_stop_inverse',smin,smax
         stop 1
      endif
      this%invm(i+im1)=dot(this%pp(0:3,i+im1),this%pp(0:3,i+im1))
      ix=ix+1
      call var_to_random(this%invm(i+im1),ip,smin,smax,this%x(ix),this%jac)
      ! Generate the momenta from the integration variables. Since there is an
      ! ambiguity in phi, get both of them and pick the one that passes the cuts
      ! (if it's only one). If both pass, simply pick one of the two at random
      ! with a flat prior.
      phi1=getphifroms(this%invm(i+im1),this%invm(ir+i),this%invm(ir),this%invm(ir+i+im1)&
           &,this%invm(ir+i+ib),V,sqrtGG,1d0)
      call gentcms2(this%pp(0,ib),this%pp(0,ib+ir+i),this%pp(0,ib+ir+i+im1),this%invm(ir+ib),phi1 &
           &,sqrt(this%invm(i)),sqrt(this%invm(ir)),pi1,ppibir1)
      pr1(0:3)=this%pp(0:3,ir+i)-pi1(0:3)
      phi2=getphifroms(this%invm(i+im1),this%invm(ir+i),this%invm(ir),this%invm(ir+i+im1)&
           &,this%invm(ir+i+ib),V,sqrtGG,0d0)
      call gentcms2(this%pp(0,ib),this%pp(0,ib+ir+i),this%pp(0,ib+ir+i+im1),this%invm(ir+ib),phi2 &
           &,sqrt(this%invm(i)),sqrt(this%invm(ir)),pi2,ppibir2)
      pr2(0:3)=this%pp(0:3,ir+i)-pi2(0:3)
      if ( pi1(0)**2-pi1(3)**2.ge.this%ETmin(i)**2 .and. pr1(0)**2-pr1(3)**2.ge.this%ETmin(ir)**2 .and. &
           pi2(0)**2-pi2(3)**2.ge.this%ETmin(i)**2 .and. pr2(0)**2-pr2(3)**2.ge.this%ETmin(ir)**2 ) then
         continue
      elseif (pi1(0)**2-pi1(3)**2.ge.this%ETmin(i)**2 .and. pr1(0)**2-pr1(3)**2.ge.this%ETmin(ir)**2) then
         this%jac=this%jac/2d0
      elseif (pi2(0)**2-pi2(3)**2.ge.this%ETmin(i)**2 .and. pr2(0)**2-pr2(3)**2.ge.this%ETmin(ir)**2) then
         this%jac=this%jac/2d0
      endif
      ! Compute the This%Jacobian
      gram4=gram_determinant4(this%invm(ir+i+im1),this%invm(ir+ib),this%invm(ir+i+ib)&
           &,this%invm(ir+i),this%invm(i+im1),this%invm(ir+ib+i+im1),this%invm(ir),this%invm(i)&
           &,this%invm(im1))
      if (gram4.ge.0d0) then 
         write (*,*) 'error, gram4 greater than or equal to zero in gen23_one_step_inverse',gram4,i,ir
         this%jac=-5d0
         return
      endif
      this%jac=this%jac/(8d0*sqrt(-gram4))
    end subroutine gen23_one_step_inverse
    subroutine fill_momentum_array
      implicit none
      integer :: i,j
      real(kind=8),dimension(0:3) :: pp
      ycm=log((p(0,1)+p(0,2)+p(3,1)+p(3,2))/(p(0,1)+p(0,2)-p(3,1)-p(3,2)))/2d0
      do i=1,maskr(this%next)
         pp(0:3)=0d0
         do j=0,this%next-1
            if (btest(i,j)) then
               if (j.le.1 .and. popcnt(i).ne.1) then
                  pp(0:3)=pp(0:3)-p(0:3,j+1)
               else
                  pp(0:3)=pp(0:3)+p(0:3,j+1)
               endif
            endif
         enddo
         call boostz(pp(0:3),ycm,this%pp(0:3,i))
      enddo
    end subroutine fill_momentum_array
    subroutine var_to_random(variable,power_in,var_min,var_max,x,jac)
      ! Given a random variable var between varmin and varmax, compute
      ! the corresponding value of x between 0 and 1.
      implicit none
      real(kind=8),intent(in) :: variable,power_in,var_min,var_max
      real(kind=8),intent(out) :: x
      real(kind=8),intent(inout) :: jac
      integer(kind=4) :: ip
      real(kind=8) :: varmin,varmax,power,var
      if (variable.lt.var_min .or. variable.gt.var_max) then
         write (*,*) 'variable not between varmin and varmax',var_min,variable,var_max
         stop 1
      endif
      if (var_min.lt.0d0 .and. var_max.le.0d0) then
         power=power_in
         varmin=-var_max
         varmax=-var_min
         var=-variable
      elseif (var_min.lt.0d0 .and. var_max.gt.0d0 .and. (abs(power_in).gt.vtiny)) then
         write (*,*) 'ERROR: in var_to_random one of the two limits '/&
              &/'is negative',var_min,var_max,power_in,jac,x
         write (*,*) 'using flat transformation'
         power=0d0
         varmin=var_min
         varmax=var_max
         var=variable
      else
         power=power_in
         varmin=var_min
         varmax=var_max
         var=variable
      endif
      ip=nint(power)
      if (dble(ip).eq.power) then
         ! integer
         if (ip.eq.-1) then
            x=log(varmin/var)/log(varmin/varmax)
            jac=jac*var*log(varmax/varmin)
         elseif (ip.eq.-2) then
            x=varmax/var * ((var-varmin)/(varmax-varmin))
            jac=jac*var**2*(varmax-varmin)/(varmax*varmin)
         elseif (ip.eq.0) then
            x=(var-varmin)/(varmax-varmin)
            jac=jac*(varmax-varmin)
         else
            x=(var**(1+ip)-varmin**(1+ip))/(varmax**(1+ip)-varmin**(1+ip))
            jac=jac*(varmax**(1+ip)-varmin**(1+ip))* &
                 (varmin**(1+ip)*(1d0-x)+varmax**(1+ip)*x)**(-power/(1d0+power))/ &
                 (1d0+power)
         endif
      else
         x=(var**(1d0+power)-varmin**(1d0+power))/(varmax**(1d0+power)-varmin**(1d0+power))
         jac=jac*(varmax**(1d0+power)-varmin**(1d0+power))* &
              (varmin**(1d0+power)*(1d0-x)+varmax**(1d0+power)*x)**(-power/(1d0+power))/&
              (1d0+power)
      endif
    end subroutine var_to_random
    subroutine generate_mass_inverse(i,shatmin,shatmax)
      implicit none
      integer :: i
      real(kind=8) :: shatmin,shatmax
      if (this%invm_min(i).ne.0d0) shatmin=max(shatmin,this%invm_min(i))
      if (this%invm_max(i).ne.0d0) shatmax=min(shatmax,this%invm_max(i))
      this%invm(i)=dot(this%pp(0:3,i),this%pp(0:3,i))
      ix=ix+1
      call var_to_random(this%invm(i),-0.5d0,shatmin,shatmax,this%x(ix),this%jac)
    end subroutine generate_mass_inverse
    
  end subroutine gen23_compute_x_from_momenta
  real(kind=8) function dot(p1,p2)
    ! Inner product between two 4-vectors
    implicit none
    real(kind=8),intent(in),dimension(0:3) :: p1,p2
    dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
  end function dot
  subroutine shatminmax(this,j1,j2,shatmin,shatmax)
    ! Determines minimum and maximum allowed s-channel invariant
    ! masses based on previously generated masses and the masses of
    ! the final-state particles.
    implicit none
    class(phase_space_gen23),intent(in) :: this
    integer(kind=4),intent(in) :: j1,j2
    real(kind=8),intent(out) :: shatmin,shatmax
    integer(kind=4) :: j
    shatmin=0d0
    do j=0,this%next-1
       if (btest(j1,j)) then
          shatmin=shatmin+sqrt(this%invm(ibset(0,j)))
       endif
    enddo
    shatmin=shatmin**2
    shatmax=(sqrt(this%invm(j1+j2))-sqrt(max(this%invm(j2),this%invm_min(j2))))**2
  end subroutine shatminmax
  real(kind=8) function lambda(s,xa2,xb2)
    ! The usual two dimensional phase-space volume factor. See, e.g.,
    ! Eq.(A2) of E.~Byckling and K.~Kajantie, ``Reductions of the
    ! phase-space integral in terms of simpler processes,'' Phys. Rev. 187
    ! (1969), 2008-2016, doi:10.1103/PhysRev.187.2008
    implicit none
    real(kind=8) :: xa2,xb2,S
    lambda=s**2-2d0*(xa2+xb2)*s+(xa2-xb2)**2
  end function lambda

  subroutine boostm(p,q,m,pboost)
    ! This subroutine performs the Lorentz boost of a four-momentum.  The
    ! momentum p is assumed to be given in the rest frame of q.  pboost is
    ! the momentum p boosted to the frame in which q is given.  q must be a
    ! timelike momentum.
    ! input:
    !       real    p(0:3)         : four-momentum p in the q rest  frame
    !       real    q(0:3)         : four-momentum q in the boosted frame
    !       real    m              : mass of q (for numerical stability)
    ! output:
    !       real    pboost(0:3)    : four-momentum p in the boosted frame
    implicit none
    real(kind=8) ,dimension(0:3) :: p,q,pboost
    real(kind=8) :: pq,qq,m,lf
    real(kind=8),parameter :: rZero=0d0
    qq = q(1)**2+q(2)**2+q(3)**2
    if ( qq.ne.rZero ) then
       pq = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
       lf = ((q(0)-m)*pq/qq+p(0))/m
       pboost(0) = (p(0)*q(0)+pq)/m
       pboost(1:3) =  p(1:3)+q(1:3)*lf
    else
       pboost(0:3) = p(0:3)
    endif
  end subroutine boostm
  subroutine boostz(p,yb,pb)
    ! boost in the z-direction with rapidity yb
    implicit none
    real(kind=8),dimension(0:3) :: p,pb
    real(kind=8) :: yb
    pb(0)=p(0)*cosh(yb)-p(3)*sinh(yb)
    pb(1:2)=p(1:2)
    pb(3)=p(3)*cosh(yb)-p(0)*sinh(yb)
  end subroutine boostz
  real(kind=8) function gram_determinant4(shat_ip1,t_im1,t_i,shat_i,s_i,t_ip1&
       &,shat_im1,m_i_2,m_ip1_2)
    ! Computes the 4x4 Gram determinant (with m_a=0) as in eq.(B6), of
    ! E.~Byckling and K.~Kajantie, ``Reductions of the phase-space
    ! integral in terms of simpler processes,'' Phys. Rev. 187 (1969),
    ! 2008-2016, doi:10.1103/PhysRev.187.2008.
    use LUPdecomposition
    implicit none
    real(kind=8) :: shat_ip1,t_im1,t_i,shat_i,s_i,t_ip1,shat_im1,m_i_2,m_ip1_2
    integer(kind=4),parameter :: n=4
    real(kind=8),dimension(n,n) :: a
    integer(kind=4),dimension(0:n) :: p
    real(kind=8),parameter :: tol=1d-8
    real(kind=8) :: deter
    logical :: success
    a(1:4,1)=(/ 0d0           , t_im1-shat_im1 , t_i-shat_i       , t_ip1-shat_ip1    /)
    a(1:4,2)=(/ t_im1-shat_im1, 2d0*t_im1      , t_i+t_im1-m_i_2  , t_im1+t_ip1-s_i   /)
    a(1:4,3)=(/ t_i-shat_i    , t_i+t_im1-m_i_2, 2d0*t_i          , t_i+t_ip1-m_ip1_2 /)
    a(1:4,4)=(/ t_ip1-shat_ip1, t_im1+t_ip1-s_i, t_i+t_ip1-m_ip1_2, 2d0*t_ip1         /)
    call LUPdecompose(a,n,tol,p,success)
    if (success) then
       call LUPdeterminant(a,p,n,deter)
       gram_determinant4=deter/16d0
    else
       gram_determinant4=1d0
    endif
  end function gram_determinant4
  subroutine tminmax(X,Z,U,V,W,tmin,tmax)
    ! Determines the integration limits for t_i, as determined by
    ! setting cos(theta) to +/- 1 in eq.(10) of E.~Byckling and
    ! K.~Kajantie, ``Reductions of the phase-space integral in terms
    ! of simpler processes,'' Phys. Rev. 187 (1969), 2008-2016,
    ! doi:10.1103/PhysRev.187.2008.
    ! x = shat_i
    ! z = t_i
    ! u = shat_{i-1}
    ! v = m_i^2
    ! w = m_a^2
    implicit none
    real(kind=8),intent(in) :: x,z,u,v,w 
    real(kind=8),intent(out):: tmin,tmax
    real(kind=8) ::  t1,t2,yr
    yr = lambda(x,u,v)*lambda(x,w,z)
    if (yr.le.0d0) then
       write (*,*) 'No allowed range for t: tmin=tmax',yr
!!$       stop 1
       yr=0d0
    endif
    yr=sqrt(yr)
    t1 = u+w - ((x+u-v)*(x+w-z) - yr)/(2d0*x)
    t2 = u+w - ((x+u-v)*(x+w-z) + yr)/(2d0*x)
    tmin = min(t1,t2)
    tmax = max(t1,t2)
  end subroutine tminmax
  real(kind=8) function computeV(shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2)
    ! Computes the determinant of V (with m_a=0) based on eq.(11) of
    ! E.~Byckling and K.~Kajantie, ``Reductions of the phase-space
    ! integral in terms of simpler processes,'' Phys. Rev. 187 (1969),
    ! 2008-2016, doi:10.1103/PhysRev.187.2008
    use LUPdecomposition
    real(kind=8),intent(in) :: shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2
    real(kind=8) :: deter
    integer(kind=4),parameter :: n=3
    real(kind=8),dimension(n,n) :: a
    integer(kind=4),dimension(0:n) :: p
    real(kind=8),parameter :: tol=1d-10
    logical :: success
    a(1:3,1)=(/2d0*shat_i             , shat_i-t_i    , shat_i+shat_im1-m_i_2/)
    a(1:3,2)=(/shat_i-t_i             , 0d0           , shat_im1-t_im1       /)
    a(1:3,3)=(/shat_ip1+shat_i-m_ip1_2, shat_ip1-t_ip1, 0d0                  /)
    call LUPdecompose(a,n,tol,p,success)
    if (success) then
       call LUPdeterminant(a,p,n,deter)
       computeV=-deter/8d0
    else
       computeV=-99d99
    endif
  end function computeV
  real(kind=8) function G(x,y,z,u,v,w)
    ! The G-function, eq.(A5) of E.~Byckling and K.~Kajantie,
    ! ``Reductions of the phase-space integral in terms of simpler
    ! processes,'' Phys. Rev. 187 (1969), 2008-2016,
    ! doi:10.1103/PhysRev.187.2008
    implicit none
    real(kind=8),intent(in) :: x,y,z,u,v,w
    G=x**2*y+x*y**2+z**2*u+z*u**2+v**2*w+v*w**2&
         & +x*z*w+x*u*v+y*z*v+y*u*w &
         & -x*y*(z+u+v+w)-z*u*(x+y+v+w)-v*w*(x+y+z+u)
  end function G
  subroutine sminmax(shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2&
       &,m_ip1_2,smin,smax,V,sqrtGG)
    ! Determines the integration limits for s_i, as determined by
    ! setting cos(phi) to +/- 1 in eq.(11) of E.~Byckling and
    ! K.~Kajantie, ``Reductions of the phase-space integral in terms
    ! of simpler processes,'' Phys. Rev. 187 (1969), 2008-2016,
    ! doi:10.1103/PhysRev.187.2008.
    implicit none
    real(kind=8),intent(in) :: shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2
    real(kind=8),intent(out) :: smin,smax,V,sqrtGG
    real(kind=8) :: GG,s1,s2
    V=computeV(shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2)
    GG = G(t_i  , shat_ip1, shat_i  , t_ip1, m_ip1_2, 0d0) &
         & *G(t_im1, shat_i  , shat_im1, t_i  , m_i_2  , 0d0)
    if (GG.le.0d0 .or. V.eq.-99d99) then
!!$       write (*,*) 'No allowed range for s: smin=smax',GG,V
!!$       stop 1
       GG=0d0
       V=0d0
    endif
    sqrtGG=sqrt(GG)
    s1=shat_im1+shat_ip1+2d0/lambda(shat_i,t_i,0d0) * (4d0*V + sqrtGG)
    s2=shat_im1+shat_ip1+2d0/lambda(shat_i,t_i,0d0) * (4d0*V - sqrtGG)
    smin=min(s1,s2)
    smax=max(s1,s2)
  end subroutine sminmax
    subroutine rotxxx(p,q,prot)
      ! This subroutine performs the spacial rotation of a four-momentum.
      ! the momentum p is assumed to be given in the frame where the spacial
      ! component of q points the positive z-axis.  prot is the momentum p
      ! rotated to the frame where q is given.
      ! input:
      !       real    p(0:3)         : four-momentum p in q(1)=q(2)=0 frame
      !!       real    q(0:3)         : four-momentum q in the rotated frame
      ! output:
      !       real    prot(0:3)      : four-momentum p in the rotated frame
      implicit none
      real(kind=8),dimension(0:3) :: p,q
      real(kind=8),dimension(0:3) :: prot
      real(kind=8) :: qt2,qt,psgn,qq
      real(kind=8),parameter :: rZero=0d0,rOne=1d0
      prot(0) = p(0)
      qt2 = q(1)**2 + q(2)**2
      if ( qt2.eq.rZero ) then
         if ( q(3).eq.rZero ) then
            prot(1:3) = p(1:3)
         else
            psgn = sign(rOne,q(3))
            prot(1:3) = p(1:3)*psgn
         endif
      else
         qq = sqrt(qt2+q(3)**2)
         qt = sqrt(qt2)
         prot(1) = q(1)*q(3)/qq/qt*p(1) -q(2)/qt*p(2) +q(1)/qq*p(3)
         prot(2) = q(2)*q(3)/qq/qt*p(1) +q(1)/qt*p(2) +q(2)/qq*p(3)
         prot(3) =          -qt/qq*p(1)               +q(3)/qq*p(3)
      endif
    end subroutine rotxxx
    subroutine rotxxx_inv(p,q,prot)
      ! Same as rotxxx, but inverse. That is, first doing
      ! rotxxx(p,q,prot) and then rotxxx_inv(prot,q,p) should give you
      ! back the original p.
      implicit none
      real(kind=8),dimension(0:3),intent(in) :: p,q
      real(kind=8),dimension(0:3),intent(out) :: prot
      real(kind=8) :: qt2,qt,psgn,qq
      prot(0) = p(0)
      qt2 = q(1)**2 + q(2)**2
      if ( qt2.lt.vtiny ) then
         if ( q(3).eq.0d0 ) then
            prot(1:3)=p(1:3)
         else
            psgn = sign(1d0,q(3))
            prot(1:3)=p(1:3)*psgn
         endif
      else
         qq = sqrt(qt2+q(3)**2)
         qt = sqrt(qt2)
         prot(1) = q(1)*q(3)/qq/qt*p(1) +q(2)*q(3)/qq/qt*p(2) -  qt/qq*p(3)
         prot(2) =        -q(2)/qt*p(1) +        q(1)/qt*p(2)
         prot(3) =   qt*q(1)/qq/qt*p(1) +        q(2)/qq*p(2) +q(3)/qq*p(3)
      endif
    end subroutine rotxxx_inv
    subroutine rotz(p,phi,prot)
      implicit none
      real(kind=8),dimension(0:3) :: p,prot
      real(kind=8) :: phi
      prot(0)=p(0)
      prot(1)=p(1)*cos(phi)-p(2)*sin(phi)
      prot(2)=p(1)*sin(phi)+p(2)*cos(phi)
      prot(3)=p(3)
    end subroutine rotz
    subroutine gentcms2(pa,pb,pc,t,phi,m1,m2,p1,pr)
      ! Generates 4 momentum for particle p1, and remainder pr given the
      ! values t, and phi in the process pa+pb -> pr+p1.  Assumes incoming
      ! particles with momenta pa, pb and outgoing particles with mass
      ! m1,m2; t=(pb-p1)^2=(pa-pr)^2. Assumes that pa is a massless
      ! momentum; phi is the azimuthal angle between pr and pc in the pa+pb
      ! rest frame, with pa aligned with the positive z-axis.
      implicit none
      real(kind=8),intent(in) :: t,phi,m1,m2
      real(kind=8),intent(in),dimension(0:3) :: pa,pb,pc
      real(kind=8),intent(out),dimension(0:3) :: p1,pr
      real(kind=8) :: E_acms,p_acms,esum,esum2,ed,pp2,md2,pt,pt2,phi_off
      real(kind=8),dimension(0:3) :: ptot,pa_cms,ptotm,Pii,pc_cms,pc_rot,pii_rot
      real(kind=8),parameter :: tiny=1d-8
      ptot(0:3)=pa(0:3)+pb(0:3)
      ptotm(0)=ptot(0)
      ptotm(1:3)=-ptot(1:3)
      ! determine magnitude of Pii in cms frame (from dhelas routine mom2cx)
      ESUM2 = dot(ptot,ptot)
      if (esum2 .le. 0d0) then
         write (*,*) "error :: must be time-like momentum in gentcms2",esum2
         stop 1
      endif
      esum=sqrt(esum2)
      MD2=(M2-M1)*(M1+M2)
      ED=MD2/ESUM
      IF (M1*M2.EQ.0.d0) THEN
         PP2=0.25d0*(ESUM-ABS(ED))**2
      ELSE
         PP2=0.25d0*((MD2/ESUM)**2-2d0*(M1**2+M2**2)+ESUM**2)
         if(pp2.lt.0d0) then
            write(*,*) 'Error #12 in genps_fks.f: magnitude^2 of '/&
                 &/'3-momentum smaller than 0',pp2
            stop 1
         endif
      ENDIF
      call boostm(pa,ptotm,esum,pa_cms)
      E_acms = pa_cms(0)
      p_acms = sqrt(pa_cms(1)**2+pa_cms(2)**2+pa_cms(3)**2)

      ! determine the offset in phi; the frame in which phi is defined
      ! is in the ptot rest-frame, with pa_cms aligned with the z-axis,
      ! and pc having zero phi angle.
      call boostm(pc,ptotm,esum,pc_cms)
      call rotxxx_inv(pc_cms,pa_cms,pc_rot)
      if (pc_rot(1).ne.0d0) then
         phi_off=atan(pc_rot(2)/pc_rot(1))
      else
         phi_off=0d0
      endif
      if (pc_rot(1).lt.0d0) then
         phi_off=phi_off+pi
      endif

      ! define Pii in the frame where pa_cms is aligned with the positive z axis
      Pii(0) = MAX((ESUM+ED)*0.5d0,0.d0)
      if (esum+ed.le.0d0) then
         write (*,*) 'Error #15 in genps_fks.f: negative energy',esum,ed
         stop 1
      endif
      Pii(3) = -(m2**2-t-2d0*Pii(0)*E_acms)/(2d0*p_acms)
      pt2=pp2-Pii(3)**2
      if (pt2/esum2.lt.-tiny) then
         write (*,*) 'Error #16 in genps_fks.f: relative pt^2 smaller than 0',pt2
         stop 1
      elseif (pt2.lt.0d0) then
         pt2=0d0
      endif
      pt = sqrt(pt2)
      Pii(1) = pt*cos(phi+phi_off)
      Pii(2) = pt*sin(phi+phi_off)
      call rotxxx(Pii,pa_cms,Pii_rot)       !Rotate Pii to the pa_cms frame
      call boostm(Pii_rot,ptot,esum,Pii)    !boost back to lab fram
      p1(0:3)=ptot(0:3)-pii(0:3)
      pr(0:3)=pb(0:3)-p1(0:3)         !Return remainder of momentum
    end subroutine gentcms2
    real(kind=8) function getphifroms(si,shat_i,shat_im1,shat_ip1,t_i,V,sqrtGG,ran)
      ! Given s_i (invariant mass of p_i and p_i+1, it transforms it
      ! into phi_i. Note that there are two possibilities for phi: need
      ! to pick one at random.
      ! Based on eq.(11) of E.~Byckling and K.~Kajantie, ``Reductions of
      ! the phase-space integral in terms of simpler processes,''
      ! Phys. Rev. 187 (1969), 2008-2016, doi:10.1103/PhysRev.187.2008
      implicit none
      real(kind=8),intent(in) :: si,shat_i,shat_im1,shat_ip1,t_i,V,sqrtGG,ran
      real(kind=8) :: cosphi,x
      cosphi=((si-shat_im1-shat_ip1)*0.5d0*lambda(shat_i,t_i,0d0)-4d0*V)/sqrtGG
      x=ran
      if (x.gt.0.5d0) then
         getphifroms=acos(cosphi)
      else
         getphifroms=-acos(cosphi)+2d0*pi
      endif
    end function getphifroms


end module phase_space_gen23_mod
