
      function ran2()
!     Wrapper for the random numbers; needed for the NLO stuff
      use random_number_interface, only: ntuple
      implicit none
      double precision ran2,x,a,b
      integer jconfig
      a=0d0                     ! min allowed value for x
      b=1d0                     ! max allowed value for x
c$$$      jconfig=iconfig           ! integration channel (for off-set)
      jconfig=1
      call ntuple(x,a,b,jconfig)
      ran2=x
      return
      end function ran2


      subroutine ntuple(x,a,b,jconfig)
c-------------------------------------------------------
c     Front to ranmar which allows user to easily
c     choose the seed.
c------------------------------------------------------
      use random_number_interface, only: get_offset,get_base,
     &     get_moffset,rmarin,ranmar
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
c
c     Arguments
c
      double precision,intent(out) :: x
      double precision,intent(in) :: a,b
      integer,intent(in) :: jconfig
c
c     Local
c
      integer init, ioffset, joffset
      integer ij, kl
      integer*8 ij_work, kl_work, joffset_work, seed_work
      double precision span

c
c     Global
c
c-------
c     18/6/2012 tjs promoted to integer*8 to avoid overflow for iseed > 60K
c------
      integer*8       iseed
      integer*8       max_seed
      parameter      (max_seed=904866561)
      common /to_seed/iseed
c
c     Data
c
      data init /1/
      save ij, kl
c-----
c  Begin Code
c-----
      if (init .eq. 1) then
         init = 0
         call get_offset(ioffset)
cRF : always read the seed from the randinit file (this file is updated
c by amcatnlo_run_interface every time a run starts). This makes sure
c that the code does not need any remcompilation when only the seed is
c changed (useful for NLO gridpack mode).
         if (iseed .eq. 0) call get_base(iseed)
         if (iseed .lt. 0 .or. iseed .gt. max_seed) then
            write(*,*) 'Random seed must be between 0 and',max_seed,
     &           ':',iseed
            stop 1
         endif
c$$$         call get_base(iseed)
c
c     TJS 3/13/2008
c     Modified to allow for more sequences 
c     iseed can be between 0 and 30081*30081
c     before pattern repeats
c
c     TJS 12/3/2010
c     multipied iseed to give larger values more likely to make change
c     get offset for multiple runs of single process
c
c     TJS 18/6/2012
c     Updated to better divide iseed among ij and kl seeds
c     Note it may still be possible to get identical ij,kl for
c     different iseed, if have exactly compensating joffset, ioffset, jconfig
c
         call get_moffset(joffset)
         joffset_work = int(joffset,8) * 3157_8
         seed_work = iseed * 31300_8
         ij_work=1802_8+int(jconfig,8) + mod(seed_work,30081_8)
         kl_work=9373_8+(seed_work/30081_8)+int(ioffset,8)
     &        +joffset_work                          !Switched to 30081
                                                     !20/6/12 to avoid
                                                     !dupes in range
                                                     !30082-31328
         write(99,'(a,i0,a3,i0,a3,i0)') 'Using random seed offsets:'
     &        ,jconfig," , ",ioffset," , ",joffset_work
         write(99,*) ' with seed', iseed
         if (ij_work .gt. 31328_8)
     &        ij_work=1_8+modulo(ij_work-1_8,31328_8)
         if (kl_work .gt. 30081_8)
     &        kl_work=1_8+modulo(kl_work-1_8,30081_8)
         if (ij_work .lt. 0_8 .or. kl_work .lt. 0_8) then
            write(*,*) 'Random seed offsets produce negative seeds:',
     &           ij_work,kl_work
            stop 1
         endif
         ij=int(ij_work,kind(ij))
         kl=int(kl_work,kind(kl))
        call rmarin(ij,kl)         
      endif
      if (.not.ieee_is_finite(a) .or. .not.ieee_is_finite(b)
     &     .or. b .lt. a) then
         write(*,*) 'Invalid random-number interval:',a,b
         stop 1
      endif
      if (a .lt. 0d0) then
         if (b .gt. huge(b)+a) then
            write(*,*) 'Random-number interval is too wide:',a,b
            stop 1
         endif
      endif
      span=b-a
      call ranmar(x)
      do while (x .lt. 1d-16)
         call ranmar(x)
      enddo
      x = a+x*span
      if (.not.ieee_is_finite(x)) then
         write(*,*) 'Random-number interval mapping overflowed:',a,b
         stop 1
      endif
      end

      subroutine get_base(iseed)
c-------------------------------------------------------
c     Looks for file iproc.dat to offset random number gen
c------------------------------------------------------
      implicit none
c
c     Constants
c
      integer*8,intent(out) :: iseed
c
c     Local
c
      character*256 fname,line
      integer lun,level,ios,close_status
c-----
c  Begin Code
c-----

      iseed = 0
      fname = 'randinit'
      do level=1,4
         open(newunit=lun,file=trim(fname),status='old',action='read',
     &        iostat=ios)
         if (ios .eq. 0) then
            read(lun,'(a)',iostat=ios) line
            close(lun,iostat=close_status)
            if (ios .eq. 0 .and. close_status .eq. 0) then
               ios=index(line,'=')
               if (ios .gt. 0) line=line(ios+1:)
               read(line,*,iostat=ios) iseed
               if (ios .eq. 0) return
               iseed=0
            endif
         endif
         if (len_trim(fname) .gt. len(fname)-3) return
         fname='../'//trim(fname)
      enddo
c      write(*,*) 'No base found using iseed=0'
      end

      subroutine get_offset(iseed)
c-------------------------------------------------------
c     Looks for file iproc.dat to offset random number gen
c------------------------------------------------------
      implicit none
c
c     Constants
c
      integer,intent(out) :: iseed
c
c     Local
c
      integer lun,ios,close_status,value
c-----
c  Begin Code
c-----

      iseed=0
      open(newunit=lun,file='./iproc.dat',status='old',action='read',
     &     iostat=ios)
      if (ios .eq. 0) then
         read(lun,*,iostat=ios) value
         close(lun,iostat=close_status)
         if (ios .eq. 0 .and. close_status .eq. 0) then
            iseed=value
            return
         endif
      endif
      open(newunit=lun,file='../iproc.dat',status='old',action='read',
     &     iostat=ios)
      if (ios .eq. 0) then
         read(lun,*,iostat=ios) value
         close(lun,iostat=close_status)
         if (ios .eq. 0 .and. close_status .eq. 0) iseed=value
      endif
      end

      subroutine get_moffset(iseed)
c-------------------------------------------------------
c     Looks for file moffset.dat to offset random number gen
c------------------------------------------------------
      implicit none
c
c     Constants
c
      integer,intent(out) :: iseed
c
c     Local
c
      integer lun,ios,close_status,value
c-----
c  Begin Code
c-----

      iseed=0
      open(newunit=lun,file='./moffset.dat',status='old',action='read',
     &     iostat=ios)
      if (ios .eq. 0) then
         read(lun,*,iostat=ios) value
         close(lun,iostat=close_status)
         if (ios .eq. 0 .and. close_status .eq. 0) then
            iseed=value
            write(99,*) "Got moffset",iseed
         endif
      endif
      end

      subroutine ranmar(rvec)
*     -----------------
* universal random number generator proposed by marsaglia and zaman
* in report fsu-scri-87-50
* in this version rvec is a double precision variable.
      use random_number_state, only: ranmar_initialized
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit real*8(a-h,o-z)
      double precision,intent(out) :: rvec
      common/ raset1 / ranu(97),ranc,rancd,rancm
      common/ raset2 / iranmr,jranmr
      save /raset1/,/raset2/
      rvec=0d0
      if (.not.ranmar_initialized) then
         write(*,*) 'RANMAR used before initialization'
         stop 1
      endif
      if (iranmr .lt. 1 .or. iranmr .gt. 97 .or.
     &     jranmr .lt. 1 .or. jranmr .gt. 97) then
         write(*,*) 'RANMAR state contains invalid indices:',
     &        iranmr,jranmr
         stop 1
      endif
      uni = ranu(iranmr) - ranu(jranmr)
      if(uni .lt. 0d0) uni = uni + 1d0
      ranu(iranmr) = uni
      iranmr = iranmr - 1
      jranmr = jranmr - 1
      if(iranmr .eq. 0) iranmr = 97
      if(jranmr .eq. 0) jranmr = 97
      ranc = ranc - rancd
      if(ranc .lt. 0d0) ranc = ranc + rancm
      uni = uni - ranc
      if(uni .lt. 0d0) uni = uni + 1d0
      rvec = uni
      if (.not.ieee_is_finite(rvec)) then
         write(*,*) 'RANMAR produced an invalid value:',rvec
         stop 1
      endif
      if (rvec .lt. 0d0 .or. rvec .ge. 1d0) then
         write(*,*) 'RANMAR produced an invalid value:',rvec
         stop 1
      endif
      end
 
      subroutine rmarin(ij,kl)
*     -----------------
* initializing routine for ranmar, must be called before generating
* any pseudorandom numbers with ranmar. the input values should be in
* the ranges 0<=ij<=31328 ; 0<=kl<=30081
      use random_number_state, only: ranmar_initialized
      implicit real*8(a-h,o-z)
      integer,intent(in) :: ij,kl
      common/ raset1 / ranu(97),ranc,rancd,rancm
      common/ raset2 / iranmr,jranmr
      save /raset1/,/raset2/
      ranmar_initialized=.false.
* this shows correspondence between the simplified input seeds ij, kl
* and the original marsaglia-zaman seeds i,j,k,l.
* to get the standard values in the marsaglia-zaman paper (i=12,j=34
* k=56,l=78) put ij=1802, kl=9373
      write(99,*) "Ranmar initialization seeds",ij,kl
c
c    18/6/2012 TJS  Added check to ensure ij and kl are in range
c      
      if (ij .lt. 0 .or. ij .gt. 31328 .or.
     $     kl .lt. 0 .or. kl .gt. 30081) then
         if (ij .lt. 0 .or. ij .gt. 31328) then
            write(*,*) 'Bad initialization value of ij in rmarin ', ij
         elseif (kl .lt. 0 .or. kl .gt. 30081) then
            write(*,*) 'Bad initialization value of kl in rmarin ', kl
         endif
         stop 1
      endif

      i = mod( ij/177 , 177 ) + 2
      j = mod( ij     , 177 ) + 2
      k = mod( kl/169 , 178 ) + 1
      l = mod( kl     , 169 )
      do 300 ii = 1 , 97
        s =  0d0
        t = .5d0
        do 200 jj = 1 , 24
          m = mod( mod(i*j,179)*k , 179 )
          i = j
          j = k
          k = m
          l = mod( 53*l+1 , 169 )
          if(mod(l*m,64) .ge. 32) s = s + t
          t = .5d0*t
  200   continue
        ranu(ii) = s
  300 continue
      ranc  =   362436d0 / 16777216d0
      rancd =  7654321d0 / 16777216d0
      rancm = 16777213d0 / 16777216d0
      iranmr = 97
      jranmr = 33
      ranmar_initialized=.true.
      end
