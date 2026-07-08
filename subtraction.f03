module subtraction
  use handling_processes
  use particles
  implicit none
  integer :: n
  type dipole
     integer,dimension(:),allocatable :: process_r
     integer,dimension(3) :: dip_ijk,dip_ijk_f
     integer,dimension(2) :: dip_r_ijk,dip_r_ijk_f
     integer :: dipole_type=0 ! 0:II, 1:IF, 2:FI, 3:FF
   contains
     final :: finalize_dipole
  end type dipole
  integer :: ndip
  type(dipole),dimension(:),allocatable :: dl
  private
  public :: initialise_subtraction
contains
  subroutine initialise_subtraction(igroup,iamp)
    implicit none
    integer,intent(in) :: igroup,iamp
    integer :: ipart,is_dipole,ipart_l,ipart_r,idip
    n=pgl(igroup)%next
    ! First pass: just count how many dipoles we need so we can
    ! allocated the right size dl
    ndip=0
    do ipart=3,pgl(igroup)%next
       call is_valid_dipole(ipart,pgl(igroup)%processes(:,iamp),pgl(igroup)%phase_space_orders(:),is_dipole,ipart_l,ipart_r)
       ndip=ndip+popcnt(is_dipole)
    enddo
    write (*,*) 'Need',ndip,'dipoles'
    allocate(dl(ndip))

    ! Second pass: now we really fill the appropriate information
    idip=0
    do ipart=3,pgl(igroup)%next
       call is_valid_dipole(ipart,pgl(igroup)%processes(:,iamp),pgl(igroup)%phase_space_orders(:),is_dipole,ipart_l,ipart_r)
       if (is_dipole.eq.0) cycle
       if (btest(is_dipole,0)) then
          ! valid dipole with particle on the left as emitter
          idip=idip+1
          call fill_dipole(idip,pgl(igroup)%processes(1:n,iamp),ipart_l,ipart,ipart_r)
       endif
       if (btest(is_dipole,1)) then
          ! valid dipole with particle on the right as emitter
          idip=idip+1
          call fill_dipole(idip,pgl(igroup)%processes(1:n,iamp),ipart_r,ipart,ipart_l)
       endif
    enddo
    call print_dipoles(pgl(igroup)%processes(:,iamp))
    do idip=1,ndip
       call finalize_dipole(dl(idip))
    enddo
    deallocate(dl)
  end subroutine initialise_subtraction
  subroutine print_dipoles(process)
    implicit none
    integer,dimension(*),intent(in) :: process
    integer :: idip
    write (*,*) 'process',process(1:n)
    do idip=1,ndip
       write (*,*) '------------------'
       write (*,*) 'dipole',idip
       write (*,*) 'i,j,k',dl(idip)%dip_ijk
       write (*,*) 'i,j,k',dl(idip)%dip_ijk_f
       write (*,*) 'process reduced',dl(idip)%process_r
    enddo
    write (*,*) '------------------'
    write (*,*) ''
    write (*,*) ''
  end subroutine print_dipoles
  subroutine fill_dipole(idip,process,dip_i,dip_j,dip_k)
    implicit none
    integer,dimension(*),intent(in) :: process
    integer,intent(in) :: idip,dip_i,dip_j,dip_k
    integer :: ipart,i
    dl(idip)%dip_ijk(1:3)=[dip_i,dip_j,dip_k]
    dl(idip)%dip_ijk_f(1:3)=[process(dip_i),process(dip_j),process(dip_k)]
    if (dip_i.gt.2) dl(idip)%dipole_type=ibset(dl(idip)%dipole_type,0)
    if (dip_k.gt.2) dl(idip)%dipole_type=ibset(dl(idip)%dipole_type,1)
    ! reduced process and dipole info
    dl(idip)%dip_r_ijk(1)=min(dip_i,dip_j)
    if (dip_j .lt. dip_k) then
       dl(idip)%dip_r_ijk(2)=dip_k-1
    else
       dl(idip)%dip_r_ijk(2)=dip_k
    endif
    if (phys_model%is_gluon(process(dip_j))) then
       dl(idip)%dip_r_ijk_f(1)=dl(idip)%dip_ijk_f(1)
    elseif (phys_model%is_gluon(process(dip_i))) then
       if (btest(dl(idip)%dipole_type,0)) then
          write (*,*) 'error in dipoles: emitter is a final-state gluon and '// &
               'emitted is a quark'
          write (*,*) dl(idip)%dip_ijk
          write (*,*) dl(idip)%dipole_type
          write (*,*) dl(idip)%dip_ijk_f
          stop 1
       endif
       dl(idip)%dip_r_ijk_f(1)=phys_model%get_antipart(dl(idip)%dip_ijk_f(2))
    else
       dl(idip)%dip_r_ijk_f(1)=21
    endif
    dl(idip)%dip_r_ijk_f(2)=dl(idip)%dip_ijk_f(3)
    allocate(dl(idip)%process_r(n-1))
    i=0
    do ipart=1,n
       if (ipart.eq.dip_j) cycle
       i=i+1
       if (ipart.eq.dip_i) then
          dl(idip)%process_r(i)=dl(idip)%dip_r_ijk_f(1)
       elseif(ipart.eq.dip_k) then
          dl(idip)%process_r(i)=dl(idip)%dip_r_ijk_f(2)
       else
          dl(idip)%process_r(i)=process(ipart)
       endif
    enddo
  end subroutine fill_dipole
  subroutine is_valid_dipole(ipart,process,order,is_dipole,ipart_l,ipart_r)
    ! Checks if the two particles next to ipart in the colour order
    ! form a valid dipole that could have radiated particle ipart
    implicit none
    integer,intent(in) :: ipart
    integer,dimension(n),intent(in) :: process,order
    integer,intent(out) :: is_dipole,ipart_l,ipart_r
    integer :: i
    is_dipole=0
    do i=1,n
       if (order(i).eq.ipart) exit
    enddo
    ipart_l=order(mod(n+i-1,n)) ! left of ipart in colour order
    ipart_r=order(mod(i,n)+1)   ! right of ipart in colour order
    
    if (phys_model%get_mass(process(ipart)).ne.0d0) return
    if (phys_model%get_colour_rep(process(ipart)).eq.1) return
    if (phys_model%get_colour_rep(process(ipart_l)).eq.1) return
    if (phys_model%get_colour_rep(process(ipart_r)).eq.1) return

    if (phys_model%is_gluon(process(ipart))) then
       is_dipole=3 ! both left and right can be emitters
    else ! must be a quark
       if (ipart_l.gt.2) then
          if (phys_model%get_antipart(process(ipart)).eq.process(ipart_l)) is_dipole=ibset(is_dipole,0)
       else
          if (process(ipart).eq.process(ipart_l)) is_dipole=ibset(is_dipole,0)
          if (phys_model%is_gluon(process(ipart_l))) is_dipole=ibset(is_dipole,0)
       endif
       if (ipart_r.gt.2) then
          if (phys_model%get_antipart(process(ipart)).eq.process(ipart_r)) is_dipole=ibset(is_dipole,1)
       else
          if (process(ipart).eq.process(ipart_r)) is_dipole=ibset(is_dipole,1)
          if (phys_model%is_gluon(process(ipart_r))) is_dipole=ibset(is_dipole,1)
       endif
    endif
  end subroutine is_valid_dipole
  subroutine finalize_dipole(di)
    type(dipole) :: di
    if (allocated(di%process_r)) deallocate(di%process_r)
  end subroutine finalize_dipole
end module subtraction
