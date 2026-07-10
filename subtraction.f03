module subtraction
  use handling_processes
  use particles
  implicit none
  integer :: n
  private
  public :: initialise_subtraction
contains
  subroutine initialise_subtraction(igroup,iamp)
    implicit none
    integer,intent(in) :: iamp,igroup
    integer :: ipart,is_dipole,ipart_l,ipart_r,idip
    n=pgl(igroup)%next
    if (.not.allocated(pgl(igroup)%dpl)) allocate(pgl(igroup)%dpl(pgl(igroup)%nproc))
    if (allocated(pgl(igroup)%dpl(iamp)%dl)) then
       call finalize_dipole_set(pgl(igroup)%dpl(iamp))
    endif
    ! First pass: just count how many dipoles we need so we can
    ! allocate the right size dl
    pgl(igroup)%dpl(iamp)%ndip=0
    do ipart=3,pgl(igroup)%next
       call is_valid_dipole(ipart,pgl(igroup)%processes(:,iamp),pgl(igroup)%phase_space_orders(:),is_dipole,ipart_l,ipart_r)
       pgl(igroup)%dpl(iamp)%ndip=pgl(igroup)%dpl(iamp)%ndip+popcnt(is_dipole)
    enddo
    write (*,*) 'Need',pgl(igroup)%dpl(iamp)%ndip,'dipoles'
    allocate(pgl(igroup)%dpl(iamp)%dl(pgl(igroup)%dpl(iamp)%ndip))

    ! Second pass: now we really fill the appropriate information
    idip=0
    do ipart=3,pgl(igroup)%next
       call is_valid_dipole(ipart,pgl(igroup)%processes(:,iamp),pgl(igroup)%phase_space_orders(:),is_dipole,ipart_l,ipart_r)
       if (is_dipole.eq.0) cycle
       if (btest(is_dipole,0)) then
          ! valid dipole with particle on the left as emitter
          idip=idip+1
          call fill_dipole(pgl(igroup)%dpl(iamp)%dl(idip),pgl(igroup)%processes(1:n,iamp),ipart_l,ipart,ipart_r,.true.)
       endif
       if (btest(is_dipole,1)) then
          ! valid dipole with particle on the right as emitter
          idip=idip+1
          call fill_dipole(pgl(igroup)%dpl(iamp)%dl(idip),pgl(igroup)%processes(1:n,iamp),ipart_r,ipart,ipart_l,.false.)
       endif
    enddo
    do idip=1,pgl(igroup)%dpl(iamp)%ndip
       allocate(pgl(igroup)%dpl(iamp)%dl(idip)%reduced_color_order(n-1))
       call build_reduced_color_order(pgl(igroup)%color_orders(1:n,iamp), &
            pgl(igroup)%dpl(iamp)%dl(idip)%dip_ijk(2),pgl(igroup)%dpl(iamp)%dl(idip)%process_r, &
            pgl(igroup)%dpl(iamp)%dl(idip)%reduced_color_order)
       pgl(igroup)%dpl(iamp)%dl(idip)%col_fac=lc_colour_factor(pgl(igroup)%dpl(iamp)%dl(idip)%process_r)
       call pgl(igroup)%dpl(iamp)%dl(idip)%amp%init(1,n-1,1,pgl(igroup)%dpl(iamp)%dl(idip)%process_r,&
            pgl(igroup)%spin(0:3,pgl(igroup)%dpl(iamp)%dl(idip)%dip_map(1:n-1)), &
            pgl(igroup)%dpl(iamp)%dl(idip)%reduced_color_order,&
            phys_model)
    enddo
    call print_dipoles(pgl(igroup)%processes(:,iamp),pgl(igroup)%color_orders(:,iamp),pgl(igroup)%dpl(iamp)%dl)
  end subroutine initialise_subtraction
  subroutine print_dipoles(process,order,dips)
    implicit none
    integer,dimension(*),intent(in) :: process,order
    type(dipole),dimension(:),intent(in) :: dips
    integer :: idip
    write (*,*) 'process',process(1:n)
    write (*,*) 'color-order',order(1:n)
    do idip=1,size(dips)
       write (*,*) '------------------'
       write (*,*) 'dipole',idip
       write (*,*) 'i,j,k',dips(idip)%dip_ijk
       write (*,*) 'i,j,k',dips(idip)%dip_ijk_f
       write (*,*) 'process reduced',dips(idip)%process_r
       write (*,*) 'color order reduced',dips(idip)%reduced_color_order
    enddo
    write (*,*) '------------------'
    write (*,*) ''
    write (*,*) ''
  end subroutine print_dipoles
  subroutine fill_dipole(dip,process,dip_i,dip_j,dip_k,reverse)
    implicit none
    integer,dimension(*),intent(in) :: process
    type(dipole),intent(inout) :: dip
    integer,intent(in) :: dip_i,dip_j,dip_k
    logical,intent(in) :: reverse
    integer :: ipart,i
    dip%dip_ijk(1:3)=[dip_i,dip_j,dip_k]
    dip%dip_ijk_f(1:3)=[process(dip_i),process(dip_j),process(dip_k)]
    if (dip_i.gt.2) dip%dipole_type=ibset(dip%dipole_type,0)
    if (dip_k.gt.2) dip%dipole_type=ibset(dip%dipole_type,1)
    ! reduced process and dipole info
    if (dip_j .lt. dip_i) then
       dip%dip_r_ijk(1)=dip_i-1
    else
       dip%dip_r_ijk(1)=dip_i
    endif
    if (dip_j .lt. dip_k) then
       dip%dip_r_ijk(2)=dip_k-1
    else
       dip%dip_r_ijk(2)=dip_k
    endif
    if (phys_model%is_gluon(process(dip_j))) then
       dip%dip_r_ijk_f(1)=dip%dip_ijk_f(1)
    elseif (phys_model%is_gluon(process(dip_i))) then
       if (btest(dip%dipole_type,0)) then
          write (*,*) 'error in dipoles: emitter is a final-state gluon and '// &
               'emitted is a quark'
          write (*,*) dip%dip_ijk
          write (*,*) dip%dipole_type
          write (*,*) dip%dip_ijk_f
          stop 1
       endif
       dip%dip_r_ijk_f(1)=phys_model%get_antipart(dip%dip_ijk_f(2))
    else
       dip%dip_r_ijk_f(1)=combined_gluon_type(dip_i,process(dip_i),dip_j,process(dip_j),reverse)
    endif
    if (phys_model%is_gluon(dip%dip_r_ijk_f(1))) dip%lc_weight=0.5d0
    dip%dip_r_ijk_f(2)=dip%dip_ijk_f(3)
    allocate(dip%process_r(n-1))
    allocate(dip%dip_map(n-1))
    i=0
    do ipart=1,n
       if (ipart.eq.dip_j) cycle
       i=i+1
       dip%dip_map(i)=ipart
       if (ipart.eq.dip_i) then
          dip%process_r(i)=dip%dip_r_ijk_f(1)
       elseif(ipart.eq.dip_k) then
          dip%process_r(i)=dip%dip_r_ijk_f(2)
       else
          dip%process_r(i)=process(ipart)
       endif
    enddo
  end subroutine fill_dipole
  integer function lc_colour_factor(process)
    implicit none
    integer,dimension(:),intent(in) :: process
    integer :: i,ifac
    real(kind=8) :: fac
    fac=0d0
    do i=1,size(process)
       if (phys_model%is_gluon(process(i))) then
          fac=fac+1d0
       elseif (phys_model%is_quark(process(i)) .or. phys_model%is_antiquark(process(i))) then
          fac=fac+0.5d0
       endif
    enddo
    ifac=nint(fac)
    if (dble(ifac).ne.fac) then
       write (*,*) 'There is some issue with the reduced LC colour factor computation: ',ifac,fac
       stop 1
    endif
    lc_colour_factor=3**ifac
  end function lc_colour_factor

  integer function combined_gluon_type(dip_i,part_i,dip_j,part_j,reverse)
    implicit none
    integer,intent(in) :: dip_i,part_i,dip_j,part_j
    logical,intent(in) :: reverse
    if (((dip_i.le.2 .and. phys_model%is_quark(part_i)) .or. &
         (dip_i.gt.2 .and. phys_model%is_antiquark(part_i))) .and. &
        ((dip_j.le.2 .and. phys_model%is_antiquark(part_j)) .or. &
         (dip_j.gt.2 .and. phys_model%is_quark(part_j)))) then
       if (reverse) then
          combined_gluon_type=21
       else
          combined_gluon_type=99
       endif
    elseif (((dip_i.le.2 .and. phys_model%is_antiquark(part_i)) .or. &
             (dip_i.gt.2 .and. phys_model%is_quark(part_i))) .and. &
            ((dip_j.le.2 .and. phys_model%is_quark(part_j)) .or. &
             (dip_j.gt.2 .and. phys_model%is_antiquark(part_j)))) then
       if (reverse) then
          combined_gluon_type=99
       else
          combined_gluon_type=21
       endif
    else
       write (*,*) 'ERROR: cannot infer combined gluon type from dipole pair'
       write (*,*) part_i,part_j
       stop 1
    endif
  end function combined_gluon_type
  subroutine build_reduced_color_order(parent_order,removed_pos,process,order)
    implicit none
    integer,dimension(:),intent(in) :: parent_order
    integer,intent(in) :: removed_pos
    integer,dimension(:),intent(in) :: process
    integer,dimension(:),intent(out) :: order
    integer :: i,ipos,insert
    logical :: valid
    ipos=0
    do i=1,n
       if (parent_order(i).eq.removed_pos) cycle
       ipos=ipos+1
       if (parent_order(i).gt.removed_pos) then
          order(ipos)=parent_order(i)-1
       else
          order(ipos)=parent_order(i)
       endif
    enddo
    do i=1,n-1
       if ((order(i).le.2 .and. phys_model%is_antiquark(process(order(i)))) .or. &
            (order(i).gt.2 .and. phys_model%is_quark(process(order(i))))) then
          ! found quark to start colour order with
          order=[order(i:),order(:i-1)]
          exit
       endif
    enddo
    do i=n-1,1,-1
       if ((order(i).le.2 .and. phys_model%is_quark(process(order(i)))) .or. &
            (order(i).gt.2 .and. phys_model%is_antiquark(process(order(i))))) then
          ! found the last antiquark
          insert=i
          exit
       endif
    enddo
    i=n-1
    do
       if (phys_model%is_singlet(process(order(i)))) then
          if (i.gt.insert) then
             order=[order(1:insert-1),order(i),order(insert:i-1),order(i+1:n-1)]
             insert=insert+1
             i=i+1
          elseif (i.lt.insert) then
             order=[order(1:i-1),order(i+1:insert-1),order(i),order(insert:n-1)]
             insert=insert-1 ! new position
          endif
       endif
       i=i-1
       if (i.eq.0) exit
    enddo
  end subroutine build_reduced_color_order
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
end module subtraction
