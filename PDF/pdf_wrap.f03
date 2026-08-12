module pdf_wrap
  use handling_processes
contains
    subroutine set_ipdgs_for_PDF(pgl)
    ! determines for which flavours the PDFs should be evolved
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: iflav,iproc
    pgl%ipdgs(-6:7,1:2)=.false.
    do iflav=-6,7
       do iproc=1,pgl%nproc
          if (iflav.eq.0) then    ! gluon
             if (any(pgl%iden_processes(1,1:pgl%iden_iproc(iproc),iproc).eq.21)) pgl%ipdgs(iflav,1)=.true.
             if (any(pgl%iden_processes(2,1:pgl%iden_iproc(iproc),iproc).eq.21)) pgl%ipdgs(iflav,2)=.true.
          elseif(iflav.eq.7) then ! photon
             if (any(pgl%iden_processes(1,1:pgl%iden_iproc(iproc),iproc).eq.22)) pgl%ipdgs(iflav,1)=.true.
             if (any(pgl%iden_processes(2,1:pgl%iden_iproc(iproc),iproc).eq.22)) pgl%ipdgs(iflav,2)=.true.
          else                    ! quarks and anti-quarks
             if (any(pgl%iden_processes(1,1:pgl%iden_iproc(iproc),iproc).eq.iflav)) pgl%ipdgs(iflav,1)=.true.
             if (any(pgl%iden_processes(2,1:pgl%iden_iproc(iproc),iproc).eq.iflav)) pgl%ipdgs(iflav,2)=.true.
          endif
       enddo
    enddo
  end subroutine set_ipdgs_for_PDF

  subroutine include_PDF_and_identical_procs(val,val_abs,pgl,iint)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    real(kind=8),intent(inout),dimension(*) :: val,val_abs
    integer,intent(in) :: iint
    integer :: iproc,ip,iip,ip_start,ip_end,iset,i
    real(kind=8), dimension(-6:7,2) :: PDF=0
    if (include_pdf) then
       ! Include the PDFs
       if (use_lhapdf) then
! The following commented-out section would evolve only the necessary partons
!!$          call getnset(iset)
!!$          do i=-6,7
!!$             if (.not.pgl%ipdgs(i,1)) cycle
!!$             call evolvepartm(iset,i,pgl%ps(1)%xbjrk(1),scale_fac,PDF(i,1))
!!$          enddo
!!$          do i=-6,7
!!$             if (.not.pgl%ipdgs(i,2)) cycle
!!$             call evolvepartm(iset,i,pgl%ps(1)%xbjrk(2),scale_fac,PDF(i,2))
!!$          enddo
          call evolvePDF(pgl%ps(1)%xbjrk(1),scale_fac,PDF(-6,1))
          call evolvePDF(pgl%ps(1)%xbjrk(2),scale_fac,PDF(-6,2))
          PDF(:,1)=PDF(:,1)/pgl%ps(1)%xbjrk(1)
          PDF(:,2)=PDF(:,2)/pgl%ps(1)%xbjrk(2)
       else
          call PDF_eval(1,pgl%ipdgs(-6,1),pgl%ps(1)%xbjrk(1),scale_fac,PDF(-6,1))
          call PDF_eval(1,pgl%ipdgs(-6,2),pgl%ps(1)%xbjrk(2),scale_fac,PDF(-6,2))
       endif
    endif
    if (iint.gt.0) then
       ip_start=iint
       ip_end=iint
    else
       ip_start=1
       ip_end=pgl%nproc
    endif
    iip=0
    do iproc=ip_start,ip_end
       iip=iip+1
       pgl%alias_factors(1:pgl%iden_iproc(iproc),iproc)=&
            pgl%idenCOandMAPfactor(1:pgl%iden_iproc(iproc),iproc)
       if (include_pdf) then
          do ip=1,pgl%iden_iproc(iproc)
             ! first incoming particle
             if (pgl%iden_processes(1,ip,iproc).eq.21) then
                pgl%alias_factors(ip,iproc)=pgl%alias_factors(ip,iproc)*PDF(0,1)
             elseif(pgl%iden_processes(1,ip,iproc).eq.22) then
                pgl%alias_factors(ip,iproc)=pgl%alias_factors(ip,iproc)*PDF(7,1)
             else
                pgl%alias_factors(ip,iproc)=pgl%alias_factors(ip,iproc)*PDF(pgl%iden_processes(1,ip,iproc),1)
             endif
             ! second incoming particle
             if (pgl%iden_processes(2,ip,iproc).eq.21) then
                pgl%alias_factors(ip,iproc)=pgl%alias_factors(ip,iproc)*PDF(0,2)
             elseif(pgl%iden_processes(2,ip,iproc).eq.22) then
                pgl%alias_factors(ip,iproc)=pgl%alias_factors(ip,iproc)*PDF(7,2)
             else
                pgl%alias_factors(ip,iproc)=pgl%alias_factors(ip,iproc)*PDF(pgl%iden_processes(2,ip,iproc),2)
             endif
          enddo
       endif
       pgl%val_procs(1:pgl%iden_iproc(iproc),iproc)=&
            val(iip)*pgl%alias_factors(1:pgl%iden_iproc(iproc),iproc)
       pgl%val_procs_abs(1:pgl%iden_iproc(iproc),iproc)=&
            val_abs(iip)*abs(pgl%alias_factors(1:pgl%iden_iproc(iproc),iproc))
       val(iip)=sum(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc))
       val_abs(iip)=sum(pgl%val_procs_abs(1:pgl%iden_iproc(iproc),iproc))
    enddo
  end subroutine include_PDF_and_identical_procs

end module pdf_wrap
