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

  subroutine include_PDF_and_identical_procs(val,val_abs,pgl,ps,iint,ivec)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    type(psv),intent(in) :: ps
    real(kind=8),intent(inout),dimension(*) :: val,val_abs
    integer,intent(in) :: iint,ivec
    integer :: iproc,ip,iip,ip_start,ip_end
    real(kind=8) :: xmu_fac
    real(kind=8), dimension(-6:7,2) :: PDF
    if (include_pdf) then
       ! Include the PDFs
       xmu_fac=91.188d0 ! factorisation scale
       call PDF_eval(1,pgl%ipdgs(-6,1),ps%xbjrk(1),xmu_fac,PDF(-6,1))
       call PDF_eval(1,pgl%ipdgs(-6,2),ps%xbjrk(2),xmu_fac,PDF(-6,2))
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
       pgl%val_procs(1:pgl%iden_iproc(iproc),iproc,ivec)=val(iip)*pgl%idenCOandMAPfactor(1:pgl%iden_iproc(iproc),iproc)
       if (include_pdf) then
          do ip=1,pgl%iden_iproc(iproc)
             ! first incoming particle
             if (pgl%iden_processes(1,ip,iproc).eq.21) then
                pgl%val_procs(ip,iproc,ivec)=pgl%val_procs(ip,iproc,ivec)*PDF(0,1)
             elseif(pgl%iden_processes(1,ip,iproc).eq.22) then
                pgl%val_procs(ip,iproc,ivec)=pgl%val_procs(ip,iproc,ivec)*PDF(7,1)
             else
                pgl%val_procs(ip,iproc,ivec)=pgl%val_procs(ip,iproc,ivec)*PDF(pgl%iden_processes(1,ip,iproc),1)
             endif
             ! second incoming particle
             if (pgl%iden_processes(2,ip,iproc).eq.21) then
                pgl%val_procs(ip,iproc,ivec)=pgl%val_procs(ip,iproc,ivec)*PDF(0,2)
             elseif(pgl%iden_processes(2,ip,iproc).eq.22) then
                pgl%val_procs(ip,iproc,ivec)=pgl%val_procs(ip,iproc,ivec)*PDF(7,2)
             else
                pgl%val_procs(ip,iproc,ivec)=pgl%val_procs(ip,iproc,ivec)*PDF(pgl%iden_processes(2,ip,iproc),2)
             endif
          enddo
       endif
       val(iip)=sum(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc,ivec),dim=1)
       val_abs(iip)=sum(abs(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc,ivec)),dim=1)
    enddo
  end subroutine include_PDF_and_identical_procs

end module pdf_wrap
