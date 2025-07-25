module pdf_help
  use common
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

  subroutine include_PDF_and_identical_procs(val,val_abs,pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    real(kind=8),intent(inout),dimension(*) :: val,val_abs
    integer :: iproc,ip
    real(kind=8) :: xmu_fac
    real(kind=8), dimension(-6:7,2) :: PDF
    if (include_pdf) then
       ! Include the PDFs
       xmu_fac=91.188d0 ! factorisation scale
       call PDF_eval(1,pgl%ipdgs(-6,1),pgl%phase_space%xbjrk(1),xmu_fac,PDF(-6,1))
       call PDF_eval(1,pgl%ipdgs(-6,2),pgl%phase_space%xbjrk(2),xmu_fac,PDF(-6,2))
    endif
    do iproc=1,pgl%nproc
       pgl%val_procs(1:pgl%iden_iproc(iproc),iproc)=val(iproc)*pgl%idenCOandMAPfactor(1:pgl%iden_iproc(iproc),iproc)
       if (include_pdf) then
          do ip=1,pgl%iden_iproc(iproc)
             ! first incoming particle
             if (pgl%iden_processes(1,ip,iproc).eq.21) then
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(0,1)
             elseif(pgl%iden_processes(1,ip,iproc).eq.22) then
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(7,1)
             else
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(pgl%iden_processes(1,ip,iproc),1)
             endif
             ! second incoming particle
             if (pgl%iden_processes(2,ip,iproc).eq.21) then
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(0,2)
             elseif(pgl%iden_processes(2,ip,iproc).eq.22) then
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(7,2)
             else
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(pgl%iden_processes(2,ip,iproc),2)
             endif
          enddo
       endif
       val(iproc)=sum(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc))
       val_abs(iproc)=sum(abs(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc)))
    enddo
  end subroutine include_PDF_and_identical_procs

end module pdf_help
