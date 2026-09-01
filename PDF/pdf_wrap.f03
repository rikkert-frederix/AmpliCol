module pdf_wrap
  use handling_processes
  use pdf_internal_interface, only: PDF_eval,PDF_initialise
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  real(kind=8),parameter :: pdf_weight_limit=0.25d0*huge(1d0)**0.25d0
  interface
     subroutine InitPDFsetbyname(name)
       implicit none
       character(len=*),intent(in) :: name
     end subroutine InitPDFsetbyname
     subroutine initPDF(member)
       implicit none
       integer,intent(in) :: member
     end subroutine initPDF
     subroutine setlhaparm(setting)
       implicit none
       character(len=*),intent(in) :: setting
     end subroutine setlhaparm
     function alphaspdf(scale) result(value)
       implicit none
       real(kind=8),intent(in) :: scale
       real(kind=8) :: value
     end function alphaspdf
     subroutine evolvePDF(x,scale,pdfs)
       implicit none
       real(kind=8),intent(in) :: x,scale
       real(kind=8),intent(out) :: pdfs(-6:6)
     end subroutine evolvePDF
  end interface
contains
  subroutine evaluate_pdf_table(x,scale,all_pdf,status)
    ! Evaluate every PDF index used by AmpliCol in one backend call.  This is
    ! important for integrated convolutions, where many splitting histories
    ! query different flavours at the same x and scale.
    implicit none
    real(kind=8),intent(in) :: x,scale
    real(kind=8),intent(out) :: all_pdf(-6:7)
    integer,intent(out),optional :: status
    logical :: requested(-6:7)

    all_pdf=0d0
    if (present(status)) status=0
    if (.not.include_pdf) return
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(scale)) then
       if (present(status)) status=-20
       return
    endif
    if (x.le.0d0 .or. x.gt.1d0 .or. scale.lt.sqrt(tiny(1d0)) .or. &
         scale.gt.pdf_weight_limit) then
       if (present(status)) status=-20
       return
    endif
    if (x.eq.1d0) return
    if (x.lt.sqrt(tiny(1d0))) then
       if (present(status)) status=-20
       return
    endif
    if (use_lhapdf) then
       call evolvePDF(x,scale,all_pdf(-6:6))
       if (.not.all(ieee_is_finite(all_pdf))) then
          all_pdf=0d0
          if (present(status)) status=-20
          return
       endif
       if (any(abs(all_pdf).gt.pdf_weight_limit*x)) then
          all_pdf=0d0
          if (present(status)) status=-20
          return
       endif
       all_pdf=all_pdf/x
    else
       requested=.true.
       call PDF_eval(1,requested,x,scale,all_pdf(-6:7))
    endif
    if (.not.all(ieee_is_finite(all_pdf))) then
       all_pdf=0d0
       if (present(status)) status=-20
       return
    endif
    if (any(abs(all_pdf).gt.pdf_weight_limit)) then
       all_pdf=0d0
       if (present(status)) status=-20
       return
    endif
    where (abs(all_pdf).lt.tiny(1d0)) all_pdf=0d0
  end subroutine evaluate_pdf_table

  pure real(kind=8) function pdf_table_flavour(all_pdf,flavour) result(value)
    implicit none
    real(kind=8),intent(in) :: all_pdf(-6:7)
    integer,intent(in) :: flavour
    integer :: index

    if (flavour.eq.21 .or. flavour.eq.99) then
       index=0
    elseif (flavour.eq.22) then
       index=7
    elseif ((flavour.ge.1 .and. flavour.le.6) .or. &
         (flavour.le.-1 .and. flavour.ge.-6)) then
       index=flavour
    else
       value=0d0
       return
    endif
    value=all_pdf(index)
  end function pdf_table_flavour

  subroutine evaluate_pdf_flavour(flavour,x,scale,value,status)
    ! Return f_flavour(x,scale), using PDG 21 for a gluon and 22 for a
    ! photon.  LHAPDF's evolvePDF returns x*f and is converted here to f.
    implicit none
    integer,intent(in) :: flavour
    real(kind=8),intent(in) :: x,scale
    real(kind=8),intent(out) :: value
    integer,intent(out),optional :: status
    real(kind=8) :: all_pdf(-6:7)
    integer :: pdf_status

    if (.not.pdf_flavour_is_supported(flavour)) then
       value=0d0
       if (present(status)) status=-4
       return
    endif
    call evaluate_pdf_table(x,scale,all_pdf,pdf_status)
    if (present(status)) status=pdf_status
    value=pdf_table_flavour(all_pdf,flavour)
  end subroutine evaluate_pdf_flavour

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

  subroutine include_PDF_and_identical_procs(val,val_abs,pgl,iint,status)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    real(kind=8),intent(inout),dimension(*) :: val,val_abs
    integer,intent(in) :: iint
    integer,intent(out),optional :: status
    integer :: iproc,ip,iip,ip_start,ip_end,pdf_status
    real(kind=8), dimension(-6:7,2) :: PDF=0
    real(kind=8) :: pdf_first,pdf_second
    if (present(status)) status=0
    if (iint.gt.0) then
       ip_start=iint
       ip_end=iint
    else
       ip_start=1
       ip_end=pgl%nproc
    endif
    if (include_pdf) then
       call evaluate_pdf_table(pgl%ps(1)%xbjrk(1),scale_fac,PDF(:,1),pdf_status)
       if (pdf_status.ne.0) then
          call pdf_failure(pdf_status)
          return
       endif
       call evaluate_pdf_table(pgl%ps(1)%xbjrk(2),scale_fac,PDF(:,2),pdf_status)
       if (pdf_status.ne.0) then
          call pdf_failure(pdf_status)
          return
       endif
    endif
    iip=0
    do iproc=ip_start,ip_end
       iip=iip+1
       if (include_pdf) then
          if (any(.not.pdf_flavour_is_supported(&
               pgl%iden_processes(1:2,1:pgl%iden_iproc(iproc),iproc)))) then
             call pdf_failure(-4)
             return
          endif
       endif
       if (.not.ieee_is_finite(val(iip)) .or. &
            .not.all(ieee_is_finite(pgl%idenCOandMAPfactor(1:pgl%iden_iproc(iproc),iproc)))) then
          call pdf_failure(-20)
          return
       endif
       do ip=1,pgl%iden_iproc(iproc)
          if (.not.product_is_representable(val(iip),pgl%idenCOandMAPfactor(ip,iproc))) then
             call pdf_failure(-20)
             return
          endif
          pgl%val_procs(ip,iproc)=val(iip)*pgl%idenCOandMAPfactor(ip,iproc)
       enddo
       if (include_pdf) then
          do ip=1,pgl%iden_iproc(iproc)
             pdf_first=pdf_table_flavour(PDF(:,1),pgl%iden_processes(1,ip,iproc))
             pdf_second=pdf_table_flavour(PDF(:,2),pgl%iden_processes(2,ip,iproc))
             if (.not.product_is_representable(pgl%val_procs(ip,iproc),pdf_first)) then
                call pdf_failure(-20)
                return
             endif
             pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*pdf_first
             if (.not.product_is_representable(pgl%val_procs(ip,iproc),pdf_second)) then
                call pdf_failure(-20)
                return
             endif
             pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*pdf_second
          enddo
       endif
       val(iip)=sum(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc))
       val_abs(iip)=sum(abs(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc)))
       if (.not.ieee_is_finite(val(iip)) .or. .not.ieee_is_finite(val_abs(iip))) then
          call pdf_failure(-20)
          return
       endif
       if (abs(val(iip)).gt.pdf_weight_limit .or. val_abs(iip).gt.pdf_weight_limit) then
          call pdf_failure(-20)
          return
       endif
    enddo
  contains
    subroutine pdf_failure(code)
      integer,intent(in) :: code
      integer :: failed_process,failed_value
      failed_value=0
      do failed_process=ip_start,ip_end
         failed_value=failed_value+1
         val(failed_value)=0d0
         val_abs(failed_value)=0d0
         pgl%val_procs(1:pgl%iden_iproc(failed_process),failed_process)=0d0
      enddo
      if (present(status)) then
         status=code
      else
         write(*,*) 'ERROR: invalid PDF-weighted point',code
         stop 1
      endif
    end subroutine pdf_failure
  end subroutine include_PDF_and_identical_procs

  pure logical function product_is_representable(first,second)
    real(kind=8),intent(in) :: first,second
    product_is_representable=.false.
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    if ((first.ne.0d0 .and. abs(first).lt.tiny(1d0)) .or. &
         (second.ne.0d0 .and. abs(second).lt.tiny(1d0))) return
    if (first.eq.0d0 .or. second.eq.0d0) then
       product_is_representable=.true.
    else
       if (abs(first).gt.pdf_weight_limit/abs(second)) return
       if (abs(second).lt.1d0) then
          if (abs(first).lt.tiny(1d0)/abs(second)) return
       elseif (abs(first).lt.1d0) then
          if (abs(second).lt.tiny(1d0)/abs(first)) return
       endif
       product_is_representable=.true.
    endif
  end function product_is_representable

  elemental pure logical function pdf_flavour_is_supported(flavour)
    integer,intent(in) :: flavour
    pdf_flavour_is_supported=(flavour.ge.1 .and. flavour.le.6) .or. &
         (flavour.le.-1 .and. flavour.ge.-6) .or. &
         flavour.eq.21 .or. flavour.eq.22 .or. flavour.eq.99
  end function pdf_flavour_is_supported

end module pdf_wrap
