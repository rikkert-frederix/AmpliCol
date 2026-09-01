module integration_histograms
  ! Accuracy-only integration histograms.
  !
  ! Every (channel,integral) pair is an independent stratum.  Contributions
  ! associated with one random point (notably R and all mapped CS dipoles) are
  ! first summed in point(:,:) and only then squared.  Leaf sample means are
  ! combined across iterations by point count and summed across leaves.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer,parameter :: curve_nlo=1,curve_born=2
  integer,parameter :: title_length=160
  real(kind=8),parameter :: histogram_value_limit=0.25d0*huge(1d0)**0.25d0

  type :: histogram_data
     integer :: label=0
     integer :: nbin=0
     real(kind=8) :: xmin=0d0,xmax=0d0,step=0d0
     character(len=title_length) :: title=''
     real(kind=8),allocatable :: point(:,:)
     logical,allocatable :: touched(:)
     integer,allocatable :: touched_bins(:)
     integer :: ntouched=0
     real(kind=8),allocatable :: accum(:,:,:),accum2(:,:,:)
     real(kind=8),allocatable :: result(:,:,:),uncertainty(:,:,:)
  end type histogram_data

  type(histogram_data),allocatable,save :: histograms(:)
  integer,allocatable,save :: leaf_offset(:)
  integer,allocatable,save :: label_index(:)
  integer,allocatable,save :: touched_histograms(:)
  integer(kind=8),allocatable,save :: npoints_iter(:),npoints_total(:)
  integer,save :: nleaves=0,current_leaf=0,nhistograms_touched=0
  logical,save :: initialized=.false.,point_open=.false.,point_invalid=.false.,write_nlo=.false.
  character(len=256),save :: output_file=''

  public :: histogram_initialize,histogram_book,histogram_begin_point
  public :: histogram_fill,histogram_commit_point,histogram_finalize_iteration
  public :: histogram_write,histogram_get_bin

contains

  subroutine histogram_initialize(nintegrals,has_nlo,filename)
    integer,intent(in) :: nintegrals(:)
    logical,intent(in) :: has_nlo
    character(len=*),intent(in) :: filename
    integer :: i

    if (allocated(histograms)) deallocate(histograms)
    if (allocated(leaf_offset)) deallocate(leaf_offset)
    if (allocated(label_index)) deallocate(label_index)
    if (allocated(touched_histograms)) deallocate(touched_histograms)
    if (allocated(npoints_iter)) deallocate(npoints_iter)
    if (allocated(npoints_total)) deallocate(npoints_total)

    allocate(histograms(0),leaf_offset(size(nintegrals)+1),label_index(0),touched_histograms(0))
    leaf_offset(1)=0
    do i=1,size(nintegrals)
       leaf_offset(i+1)=leaf_offset(i)+nintegrals(i)
    enddo
    nleaves=leaf_offset(size(leaf_offset))
    allocate(npoints_iter(nleaves),npoints_total(nleaves))
    npoints_iter=0_8
    npoints_total=0_8
    write_nlo=has_nlo
    output_file=trim(filename)
    current_leaf=0
    nhistograms_touched=0
    point_open=.false.
    point_invalid=.false.
    initialized=.true.
  end subroutine histogram_initialize

  subroutine histogram_book(label,title,nbin,xmin,xmax)
    integer,intent(in) :: label,nbin
    character(len=*),intent(in) :: title
    real(kind=8),intent(in) :: xmin,xmax
    type(histogram_data),allocatable :: tmp(:)
    integer,allocatable :: new_label_index(:),new_touched_histograms(:)
    integer :: i,nold

    call require_initialized('histogram_book')
    if (label.le.0 .or. nbin.le.0 .or. .not.ieee_is_finite(xmin) .or. &
         .not.ieee_is_finite(xmax)) then
       write(*,*) 'ERROR: invalid histogram booking:',label,nbin,xmin,xmax
       stop 1
    endif
    if (abs(xmin).gt.histogram_value_limit .or. &
         abs(xmax).gt.histogram_value_limit .or. xmax.le.xmin) then
       write(*,*) 'ERROR: invalid histogram booking:',label,nbin,xmin,xmax
       stop 1
    endif
    do i=1,size(histograms)
       if (histograms(i)%label.eq.label) then
          write(*,*) 'ERROR: duplicate histogram label:',label
          stop 1
       endif
    enddo

    nold=size(histograms)
    allocate(tmp(nold+1))
    if (nold.gt.0) tmp(1:nold)=histograms
    tmp(nold+1)%label=label
    tmp(nold+1)%title=title
    tmp(nold+1)%nbin=nbin
    tmp(nold+1)%xmin=xmin
    tmp(nold+1)%xmax=xmax
    tmp(nold+1)%step=(xmax-xmin)/dble(nbin)
    allocate(tmp(nold+1)%point(2,nbin))
    allocate(tmp(nold+1)%touched(nbin),tmp(nold+1)%touched_bins(nbin))
    allocate(tmp(nold+1)%accum(2,nleaves,nbin),tmp(nold+1)%accum2(2,nleaves,nbin))
    allocate(tmp(nold+1)%result(2,nleaves,nbin),tmp(nold+1)%uncertainty(2,nleaves,nbin))
    tmp(nold+1)%point=0d0
    tmp(nold+1)%touched=.false.
    tmp(nold+1)%touched_bins=0
    tmp(nold+1)%ntouched=0
    tmp(nold+1)%accum=0d0
    tmp(nold+1)%accum2=0d0
    tmp(nold+1)%result=0d0
    tmp(nold+1)%uncertainty=0d0
    call move_alloc(tmp,histograms)
    allocate(new_touched_histograms(nold+1))
    new_touched_histograms=0
    if (nold.gt.0) new_touched_histograms(1:nold)=touched_histograms
    call move_alloc(new_touched_histograms,touched_histograms)
    if (label.gt.size(label_index)) then
       allocate(new_label_index(label))
       new_label_index=0
       if (size(label_index).gt.0) new_label_index(1:size(label_index))=label_index
       call move_alloc(new_label_index,label_index)
    endif
    label_index(label)=nold+1
  end subroutine histogram_book

  subroutine histogram_begin_point(ichan,iint)
    integer,intent(in) :: ichan,iint
    call require_initialized('histogram_begin_point')
    if (point_open) then
       write(*,*) 'ERROR: histogram point already open'
       stop 1
    endif
    if (ichan.lt.1 .or. ichan.ge.size(leaf_offset)) then
       write(*,*) 'ERROR: invalid histogram channel:',ichan
       stop 1
    endif
    current_leaf=leaf_offset(ichan)+iint
    if (current_leaf.le.leaf_offset(ichan) .or. current_leaf.gt.leaf_offset(ichan+1)) then
       write(*,*) 'ERROR: invalid histogram integral:',ichan,iint
       stop 1
    endif
    point_open=.true.
    point_invalid=.false.
  end subroutine histogram_begin_point

  subroutine histogram_fill(label,x,wgt_nlo,wgt_born)
    integer,intent(in) :: label
    real(kind=8),intent(in) :: x,wgt_nlo,wgt_born
    integer :: ih,ibin
    if (.not.point_open) then
       write(*,*) 'ERROR: histogram_fill called outside a point'
       stop 1
    endif
    if (point_invalid) return
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(wgt_nlo) .or. &
         .not.ieee_is_finite(wgt_born)) then
       point_invalid=.true.
       call report_invalid_histogram_point(label)
       return
    endif
    if (abs(wgt_nlo).gt.histogram_value_limit .or. &
         abs(wgt_born).gt.histogram_value_limit) then
       point_invalid=.true.
       call report_invalid_histogram_point(label)
       return
    endif
    if (wgt_nlo.eq.0d0 .and. wgt_born.eq.0d0) return
    ih=find_histogram(label)
    if (x.lt.histograms(ih)%xmin .or. x.ge.histograms(ih)%xmax) return
    ibin=int((x-histograms(ih)%xmin)/histograms(ih)%step)+1
    if (ibin.lt.1 .or. ibin.gt.histograms(ih)%nbin) return
    if (.not.histograms(ih)%touched(ibin)) then
       if (histograms(ih)%ntouched.eq.0) then
          nhistograms_touched=nhistograms_touched+1
          touched_histograms(nhistograms_touched)=ih
       endif
       histograms(ih)%ntouched=histograms(ih)%ntouched+1
       histograms(ih)%touched_bins(histograms(ih)%ntouched)=ibin
       histograms(ih)%touched(ibin)=.true.
    endif
    if (abs(histograms(ih)%point(curve_nlo,ibin)+wgt_nlo).gt.histogram_value_limit .or. &
         abs(histograms(ih)%point(curve_born,ibin)+wgt_born).gt.histogram_value_limit) then
       point_invalid=.true.
       call report_invalid_histogram_point(label)
       return
    endif
    histograms(ih)%point(curve_nlo,ibin)=histograms(ih)%point(curve_nlo,ibin)+wgt_nlo
    histograms(ih)%point(curve_born,ibin)=histograms(ih)%point(curve_born,ibin)+wgt_born
  end subroutine histogram_fill

  subroutine histogram_commit_point()
    integer :: iactive,ih,itouched,ibin
    if (.not.point_open) then
       write(*,*) 'ERROR: histogram_commit_point called without begin'
       stop 1
    endif
    npoints_iter(current_leaf)=npoints_iter(current_leaf)+1_8
    do iactive=1,nhistograms_touched
       ih=touched_histograms(iactive)
       do itouched=1,histograms(ih)%ntouched
          ibin=histograms(ih)%touched_bins(itouched)
          if (.not.point_invalid) then
             histograms(ih)%accum(:,current_leaf,ibin)=&
                  histograms(ih)%accum(:,current_leaf,ibin)+histograms(ih)%point(:,ibin)
             histograms(ih)%accum2(:,current_leaf,ibin)=&
                  histograms(ih)%accum2(:,current_leaf,ibin)+histograms(ih)%point(:,ibin)**2
          endif
          histograms(ih)%point(:,ibin)=0d0
          histograms(ih)%touched(ibin)=.false.
       enddo
       histograms(ih)%ntouched=0
    enddo
    nhistograms_touched=0
    current_leaf=0
    point_open=.false.
    point_invalid=.false.
  end subroutine histogram_commit_point

  subroutine report_invalid_histogram_point(label)
    integer,intent(in) :: label
    integer(kind=8),save :: invalid_points=0_8
    invalid_points=invalid_points+1_8
    if (invalid_points.le.10_8 .or. mod(invalid_points,1000_8).eq.0_8) then
       write(*,'(a,i0,a,i0)') 'WARNING: discarded invalid histogram point, count=',&
            invalid_points,' label=',label
       write(99,'(a,i0,a,i0)') 'WARNING: discarded invalid histogram point, count=',&
            invalid_points,' label=',label
    elseif (invalid_points.eq.11_8) then
       write(*,'(a)') 'WARNING: further invalid-histogram-point messages suppressed'
       write(99,'(a)') 'WARNING: further invalid-histogram-point messages suppressed'
    endif
  end subroutine report_invalid_histogram_point

  subroutine histogram_finalize_iteration()
    integer :: ih,ileaf,icurve,ibin
    integer(kind=8) :: nold,nnew
    real(kind=8) :: mean_iter,mean2_iter,unc_iter
    if (point_open) then
       write(*,*) 'ERROR: cannot finalize histograms with an open point'
       stop 1
    endif
    do ileaf=1,nleaves
       nnew=npoints_iter(ileaf)
       if (nnew.eq.0_8) cycle
       nold=npoints_total(ileaf)
       do ih=1,size(histograms)
          do ibin=1,histograms(ih)%nbin
             do icurve=1,2
                mean_iter=histograms(ih)%accum(icurve,ileaf,ibin)/dble(nnew)
                mean2_iter=histograms(ih)%accum2(icurve,ileaf,ibin)/dble(nnew)
                unc_iter=sqrt(abs(mean2_iter-mean_iter**2)/dble(nnew))
                if (nold.eq.0_8) then
                   histograms(ih)%result(icurve,ileaf,ibin)=mean_iter
                   histograms(ih)%uncertainty(icurve,ileaf,ibin)=unc_iter
                else
                   call combine_mean_and_uncertainty(&
                        histograms(ih)%result(icurve,ileaf,ibin),&
                        histograms(ih)%uncertainty(icurve,ileaf,ibin),nold,&
                        mean_iter,unc_iter,nnew)
                endif
             enddo
          enddo
          histograms(ih)%accum(:,ileaf,:)=0d0
          histograms(ih)%accum2(:,ileaf,:)=0d0
       enddo
       npoints_total(ileaf)=nold+nnew
       npoints_iter(ileaf)=0_8
    enddo
    call histogram_write()
  end subroutine histogram_finalize_iteration

  subroutine combine_mean_and_uncertainty(mean,unc,nold,mean_new,unc_new,nnew)
    real(kind=8),intent(inout) :: mean,unc
    integer(kind=8),intent(in) :: nold,nnew
    real(kind=8),intent(in) :: mean_new,unc_new
    integer(kind=8) :: ntotal
    ntotal=nold+nnew
    unc=sqrt((unc**2*dble(nold)**2+unc_new**2*dble(nnew)**2)/dble(ntotal)**2+&
         dble(nold)*dble(nnew)*(mean-mean_new)**2/dble(ntotal)**3)
    mean=(dble(nold)*mean+dble(nnew)*mean_new)/dble(ntotal)
  end subroutine combine_mean_and_uncertainty

  subroutine histogram_write()
    integer :: unit,ih
    call require_initialized('histogram_write')
    open(newunit=unit,file=trim(output_file),status='replace',action='write')
    write(unit,'(a)') '##& xmin & xmax & central value & dy'
    write(unit,'(a)') ''
    do ih=1,size(histograms)
       if (write_nlo) call write_curve(unit,histograms(ih),curve_nlo,'NLO')
       call write_curve(unit,histograms(ih),curve_born,'LO')
    enddo
    close(unit)
  end subroutine histogram_write

  subroutine histogram_get_bin(label,ibin,nlo,value,error)
    integer,intent(in) :: label,ibin
    logical,intent(in) :: nlo
    real(kind=8),intent(out) :: value,error
    integer :: ih,curve
    ih=find_histogram(label)
    if (ibin.lt.1 .or. ibin.gt.histograms(ih)%nbin) then
       write(*,*) 'ERROR: invalid histogram bin:',label,ibin
       stop 1
    endif
    if (nlo) then
       curve=curve_nlo
    else
       curve=curve_born
    endif
    value=sum(histograms(ih)%result(curve,:,ibin))
    error=sqrt(sum(histograms(ih)%uncertainty(curve,:,ibin)**2))
  end subroutine histogram_get_bin

  subroutine write_curve(unit,hist,curve,curve_name)
    integer,intent(in) :: unit,curve
    type(histogram_data),intent(in) :: hist
    character(len=*),intent(in) :: curve_name
    integer :: ibin
    real(kind=8) :: xlow,xhigh,value,error
    write(unit,'(a,i4.4,a,a,a,a,a)') '<histogram> ',hist%nbin,' " ',&
         trim(hist%title),' |T@',trim(curve_name),'"'
    do ibin=1,hist%nbin
       xlow=hist%xmin+dble(ibin-1)*hist%step
       xhigh=xlow+hist%step
       value=sum(hist%result(curve,:,ibin))
       error=sqrt(sum(hist%uncertainty(curve,:,ibin)**2))
       write(unit,'(4(2x,e14.7))') xlow,xhigh,value,error
    enddo
    write(unit,'(a)') '<\histogram>'
    write(unit,'(a)') ''
  end subroutine write_curve

  integer function find_histogram(label)
    integer,intent(in) :: label
    if (label.ge.1 .and. label.le.size(label_index)) then
       find_histogram=label_index(label)
       if (find_histogram.ne.0) return
    endif
    write(*,*) 'ERROR: histogram label was not booked:',label
    stop 1
  end function find_histogram

  subroutine require_initialized(caller)
    character(len=*),intent(in) :: caller
    if (.not.initialized) then
       write(*,*) 'ERROR: ',trim(caller),' called before histogram_initialize'
       stop 1
    endif
  end subroutine require_initialized

end module integration_histograms
