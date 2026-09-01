program color_nlc_regression
  use color_algebra, only: check_NLC,check_NLC_1qqbar,check_NLC_2qqbar,&
       check_NLC_2qqbar_SF,ipnext
  implicit none
  integer :: ng,next,ri,rj,ii,jj,acc
  integer,allocatable :: first(:),second(:),identity(:)

  do ng=1,5
     allocate(first(ng),second(ng),identity(ng))
     identity=[(ii,ii=1,ng)]
     first=identity
     do
        second=identity
        do
           acc=123456789
           call check_NLC(ng,first,second,acc)
           call require_coefficient(acc,'all-gluon')
           acc=123456789
           call check_NLC_1qqbar(ng+2,first,second,acc)
           call require_coefficient(acc,'one-quark-line')
           if (ng.le.4) then
              next=ng+4
              do ri=0,ng
                 do rj=0,ng
                    do ii=1,2
                       do jj=1,2
                          acc=123456789
                          call check_NLC_2qqbar(next,first,second,ri,rj,ii,jj,acc)
                          call require_coefficient(acc,'two-quark-line different-flavour')
                          acc=123456789
                          call check_NLC_2qqbar_SF(next,first,second,ri,rj,ii,jj,acc)
                          call require_coefficient(acc,'two-quark-line same-flavour')
                       enddo
                    enddo
                 enddo
              enddo
           endif
           call ipnext(second,ng)
           if (all(second.eq.identity)) exit
        enddo
        call ipnext(first,ng)
        if (all(first.eq.identity)) exit
     enddo
     deallocate(first,second,identity)
  enddo
  write(*,'(a)') 'Colour-NLC regression: PASS'

contains
  subroutine require_coefficient(value,label)
    integer,intent(in) :: value
    character(len=*),intent(in) :: label
    if (value.ne.-1 .and. value.ne.0 .and. value.ne.1 .and. value.ne.99) then
       write(*,*) 'Invalid ',trim(label),' NLC coefficient:',value
       error stop 1
    endif
  end subroutine require_coefficient
end program color_nlc_regression
