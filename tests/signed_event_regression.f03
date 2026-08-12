program signed_event_regression
  use handling_events, only: event_update_wgt
  implicit none

  integer :: input_unit,output_unit,nup,idprup
  real(kind=8) :: xwgtup,scalup,aqedup,aqcdup
  real(kind=8),dimension(3) :: sampling_weights
  character(len=1024) :: line

  open(newunit=input_unit,status='scratch',action='readwrite')
  open(newunit=output_unit,status='scratch',action='readwrite')
  write(input_unit,'(a)') '<event>'
  write(input_unit,503) 4,1,-1d0,100d0,0.0075d0,0.118d0
  write(input_unit,'(a)') '# signed-interference payload'
  write(input_unit,'(a)') '</event>'
  rewind(input_unit)

  sampling_weights=[2d0,3d0,4d0]
  call event_update_wgt(input_unit,output_unit,sampling_weights)
  rewind(output_unit)
  read(output_unit,'(a)') line
  if (trim(line).ne.'<event>') then
     write (*,*) 'Signed event record lost its opening tag:',trim(line)
     stop 1
  endif
  read(output_unit,503) nup,idprup,xwgtup,scalup,aqedup,aqcdup
  if (nup.ne.4 .or. idprup.ne.1 .or. abs(xwgtup+2d0).gt.1d-12) then
     write (*,*) 'Final LHE weight did not retain the negative event sign:',&
          nup,idprup,xwgtup
     stop 1
  endif
  write (*,'(a,es12.4)') 'Signed-event regression passed; final weight=',xwgtup

503 format(1x,i2,1x,i6,4(1x,e14.8))
end program signed_event_regression
