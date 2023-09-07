! gfortran -fbounds-check -o test_QCD feynmanrules.f03 amplitude_QCD.f03 test_QCD.f03 

program test_QCD
  use amplitude_mod
  implicit none
  type(amplitude) :: amplitudes
  integer :: n
  integer,dimension(:),allocatable :: part,order
  integer :: i
  real(kind=8),dimension(:,:),allocatable :: p
  
  n=5
  allocate(part(n))
  allocate(order(n))

  part(1:n)=[-1,21,21,21,-1]
  order(1:n)=[1,2,3,4,5]


  call amplitudes%init_OneOrder(n,part,order)

  do i=1,amplitudes%n_vert
     write (*,*) i,':',amplitudes%interaction_list(i)%currents(1:2),':',amplitudes%interaction_list(i)%type
  enddo
  
  do i=1,amplitudes%n_cur
     write (*,*) i,':',amplitudes%current_list(i)%order(1:popcnt(amplitudes%current_list(i)%bin)),':',&
          amplitudes%current_list(i)%type,':',amplitudes%current_list(i)%bin
     write (*,*) amplitudes%current_list(i)%vertices(1:amplitudes%current_list(i)%n_vert)
  enddo

  
  allocate(p(0:3,n))
  
  if (n.eq.5) then
     p(0:3,1)=(/   500.00000000000000      ,   0.0000000000000000      ,   0.0000000000000000      ,   500.00000000000000      /)
     p(0:3,2)=(/   500.00000000000000      ,   0.0000000000000000      ,   0.0000000000000000      ,  -500.00000000000000      /)
     p(0:3,3)=(/   362.35913008614813      ,   211.67753703547245      ,  -165.00336466515242      ,   243.45564097092699      /)
     p(0:3,4)=(/   248.22961274839184      ,   68.500483807498668      ,   237.98757697004567      ,   16.956932838275925      /)
     p(0:3,5)=(/   389.41125716545980      ,  -280.17802084297108      ,  -72.984212304893376      ,  -260.41257380920291      /)
  endif
  
  call amplitudes%evaluate_OneOrder(n,p)

  
end program test_QCD
