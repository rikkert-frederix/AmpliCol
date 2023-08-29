! gfortran -fbounds-check -o test_QCD test_QCD.f03 amplitude_real_QCD.f03

program test_QCD
  use amplitude_mod
  implicit none
  type(amplitude) :: amplitudes
  integer :: n
  integer,dimension(:),allocatable :: part,order
  integer :: i

  n=5
  allocate(part(n))
  allocate(order(n))

  part(1:n)=[21,21,1,-1,21]
  order(1:n)=[3,1,2,5,4]
  
  call amplitudes%init_OneOrder(n,part,order)

  do i=1,amplitudes%n_vert
     write (*,*) i,':',amplitudes%interaction_list(i)%currents(1:2),':',amplitudes%interaction_list(i)%type
  enddo
  
  do i=1,amplitudes%n_cur
     write (*,*) i,':',amplitudes%current_list(i)%order(1:popcnt(amplitudes%current_list(i)%bin)),':',&
          amplitudes%current_list(i)%type,':',amplitudes%current_list(i)%bin
     write (*,*) amplitudes%current_list(i)%vertices(1:amplitudes%current_list(i)%n_vert)
  enddo

end program test_QCD
