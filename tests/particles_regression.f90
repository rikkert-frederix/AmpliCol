program particles_regression
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use particles, only: physics_model
  implicit none
  type(physics_model) :: model
  character(len=32) :: mode
  integer :: minimum_integer,value

  mode='success'
  if (command_argument_count().ge.1) call get_command_argument(1,mode)
  open(unit=99,status='scratch',action='write')
  minimum_integer=-huge(0)-1

  select case (trim(mode))
  case ('success')
     if (model%is_lepton_any(minimum_integer)) &
          error stop 'minimum integer was classified as a lepton'
     if (model%is_w_aux_tensor(minimum_integer)) &
          error stop 'minimum integer was classified as an auxiliary tensor'
     if (model%is_massive_vector(minimum_integer)) &
          error stop 'minimum integer was classified as a massive vector'
     if (model%is_massless_fermion(minimum_integer)) &
          error stop 'minimum integer was classified as a massless fermion'
     if (model%is_jet(minimum_integer)) &
          error stop 'minimum integer was classified as a jet'
     if (.not.model%is_singlet(minimum_integer)) &
          error stop 'singlet classification changed for a non-coloured sentinel'
     if (model%is_massless_fermion(22)) &
          error stop 'a photon was classified as a massless fermion'
     if (model%is_jet(25)) error stop 'a Higgs was classified as a jet'

     call model%init_part()
     call model%init_vert()
     if (.not.model%is_massless_fermion(1)) &
          error stop 'a light quark was not classified as massless'
     if (.not.model%is_jet(21)) error stop 'a gluon was not classified as a jet'
     if (model%is_chiral_eligible(6)) &
          error stop 'a massive top quark was classified as chiral eligible'
     call model%set_width(23,0d0)
     if (model%get_width(23).ne.0d0) error stop 'valid width update failed'
     write(*,'(a)') 'Particles regression: PASS'
  case ('nan-width')
     call model%init_part()
     call model%set_width(23,ieee_value(0d0,ieee_quiet_nan))
  case ('huge-width')
     call model%init_part()
     call model%set_width(23,huge(1d0))
  case ('tiny-width')
     call model%init_part()
     call model%set_width(23,tiny(1d0))
  case ('vertices-before-particles')
     call model%init_vert()
  case ('stale-vertices')
     call model%init_part()
     call model%init_vert()
     call model%init_part()
     value=model%get_inter_dim(0)
  case default
     error stop 'unknown particles regression mode'
  end select
  if (value.eq.-huge(0)) write(*,*) 'unreachable'
  close(99)
end program particles_regression
