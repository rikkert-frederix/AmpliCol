program run_parameters_regression
  use run_parameters
  use particles
  implicit none
  integer,parameter :: dp=kind(1d0),n=6,nprocs=3
  type(physics_model) :: model,ignore_model,higgs_model
  integer,dimension(n,nprocs) :: processes
  integer,dimension(n,1) :: higgs_process
  character(len=256) :: default_card,custom_card,ignore_card
  integer :: i
  real(kind=dp) :: expected

  if (command_argument_count().ne.3) then
     write (*,*) 'Usage: run_parameters_regression DEFAULT CUSTOM IGNORE'
     stop 1
  endif
  call get_command_argument(1,default_card)
  call get_command_argument(2,custom_card)
  call get_command_argument(3,ignore_card)
  open(unit=99,file='/dev/null',status='unknown',action='write')

  call read_run_parameters(trim(default_card))
  call assert_close(top_mass,173d0,'default top mass')
  call assert_close(top_width,1.4915d0,'default top width')
  call assert_close(z_mass,91.188d0,'default Z mass')
  call assert_close(z_width,2.441404d0,'default Z width')
  call assert_close(w_mass,80.419002445756163d0,'default W mass')
  call assert_close(w_width,2.0476d0,'default W width')
  call assert_close(higgs_mass,125d0,'default Higgs mass')
  call assert_close(higgs_width,0.0063823389999999999d0,'default Higgs width')
  call assert_close(sw,0.47143025548407230d0,'default weak mixing sine')
  call assert_close(alphaS_MZ,0.119d0,'default alphaS')
  call assert_close(alphaEW,0.007546771114d0,'default alphaEW')
  call assert_close(sqrts,14000d0,'default collider energy')
  if (scale_choice.ne.2 .or. (.not.include_pdf) .or. (.not.use_lhapdf) .or.&
       lhapdf_member.ne.0 .or. pdf_lhaid.ne.244800) then
     write (*,*) 'Default run/PDF configuration changed unexpectedly'
     stop 1
  endif
  if (trim(lhapdfset).ne.'NNPDF23_nlo_as_0119_qed') then
     write (*,*) 'Default LHAPDF set changed unexpectedly: ',trim(lhapdfset)
     stop 1
  endif

  call read_run_parameters(trim(custom_card))
  call assert_close(top_mass,181d0,'custom top mass')
  call assert_close(top_width,2.3d0,'custom top width')
  call assert_close(z_mass,92d0,'custom Z mass')
  call assert_close(z_width,3.1d0,'custom Z width')
  call assert_close(w_mass,81d0,'custom W mass')
  call assert_close(w_width,2.7d0,'custom W width')
  call assert_close(higgs_mass,126d0,'custom Higgs mass')
  call assert_close(higgs_width,0.02d0,'custom Higgs width')
  call assert_close(sw,0.5d0,'custom weak mixing sine')
  call assert_close(alphaS_MZ,0.121d0,'custom alphaS')
  call assert_close(alphaEW,0.008d0,'custom alphaEW')
  call assert_close(sqrts,13000d0,'custom collider energy')
  if (scale_choice.ne.5 .or. include_pdf .or. use_lhapdf .or.&
       lhapdf_member.ne.3 .or. pdf_lhaid.ne.999999) then
     write (*,*) 'Custom run/PDF configuration was not read'
     stop 1
  endif
  if (trim(lhapdfset).ne.'CUSTOM_LHAPDF_SET' .or.&
       trim(internal_pdf_grid).ne.'custom_internal.grid') then
     write (*,*) 'Custom PDF names were not read'
     stop 1
  endif
  call assert_close(pTj_min,41d0,'custom jet pT cut')
  call assert_close(DRjj_min,0.51d0,'custom jet DR cut')
  call assert_close(etaj_max,4.1d0,'custom jet eta cut')
  call assert_close(sqrt_sjj_min,101d0,'custom jet invariant-mass cut')
  call assert_close(pTa_min,42d0,'custom photon pT cut')
  call assert_close(DRaa_min,0.52d0,'custom photon DR cut')
  call assert_close(etaa_max,4.2d0,'custom photon eta cut')
  call assert_close(sqrt_saa_min,102d0,'custom photon invariant-mass cut')
  call assert_close(pTl_min,43d0,'custom lepton pT cut')
  call assert_close(DRll_min,0.53d0,'custom lepton DR cut')
  call assert_close(etal_max,4.3d0,'custom lepton eta cut')
  call assert_close(sqrt_sll_min,103d0,'custom lepton invariant-mass cut')
  call assert_close(DRja_min,0.54d0,'custom jet-photon DR cut')
  call assert_close(sqrt_sja_min,104d0,'custom jet-photon invariant-mass cut')
  call assert_close(DRjl_min,0.55d0,'custom jet-lepton DR cut')
  call assert_close(sqrt_sjl_min,105d0,'custom jet-lepton invariant-mass cut')
  call assert_close(DRla_min,0.56d0,'custom lepton-photon DR cut')
  call assert_close(sqrt_sla_min,106d0,'custom lepton-photon invariant-mass cut')

  call model%init_part()
  call model%init_vert()
  call assert_close(model%get_mass(6),181d0,'model top mass')
  call assert_close(model%get_width(-6),2.3d0,'model anti-top width')
  call assert_close(model%get_mass(24),81d0,'model W mass')
  expected=81d0/0.5d0
  do i=1,model%nint
     if (model%vertex_list(i)%type.eq.17 .and.&
          all(model%vertex_list(i)%particles.eq.[24,-24,25])) then
        call assert_close(model%vertex_list(i)%coupl(1),expected,&
             'configured HWW coupling')
        exit
     endif
  enddo
  if (i.gt.model%nint) then
     write (*,*) 'Could not find HWW vertex'
     stop 1
  endif

  ! Top and W occur only in subprocess 2, and Z only in subprocess 3.
  ! A Higgs in an incoming slot must retain its nominal width.
  processes(:,1)=[25,21,1,-1,21,21]
  processes(:,2)=[1,21,21,6,-24,21]
  processes(:,3)=[1,21,21,23,21,21]
  call model%apply_final_state_widths(n,nprocs,processes)
  call assert_close(model%get_width(6),0d0,'final top width')
  call assert_close(model%get_width(-6),0d0,'final anti-top width')
  call assert_close(model%get_width(24),0d0,'final W width')
  call assert_close(model%get_width(-24),0d0,'final anti-W width')
  call assert_close(model%get_width(23),0d0,'final Z width')
  call assert_close(model%get_width(25),0.02d0,'initial-only Higgs width')

  call read_run_parameters(trim(custom_card))
  call higgs_model%init_part()
  higgs_process(:,1)=[1,21,21,21,21,25]
  call higgs_model%apply_final_state_widths(n,1,higgs_process)
  call assert_close(higgs_model%get_width(25),0d0,'final Higgs width')

  call read_run_parameters(trim(ignore_card))
  if (.not.ignore_final_state_width_fix) then
     write (*,*) 'Width-fix opt-out flag was not read'
     stop 1
  endif
  call ignore_model%init_part()
  call ignore_model%apply_final_state_widths(n,nprocs,processes)
  call assert_close(ignore_model%get_width(6),1.4915d0,'opt-out top width')
  call assert_close(ignore_model%get_width(23),2.441404d0,'opt-out Z width')
  call assert_close(ignore_model%get_width(24),2.0476d0,'opt-out W width')
  call assert_close(ignore_model%get_width(25),0.0063823389999999999d0,&
       'opt-out Higgs width')

  write (*,'(a)') 'Run-parameter and final-state-width regression passed'

contains

  subroutine assert_close(value,reference,label)
    implicit none
    real(kind=dp),intent(in) :: value,reference
    character(len=*),intent(in) :: label
    if (abs(value-reference).gt.1d-13*max(1d0,abs(reference))) then
       write (*,*) trim(label),' mismatch:',value,reference
       stop 1
    endif
  end subroutine assert_close

end program run_parameters_regression
