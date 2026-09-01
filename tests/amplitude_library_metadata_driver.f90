program amplitude_library_metadata_driver
  use common, only: phys_model
  use run_parameters, only: reset_run_parameters
  use amplitude_library, only: read_amplitude_lib
  implicit none
  character(len=1024) :: filename

  call get_command_argument(1,filename)
  if (len_trim(filename).eq.0) error stop 'missing metadata filename'
  open(unit=99,status='scratch',action='readwrite')
  call reset_run_parameters()
  call phys_model%init_part()
  call phys_model%init_vert()
  call read_amplitude_lib(trim(filename))
  error stop 'corrupt amplitude metadata was accepted'
end program amplitude_library_metadata_driver
