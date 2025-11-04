.DEFAULT_GOAL := matrix_integrate_QCD

AMPLIB = libamp.a

FILES_M_INT_QCD=bitset.o pdf.o NNPDFDriver.o mint_module.o ranmar.o HwU.o phase_space.o	\
LUPdecompose.o phase_space_gen23.o color_algebra.o math_functions.o	\
feynmanrules.o particles.o amplitude_QCD.o matrix_integrate_QCD.o common.o	\
phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
read_process_file.o multichannel.o handling_processes.o simple_integrator.o \
helper_modules.o amplitude_library.o command_line_parser.o mg_checks.o

FILES_M_RWGT_QCD=bitset.o color_algebra.o math_functions.o feynmanrules.o particles.o	\
amplitude_QCD.o matrix_reweight_QCD.o

FILES_M_UNWGT_QCD=color_algebra.o math_functions.o feynmanrules.o particles.o    \
amplitude_QCD.o matrix_unweight_QCD.o

FILES_M_COMBINE_QCD=color_algebra.o math_functions.o feynmanrules.o particles.o    \
amplitude_QCD.o matrix_combine_QCD.o

FC=gfortran
FFLAGS=-ffast-math -O3 -ffree-line-length-0
#FFLAGS=-fbounds-check -g -ffpe-trap=invalid,zero,overflow,underflow,denormal

# Files for all executables

%.o: %.f03
	$(FC) $(FFLAGS) -c -I. $<
%.o: %.f95
	$(FC) $(FFLAGS) -c -I. $<
%.o: %.f
	$(FC) $(FFLAGS) -c -I. $<
%.o: PDF/%.f
	$(FC) $(FFLAGS) -c -I. -IPDF $<
%.o: SimpleMint/%.f
	$(FC) $(FFLAGS) -c -I. -ISimpleMint $<
%.o: SimpleMint/%.f90
	$(FC) $(FFLAGS) -c -I. -ISimpleMint $<
%.o: SimpleMint/%.f03
	$(FC) $(FFLAGS) -c -I. -ISimpleMint $<
%.o: PhaseSpace/%.f90
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace $<
%.o: PhaseSpace/%.f03
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace $<

#.PHONY: $(AMPLIB)
$(AMPLIB):
	@rm -f $(AMPLIB)
	@for src in $(shell find library/ -name 'amp*_lib.f03'); do $(FC) $(FFLAGS) -c $$src ; done
	$(FC) $(FFLAGS) -c library/amplib.f03
	ar rcs $(AMPLIB) amp*_lib.o amplib.o

matrix_integrate_QCD: cleanlib $(FILES_M_INT_QCD) $(AMPLIB)
	$(FC) $(FFLAGS) -o matrix_integrate_QCD $(FILES_M_INT_QCD) $(AMPLIB) `lhapdf-config --ldflags` -lstdc++
	@rm -f $(AMPLIB)

matrix_integrate_QCD_library: $(FILES_M_INT_QCD) $(AMPLIB)
	$(FC) $(FFLAGS) -o matrix_integrate_QCD $(FILES_M_INT_QCD) $(AMPLIB) `lhapdf-config --ldflags` -lstdc++
	@rm -f $(AMPLIB)

matrix_reweight_QCD: $(FILES_M_RWGT_QCD) 
	$(FC) $(FFLAGS) -o matrix_reweight_QCD $(FILES_M_RWGT_QCD)

matrix_unweight_QCD: $(FILES_M_UNWGT_QCD)
	$(FC) $(FFLAGS) -o matrix_unweight_QCD $(FILES_M_UNWGT_QCD) `lhapdf-config --ldflags` -lstdc++

matrix_combine_QCD: $(FILES_M_COMBINE_QCD)
	$(FC) $(FFLAGS) -o matrix_combine_QCD $(FILES_M_COMBINE_QCD)

clean:
	rm -f *.o *.mod library/amp*.f03 library/amp*.data library/amplitudes.bin

cleanlib:
	rm -f amp*lib.o amp*lib.mod library/amp*lib*
	$(FC) $(FFLAGS) -c library/dummy.f03
	ar rcs $(AMPLIB) dummy.o

matrix_reweight_QCD.o : amplitude_QCD.o math_functions.o particles.o 
ranmar.o : mint_module.o
phase_space_gen23.o : phase_space.o LUPdecompose.o
phase_space_genpt.o : phase_space.o
haag.o : phase_space.o
amplitude_QCD.o : bitset.o math_functions.o feynmanrules.o color_algebra.o particles.o
matrix_integrate_QCD.o : amplitude_QCD.o phase_space_gen23.o mint_module.o common.o math_functions.o particles.o phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o read_process_file.o multichannel.o handling_processes.o simple_integrator.o amplitude_library.o command_line_parser.o mg_checks.o
common.o : particles.o simple_integrator.o
handling_events.o : common.o mint_module.o handling_processes.o simple_integrator.o
read_process_file.o : mint_module.o phase_space_gen23.o cuts.o handling_processes.o simple_integrator.o
multichannel.o : handling_processes.o mint_module.o math_functions.o simple_integrator.o
handling_processes.o : math_functions.o common.o phase_space.o amplitude_QCD.o
cuts.o : common.o particles.o handling_processes.o
pdf_wrap.o : handling_processes.o
simple_integrator.o : helper_modules.o
amplitude_library.o : handling_processes.o read_process_file.o
mg_checks.o : handling_processes.o
