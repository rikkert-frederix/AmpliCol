.DEFAULT_GOAL := matrix_integrate_QCD

FILES_M_INT_QCD=pdf.o NNPDFDriver.o mint_module.o ranmar.o HwU.o phase_space.o	\
LUPdecompose.o phase_space_gen23.o color_algebra.o math_functions.o	\
feynmanrules.o particles.o amplitude_QCD.o matrix_integrate_QCD.o common.o	\
phase_space_genpt.o phase_space_haag.o

FILES_M_RWGT_QCD=color_algebra.o math_functions.o feynmanrules.o particles.o	\
amplitude_QCD.o matrix_reweight_QCD.o

FILES_M_UNWGT_QCD=color_algebra.o math_functions.o feynmanrules.o particles.o    \
amplitude_QCD.o matrix_unweight_QCD.o

FILES_M_COMBINE_QCD=color_algebra.o math_functions.o feynmanrules.o particles.o    \
amplitude_QCD.o matrix_combine_QCD.o

FC=gfortran
FFLAGS=-ffast-math -O3
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
%.o: PhaseSpace/%.f90
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace $<
%.o: PhaseSpace/%.f03
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace $<


matrix_integrate_QCD:  $(FILES_M_INT_QCD) 
	$(FC) $(FFLAGS) -o matrix_integrate_QCD $(FILES_M_INT_QCD) `lhapdf-config --ldflags` -lstdc++

matrix_reweight_QCD: $(FILES_M_RWGT_QCD) 
	$(FC) $(FFLAGS) -o matrix_reweight_QCD $(FILES_M_RWGT_QCD)

matrix_unweight_QCD: $(FILES_M_UNWGT_QCD)
	$(FC) $(FFLAGS) -o matrix_unweight_QCD $(FILES_M_UNWGT_QCD) `lhapdf-config --ldflags` -lstdc++

matrix_combine_QCD: $(FILES_M_COMBINE_QCD)
	$(FC) $(FFLAGS) -o matrix_combine_QCD $(FILES_M_COMBINE_QCD)

clean:
	rm -f *.o *.mod

matrix_reweight_QCD.o : amplitude_QCD.o math_functions.o particles.o 
ranmar.o : mint_module.o
phase_space_gen23.o : phase_space.o LUPdecompose.o
phase_space_genpt.o : phase_space.o
haag.o : phase_space.o
amplitude_QCD.o : math_functions.o feynmanrules.o color_algebra.o particles.o
matrix_integrate_QCD.o : amplitude_QCD.o phase_space_gen23.o mint_module.o common.o math_functions.o particles.o phase_space_genpt.o phase_space_haag.o
common.o : amplitude_QCD.o
