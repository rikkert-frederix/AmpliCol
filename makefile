.DEFAULT_GOAL := matrix_integrate_QCD

FILES_M_INT=mint_module.o MC_integer.o ranmar.o HwU.o LUPdecompose.o		\
phase_space_gen23.o haag.o color_algebra.o math_functions.o feynmanrules.o	\
amplitude_QCD.o amplitude_real.o matrix_integrate.o common.o

FILES_M_INT_QCD=pdf.o NNPDFDriver.o mint_module.o MC_integer.o ranmar.o HwU.o	\
LUPdecompose.o phase_space_gen23.o haag.o color_algebra.o math_functions.o	\
feynmanrules.o amplitude_QCD.o amplitude_real.o matrix_integrate_QCD.o		\
common.o

FILES_INT_QCD=mint_module.o MC_integer.o ranmar.o HwU.o LUPdecompose.o		\
phase_space_gen23.o haag.o math_functions.o color_algebra.o feynmanrules.o	\
amplitude_QCD.o amplitude_real.o integrate_QCD.o common.o

FILES_M_RWGT=random.o color_algebra.o amplitude_real.o math_functions.o	\
feynmanrules.o amplitude_QCD.o matrix_reweight.o

FILES_M_RWGT_QCD=random.o color_algebra.o math_functions.o	\
feynmanrules.o amplitude_QCD.o matrix_reweight_QCD.o

FC=gfortran
FFLAGS=-ffast-math -O3
#FFLAGS=-fbounds-check

# Files for all executables

%.o: %.f03
	$(FC) $(FFLAGS) -c -I. $<
%.o: %.f95
	$(FC) $(FFLAGS) -c -I. $<
%.o: %.f
	$(FC) $(FFLAGS) -c -I. $<
%.o: PDF/%.f
	$(FC) $(FFLAGS) -c -I. -IPDF $<
%.o: simple_mint/%.f
	$(FC) $(FFLAGS) -c -I. -Isimple_mint $<
%.o: simple_mint/%.f90
	$(FC) $(FFLAGS) -c -I. -Isimple_mint $<
%.o: PhaseSpace_BycklingKajantie/%.f90
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace_BycklingKajantie $<
%.o: PhaseSpace_haag/%.f90
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace_haag $<


matrix_integrate: $(FILES_M_INT)
	$(FC) $(FFLAGS) -o matrix_integrate $(FILES_M_INT)

matrix_integrate_QCD:  $(FILES_M_INT_QCD)
	$(FC) $(FFLAGS) -o matrix_integrate_QCD $(FILES_M_INT_QCD)

integrate_QCD: $(FILES_INT_QCD)
	$(FC) $(FFLAGS) -o integrate_QCD $(FILES_INT_QCD)

matrix_reweight: $(FILES_M_RWGT)
	$(FC) $(FFLAGS) -o matrix_reweight $(FILES_M_RWGT)

matrix_reweight_QCD: $(FILES_M_RWGT_QCD)
	$(FC) $(FFLAGS) -o matrix_reweight_QCD $(FILES_M_RWGT_QCD)

clean:
	rm *.o *.mod

matrix_reweight.o : amplitude_real.o amplitude_QCD.o
matrix_reweight_QCD.o : amplitude_QCD.o math_functions.o
ranmar.o : mint_module.o
phase_space_gen23.o : common.o
haag.o : common.o
amplitude_QCD.o : math_functions.o feynmanrules.o color_algebra.o
amplitude_real.o : color_algebra.o
matrix_integrate.o : amplitude_real.o amplitude_QCD.o phase_space_gen23.o haag.o mint_module.o common.o
matrix_integrate_QCD.o : amplitude_real.o amplitude_QCD.o phase_space_gen23.o haag.o mint_module.o common.o
integrate_QCD.o : amplitude_real.o amplitude_QCD.o phase_space_gen23.o haag.o mint_module.o common.o
common.o : amplitude_real.o amplitude_QCD.o
