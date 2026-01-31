.DEFAULT_GOAL := matrix_integrate_QCD

FC = gfortran
#FFLAGS= -fbounds-check -g -ffpe-trap=invalid,zero,overflow,underflow,denormal
#FFLAGS = -ffast-math -O3 -mcmodel=large
FFLAGS = -ffast-math -O3

# ----------------------------------------------------------------------
# 1. Detect amplitude sources and group them
# ----------------------------------------------------------------------

# All amplitude source files: amp<GROUP>_<ID>_lib.f03
AMPSRC := $(shell find library/ -name 'amp*_lib.f03')

# All amplitude object files
AMPOBJ := $(notdir $(AMPSRC:.f03=.o))

# Explicit rule so Make knows how to build amplitude objects
$(AMPOBJ): %.o : library/%.f03
	$(FC) $(FFLAGS) -fPIC -c -I. -Ilibrary $<

# Extract GROUP names (amp32, amp81, amp252, ...)
AMPGROUPS := $(sort $(foreach f,$(AMPOBJ),$(word 1,$(subst _, ,$(notdir $(basename $f))))))

# One shared library per group
AMPLIBS := $(foreach g,$(AMPGROUPS),lib$(g).so)

# ----------------------------------------------------------------------
# 2. Generic compilation rules
# ----------------------------------------------------------------------

%.o: %.f03
	$(FC) $(FFLAGS) -c -I. $<

%.o: %.f95
	$(FC) $(FFLAGS) -c -I. $<

%.o: %.f
	$(FC) $(FFLAGS) -c -I. $<

%.o: PDF/%.f
	$(FC) $(FFLAGS) -c -I. -IPDF $<

%.o: PDF/%.cc
	$(CXX) $(CXXFLAGS) -c -I. -IPDF $< -std=c++11 -stdlib=libc++

%.o: PDF/%.f90
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

%.o: library/%.f03
	$(FC) $(FFLAGS) -fPIC -c -I. -Ilibrary $<

# ----------------------------------------------------------------------
# 3. Build one shared library per amplitude group
# ----------------------------------------------------------------------

define one_lib_template
lib$(1).so: $(filter $(1)_%_lib.o,$(AMPOBJ))
	$$(FC) -shared -o $$@ $$^
endef

$(foreach g,$(AMPGROUPS),$(eval $(call one_lib_template,$(g))))

# ----------------------------------------------------------------------
# 4. Main program object lists
# ----------------------------------------------------------------------

FILES_M_INT_QCD = bitset.o pdf.o NNPDFDriver.o mint_module.o ranmar.o HwU.o phase_space.o \
LUPdecompose.o phase_space_gen23.o color_algebra.o math_functions.o \
feynmanrules.o particles.o amplitude_QCD.o matrix_integrate_QCD.o common.o \
phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
read_process_file.o multichannel.o handling_processes.o simple_integrator.o \
helper_modules.o amplitude_library.o command_line_parser.o mg_checks.o scales.o \
pdf_lhapdf62.o

FILES_M_RWGT_QCD = bitset.o color_algebra.o math_functions.o feynmanrules.o particles.o \
amplitude_QCD.o matrix_reweight_QCD.o ranmar.o

FILES_M_UNWGT_QCD = color_algebra.o math_functions.o feynmanrules.o particles.o \
amplitude_QCD.o matrix_unweight_QCD.o

FILES_M_COMBINE_QCD = color_algebra.o math_functions.o feynmanrules.o particles.o \
amplitude_QCD.o matrix_combine_QCD.o

# ----------------------------------------------------------------------
# 5. Build executables
# ----------------------------------------------------------------------

matrix_integrate_QCD: cleanlib $(FILES_M_INT_QCD) dummy.o
	$(FC) $(FFLAGS) -o $@ $(FILES_M_INT_QCD) dummy.o `lhapdf-config --ldflags` -lstdc++ -lc++

matrix_integrate_QCD_library: $(FILES_M_INT_QCD) amplib.o $(AMPLIBS)
	$(FC) $(FFLAGS) -o matrix_integrate_QCD $(FILES_M_INT_QCD) amplib.o $(AMPLIBS) \
	`lhapdf-config --ldflags` -lstdc++ -Wl,-rpath,$(PWD)

matrix_reweight_QCD: $(FILES_M_RWGT_QCD)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_RWGT_QCD)

matrix_unweight_QCD: $(FILES_M_UNWGT_QCD)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_UNWGT_QCD) `lhapdf-config --ldflags` -lstdc++

matrix_combine_QCD: $(FILES_M_COMBINE_QCD)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_COMBINE_QCD)

# ----------------------------------------------------------------------
# 6. Manual dependency rules
# ----------------------------------------------------------------------

matrix_reweight_QCD.o : amplitude_QCD.o math_functions.o particles.o
ranmar.o : mint_module.o
phase_space_gen23.o : phase_space.o LUPdecompose.o particles.o
phase_space_genpt.o : phase_space.o particles.o
phase_space.o : particles.o
haag.o : phase_space.o
amplitude_QCD.o : bitset.o math_functions.o feynmanrules.o color_algebra.o particles.o
matrix_integrate_QCD.o : amplitude_QCD.o phase_space_gen23.o mint_module.o common.o math_functions.o \
	particles.o phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
	read_process_file.o multichannel.o handling_processes.o simple_integrator.o amplitude_library.o \
	command_line_parser.o mg_checks.o scales.o
common.o : particles.o simple_integrator.o
handling_events.o : common.o mint_module.o handling_processes.o simple_integrator.o
read_process_file.o : mint_module.o phase_space_gen23.o cuts.o handling_processes.o simple_integrator.o
multichannel.o : handling_processes.o mint_module.o math_functions.o simple_integrator.o
handling_processes.o : math_functions.o common.o phase_space.o amplitude_QCD.o
cuts.o : common.o particles.o handling_processes.o
pdf_wrap.o : handling_processes.o
simple_integrator.o : helper_modules.o
amplitude_library.o : handling_processes.o read_process_file.o
mg_checks.o : common.o amplitude_QCD.o command_line_parser.o handling_processes.o
scales.o : common.o particles.o cuts.o
amplib.o: $(notdir $(AMPSRC:.f03=.o))

# ----------------------------------------------------------------------
# 7. Cleanup
# ----------------------------------------------------------------------

clean:
	rm -f *.o *.mod library/amp*.f03 library/amp*.data library/amplitudes.bin lib*.so

cleanlib:
	rm -f libamp*.so amp*lib.o amp*lib.mod
	$(FC) $(FFLAGS) -c library/dummy.f03
