.DEFAULT_GOAL := amplicol_generate

.PHONY: test_matrix_elements test_integrated_kernels test_massive_dipole_kernel \
	test_simple_integrator_mixed_dims test_simple_integrator_uncertainty_sampling \
	test_histograms_cli test_integration_histograms test_integration_histograms_physics \
	update_matrix_cases update_matrix_goldens

FC = gfortran
#FFLAGS= -fbounds-check -g -ffpe-trap=invalid,zero,overflow,underflow,denormal
#FFLAGS = -ffast-math -O3 -mcmodel=large
FFLAGS = -ffast-math -O3
PYTHON ?= python

CXX_ORIGIN := $(origin CXX)
UNAME_S := $(shell uname -s)
ifeq ($(CXX_ORIGIN),default)
  ifeq ($(UNAME_S),Darwin)
    CXX = clang++
  else
    CXX = g++
  endif
endif

ifeq ($(shell $(CXX) --version | grep -c clang),1)
  STDLIB_FLAG = -stdlib=libc++
  STDLIB_LDLIBS   = -lc++
else
  STDLIB_FLAG =
  STDLIB_LDLIBS   = -lstdc++
endif

LHAPDF_CFLAGS  := $(shell lhapdf-config --cflags)


# ----------------------------------------------------------------------
# 1. Detect amplitude sources and group them
# ----------------------------------------------------------------------

# All amplitude source files: amp<GROUP>_<ID>_lib.f03
AMPSRC := $(shell find Library/ -name 'amp*_lib.f03')

# All amplitude object files
AMPOBJ := $(notdir $(AMPSRC:.f03=.o))

# Explicit rule so Make knows how to build amplitude objects
$(AMPOBJ): %.o : Library/%.f03
	$(FC) $(FFLAGS) -fPIC -c -I. -ILibrary $<

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
	$(CXX) $(CXXFLAGS) $(STDLIB_FLAG) $(LHAPDF_CFLAGS)  -c -I. -IPDF $< -std=c++11

%.o: PDF/%.f90
	$(FC) $(FFLAGS) -c -I. -IPDF $<

%.o: PDF/%.f03
	$(FC) $(FFLAGS) -c -I. -IPDF $<

%.o: SimpleIntegrator/%.f
	$(FC) $(FFLAGS) -c -I. -ISimpleIntegrator $<

%.o: SimpleIntegrator/%.f03
	$(FC) $(FFLAGS) -c -I. -ISimpleIntegrator $<

%.o: PhaseSpace/%.f90
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace $<

%.o: PhaseSpace/%.f03
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace $<

%.o: CS_Dipoles/%.f03
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace $<

%.o: CS_Dipoles/%.f90
	$(FC) $(FFLAGS) -c -I. -IPhaseSpace $<

%.o: Library/%.f03
	$(FC) $(FFLAGS) -fPIC -c -I. -ILibrary $<

%.o: Utilities/plot_events/%.f03
	$(FC) $(FFLAGS) -c -I. -IUtilities/plot_events $<

# ----------------------------------------------------------------------
# 3. Build one shared library per amplitude group
# ----------------------------------------------------------------------

define one_lib_template
lib$(1).so: $(filter $(1)_%_lib.o,$(AMPOBJ))
	$$(FC) feynmanrules.o -shared -o $$@ $$^
endef

$(foreach g,$(AMPGROUPS),$(eval $(call one_lib_template,$(g))))

# ----------------------------------------------------------------------
# 4. Main program object lists
# ----------------------------------------------------------------------

FILES_M_INT_QCD = bitset.o pdf.o NNPDFDriver.o ranmar.o phase_space.o \
LUPdecompose.o phase_space_gen23.o color_algebra.o math_functions.o \
feynmanrules.o particles.o amplitude_QCD.o amplicol_generate.o common.o \
phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
read_process_file.o multichannel.o handling_processes.o 	\
simple_integrator.o helper_modules.o amplitude_library.o command_line_parser.o \
mg_checks.o scales.o pdf_lhapdf62.o phase_space_module.o cs_dipole_mappings.o \
cs_lc_dipoles.o cs_integrated_kernels.o subtraction.o integrated_dipoles.o
FILES_M_INT_QCD += integration_histograms.o integration_analysis.o

FILES_M_RWGT_QCD = bitset.o color_algebra.o math_functions.o feynmanrules.o particles.o \
amplitude_QCD.o amplicol_reweight.o ranmar.o

FILES_M_TEST_ME = bitset.o color_algebra.o math_functions.o feynmanrules.o particles.o \
amplitude_QCD.o matrix_element_regression.o

test_integrated_kernels: cs_dipole_mappings.o cs_integrated_kernels.o
	$(FC) $(FFLAGS) -I. -IPhaseSpace -o tests/integrated_kernels.exe \
		tests/integrated_kernels.f90 cs_integrated_kernels.o cs_dipole_mappings.o
	./tests/integrated_kernels.exe

test_massive_dipole_kernel: cs_dipole_mappings.o cs_lc_dipoles.o
	$(FC) $(FFLAGS) -I. -IPhaseSpace -o tests/massive_dipole_kernel.exe \
		tests/massive_dipole_kernel.f90 cs_lc_dipoles.o cs_dipole_mappings.o
	./tests/massive_dipole_kernel.exe

test_simple_integrator_mixed_dims: simple_integrator.o helper_modules.o ranmar.o
	$(FC) $(FFLAGS) -I. -o tests/simple_integrator_mixed_dims.exe \
		tests/simple_integrator_mixed_dims.f90 simple_integrator.o helper_modules.o ranmar.o
	./tests/simple_integrator_mixed_dims.exe

test_simple_integrator_uncertainty_sampling: simple_integrator.o helper_modules.o ranmar.o
	$(FC) $(FFLAGS) -I. -o tests/simple_integrator_uncertainty_sampling.exe \
		tests/simple_integrator_uncertainty_sampling.f90 simple_integrator.o helper_modules.o ranmar.o
	./tests/simple_integrator_uncertainty_sampling.exe

test_histograms_cli:
	python3 tests/histograms_cli.py

test_integration_histograms: integration_histograms.o
	$(FC) $(FFLAGS) -I. -IUtilities/plot_events -o tests/integration_histograms.exe \
		tests/integration_histograms.f90 integration_histograms.o
	./tests/integration_histograms.exe
	python3 Utilities/plot_events/internal/histograms.py tests/integration_histograms.HwU \
		--gnuplot --out=tests/integration_histograms_plot --no_open

test_integration_histograms_physics: amplicol_generate
	python3 tests/integration_histograms_physics.py

# ----------------------------------------------------------------------
# 5. Build executables
# ----------------------------------------------------------------------

amplicol_generate: cleanlib $(FILES_M_INT_QCD) dummy.o
	$(FC) $(FFLAGS) -o $@ $(FILES_M_INT_QCD) dummy.o `lhapdf-config --ldflags` $(STDLIB_LDLIBS)

amplicol_generate_library: $(FILES_M_INT_QCD) amplib.o $(AMPLIBS)
	$(FC) $(FFLAGS) $(STDLIB_LDLIBS) -o amplicol_generate $(FILES_M_INT_QCD) amplib.o $(AMPLIBS) \
	`lhapdf-config --ldflags` $(STDLIB_LDLIBS) -Wl,-rpath,$(PWD)

amplicol_reweight: $(FILES_M_RWGT_QCD)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_RWGT_QCD)

matrix_element_regression: $(FILES_M_TEST_ME)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_ME)

matrix_element_regression.o: tests/matrix_elements/matrix_element_regression.f03 amplitude_QCD.o particles.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

tests/matrix_elements/cases.dat: tests/matrix_elements/generate_matrix_cases.py process_list.py
	$(PYTHON) tests/matrix_elements/generate_matrix_cases.py --output $@

tests/matrix_elements/golden.dat: matrix_element_regression tests/matrix_elements/cases.dat
	./matrix_element_regression --write tests/matrix_elements/cases.dat $@

test_matrix_elements: matrix_element_regression
	./matrix_element_regression --check tests/matrix_elements/cases.dat tests/matrix_elements/golden.dat

update_matrix_cases: tests/matrix_elements/generate_matrix_cases.py process_list.py
	$(PYTHON) tests/matrix_elements/generate_matrix_cases.py --output tests/matrix_elements/cases.dat

update_matrix_goldens: matrix_element_regression update_matrix_cases
	./matrix_element_regression --write tests/matrix_elements/cases.dat tests/matrix_elements/golden.dat

# ----------------------------------------------------------------------
# 6. Manual dependency rules
# ----------------------------------------------------------------------

amplicol_reweight.o : amplitude_QCD.o math_functions.o particles.o
phase_space_gen23.o : phase_space.o LUPdecompose.o particles.o
phase_space_genpt.o : phase_space.o particles.o
phase_space.o : particles.o
phase_space_haag.o : phase_space.o
amplitude_QCD.o : bitset.o math_functions.o feynmanrules.o color_algebra.o particles.o
amplicol_generate.o : amplitude_QCD.o phase_space_gen23.o common.o math_functions.o \
	particles.o phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
	read_process_file.o multichannel.o handling_processes.o simple_integrator.o amplitude_library.o \
	command_line_parser.o mg_checks.o scales.o subtraction.o integrated_dipoles.o phase_space_module.o \
	cs_dipole_mappings.o cs_lc_dipoles.o cs_integrated_kernels.o feynmanrules.o \
	integration_histograms.o integration_analysis.o
integration_analysis.o : integration_histograms.o common.o cuts.o
common.o : particles.o simple_integrator.o
handling_events.o : common.o handling_processes.o simple_integrator.o
read_process_file.o : phase_space_gen23.o cuts.o handling_processes.o simple_integrator.o
multichannel.o : handling_processes.o math_functions.o simple_integrator.o
handling_processes.o : math_functions.o common.o phase_space.o amplitude_QCD.o
cuts.o : common.o particles.o handling_processes.o
pdf_wrap.o : handling_processes.o
simple_integrator.o : helper_modules.o
amplitude_library.o : handling_processes.o read_process_file.o
mg_checks.o : common.o amplitude_QCD.o command_line_parser.o handling_processes.o
scales.o : common.o particles.o cuts.o
amplib.o : $(notdir $(AMPSRC:.f03=.o))
subtraction.o : particles.o handling_processes.o amplitude_QCD.o cuts.o phase_space_module.o \
	cs_dipole_mappings.o cs_lc_dipoles.o feynmanrules.o
cs_lc_dipoles.o : cs_dipole_mappings.o
cs_integrated_kernels.o : cs_dipole_mappings.o
integrated_dipoles.o : handling_processes.o cs_dipole_mappings.o cs_integrated_kernels.o pdf_wrap.o

# ----------------------------------------------------------------------
# 7. Cleanup
# ----------------------------------------------------------------------

clean:
	rm -f *.o *.mod Library/amp*.f03 Library/amp*.data Library/amplitudes*.bin lib*.so

cleanlib:
	rm -f libamp*.so amp*lib.o amp*lib.mod Library/amp*.f03 Library/amp*.data Library/amplitudes*.bin
	$(FC) $(FFLAGS) -c Library/dummy.f03
