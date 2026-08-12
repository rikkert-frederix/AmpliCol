.DEFAULT_GOAL := amplicol_generate

.PHONY: test_matrix_elements test_fermi_statistics test_mixed_spinors \
	test_three_quark_line_reweight test_three_quark_line_multichannel \
	test_coupling_order_parser test_coupling_order_sectors test_coupling_order_pruning \
	test_q3_maxas_compaction test_q3_maxas_mixed_compaction \
	test_lower_line_maxas_compaction \
	test_coupling_order_lhe_roundtrip \
	test_find_unique_coupling_sectors test_signed_event test_run_parameters \
	update_matrix_cases update_matrix_goldens

FC = gfortran
#FFLAGS= -fbounds-check -g -ffpe-trap=invalid,zero,overflow,underflow,denormal
#FFLAGS = -ffast-math -O3 -mcmodel=large
FFLAGS = -ffast-math -O3
PYTHON ?= python

CXX ?= g++

ifeq ($(shell $(CXX) --version | grep -c clang),1)
  STDLIB_FLAG = -stdlib=libc++
  STDLIB_LDLIBS   = -lc++
else
  STDLIB_FLAG =
  STDLIB_LDLIBS   =
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

%.o: Library/%.f03
	$(FC) $(FFLAGS) -fPIC -c -I. -ILibrary $<

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
feynmanrules.o run_parameters.o particles.o amplitude_QCD.o amplicol_generate.o common.o \
phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
coupling_orders.o read_process_file.o multichannel.o handling_processes.o simple_integrator.o \
helper_modules.o amplitude_library.o command_line_parser.o mg_checks.o scales.o \
pdf_lhapdf62.o

FILES_M_RWGT_QCD = bitset.o color_algebra.o math_functions.o feynmanrules.o run_parameters.o particles.o \
coupling_orders.o amplitude_QCD.o amplicol_reweight.o ranmar.o

FILES_M_TEST_ME = bitset.o color_algebra.o math_functions.o feynmanrules.o run_parameters.o particles.o \
	amplitude_QCD.o matrix_element_regression.o

FILES_M_TEST_FERMI = bitset.o color_algebra.o math_functions.o feynmanrules.o run_parameters.o particles.o \
	amplitude_QCD.o fermi_statistics_regression.o

FILES_M_TEST_THREE_QUARK = bitset.o color_algebra.o math_functions.o feynmanrules.o run_parameters.o particles.o \
	amplitude_QCD.o three_quark_line_reweight_regression.o

FILES_M_TEST_MIXED_SPINOR = bitset.o color_algebra.o math_functions.o feynmanrules.o run_parameters.o particles.o \
	amplitude_QCD.o mixed_spinor_regression.o

FILES_M_TEST_COUPLING_ORDERS = bitset.o color_algebra.o math_functions.o feynmanrules.o run_parameters.o \
	particles.o amplitude_QCD.o coupling_order_sector_regression.o

FILES_M_TEST_COUPLING_PRUNING = bitset.o color_algebra.o math_functions.o feynmanrules.o run_parameters.o \
	particles.o amplitude_QCD.o coupling_order_pruning_regression.o

FILES_M_TEST_SIGNED_EVENT = ranmar.o helper_modules.o simple_integrator.o run_parameters.o particles.o \
	common.o coupling_orders.o handling_events.o signed_event_regression.o

FILES_M_TEST_UNIQUE_COUPLING = $(filter-out amplicol_generate.o,$(FILES_M_INT_QCD)) dummy.o \
	find_unique_coupling_sector_regression.o

FILES_M_TEST_RUN_PARAMETERS = run_parameters.o particles.o run_parameters_regression.o

# ----------------------------------------------------------------------
# 5. Build executables
# ----------------------------------------------------------------------

amplicol_generate: cleanlib $(FILES_M_INT_QCD) dummy.o
	$(FC) $(FFLAGS) -o $@ $(FILES_M_INT_QCD) $(STDLIB_LDLIBS) dummy.o `lhapdf-config --ldflags` -lstdc++ 

amplicol_generate_library: $(FILES_M_INT_QCD) amplib.o $(AMPLIBS)
	$(FC) $(FFLAGS) $(STDLIB_LDLIBS) -o amplicol_generate $(FILES_M_INT_QCD) amplib.o $(AMPLIBS) \
	`lhapdf-config --ldflags` -lstdc++ -Wl,-rpath,$(PWD)

amplicol_reweight: $(FILES_M_RWGT_QCD)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_RWGT_QCD)

matrix_element_regression: $(FILES_M_TEST_ME)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_ME)

fermi_statistics_regression: $(FILES_M_TEST_FERMI)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_FERMI)

three_quark_line_reweight_regression: $(FILES_M_TEST_THREE_QUARK)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_THREE_QUARK)

mixed_spinor_regression: $(FILES_M_TEST_MIXED_SPINOR)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_MIXED_SPINOR)

coupling_order_sector_regression: $(FILES_M_TEST_COUPLING_ORDERS)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_COUPLING_ORDERS)

coupling_order_pruning_regression: $(FILES_M_TEST_COUPLING_PRUNING)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_COUPLING_PRUNING)

signed_event_regression: $(FILES_M_TEST_SIGNED_EVENT)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_SIGNED_EVENT)

find_unique_coupling_sector_regression: $(FILES_M_TEST_UNIQUE_COUPLING)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_UNIQUE_COUPLING) $(STDLIB_LDLIBS) \
		`lhapdf-config --ldflags` -lstdc++

run_parameters_regression: $(FILES_M_TEST_RUN_PARAMETERS)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_TEST_RUN_PARAMETERS)

matrix_element_regression.o: tests/matrix_elements/matrix_element_regression.f03 amplitude_QCD.o particles.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

fermi_statistics_regression.o: tests/matrix_elements/fermi_statistics_regression.f03 amplitude_QCD.o particles.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

three_quark_line_reweight_regression.o: \
		tests/matrix_elements/three_quark_line_reweight_regression.f03 amplitude_QCD.o particles.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

mixed_spinor_regression.o: \
		tests/matrix_elements/mixed_spinor_regression.f03 amplitude_QCD.o particles.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

coupling_order_sector_regression.o: \
		tests/matrix_elements/coupling_order_sector_regression.f03 amplitude_QCD.o particles.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

coupling_order_pruning_regression.o: \
		tests/matrix_elements/coupling_order_pruning_regression.f03 amplitude_QCD.o particles.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

signed_event_regression.o: tests/signed_event_regression.f03 handling_events.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

find_unique_coupling_sector_regression.o: \
		tests/find_unique_coupling_sector_regression.f03 read_process_file.o handling_processes.o common.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

run_parameters_regression.o: tests/run_parameters_regression.f03 run_parameters.o particles.o
	$(FC) $(FFLAGS) -c -I. $< -o $@

tests/matrix_elements/cases.dat: tests/matrix_elements/generate_matrix_cases.py process_list.py
	$(PYTHON) tests/matrix_elements/generate_matrix_cases.py --output $@

tests/matrix_elements/golden.dat: matrix_element_regression tests/matrix_elements/cases.dat
	./matrix_element_regression --write tests/matrix_elements/cases.dat $@

test_matrix_elements: matrix_element_regression
	./matrix_element_regression --check tests/matrix_elements/cases.dat tests/matrix_elements/golden.dat

test_fermi_statistics: fermi_statistics_regression \
		tests/matrix_elements/run_fermi_library_regression.py \
		tests/matrix_elements/fermi_statistics_library_regression.f03
	./fermi_statistics_regression
	$(PYTHON) tests/matrix_elements/run_fermi_library_regression.py \
		--generator $(CURDIR)/fermi_statistics_regression \
		--compiler "$(FC)" --fflags="$(FFLAGS)"

test_mixed_spinors: mixed_spinor_regression \
		tests/matrix_elements/run_mixed_spinor_library_regression.py \
		tests/matrix_elements/mixed_spinor_library_regression.f03
	./mixed_spinor_regression
	$(PYTHON) tests/matrix_elements/run_mixed_spinor_library_regression.py \
		--generator $(CURDIR)/mixed_spinor_regression \
		--compiler "$(FC)" --fflags="$(FFLAGS)"

test_three_quark_line_reweight: three_quark_line_reweight_regression \
		amplicol_reweight tests/matrix_elements/run_three_quark_line_reweight_regression.py
	./three_quark_line_reweight_regression
	$(PYTHON) tests/matrix_elements/run_three_quark_line_reweight_regression.py \
		--reweighter $(CURDIR)/amplicol_reweight --input-card $(CURDIR)/run_card.dat

test_three_quark_line_multichannel: \
		process_list.py tests/process_list_three_quark_multichannel_regression.py
	$(PYTHON) tests/process_list_three_quark_multichannel_regression.py

test_coupling_order_parser: process_list.py tests/process_list_coupling_orders_regression.py
	$(PYTHON) tests/process_list_coupling_orders_regression.py

test_coupling_order_sectors: coupling_order_sector_regression coupling_order_pruning_regression \
		tests/matrix_elements/run_coupling_order_library_regression.py \
		tests/matrix_elements/coupling_order_library_regression.f03 \
		tests/matrix_elements/run_coupling_order_pruning_library_regression.py \
		tests/matrix_elements/coupling_order_pruning_library_regression.f03 \
		find_unique_coupling_sector_regression amplicol_reweight \
		tests/coupling_order_lhe_roundtrip_regression.py
	./coupling_order_sector_regression
	./coupling_order_pruning_regression
	$(PYTHON) tests/matrix_elements/run_coupling_order_library_regression.py \
		--generator $(CURDIR)/coupling_order_sector_regression \
		--compiler "$(FC)" --fflags="$(FFLAGS)"
	$(PYTHON) tests/matrix_elements/run_coupling_order_pruning_library_regression.py \
		--generator $(CURDIR)/coupling_order_pruning_regression \
		--compiler "$(FC)" --fflags="$(FFLAGS)"
	./find_unique_coupling_sector_regression
	$(PYTHON) tests/coupling_order_lhe_roundtrip_regression.py \
		--reweighter $(CURDIR)/amplicol_reweight --input-card $(CURDIR)/run_card.dat

test_find_unique_coupling_sectors: find_unique_coupling_sector_regression
	./find_unique_coupling_sector_regression

test_coupling_order_pruning: coupling_order_pruning_regression \
		tests/matrix_elements/run_coupling_order_pruning_library_regression.py \
		tests/matrix_elements/coupling_order_pruning_library_regression.f03
	./coupling_order_pruning_regression
	$(PYTHON) tests/matrix_elements/run_coupling_order_pruning_library_regression.py \
		--generator $(CURDIR)/coupling_order_pruning_regression \
		--compiler "$(FC)" --fflags="$(FFLAGS)"

test_q3_maxas_compaction: coupling_order_pruning_regression
	./coupling_order_pruning_regression --maxas-q3

test_q3_maxas_mixed_compaction: coupling_order_pruning_regression
	./coupling_order_pruning_regression --maxas-q3-mixed

test_lower_line_maxas_compaction: coupling_order_pruning_regression
	./coupling_order_pruning_regression --maxas-lower-lines

test_coupling_order_lhe_roundtrip: amplicol_reweight \
		tests/coupling_order_lhe_roundtrip_regression.py
	$(PYTHON) tests/coupling_order_lhe_roundtrip_regression.py \
		--reweighter $(CURDIR)/amplicol_reweight --input-card $(CURDIR)/run_card.dat

test_signed_event: signed_event_regression
	./signed_event_regression

test_run_parameters: run_parameters_regression
	./run_parameters_regression run_card.dat tests/input/custom_run_card.dat \
		tests/input/ignore_final_width_run_card.dat

update_matrix_cases: tests/matrix_elements/generate_matrix_cases.py process_list.py
	$(PYTHON) tests/matrix_elements/generate_matrix_cases.py --output tests/matrix_elements/cases.dat

update_matrix_goldens: matrix_element_regression update_matrix_cases
	./matrix_element_regression --write tests/matrix_elements/cases.dat tests/matrix_elements/golden.dat

# ----------------------------------------------------------------------
# 6. Manual dependency rules
# ----------------------------------------------------------------------

amplicol_reweight.o : amplicol_reweight.f03 amplitude_QCD.o coupling_orders.o math_functions.o particles.o run_parameters.o
	$(FC) $(FFLAGS) -fno-finite-math-only -c -I. $<
phase_space_gen23.o : phase_space.o LUPdecompose.o particles.o run_parameters.o
phase_space_genpt.o : phase_space.o particles.o
phase_space.o : particles.o
phase_space_haag.o : phase_space.o
amplitude_QCD.o : bitset.o math_functions.o feynmanrules.o color_algebra.o particles.o
amplicol_generate.o : amplitude_QCD.o phase_space_gen23.o common.o math_functions.o \
	particles.o phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
	read_process_file.o multichannel.o handling_processes.o simple_integrator.o amplitude_library.o \
	command_line_parser.o coupling_orders.o mg_checks.o scales.o
particles.o : run_parameters.o
run_parameters.o : run_parameters.f03
	$(FC) $(FFLAGS) -fno-finite-math-only -c -I. $<
common.o : particles.o run_parameters.o simple_integrator.o
handling_events.o : common.o coupling_orders.o handling_processes.o simple_integrator.o
read_process_file.o : coupling_orders.o phase_space_gen23.o cuts.o handling_processes.o simple_integrator.o
amplitude_library.o : coupling_orders.o pdf_wrap.o
# The inverse-map guard must distinguish NaN/Inf values.  The finite-math
# part of -ffast-math otherwise makes ieee_is_finite fold to true.
multichannel.o : multichannel.f03 handling_processes.o math_functions.o simple_integrator.o
	$(FC) $(FFLAGS) -fno-finite-math-only -c -I. $<
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
	rm -f *.o *.mod Library/amp*.f03 Library/amp*.data Library/amplitudes*.bin lib*.so

cleanlib:
	rm -f libamp*.so amp*lib.o amp*lib.mod Library/amp*.f03 Library/amp*.data Library/amplitudes*.bin
	$(FC) $(FFLAGS) -c Library/dummy.f03
