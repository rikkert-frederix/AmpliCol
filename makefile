.DEFAULT_GOAL := amplicol_generate

FC   = gfortran
CC   = gcc
CXX  ?= g++
NVCC ?= nvcc

FP ?= d
ifeq ($(FP),f)
  FPFLAG := -DAC_FP_SINGLE
else ifeq ($(FP),d)
  FPFLAG :=
else
  $(error FP must be 'd' (double, default) or 'f' (float), got '$(FP)')
endif

# FFLAGS    = -O0 -g
# CFLAGS    = -O0 -g $(FPFLAG)
# NVCCFLAGS = -O0 -std=c++14 -rdc=true -Xcompiler -fPIC -I. -ILibrary -diag-suppress 177 --fmad=false $(FPFLAG)
# FFLAGS = -O1
# CFLAGS = -O1
# NVCCFLAGS = -O1 -std=c++14 -rdc=true -Xcompiler -fPIC -I. -ILibrary -diag-suppress 177
# FFLAGS = -O2
# CFLAGS = -O2
# NVCCFLAGS = -O2 -std=c++14 -rdc=true -Xcompiler -fPIC -I. -ILibrary -diag-suppress 177
# FFLAGS = -O3
# CFLAGS = -O3
# NVCCFLAGS = -O3 -std=c++14 -rdc=true -Xcompiler -fPIC -I. -ILibrary -diag-suppress 177
FFLAGS = -O3 -ffast-math
CFLAGS = -O3 -ffast-math
NVCCFLAGS = -O3 -std=c++14 -rdc=true -Xcompiler -fPIC -I. -ILibrary -diag-suppress 177 --use_fast_math


ifeq ($(shell $(CXX) --version | grep -c clang),1)
  STDLIB_FLAG   = -stdlib=libc++
  STDLIB_LDLIBS = -lc++
else
  STDLIB_FLAG   =
  STDLIB_LDLIBS =
endif

#LHAPDF_CFLAGS := $(shell lhapdf-config --cflags)

# ----------------------------------------------------------------------
# 1. Detect amplitude sources and group them
# ----------------------------------------------------------------------

# All amplitude source files: amp<GROUP>_<ID>_lib.f03
AMPSRC    := $(shell find Library/ -name 'amp*_lib.f03')
AMPSRCC   := $(shell find Library/ -name 'amp*_libc.c')
AMPSRCCU  := $(shell find Library/ -name 'amp*_libcu.cu')

# All amplitude object files
AMPOBJ    := $(notdir $(AMPSRC:.f03=.o))
AMPOBJC   := $(notdir $(AMPSRCC:.c=.o))
AMPOBJCU  := $(notdir $(AMPSRCCU:.cu=.o))

# Explicit rule so Make knows how to build amplitude objects
$(AMPOBJ): %.o : Library/%.f03
	$(FC) -frecursive $(FFLAGS) -fPIC -c -I. -ILibrary $<
$(AMPOBJC): %.o : Library/%.c
	$(CC) $(CFLAGS) -fPIC -c -I. -ILibrary $<
$(AMPOBJCU): %.o : Library/%.cu
	$(NVCC) $(NVCCFLAGS) -dc $< -o $@

# Extract GROUP names (amp32, amp81, amp252, ...)
AMPGROUPS    := $(sort $(foreach f,$(AMPOBJ),$(word 1,$(subst _, ,$(notdir $(basename $f))))))
AMPGROUPSC := $(sort $(foreach f,$(AMPOBJC),$(word 1,$(subst _, ,$(notdir $(basename $f))))))
AMPGROUPSCU := $(sort $(foreach f,$(AMPOBJCU),$(word 1,$(subst _, ,$(notdir $(basename $f))))))

# One shared library per group
AMPLIBS    := $(foreach g,$(AMPGROUPS),lib$(g).so)
AMPLIBSC := $(foreach g,$(AMPGROUPSC),lib$(g)_c.so)
AMPLIBSCU := $(foreach g,$(AMPGROUPSCU),lib$(g)_cu.a)

# For MadSpace interface
AMPSPACELIBF   := libamplicolmadspace_f.so
AMPSPACELIBC   := libamplicolmadspace_c.so
AMPSPACELIBCU  := libamplicolmadspace_cu.so

$(AMPSPACELIBF): amplib.o umami_impl.o umami.o $(AMPLIBS) feynmanrules.o
	$(FC) -frecursive $(FFLAGS) -fPIC -shared -o $@ $^

$(AMPSPACELIBC): amplibc.o umamic.o $(AMPLIBSC) FeynmanRulesC.o
	$(CXX) $(CFLAGS) $(STDLIB_FLAG) -fPIC -shared -o $@ $^

FeynmanRulesC.o: FeynmanRules.c FeynmanRules.h AmpliColTypes.h
	$(CC) $(CFLAGS) -fPIC -c FeynmanRules.c -o FeynmanRulesC.o

FeynmanRules_device.o: FeynmanRules_device.cu
	$(NVCC) $(NVCCFLAGS) -dc $< -o $@

amplibcu.o: Library/amplibcu.cu
	$(NVCC) $(NVCCFLAGS) -dc $< -o $@

amplibcu_dlink.o: amplibcu.o $(AMPLIBSCU) FeynmanRules_device.o umamicu.o
	$(NVCC) -dlink -Xcompiler -fPIC $^ -o $@

umamicu.o: umamicu.cu
	$(NVCC) $(NVCCFLAGS) -dc $< -o $@

# $(AMPSPACELIBCU): amplibcu.o amplibcu_dlink.o $(AMPLIBSCU) FeynmanRules_device.o
# 	$(NVCC) -shared -Xcompiler -fPIC --device-link=false -o $@ $^

$(AMPSPACELIBCU): amplibcu.o amplibcu_dlink.o $(AMPOBJCU) FeynmanRules_device.o umamicu.o
	$(CXX) $(CFLAGS) -fPIC -shared -o $@ $^ -L$(dir $(shell which $(NVCC)))../lib64 -lcudart

# ----------------------------------------------------------------------
# 2. Generic compilation rules
# ----------------------------------------------------------------------

%.o: %.f03
	$(FC) $(FFLAGS) -fPIC -c -I. $<

%.o: %.cpp
	$(CXX) $(CFLAGS) $(STDLIB_FLAG) -fPIC -c $< -std=c++17

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
	$(FC) -frecursive $(FFLAGS) -fPIC -c -I. -ILibrary $<

%.o: Library/%.cpp
	$(CXX) $(CFLAGS) $(STDLIB_FLAG) -fPIC -c -I. -ILibrary $< -std=c++17

%.o: Library/%.c
	$(CC) $(CFLAGS) -fPIC -c -I. -ILibrary $<

%.o: Library/%.cu
	$(NVCC) $(NVCCFLAGS) -dc $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -fPIC -c $<

# ----------------------------------------------------------------------
# 3. Build one shared library per amplitude group
# ----------------------------------------------------------------------

define one_lib_template
lib$(1).so: feynmanrules.o $(filter $(1)_%_lib.o,$(AMPOBJ))
	$$(FC) -shared -o $$@ $$^
endef

define one_lib_template_c
lib$(1)_c.so: FeynmanRulesC.o $(filter $(1)_%_libc.o,$(AMPOBJC))
	$$(CC) -shared -o $$@ $$^
endef

define one_lib_template_cu
lib$(1)_cu.a: $(1)_libcu.o $(filter $(1)_%_libcu.o,$(AMPOBJCU))
	$$(NVCC) -lib -o $$@ $$^
endef

$(foreach g,$(AMPGROUPS),$(eval $(call one_lib_template,$(g))))
$(foreach g,$(AMPGROUPSC),$(eval $(call one_lib_template_c,$(g))))
$(foreach g,$(AMPGROUPSCU),$(eval $(call one_lib_template_cu,$(g))))

# ----------------------------------------------------------------------
# 4. Main program object lists
# ----------------------------------------------------------------------

FILES_M_INT_QCD = bitset.o pdf.o NNPDFDriver.o ranmar.o phase_space.o \
LUPdecompose.o phase_space_gen23.o color_algebra.o math_functions.o \
feynmanrules.o FeynmanRulesC.o particles.o amplitude_QCD.o amplicol_generate.o common.o \
phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
read_process_file.o multichannel.o handling_processes.o simple_integrator.o \
helper_modules.o amplitude_library.o command_line_parser.o mg_checks.o scales.o \
pdf_lhapdf62.o

FILES_M_RWGT_QCD = bitset.o color_algebra.o math_functions.o feynmanrules.o FeynmanRulesC.o particles.o \
amplitude_QCD.o amplicol_reweight.o ranmar.o

# ----------------------------------------------------------------------
# 5. Build executables
# ----------------------------------------------------------------------

amplicol_generate: cleanlib $(FILES_M_INT_QCD) dummy.o
	$(FC) $(FFLAGS) -o $@ $(FILES_M_INT_QCD) $(STDLIB_LDLIBS) dummy.o `lhapdf-config --ldflags` -lstdc++ 

amplicol_generate_library_f: $(FILES_M_INT_QCD) amplib.o $(AMPLIBS) $(AMPSPACELIBF)
	$(FC) -frecursive $(FFLAGS) $(STDLIB_LDLIBS) -o amplicol_generate $(FILES_M_INT_QCD) amplib.o $(AMPLIBS) \
	`lhapdf-config --ldflags` -lstdc++ -Wl,-rpath,$(PWD)

amplicol_generate_library_c: $(FILES_M_INT_QCD) amplib.o $(AMPLIBS) $(AMPSPACELIB) amplibc.o $(AMPLIBSC) $(AMPSPACELIBC)
	$(CXX) $(CFLAGS) $(STDLIB_LDLIBS) -o amplicol_generate $(FILES_M_INT_QCD) amplib.o $(AMPLIBS) amplibc.o $(AMPLIBSC) \
	-lstdc++ -lgfortran -Wl,-rpath,$(PWD)

amplicol_reweight: $(FILES_M_RWGT_QCD)
	$(FC) $(FFLAGS) -o $@ $(FILES_M_RWGT_QCD)

# ----------------------------------------------------------------------
# 6. From-scratch builds for the MadSpace interface libraries
# ----------------------------------------------------------------------

.PHONY: madspace_lib madspace_lib_c madspace_lib_cu

madspace_lib:
	$(MAKE) $(AMPSPACELIBF)

madspace_lib_c:
	rm -f $(AMPSPACELIBC) $(AMPLIBSC) $(AMPOBJC) amplibc.o umamic.o FeynmanRulesC.o
	$(MAKE) $(AMPSPACELIBC)

madspace_lib_cu:
	$(MAKE) $(AMPSPACELIBCU)

# ----------------------------------------------------------------------
# 7. Manual dependency rules
# ----------------------------------------------------------------------

amplicol_reweight.o : amplitude_QCD.o math_functions.o particles.o
phase_space_gen23.o : phase_space.o LUPdecompose.o particles.o
phase_space_genpt.o : phase_space.o particles.o
phase_space.o : particles.o
phase_space_haag.o : phase_space.o
amplitude_QCD.o : bitset.o math_functions.o feynmanrules.o FeynmanRulesC.o color_algebra.o particles.o
amplicol_generate.o : amplitude_QCD.o phase_space_gen23.o common.o math_functions.o \
	particles.o phase_space_genpt.o phase_space_haag.o cuts.o pdf_wrap.o handling_events.o \
	read_process_file.o multichannel.o handling_processes.o simple_integrator.o amplitude_library.o \
	command_line_parser.o mg_checks.o scales.o
common.o : particles.o simple_integrator.o
handling_events.o : common.o handling_processes.o simple_integrator.o
read_process_file.o : phase_space_gen23.o cuts.o handling_processes.o simple_integrator.o
multichannel.o : handling_processes.o math_functions.o simple_integrator.o
handling_processes.o : math_functions.o common.o phase_space.o amplitude_QCD.o
cuts.o : common.o particles.o handling_processes.o
pdf_wrap.o : handling_processes.o
simple_integrator.o : helper_modules.o
amplitude_library.o : handling_processes.o read_process_file.o common.o
mg_checks.o : common.o amplitude_QCD.o command_line_parser.o handling_processes.o
scales.o : common.o particles.o cuts.o
amplib.o: $(notdir $(AMPSRC:.f03=.o))
amplibc.o: $(notdir $(AMPSRCC:.c=.o))
umami.o : amplib.o
umamic.o : amplibc.o
umamicu.o : amplibcu.o Library/amplib.cuh

# ----------------------------------------------------------------------
# 8. Cleanup
# ----------------------------------------------------------------------

clean:
	rm -f *.o *.mod Library/amp*.f03 Library/amp*.cu Library/amp*.cuh Library/amp*.c Library/amp*.h Library/amp*.data Library/amplitudes*.bin lib*.so lib*.a amplicol_generate amplicol_reweight madspace.json

.PHONY: clean_madspace_lib clean_madspace_lib_c clean_madspace_lib_cu
clean_madspace_lib:
	rm -f $(AMPSPACELIBF) $(AMPLIBS) $(AMPOBJ) amplib.o umami_impl.o umami.o

clean_madspace_lib_c:
	rm -f $(AMPSPACELIBC) $(AMPLIBSC) $(AMPOBJC) amplibc.o umamic.o

clean_madspace_lib_cu:
	rm -f $(AMPSPACELIBCU) $(AMPLIBSCU) $(AMPOBJCU) amplibcu.o amplibcu_dlink.o umamicu.o
# 	rm -f $(AMPSPACELIBCU) $(AMPOBJCU) FeynmanRules_device.o amplibcu.o amplibcu_dlink.o

cleanlib:
	rm -f libamp*.so amp*lib.o amp*lib.mod Library/amp*.f03 Library/amp*.cpp Library/amp*.hpp Library/amp*.c Library/amp*.h Library/amp*.data Library/amplitudes*.bin madspace.json
	$(FC) $(FFLAGS) -c Library/dummy.f03
