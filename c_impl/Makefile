# Makefile for FeynmanRules library and tests
#
# Targets:
#   all       – build everything (default)
#   test      – build and run the C++ self-consistency tests
#   ref       – build and run the Fortran reference driver (requires gfortran)
#   fulltest  – run ref then test (cross-language diff included)
#   clean     – remove all build artefacts

CXX      := g++
CC       := gcc
FC       := gfortran

FP ?= d
ifeq ($(FP),f)
  FPFLAG := -DAC_FP_SINGLE
else ifeq ($(FP),d)
  FPFLAG :=
else
  $(error FP must be 'd' (double) or 'f' (float))
endif

CXXFLAGS := -std=c++17 -O2 -Wall -Wextra $(FPFLAG)
CFLAGS   := -std=c11 -O2 -Wall -Wextra $(FPFLAG)
FCFLAGS  := -O2

# Sources
LIB_SRC  := FeynmanRules.c
LIB_OBJ  := FeynmanRulesC.o
TEST_SRC := test_feynman_rules.cpp
REF_SRC  := test_feynman_rules_reference.f03
MOD_SRC  := feynmanrules.f03

# Outputs
TEST_BIN := test_feynman_rules
REF_BIN  := fortran_ref
REF_DAT  := fortran_reference.dat
MOD_OBJ  := feynmanrules.o

# ── default ──────────────────────────────────────────────────────────────────
.PHONY: all
all: $(TEST_BIN) $(REF_BIN)

# ── C library object ─────────────────────────────────────────────────────────
$(LIB_OBJ): $(LIB_SRC) FeynmanRules.h AmpliColTypes.h
	$(CC) $(CFLAGS) -c $< -o $@

# ── C++ test binary ──────────────────────────────────────────────────────────
$(TEST_BIN): $(TEST_SRC) $(LIB_OBJ) FeynmanRules.h AmpliColTypes.h
	$(CXX) $(CXXFLAGS) $(TEST_SRC) $(LIB_OBJ) -o $@

# ── Fortran module object + .mod file ────────────────────────────────────────
# feynmanrules.f03 must be compiled first so gfortran generates feynmanrules.mod
# before the reference driver tries to `use FeynmanRules`.
$(MOD_OBJ) feynmanrules.mod &: $(MOD_SRC)
	$(FC) $(FCFLAGS) -c $(MOD_SRC) -o $(MOD_OBJ)

# ── Fortran reference binary ─────────────────────────────────────────────────
$(REF_BIN): $(REF_SRC) $(MOD_OBJ) feynmanrules.mod
	$(FC) $(FCFLAGS) $(REF_SRC) $(MOD_OBJ) -o $@

# ── run targets ──────────────────────────────────────────────────────────────
.PHONY: test
test: $(TEST_BIN)
	./$(TEST_BIN)

.PHONY: ref
ref: $(REF_BIN)
	./$(REF_BIN)

.PHONY: fulltest
fulltest: ref test

# ── clean ─────────────────────────────────────────────────────────────────────
.PHONY: clean
clean:
	rm -f $(LIB_OBJ) $(MOD_OBJ) $(TEST_BIN) $(REF_BIN) $(REF_DAT) *.mod
