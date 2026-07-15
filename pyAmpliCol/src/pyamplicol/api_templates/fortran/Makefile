RUSTICOL_CONFIG ?= rusticol-config
FC = gfortran
RUSTICOL_FORTRAN := $(shell $(RUSTICOL_CONFIG) --fortran-module)
FFLAGS ?= -O2
FFLAGS += -std=f2008 -ffree-line-length-none -Wall -Wextra
LDLIBS += $(shell $(RUSTICOL_CONFIG) --libs)

.PHONY: all clean run

all: check_standalone

check_standalone: rusticol.o check_standalone.o
	$(FC) $(FFLAGS) $^ -o $@ $(LDLIBS)

rusticol.o rusticol.mod: $(RUSTICOL_FORTRAN)
	$(FC) $(FFLAGS) -c $< -o rusticol.o

check_standalone.o: check_standalone.f90 rusticol.mod
	$(FC) $(FFLAGS) -c $< -o $@

run: check_standalone
	./check_standalone

clean:
	$(RM) check_standalone check_standalone.o rusticol.o rusticol.mod
