RUSTICOL_CONFIG ?= rusticol-config
CXX ?= c++
CPPFLAGS += $(shell $(RUSTICOL_CONFIG) --cflags)
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -Wall -Wextra -pedantic
LDLIBS += $(shell $(RUSTICOL_CONFIG) --libs)

.PHONY: all clean run

all: check_standalone

check_standalone: check_standalone.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

run: check_standalone
	./check_standalone

clean:
	$(RM) check_standalone
