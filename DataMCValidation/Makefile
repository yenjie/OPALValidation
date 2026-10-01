CXX       ?= g++
ROOTCFG   ?= root-config
CXXFLAGS  += -O2 -Wall -Wextra $(shell $(ROOTCFG) --cflags)
LDLIBS    += $(shell $(ROOTCFG) --libs)

all: CompareDataMC

CompareDataMC: CompareDataMC.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f CompareDataMC

.PHONY: all clean
