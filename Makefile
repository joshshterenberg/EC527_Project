# Flags that can be overridden at make time
CXXFLAGS=-Wall -g -O1
EIGEN_INCLUDE=/usr/include/eigen3
# Flags that are needed for correct compilation
ALL_CXXFLAGS=-std=c++17 -I./cmssw_include -I$(EIGEN_INCLUDE) $(CXXFLAGS)
# Enable the ALADDIN libmd5 interface, used by cmssw
ALL_CPPFLAGS=-DLIBMD_MD5_ALADDIN $(CPPFLAGS)
# Linked dependencies
LDLIBS=-lm -lrt -lcudart -ltinyxml2 -ltbb -lmd -luuid
LINK.o=$(LINK.cc)

# C++ Executable outputs
CXXBINS=\
     test_vertex_fitter \
     # End of CXXBINS

CMSSW_CPP_SOURCES=$(wildcard cmssw_include/HeterogeneousCore/CUDAUtilities/src/*.cc cmssw_include/FWCore/MessageLogger/src/*.cc  cmssw_include/FWCore/Utilities/src/*.cc)
CMSSW_CPP_OBJS=$(CMSSW_CPP_SOURCES:.cc=.o)

.PHONY: all clean check
all: $(CXXBINS)

# Default compilation rule for C++ source files
%.o: %.cc
	$(CXX) -c $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $< -o $@

$(CMSSW_CPP_OBJS): %.o : %.cc
	$(CXX) -c $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $< -o $@

$(CXXBINS): % : %.o $(CMSSW_CPP_SOURCES:.cc=.o)

# Run any self-tests we have
check: test_vertex_fitter
	./test_vertex_fitter

clean:
	rm -f $(CXXBINS) $(CXXBINS:=.o) $(CMSSW_CPP_OBJS)
