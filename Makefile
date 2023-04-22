# Flags that can be overridden at make time
CXXFLAGS=-Wall -g -O1
# Flags that are needed for correct compilation
ALL_CXXFLAGS=-std=c++11 $(CXXFLAGS)
# Linked dependencies
LDLIBS=-lm -lrt

# C++ Executable outputs
CXXBINS=\
     test_vertex_fitter \
     # End of CXXBINS

.PHONY: all clean check
all: $(CXXBINS)

# Default compilation rule for C++ source files
.cc.o:
	$(CXX) -c $(CPPFLAGS) $(ALL_CXXFLAGS) $<

# Run any self-tests we have
check: test_vertex_fitter
	./test_vertex_fitter

clean:
	rm -f $(CXXBINS) $(CXXBINS:=.o)
