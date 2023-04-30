# Flags that can be overridden at make time
CXXFLAGS=
CUDART_CFLAGS?=
CUDART_LDFLAGS?=
CUDART_LDLIBS=-lcudart
# Flags that are needed for correct compilation
ALL_CXXFLAGS=-std=c++11 -Wall -g -O3 -march=native $(CUDART_CFLAGS) $(CXXFLAGS)
ALL_CPPFLAGS=$(CPPFLAGS)
# Linked dependencies
LDLIBS=-lm -lrt $(CUDART_LDLIBS) -pthread
LDFLAGS=$(CUDART_LDFLAGS)
LINK.o=$(LINK.cc)

# C++ Executable outputs
CXXBINS=\
     test_vertex_fitter \
     # End of CXXBINS
CUDABINS=\
     test_vertex_fitter_CUDA \
     test_vertex_fitter_CUDAv2 \
     test_vertex_fitter_CUDAv3 \
     # End of CUDABINS

.PHONY: all clean check
all: $(CXXBINS) $(CUDABINS)

# Default compilation rule for C++ source files
%.o: %.cc
	$(CXX) -c $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $< -o $@

$(CUDABINS): %: %.cu
	nvcc -arch sm_35 $< -o $@  

$(CXXBINS): % : %.o

# Run any self-tests we have
check: $(CXXBINS) $(CUDABINS)
	./test_vertex_fitter
	./test_vertex_fitter_CUDA
	./test_vertex_fitter_CUDAv2
	./test_vertex_fitter_CUDAv3

clean:
	rm -f $(CXXBINS) $(CXXBINS:=.o) $(CUDABINS)
