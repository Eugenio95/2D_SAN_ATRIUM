EXE:=2D_parallel_SAN_atrium
SRC:=$(wildcard *.cpp) $(wildcard *.cu)
INC:=$(wildcard *.hpp)
OBJ:=${SRC:%.cpp=%.o}
OBJ:=${OBJ:%.cu=%.o} # [FIXME] horrible hack
NVCC:=nvcc
#NVCFLAGS+=--gpu-architecture compute_70 --gpu-code=compute_70 -Wno-deprecated-gpu-target
NVCFLAGS+=--gpu-architecture=compute_70 --gpu-code=compute_70
#NVCFLAGS+=-DNO_CUDA_CHECK_ERROR

## compiler flags
CXXFLAGS += -std=c++11 -Wall -Wpedantic
##LDLIBS += -L/usr/local/cuda-10.2/lib64/ -lm -lcudart
##LDLIBS += -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/cuda/11.0/lib64 -lm -lcudart
LDLIBS += -L$(LD_LIBRARY_PATH) -lm -lcudart

.PHONY: test clean

$(EXE): $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

%.o : %.cu
	$(NVCC) $(NVCFLAGS) -c $<

test: $(EXE)
	cd ../m_scripts/Sim_Param_folder/ && ../../CUDA/$(EXE)

clean:
	\rm -f $(OBJ) $(EXE) *~ gmon.out

distclean: clean
