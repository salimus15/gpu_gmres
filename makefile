ALL:  blib test exec



#NVIDIA = /apps/cuda/cuda-6.5/
CUDA = /apps/cuda/cuda-6.5
CUSP = /apps/cusplibrary
NVIDINCADD = -I$(NVIDIA)/common/inc
CUDAINCADD = -I$(CUDA)/include 
CC =  -lstdc++
GPPFLAGS = -c

#GCCOPT = -O2 -fno-rtti -fno-exceptions 
#INTELOPT = -O3 -fno-rtti -xW -restrict -fno-alias
#SPMV = $(CUSP)/cusp/detail/device/spmv
#DEB = -g
#NVCC = -G
#ARCH = -arch=sm_20
ARCH = -arch=sm_35

blib:
#    nvcc $(DEB) -c $(CC) -lm -I${CUDAINCADD} -I$(CUSP) -I$(SPMV) -I$(THRUST) -lcusparse -I./ gmres_cuda.cu
	nvcc -c -I{CUDAINCADD} -I${CUSP} -I ./ gmres_cuda.cu -o gmres.o
test:
	g++ -c test_cuda.cpp -o test.o	
exec:
	nvcc gmres.o test.o -o runit