# Tested with 
#  - gcc ver. 4.6.7
#  - SHMEM-versions: UH Reference implementation (gasnet-ibv)
#  - Prebuilt HDF5 binary (Linux 2.6 CentOS 6 x86_64, gcc-4.4.7)
#  - ACML installation (version 4.4.0, 64bit, gfortran)


#extra
export USE_OPENMP=0

ACML_DIR=/path/to/acml_install/gfortran64

export LIBS += -L$(HDF5_DIR)/lib -Wl,rpath,$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -L/usr/lib64 -L$(ACML_DIR)/lib -Wl,rpath,$(ACML_DIR)/lib -lacml  /usr/lib64/libppl_c.so.4  -pthread  -ldl -lm -lnuma -lgfortran -Wl,--export-dynamic -lrt -lnsl -lutil

export ADD_LIBS += -L$(TOP_DIR)/CBLAS/lib -lcblas_LINUX -L$(ACML_DIR) -lacml -lgfortran $(PAPI_POST_LINK_OPTS)

export INC_PATH += -I$(HDF5_DIR)/include -I$(TOP_DIR)/CBLAS/include -I$(ACML_DIR)/include -I/opt/openmpi/gnu/1.8.8/include

#export ADDITIONAL_TARGETS = CBLAS_target

export BOOST_ROOT=$(TOP_DIR)

## For UH-shmem reference implementation, the compiler wrappers are named as oshcc/oshcxx/oshfort
export F77 = oshfort   
export CXX = oshcxx  -std=c++0x
export CC = oshcc 


export LUACXX = $(CXX)
