# Tested with 
# gcc/4.3.2
# hdf5/1.6.10
# ompi/1.4.2-gnu4.3.2


#extra
export USE_OPENMP=0

ACML_DIR=/home/siddhart/opt/acml_install_4_4_0/gfortran64

export LIBS += -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -L/usr/lib64 -lsz -lz -L$(ACML_DIR)/lib -lacml -lacml_mv /usr/lib64/libppl_c.so.2 -L/sw/global/compilers/gcc/4.8.2/lib64/ -pthread -L/sw/taurus/libraries/bullxmpi/1.2.4.3/lib/ -lmpi -lmpi_cxx -ldl -lm -lnuma -lgfortran -Wl,--export-dynamic -lrt -lnsl -lutil

export ADD_LIBS += -L$(TOP_DIR)/CBLAS/lib -lcblas_LINUX -L$(ACML_DIR) -lacml -lgfortran $(PAPI_POST_LINK_OPTS)

export INC_PATH += -I$(HDF5_DIR)/include -I$(TOP_DIR)/CBLAS/include -I$(ACML_DIR)/include -I/sw/taurus/libraries/bullxmpi/1.2.4.3/include

#export ADDITIONAL_TARGETS = CBLAS_target

export BOOST_ROOT=$(TOP_DIR)

## For UH-shmem, oshcc/oshcxx is used for C/C++  
export F77 = oshfort -g  
export CXX = oshcxx -std=c++0x -g -DPE_ON_FPGA=$(FPGA_PE)
export CC = oshcc -g 

## For mellanox, shmemcc is used for both c and c++. There is no wrapper for Fortran
#export F77 = shmemf77 -g 
#export CXX = shmemcc -std=c++0x -g  
#export CC =  shmemcc -g 
#
export LUACXX = $(CXX)
