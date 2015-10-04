// -*- mode: c++; -*-

#ifndef BUILDKKRMATRIX_GPU_H
#define BUILDKKRMATRIX_GPU_H

#include <stdlib.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"


// #include "TestStructures.hpp"
#include "Misc/Indices.hpp"
#include "Main/SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"
#include "Misc/Coeficients.hpp"



void *allocateDConst(void);
void freeDConst(void * d_store);

void setupForBuildKKRMatrix_gpu_opaque(LSMSSystemParameters &lsms, AtomData &atom, void *d_const);

void buildKKRMatrix_gpu_opaque(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom, Complex energy, Complex prel, int iie, Matrix<Complex> &m, void *d_const);

#endif
