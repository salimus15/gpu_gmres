// AUTHOR   NAHI SALIM 2015 

#ifndef __GMRES_CUDA_H
#define __GMRES_CUDA_H

#include <cusp/hyb_matrix.h>
//#include <cusp/gallery/poisson.h>
#include <cusp/krylov/gmres.h>
#include <cusp/io/matrix_market.h>
#include <cusp/csr_matrix.h>
#include <cusp/monitor.h>
#include <cusp/array1d.h>
#include <cusp/array2d.h>
#include <cusp/blas.h>
#include <cusp/multiply.h>
#include <cusp/linear_operator.h>
#include <iostream>  
//#include <stdlib.h>	
#include <vector> 
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

// where to perform the computation
typedef cusp::device_memory MemorySpace;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ici travail à faire pour tunner le type des données via la commande d'execution //
// which floating point type to use
//typedef float ValueType;
////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//typedef int IndexType;


// We define the type CudaMatrix which is in fact a csr_matrix class
typedef struct cusp::csr_matrix<IndexType, ValueType, MemorySpace> CudaMatrix;

// We define a CudaVector which is in fact an array1d class
//typedef struct cusp::array1d<ValueType, MemorySpace> CudaVector;







// reading a matrix from a matrix market file 
int read_Operator_A_mm(CudaMatrix& mtx, const std::string & filename);

//
int initialize_problem(CudaMatrix& mtx, const std::string& filename, CudaVector& b, CudaVector& x, int& mGmres, int& tolerance);  

// calling the GMRES function implemented in CUSP
int call_cusp_GMRES(CudaMatrix& A, CudaVector& x, CudaVector b, int restart);

//
int cusp_GMRES(int argc, char ** argv);



#ifdef __cplusplus
}
#endif
#endif
