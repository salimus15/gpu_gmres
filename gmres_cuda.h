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
#include <iostream>  
#include <vector> 

#ifdef __cplusplus
extern "C" {
#endif


// We define the type CudaMatrix which is in fact a csr_matrix class
typedef struct cusp::csr_matrix CudaMatrix;

// We define a CudaVector which is in fact an array1d class
typedef struct cusp::array1d CudaVector;

// where to perform the computation
typedef cusp::device_memory MemorySpace;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ici travail à faire pour tunner le type des données via la commande d'execution //
// which floating point type to use
typedef float ValueType;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



// reading a matrix from a matrix market file 
int read_Operator_A_mm(CudaMatrix& mtx, const std::string& filename);

//
int initialize_problem(CudaMatrix& mtx, const std::string& filename, Vector& x, cusp::cusp::default_monitor& monitor, int& mGmres, int& tolerance);  

// calling the GMRES function implemented in CUSP
int call_cusp_GMRES(CudaMatrix& A, CudaVector& x, CudaVector b, int restart, cusp::cusp::default_monitor& monitor);



//
int cusp_GMRES(int argc, char ** argv);



#ifdef __cplusplus
}
#endif
#endif
