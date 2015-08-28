#ifndef __UTILS_H
#define __UTILS_H

#ifdef __cplusplus
extern "C"{
#endif

//#include <cusp/gallery/poisson.h>
//#include <cusp/krylov/gmres.h>
//#include <cusp/io/matrix_market.h>
#include <cusp/csr_matrix.h>
//#include <cusp/monitor.h>
//#include <cusp/array1d.h>
#include <cusp/blas.h>
#include <cusp/multiply.h>
#include <cusp/linear_operator.h>
#include <iostream> 

// where to perform the computation
typedef cusp::device_memory MemorySpace;


typedef float ValueType;

void rotationplan(ValueType& dx, ValueType& dy, ValueType& cs, ValueType& sn);

void genererrotaionplan(ValueType& dx, ValueType& dy, ValueType& cs, ValueType& sn);

void applyrotationplan(cusp::csr_matrix<int, ValueType, MemorySpace>& H, ValueType& cs, ValueType& sn, ValueType& s, int i);


#ifdef __cplusplus
}
#endif
#endif

