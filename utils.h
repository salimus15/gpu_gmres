#ifndef __UTILS_H
#define __UTILS_H

//#include <cusp/gallery/poisson.h>
//#include <cusp/krylov/gmres.h>
//#include <cusp/io/matrix_market.h>
#include <cusp/csr_matrix.h>
//#include <cusp/monitor.h>
//#include <cusp/array1d.h>
#include <cusp/array2d.h>
#include <cusp/blas.h>
#include <cusp/multiply.h>
#include <cusp/linear_operator.h>
#include <iostream> 


#ifdef __cplusplus
extern "C"{
#endif


// where to perform the computation
typedef cusp::device_memory MemorySpace;
typedef cusp::host_memory LocalSpace;

typedef float ValueType;
typedef cusp::array1d<ValueType, cusp::host_memory> CuspArray;
void rotationplan(ValueType& dx, ValueType& dy, ValueType& cs, ValueType& sn);

void genererrotaionplan(ValueType& dx, ValueType& dy, ValueType& cs, ValueType& sn);

void applyrotationplan(cusp::array2d<ValueType, LocalSpace, cusp::column_major>& H, CuspArray& cs, CuspArray& sn, CuspArray& s, int i);

int my_GMRES(cusp::csr_matrix<IndexType, ValueType, MemorySpace>& A, CudaVector& x, CudaVector& b, int restart, cusp::default_monitor<ValueType>& monitor);
//	       Preconditioner& M)

#ifdef __cplusplus
}
#endif
#endif

