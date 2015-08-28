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
void rotationplan(CuspArray& dx, CuspArray& dy, CuspArray& cs, CuspArray& sn);

void genererrotaionplan(CuspArray& dx, CuspArray& dy, CuspArray& cs, CuspArray& sn);

void applyrotationplan(cusp::array2d<ValueType, LocalSpace, cusp::column_major>& H, CuspArray& cs, CuspArray& sn, CuspArray& s, int i);


#ifdef __cplusplus
}
#endif
#endif

