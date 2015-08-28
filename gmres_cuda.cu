#include "gmres_cuda.h"
extern "C" {

// reading a matrix from a matrix market file 
int read_Operator_A_mm(CudaMatrix& mtx, const std::string& filename){
	std::cout << " Going to make read of the matrix \n";
	cusp::io::read_matrix_market_file(mtx, "rdb968.mtx");
	std::cout << " Matrix reading done \n";
	return 0;
}


//
int initialize_problem(CudaMatrix& mtx, const std::string& filename, CudaVector& b, CudaVector& x, int& mGmres, int& tolerance){
	//cusp::csr_matrix<int, double, cusp::device_memory> A;
	// allocate storage for solution (x) and right hand side (b)
	//cusp::array1d<ValueType, MemorySpace> x(A.num_rows, ValueType(1));
	//cusp::array1d<ValueType, MemorySpace> b(A.num_rows);
	
	read_Operator_A_mm( mtx, filename);
	std::cout << " Matrix read and has : " << mtx.num_rows << "rows " << mtx.num_cols << "cols " << mtx.num_entries << " entries \n";
	// here we gonna set the vectors sizes
	x.resize(mtx.num_rows);
	b.resize(mtx.num_rows);	

	// set initial guess
	thrust::fill( x.begin(), x.end(), ValueType(1) );	
	std::cout << " vector x set to size of : " << x.size() << "\n";
	thrust::fill( b.begin(), b.end(), ValueType(2) );
	std::cout << " vector b set to size of : " << b.size() << "\n";
	// set stopping criteria:
	//  iteration_limit    = 100
	//  relative_tolerance = 1e-6
	//	cusp::verbose_monitor<ValueType> monitor(b, 100, 1e-6);
	//	int restart = 50;
	
//	on initialise le moniteur de convergence
	
	return 0;
}

// calling the GMRES function implemented in CUSP
int call_cusp_GMRES(CudaMatrix& A, CudaVector& x, CudaVector b, int restart){
	 // solve the linear system A * x = b with the GMRES
	 
    cusp::krylov::gmres(A, x, b,restart);

	return 0;
}



//
int cusp_GMRES(int argc, char ** argv){
	int i;
	std::string filename;
	int tolerance, mGmres;
	
	CudaMatrix mtx;
	CudaVector x,b;

	if(argc > 10){
		printf("\nje sais pas trop !!!!");
		return 1;
	}
		
	for(i = 0; i < argc; ++i){
		// we check if the matrix is contained in a matrix market file 
		if (strcmp(argv[i], " --matrix-from-file") == 0){
			// we get the name of the file from where to get the matrix 
			filename.assign("./rdb968.mtx");
		}
	
		// we check if tolerance was specified 
		if (strcmp(argv[i], " --tolerance") == 0){
			// we get the value of the tolerance 
			tolerance = atoi(argv[i+1]);
		}
	
		// we check if number of iterations was specified
		if (strcmp(argv[i], " --restart") == 0){
			// we get the number of iterations before a restart 
			mGmres = atoi(argv[i+1]);
		}
	}
	
	//read_Operator_A_mm(mtx, filename);
	initialize_problem(mtx, filename, b, x, mGmres, tolerance);
	std::cout << "problem initialization done !\n ";
	std::cout << " now follow the data states before callin gmres :\n";
	std::cout << " Matrix read and has : " << mtx.num_rows << "rows " << mtx.num_cols << "cols " << mtx.num_entries << " entries \n";
	std::cout << " vector x set to size of : " << x.size() << "\n";
	std::cout << " vector b set to size of : " << b.size() << "\n";
	cusp::default_monitor<ValueType> monitor(b, 100, 1e-6);
//	call_cusp_GMRES( mtx, x, b, mGmres);
	my_GMRES( mtx, x, b, mGmres, monitor );
	std::cout << " gmres solving done !!!\n";
	return 0;
}


// cusp gmres modified. it runs on one gpu 
// coming a version running on multiple gpu(s) i guess

int my_GMRES(CudaMatrix& A, CudaVector& x,  CudaVector& b, int restart, cusp::default_monitor<ValueType>& monitor)
//	       Preconditioner& M)
{
//      typedef typename LinearOperator::value_type   ValueType;
//      typedef typename LinearOperator::memory_space MemorySpace;
 //     typedef typename norm_type<ValueType>::type NormType;
      // here we check that it's a squar matrix
      assert(A.num_rows == A.num_cols);        // sanity check
      const size_t N = A.num_rows;
      const int R = restart;
      int i, j, k;
      ValueType beta = 0;
      ValueType resid0 = 0;
      cusp::array1d<ValueType,cusp::host_memory> rel_resid(1);
      //allocate workspace
      cusp::array1d<ValueType,MemorySpace> w(N);
      cusp::array1d<ValueType,MemorySpace> V0(N); //Arnoldi matrix pos 0
      cusp::array2d<ValueType,MemorySpace,cusp::column_major> V(N,R+1,ValueType(0.0)); //Arnoldi matrix
      //duplicate copy of s on GPU
      cusp::array1d<ValueType,MemorySpace> sDev(R+1);
      //HOST WORKSPACE
      cusp::array2d<ValueType,cusp::host_memory,cusp::column_major> H(R+1, R); //Hessenberg matrix
      cusp::array1d<ValueType,cusp::host_memory> s(R+1);
      cusp::array1d<ValueType,cusp::host_memory> cs(R);
      cusp::array1d<ValueType,cusp::host_memory> sn(R);
      ValueType b_norm = blas::nrm2(b);
      
      do{
			// compute initial residual and its norm //
			cusp::multiply(A, x, w);                     // V(0) = A*x        //
			blas::axpy(b,w,ValueType(-1));               // V(0) = V(0) - b   //
		//	cusp::multiply(M,w,w);                       // V(0) = M*V(0)     //
			beta = blas::nrm2(w);                        // beta = norm(V(0)) //
			blas::scal(w, ValueType(-1.0/beta));         // V(0) = -V(0)/beta //
			blas::copy(w,V.column(0));
			// save very first residual norm //
			if (monitor.iteration_count()== 0){
			  //resid0 = beta;
		//	  cusp::multiply(M,b,V0);
			  resid0 = blas::nrm2(V0)/b_norm;
			}
			//s = 0 //
			blas::fill(s,ValueType(0.0));
			s[0] = beta;
			i = -1;
	
			do{
			  ++i;
			  ++monitor;
			  
			  //apply preconditioner
			  //can't pass in ref to column in V so need to use copy (w)
			  cusp::multiply(A,w,V0);
			  //V(i+1) = A*w = M*A*V(i)    //
		//	  cusp::multiply(M,V0,w);
			  
			  for (k = 0; k <= i; k++){
				 //  H(k,i) = <V(i+1),V(k)>    //
				 H(k, i) = blas::dotc(w, V.column(k));
				 // V(i+1) -= H(k, i) * V(k)  //
				 blas::axpy(V.column(k),w,-H(k,i));
			  }
			  
			  H(i+1,i) = blas::nrm2(w);   
			  // V(i+1) = V(i+1) / H(i+1, i) //
			  blas::scal(w,ValueType(1.0)/H(i+1,i));
			  blas::copy(w,V.column(i+1));
			  
			  applyrotationplan(H,cs,sn,s,i);
			  
			  rel_resid[0] = abs(s[i+1]) / resid0 + monitor.absolute_tolerance();
			  
			  //check convergence condition
			  //if (rel_resid < monitor.relative_tolerance())
			  if (monitor.finished(rel_resid)){
				 break;
			  }
			}while (i+1 < R && monitor.iteration_count()+1 <= monitor.iteration_limit());
	

			// solve upper triangular system in place //
			for (j = i; j >= 0; j--){
			  s[j] /= H(j,j);
			  //S(0:j) = s(0:j) - s[j] H(0:j,j)
			  for (k = j-1; k >= 0; k--){
				 s[k] -= H(k,j) * s[j];
			  }
			}
	
			// update the solution //
	
			//copy s to gpu 
			blas::copy(s,sDev);
			// x= V(1:N,0:i)*s(0:i)+x //
			for (j = 0; j <= i; j++){
			  // x = x + s[j] * V(j) //
			  blas::axpy(V.column(j),x,s[j]);
			}
		} while (rel_resid[0] >= monitor.tolerance() &&  monitor.iteration_count()+1 <= monitor.iteration_limit());
	 	return 0;
}

}// for the extern 
