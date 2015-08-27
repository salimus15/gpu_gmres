#include "gmres_cuda.h"


// reading a matrix from a matrix market file 
int read_Operator_A_mm(CudaMatrix& mtx, const std::string& filename){

	cusp::io::read_matrix_market_file(mtx, filename);
	return 0;
}


//
int initialize_problem(CudaMatrix& mtx, const std::string& filename, CudaVector& b, CudaVector& x, int& mGmres, int& tolerance){
	//cusp::csr_matrix<int, double, cusp::device_memory> A;
	// allocate storage for solution (x) and right hand side (b)
	//cusp::array1d<ValueType, MemorySpace> x(A.num_rows, ValueType(1));
	//cusp::array1d<ValueType, MemorySpace> b(A.num_rows);
	
	read_Operator_A_mm( mtx, filename);
	// set initial guess
   thrust::fill( x.begin(), x.end(), ValueType(1) );	
	// set stopping criteria:
	//  iteration_limit    = 100
	//  relative_tolerance = 1e-6
	//	cusp::verbose_monitor<ValueType> monitor(b, 100, 1e-6);
	//	int restart = 50;
	
//	on initialise le moniteur de convergence
	
	return 0;
}

// calling the GMRES function implemented in CUSP
int call_cusp_GMRES(CudaMatrix& A, CudaVector& x, CudaVector b, int restart, cusp::default_monitor<ValueType>& monitor){
	 // solve the linear system A * x = b with the GMRES
    cusp::krylov::gmres(A, x, b,restart, monitor);

	return 0;
}



//
int cusp_GMRES(int argc, char ** argv){
	int i;
	char * filename;
	int tolerance, mGmres;
	
	CudaMatrix mtx;
	CudaVector x,b;

	if(argc < 10){
		printf("\nje sais pas trop !!!!");
		return 1;
	}
		
	for(i = 0; i < argc; ++i){
		// we check if the matrix is contained in a matrix market file 
		if (strcmp(argv[i], " --matrix-from-file") == 0){
			// we get the name of the file from where to get the matrix 
			filename = argv[i+1];
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
	
	read_Operator_A_mm(mtx, filename);
	initialize_problem(mtx, filename, b, x, mGmres, tolerance);
	cusp::default_monitor<ValueType> monitor(b, mGmres, tolerance);
	call_cusp_GMRES( mtx, x, b, mGmres, monitor);
	return 0;
}

