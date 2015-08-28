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
	call_cusp_GMRES( mtx, x, b, 100);
//	my_GMRES( mtx, x, b, 100, monitor );
	std::cout << " gmres solving done !!!\n";
	return 0;
}



}// for the extern 
