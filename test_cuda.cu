
#include <iostream>
#include "gmres_cuda.h"

using namespace std;

int main(int argc, char ** argv){
        int i;
	std::string filename;
	int tolerance, mGmres;
	
	if(argc > 10){
		printf("\nje sais pas trop !!!!");
		return 1;
	}
		
	for(i = 0; i < argc; ++i){
		// we check if the matrix is contained in a matrix market file 
		if (strcmp(argv[i], " --matrix-from-file") == 0){
			// we get the name of the file from where to get the matrix 
			//filename.assign("./rdb968.mtx");
			filename.assign(argv[i+1];
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

        if(cusp_GMRES(filename, tolerance, mGmres)== 0){
          cout << " it seems to be ok\n";
        }



	return 0;
}


