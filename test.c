#include <stdlib.h>
#include <stdio.h>
#include "gmres_cuda.h"

int main(int argc, char ** argv){

	if(cusp_GMRES(argc, argv)== 0){
		printf(" it seems to be ok\n");
	}
	
	return 0;
}
