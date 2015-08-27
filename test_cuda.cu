
#include <iostream>
#include "gmres_cuda.h"

using namespace std;

int main(int argc, char ** argv){

        if(cusp_GMRES(argc, argv)== 0){
          cout << " it seems to be ok\n";
        }

	return 0;
}


