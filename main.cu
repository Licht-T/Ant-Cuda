#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

#include <sys/stat.h>
#include <sys/types.h>

// includes CUDA Runtime
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

#include "Constants.h"
#include "DataStructures.h"
#include "Variables.h"
#include "kernel.h"
#include "Display.h"
#include "IO.h"

int main(int argc, char *argv[]){

	getHoming homeOp;
	thrust::plus<int> binary_op;

	IOInit();

	initialize();


	double normalSum = 0.0;

	for(unsigned long long int dummy=1; dummy<=MAX_STEP; dummy++){
		reset(1,0,dummy);

		//cudaDeviceSynchronize();
		//display(argc,argv);

		for(int t=0; t<MAX_TIME; t++){
			calculation();
			if(t%500==0){
				IOEffPoll(0,500,dummy,t);
			}
		}

		normalSum += thrust::transform_reduce(ants_d_ptr, ants_d_ptr+NMAX, homeOp, 0, binary_op);
	}

	for (int n=0; n<=500; n+=50){
		IOEffWrite(0,n,normalSum);
	}

	for (int pw=1; pw<=7; pw++){
		for (int n=0; n<=450; n+=50){
			double sensor = pow(10,-pw);
			int naho = (NMAX - n);

			double sum = 0.0;
			for(unsigned long long int dummy=1; dummy<=MAX_STEP; dummy++){
				reset(sensor,naho,dummy);

				//display(argc,argv);
				for(int t=1; t<=MAX_TIME; t++){
					calculation();
					if(t%500==0){
						IOEffPoll(pw,n,dummy,t);
					}
				}

				sum += thrust::transform_reduce(ants_d_ptr, ants_d_ptr+NMAX, homeOp, 0, binary_op);
			}
			IOCellWrite(pw,n);
			IOEffWrite(pw,n,sum);
		}
		IOEffWrite(pw,500,normalSum);
	}

}

