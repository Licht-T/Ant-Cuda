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
#include "kernel7.h"
#include "Display.h"
#include "IO7.h"

int main(int argc, char *argv[]){

    getHoming homeOp;
    thrust::plus<int> binary_op;

    IOInit();

    initialize();

    double normalSum = 0.0;

    int access[MACRO_NUM_FOODS] = {0, 0};
    double prob[5] = {0.0, 0.0, 0.0, 0.0, 0.0}; // 0:どちらにもアクセスしていない, 1:片方にアクセス, 2:両方にアクセス, 3:id=0にアクセスしている状態, 4:id=1にアクセスしている状態
    double probNormal[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double foodspre[MACRO_NUM_FOODS] = {0.0, 0.0};

    /* for(unsigned long long int dummy=1; dummy<=MACRO_MAX_STEP; dummy++){ */
    /*     reset(pow(10,-3),50,dummy); */
    /*     display(argc,argv); */
    /*  */
    /*     for(int t=0; t<MACRO_MAX_TIME; t++){ */
    /*         calculation(); */
    /*     } */
    /* } */

    // 馬鹿ありのいない場合
    for(unsigned long long int dummy=1; dummy<=MACRO_MAX_STEP; dummy++){
        reset(1,0,dummy);

        // display(argc,argv);

        for (int id=0; id<MACRO_NUM_FOODS; id++)
            foodspre[id] = MACRO_FOODSOURCE;
        for(int t=0; t<MACRO_MAX_TIME; t++){
            calculation();
            cudaMemcpyFromSymbol(&foods, foods_d, sizeof(Food)*MACRO_NUM_FOODS);
            for (int id=0; id<MACRO_NUM_FOODS; id++){
                if( foods[id].vol < foodspre[id] + MACRO_REC ){
                    access[id] = 1;
                    probNormal[3+id] += 1;
                }
                else{
                    access[id] = 0;
                }
                foodspre[id] = foods[id].vol;
            }
            probNormal[access[0]+access[1]] += 1;
            // if(t%500==0){
            //     IOEffPoll(0,500,dummy,t);
            // }

        }

        normalSum += thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, homeOp, 0, binary_op);
    }

    for (int n=0; n<=500; n+=50){
        IOEffWrite(0,n,normalSum);
        IOProbWrite(0,n,probNormal);
    }

    // 馬鹿アリの感受性
    for (int pw=1; pw<=7; pw++){
        // 正常アリの数
        for (int n=0; n<=450; n+=50){
            double sensor = pow(10,-pw);
            int naho = (MACRO_NMAX - n);

            double sum = 0.0;
            for (int p=0; p<5; p++)
                prob[p] = 0;

            for(unsigned long long int dummy=1; dummy<=MACRO_MAX_STEP; dummy++){
                reset(sensor,naho,dummy);

                for (int id=0; id<MACRO_NUM_FOODS; id++)
                    foodspre[id] = MACRO_FOODSOURCE;

                for(int t=1; t<=MACRO_MAX_TIME; t++){
                    calculation();
                    cudaMemcpyFromSymbol(&foods, foods_d, sizeof(Food)*MACRO_NUM_FOODS);

                    for (int id=0; id<MACRO_NUM_FOODS; id++){
                        if( foods[id].vol < foodspre[id] + MACRO_REC ){
                            access[id] = 1;
                            prob[3+id] += 1;
                        }
                        else{
                            access[id] = 0;
                        }
                        foodspre[id] = foods[id].vol;
                    }
                    prob[access[0]+access[1]] += 1;
                    // if(t%500==0){
                    //     IOEffPoll(pw,n,dummy,t);
                    // }

                }

                sum += thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, homeOp, 0, binary_op);
            }
            // IOCellWrite(pw,n);
           IOEffWrite(pw,n,sum);
           IOProbWrite(pw,n,prob);
        }
        IOEffWrite(pw,500,normalSum);
        IOProbWrite(pw,500,probNormal);
    }

}

