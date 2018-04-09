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
#include "kernel4.h"
#include "Display.h"
#include "IO4.h"

int main(int argc, char *argv[]){

    getHoming homeOp;
    thrust::plus<int> binary_op;

    IOInit();

    initialize();

    double normalEndTime[MACRO_NUM_FOODS]; // 有限量の餌を採り尽くす時間
    double normalFindTime[MACRO_NUM_FOODS]; // 有限量の餌を採り尽くす時間
    for (int id=0; id<MACRO_NUM_FOODS; id++){
        normalEndTime[id] = 0;
        normalFindTime[id] = 0;
    }
    double normalSum = 0.0;

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
        int FindTime[MACRO_NUM_FOODS]; // 有限量の餌を見つける時間
        int FindFlag[MACRO_NUM_FOODS];
        int EndTime[MACRO_NUM_FOODS]; // 有限量の餌を採り尽くす時間
        int EndFlag[MACRO_NUM_FOODS];
        for (int id=0; id<MACRO_NUM_FOODS; id++){
            EndTime[id] = 0;
            EndFlag[id] = 0;
            FindTime[id] = 0;
            FindFlag[id] = 0;
        }
        reset(1,0,dummy);
        foods[0].vol = MACRO_FOODSOURCE;
        foods[1].vol = 0;
        cudaMemcpyFromSymbol(foods_d, foods, sizeof(Food)*MACRO_NUM_FOODS);

        // display(argc,argv);

        for(int t=0; t<MACRO_MAX_TIME; t++){
            calculation();
            cudaMemcpyFromSymbol(&foods, foods_d, sizeof(Food)*MACRO_NUM_FOODS);
            if (foods[0].vol < MACRO_FOODSOURCE && FindFlag[0] == 0){
                FindTime[0] = t;
                FindFlag[0] = 1;
            }
            if (foods[1].vol < MACRO_FOODSOURCE && FindFlag[1] == 0 && EndFlag[0] == 1){
                FindTime[1] = t;
                FindFlag[1] = 1;
            }
            if (foods[0].vol < 0.04 && EndFlag[0]==0){
                foods[0].vol = 0;
                foods[1].vol = MACRO_FOODSOURCE;
                cudaMemcpyFromSymbol(foods_d, foods, sizeof(Food)*MACRO_NUM_FOODS);
                EndTime[0] = t;
                EndFlag[0] = 1;
            }
            if (EndFlag[0]==1 && foods[1].vol < 0.04){
                EndTime[1] = t;
                EndFlag[1] = 1;
                break;
            }
            // if(t%500==0){
            //     IOEffPoll(0,500,dummy,t);
            // }

        }
        for (int id=0; id<MACRO_NUM_FOODS; id++){
            normalEndTime[id] += EndTime[id];
            normalFindTime[id] += FindTime[id];
        }

        normalSum += thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, homeOp, 0, binary_op);
    }

    for (int n=0; n<=500; n+=50){
        IOEffWrite(0,n,normalSum);
        IOFoodWrite(0,n, normalFindTime, normalEndTime);
    }

    // 馬鹿アリの感受性
    for (int pw=1; pw<=7; pw++){
        // 正常アリの数
        for (int n=0; n<=450; n+=50){
            double sensor = pow(10,-pw);
            int naho = (MACRO_NMAX - n);

            double sum = 0.0;
            double EndTimeAve[MACRO_NUM_FOODS];
            double FindTimeAve[MACRO_NUM_FOODS];
            for (int id=0; id<MACRO_NUM_FOODS; id++){
                EndTimeAve[id] = 0;
                FindTimeAve[id] = 0;
            }
            for(unsigned long long int dummy=1; dummy<=MACRO_MAX_STEP; dummy++){
                int FindTime[MACRO_NUM_FOODS]; // 有限量の餌を見つける時間
                int FindFlag[MACRO_NUM_FOODS];
                int EndTime[MACRO_NUM_FOODS]; // 有限量の餌を採り尽くす時間
                int EndFlag[MACRO_NUM_FOODS];
                for (int id=0; id<MACRO_NUM_FOODS; id++){
                    EndTime[id] = 0;
                    EndFlag[id] = 0;
                    FindTime[id] = 0;
                    FindFlag[id] = 0;
                }
                reset(sensor,naho,dummy);
                foods[0].vol = MACRO_FOODSOURCE;
                foods[1].vol = 0;
                cudaMemcpyFromSymbol(foods_d, foods, sizeof(Food)*MACRO_NUM_FOODS);

                for(int t=1; t<=MACRO_MAX_TIME; t++){
                    calculation();
                    cudaMemcpyFromSymbol(&foods, foods_d, sizeof(Food)*MACRO_NUM_FOODS);
                    if (foods[0].vol < MACRO_FOODSOURCE && FindFlag[0] == 0){
                        FindTime[0] = t;
                        FindFlag[0] = 1;
                    }
                    if (foods[1].vol < MACRO_FOODSOURCE && FindFlag[1] == 0 && EndFlag[0] == 1){
                        FindTime[1] = t;
                        FindFlag[1] = 1;
                    }
                    if (foods[0].vol < 0.04 && EndFlag[0]==0){
                        foods[0].vol = 0;
                        foods[1].vol = MACRO_FOODSOURCE;
                        cudaMemcpyFromSymbol(foods_d, foods, sizeof(Food)*MACRO_NUM_FOODS);
                        EndTime[0] = t;
                        EndFlag[0] = 1;
                    }
                    if (EndFlag[0]==1 && foods[1].vol < 0.04){
                        EndTime[1] = t;
                        EndFlag[1] = 1;
                        break;
                    }
                    // if(t%500==0){
                    //     IOEffPoll(pw,n,dummy,t);
                    // }

                }

                for (int id=0; id<MACRO_NUM_FOODS; id++){
                    EndTimeAve[id] += EndTime[id];
                    FindTimeAve[id] += FindTime[id];
                }
                sum += thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, homeOp, 0, binary_op);
            }
            // IOCellWrite(pw,n);
            IOEffWrite(pw,n,sum);
            IOFoodWrite(pw,n,FindTimeAve,EndTimeAve);
        }
        IOEffWrite(pw,500,normalSum);
        IOFoodWrite(pw,500, normalFindTime, normalEndTime);
    }

}

