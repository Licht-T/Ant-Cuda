#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

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
#include "kernel2.h"
//#include "kernel.h"
#include "Display.h"
#include "IO2.h"

int main(int argc, char *argv[]){


    getHoming homeOp;
    thrust::plus<int> binary_op;

    IOInit();
    std::string step = toString(MACRO_MAX_STEP);
    std::string angle = toString(MACRO_FOOD_ANGLE);

    initialize();

    double normalEndTime[MACRO_NUM_FOODS]; // 有限量の餌を採り尽くす時間
    double normalFindTime[MACRO_NUM_FOODS]; // 有限量の餌を採り尽くす時間
    for (int id=0; id<MACRO_NUM_FOODS; id++){
        normalEndTime[id] = 0;
        normalFindTime[id] = 0;
    }
    double normalSum = 0.0;
    // std::vector< std::vector<double> > FoodsAmount(MACRO_MAX_TIME, std::vector<double>(MACRO_NUM_FOODS, 0));

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
        // std::string pwstr = toString(0);
        // std::string nstr = toString(0);
        // std::string estr = toString(dummy);
        // std::string fool_num_foodamount = std::string(step+"steps_"+angle+"deg_"+pwstr+"_"+nstr+"_"+estr+"_"+"_foodamount.dat");
        // std::ofstream ofsfoodamount(fool_num_foodamount.c_str());
        // FoodsAmount = std::vector< std::vector<double> >(MACRO_MAX_TIME, std::vector<double>(MACRO_NUM_FOODS, 0));
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
        int FinFoodNum = 0;
        reset(1,0,dummy);

        // display(argc,argv);

        for(int t=0; t<MACRO_MAX_TIME; t++){
            calculation();
            // if(t%500==0){
            //     IOEffPoll(0,500,dummy,t);
            // }

            cudaMemcpyFromSymbol(&foods, foods_d, sizeof(Food)*MACRO_NUM_FOODS);
            for (int id=0; id<MACRO_NUM_FOODS; id++){
                // FoodsAmount[t][id] = foods[id].vol;
                if (foods[id].vol < MACRO_FOODSOURCE && FindFlag[id]==0)
                {
                    FindTime[id] = t;
                    FindFlag[id] = 1;
                }
                if (foods[id].vol < 0.04 && EndFlag[id] == 0){
                    EndTime[id] = t;
                    EndFlag[id] = 1;
                    FinFoodNum++;
                }
            }
            if (FinFoodNum == MACRO_NUM_FOODS)
                break;
        }
        for (int id=0; id<MACRO_NUM_FOODS; id++){
            normalEndTime[id] += EndTime[id];
            normalFindTime[id] += FindTime[id];
        }

        // IOFoodAmountWrite(ofsfoodamount, FoodsAmount);
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
                // std::string pwstr = toString(pw);
                // std::string nstr = toString(n);
                // std::string estr = toString(dummy);
                // std::string fool_num_foodamount = std::string(step+"steps_"+angle+"deg_"+pwstr+"_"+nstr+"_"+estr+"_"+"_foodamount.dat");
                // std::ofstream ofsfoodamount(fool_num_foodamount.c_str());
                // FoodsAmount = std::vector< std::vector<double> >(MACRO_MAX_TIME, std::vector<double>(MACRO_NUM_FOODS, 0));
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
                int FinFoodNum = 0;
                reset(sensor,naho,dummy);

                for(int t=1; t<=MACRO_MAX_TIME; t++){
                    calculation();
                    // if(t%500==0){
                    //     IOEffPoll(pw,n,dummy,t);
                    // }

                    cudaMemcpyFromSymbol(&foods, foods_d, sizeof(Food)*MACRO_NUM_FOODS);
                    for (int id=0; id<MACRO_NUM_FOODS; id++){
                        // FoodsAmount[t][id] = foods[id].vol;
                        if (foods[id].vol < MACRO_FOODSOURCE && FindFlag[id]==0)
                        {
                            FindTime[id] = t;
                            FindFlag[id] = 1;
                        }
                        if (foods[id].vol < 0.04 && EndFlag[id] == 0){
                            EndTime[id] = t;
                            EndFlag[id] = 1;
                            FinFoodNum++;
                        }
                    }
                    if (FinFoodNum == MACRO_NUM_FOODS)
                        break;
                }

                for (int id=0; id<MACRO_NUM_FOODS; id++){
                    EndTimeAve[id] += EndTime[id];
                    FindTimeAve[id] += FindTime[id];
                }
                // IOFoodAmountWrite(ofsfoodamount, FoodsAmount);
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

