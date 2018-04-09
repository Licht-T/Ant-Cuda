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
#include <thrust/inner_product.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

#include "Constants.h"
#include "DataStructures2.h"
#include "Variables3.h"
#include "kernel11.h"
#include "Display11.h"
#include "IO11.h"

int main(int argc, char *argv[]){

    getHoming homeOp;
    thrust::plus<int> binary_op;
    thrust::plus<double> binary_op_add;

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

    double gp1norm = sqrt(thrust::inner_product(min_path1_d_ptr, min_path1_d_ptr + MACRO_MAX*MACRO_MAX, min_path1_d_ptr, 0.0, binary_op_add, cell_phero_mul()));
    double gp2norm = sqrt(thrust::inner_product(min_path2_d_ptr, min_path2_d_ptr + MACRO_MAX*MACRO_MAX, min_path2_d_ptr, 0.0, binary_op_add, cell_phero_mul()));

    cudaMemcpyFromSymbol(&cells, cells_d, sizeof(Cell)*MACRO_MAX*MACRO_MAX);
    // 個々のセルがどの最短経路上のセルに対応しているかの情報を持つ
    cudaMemcpyFromSymbol(&min_path1, min_path1_d, sizeof(Cell)*MACRO_MAX*MACRO_MAX);
    cudaMemcpyFromSymbol(&min_path2, min_path2_d, sizeof(Cell)*MACRO_MAX*MACRO_MAX);

    size_t mp_len[MACRO_NUM_FOODS] = {0, 0};
    std::vector< std::vector< std::vector<int> > > gp;
    gp.resize(MACRO_NUM_FOODS);

    Cell (*p_mp[2])[MACRO_MAX];
    p_mp[0] = min_path1;
    p_mp[1] = min_path2;

    for (int fi=0; fi<MACRO_NUM_FOODS; fi++){
        for (int ci=0; ci<MACRO_MAX; ci++){
            for (int cj=0; cj<MACRO_MAX; cj++){
                Cell * c1 = &p_mp[fi][ci][cj];
                // Cell * c1 = &p_mp[ci][cj];
                std::cout << fi << " : " << ci*MACRO_MAX+(cj+1) << " / " << MACRO_MAX*MACRO_MAX << std::endl;
                if (c1->i == ci & c1->j == cj){
                    std::cout << "on" << std::endl;
                    std::vector<int> v(ci, cj);
                    gp[fi].push_back(v);
                }
            }
        }
        mp_len[fi] = gp[fi].size();
    }

    std::vector< std::vector< std::vector<double> > > pathlen(MACRO_NUM_FOODS, std::vector< std::vector<double> >(int(MACRO_MAX_TIME/1000), std::vector<double>(MACRO_MAX_STEP, 0.0)));

    // // 馬鹿ありのいない場合
    // for(unsigned long long int dummy=1; dummy<=MACRO_MAX_STEP; dummy++){
    //     reset(1,0,dummy);

    //     // display(argc,argv);

    //     int ti = 0;
    //     for (int id=0; id<MACRO_NUM_FOODS; id++)
    //         foodspre[id] = MACRO_FOODSOURCE;
    //     for(int t=0; t<MACRO_MAX_TIME; t++){
    //         calculation();
    //         cudaMemcpyFromSymbol(&foods, foods_d, sizeof(Food)*MACRO_NUM_FOODS);
    //         for (int id=0; id<MACRO_NUM_FOODS; id++){
    //             if( foods[id].vol < foodspre[id] + MACRO_REC ){
    //                 access[id] = 1;
    //                 probNormal[3+id] += 1;
    //             }
    //             else{
    //                 access[id] = 0;
    //             }
    //             foodspre[id] = foods[id].vol;
    //         }
    //         probNormal[access[0]+access[1]] += 1;
    //         // if(t%500==0){
    //         //     IOEffPoll(0,500,dummy,t);
    //         // }
    //         if(t%1000==0){
    //             cudaMemcpyFromSymbol(&cells, cells_d, sizeof(Cell)*MACRO_MAX*MACRO_MAX);
    //             std::vector< std::vector< std::vector<double> > > mpxyp;//(2, std::vector< std::vector< int> > (mp_len[0], std::vector<int> (3, 0)) );
    //             mpxyp.resize(MACRO_NUM_FOODS);
    //             for(int fi=0; fi<MACRO_NUM_FOODS; fi++){
    //                 mpxyp[fi].resize(mp_len[fi]);
    //                 for(int l=0; l<mp_len[fi]; l++){
    //                     mpxyp[fi][l].resize(4, 0.0);
    //                 }
    //             }
    //             // パスを得るアルゴリズム!!
    //             for (int fi=0; fi<MACRO_NUM_FOODS; fi++){
    //                 // パスを得る
    //                 for (int l_i=0; l_i<mp_len[fi]; l_i++){
    //                     for (int ci=0; ci<MACRO_MAX; ci++){
    //                         for (int cj=0; cj<MACRO_MAX; cj++){
    //                             Cell * c = &p_mp[fi][ci][cj];
    //                             if (c->i == gp[fi][l_i][0] & c->j == gp[fi][l_i][1])
    //                                 if ( cells[c->i][c->j].phero > mpxyp[fi][l_i][2] ){
    //                                     mpxyp[fi][l_i][0] = c->cart.x;
    //                                     mpxyp[fi][l_i][1] = c->cart.x;
    //                                     mpxyp[fi][l_i][2] = cells[c->i][c->j].phero;
    //                                 }
    //                         }
    //                     }
    //                 }
    //                 // 長さを測る
    //                 for (int l_i=1; l_i<mp_len[fi]; l_i++){
    //                     mpxyp[fi][l_i][3] += sqrt( (mpxyp[fi][l_i-1][0] - mpxyp[fi][l_i][0])*(mpxyp[fi][l_i-1][0] - mpxyp[fi][l_i][0]) + (mpxyp[fi][l_i-1][1] - mpxyp[fi][l_i][1])*(mpxyp[fi][l_i-1][1] - mpxyp[fi][l_i][1]) );
    //                 }
    //                 pathlen[fi][t%1000][dummy-1] = mpxyp[fi][mp_len[fi]-1][3];
    //             }

    //             IOCellWrite(0, 500, dummy, t);
    //             IOPheroPathWrite(0, 500, dummy, t, mpxyp);
    //         }

    //     }

    //     normalSum += thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, homeOp, 0, binary_op);
    // }
    // IOPheroPathLengthWrite(0, 500, pathlen);
    // for (int n=0; n<=500; n+=50){
    //     IOEffWrite(0,n,normalSum);
    //     // IOProbWrite(0,n,probNormal);
    // }

    // 馬鹿アリの感受性
    // for (int pw=1; pw<=7; pw++)
    {
        int pw = 2;
        // 正常アリの数
        // for (int n=0; n<=450; n+=50)
        {
            int n = 200;
            double sensor = pow(10,-pw);
            int naho = (MACRO_NMAX - n);

            double sum = 0.0;
            for (int p=0; p<5; p++)
                prob[p] = 0;

            std::vector< std::vector< std::vector<double> > > pathlen(MACRO_NUM_FOODS, std::vector< std::vector<double> >(int(MACRO_MAX_TIME/1000), std::vector<double>(MACRO_MAX_STEP, 0.0)));

            for(unsigned long long int dummy=1; dummy<=10; dummy++){
            // for(unsigned long long int dummy=1; dummy<=MACRO_MAX_STEP; dummy++){
                reset(sensor,naho,dummy);

                // std::string pwstr = toString(pw);
                // std::string nstr = toString(n);
                // std::string anglestr = toString(MACRO_FOOD_ANGLE);
                // std::string samplestr = toString(dummy);

                // std::string pherostatedata(path+"cell_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_sampleNo"+samplestr+"_phero_state.dat");
                // std::ofstream pherostate_fs(pherostatedata.c_str());

                for (int id=0; id<MACRO_NUM_FOODS; id++)
                    foodspre[id] = MACRO_FOODSOURCE;

                int ti = 0;

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
                    if(t%1000==0){
                        cudaMemcpyFromSymbol(&cells, cells_d, sizeof(Cell)*MACRO_MAX*MACRO_MAX);
                        std::vector< std::vector< std::vector<double> > > mpxyp;//(2, std::vector< std::vector< int> > (mp_len[0], std::vector<int> (3, 0)) );
                        mpxyp.resize(MACRO_NUM_FOODS);
                        for(int fi=0; fi<MACRO_NUM_FOODS; fi++){
                            mpxyp[fi].resize(mp_len[fi]);
                            for(int l=0; l<mp_len[fi]; l++){
                                mpxyp[fi][l].resize(4, 0.0);
                            }
                        }
                        // パスを得るアルゴリズム!!
                        for (int fi=0; fi<MACRO_NUM_FOODS; fi++){
                            // パスを得る
                            for (int l_i=0; l_i<mp_len[fi]; l_i++){
                                for (int ci=0; ci<MACRO_MAX; ci++){
                                    for (int cj=0; cj<MACRO_MAX; cj++){
                                        Cell * c = &p_mp[fi][ci][cj];
                                        if (c->i == gp[fi][l_i][0] & c->j == gp[fi][l_i][1])
                                            if ( cells[c->i][c->j].phero > mpxyp[fi][l_i][2] ){
                                                mpxyp[fi][l_i][0] = c->cart.x;
                                                mpxyp[fi][l_i][1] = c->cart.x;
                                                mpxyp[fi][l_i][2] = cells[c->i][c->j].phero*10;
                                            }
                                    }
                                }
                            }
                            // 長さを測る
                            for (int l_i=1; l_i<mp_len[fi]; l_i++){
                                mpxyp[fi][l_i][3] += sqrt( (mpxyp[fi][l_i-1][0] - mpxyp[fi][l_i][0])*(mpxyp[fi][l_i-1][0] - mpxyp[fi][l_i][0]) + (mpxyp[fi][l_i-1][1] - mpxyp[fi][l_i][1])*(mpxyp[fi][l_i-1][1] - mpxyp[fi][l_i][1]) );
                            }
                            pathlen[fi][t%1000][dummy-1] = mpxyp[fi][mp_len[fi]-1][3];
                        }

                        IOCellWrite(pw, n, dummy, t);
                        IOPheroPathWrite(pw, n, dummy, t, mpxyp);
                    }

                }

                IOPheroPathLengthWrite(pw, n, pathlen);
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
