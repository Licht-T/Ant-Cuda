#include "Variables2.h"

Food foods[MACRO_NUM_FOODS];
__device__ Food foods_d[MACRO_NUM_FOODS];
Food *foods_d_ptr_raw;
thrust::device_ptr<Food> foods_d_ptr;

Ant ants[MACRO_NMAX];
__device__ Ant ants_d[MACRO_NMAX];
Ant *ants_d_ptr_raw;
thrust::device_ptr<Ant> ants_d_ptr;

Cell cells[MACRO_MAX][MACRO_MAX];
__device__ Cell cells_d[MACRO_MAX][MACRO_MAX];
Cell *cells_d_ptr_raw;
thrust::device_ptr<Cell> cells_d_ptr;

Cell gauss_path1[MACRO_MAX][MACRO_MAX];
__device__ Cell gauss_path1_d[MACRO_MAX][MACRO_MAX];
Cell *gauss_path1_d_ptr_raw;
thrust::device_ptr<Cell> gauss_path1_d_ptr;

Cell gauss_path2[MACRO_MAX][MACRO_MAX];
__device__ Cell gauss_path2_d[MACRO_MAX][MACRO_MAX];
Cell *gauss_path2_d_ptr_raw;
thrust::device_ptr<Cell> gauss_path2_d_ptr;

__device__ unsigned int sort_key_d[MACRO_NMAX];
unsigned int *sort_key_d_ptr_raw;
thrust::device_ptr<unsigned int> sort_key_d_ptr;

__device__ unsigned long long int seeds_d[MACRO_NMAX];
unsigned long long int *seeds_d_ptr_raw;
thrust::device_ptr<unsigned long long int> seeds_d_ptr;

int homing = 0;
__device__ int homing_d = 0;
double x_Nest2Food[2] = {0.0, 0.0};
double y_Nest2Food[2] = {0.0, 0.0};
__device__ double x_Nest2Food_d[2];
__device__ double y_Nest2Food_d[2];

