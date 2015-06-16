#include "Variables.h"

Food foods[NUM_FOODS];
__device__ Food foods_d[NUM_FOODS];
Food *foods_d_ptr_raw;
thrust::device_ptr<Food> foods_d_ptr;

Ant ants[NMAX];
__device__ Ant ants_d[NMAX];
Ant *ants_d_ptr_raw;
thrust::device_ptr<Ant> ants_d_ptr;

Cell cells[MAX][MAX];
__device__ Cell cells_d[MAX][MAX];
Cell *cells_d_ptr_raw;
thrust::device_ptr<Cell> cells_d_ptr;

__device__ unsigned int sort_key_d[NMAX];
unsigned int *sort_key_d_ptr_raw;
thrust::device_ptr<unsigned int> sort_key_d_ptr;

__device__ unsigned long long int seeds_d[NMAX];
unsigned long long int *seeds_d_ptr_raw;
thrust::device_ptr<unsigned long long int> seeds_d_ptr;

int homing = 0;
__device__ int homing_d = 0;

