#pragma once

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "Constants.h"
#include "DataStructures.h"


extern Food foods[NUM_FOODS];
extern __device__ Food foods_d[NUM_FOODS];
extern Food *foods_d_ptr_raw;
extern thrust::device_ptr<Food> foods_d_ptr;

extern Ant ants[NMAX];
extern __device__ Ant ants_d[NMAX];
extern Ant *ants_d_ptr_raw;
extern thrust::device_ptr<Ant> ants_d_ptr;

extern Cell cells[MAX][MAX];
extern __device__ Cell cells_d[MAX][MAX];
extern Cell *cells_d_ptr_raw;
extern thrust::device_ptr<Cell> cells_d_ptr;

extern __device__ unsigned int sort_key_d[NMAX];
extern unsigned int *sort_key_d_ptr_raw;
extern thrust::device_ptr<unsigned int> sort_key_d_ptr;

extern __device__ unsigned long long int seeds_d[NMAX];
extern unsigned long long int *seeds_d_ptr_raw;
extern thrust::device_ptr<unsigned long long int> seeds_d_ptr;

extern int homing;
extern __device__ int homing_d;


