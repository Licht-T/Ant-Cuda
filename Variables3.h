#pragma once

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "Constants.h"
#include "DataStructures2.h"


extern Food foods[MACRO_NUM_FOODS];
extern __device__ Food foods_d[MACRO_NUM_FOODS];
extern Food *foods_d_ptr_raw;
extern thrust::device_ptr<Food> foods_d_ptr;

extern Ant ants[MACRO_NMAX];
extern __device__ Ant ants_d[MACRO_NMAX];
extern Ant *ants_d_ptr_raw;
extern thrust::device_ptr<Ant> ants_d_ptr;

extern Cell cells[MACRO_MAX][MACRO_MAX];
extern __device__ Cell cells_d[MACRO_MAX][MACRO_MAX];
extern Cell *cells_d_ptr_raw;
extern thrust::device_ptr<Cell> cells_d_ptr;

extern Cell min_path1[MACRO_MAX][MACRO_MAX];
extern __device__ Cell min_path1_d[MACRO_MAX][MACRO_MAX];
extern Cell *min_path1_d_ptr_raw;
extern thrust::device_ptr<Cell> min_path1_d_ptr;

extern Cell min_path2[MACRO_MAX][MACRO_MAX];
extern __device__ Cell min_path2_d[MACRO_MAX][MACRO_MAX];
extern Cell *min_path2_d_ptr_raw;
extern thrust::device_ptr<Cell> min_path2_d_ptr;

extern __device__ unsigned int sort_key_d[MACRO_NMAX];
extern unsigned int *sort_key_d_ptr_raw;
extern thrust::device_ptr<unsigned int> sort_key_d_ptr;

extern __device__ unsigned long long int seeds_d[MACRO_NMAX];
extern unsigned long long int *seeds_d_ptr_raw;
extern thrust::device_ptr<unsigned long long int> seeds_d_ptr;

extern int homing;
extern __device__ int homing_d;
extern double x_Nest2Food[2];
extern double y_Nest2Food[2];
extern __device__ double x_Nest2Food_d[2];
extern __device__ double y_Nest2Food_d[2];
