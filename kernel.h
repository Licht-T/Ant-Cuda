#include <cmath>
#include <algorithm>
#include <cstdlib>

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/sort.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/copy.h>

#include "Constants.h"
#include "DataStructures.h"
#include "Variables.h"


__host__ void initialize();
__host__ void reset(double sensor, int naho, unsigned long long int step);
__host__ void calculation();

__global__ void sortKeyInit();

__device__ __host__ enum Direction operator<<(enum Direction d, int i);
__device__ __host__ enum Direction operator>>(enum Direction d, int i);
__device__ __host__ enum Direction operator|(enum Direction d1, enum Direction d2);
__device__ __host__ enum Direction operator&(enum Direction d1, enum Direction d2);
__device__ __host__ enum Direction& operator|=(enum Direction& d1, enum Direction d2);
__device__ __host__ enum Direction& operator&=(enum Direction& d1, enum Direction d2);
__device__ __host__ enum Direction& operator<<=(enum Direction& d1, int i);
__device__ __host__ enum Direction& operator>>=(enum Direction& d1, int i);
__device__ __host__ bool operator<=(enum Direction d1, enum Direction d2);
__device__ __host__ Cell* getCell(Cell cells[MAX][MAX],int i,int j, enum Direction dir);



__device__ __host__  enum CELLStatus operator<<(enum CELLStatus d, int i);
__device__ __host__  enum CELLStatus operator>>(enum CELLStatus d, int i);
__device__ __host__  enum CELLStatus operator|(enum CELLStatus d1, enum CELLStatus d2);
__device__ __host__  enum CELLStatus operator&(enum CELLStatus d1, enum CELLStatus d2);
__device__ __host__  enum CELLStatus& operator|=(enum CELLStatus& d1, enum CELLStatus d2);
__device__ __host__  enum CELLStatus& operator&=(enum CELLStatus& d1, enum CELLStatus d2);

__device__ __host__ enum Direction left(enum Direction dir);
__device__ __host__ enum Direction right(enum Direction dir);
__device__ __host__ Cell* getCell(Cell cells[MAX][MAX],int i,int j, enum Direction dir);

