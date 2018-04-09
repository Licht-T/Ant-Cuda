#pragma once

#include "Constants.h"
#include <functional>

enum AntCharacter{
    NORMAL_CH,
    FOOL_CH
};

enum AntStatus{
    FORAGE,
    GOHOME,
    RANDOM_SEARCH,
    EMERGENCY
};

enum CELLStatus{
    NORMAL_CELL         = 0,
    FOOD_CELL           = 1,
    FOOD_NEIGHBOUR_CELL = 1 << 1,
    NEST_CELL           = 1 << 2,
    NEST_NEIGHBOUR_CELL = 1 << 3,

    NEAR_FOOD           = FOOD_CELL | FOOD_NEIGHBOUR_CELL,
    NEAR_NEST           = NEST_CELL | NEST_NEIGHBOUR_CELL
};

enum Direction{
    NONE      = 0,
    UP        = 1 ,
    UPRIGHT   = 1 << 1,
    LOWRIGHT  = 1 << 2,
    LOW       = 1 << 3,
    LOWLEFT   = 1 << 4,
    UPLEFT    = 1 << 5,
    UP_AND_UR = UP | UPRIGHT,
    RIGHT     = UPRIGHT | LOWRIGHT,
    LR_AND_LO = LOWRIGHT | LOW,
    LO_AND_LL = LOW | LOWLEFT,
    LEFT      = LOWLEFT | UPLEFT,
    UL_AND_UP = UPLEFT | UP,

    ALL_DIR   = UP | UPRIGHT | LOWRIGHT | LOW | LOWLEFT | UPLEFT,
    EDGE      = ALL_DIR
};

enum BitCalcAllowedTypes {
    CELLStatus,
    Direction
};

typedef struct{
    double x;
    double y;
} Cartesian;

typedef struct {
    enum AntStatus status;
    int i,j;
    enum Direction dir;
    int searchTime;
    enum AntCharacter ch;
    int homing[MACRO_NUM_FOODS];
    int _foodNo;
} Ant;


struct getHoming{
    __host__ __device__ int operator()(const Ant& x) const { 
        int sum = 0;
        for (int i=0; i<MACRO_NUM_FOODS; i++){
            sum += x.homing[i];
        }
        return sum;
    }
};

template <AntCharacter t> struct getHomingWithType{
    __host__ __device__ int operator()(const Ant& x) const {
        if(x.ch==t){
            int sum = 0;
            for (int i=0; i<MACRO_NUM_FOODS; i++){
                sum += x.homing[i];
            }
            return sum;
        }
        else{
            return 0;
        }
    }
};

template <int t> struct getHomingWithFoodNum{
    __host__ __device__ int operator()(const Ant& x) const {
        return x.homing[t];
    }
};

template <AntCharacter t,int i> struct getHomingWithTypeAndFoodNum{
    __host__ __device__ int operator()(const Ant& x) const {
        if(x.ch==t){
            return x.homing[i];
        }
        else{
            return 0;
        }
    }
};

typedef struct {
    double vol;
    int i,j;
} Food;

typedef struct {
    enum CELLStatus status;
    enum Direction edge;
    enum Direction nestDir;
    enum Direction criticalAngle;
    double phero;
    int foodNo;
    int i,j;
    Cartesian cart;
    double distFromNest;
    enum Direction nearestDirFromNestList[6];
} Cell;

struct cell_phero_mul : public std::binary_function<Cell, Cell, double> {
    double operator()(Cell a, Cell b) { return a.phero * b.phero; }
};

struct cell_phero_plus : public std::binary_function<Cell, Cell, double> {
    double operator()(Cell a, Cell b) { return a.phero + b.phero; }
};
