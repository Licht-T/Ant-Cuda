#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <cuda_runtime.h>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef _WIN32
#include <direct.h>
#endif

#include "Constants.h"
#include "DataStructures.h"
#include "Variables.h"
#include "kernel12.h"

/* const std::string constHeader("Constants.h"); */
/* std::string path; */
/* std::string fool_num_plot; */
/* std::ofstream *ofs; */
/* std::string fool_num_plot_prob; */
/* std::ofstream *ofsProb; */
extern std::string fool_num_food;
/* std::ofstream *ofsfood; */


template <class T> std::string toString(T t){
    static std::ostringstream ss;
    ss << (T)t;

    std::string str(ss.str());
    ss.str("");
    ss.clear(std::stringstream::goodbit);

    return str;
}


void IOInit();

void IOEffWrite(int pw, int n, double sum);
void IOCellWrite(int pw, int n);
void IOEffPoll(int pw, int n, int step, int t);
void IOFinWrite(int pw, int n, double sum);
void IOFoodWrite(int pw, int n, double ft[], double et[]);
void IOFoodAmountWrite(std::ofstream & ofs, double ft[][MACRO_NUM_FOODS]);// std::vector< std::vector<double> > & ft);
void IOProbWrite(int pw, int n, double prob[]);
