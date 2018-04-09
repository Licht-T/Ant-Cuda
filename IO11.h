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
#include "DataStructures2.h"
#include "Variables3.h"
#include "kernel11.h"



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
//void IOCellWrite(int pw, int n);
void IOCellWrite(int pw, int n, int sample, int t);
void IOEffPoll(int pw, int n, int step, int t);
void IOFinWrite(int pw, int n, double sum);
void IOFoodWrite(int pw, int n, double ft[], double et[]);
void IOFoodAmountWrite(std::ofstream & ofs, std::vector< std::vector<double> > & ft);
void IOProbWrite(int pw, int n, double prob[]);
void IOPheroStateWrite(int pw, int n, int sample, double conv [][3]);
void IOPheroPathWrite(int pw, int n, int sample, int t, std::vector< std::vector< std::vector<double> > > & pathxyp);
void IOPheroPathLengthWrite(int pw, int n, std::vector< std::vector< std::vector<double> > > & pathlen);