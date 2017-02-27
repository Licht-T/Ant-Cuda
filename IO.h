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
#include "kernel.h"

void IOInit();

void IOEffWrite(int pw, int n, double sum);
void IOCellWrite(int pw, int n);
void IOEffPoll(int pw, int n, int step, int t);
