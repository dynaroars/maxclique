#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif



#include "MT_rand.h"
#include "Dimacs.h"
#include "config.h"
#include "tvnFuncs.h"
#include "AntsOps.h" //header for Ants operations s.a moving etc
#include "localOptimize.h"
#include "localOptimizePAR.h"
//#include "Initialize.h"


