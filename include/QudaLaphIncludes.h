#ifndef QUDA_LAPH_INCLUDES_H
#define QUDA_LAPH_INCLUDES_H

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <map>

#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

#if defined(USE_GSL_CBLAS)
#include <gsl_cblas.h>
#elif defined(USE_OPENBLAS)
#include <cblas.h>
#endif

#include <quda.h>
#include <quda_info.h>

#include "stop_watch.h"

#endif
