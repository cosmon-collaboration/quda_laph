#ifndef QUDA_LAPH_BLAS_H
#define QUDA_LAPH_BLAS_H

#if defined(USE_GSL_CBLAS)
#include "gsl_cblas.h"
#elif defined(USE_OPENBLAS)
#include "cblas.h"
#endif

#endif
