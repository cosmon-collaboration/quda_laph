#ifndef INIT_QUDA_LAPH
#define INIT_QUDA_LAPH

#include "xml_handler.h"

// returns "ntasks"
int init_quda_laph( int argc , char *argv[] , XMLHandler &xml_in ) ;

void run_tasks( XMLHandler &xml_info , const int ntasks ) ;

void finalize() ;

#endif
