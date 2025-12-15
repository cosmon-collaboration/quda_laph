/**
   Test that color color contraction is sensible and then do a stress test
 **/
#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

#include <quda.h>
#include <omp.h>
#include <timer.h>
#include <blas_lapack.h>
#include <blas_quda.h>
#include <tune_quda.h>
#include <color_spinor_field.h>
#include <contract_quda.h>

#include <cassert>

#include "DWF.h"

using namespace LaphEnv ;
using namespace quda ;

int main(int argc, char *argv[]) {
#ifdef LAPH_DOMAIN_WALL
  
  XMLHandler xml_in;

  // just reuse tasks to do smearing + an inversion
  const int ntasks = init_quda_laph(argc, argv, xml_in) ;
  run_tasks( xml_in , ntasks ) ;

  int global = 1 ;
  #ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
  #endif
  assert( global == 1 ) ;

  const double m5 = 1.0 ;
  const double mass = 0.01 ;

  QudaInvertParam inv_param = newQudaInvertParam();
  inv_param.dslash_type = QUDA_WILSON_DSLASH;
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  inv_param.solve_type = QUDA_DIRECT_SOLVE;
  inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
  inv_param.cuda_prec = QUDA_DOUBLE_PRECISION;
  inv_param.dirac_order = QUDA_DIRAC_ORDER;
  inv_param.gamma_basis = QUDA_DIRAC_PAULI_GAMMA_BASIS ;
  inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.m5 = -1.0 ;
  inv_param.Ls = 8 ;
  inv_param.kappa = 1./(2.*(4.-m5)+1.) ;
  inv_param.mass = mass ;
  inv_param.tol = 1E-8 ;
  inv_param.reliable_delta = inv_param.reliable_delta_refinement = 0.1 ;
  inv_param.precondition_cycle = 1 ;
  inv_param.tol_precondition = 0.1 ;
  inv_param.maxiter_precondition = 1 ;
  inv_param.maxiter = 4000 ;
  inv_param.cuda_prec_sloppy = QUDA_SINGLE_PRECISION ;
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES ;
  inv_param.dirac_order = QUDA_DIRAC_ORDER ;
  inv_param.input_location = inv_param.output_location = QUDA_CPU_FIELD_LOCATION ;
  inv_param.dagger = QUDA_DAG_NO ;
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION ;
  inv_param.solution_type = QUDA_MAT_SOLUTION ;
  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN ;
  inv_param.verbosity = QUDA_SUMMARIZE ;
  inv_param.compute_true_res = QUDA_BOOLEAN_FALSE ;
  
  for( int k = 0 ; k < inv_param.Ls ; k++ ) {
    inv_param.b_5[k] = 1.5 ; inv_param.c_5[k] = 0.5 ;
  }
  inv_param.dslash_type = QUDA_MOBIUS_DWF_DSLASH ;
  inv_param.inv_type = QUDA_CG_INVERTER ;
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE ;
  printQudaInvertParam( &inv_param ) ;
  
  // create a point source at (0,0,0,0) in spin-color space
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int V =X[0]*X[1]*X[2]*X[3] ;

  // make a correlator
  std::vector<std::complex<double>> C( X[3] ) ;
  for( int t = 0 ; t < X[3] ; t++ ) {
    C[t] = 0.0 ;
  }
  
  // loop solution indices
  const int Nd = 4 , Nc = 3 ;
  for( int nspin = 0 ; nspin < Nd ; nspin++ ) {
    for( int nc = 0 ; nc < Nc ; nc++ ) {

      StopWatch Casio ;
      Casio.start() ;
      
      // latticefermion type objects
      LattField src( FieldSiteType::ColorSpinVector);
      LattField snk( FieldSiteType::ColorSpinVector);

      // point source at 0,0,0 should be the first element
      double *srcPtr = (double*)src.getDataPtr() ;
      double *snkPtr = (double*)snk.getDataPtr() ;

      memset( srcPtr , 0.0 , 2*12*V*sizeof(double) ) ;
      memset( snkPtr , 0.0 , 2*12*V*sizeof(double) ) ;

      // point source at origin? I think so...
      srcPtr[ 2*(nc + 3*nspin) ] = 1.0 ;

      // do solve
      solve_DWF<double>( (void*)snkPtr , (void*)srcPtr , &inv_param) ;

      Casio.stop() ;

      std::cout<<"Time for (" << nspin << "," << nc << ") -> " << Casio.getTimeInSeconds() << std::endl ;
      
      // contract < \psi | \psi > ~ pion
      for( int cb = 0 ; cb < 2 ; cb++ ) {
	const size_t nSp = 12*X[0]*X[1]*X[2] ;
	for( int t = 0 ; t < X[3] ; t++ ) {
	  double sumre = 0 , sumim = 0 ;
	  
	  for( size_t j = 0 ; j < nSp ; j+=2 ) {
	    sumre += snkPtr[j+0+nSp*t]*snkPtr[j+0+nSp*t] ;
	    sumre += snkPtr[j+1+nSp*t]*snkPtr[j+1+nSp*t] ;
	    
	    sumim += snkPtr[j+0+nSp*t]*snkPtr[j+1+nSp*t] ;
	    sumim -= snkPtr[j+1+nSp*t]*snkPtr[j+0+nSp*t] ;
	  }
	  C[t] += std::complex( sumre , sumim ) ;
	}
	snkPtr+=nSp*X[3] ;
      }
    }
  }

  // should be even-odd
  for( int t = 0 ; t < X[3] ; t++ ) {
    std::cout<< "CORR_" << t << " " << C[t].real() << " " << C[t].imag() << std::endl ;
  }
  
  finalize();

#else
  std::cout<<"quda_laph not set up for DWF"<<std::endl ;
#endif

  return 0;
}
