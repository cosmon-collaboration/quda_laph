/**
   Test the feasibility of Meson contracton on the device
 */
#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

#include <random>
#include <cassert>
#include <complex.h>

#include <quda.h>
#include <timer.h>
#include <blas_lapack.h>
#include <blas_quda.h>
#include <tune_quda.h>
#include <color_spinor_field.h>
#include <contract_quda.h>

using namespace quda ;
using namespace LaphEnv ;

// Cpu code

// does A_{ij} = B_{ik}C_{jk}
static inline void
innerevprod( const double _Complex *B ,
	     const double _Complex *C ,
	     const size_t nEv ,
	     double _Complex *A )
{
#if 1
  double _Complex alpha = 1.0 ;
  double _Complex beta  = 0.0 ;
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans ,
	      nEv, nEv, nEv , &alpha , B , nEv , C , nEv , &beta , A , nEv ) ;
#else
  // slow loopy version
  for( size_t i = 0 ; i < nEv ; i++ ) {
    for( size_t j = 0 ; j < nEv ; j++ ) {
      double _Complex sum = 0. ;
      for( size_t k = 0 ; k < nEv ; k++ ) {
	sum += B[ k + nEv*i ]*C[ j + nEv*k ] ;
      }
      A[ j + i*nEv ] = sum ;
    }
  }
#endif
}

// fwd guy is \tau^{\alpha\beta}_{jk}(t).\phi(p,tsrc)_{ij}
static void
compute_Mfwd( const double _Complex *phi ,
	      const double _Complex *tau ,
	      const size_t nEv ,
	      const size_t nmom ,
	      const size_t tsrc ,
	      const size_t LT ,
	      double _Complex *M )
{
  const double _Complex *B = tau ;
  for( size_t p = 0 ; p < nmom ; p++ ) {
    const double _Complex *phiP = phi + nEv*nEv*( tsrc + LT*p ) ;
    double _Complex *C = M + nEv*nEv*4*4*LT*p ;
    double _Complex alpha = 1.0 ;
    double _Complex beta  = 0.0 ;
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans ,
		LT*4*4*nEv, nEv, nEv , &alpha , B , nEv , phiP , nEv , &beta , C , nEv ) ;
  }
}

// backward prop is (\tau^{\alpha\beta}_{ik})(t) \phi(p,t)_{kj}
// note there is no gamma business happening here yet and no conjugation
static void
compute_Mbwd( const double _Complex *phi ,
	      const double _Complex *tau ,
	      const size_t nEv ,
	      const size_t nmom ,
	      const size_t tsrc ,
	      const size_t LT ,
	      double _Complex *M )
{
  for( size_t p = 0 ; p < nmom ; p++ ) {
    for( size_t t = 0 ; t < LT ; t++ ) {      
      const double _Complex *phiP = phi + nEv*nEv*( t + LT*p ) ;
      const double _Complex *B = tau + nEv*nEv*4*4*t ;
      double _Complex *C = M + nEv*nEv*4*4*(t+LT*p) ;
      double _Complex alpha = 1.0 ;
      double _Complex beta  = 0.0 ;
      cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans ,
		  4*4*nEv, nEv, nEv , &alpha , B , nEv , phiP , nEv , &beta , C , nEv ) ;
    }
  }
}

// trace of the product of ev indices Tr[ B.C^\dagger ]
static inline double _Complex
ev_traceprod_dag( const double _Complex *B ,
		  const double _Complex *C ,
		  const size_t nEv )
{
#if 1 
  return cblas_zdotc( nEv , C , 1 , B , 1 ) ;
#else
  double _Complex sum = 0.0 ;
  for( size_t i = 0 ; i < nEv ; i++ ) {
    sum += B[i] * conj( C[i] ) ;
  }
  return sum ;
#endif
}

// contract C^{\alpha\beta\kappa\delta}(psrc,psnk,t) = Mfwd.Mbwd
static void
contract_meson( const double _Complex *Mfwd ,
		const double _Complex *Mbwd ,
		const size_t nEv ,
		const size_t nmom ,
		const size_t LT ,
		double _Complex *C )
{
  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
      for( size_t t = 0 ; t < LT ; t++ ) {
	// open spin indices
	for( size_t osi = 0 ; osi < 4*4*4*4 ; osi++ ) {
	  const size_t alpha = (osi/64) ;
	  const size_t beta  = (osi/16)%4 ;
	  const size_t kappa = (osi/8)%4 ;
	  const size_t delta = (osi)%4 ;
	  const double _Complex *Mp1 = Mfwd + nEv*nEv*( alpha + 4*( beta + 4*( t + LT*psrc ) ) ) ;
	  const double _Complex *Mp2 = Mbwd + nEv*nEv*( delta + 4*( kappa + 4*( t + LT*psnk ) ) ) ;

	  *C = ev_traceprod_dag( Mp1 , Mp2 , nEv ) ; C++ ;
	}
      }
    }
  }
}

// Quda interface

//#define VERBOSE_COMPARISON
//#define GPU_STRESS
#define CPU_CROSSCHECK


int main(int argc, char *argv[]) {
  XMLHandler xml_in;

  if( init_quda_laph(argc, argv, xml_in) != 0 ) {
    exit(1) ;
  }

  int global = 1 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
#endif
  assert( global == 1 ) ;
  setVerbosityQuda(QUDA_VERBOSE, "#" , stdout ) ;

  // test for this many EVs
  const int nEv = 96 ;
  const int nmom = 64 ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  
  // create a "perambulator" - tau - of noise with these indices
  double _Complex *tau = (double _Complex*)malloc( X[3]*4*4*nEv*nEv*sizeof(double _Complex));
  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt ;
  for( size_t i = 0 ; i < X[3]*4*4*nEv*nEv ; i++ ) {
    const std::complex z( unif(mt) , unif(mt) ) ;
    tau[i] = z.real() + I*z.imag() ;
  }

  // create our "phi" matrix out of noise too
  double _Complex *phi = (double _Complex*)malloc( nmom*X[3]*nEv*nEv*sizeof(double _Complex));
  for( size_t i = 0 ; i < nmom*X[3]*nEv*nEv ; i++ ) {
    const std::complex z( unif(mt) , unif(mt) ) ;
    phi[i] = z.real() + I*z.imag() ;
  }

  // create a meson field or whatever bullshit people call this object for a specific psrc,psnk pair
  // M^{\alpha\beta}_{ij}(p,t)
  double _Complex *Mfwd = (double _Complex*)malloc( 4*4*nmom*X[3]*nEv*nEv*sizeof( double _Complex ) ) ;
  double _Complex *Mbwd = (double _Complex*)malloc( 4*4*nmom*X[3]*nEv*nEv*sizeof( double _Complex ) ) ;

  StopWatch fwd ; fwd.start() ;
  compute_Mfwd( phi , tau , nEv , nmom , 0 , X[3] , Mfwd ) ;
  fwd.stop() ;
  double CPUtime = fwd.getTimeInSeconds() ;
  printLaph( make_strf( "\n Mfwd in %g seconds\n" , CPUtime ) ) ;
  
  StopWatch bwd ; bwd.start() ;
  compute_Mbwd( phi , tau , nEv , nmom , 0 , X[3] , Mbwd ) ;
  bwd.stop() ;
  CPUtime = bwd.getTimeInSeconds() ;
  printLaph( make_strf( "\n Mbwd in %g seconds\n" , CPUtime ) ) ;

  // so now the contraction is -tr[ G1.mfwd.(gt.G2^\dagger.gt).G5.Mbwd.G5 ] over all indices but we can cheat by closing the ev indices first T^{alpha,beta,kappa,delta} = M^{alpha,beta} \bar{M^{kappa,delta} hmmmm
  free( tau ) ;
  free( phi ) ;

  StopWatch traces ; traces.start() ;
  // KISS look at the gamma_5 first Tr[ M M\bar ] all indices closed
  double _Complex *C = (double _Complex*)malloc( nmom*nmom*X[3]*256*sizeof(double _Complex) );

  contract_meson( Mfwd , Mbwd , nEv , nmom , X[3] , C ) ;

  traces.stop() ;

  CPUtime = traces.getTimeInSeconds() ;
  printLaph( make_strf( "\n CPU contraction time in %g seconds\n" , CPUtime ) ) ;
  free( C ) ;

  free( Mfwd ) ;
  free( Mbwd ) ;
  
  finalize( ) ;
  return 0;
}
