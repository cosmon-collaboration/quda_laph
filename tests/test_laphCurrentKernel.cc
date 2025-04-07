#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

#include <random>
#include <cassert>
#include <complex.h>

using namespace LaphEnv ;

//#define VERBOSE

static void cpu_code( const int n1,
		      const int n2,
		      const int nMom,
		      const int block_size_mom_proj,
		      void **host_quark,
		      void **host_quark_bar,
		      const double _Complex *host_mom,
		      void *ret_arr,
		      const int X[4])
{
  const size_t nSp = X[0]*X[1]*X[2] ;
  const size_t nSites = nSp*X[3] ;
  double _Complex loc_sum[ X[3]*nMom*n1*n2 ] = {} ;
  
  // loop over all sites
#pragma omp parallel for reduction(+:loc_sum[:X[3]*nMom*n1*n2])
  for( size_t i = 0 ; i < nSites ; i++ ) {
    const size_t spidx = i%(nSp) , T = i/nSp ;
    const double _Complex *q1[n1] , *q2[n2] ;
    for( int dil1 = 0 ; dil1 < n1 ; dil1++ ) {
      q1[dil1] = (const double _Complex*)host_quark_bar[dil1] + i*3 ;
    }
    for( int dil2 = 0 ; dil2 < n2 ; dil2++ ) {
      q2[dil2] = (const double _Complex*)host_quark[dil2] + i*3 ;
    }
    const double _Complex *Mom = host_mom + spidx ;
    for( int dil1 = 0 ; dil1 < n1 ; dil1++ ) {

      // diagonal guy
      const double _Complex iP =		\
	conj(q1[dil1][0])*(q2[dil1][0])+
	conj(q1[dil1][1])*(q2[dil1][1])+
	conj(q1[dil1][2])*(q2[dil1][2]) ;
      // momentum project
      for( int p = 0 ; p < nMom ; p++ ) {
	loc_sum[ T + X[3]*( p + nMom*( dil1 + n2*dil1 ) ) ] += Mom[ p*nSp ]*iP ;
      }

      for( int dil2 = dil1+1 ; dil2 < n2 ; dil2++ ) {
	// do the inner product over color
	const double _Complex iP =		\
	  conj(q1[dil1][0])*(q2[dil2][0])+
	  conj(q1[dil1][1])*(q2[dil2][1])+
	  conj(q1[dil1][2])*(q2[dil2][2]) ;
	// momentum project
	for( int p = 0 ; p < nMom ; p++ ) {
	  loc_sum[ T + X[3]*( p + nMom*( dil2 + n2*dil1 ) ) ] += Mom[ p*nSp ]*iP ;
	  loc_sum[ T + X[3]*( p + nMom*( dil1 + n2*dil2 ) ) ] += Mom[ p*nSp ]*conj(iP) ;
	}
      }
    }
  }  
  memcpy( ret_arr , loc_sum , X[3]*n1*n2*nMom*sizeof(double _Complex) ) ;
}

static void
set_constant( std::vector<LattField> &laphEigvecs )
{
  int myrank = 0 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_rank( MPI_COMM_WORLD , &myrank ) ;
#endif
  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt ;
  // give them some bullshit values
  for( size_t n = 0 ; n < laphEigvecs.size() ; n++ ) {
    #ifdef ALL_CONSTANT
    std::complex<double> z( (n+1) , (n+1) ) ;
    setConstantField( laphEigvecs[n], z );
    #else
    std::complex<double> *ptr = (std::complex<double>*)laphEigvecs[n].getDataPtr() ;
    const size_t V = (size_t)LayoutInfo::getRankLatticeNumSites() ;
    for( size_t i = 0 ; i < V ; i++ ) {
      for( size_t c = 0 ; c < (size_t)FieldNcolor ; c++ ) {
	#ifdef PSEUDOCONSTANT
	const std::complex<double> z( n+laphEigvecs.size()*(c+FieldNcolor*(i+V*myrank))+1 ,
				      n+laphEigvecs.size()*(c+FieldNcolor*(i+V*myrank))+1 ) ;
	#else
	const std::complex<double> z( unif(mt) , unif(mt) ) ;
	#endif
	*ptr = z ;
	ptr++ ;
      }
    }
    #endif
  }
}

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
  
  // call rephase here
  const int Nev = 8 ;
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);
  set_constant( laphEigvecs ) ;

  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 2 , blockSizeMomProj = 1 , n1 = 2 , n2 = 2 ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int nspat  = X[0]*X[1]*X[2] ;

  // host_mom should be complex
  double _Complex host_mom[ nmom*nspat ] = {} ;
  // host_mom should be complex                                                                                                           
  for( size_t p = 0 ; p < nmom ; p++ ) {
    const size_t mom[3] = { p+1, p+2, p+3 } ;
    for( size_t i = 0 ; i < (size_t)nspat ; i++ ) {
      const size_t x = i%X[0] ;
      const size_t y = (i/X[0])%X[1] ;
      const size_t z = (i/(X[0]*X[1]))%X[2] ;
      const double arg = 2*M_PI*( mom[0]*x/(double)X[0] +
                                  mom[1]*y/(double)X[1] +
                                  mom[2]*z/(double)X[2] ) ;
      double s = 1 ,c = 1 ;
      sincos( arg , &s , &c ) ;
      host_mom[ i + nspat*p ] = c+I*s ;
    }
  }
  
  double _Complex GPU_ret[ n1*n2*nmom*X[3] ] = {} ;
  laphCurrentKernel( n1, n2,
		     nmom,
		     blockSizeMomProj,
		     evList.data() , 
		     evList.data() ,
		     host_mom ,
		     GPU_ret ,
		     X ) ;
  // GPU version
  StopWatch gpu ;
  gpu.start() ;
  laphCurrentKernel( n1, n2,
		     nmom,
		     blockSizeMomProj,
		     evList.data() , 
		     evList.data() ,
		     host_mom ,
		     GPU_ret ,
		     X ) ;
  gpu.stop();
  printLaph(make_strf("\nGPU current kernel in = %g seconds\n", gpu.getTimeInSeconds()));

  // CPU version
  StopWatch cpu ;
  cpu.start() ;
  double _Complex CPU_ret[ n1*n2*nmom*X[3] ] = {} ;
  cpu_code( n1, n2,
	    nmom,
	    blockSizeMomProj,
	    evList.data() , 
	    evList.data() ,
	    host_mom ,
	    CPU_ret ,
	    X ) ;
  cpu.stop() ;
  printLaph(make_strf("\nCPU current kernel in = %g seconds\n", cpu.getTimeInSeconds()));

  for( int dil1 = 0 ; dil1 < n1 ; dil1++ ) {
    for( int dil2 = 0 ; dil2 < n2 ; dil2++ ) {
      for( int p = 0 ; p < nmom ; p++ ) {
	double sum = 0.0 ;
	#ifdef VERBOSE
	printf( "(di1,dil2,p) %d,%d,%d\n" , dil1, dil2, p) ;
	#endif
	for( int t = 0 ; t < X[3] ; t++ ) {
	  #ifdef VERBOSE
	  printf( "(%f %f) == (%f %f)\n" ,
		  creal( CPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] ) ,
		  cimag( CPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] ) ,
		  creal( GPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] ) ,
		  cimag( GPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] )
		  ) ;
	  #endif
	  sum += cabs( CPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] - GPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] ) ; 
	}
	printf( "diff %e\n" , sum ) ;
      }
    }
  }
  finalize();

  return 0;
}
