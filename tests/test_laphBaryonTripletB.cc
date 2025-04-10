/**
   Test the laphBaryonTripletB code in my fork of quda
 */
#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

#include <random>
#include <cassert>
#include <complex.h>

using namespace LaphEnv ;
using namespace quda ;

// ok so what is this doing?
void
cpu_code( const int n1,
	  const int n2,
	  const int n3,
	  const int nMom,
	  const int nEv,
	  const double _Complex *host_coeffs1,
	  const double _Complex *host_coeffs2,
	  const double _Complex *host_coeffs3,
	  const double _Complex *host_mode_trip_buf,
	  double _Complex *host_ret_arr )
{
  int global = 1 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
#endif
  const int nSubEv = nEv/global ;

  // zgemm d_mtb*d_coeffs = d_q3 
  double _Complex mtb[ nMom*nSubEv*nEv*nEv ] = {} ;
  memcpy( mtb , host_mode_trip_buf , nMom*nSubEv*nEv*nEv*sizeof( double _Complex ) ) ;

  // is the first product
  double _Complex q3[ nMom*nSubEv*nEv*n3 ] = {} ;
  
  // eg this is doing C_ij \sum_k A_ik Bjk strided for nMom*nSubEv
  // C_ij is size Nmom
  // nmom*nSubEv*Ev*Ev | n3*nEv
  for( int p = 0 ; p < nMom ; p++ ) {
    for( int n = 0 ; n < nSubEv ; n++ ) {
      for( int evo = 0 ; evo < nEv ; evo++ ) { 
	double _Complex *pt1 = (double _Complex*)mtb + nEv*( evo + nEv*( n + nSubEv*p ) ) ;
	for( int nin = 0 ; nin < n3 ; nin++ ) {
	  double _Complex *pt2 = (double _Complex*)host_coeffs3+nin*nEv ;
	  double _Complex sum = 0.0 ;
	  for( int ev = 0 ; ev < nEv ; ev++ ) {
	    sum += pt1[ev]*pt2[ev] ; 
	  }
	  q3[nin+n3*( evo + nEv*( n + nSubEv*p ) ) ] = sum ; 
	}
      }
    }
  }

  // second set of gemms
  // host_coeffs2*d_q3 -> dtmp
  // n1*nEv . nMom*nSubEv*nEv*n3 -> nSubEv*n1*n3 
  for( int p = 0 ; p < nMom ; p++ ) {
    // again just a product over nEv
    double _Complex tmp[ nSubEv*n2*n3 ] = {} ;
    for( int nsub = 0 ; nsub < nSubEv ; nsub++ ) {
      for( int nout = 0 ; nout < n2 ; nout++ ) {
	for( int nin = 0 ; nin < n3 ; nin++ ) { 
	  double _Complex sum = 0 ;
	  for( int ev = 0 ; ev < nEv ; ev++ ) {
	    sum += host_coeffs2[ ev + nout*nEv ] * q3[ nin + n3*( ev + nEv*( nsub + nSubEv*p ) ) ] ; 
	  }
	  tmp[ nin + n3*( nout + n2*nsub ) ] = sum ;
	}
      }
    }
    // host_coeffs1*d_tmp -> d_ret
    // n1*nEv . nSubEv*n2*n3 -> nmom.n1.n2.n3
    for( int i = 0 ; i < n1 ; i++ ) {
      for( int j = 0 ; j < n2 ; j++ ) {
	for( int k = 0 ; k < n3 ; k++ ) {
	  // only works for 1 rank so that nSubEv == nEv 
	  double _Complex sum = 0.0 ;
	  for( int ev = 0 ; ev < nEv ; ev++ ) {
	    sum += host_coeffs1[ev+i*nEv]*tmp[ k + n3*( j + n2*ev ) ];
	  }
	  host_ret_arr[ k + n3*( j + n2*( i + n1*p ) ) ] = sum ;
	}
      }
    }
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
  assert( global == 1 ) ; // for now

  const int Nev = 16 ;
  const int NsubEv = 16/global ;
  const int nmom = 4 , n1 = 1 , n2 = 2 , n3 = 3 ;

  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt ;
  
  double _Complex host_coeffs1[ n1*Nev ] = {} ;
  for( int i = 0 ; i < n1*Nev ; i++ ) {
    host_coeffs1[ i ] = 1. ; //unif(mt) + I*unif(mt) ;
  }  
  double _Complex host_coeffs2[ n2*Nev ] = {} ;
  for( int i = 0 ; i < n2*Nev ; i++ ) {
    host_coeffs2[ i ] = 1. ; //unif(mt) + I*unif(mt) ;
  }
  double _Complex host_coeffs3[ n3*Nev ] = {} ;
  for( int i = 0 ; i < n3*Nev ; i++ ) {
    host_coeffs3[ i ] = 1. ; //unif(mt) + I*unif(mt) ;
  }
  double _Complex host_mode_trip_buf[ NsubEv*Nev*Nev*nmom ] = {} ;
  for( int i = 0 ; i < NsubEv*Nev*Nev*nmom ; i++ ) {
    host_mode_trip_buf[ i ] = 1. ; //unif(mt) + I*unif(mt) ;
  }

  double _Complex host_ret_arr[ nmom*n1*n2*n3 ] = {} ;
  
  laphBaryonKernelComputeModeTripletB( n1, n2, n3,
				       nmom,
				       Nev ,
				       host_coeffs1 ,
				       host_coeffs2 ,
				       host_coeffs3 ,
				       host_mode_trip_buf ,
				       host_ret_arr ) ;

  double _Complex cpu_ret_arr[ nmom*n1*n2*n3 ] = {} ;

  cpu_code( n1, n2, n3,
	    nmom,
	    Nev ,
	    host_coeffs1 ,
	    host_coeffs2 ,
	    host_coeffs3 ,
	    host_mode_trip_buf ,
	    cpu_ret_arr ) ;

  for( int p = 0 ; p < nmom ; p++ ) {
    double Sum = 0. ;
    for( int dil = 0 ; dil < n1*n2*n3 ; dil++ ) {
      Sum += cabs( cpu_ret_arr[ dil + n1*n2*n3*p ] - host_ret_arr[ dil + n1*n2*n3*p ] ) ;

      printf( "(%f,%f) == (%f,%f)\n" ,
	      creal( cpu_ret_arr[dil+n1*n2*n3*p] ) , cimag( cpu_ret_arr[dil+n1*n2*n3*p] ) ,
	      creal( host_ret_arr[dil+n1*n2*n3*p] ) , cimag( host_ret_arr[dil+n1*n2*n3*p] ) ) ;
    }
    printf( "%d Sub %e\n" , p , Sum ) ;
  }
  
  finalize();

  return 0;
}
