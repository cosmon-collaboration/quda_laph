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
#define BLASFUNCS

// does A_{ij} = B_{ik}C_{jk}
static inline void
innerevprod( const double _Complex *B ,
	     const double _Complex *C ,
	     const size_t nEv ,
	     double _Complex *A )
{
#ifdef BLASFUNCS
  double _Complex alpha = 1.0 , beta  = 0.0 ;
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

// does A_{ij} = (B_{ik}C_{jk}^\dagger )^\dagger = C B^\dagger
// the reasoning is that when we do a lot of traces we want to cache-cohere
static inline void
innerevprod_dag( const double _Complex *B ,
		 const double _Complex *C ,
		 const size_t nEv ,
		 double _Complex *A )
{
#ifdef BLASFUNCS
  double _Complex alpha = 1.0 , beta  = 0.0 ;
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans ,
	      nEv, nEv, nEv , &alpha , C , nEv , B , nEv , &beta , A , nEv ) ;
#else
  // slow loopy version
  for( size_t i = 0 ; i < nEv ; i++ ) {
    for( size_t j = 0 ; j < nEv ; j++ ) {
      double _Complex sum = 0. ;
      for( size_t k = 0 ; k < nEv ; k++ ) {
	sum += C[ k + nEv*i ]*conj( B[ k + nEv*j ] ) ;
      }
      A[ j + i*nEv ] = sum ;
    }
  }
#endif
}

// fwd guy is \phi(p,tsrc).\tau^{\alpha\beta}_{jk}(t)
static void
compute_Mfwd1( const double _Complex *phi ,
	       const double _Complex *tau ,
	       const size_t nEv ,
	       const size_t nmom ,
	       const size_t tsrc ,
	       const size_t LT ,
	       double _Complex *M )
{
#pragma omp parallel for
  for( size_t p = 0 ; p < nmom ; p++ ) {
    const double _Complex *PhiP = phi + nEv*nEv*( tsrc + LT*p ) ;
    for( size_t t = 0 ; t < LT ; t++ ) {
      for( size_t ab = 0 ; ab < 16 ; ab++ ) {
	const double _Complex *Tp = tau + nEv*nEv*( ab + 16*t ) ;
	double _Complex *Mp = M + nEv*nEv*( ab + 16*( t + LT*p ) ) ;
	innerevprod( PhiP , Tp , nEv , Mp ) ;
      }
    }
  }
}

// fwd guy is \phi(p,t).(\gamma_5 (\tau^{\alpha\beta}_{jk}(t))^\dagger \gamma_5)
static void
compute_Mbwd1( const double _Complex *phi ,
	       const double _Complex *tau ,
	       const size_t nEv ,
	       const size_t nmom ,
	       const size_t tsrc ,
	       const size_t LT ,
	       double _Complex *M )
{
  // shuffle matrix \gamma_5 \tau^\dagger \gamma_5 == { D* , B* , C* , A* } block shuffle just a swap of D and A and 
  const int shuf[16] = { 10, 11, 2, 3,
                         14, 15, 6, 7,
			  8,  9, 0, 1,
			 12, 13, 4, 5 } ;
  const double _Complex *B = tau ;
#pragma omp parallel for
  for( size_t p = 0 ; p < nmom ; p++ ) {
    const double _Complex *PhiP = phi + nEv*nEv*( tsrc + LT*p ) ;
    for( size_t t = 0 ; t < LT ; t++ ) {
      for( size_t ab = 0 ; ab < 16 ; ab++ ) {
	const double _Complex *Tp = tau + nEv*nEv*( shuf[ab] + 16*t ) ;
	double _Complex *Mp = M + nEv*nEv*( shuf[ab] + 16*( t + LT*p ) ) ;
	innerevprod_dag( PhiP , Tp , nEv , Mp ) ;
      }
    }
  }
}

// trace of the product of ev indices Tr[ B.C ]
static inline double _Complex
ev_traceprod( const double _Complex *B ,
	      const double _Complex *C ,
	      const size_t nEv )
{
  double _Complex sum = 0.0 ;
  for( size_t i = 0 ; i < nEv ; i++ ) {
    #ifdef BLASFUNCS
    sum += cblas_zdotc( nEv , B+i*nEv , 1 , C+i*nEv , 1 ) ; 
    #else
    for( size_t j = 0 ; j < nEv ; j++ ) {
      sum += B[j+i*nEv] * conj( C[j+i*nEv] ) ;
    }
    #endif
  }
  return sum ;
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
#pragma omp parallel for collapse(4)
  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
      for( size_t t = 0 ; t < LT ; t++ ) {
	for( size_t osi = 0 ; osi < 4*4*4*4 ; osi++ ) {
	  #if 1
	  const size_t ab = osi/16 , kd = osi%16 ;
	  const double _Complex *Mp1 = Mfwd + nEv*nEv*( ab + 16*( t + LT*psrc ) ) ;
	  const double _Complex *Mp2 = Mbwd + nEv*nEv*( kd + 16*( t + LT*psnk ) ) ;
	  C[ osi + 256*( t + LT*( psnk + nmom*psrc ) ) ] = ev_traceprod( Mp1 , Mp2 , nEv ) ;
	  #else
	  const size_t alpha = (osi/64) ;
	  const size_t beta  = (osi/16)%4 ;
	  const size_t kappa = (osi/8)%4 ;
	  const size_t delta = (osi)%4 ;
	  const double _Complex *Mp1 = Mfwd + nEv*nEv*( alpha + 4*( beta  + 4*( t + LT*psrc ) ) ) ;
	  const double _Complex *Mp2 = Mbwd + nEv*nEv*( delta + 4*( kappa + 4*( t + LT*psnk ) ) ) ;

	  C[ osi + 256*( t + LT*( psnk + nmom*psrc ) ) ] = ev_traceprod( Mp1 , Mp2 , nEv ) ;
	  #endif
	}
      }
    }
  }
}

// Quda interface

//#define VERBOSE_COMPARISON
//#define GPU_STRESS
#define CPU_CROSSCHECK

static void
H2Dwrap( void *d_buf , const double _Complex *buf , const size_t arr_size , const int precision )
{
  if( precision == QUDA_SINGLE_PRECISION ) {
    float _Complex *fbuf = (float _Complex*)malloc( 2*arr_size*precision ) ;
    for( size_t i = 0 ; i < arr_size ; i++ ) {
      fbuf[i] = (float _Complex)buf[i] ;
    }
    qudaMemcpy(d_buf, fbuf, 2*arr_size*precision, qudaMemcpyHostToDevice);
    free( fbuf ) ;
  } else {
    qudaMemcpy(d_buf, buf , 2*arr_size*precision, qudaMemcpyHostToDevice);
  }
}

static void
D2Hwrap( double _Complex *buf , const void *d_buf , const size_t arr_size , const int precision )
{
  if( precision == QUDA_SINGLE_PRECISION ) {
    float _Complex *fbuf = (float _Complex*)malloc( 2*arr_size*precision ) ;
    qudaMemcpy(fbuf, d_buf, 2*arr_size*precision, qudaMemcpyDeviceToHost);
    for( size_t i = 0 ; i < arr_size ; i++ ) {
      buf[i] = (double _Complex)fbuf[i] ;
    }
    free( fbuf ) ;
  } else {
    qudaMemcpy(buf, d_buf , 2*arr_size*precision, qudaMemcpyDeviceToHost);
  }
}

// GPU code
// Does T_ij[0] S_jk(t) ( S_{kl} T_{li} )^\dagger or something
static void
McontractGPU( double _Complex *sum , 
	      const double _Complex *T ,
	      const std::array<const double _Complex*,2> &S ,
	      const size_t LT ,
	      const size_t nEv ,
	      const size_t nmom ,
	      const int precision = QUDA_DOUBLE_PRECISION )
{
  // sbytes
  const size_t sbytes = LT*4*4*nEv*nEv*2*precision ;
  void *d_S = pool_device_malloc( sbytes ) ;
  
  // compute the fwd "elemental"
  const size_t mbytes = LT*nmom*nEv*nEv*16*2*precision ;
  void *d_fwd = pool_device_malloc( mbytes ) ;
  void *d_bwd = pool_device_malloc( mbytes ) ;

  // device return bytes
  const size_t rbytes = LT*nmom*nmom*2*precision ;
  void *d_ret = pool_device_malloc( rbytes ) ;
  
  // tbytes
  const size_t tbytes = LT*nmom*nEv*nEv*2*precision ;
  void *d_T = pool_device_malloc( tbytes ) ;
  H2Dwrap( d_T , T , LT*nmom*nEv*nEv , precision ) ;
  
  // fwd is T_{ij}(0,p) S_{jk}^{\alpha\beta}(t)
  for( size_t t = 0 ; t < LT ; t++ ) {
    H2Dwrap( (char*)d_S + nEv*nEv*16*t*2*precision ,
	     S[0] + nEv*nEv*16*t ,
	     nEv*nEv*16, precision) ;
  }

  // BLAS function1 creating fwd elemental - CPU GPU agree
  // T_{p,i,j}S_{t,alpha,beta,j,k} = M_{p,t,alpha,beta,i,k}
  QudaBLASParam cublas_param1 = newQudaBLASParam();
  cublas_param1.trans_a = QUDA_BLAS_OP_N ;
  cublas_param1.trans_b = QUDA_BLAS_OP_N ;
  cublas_param1.m = nEv ;
  cublas_param1.n = nEv ;
  cublas_param1.k = nEv ;
  cublas_param1.lda = nEv ;
  cublas_param1.ldb = nEv ;
  cublas_param1.ldc = nEv ;
  cublas_param1.a_stride = 0 ;
  cublas_param1.b_stride = nEv*nEv ;
  cublas_param1.c_stride = nEv*nEv ;
  cublas_param1.batch_count = LT*16 ; // do "T x spin x spin" batches
  cublas_param1.alpha = 1.0 ; cublas_param1.beta = 0.0 ;
  cublas_param1.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param1.data_type = ( precision == QUDA_SINGLE_PRECISION ) ?	\
    QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;

  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    blas_lapack::native::stridedBatchGEMM( (char*)d_T + 2*precision*nEv*nEv*psrc*LT ,
					   d_S ,
					   (char*)d_fwd + 2*precision*nEv*nEv*16*psrc*LT ,
					   cublas_param1, QUDA_CUDA_FIELD_LOCATION);
  }

#ifdef DIRECT_COMPARISON
  double _Complex *Mfwd1 = (double _Complex*)malloc( 16*nmom*LT*nEv*nEv*sizeof( double _Complex ) ) ;
  double _Complex *Mfwd2 = (double _Complex*)malloc( 16*nmom*LT*nEv*nEv*sizeof( double _Complex ) ) ;

  compute_Mfwd1( T , S[0] , nEv , nmom , 0 , LT , Mfwd1 ) ;
  D2Hwrap( Mfwd2 , d_fwd , nEv*nEv*16*nmom*LT , precision ) ;

  // comparator
  #if 0
  for( size_t t = 0 ; t < LT ; t++ ) {
    for( size_t p = 0 ; p < nmom ; p++ ) {
      for( size_t alpha = 0 ; alpha < 4 ; alpha++ ) {
	for( size_t beta = 0 ; beta < 4 ; beta++ ) {
	  for( size_t i = 0 ; i < nEv ; i++ ) {
	    for( size_t j = 0 ; j < nEv ; j++ ) {

	      const size_t idx1 = j + nEv*( i + nEv*( beta + 4*( alpha + 4*( t + LT*p ) ) ) ) ;
	      const size_t idx2 = j + nEv*( i + nEv*( beta + 4*( alpha + 4*( p + nmom*t ) ) ) ) ;
	      const double _Complex a = Mfwd1[idx1] , b = Mfwd2[idx1] ;
	      printf( "[t,p,alpha,beta,i,j] (%zu %zu %zu %zu %zu %zu) :: CPU GPU ( %e %e ) == ( %e %e ) || %e\n" ,
		      t , p , alpha , beta , i , j ,
		      creal( a ) , cimag( a ) , creal( b ) , cimag( b ) ,
		      cabs( a - b ) ) ;

	    }
	  }
	}
      }
    }
  } 
  #endif
#endif
  
  // bwd is S2_{jk}^{\alpha\beta}(t) T_{t,p,k,i}
  for( size_t t = 0 ; t < LT ; t++ ) {
    H2Dwrap( (char*)d_S + nEv*nEv*16*t*2*precision ,
	     S[1] + nEv*nEv*16*t ,
	     nEv*nEv*16, precision) ;
  }

  // OK now we do S2_{t,alpha,beta,j,k} T_{p,t,k,i}
  for( size_t ab = 0 ; ab < 16 ; ab++ ) {
    QudaBLASParam cublas_param2 = newQudaBLASParam();
    cublas_param2.trans_a = QUDA_BLAS_OP_N ;
    cublas_param2.trans_b = QUDA_BLAS_OP_N ;
    cublas_param2.m = nEv ;
    cublas_param2.n = nEv ;
    cublas_param2.k = nEv ;
    cublas_param2.lda = nEv ;
    cublas_param2.ldb = nEv ;
    cublas_param2.ldc = nEv ;
    cublas_param2.a_stride = nEv*nEv*16 ;
    cublas_param2.b_stride = nEv*nEv ;
    cublas_param2.c_stride = nEv*nEv*16 ;
    cublas_param2.batch_count = LT ; // do "T" batches
    cublas_param2.alpha = 1.0 ; cublas_param2.beta = 0.0 ;
    cublas_param2.data_order = QUDA_BLAS_DATAORDER_ROW;
    cublas_param2.data_type = ( precision == QUDA_SINGLE_PRECISION ) ?	\
      QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
  
    for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
      blas_lapack::native::stridedBatchGEMM( (char*)d_S + 2*precision*nEv*nEv*ab ,
					     (char*)d_T + 2*precision*nEv*nEv*LT*psnk ,
					     (char*)d_bwd + 2*precision*nEv*nEv*(ab+16*psnk*LT) ,
					     cublas_param2, QUDA_CUDA_FIELD_LOCATION);
    }
  }

#ifdef DIRECT_COMPARISON
  double _Complex *Mbwd1 = (double _Complex*)malloc( 16*nmom*LT*nEv*nEv*sizeof( double _Complex ) ) ;
  double _Complex *Mbwd2 = (double _Complex*)malloc( 16*nmom*LT*nEv*nEv*sizeof( double _Complex ) ) ;

  for( size_t p = 0 ; p < nmom ; p++ ) {
    for( size_t t = 0 ; t < LT ; t++ ) {
      const double _Complex *PhiP = T + nEv*nEv*( t + LT*p ) ;
      for( size_t ab = 0 ; ab < 16 ; ab++ ) {
	const double _Complex *Tp = S[1] + nEv*nEv*( ab + 16*t ) ;
	double _Complex *Mp = Mbwd1 + nEv*nEv*( ab + 16*( t + LT*p ) ) ;
	innerevprod( Tp , PhiP , nEv , Mp ) ;
      }
    }
  }

  D2Hwrap( Mbwd2 , d_bwd , nEv*nEv*16*nmom*LT , precision ) ;

  // comparator
  #if 0
  for( size_t t = 0 ; t < LT ; t++ ) {
    for( size_t p = 0 ; p < nmom ; p++ ) {
      for( size_t alpha = 0 ; alpha < 4 ; alpha++ ) {
	for( size_t beta = 0 ; beta < 4 ; beta++ ) {
	  for( size_t i = 0 ; i < nEv ; i++ ) {
	    for( size_t j = 0 ; j < nEv ; j++ ) {

	      const size_t idx1 = j + nEv*( i + nEv*( beta + 4*( alpha + 4*( t + LT*p ) ) ) ) ;
	      const size_t idx2 = j + nEv*( i + nEv*( beta + 4*( alpha + 4*( p + nmom*t ) ) ) ) ;
	      const double _Complex a = Mbwd1[idx1] , b = Mbwd2[idx1] ;
	      printf( "[t,p,alpha,beta,i,j] (%zu %zu %zu %zu %zu %zu) :: CPU GPU ( %e %e ) == ( %e %e ) || %e\n" ,
		      t , p , alpha , beta , i , j ,
		      creal( a ) , cimag( a ) , creal( b ) , cimag( b ) ,
		      cabs( a - b ) ) ;

	    }
	  }
	}
      }
    }
  }
  #endif
#endif
  
  // innerproduct trace over all psrc psnk
  // fwd_{t,psrc,alpha,beta,i,j} bwd_{t,psnk,beta,alpha,j,i}
  QudaBLASParam cublas_param3 = newQudaBLASParam();
  cublas_param3.trans_a = QUDA_BLAS_OP_N ;
  cublas_param3.trans_b = QUDA_BLAS_OP_C ;
  cublas_param3.m = nmom ;
  cublas_param3.n = nmom ;
  cublas_param3.k = nEv*nEv*4*4 ;
  cublas_param3.lda = nEv*nEv*4*4 ;
  cublas_param3.ldb = nEv*nEv*4*4 ;
  cublas_param3.ldc = nmom ;
  cublas_param3.a_stride = nmom*4*4*nEv*nEv ;
  cublas_param3.b_stride = nmom*4*4*nEv*nEv ;
  cublas_param3.c_stride = nmom*nmom ;
  cublas_param3.batch_count = LT ; // do "T x spin x spin" batches
  cublas_param3.alpha = 1.0 ; cublas_param3.beta = 0.0 ;
  cublas_param3.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param3.data_type = ( precision == QUDA_SINGLE_PRECISION ) ?	\
    QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;

  blas_lapack::native::stridedBatchGEMM( d_fwd, d_bwd, d_ret,
					 cublas_param3,
					 QUDA_CUDA_FIELD_LOCATION);
  D2Hwrap( sum , d_ret , nmom*nmom*LT , precision ) ;

#ifdef DIRECT_COMPARISON
  double _Complex *C = (double _Complex*)malloc( nmom*nmom*LT*sizeof( double _Complex ) ) ;
  // test inner product code
  for( size_t t = 0 ; t < LT ; t++ ) {
    for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
      for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
	// inner
	double _Complex sm = 0 ;
	for( size_t i = 0 ; i < nEv*nEv*16 ; i++ ) {
	  const double _Complex a = Mfwd1[ i + nEv*nEv*16*( psrc + nmom*t ) ] ;
	  const double _Complex b = Mbwd1[ i + nEv*nEv*16*( psnk + nmom*t ) ] ;
	  sm += a*conj(b) ;
	}
	C[ t + LT*( psnk + nmom*psrc ) ] = sm ;

	const double _Complex a = sum[ psnk + nmom*( psrc + nmom*t ) ] ; 
	printf( "(t,psrc,psnk) (%zu %zu %zu) CPU GPU (%e %e) || (%e %e) :: %e\n" ,
		t , psnk , psrc ,
		creal( sm ) , cimag( sm ) ,
		creal( a ) , cimag( a ) ,
		cabs( a - sm ) ) ;
      }
    }
  }

  free( Mfwd1 ) ;
  free( Mfwd2 ) ;
  free( Mbwd1 ) ;
  free( Mbwd2 ) ;
  free( C ) ;
#endif
  
  pool_device_free( d_ret ) ;
  pool_device_free( d_fwd ) ;
  pool_device_free( d_bwd ) ;
  pool_device_free( d_S ) ;
  pool_device_free( d_T ) ;
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
  setVerbosityQuda(QUDA_VERBOSE, "#" , stdout ) ;

  // test for this many EVs
  const int nEv = 64 ;
  const int nmom = 42 ;
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

  const std::array<const double _Complex*,2> per = { tau , tau  } ;
  // GPU
  double _Complex C2[ nmom*nmom*X[3] ] ;
  McontractGPU( C2 , phi , per , X[3] , nEv , nmom , QUDA_DOUBLE_PRECISION ) ;

  StopWatch GPU ; GPU.start() ;
  McontractGPU( C2 , phi , per , X[3] , nEv , nmom , QUDA_DOUBLE_PRECISION ) ;
  GPU.stop() ;
  double GPUtime = GPU.getTimeInSeconds() ;
  printLaph( make_strf( "\n GPU in %g seconds\n" , GPUtime ) ) ;
  
  // create a meson field or whatever bullshit people call this object for a specific psrc,psnk pair
  // M^{\alpha\beta}_{ij}(p,t)
  double _Complex *Mfwd = (double _Complex*)malloc( 16*nmom*X[3]*nEv*nEv*sizeof( double _Complex ) ) ;
  double _Complex *Mbwd = (double _Complex*)malloc( 16*nmom*X[3]*nEv*nEv*sizeof( double _Complex ) ) ;

  StopWatch fwd ; fwd.start() ;
  compute_Mfwd1( phi , tau , nEv , nmom , 0 , X[3] , Mfwd ) ;
  fwd.stop() ;
  double CPUtime = fwd.getTimeInSeconds() ;
  printLaph( make_strf( "\n Mfwd in %g seconds\n" , CPUtime ) ) ;
  
  StopWatch bwd ; bwd.start() ;
  compute_Mbwd1( phi , tau , nEv , nmom , 0 , X[3] , Mbwd ) ;
  bwd.stop() ;
  CPUtime += bwd.getTimeInSeconds() ;
  printLaph( make_strf( "\n Mbwd in %g seconds\n" , bwd.getTimeInSeconds() ) ) ;

  // so now the contraction is -tr[ G1.mfwd.(gt.G2^\dagger.gt).G5.Mbwd.G5 ] over all indices but we can cheat by closing the ev indices first T^{alpha,beta,kappa,delta} = M^{alpha,beta} \bar{M^{kappa,delta} hmmmm
  free( tau ) ;
  free( phi ) ;

  StopWatch traces ; traces.start() ;
  // KISS look at the gamma_5 first Tr[ M M\bar ] all indices closed
  double _Complex *C = (double _Complex*)malloc( nmom*nmom*X[3]*256*sizeof(double _Complex) );

  contract_meson( Mfwd , Mbwd , nEv , nmom , X[3] , C ) ;

  for( size_t t = 0 ; t < X[3] ; t++ ) {
    for( size_t osi = 0 ; osi < 1 ; osi++ ) {
      printf( "%zu %zu %e %e\n" , t , osi ,
	      creal(C[osi + 256*t]) , cimag(C[osi + 256*t]) ) ;
    }
  }

  traces.stop() ;

  CPUtime += traces.getTimeInSeconds() ;
  printLaph( make_strf( "\n CPU contraction time in %g seconds\n" , CPUtime ) ) ;
  free( C ) ;

  free( Mfwd ) ;
  free( Mbwd ) ;
  
  finalize( ) ;
  return 0;
}
