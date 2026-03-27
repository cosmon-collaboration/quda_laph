/**
   Test the feasibility of Baryon contractions on the device


   TODO spin indices, ... etc
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

#include <cuComplex.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>

using namespace quda ;
using namespace LaphEnv ;

// Quda interface

//#define VERBOSE_COMPARISON
//#define GPU_STRESS
#define CPU_CROSSCHECK

// Cpu code
#define BLASFUNCS

// does TL_{ijk} S1_{ii'} S2_{jj'} S3_{kk'} TR_{i'j'k'}
double _Complex
Tcontract1( const double _Complex *TL ,
	    const double _Complex *S1 ,
	    const double _Complex *S2 ,
	    const double _Complex *S3 ,
	    const double _Complex *TR ,
	    const size_t nEv )
{
  double _Complex sum = 0. ;
  for( size_t i = 0 ; i < nEv ; i++ ) {
    for( size_t ip = 0 ; ip < nEv ; ip++ ) {

      double _Complex sumj = 0 ;
      for( size_t j = 0 ; j < nEv ; j++ ) {
	for( size_t jp = 0 ; jp < nEv ; jp++ ) {

	  double _Complex sumk = 0 ;
	  for( size_t k = 0 ; k < nEv ; k++ ) {
	    for( size_t kp = 0 ; kp < nEv ; kp++ ) {
	      sumk += TL[ k + nEv*( j + nEv*i ) ]*S3[ kp + nEv*k ]*TR[ kp + nEv*( jp + nEv*ip ) ] ; 
	    }
	  }

	  sumj += S2[ jp + nEv*j ]*sumk ;

	  // jjp
	}
      }
      sum += S1[ ip + nEv*i ]*sumj ;
      // iip
    }
  }	  
  return sum ;
}


// does TL_{ijk} S1_{ii'} S2_{jj'} S3_{kk'} TR_{i'j'k'}
double _Complex
Tcontract2( const double _Complex *TL ,
	    const double _Complex *S1 ,
	    const double _Complex *S2 ,
	    const double _Complex *S3 ,
	    const double _Complex *TR ,
	    const size_t nEv )
{
  double _Complex sum = 0. ;
  for( size_t i = 0 ; i < nEv ; i++ ) {
    for( size_t ip = 0 ; ip < nEv ; ip++ ) {

      double _Complex sumj = 0. ;
      for( size_t j = 0 ; j < nEv ; j++ ) {
	for( size_t jp = 0 ; jp < nEv ; jp++ ) {

	  double _Complex sumk = 0 ;
	  for( size_t k = 0 ; k < nEv ; k++ ) {
	    double _Complex innerK = 0. ;
	    innerK = cblas_zdotu( nEv , S3 + nEv*k , 1 , TR + nEv*( jp + nEv*ip ) , 1 ) ;
	    sumk += TL[ k + nEv*(j + nEv*i) ]*innerK ;
	  }

	  sumj += S2[ jp + nEv*j ]*sumk ;

	  // jjp
	}
      }
      sum += S1[ ip + nEv*i ]*sumj ;
      // iip
    }
  }
  

	  
  return sum ;
}

// does TL_{ijk} S1_{ii'} S2_{jj'} S3_{kk'} TR_{i'j'k'}
double _Complex
Tcontract3( const double _Complex *TL ,
	    const double _Complex *S1 ,
	    const double _Complex *S2 ,
	    const double _Complex *S3 ,
	    const double _Complex *TR ,
	    const size_t nEv )
{
  double _Complex ini[ nEv ] , inj[ nEv ] , ink[ nEv ] , sum = 0. ;
  for( size_t i = 0 ; i < nEv ; i++ ) {
    for( size_t ip = 0 ; ip < nEv ; ip++ ) {
      ini[ip] = 0.0 ;
      for( size_t j = 0 ; j < nEv ; j++ ) {
	for( size_t jp = 0 ; jp < nEv ; jp++ ) {
	  for( size_t k = 0 ; k < nEv ; k++ ) {
	    ink[ k ] = cblas_zdotu( nEv , S3 + nEv*k , 1 , TR + nEv*( jp + nEv*ip ) , 1 ) ;
	  }
	  inj[jp] = cblas_zdotu( nEv , TL + nEv*(j+nEv*i) , 1 , ink , 1 ) ;
	}
	ini[ip] += cblas_zdotu( nEv , S2 +  nEv*j , 1 , inj , 1 ) ;
      }
    }
    sum += cblas_zdotu( nEv , S1 +  nEv*i , 1 , ini , 1 ) ;
  }
  
  return sum ;
}

typedef enum { ijk , jik , kji , kij , ikj , jki } bartype ;


// does TL_{ijk} S1_{ii'} S2_{jj'} S3_{kk'} TR_{i'j'k'}
double _Complex
Tcontract4( const double _Complex *TL ,
	    const std::array<const double _Complex*,3> &S ,
	    const double _Complex *TR ,
	    const int nEv ,
	    const bartype btype )
{
  double _Complex inj[nEv] , ini[nEv] , sum = 0. ;
  for( int i = 0 ; i < nEv ; i++ ) {
    for( int j = 0 ; j < nEv ; j++ ) {
      for( int k = 0 ; k < nEv ; k++ ) {

	// map of inner indices for perambulators
	std::array<const double _Complex*,3> buf = { S[0] , S[1] , S[2] } ;
	std::array<int,3> pmap = { i , j , k } ;
	switch( btype ){
	case ijk : break ;
	case jik :
	  buf[0] = S[1] ; buf[1] = S[0] ; buf[2] = S[2] ;
	  pmap[0] = j ; pmap[1] = i ; pmap[2] = k ; 
	  break ;
	case kji :
	  buf[0] = S[2] ; buf[1] = S[1] ; buf[2] = S[0] ;
	  pmap[0] = k ; pmap[1] = j ; pmap[2] = i ; 
	  break ;
	case kij :
	  buf[0] = S[2] ; buf[1] = S[0] ; buf[2] = S[1] ;
	  pmap[0] = k ; pmap[1] = i ; pmap[2] = j ; 
	  break ;
	case ikj :
	  buf[0] = S[0] ; buf[1] = S[2] ; buf[2] = S[1] ;
	  pmap[0] = i ; pmap[1] = k ; pmap[2] = j ; 
	  break ;
	case jki :
	  buf[0] = S[1] ; buf[1] = S[2] ; buf[2] = S[0] ;
	  pmap[0] = j ; pmap[1] = k ; pmap[2] = i ; 
	  break ;
	}

	// hot inner products on sink idx
	const double _Complex *pTr = TR ;
	for( int ip = 0 ; ip < nEv ; ip++ ) {
	  for( int jp = 0 ; jp < nEv ; jp++ ) {
	    // kp innerprod
	    inj[jp] = cblas_zdotu( nEv , buf[2] + nEv*pmap[2] , 1 , pTr , 1 ) ;
	    pTr += nEv ;
	  }
	  ini[ip] = cblas_zdotu( nEv , buf[1] + nEv*pmap[1] , 1 , inj , 1 ) ;
	}
	// ijk
	sum += TL[ k + nEv*(j + nEv*i) ]*\
	  cblas_zdotu( nEv , buf[0] + nEv*pmap[0] , 1 , ini , 1 ) ;
      }
    }
  }	  
  return sum ;
}

// does TL_{ijk}(T) S1_{ii'}(t) S2_{jj'}(t) S3_{kk'}(t) TR_{i'j'k'}(t)
static void
Tcontract5( double _Complex *sum , 
	    const double _Complex *TL ,
	    const std::array<const double _Complex*,3> &S ,
	    const double _Complex *TR ,
	    const int LT ,
	    const int nEv ,
	    const bartype btype )
{
  double _Complex inj[LT][nEv] , ini[LT][nEv] ;
  memset( sum , 0. , LT*sizeof( double _Complex ) ) ; 
  for( int i = 0 ; i < nEv ; i++ ) {
    for( int j = 0 ; j < nEv ; j++ ) {
      for( int k = 0 ; k < nEv ; k++ ) {
	  
	// map of inner indices for perambulators
	std::array<const double _Complex*,3> buf =	\
	  { S[0] , S[1] , S[2] } ;
	  std::array<int,3> pmap = { i , j , k } ;
	  switch( btype ){
	  case ijk : break ;
	  case jik :
	    buf[0] = S[1] ; buf[1] = S[0] ; buf[2] = S[2] ;
	    pmap[0] = j ; pmap[1] = i ; pmap[2] = k ; 
	    break ;
	  case kji :
	    buf[0] = S[2] ; buf[1] = S[1] ; buf[2] = S[0] ;
	    pmap[0] = k ; pmap[1] = j ; pmap[2] = i ; 
	    break ;
	  case kij :
	    buf[0] = S[2] ; buf[1] = S[0] ; buf[2] = S[1] ;
	    pmap[0] = k ; pmap[1] = i ; pmap[2] = j ; 
	    break ;
	  case ikj :
	    buf[0] = S[0] ; buf[1] = S[2] ; buf[2] = S[1] ;
	    pmap[0] = i ; pmap[1] = k ; pmap[2] = j ; 
	    break ;
	  case jki :
	    buf[0] = S[1] ; buf[1] = S[2] ; buf[2] = S[0] ;
	    pmap[0] = j ; pmap[1] = k ; pmap[2] = i ; 
	    break ;
	  }

	  // hot inner products on sink idx
	  const double _Complex *pTr = TR ;
	  for( int ip = 0 ; ip < nEv ; ip++ ) {
	    for( int jp = 0 ; jp < nEv ; jp++ ) {

	      // would be a batched blas (1xN) (Nxp) call over all t indices to create an object in[t][p][j']
	      
	      for( int t = 0 ; t < LT ; t++ ) {
		// kp innerprod
		inj[t][jp] = cblas_zdotu( nEv , buf[2] + nEv*( pmap[2] + nEv*16*t ) , 1 , pTr + nEv*nEv*nEv*t , 1 ) ;
	      }
	      pTr += nEv ;
	    }
	    // batched over t (1xN)(Nxp) again
	    for( int t = 0 ; t < LT ; t++ ) {
	      ini[t][ip] = cblas_zdotu( nEv , buf[1] + nEv*(pmap[1] + nEv*16*t) , 1 , inj[t] , 1 ) ;
	    }
	  }
	  for( int t = 0 ; t < LT ; t++ ) {
	    // ijk
	    sum[t] += TL[ k + nEv*(j + nEv*i) ]*			\
	      cblas_zdotu( nEv , buf[0] + nEv*(pmap[0] + nEv*16*t) , 1 , ini[t] , 1 ) ;

	  }
	  // T loop
      }
    }
  }
  return ;
}

// CPU version of GPU code
static void
Tcontract6( double _Complex *sum ,
	    const double _Complex *T ,
	    const std::array<const double _Complex*,3> &S ,
	    const size_t LT ,
	    const size_t nEv ,
	    const size_t nmom ,
	    const bartype btype )
{
  //#pragma omp parallel for
  for( size_t t = 0 ; t < LT ; t++ ) {

    double _Complex *A = (double _Complex*)malloc( nmom*nEv*nEv*nEv*sizeof(double _Complex) ) ;
    double _Complex *B = (double _Complex*)malloc( nmom*nEv*nEv*nEv*sizeof(double _Complex) ) ;
    
    memset( A , 0. , nmom*nEv*nEv*nEv*sizeof( double _Complex ) ) ;
    for( size_t ip = 0 ; ip < nEv ; ip++ ) {
      for( size_t j = 0 ; j < nEv ; j++ ) {
	for( size_t k = 0 ; k < nEv ; k++ ) {
	  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
	    A[ k + nEv*( j + nEv*( ip + nEv*psrc )) ] =		\
	      cblas_zdotu( nEv ,
			   T + k + nEv*j + nEv*nEv*nEv*LT*psrc , nEv*nEv ,
			   S[0] + ip + nEv*(nEv*16*t) , nEv ) ;
	  }
	}
      }
    }

    // compute B
    memset( B , 0. , nmom*nEv*nEv*nEv*sizeof( double _Complex ) ) ;
    for( size_t ip = 0 ; ip < nEv ; ip++ ) {
      for( size_t k = 0 ; k < nEv ; k++ ) {
	for( size_t jp = 0 ; jp < nEv ; jp++ ) {
	  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
	    B[ k + nEv*( jp + nEv*( ip + nEv*psrc )) ] =		\
	      cblas_zdotu( nEv , A + k + nEv*( nEv*(ip + nEv*psrc )) , nEv ,
			   S[1] + jp + nEv*(nEv*16*t) , nEv ) ;
	  }
	}
      }
    }
    
    // get A again
    memset( A , 0. , nmom*nEv*nEv*nEv*sizeof( double _Complex ) ) ;
    for( size_t ip = 0 ; ip < nEv ; ip++ ) {
      for( size_t jp = 0 ; jp < nEv ; jp++ ) {
	for( size_t kp = 0 ; kp < nEv ; kp++ ) {
	  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
	    A[ kp + nEv*( jp + nEv*(ip + nEv*psrc ) ) ] =				\
	      cblas_zdotu( nEv , B + nEv*( jp + nEv*( ip + nEv*psrc ) ) , 1 ,
			   S[2] + kp + nEv*nEv*16*t , nEv ) ;
	  }
	}
      }
    }
    
    // compute the sum
    for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
      for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
	sum[ t + LT*( psnk + nmom*psrc ) ] = \
	  cblas_zdotu( nEv*nEv*nEv , A + nEv*nEv*nEv*psrc , 1 ,
		       T + nEv*nEv*nEv*(t+LT*psnk) , 1 ) ;
      }
    }
    
    free( A ) ;
    free( B ) ;
  }
}

// need a better idea T_{ijk} D_{ii'} -> A_{i'jk}
// A_{i'jk}D_{jj'} -> B_{i'j'k}
// B_{i'j'k}D_{kk'} -> A_{i'j'k'}
// T_{i'j'k'}(t) C_{i'j'k'} = sum

// GPU code
// does TL_{ijk}(T) S1_{ii'}(t) S2_{jj'}(t) S3_{kk'}(t) TR_{i'j'k'}(t)
static void
TcontractGPU( double _Complex *sum , 
	      const double _Complex *T ,
	      const std::array<const double _Complex*,3> &S ,
	      const size_t LT ,
	      const size_t nEv ,
	      const size_t nmom ,
	      const bartype btype )
{
  // device "T" object
  //const size_t tbytes = LT*nmom*nEv*nEv*nEv*2*sizeof(double) ;
  //void *d_T = pool_device_malloc( tbytes ) ;
  // device buffers for perambulators
  const size_t sbytes = LT*16*nEv*nEv*2*sizeof(double) ;
  void *d_buf0 = pool_device_malloc( sbytes ) ;

  const size_t tbytes = LT*nmom*nEv*nEv*nEv*2*sizeof(double) ;
  void *d_T = pool_device_malloc( tbytes ) ;
  void *d_A = pool_device_malloc( tbytes ) ;
  void *d_B = pool_device_malloc( tbytes ) ;

  const size_t retbytes = LT*nmom*nmom*2*sizeof(double) ;
  void *d_ret = pool_device_malloc( retbytes ) ;

  const double OneGB = 1024.*1024.*1024. ;
  printfQuda( "Trip %g GB | buf0 %g GB | ret %g GB | A %g GB | B %g GB \n" ,
	      tbytes/OneGB , sbytes/OneGB , retbytes/OneGB , tbytes/OneGB , tbytes/OneGB ) ;
  
  // quda memcpys
  qudaMemcpy(d_T,    T   , tbytes, qudaMemcpyHostToDevice);    
  qudaMemcpy(d_buf0, S[0], sbytes, qudaMemcpyHostToDevice);

  // D_{ii'}^T(t) T_{ijk}(p,t) -> A_{i'jk}(p,t)
  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    QudaBLASParam cublas_param1 = newQudaBLASParam();
    cublas_param1.trans_a = QUDA_BLAS_OP_T ;
    cublas_param1.trans_b = QUDA_BLAS_OP_N ;
    cublas_param1.m = nEv ;
    cublas_param1.n = nEv*nEv ;
    cublas_param1.k = nEv ;
    cublas_param1.lda = nEv;
    cublas_param1.ldb = nEv*nEv;
    cublas_param1.ldc = nEv*nEv;
    cublas_param1.a_stride = 16*nEv*nEv ;
    cublas_param1.b_stride = 0 ;
    cublas_param1.c_stride = nEv*nEv*nEv ;
    cublas_param1.batch_count = LT ; // do "T" batches
    cublas_param1.alpha = 1.0 ; cublas_param1.beta = 0.0 ;
    cublas_param1.data_order = QUDA_BLAS_DATAORDER_ROW;
    cublas_param1.data_type = QUDA_BLAS_DATATYPE_Z;

    // checked and is ok
    blas_lapack::native::stridedBatchGEMM( d_buf0,
					   (double _Complex*)d_T + nEv*nEv*nEv*LT*psrc ,
					   (double _Complex*)d_A + nEv*nEv*nEv*LT*psrc ,
					   cublas_param1, QUDA_CUDA_FIELD_LOCATION);
  }
  
  // second part break up into t i' batches
  qudaMemcpy(d_buf0, S[1], sbytes, qudaMemcpyHostToDevice);  
  
  // D_{jj'}^T(t) A_{i'jk}(p,t) -> B_{i'jk}
  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    for( size_t t = 0 ; t < LT ; t++ ) {
      QudaBLASParam cublas_param2 = newQudaBLASParam();
      cublas_param2.trans_a = QUDA_BLAS_OP_T ;
      cublas_param2.trans_b = QUDA_BLAS_OP_N ;
      cublas_param2.m = nEv ;
      cublas_param2.n = nEv ;
      cublas_param2.k = nEv ;
      cublas_param2.lda = nEv;
      cublas_param2.ldb = nEv;
      cublas_param2.ldc = nEv;
      cublas_param2.a_stride = 0 ;
      cublas_param2.b_stride = nEv*nEv ;
      cublas_param2.c_stride = nEv*nEv ;
      cublas_param2.batch_count = nEv ; // do "nEv" batches
      cublas_param2.alpha = 1.0 ; cublas_param2.beta = 0.0 ;
      cublas_param2.data_order = QUDA_BLAS_DATAORDER_ROW;
      cublas_param2.data_type = QUDA_BLAS_DATATYPE_Z;

      blas_lapack::native::stridedBatchGEMM( (double _Complex*)d_buf0 + nEv*nEv*16*t ,
					     (double _Complex*)d_A + nEv*nEv*nEv*(t + LT*psrc),
					     (double _Complex*)d_B + nEv*nEv*nEv*(t + LT*psrc),
					     cublas_param2, QUDA_CUDA_FIELD_LOCATION);
    }
  }

  // third product
  qudaMemcpy(d_buf0, S[2], sbytes, qudaMemcpyHostToDevice);
  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    // A_{i'jk'}(p,t) = B_{i'j'k}(p,t)D_{kk'}(t)
    QudaBLASParam cublas_param3 = newQudaBLASParam();
    cublas_param3.trans_a = QUDA_BLAS_OP_N ;
    cublas_param3.trans_b = QUDA_BLAS_OP_N ;
    cublas_param3.m = nEv*nEv ;
    cublas_param3.n = nEv ;
    cublas_param3.k = nEv ;
    cublas_param3.lda = nEv;
    cublas_param3.ldb = nEv;
    cublas_param3.ldc = nEv;
    cublas_param3.a_stride = nEv*nEv*nEv ;
    cublas_param3.b_stride = 16*nEv*nEv ;
    cublas_param3.c_stride = nEv*nEv*nEv ;
    cublas_param3.batch_count = LT ; // do "T" batches
    cublas_param3.alpha = 1.0 ; cublas_param3.beta = 0.0 ;
    cublas_param3.data_order = QUDA_BLAS_DATAORDER_ROW;
    cublas_param3.data_type = QUDA_BLAS_DATATYPE_Z;
    
    // checked and is ok
    blas_lapack::native::stridedBatchGEMM( (double _Complex*)d_B + nEv*nEv*nEv*LT*psrc ,
					   d_buf0 ,
					   (double _Complex*)d_A + nEv*nEv*nEv*LT*psrc ,
					   cublas_param3, QUDA_CUDA_FIELD_LOCATION);
  }
  
  // batched inner product A_{i'j'k'}(psrc,t)T_{i'j'k'}(psnk,t)
  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
      QudaBLASParam cublas_param4 = newQudaBLASParam();
      cublas_param4.trans_a = QUDA_BLAS_OP_N ;
      cublas_param4.trans_b = QUDA_BLAS_OP_T ;
      cublas_param4.m = 1 ;
      cublas_param4.n = 1 ;
      cublas_param4.k = nEv*nEv*nEv ;
      cublas_param4.lda = nEv*nEv*nEv;
      cublas_param4.ldb = nEv*nEv*nEv;
      cublas_param4.ldc = 1;
      cublas_param4.a_stride = nEv*nEv*nEv ;
      cublas_param4.b_stride = nEv*nEv*nEv ;
      cublas_param4.c_stride = 1 ;
      cublas_param4.batch_count = LT ;
      cublas_param4.alpha = 1.0 ; cublas_param4.beta = 0.0 ;
      cublas_param4.data_order = QUDA_BLAS_DATAORDER_ROW;
      cublas_param4.data_type = QUDA_BLAS_DATATYPE_Z;   
      blas_lapack::native::stridedBatchGEMM( (double _Complex*)d_A + nEv*nEv*nEv*LT*psrc ,
					     (double _Complex*)d_T + nEv*nEv*nEv*LT*psnk ,
					     (double _Complex*)d_ret + LT*( psnk + nmom*psrc ) ,
					     cublas_param4, QUDA_CUDA_FIELD_LOCATION);
    }
  }
    
  // return function
  qudaMemcpy( sum, d_ret , retbytes, qudaMemcpyDeviceToHost);
  
#ifdef direct_compare
  double _Complex *A  = (double _Complex*)malloc( nEv*nEv*nEv*LT*sizeof(double _Complex)) ;
  double _Complex *B  = (double _Complex*)malloc( nEv*nEv*nEv*LT*sizeof(double _Complex)) ;
  double _Complex *tT = (double _Complex*)malloc( nEv*nEv*nEv*LT*sizeof(double _Complex)) ;

  double _Complex sum2[ LT ] ;
  for( int t = 0 ; t < LT ; t++ ) {
    memset( A + nEv*nEv*nEv*t , 0. , nEv*nEv*nEv*sizeof( double _Complex ) ) ;
    for( int i = 0 ; i < nEv ; i++ ) {
      for( int j = 0 ; j < nEv ; j++ ) {
	for( int k = 0 ; k < nEv ; k++ ) {
	  for( int ip = 0 ; ip < nEv ; ip++ ) {
	    A[ k + nEv*( j + nEv*(ip + nEv*t ) ) ] +=			\
	      T[ k + nEv*( j + nEv*(i ) ) ] * S[0][ ip + nEv*(i + nEv*16*t) ] ;
	  }
	}
      }
    }
    // compute B
    memset( B + nEv*nEv*nEv*t , 0. , nEv*nEv*nEv*sizeof( double _Complex ) ) ;
    for( int ip = 0 ; ip < nEv ; ip++ ) {
      for( int k = 0 ; k < nEv ; k++ ) {
	for( int j = 0 ; j < nEv ; j++ ) {
	  for( int jp = 0 ; jp < nEv ; jp++ ) {
	    B[ k + nEv*( jp + nEv*( ip + nEv*t )) ] +=			\
	      A[ k + nEv*( j + nEv*(ip + nEv*t )) ] * S[1][ jp + nEv*(j + nEv*16*t) ] ; // + nEv*16*t) ] ;
	  }
	}
      }
    }
    // get A again
    memset( A + nEv*nEv*nEv*t, 0. , nEv*nEv*nEv*sizeof( double _Complex ) ) ;
    for( int ip = 0 ; ip < nEv ; ip++ ) {
      for( int jp = 0 ; jp < nEv ; jp++ ) {
	for( int kp = 0 ; kp < nEv ; kp++ ) {
	  for( int k = 0 ; k < nEv ; k++ ) {
	    A[ kp + nEv*( jp + nEv*(ip + nEv*t) ) ] +=			\
	      B[ k + nEv*( jp + nEv*(ip + nEv*t) ) ] * S[2][ kp + nEv*( k + nEv*16*t) ] ;
	  }
	}
      }
    }
    // compute the sum
    sum2[t] = 0. ;
    for( int ip = 0 ; ip < nEv ; ip++ ) {
      for( int jp = 0 ; jp < nEv ; jp++ ) {
	for( int kp = 0 ; kp < nEv ; kp++ ) {
	  //sum[t] += A[kp+nEv*(jp+nEv*ip)]*TR[kp+nEv*(jp+nEv*(ip) ) ] ; //+nEv*t))] ;
	  sum2[t] += A[kp+nEv*(jp+nEv*(ip+nEv*t))]*T[kp+nEv*(jp+nEv*(ip+nEv*t) ) ] ; //+nEv*t))] ;
	}
      }
    }

  }

  //qudaMemcpy( tT , d_A , inbytes, qudaMemcpyDeviceToHost);

  for( int t = 0 ; t < LT ; t++ ) {
    double _Complex a = sum2[ t ] ;
    double _Complex b = sum[ t ] ;
    double diff = cabs( a - b ) ;
    printf( "(%d %d %d %d) %e | (%e %e) (%e %e)\n" , t , diff ,
	    creal( a ) , cimag( a ) ,
	    creal( b ) , cimag( b ) ) ;
  }
  
  for( int t = 0 ; t < LT ; t++ ) {
    for( int ip = 0 ; ip < nEv ; ip++ ) {
      for( int j = 0 ; j < nEv ; j++ ) {
	for( int k = 0 ; k < nEv ; k++ ) {
	  double _Complex a =  A[ k + nEv*(j+nEv*(ip+nEv*t) ) ] ;
	  double _Complex b = tT[ k + nEv*(j+nEv*(ip+nEv*t) ) ] ;
	  double diff = cabs( a - b ) ;
	  printf( "(%d %d %d %d) %e | (%e %e) (%e %e)\n" , t , ip , j , k , diff ,
		  creal( a ) , cimag( a ) ,
		  creal( b ) , cimag( b ) ) ;
	}
      }
    }
  }

  free( A ) ;
  free( B ) ;
  free( tT ) ;
#endif
  
  pool_device_free( d_buf0 ) ;
  pool_device_free( d_ret ) ;
  pool_device_free( d_A ) ;
  pool_device_free( d_B ) ;

  return ;
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
  const size_t nEv = 32 ;
  const size_t nmom = 16 ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;

  printf( "nEv %d | nmom %d | LT %d\n" , nEv , nmom , X[3] ) ;
  
  // create a "perambulator" - tau - of noise with these indices
  double _Complex *tau = (double _Complex*)malloc( X[3]*4*4*nEv*nEv*sizeof(double _Complex));
  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt ;
  for( size_t i = 0 ; i < X[3]*4*4*nEv*nEv ; i++ ) {
    const std::complex z( unif(mt) , unif(mt) ) ;
    tau[i] = z.real() + I*z.imag() ;
  }

  // create our "phi" matrix out of noise too
  double _Complex *T = (double _Complex*)calloc( nmom*((size_t)X[3])*nEv*nEv*nEv , sizeof(double _Complex));
  for( size_t i = 0 ; i < nmom*X[3]*nEv*nEv*nEv ; i++ ) {
    const std::complex z( unif(mt) , unif(mt) ) ;
    T[i] = z.real() + I*z.imag() ;
  }
  // enforce some symmetry in T??

  // bad choice of local dot products - performs very shittily
#ifdef SLOW_CPU
  // baryon correlator
  double _Complex C[ X[3] ] ;
  {
    StopWatch dumbCPU ; dumbCPU.start() ;
#pragma omp parallel for
    for( size_t t = 0 ; t < (size_t)X[3] ; t++ ) {
      const std::array<const double _Complex*,3> per = {
	tau + nEv*nEv*16*t ,
	tau + nEv*nEv*16*t ,
	tau + nEv*nEv*16*t } ;
      //C[t] = Tcontract1( T , per[0] , per[1] , per[2] , T + nEv*nEv*nEv*t , nEv ) ;
      //C[t] = Tcontract2( T , per[0] , per[1] , per[2] , T + nEv*nEv*nEv*t , nEv ) ;
      C[t] = Tcontract3( T , per[0] , per[1] , per[2] , T + nEv*nEv*nEv*t , nEv ) ;
      //C[t] = Tcontract4( T , per , T + nEv*nEv*nEv*t , nEv , ijk ) ;
      //printf( "C[%zu] %e %e \n" , t , creal( C[t] ) , cimag( C[t] ) ) ; 
    }
    dumbCPU.stop() ;
    double CPUtime = dumbCPU.getTimeInSeconds() ;
    printLaph( make_strf( "\n CPU contraction time in %g seconds\n" , CPUtime ) ) ;  
  }
#endif

  printf( "Dropping into Tcontract6\n" ) ;
  double _Complex C2[ nmom*nmom*X[3] ] ;
  StopWatch dumbCPU ; dumbCPU.start() ;
  const std::array<const double _Complex*,3> per = {
    tau , tau , tau } ;

  Tcontract6( C2 ,  T , per , X[3] , nEv , nmom , ijk ) ;
  /*
#pragma omp parallel for
  for( size_t t = 0 ; t < (size_t)X[3] ; t++ ) {
    const std::array<const double _Complex*,3> per = {
      tau + nEv*nEv*16*t ,
      tau + nEv*nEv*16*t ,
      tau + nEv*nEv*16*t } ;
    Tcontract6( &C2[t] ,  T , per , T + nEv*nEv*nEv*t , X[3] , nEv , ijk ) ;
  }
  */
  dumbCPU.stop() ;
  double CPUtime = dumbCPU.getTimeInSeconds() ;
  printLaph( make_strf( "\n CPU contraction time in %g seconds\n" , CPUtime ) ) ;  

#ifdef SLOW_CPU
  for( size_t t = 0 ; t < (size_t)X[3] ; t++ ) {
    printf( "DIFF CPU2-CPU :: %zu ( %e , %e ) || ( %e %e ) :: %e \n" , t ,
	    creal( C2[t] ) , cimag( C2[t] ) ,
	    creal( C[t] ) , cimag( C[t] ) ,
	    cabs( C2[t] - C[t] ) ) ; 
  }
#endif
  double _Complex C3[ nmom*nmom*X[3] ] ;
  /*
  const std::array<const double _Complex*,3> per = {
    tau , tau , tau } ;
  */
  TcontractGPU( C3 ,  T , per , X[3] , nEv , nmom , ijk ) ;
  
  memset( C3 , 0. , nmom*nmom*X[3]*sizeof(double _Complex));
  StopWatch GPU ; GPU.start() ;
  TcontractGPU( C3 ,  T , per , X[3] , nEv , nmom , ijk ) ;
  GPU.stop() ;
  double GPUtime = GPU.getTimeInSeconds() ;
  printLaph( make_strf( "\n GPU contraction time in %g seconds\n" , GPUtime ) ) ;

  printLaph( make_strf( "\nGPU speedup factor %gx\n" , CPUtime/GPUtime )  ) ;

  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
      for( size_t t = 0 ; t < (size_t)X[3] ; t++ ) {
	const double _Complex a = C3[t + X[3]*( psnk + nmom*psrc )] ;	
	const double _Complex b = C2[t + X[3]*( psnk + nmom*psrc )] ; 
	printf( "DIFF GPU-CPU :: (%d,%d,%zu) ( %e , %e ) || ( %e %e ) :: %e \n" ,
		psrc , psnk , t ,
		creal( a ) , cimag( a ) ,
		creal( b ) , cimag( b ) ,
		cabs( a - b ) ) ; 
      }
    }
  }
  
  free( tau ) ;
  free( T ) ;
  finalize( ) ;
  return 0;
}
 
