/**
   Test the feasibility of Baryon contractions on the device

   Assumes that source position is at T = 0 and batches in t

   Requires T of T_{t,p,i,j,k} order of indices
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

//#define GPU_STRESS
#define CPU_CROSSCHECK
//#define VERBOSE_COMPARISON
//#define SLOW_CPU

typedef enum bartype { ijk = 0 } ;

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

// CPU version of GPU code
static void
Tcontract7( double _Complex *sum ,
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
			   T + k + nEv*j + nEv*nEv*nEv*(psrc) , nEv*nEv ,
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
	sum[ psnk + nmom*(psrc + nmom*t) ] =			\
	  cblas_zdotu( nEv*nEv*nEv , A + nEv*nEv*nEv*psrc , 1 ,
		       T + nEv*nEv*nEv*(psnk+nmom*t) , 1 ) ;
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
// does TL_{ijk}(T) S1_{ii'}(t) S2_{jj'}(t) S3_{kk'}(t) TR_{i'j'k'}(t)
static void
TcontractGPU3( double _Complex *sum , 
	       const double _Complex *T ,
	       const std::array<const double _Complex*,3> &S ,
	       const size_t LT ,
	       const size_t nEv ,
	       const size_t nmom ,
	       const bartype btype ,
	       const std::vector< std::array<const std::array<int,2> , 3>> &sidx ,
	       const int precision = QUDA_DOUBLE_PRECISION )
{
  // device "T" object
  const size_t sbytes = LT*nEv*nEv*2*precision ;
  void *d_buf0 = pool_device_malloc( sbytes ) ;

  // device buffers for perambulators
  const size_t tbytes = LT*nmom*nEv*nEv*nEv*2*precision ;
  void *d_T = pool_device_malloc( tbytes ) ;
  void *d_A = pool_device_malloc( tbytes ) ;
  void *d_B = pool_device_malloc( tbytes ) ;

  const size_t retbytes = LT*nmom*nmom*2*precision ;
  void *d_ret = pool_device_malloc( retbytes ) ;

  const double OneGB = 1024.*1024.*1024. ;
  printfQuda( "Trip %g GB | buf0 %g GB | ret %g GB | A %g GB | B %g GB \n" ,
	      tbytes/OneGB , sbytes/OneGB , retbytes/OneGB , tbytes/OneGB , tbytes/OneGB ) ;

  StopWatch Mem ; Mem.start() ;
 
  // quda memcpy the whole baryon kernel to the device once only as it is expensive
  H2Dwrap(d_T, T, LT*nmom*nEv*nEv*nEv, precision) ;

  Mem.stop() ;
  double Memtime = Mem.getTimeInSeconds() ;
  printLaph( make_strf( "\n Mem in %g seconds\n" , Memtime ) ) ;

  printfQuda( "Picking %zu spin indices\n" , sidx.size() ) ;
  
  for( size_t nspin = 0 ; nspin < sidx.size() ; nspin++ ) {

    const std::array<const std::array<int,2> , 3> spin_idx = sidx[nspin] ;

    printfQuda( "Contracting tau_0^{%d%d} tau_1^{%d%d} tau_2^{%d%d}\n" ,
		spin_idx[0][0] , spin_idx[0][1] ,
		spin_idx[1][0] , spin_idx[1][1] ,
		spin_idx[2][0] , spin_idx[2][1] ) ;
    
    for( size_t t = 0 ; t < LT ; t++ ) {
      H2Dwrap( (char*)d_buf0 + nEv*nEv*t*2*precision ,
	       S[0] + nEv*nEv*( spin_idx[0][1] + 4*(spin_idx[0][0] + 4*t ) ) ,
	       nEv*nEv, precision) ;
    }
      
    StopWatch MMuls ; MMuls.start() ;
    QudaBLASParam cublas_param1 = newQudaBLASParam();
    cublas_param1.trans_a = QUDA_BLAS_OP_T ;
    cublas_param1.trans_b = QUDA_BLAS_OP_N ;
    cublas_param1.m = nEv ;
    cublas_param1.n = nEv*nEv ;
    cublas_param1.k = nEv ;
    cublas_param1.lda = nEv;
    cublas_param1.ldb = nEv*nEv;
    cublas_param1.ldc = nEv*nEv;
    cublas_param1.a_stride = nEv*nEv ;
    cublas_param1.b_stride = 0 ;
    cublas_param1.c_stride = nEv*nEv*nEv*nmom ;
    cublas_param1.batch_count = LT ; // do "T" batches
    cublas_param1.alpha = 1.0 ; cublas_param1.beta = 0.0 ;
    cublas_param1.data_order = QUDA_BLAS_DATAORDER_ROW;
    cublas_param1.data_type = ( precision == QUDA_SINGLE_PRECISION ) ?	\
      QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
    for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {      
      // checked and is ok
      blas_lapack::native::stridedBatchGEMM( d_buf0,
					     (char*)d_T + 2*precision*nEv*nEv*nEv*psrc ,
					     (char*)d_A + 2*precision*nEv*nEv*nEv*psrc ,
					     cublas_param1, QUDA_CUDA_FIELD_LOCATION);
    }
    
    // second part break up into t i' batches
    for( size_t t = 0 ; t < LT ; t++ ) {
      H2Dwrap( (char*)d_buf0 + nEv*nEv*t*2*precision ,
	       S[1] + nEv*nEv*( spin_idx[1][1] + 4*(spin_idx[1][0] + 4*t ) ) ,
	       nEv*nEv, precision) ;
    }
    
    // D_{jj'}^T(t) A_{i'jk}(t,p) -> B_{i'j'k}(t,p)
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
    cublas_param2.data_type = ( precision == QUDA_SINGLE_PRECISION ) ?	\
      QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
    for( size_t t = 0 ; t < LT ; t++ ) {      
      for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
	blas_lapack::native::stridedBatchGEMM( (char*)d_buf0 + nEv*nEv*t*2*precision ,
					       (char*)d_A + 2*precision*nEv*nEv*nEv*(psrc + nmom*t),
					       (char*)d_B + 2*precision*nEv*nEv*nEv*(psrc + nmom*t),
					       cublas_param2, QUDA_CUDA_FIELD_LOCATION);
      }
    }
    
    // third product
    for( size_t t = 0 ; t < LT ; t++ ) {
      H2Dwrap( (char*)d_buf0 + nEv*nEv*t*2*precision ,
	       S[2] + nEv*nEv*( spin_idx[2][1] + 4*(spin_idx[2][0] + 4*t ) ) ,
	       nEv*nEv, precision) ;
    }

    // A_{i'jk'}(t,p) = B_{i'j'k}(t,p)D_{kk'}(t)
    QudaBLASParam cublas_param3 = newQudaBLASParam();
    cublas_param3.trans_a = QUDA_BLAS_OP_N ;
    cublas_param3.trans_b = QUDA_BLAS_OP_N ;
    cublas_param3.m = nEv*nEv ;
    cublas_param3.n = nEv ;
    cublas_param3.k = nEv ;
    cublas_param3.lda = nEv;
    cublas_param3.ldb = nEv;
    cublas_param3.ldc = nEv;
    cublas_param3.a_stride = nEv*nEv*nEv*nmom ;
    cublas_param3.b_stride = nEv*nEv ;
    cublas_param3.c_stride = nEv*nEv*nEv*nmom ;
    cublas_param3.batch_count = LT ; // do "T" batches
    cublas_param3.alpha = 1.0 ; cublas_param3.beta = 0.0 ;
    cublas_param3.data_order = QUDA_BLAS_DATAORDER_ROW;
    cublas_param3.data_type = ( precision == QUDA_SINGLE_PRECISION ) ?	\
      QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
    for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {     
      blas_lapack::native::stridedBatchGEMM( (char*)d_B + nEv*nEv*nEv*psrc*2*precision ,
					     d_buf0 ,
					     (char*)d_A + nEv*nEv*nEv*psrc*2*precision ,
					     cublas_param3, QUDA_CUDA_FIELD_LOCATION);
    }
    
    MMuls.stop() ;
    double MMulstime = MMuls.getTimeInSeconds() ;
    printLaph( make_strf( "\n MMULs in %g seconds\n" , MMulstime ) ) ;
    
    StopWatch Inner ; Inner.start() ;

    // pinner
    QudaBLASParam cublas_param4 = newQudaBLASParam();
    cublas_param4.trans_a = QUDA_BLAS_OP_N ;
    cublas_param4.trans_b = QUDA_BLAS_OP_T ;
    cublas_param4.m = nmom ;
    cublas_param4.n = nmom ;
    cublas_param4.k = nEv*nEv*nEv ;
    cublas_param4.lda = nEv*nEv*nEv;
    cublas_param4.ldb = nEv*nEv*nEv;
    cublas_param4.ldc = nmom;
    cublas_param4.a_stride = nmom*nEv*nEv*nEv ;
    cublas_param4.b_stride = nmom*nEv*nEv*nEv ;
    cublas_param4.c_stride = nmom*nmom ;
    cublas_param4.batch_count = LT ;
    cublas_param4.alpha = 1.0 ; cublas_param4.beta = 0.0 ;
    cublas_param4.data_order = QUDA_BLAS_DATAORDER_ROW;
    cublas_param4.data_type = ( precision == QUDA_SINGLE_PRECISION ) ?	\
      QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
    blas_lapack::native::stridedBatchGEMM( d_A, d_T, d_ret,
					   cublas_param4, QUDA_CUDA_FIELD_LOCATION);
    
    Inner.stop() ;
    double Innertime = Inner.getTimeInSeconds() ;
    printLaph( make_strf( "\n Inner in %g seconds\n" , Innertime ) ) ;
    
    // copy back to host
    D2Hwrap( sum + nmom*nmom*LT*nspin , d_ret , nmom*nmom*LT , precision ) ;
  }
  
  pool_device_free( d_buf0 ) ;
  pool_device_free( d_ret ) ;
  pool_device_free( d_A ) ;
  pool_device_free( d_B ) ;
  pool_device_free( d_T ) ;

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
  const size_t nEv = 42 ;
  const size_t nmom = 24 ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;

  printf( "nEv %zu | nmom %zu | LT %d\n" , nEv , nmom , X[3] ) ;
  
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
  double _Complex C[ nmom*nmom*X[3] ] ;
  {
    StopWatch dumbCPU ; dumbCPU.start() ;
#pragma omp parallel for
    for( size_t t = 0 ; t < (size_t)X[3] ; t++ ) {
      const std::array<const double _Complex*,3> per = {
	tau + nEv*nEv*16*t ,
	tau + nEv*nEv*16*t ,
	tau + nEv*nEv*16*t } ;
      for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
	for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
	  C[psnk+nmom*(psrc+nmom*t)] = Tcontract1( T + nEv*nEv*nEv*psrc ,
						   per[0] , per[1] , per[2] ,
						   T + nEv*nEv*nEv*( psnk + nmom*t ) , nEv ) ;
	}
      }
    }
    dumbCPU.stop() ;
    double CPUtime = dumbCPU.getTimeInSeconds() ;
    printLaph( make_strf( "\n CPU contraction time in %g seconds\n" , CPUtime ) ) ;  
  }
#endif

  const std::array<const double _Complex*,3> per = {
    tau , tau , tau } ;
  
#ifdef CPU_CROSSCHECK
  printf( "Dropping into Tcontract6\n" ) ;
  double _Complex C2[ nmom*nmom*X[3] ] ;
  StopWatch dumbCPU ; dumbCPU.start() ;
  Tcontract7( C2 ,  T , per , X[3] , nEv , nmom , ijk ) ;

  dumbCPU.stop() ;
  double CPUtime = dumbCPU.getTimeInSeconds() ;
  printLaph( make_strf( "\n CPU contraction time in %g seconds\n" , CPUtime ) ) ;  

#ifdef SLOW_CPU
  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
      
      for( size_t t = 0 ; t < (size_t)X[3] ; t++ ) {
	const double _Complex a = C[ psnk + nmom*( psrc + nmom*t) ] ;	
	const double _Complex b = C2[ psnk + nmom*( psrc + nmom*t) ] ;
	printf( "DIFF CPU2-CPU :: (%zu,%zu,%zu) ( %e , %e ) || ( %e %e ) :: %e \n" ,
		psrc, psnk , t ,
		creal( a ) , cimag( a ) ,
		creal( b ) , cimag( b ) ,
		cabs( a-b ) ) ; 
      }
    }
  }
#endif
#endif

  // so many parentheses it's like I'm programming in LISP
  const std::vector< std::array<const std::array<int,2> , 3>> sidx = { {{{0,0},{0,0},{0,0}}} ,
								       {{{1,1},{1,1},{1,1}}} ,
								       {{{2,2},{2,2},{2,2}}} ,
								       {{{3,3},{3,3},{3,3}}} } ;
  double _Complex C3[ nmom*nmom*X[3]*sidx.size() ] ;
  TcontractGPU3( C3 ,  T , per , X[3] , nEv , nmom , ijk , sidx , QUDA_SINGLE_PRECISION ) ;
  
  memset( C3 , 0. , sidx.size()*nmom*nmom*X[3]*sizeof(double _Complex));
  StopWatch GPU ; GPU.start() ;
  TcontractGPU3( C3 ,  T , per , X[3] , nEv , nmom , ijk , sidx , QUDA_SINGLE_PRECISION ) ;
  GPU.stop() ;
  double GPUtime = GPU.getTimeInSeconds() ;
  printLaph( make_strf( "\n GPU contraction time in %g seconds\n" , GPUtime ) ) ;

#ifdef CPU_CROSSCHECK
  printLaph( make_strf( "\nGPU speedup factor %gx\n" , sidx.size()*CPUtime/GPUtime )  ) ;

  #ifdef VERBOSE_COMPARISON
  for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
    for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
      for( size_t t = 0 ; t < (size_t)X[3] ; t++ ) {
	const double nrm = cabs( C2[ psnk + nmom*( psrc + nmom*t) ] ) ;
	const double _Complex a = C3[ psnk + nmom*( psrc + nmom*t) ] ;	
	// tfastest
	const double _Complex b = C2[ psnk + nmom*( psrc + nmom*t) ] ; //C2[t + X[3]*( psnk + nmom*psrc )] ; 
	printf( "DIFF GPU-CPU :: (%zu,%zu,%zu) ( %e , %e ) || ( %e %e ) :: %e \n" ,
		psrc , psnk , t ,
		creal( a ) , cimag( a ) ,
		creal( b ) , cimag( b ) ,
		cabs( a - b )/nrm ) ; 
      }
    }
  }
  #else
  for( size_t t = 0 ; t < (size_t)X[3] ; t++ ) {
    double sum = 0. ;
    for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
      for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
	const double nrm = cabs( C2[ psnk + nmom*( psrc + nmom*t) ] ) ;
	const double _Complex a = C3[ psnk + nmom*( psrc + nmom*t) ] ;	
	const double _Complex b = C2[ psnk + nmom*( psrc + nmom*t) ] ;
	sum += cabs( a - b )/nrm ;
      }
    }
    printf( "Relative dev CPU-GPU %zu %e \n" , t , sum/(nmom*nmom) ) ;
  }
  #endif
#endif
  
  free( tau ) ;
  free( T ) ;
  finalize( ) ;
  return 0;
}
 
