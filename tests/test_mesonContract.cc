/**
   Test the feasibility of Meson contracton on the device

   The idea is to compute

   FWD(t,p)^{\alpha\kappa}_{ik} = \tau(0,p)_ij gl^{\alpha\beta} \phi1(t)^{\beta\kappa}_jk

   BWD(t,p)^{\alpha\kappa}_{ik} = g5.\phi2(t)^{\beta\kappa}_ij.g5. gr^{\alpha\beta}\tau(t,p)_jk

   Such that we can form the inner product 

   C(p1,p2,t) = FWD(t,p)^{\alpha\kappa}_ik (BWD(t,p2)^{\alpha\kappa}_ik)^\dagger

   Ok we need a g5.\phi.g5

   TODO

   Might as well compute disco while we are here, is just a trace of Mfwd and Mbwd
 */
#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"
#include "gammas.h"

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

// Quda interface

//#define VERBOSE_COMPARISON
//#define GPU_STRESS
#define CPU_CROSSCHECK
//#define VERBOSE_CPU

static inline
double _Complex getfac( const Z4 z4 )
{
  switch( z4 ) {
  case p1 : return +1 ;
  case pi : return +I ;
  case m1 : return -1 ;
  case mi : return -I ;
  }
  assert( true ) ;
  return 1 ;
}

// left-multpily a perambulor tau by a gamma matrix
// tau'^{\alpha\kappa} = L^{\alpha\beta}(tau(t)^{\beta\kappa}_{ij})
void
gamL( double _Complex *res ,
      const NRgamma L ,
      const double _Complex *tau ,
      const size_t nEv ,
      const size_t LT )
{
#pragma omp for collapse(3)
  for( size_t t = 0 ; t < LT ; t++ ) {
    for( int alpha = 0 ; alpha < 4 ; alpha++ ) {
      for( int beta = 0 ; beta < 4 ; beta++ ) {
	const int col = L.ig[alpha] ;
	const double _Complex fac = getfac( L.g[alpha]) ;
	#if 1
	memcpy( res + nEv*nEv*( beta + 4*( alpha + 4*t )) ,
		tau + nEv*nEv*( beta + 4*( col + 4*t )) ,
		nEv*nEv*sizeof(double _Complex)) ;
	cblas_zscal( nEv*nEv , &fac , res + nEv*nEv*( beta + 4*( alpha + 4*t )) , 1 ) ;
	#else
	// could this be a copy and a scal? yeah probably
	for( size_t ij = 0 ; ij < nEv*nEv ; ij++ ) {
	  res[ ij + nEv*nEv*( beta + 4*( alpha + 4*t )) ] = \
	    fac*tau[ ij + nEv*nEv*( beta + 4*( col + 4*t )) ] ;
	}
	#endif
      }
    }
  }
}

// left-and right multpily a (transposed in dirac indices) perambulor tau by a gamma matrix
// tau'^{\alpha\kappa} = L^{\alpha\beta}(tau(t)^{\beta\kappa}_{ij})
void
gamLR( double _Complex *res ,
       const NRgamma L ,
       const double _Complex *tau ,
       const NRgamma R ,
       const size_t nEv ,
       const size_t LT )
{
  #pragma omp for collapse(3)
  for( size_t t = 0 ; t < LT ; t++ ) {
    for( int alpha = 0 ; alpha < 4 ; alpha++ ) {
      for( int beta = 0 ; beta < 4 ; beta++ ) {
	const int col1 = L.ig[alpha] , col2 = R.ig[beta] ;
	const double _Complex fac = \
	  getfac( L.g[alpha] )*getfac( R.g[col2] ) ;
	// sneaky transpose in here
	memcpy( res + nEv*nEv*( beta + 4*( alpha + 4*t )) ,
		tau + nEv*nEv*( col1 + 4*( col2 + 4*t )) ,
		nEv*nEv*sizeof(double _Complex)) ;
	cblas_zscal( nEv*nEv , &fac , res + nEv*nEv*( beta + 4*( alpha + 4*t )) , 1 ) ;
      }
    }
  }
}

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

// does A_{ij} = (B_{ik}C_{jk}^\dagger 
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
	      nEv, nEv, nEv , &alpha , B , nEv , C , nEv , &beta , A , nEv ) ;
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

// back to basics....
// Implement the dumb way Eq 11 of https://arxiv.org/pdf/0905.2160 
void
basic_contraction( double _Complex *C ,
		   const double _Complex *phi ,
		   const std::array<const double _Complex*,2> &per,
		   const std::vector<std::array< NRgamma, 2>> &GlGr ,
		   const size_t nEv ,
		   const size_t nmom ,
		   const size_t tsrc ,
		   const size_t LT )
{
  // \phi(t')_ij \gsnk \tau_jk(t',t)
  double _Complex *Msnk = (double _Complex*)malloc( 16*nmom*LT*nEv*nEv*sizeof( double _Complex ) ) ;
  double _Complex *buf1 = (double _Complex*)malloc( 16*LT*nEv*nEv*sizeof( double _Complex) ) ;
  double _Complex *buf2 = (double _Complex*)malloc( 16*LT*nEv*nEv*sizeof( double _Complex) ) ;
  double _Complex *Msrc = (double _Complex*)malloc( 16*nmom*LT*nEv*nEv*sizeof( double _Complex ) ) ;
  NRgamma g5( Gamma_5 ) ;
  
  for( size_t ng = 0 ; ng < GlGr.size() ; ng++ ) {

    // left and right gammas in a parallel region to avoid opening and closing one
#pragma omp parallel
    {
      gamL( buf1 , GlGr[0][0] , per[0] , nEv , LT ) ;

      gamLR( buf2 , GlGr[0][1]*g5 , per[1] , g5 , nEv , LT ) ;
      
      #pragma omp for collapse(3)
      for( size_t p = 0 ; p < nmom ; p++ ) {
	for( size_t t = 0 ; t < LT ; t++ ) {
	  for( size_t ab = 0 ; ab < 16 ; ab++ ) {
	    // Does Msnk^{\alpha\beta}_{ij}(p,t) = phi_{ik}(p,t) Gl^{\alpha\kappa}\tau^{\kappa\beta}_{kj}(t)
	    const double _Complex *PhiP = phi + nEv*nEv*( p + nmom*t ) ;
	    const double _Complex *Tp = buf1 + nEv*nEv*( ab + 16*t ) ;
	    double _Complex *Mp = Msnk + nEv*nEv*( ab + 16*( p + nmom*t ) ) ;
	    innerevprod( PhiP , Tp , nEv , Mp ) ;
	    // Does Msrc^{\alpha\beta}_{ij}(p,t) = ( phi_{ik}(p,t) Gl^{\alpha\kappa}(\gamma_5 \tau^* \gamma_5 )^{\kappa\beta}_{jk}(t) )
	    PhiP = phi + nEv*nEv*( p + nmom*tsrc ) ;
	    Tp = buf2 + nEv*nEv*( ab + 16*t ) ;
	    Mp = Msrc + nEv*nEv*( ab + 16*( p + nmom*t ) ) ;
	    innerevprod_dag( PhiP , Tp , nEv , Mp ) ;
	  }
	}
      }
      
      // dot product Msnk^{\alpha\beta}_{ij} Msrc^{\beta\alpha}_ji
      #pragma omp for collapse(3)
      for( int psrc = 0 ; psrc < nmom ; psrc++ ) {
	for( int psnk = 0 ; psnk < nmom ; psnk++ ) {
	  for( int t = 0 ; t < LT ; t++ ) {
	    C[ psnk + nmom*( psrc + nmom*( t + LT*ng )) ] =	\
	      cblas_zdotu( nEv*nEv*16 ,
			   Msnk + nEv*nEv*16*( psnk + nmom*t ) , 1 ,
			   Msrc + nEv*nEv*16*( psrc + nmom*t ) , 1 ) ;
	  }
	}
      }
    } // close parallel region
  }

  free( Msnk ) ; free( Msrc ) ; free( buf1 ) ; free( buf2 ) ;
}

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
	      std::vector<std::array< NRgamma,2>> &GlGr ,
	      const size_t LT ,
	      const size_t nEv ,
	      const size_t nmom ,
	      const int precision = QUDA_DOUBLE_PRECISION )
{
  const NRgamma g5( Gamma_5 ) ;

  // host gamma'd versions of the perambulators
  double _Complex *gtmp1 = (double _Complex*)malloc( LT*4*4*nEv*nEv*sizeof( double _Complex ) ) ;
  double _Complex *gtmp2 = (double _Complex*)malloc( LT*4*4*nEv*nEv*sizeof( double _Complex ) ) ;

  const size_t mbytes = LT*nmom*nEv*nEv*16*2*precision ;
  const size_t rbytes = LT*nmom*nmom*2*precision ;
  const size_t tbytes = LT*nmom*nEv*nEv*2*precision ;
  if( getVerbosity() >= QUDA_SUMMARIZE ) {
    const size_t OneGB = 1024*1024*1024;
    const size_t total_bytes = 2*mbytes + rbytes + 2*tbytes ;
    printfQuda("d_ret %fGB | dT %fGB | dT2 %fGB| d_fwd %fGB | d_bwd %fGB | total = %fGB\n",
	       (double)rbytes/(OneGB) , (double)tbytes/(OneGB), (double)tbytes/OneGB ,
	       (double)mbytes/(OneGB), (double)mbytes/(OneGB) , (double)total_bytes/(OneGB)); 
  }
  void *d_fwd = pool_device_malloc( mbytes ) ;
  void *d_bwd = pool_device_malloc( mbytes ) ;
  void *d_ret = pool_device_malloc( rbytes ) ;
  void *d_T  = pool_device_malloc( tbytes ) ;
  void *d_T2 = pool_device_malloc( tbytes ) ;
  
  // copy "T" to device
  H2Dwrap( d_T , T , LT*nmom*nEv*nEv , precision ) ;

  const int tsrc = 0 ;
  for( size_t t = 0 ; t < LT ; t++ ) {
    qudaMemcpy( (char*)d_T2 + 2*precision*nEv*nEv*nmom*t,
		(char*)d_T  + 2*precision*nEv*nEv*nmom*tsrc ,
		2*precision*nEv*nEv*nmom , qudaMemcpyDeviceToDevice);
  }
  
  // loop gamma src and sink combinations
  for( size_t g = 0 ; g < GlGr.size() ; g++ ) {
  
    // fwd is T_{ij}(0,p) S_{jk}^{\alpha\beta}(t)
#pragma omp parallel
    {
      gamL( gtmp1 , GlGr[g][0] , S[0] , nEv , LT ) ;
      gamLR( gtmp2 , GlGr[g][1]*g5 , S[1] , g5 , nEv , LT ) ;
    }
    for( size_t t = 0 ; t < LT ; t++ ) {
      H2Dwrap( (char*)d_fwd + 2*precision*nEv*nEv*16*nmom*t ,
	       gtmp1 + nEv*nEv*16*t ,
	       nEv*nEv*16, precision) ;
      for( size_t p = 1 ; p < nmom ; p++ ) {
      	qudaMemcpy( (char*)d_fwd + 2*precision*16*nEv*nEv*( p + nmom*t ) ,
		    (char*)d_fwd + 2*precision*16*nEv*nEv*nmom*t ,
		    2*precision*16*nEv*nEv , qudaMemcpyDeviceToDevice);
      }
    }

    // BLAS function1 creating fwd elemental - CPU GPU agree
    // T_{p,t,i,j}S_{t,alpha,beta,j,k} = M_{p,t,alpha,beta,i,k}
    for( size_t ab = 0 ; ab < 16 ; ab++ ) {
      QudaBLASParam cublas_param1 = newQudaBLASParam();
      cublas_param1.trans_a = QUDA_BLAS_OP_N ;
      cublas_param1.trans_b = QUDA_BLAS_OP_N ;
      cublas_param1.m = nEv ;
      cublas_param1.n = nEv ;
      cublas_param1.k = nEv ;
      cublas_param1.lda = nEv ;
      cublas_param1.ldb = nEv ;
      cublas_param1.ldc = nEv ;
      cublas_param1.a_stride = nEv*nEv ;
      cublas_param1.b_stride = nEv*nEv*16 ;
      cublas_param1.c_stride = nEv*nEv*16 ; // tslowest
      cublas_param1.batch_count = LT*nmom ; // do "T" batches
      cublas_param1.alpha = 1.0 ; cublas_param1.beta = 0.0 ;
      cublas_param1.data_order = QUDA_BLAS_DATAORDER_ROW;
      cublas_param1.data_type = ( precision == QUDA_SINGLE_PRECISION ) ? \
	QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
      blas_lapack::native::stridedBatchGEMM( (char*)d_T ,
					     (char*)d_fwd + 2*precision*nEv*nEv*ab ,
					     (char*)d_fwd + 2*precision*nEv*nEv*ab ,
					     cublas_param1, QUDA_CUDA_FIELD_LOCATION);
    }
    
    // bwd
    const int tsrc = 0 ;
    for( size_t t = 0 ; t < LT ; t++ ) {
      H2Dwrap( (char*)d_bwd + nEv*nEv*16*t*nmom*2*precision ,
	       gtmp2 + nEv*nEv*16*t ,
	       nEv*nEv*16, precision) ;
      for( size_t p = 1 ; p < nmom ; p++ ) {
      	qudaMemcpy( (char*)d_bwd + 2*precision*16*nEv*nEv*( p + nmom*t ) ,
		    (char*)d_bwd + 2*precision*16*nEv*nEv*nmom*t ,
		    2*precision*16*nEv*nEv , qudaMemcpyDeviceToDevice);
      }
    }

    // BLAS function2 creating second elemental
    // T_{p,0,i,j}S1_{t,alpha,beta,j,k}^\dagger = M_{p,t,alpha,beta,i,k}
    for( size_t ab = 0 ; ab < 16 ; ab++ ) {
      QudaBLASParam cublas_param1 = newQudaBLASParam();
      cublas_param1.trans_a = QUDA_BLAS_OP_N ;
      cublas_param1.trans_b = QUDA_BLAS_OP_C ;
      cublas_param1.m = nEv ;
      cublas_param1.n = nEv ;
      cublas_param1.k = nEv ;
      cublas_param1.lda = nEv ;
      cublas_param1.ldb = nEv ;
      cublas_param1.ldc = nEv ;
      cublas_param1.a_stride = nEv*nEv ;
      cublas_param1.b_stride = nEv*nEv*16 ;
      cublas_param1.c_stride = nEv*nEv*16 ; // tslowest
      cublas_param1.batch_count = LT*nmom ; // do "T" batches for now
      cublas_param1.alpha = 1.0 ; cublas_param1.beta = 0.0 ;
      cublas_param1.data_order = QUDA_BLAS_DATAORDER_ROW;
      cublas_param1.data_type = ( precision == QUDA_SINGLE_PRECISION ) ? \
	QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
      blas_lapack::native::stridedBatchGEMM( (char*)d_T2 ,
					     (char*)d_bwd + 2*precision*nEv*nEv*ab ,
					     (char*)d_bwd + 2*precision*nEv*nEv*ab ,
					     cublas_param1, QUDA_CUDA_FIELD_LOCATION);
    }
    
    // innerproduct trace over all psrc psnk
    // fwd_{t,psrc,alpha,beta,i,j} bwd_{t,psnk,beta,alpha,j,i}
    QudaBLASParam cublas_param3 = newQudaBLASParam();
    cublas_param3.trans_a = QUDA_BLAS_OP_N ;
    cublas_param3.trans_b = QUDA_BLAS_OP_T ;
    cublas_param3.m = nmom ;
    cublas_param3.n = nmom ;
    cublas_param3.k = nEv*nEv*4*4 ;
    cublas_param3.lda = nEv*nEv*4*4 ;
    cublas_param3.ldb = nEv*nEv*4*4 ;
    cublas_param3.ldc = nmom ;
    cublas_param3.a_stride = nmom*4*4*nEv*nEv ;
    cublas_param3.b_stride = nmom*4*4*nEv*nEv ;
    cublas_param3.c_stride = nmom*nmom ;
    cublas_param3.batch_count = LT ; // do "T" batches
    cublas_param3.alpha = 1.0 ; cublas_param3.beta = 0.0 ;
    cublas_param3.data_order = QUDA_BLAS_DATAORDER_ROW;
    cublas_param3.data_type = ( precision == QUDA_SINGLE_PRECISION ) ?	\
      QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
    blas_lapack::native::stridedBatchGEMM( d_bwd, d_fwd, d_ret, cublas_param3,
					   QUDA_CUDA_FIELD_LOCATION);
    // copy back to "sum" on the host
    D2Hwrap( sum + g*nmom*nmom*LT , d_ret , nmom*nmom*LT , precision ) ;
  }
  
  free( gtmp1 ) ; free( gtmp2 ) ;
  
  pool_device_free( d_ret ) ;
  pool_device_free( d_fwd ) ;
  pool_device_free( d_bwd ) ;
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
  const size_t nEv = 96 ;
  const size_t nmom = 42 ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  
  // create a "perambulator" - tau - of noise with these indices
  double _Complex *tau = (double _Complex*)malloc( X[3]*4*4*nEv*nEv*sizeof(double _Complex));
  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt(12345) ;
  for( size_t i = 0 ; i < X[3]*4*4*nEv*nEv ; i++ ) {
    const std::complex z( unif(mt) , unif(mt) ) ;
    tau[i] = z.real() + I*z.imag() ;
  }

  // create our "phi" matrix out of noise too
  double _Complex *phi = (double _Complex*)malloc( X[3]*nmom*nEv*nEv*sizeof(double _Complex));
  for( size_t i = 0 ; i < X[3]*nmom*nEv*nEv ; i++ ) {
    const std::complex z( unif(mt) , unif(mt) ) ;
    phi[i] = z.real() + I*z.imag() ;
  }

  const std::array<const double _Complex*,2> per = { tau , tau  } ;

  Gbasis g( false ) ;
  std::vector<std::array< NRgamma, 2>> GlGr = { { g.G[Gamma_5] , g.G[Gamma_5].dagger() } ,
						{ g.G[Gamma_X] , g.G[Gamma_X].dagger() } ,
						{ g.G[Gamma_Y] , g.G[Gamma_Y].dagger() } ,
						{ g.G[Gamma_Z] , g.G[Gamma_Z].dagger() } } ;
    

  double _Complex *C1 = (double _Complex*)calloc( GlGr.size()*nmom*nmom*X[3] , sizeof( double _Complex) ) ;

  StopWatch CPU ; CPU.start() ;
  basic_contraction( C1 , phi , per , GlGr , nEv , nmom , 0 , X[3] ) ;
  CPU.stop() ;
  const double CPUtime = CPU.getTimeInSeconds() ;
  printLaph( make_strf( "\n CPU in %g seconds\n" , CPUtime ) ) ;

  // GPU
  double _Complex *C2 = (double _Complex*)calloc( GlGr.size()*nmom*nmom*X[3] , sizeof( double _Complex) ) ;
  McontractGPU( C2 , phi , per , GlGr , X[3] , nEv , nmom , QUDA_SINGLE_PRECISION ) ;

  StopWatch GPU ; GPU.start() ;
  McontractGPU( C2 , phi , per , GlGr , X[3] , nEv , nmom , QUDA_SINGLE_PRECISION ) ;
  GPU.stop() ;
  double GPUtime = GPU.getTimeInSeconds() ;
  printLaph( make_strf( "\n GPU in %g seconds\n" , GPUtime ) ) ;

  printLaph( make_strf( "GPU supremacy %gx\n" , CPUtime/GPUtime ) ) ;

#ifdef VERBOSE_CPU
  for( size_t ng = 0 ; ng < GlGr.size() ; ng++ ) {
    GlGr[ng][0].print() ; GlGr[ng][1].print() ;
    for( size_t psrc = 0 ; psrc < nmom ; psrc++ ) {
      for( size_t psnk = 0 ; psnk < nmom ; psnk++ ) {
	for( size_t t = 0 ; t < X[3] ; t++ ) {
	  const double _Complex z1 = C1[ psnk + nmom*( psrc + nmom*(t + X[3]*ng) ) ] ;
	  const double _Complex z2 = C2[ psnk + nmom*( psrc + nmom*(t + X[3]*ng) ) ] ;
	  printf( "Psrc,Psnk,t : (%zu %zu %zu) [%e,%e] [%e,%e] :: %e \n" , psrc , psnk , t ,
		  creal( z1 ) , cimag( z1 ) ,
		  creal( z2 ) , cimag( z2 ) ,
		  cabs( z1 - z2 )/cabs(z1) ) ;
	}
      }
    }
  }
#endif

  free( phi ) ; free( tau ) ;
  free( C2 ) ; free( C1 ) ;
  
  finalize( ) ;
  return 0;
}
