#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

#include <quda.h>
#include <quda_internal.h>
#include <timer.h>
#include <blas_lapack.h>
#include <blas_quda.h>
#include <tune_quda.h>
#include <color_spinor_field.h>
#include <contract_quda.h>

#include <complex>

#include <random>
#include <cassert>
#include <complex.h>

using namespace LaphEnv ;
using namespace quda ;

//#define VERBOSE_COMPARISON
//#define GPU_STRESS
#define CPUCROSSCHECK
//#define CPU_STRESS

static const double OneGB = 1024.*1024.*1024.;

static inline void
evprod( const double _Complex *coeffs ,
	const void *const *host_evec ,
	const size_t n ,
	const size_t nEv ,
	const size_t nsites ,
	double _Complex *q )
{
#pragma omp for
  for( size_t i = 0 ; i < nsites ; i++ ) {
    const double _Complex *pt2[nEv] ;
    for( size_t ev = 0 ; ev < nEv ; ev++ ) {
      pt2[ev] = (double _Complex*)host_evec[ev]+3*i ; 
    }
    for( size_t dil = 0 ; dil < n ; dil++ ) {      
      for( size_t c = 0 ; c < 3 ; c++ ) {
	double _Complex sum = 0.0 ;
	for( size_t ev = 0 ; ev < nEv ; ev++ ) {	  
	  sum += (*(pt2[ev]+c) ) * coeffs[ev+dil*nEv] ;
	}
	q[c+3*(i+nsites*dil)] = sum ;
      }
    }
  }
}

static inline void
evprodv2( const std::vector<const double _Complex *>coeffs ,
	  const void *const *host_evec ,
	  const std::vector<size_t> ndil ,
	  const size_t nEv ,
	  const size_t nsites ,
	  std::vector<double _Complex *>q )
{
#pragma omp for
  for( size_t i = 0 ; i < nsites ; i++ ) {
    const double _Complex alpha = 1. , beta = 0. ;
    double _Complex pt2[3][nEv] ;
    for( size_t ev = 0 ; ev < nEv ; ev++ ) {
      double _Complex *pt = (double _Complex*)host_evec[ev]+3*i ;
      pt2[0][ev] = *(pt+0) ;
      pt2[1][ev] = *(pt+1) ;
      pt2[2][ev] = *(pt+2) ;
    }
    for( size_t nq = 0 ; nq < q.size() ; nq++ ) {
      for( size_t dil = 0 ; dil < ndil[nq] ; dil++ ) {
	double _Complex *pc = (double _Complex*)coeffs[nq]+dil*nEv ;
	// unrolled color loop should be a matvec
        #ifdef USE_OPENBLAS
	// unrolled inner products turn out to be faster 
	//cblas_zgemv( CblasRowMajor , CblasNoTrans , 3 , nEv , &alpha , &pt2[0][0] , nEv , pc , 1 , &beta , q+3*(i+nsites*dil) , 1 ) ;
	q[nq][0+3*(i+nsites*dil)] = cblas_zdotu( nEv , pt2[0] , 1 , pc , 1 ) ;
	q[nq][1+3*(i+nsites*dil)] = cblas_zdotu( nEv , pt2[1] , 1 , pc , 1 ) ;
	q[nq][2+3*(i+nsites*dil)] = cblas_zdotu( nEv , pt2[2] , 1 , pc , 1 ) ;
        #elif (defined USE_GSL)
	cblas_zdotu( nEv , pt2[0] , 1 , pc , 1 , q[nq]+0+3*(i+nsites*dil)  ) ;
	cblas_zdotu( nEv , pt2[1] , 1 , pc , 1 , q[n1]+1+3*(i+nsites*dil) ) ;
	cblas_zdotu( nEv , pt2[2] , 1 , pc , 1 , q[nq]+2+3*(i+nsites*dil) ) ;
        #else
	for( size_t c = 0 ; c < 3 ; c++ ) {
	  double _Complex sum = 0.0 ;
	  for( size_t ev = 0 ; ev < nEv ; ev++ ) {	  
	    sum += pt2[c][ev]*pc[ev] ;
	  }
	  q[nq][c+3*(i+nsites*dil)] = sum ;
	}
        #endif
      }
    }
  }
}

// cpu color cross
static void
cpuColorCross( void *A , void *B , void *result , const int X[4] )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ;
  std::complex<double> *ptA = (std::complex<double>*)A ;
  std::complex<double> *ptB = (std::complex<double>*)B ;
  std::complex<double> *ptC = (std::complex<double>*)result ;
  for( size_t i = 0 ; i < (size_t)Nsites ; i++ ) {
    ptC[ 3*i + 0 ] =  ptA[ 3*i + 1 ]*ptB[ 3*i + 2 ] - ptA[ 3*i + 2 ]*ptB[ 3*i + 1 ] ;
    ptC[ 3*i + 1 ] = -ptA[ 3*i + 0 ]*ptB[ 3*i + 2 ] + ptA[ 3*i + 2 ]*ptB[ 3*i + 0 ] ;
    ptC[ 3*i + 2 ] =  ptA[ 3*i + 0 ]*ptB[ 3*i + 1 ] - ptA[ 3*i + 1 ]*ptB[ 3*i + 0 ] ;
  }
}

// cpu color contract
static void
cpuColorContract( void *A , void *B , void *result , const int X[4] )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ;
  const std::complex<double> *ptA = (const std::complex<double>*)A ;
  const std::complex<double> *ptB = (const std::complex<double>*)B ;
  std::complex<double> *ptC = (std::complex<double>*)result ;
  for( size_t i = 0 ; i < (size_t)Nsites ; i++ ) {
    #ifdef USE_OPENBLAS
    ptC[i] = cblas_zdotu( 3 , ptA+3*i , 1 , ptB+3*i , 1 ) ;
    #elif (defined USE_GSL_CBLAS)
    cblas_zdotu( 3 , ptA+3*i , 1 , ptB+3*i , 1 , ptC+i ) ;
    #else
    ptC[i]  = ptA[0+3*i]*(ptB[0+3*i]) ;
    ptC[i] += ptA[1+3*i]*(ptB[1+3*i]) ;
    ptC[i] += ptA[2+3*i]*(ptB[2+3*i]) ;
    #endif    
  }
}

void cpu_code_v1( const int n1, const int n2, const int n3,
		  const double _Complex *host_coeffs1, 
		  const double _Complex *host_coeffs2, 
		  const double _Complex *host_coeffs3,
		  const int nMom,
		  const double _Complex *host_mom, 
		  const int nEv,
		  const void *const *host_evec, 
		  double _Complex *return_array,
		  const int blockSizeMomProj,
		  const int X[4] )
{
  const size_t nSp = X[0]*X[1]*X[2] ;
  const size_t nsites = nSp*X[3] ;
  double _Complex *q1 = (double _Complex*)calloc( n1*nsites*3 , sizeof(double _Complex) ) ;
  double _Complex *q2 = (double _Complex*)calloc( n2*nsites*3 , sizeof(double _Complex) ) ;
  double _Complex *q3 = (double _Complex*)calloc( n3*nsites*3 , sizeof(double _Complex) ) ;

  std::vector<size_t> ndil = { (size_t)n1 , (size_t)n2 , (size_t)n3 } ;
  std::vector<double _Complex*> q = { q1 , q2 , q3 } ;
  std::vector<const double _Complex*> coeffs = { host_coeffs1 , host_coeffs2 , host_coeffs3 } ;
#pragma omp parallel
  {
    //evprod( host_coeffs1 , host_evec , n1 , nEv , nsites , q1 ) ;
    //evprod( host_coeffs2 , host_evec , n2 , nEv , nsites , q2 ) ;
    //evprod( host_coeffs3 , host_evec , n3 , nEv , nsites , q3 ) ;
    evprodv2( coeffs , host_evec , ndil , nEv , nsites , q ) ;

  #pragma omp for collapse(2)
  for( int dil1 = 0 ; dil1 < n1 ; dil1++ ) {
    for( int dil2 = 0 ; dil2 < n2 ; dil2++ ) {
      LattField Diq( FieldSiteType::ColorVector);
      LattField tmp( FieldSiteType::Complex);      
      cpuColorCross( (void*)&q1[dil1*nsites*3] , (void*)&q2[dil2*nsites*3] , (void*)Diq.getDataPtr() , X ) ;
      for( int dil3 = 0 ; dil3 < n3 ; dil3++ ) {
	cpuColorContract( (void*)Diq.getDataPtr() , (void*)&q3[dil3*nsites*3] , (void*)tmp.getDataPtr() , X ) ;
	// DFT
	for( int cb = 0 ; cb < 2 ; cb++ ) {
	  for( int p = 0 ; p < nMom ; p++ ) {
	    const double _Complex *p2 = (const double _Complex*)host_mom+nSp*p+cb*nSp/2 ;
	    for( int T = 0 ; T < X[3] ; T++ ) {
	      const double _Complex *p1 = (const double _Complex*)tmp.getDataPtr() + nSp*T + cb*nSp/2;
	      double _Complex sum = 0.0 ;
            #ifdef USE_OPENBLAS
	      sum = cblas_zdotu( nSp/2 , p2 , 1 , p1 , 1 ) ;
            #elif (defined USE_GSL_CBLAS)
	      cblas_zdotu( nSp/2 , p2 , 1 , p1 , 1 , &sum ) ;
            #else
	      for( size_t i = 0 ; i < nSp/2 ; i++ ) {
		sum += p2[i]*p1[i] ;
	      }
	    #endif
	      return_array[ T + X[3]*( p + nMom*(size_t)( dil3 + n3*( dil2 + n2*dil1 ) )) ] += sum ;
	    }
	  }
	}
      }
    }
  }
  }
  free(q1) ;
  free(q2) ;
  free(q3) ;
}

// copy Fourier twiddles to the device
static inline void
device_hostmom( const double _Complex *host_mom ,
		void *d_mom ,
		const size_t size ,
		const int precision )
{
  if( precision == QUDA_SINGLE_PRECISION ) {
    float _Complex *tmp = (float _Complex*)calloc( size , sizeof( float _Complex ) ) ;
    for( size_t i = 0 ; i < size ; i++ ) {
      tmp[i] = (float _Complex)host_mom[i] ;
    }
    qudaMemcpy(d_mom, tmp , size*2*precision, qudaMemcpyHostToDevice);  
    free( tmp ) ;
  } else {
    qudaMemcpy(d_mom, host_mom, size*2*precision, qudaMemcpyHostToDevice);  
  }
}

// return d_ret to the host handling different precisions
static inline void
hostreturn( const void *d_ret ,
	    double _Complex *return_array ,
	    const size_t size ,
	    const int precision )
{
  if( precision == QUDA_SINGLE_PRECISION ) {
    float _Complex *tmp = (float _Complex*)calloc( size , sizeof( float _Complex ) ) ;
    qudaMemcpy(tmp, d_ret, size*2*precision, qudaMemcpyDeviceToHost);
    for( size_t i = 0 ; i < size ; i++ ) {
      return_array[i] = (double _Complex)tmp[i] ;
    }
    free( tmp ) ;
  } else {
    qudaMemcpy(return_array, d_ret, size*2*precision, qudaMemcpyDeviceToHost);  
  }
}

// device-side
static void
apply_noises( const std::vector<ColorSpinorField> &evec ,
	      const ColorSpinorParam cuda_evec_param ,
	      std::vector<std::vector<ColorSpinorField>> &q ,
	      const std::vector<std::vector<std::complex<double>>> &coeffs )
{
  const size_t nEv = evec.size() ;
  std::vector<ColorSpinorField> quda_evec(1) ;
  quda_evec[0] = ColorSpinorField(cuda_evec_param);
  for (size_t i=0; i<nEv; i++) {
    quda_evec[0] = evec[i] ;
    #pragma unroll
    for( size_t n = 0 ; n < q.size() ; n++ ) {
      const size_t n1 = q[n].size() ;
      blas::block::caxpy( {coeffs[n].begin()+n1*i,coeffs[n].begin()+n1*(i+1)},
			  {quda_evec[0]}, {q[n].begin(),q[n].end()} ) ;
    }
  }
}

// because we now have the same indexing pattern all our BLAS DFT calls are the same
static QudaBLASParam
default_BLAS( const int nMom , const int X[4] , const int blockSizeMomProj , const int precision )
{
  const int nSp = X[0]*X[1]*X[2] ;
  const int nSites = nSp*X[3] ;
  QudaBLASParam cublas_param = newQudaBLASParam() ;
  cublas_param.trans_a = QUDA_BLAS_OP_N;
  cublas_param.trans_b = QUDA_BLAS_OP_T;
  cublas_param.m = nMom ;
  cublas_param.n = X[3] ;
  cublas_param.k = nSp ;
  cublas_param.lda = nSp ;
  cublas_param.ldb = nSp ;
  cublas_param.ldc = X[3] ;
  cublas_param.a_stride = 0 ;
  cublas_param.b_stride = nSites ;
  cublas_param.c_stride = X[3]*nMom ;
  cublas_param.batch_count = blockSizeMomProj;
  cublas_param.alpha = 1. ; cublas_param.beta = 0. ;
  cublas_param.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param.data_type = ( precision == QUDA_SINGLE_PRECISION ) ? \
    QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
  cublas_param.blas_type = QUDA_BLAS_GEMM ;
  return cublas_param ;
}

static void
alamode2( const int n1, const int n2, const int n3,
	  const double _Complex *host_coeffs1, 
	  const double _Complex *host_coeffs2, 
	  const double _Complex *host_coeffs3,
	  const int nMom,
	  const double _Complex *host_mom, 
	  const int nEv,
	  void **host_evec,
	  QudaInvertParam inv_param,
	  double _Complex *return_array,
	  const int blockSizeMomProj,
	  const int X[4] )
{
  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_TOTAL);
  // appropriate checks and balances
  //if( sizeof(Complex) != sizeof(double _Complex) ) {
  //  errorQuda("Irreconcilable difference between interface and internal complex number conventions");
  // }
  if( blockSizeMomProj > (n1*n2*n3) ) {
    errorQuda( "Block size mom proj %d > %d\n", blockSizeMomProj, n1*n2*n3 ) ;
  }
  if( inv_param.cuda_prec != QUDA_DOUBLE_PRECISION &&
      inv_param.cuda_prec != QUDA_SINGLE_PRECISION ) {
    errorQuda( "Unsupported device precision %d" , inv_param.cuda_prec ) ;
  }
  const size_t nSp    = X[0]*X[1]*X[2] ;
  const size_t nSites = nSp*X[3] ;
  const int precision = inv_param.cuda_prec ;
  const int nRHS = 16 , nDiq = 8 ;
  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_INIT);
  const lat_dim_t x = { X[0] , X[1] , X[2] , X[3] } ;
  ColorSpinorParam cpu_evec_param( host_evec, inv_param, x, false, QUDA_CPU_FIELD_LOCATION );
  cpu_evec_param.nSpin = 1;
  std::vector<ColorSpinorField> evec(nEv);
  for (int iEv=0; iEv<nEv; ++iEv) {
    cpu_evec_param.v = host_evec[iEv];
    evec[iEv] = ColorSpinorField(cpu_evec_param) ;
  }
  // evec parameters
  ColorSpinorParam cuda_evec_param( cpu_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION );
  cuda_evec_param.setPrecision( inv_param.cuda_prec, inv_param.cuda_prec, true );
  // Create q1, q2, and q3 temporaries
  ColorSpinorParam cuda_q1_param( cuda_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION );
  ColorSpinorParam cuda_q2_param( cuda_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION );
  ColorSpinorParam cuda_q3_param( cuda_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION );
  cuda_q1_param.create = cuda_q2_param.create = cuda_q3_param.create = QUDA_ZERO_FIELD_CREATE;
  std::vector<std::vector<std::complex<double>>> coeffs(3) ;
  std::vector<std::vector<ColorSpinorField>> quda_q(3) ;
  quda_q[0].resize(n1) ; coeffs[0].resize( n1*nEv) ;
  for(int i=0; i<n1; i++) {
    quda_q[0][i] = ColorSpinorField(cuda_q1_param) ;
    for( int j = 0 ; j < nEv ; j++ ) { coeffs[0][j*n1+i] = (std::complex<double>)host_coeffs1[j+i*nEv] ; }
  }
  quda_q[1].resize(n2) ; coeffs[1].resize( n2*nEv) ;
  for(int i=0; i<n2; i++) {
    quda_q[1][i] = ColorSpinorField(cuda_q2_param) ;
    for( int j = 0 ; j < nEv ; j++ ) { coeffs[1][j*n2+i] = (std::complex<double>)host_coeffs2[j+i*nEv] ; }
  }
  quda_q[2].resize(n3) ; coeffs[2].resize( n3*nEv) ;
  for(int i=0; i<n3; i++) {
    quda_q[2][i] = ColorSpinorField(cuda_q3_param) ;
    for( int j = 0 ; j < nEv ; j++ ) { coeffs[2][j*n3+i] = (std::complex<double>)host_coeffs3[j+i*nEv] ; }
  }
  // device temporaries, momentum, and return buffers. All pretty small
  const size_t data_tmp_bytes = std::max( nRHS , blockSizeMomProj )*nSites*2*precision ;
  const size_t data_ret_bytes = blockSizeMomProj*nMom*X[3]*2*precision ;
  const size_t data_mom_bytes = nMom*nSp*2*precision ;
  void *d_tmp = pool_device_malloc(data_tmp_bytes);
  void *d_ret = pool_device_malloc(data_ret_bytes);
  void *d_mom = pool_device_malloc(data_mom_bytes);
  if( getVerbosity() >= QUDA_SUMMARIZE ) {
    printfQuda( "Tmp %f | ret %f | mom %f [GB]\n" ,
		data_tmp_bytes/OneGB , data_ret_bytes/OneGB , data_mom_bytes/OneGB ) ;
  }
  //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_INIT);  
  // Copy host_mom data to device
  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_H2D);
  device_hostmom( host_mom , d_mom , nMom*nSp , precision ) ;
  //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_H2D);
  // apply the noises to the evecs
  //getProfileApplyNoise().TPSTART(QUDA_PROFILE_COMPUTE);
  apply_noises( evec , cuda_evec_param , quda_q , coeffs ) ;
  //getProfileApplyNoise().TPSTOP(QUDA_PROFILE_COMPUTE);
  // usual momentum contraction, strided and blocked
  QudaBLASParam cublas_param_mom_sum =				\
    default_BLAS( nMom , X , blockSizeMomProj , precision ) ;
  // Create device diquark vector
  ColorSpinorParam cuda_diq_param( cuda_evec_param , inv_param , QUDA_CUDA_FIELD_LOCATION ) ;
  std::vector< ColorSpinorField > quda_diq( nDiq ) ;
  for( size_t i = 0 ; i < quda_diq.size() ; i++ ) {
    quda_diq[i] = ColorSpinorField( cuda_diq_param ) ;
  }
  int nInBlock = 0 , blockStart = 0 ;
  for( int dil1=0; dil1<n1; dil1++ ) {
    int dil2 = 0 ;
    while( dil2 < n2 ) {
      const int diqblk = std::min( nDiq , n2-dil2 ) ;
      //getProfileColorCross().TPSTART(QUDA_PROFILE_COMPUTE);
      colorCrossQudaV( quda_q[0][dil1], { quda_q[1].begin() + dil2 , quda_q[1].begin() + dil2 + diqblk } ,
		       { quda_diq.begin() , quda_diq.begin()+diqblk } ) ;
      //getProfileColorCross().TPSTOP(QUDA_PROFILE_COMPUTE);
      for( int b = 0 ; b < diqblk ; b++ ) {
	int dil3 = 0 ;
	while( dil3 < n3 ) {
	  //getProfileColorContract().TPSTART(QUDA_PROFILE_COMPUTE);
	  const int blk = std::min( std::min( n3 - dil3 , nRHS ) , (blockSizeMomProj-nInBlock ) ) ;
	  //getProfileColorContract().TPSTART(QUDA_PROFILE_COMPUTE);
	  colorContractQudaV( quda_diq[b], { quda_q[2].begin() + dil3 , quda_q[2].begin() + dil3 + blk } ,
			      (char*)d_tmp + nSites*nInBlock*2*precision);
	  //getProfileColorContract().TPSTOP(QUDA_PROFILE_COMPUTE);
	  nInBlock+=blk ; dil3 += blk ;
	  //getProfileColorContract().TPSTOP(QUDA_PROFILE_COMPUTE);
	  if (nInBlock == blockSizeMomProj ) {
	    //getProfileBLAS().TPSTART(QUDA_PROFILE_COMPUTE);
	    blas_lapack::native::stridedBatchGEMM(d_mom, d_tmp, (char*)d_ret,
						  cublas_param_mom_sum, QUDA_CUDA_FIELD_LOCATION);
	    //getProfileBLAS().TPSTOP(QUDA_PROFILE_COMPUTE);
	    //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_D2H);
	    hostreturn( d_ret , return_array + X[3]*nMom*blockStart ,
			(size_t)nInBlock*X[3]*nMom , precision ) ;
	    //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_D2H);
	    blockStart += nInBlock ;
	    nInBlock = 0;
	  }
	}
      }// diqblk
      dil2 += diqblk ;
    }
  }
  // overspill is less efficient than exact division but more flexible
  if( nInBlock > 0 ) {
    //getProfileBLAS().TPSTART(QUDA_PROFILE_COMPUTE);
    cublas_param_mom_sum.batch_count = nInBlock;
    blas_lapack::native::stridedBatchGEMM(d_mom, d_tmp, (char*)d_ret,
					  cublas_param_mom_sum, QUDA_CUDA_FIELD_LOCATION);
    //getProfileBLAS().TPSTOP(QUDA_PROFILE_COMPUTE);
    //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_D2H);
    hostreturn( d_ret , return_array + X[3]*nMom*blockStart ,
		(size_t)nInBlock*X[3]*nMom , precision ) ;
    //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_D2H);
  }
  // Clean up memory allocations
  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_FREE);
  pool_device_free(d_tmp);
  pool_device_free(d_mom);
  pool_device_free(d_ret);
  //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_FREE);
  //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_TOTAL);
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
	// will just give zero
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
  setVerbosityQuda(QUDA_VERBOSE, "#" , stdout ) ;
  
  // init Evs
#ifdef GPU_STRESS
  const int Nev = 256 , n1 = 64 , n2 = 64 , n3 = 64 ;
#else
  //const int Nev = 32 , n1 = 8 , n2 = 8 , n3 = 8 ;
  const int Nev = 64 , n1 = 32 , n2 = 32 , n3 = 32 ;
#endif
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);
  std::cout<<"Constant Eigvecs"<<std::endl ;
  set_constant( laphEigvecs ) ;
  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 32 ;
  const int X[4] = {
    LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int nsp = X[0]*X[1]*X[2] ;
  printf( "Nmom %d | (n1,n2,n3) %d,%d,%d\n" , nmom , n1 , n2 , n3 ) ;
  
  double _Complex coeffs1[ Nev*n1 ] = {} ;
  double _Complex coeffs2[ Nev*n2 ] = {} ;
  double _Complex coeffs3[ Nev*n3 ] = {} ;
  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt ;
  for( size_t i = 0 ; i < Nev*n1 ; i++ ) {
    coeffs1[i] = unif(mt) + I*unif(mt) ;
  }
  for( size_t i = 0 ; i < Nev*n2 ; i++ ) {
    coeffs2[i] = unif(mt) + I*unif(mt) ;
  }
  for( size_t i = 0 ; i < Nev*n3 ; i++ ) {
    coeffs3[i] = unif(mt) + I*unif(mt) ;
  }  
  
  double _Complex *host_mom = (double _Complex*)calloc(nsp*nmom,sizeof(double _Complex) ) ;
    // host_mom should be complex
  for( size_t p = 0 ; p < nmom ; p++ ) {
    const size_t mom[3] = { p+1, p+2, p+3 } ;
    for( size_t i = 0 ; i < (size_t)nsp ; i++ ) {
      const size_t x = i%X[0] ;
      const size_t y = (i/X[0])%X[1] ;
      const size_t z = (i/(X[0]*X[1]))%X[2] ;
      const double arg = 2*M_PI*( mom[0]*x/(double)X[0] +
				  mom[1]*y/(double)X[1] +
				  mom[2]*z/(double)X[2] ) ;
      double s = 1 ,c = 1 ;
      sincos( arg , &s , &c ) ;
      host_mom[ i + nsp*p ] = c+I*s ;
    }
  }
  
  QudaInvertParam inv_param = newQudaInvertParam();
  inv_param.dslash_type = QUDA_WILSON_DSLASH;
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  inv_param.solve_type = QUDA_DIRECT_SOLVE;
  inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
  inv_param.cuda_prec = QUDA_SINGLE_PRECISION;
  inv_param.dirac_order = QUDA_DIRAC_ORDER;
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

  double _Complex *retGPU = (double _Complex*)calloc( X[3]*n1*n2*n3*nmom , sizeof( double _Complex) ) ;

#ifdef GPU_STRESS  
  for( int blockSizeMomProj = 2 ; blockSizeMomProj < 8192 ; blockSizeMomProj *= 2 ) {
    memset( retGPU , 0.0 , X[3]*nmom*n1*n2*n3*sizeof(double _Complex)) ;
#else
    const int blockSizeMomProj = 2048 ;
#endif

    const size_t ndil[3] = { n1 , n2 , n3 } ;
    const double _Complex *host[3] = { coeffs1, coeffs2 , coeffs3 } ;
    nKernel(
	    ndil , host,
	    nmom,
	    host_mom ,
	    Nev ,
	    evList.data(),
	    inv_param ,
	    retGPU,
	    blockSizeMomProj,
	    X , 3 ) ;
    
    double GPUtime ;
    printf( "blockSizeMomProj %d\n" , blockSizeMomProj ) ;
    //for( int Np = 1 ; Np <= 512 ; Np *=2 ) {
      StopWatch GPU ;
      GPU.start() ;
      nKernel(
	      ndil , host ,
	      nmom,
	      host_mom ,
	      Nev ,
	      evList.data(),
	      inv_param ,
	      retGPU,
	      blockSizeMomProj,
	      X , 3 ) ;
      GPU.stop() ;
      GPUtime = GPU.getTimeInSeconds() ;
      printLaph(make_strf("\nGPU baryonkernel (%d) in = %g seconds\n", 1 , GPUtime )) ;
      //}
#ifdef GPU_STRESS
  }
#elif (defined CPUCROSSCHECK)
  double _Complex *retCPU = (double _Complex*)calloc( X[3]*n1*n2*n3*nmom , sizeof( double _Complex) ) ;
  
  memset( retCPU , 0 , X[3]*n1*n2*n3*nmom*sizeof(double _Complex) ) ;
  StopWatch CPU ;
  CPU.start() ;
  cpu_code_v1( n1 , n2 , n3 ,
	       coeffs1, 
	       coeffs2, 
	       coeffs3,
	       nmom,
	       host_mom, 
	       Nev,
	       evList.data(),
	       retCPU,
	       blockSizeMomProj,
	       X ) ;
  CPU.stop() ;
  const double CPUtime = CPU.getTimeInSeconds() ;
  printLaph(make_strf("\nCPU baryonkernel in = %g seconds\n", CPUtime)) ;

  
  printf( "\n*************************************\n" ) ;
  printf( "-----> GPU speedup factor %gx\n" , CPUtime/GPUtime ) ;
  printf( "*************************************\n\n" ) ;

  // test outputs
  printf( "CPU == GPU\n" ) ;
  for( size_t p = 0 ; p < nmom ; p++ ) {
    double sum = 0 ;
    for( size_t T = 0 ; T < (size_t)X[3] ; T++ ) {
      for( size_t i = 0 ; i < n1*n2*n3 ; i++ ) {
	const size_t idx = T + X[3]*(i + (n1*n2*n3)*p) ;
	sum += cabs(( retCPU[idx] - retGPU[idx] )/retCPU[idx] );
        #ifdef VERBOSE_COMPARISON
	printf( " (%f %f) == (%f %f)\n" ,
		creal(retCPU[idx]) , cimag(retCPU[idx]) ,
		creal(retGPU[idx]) , cimag(retGPU[idx]) ) ;
        #endif
      }
    }
    std::cout<<"Summed diff p="<<p<<" "<<sum/(n1*n2*n3*X[3])<<std::endl ;
  }
  free( retCPU ) ;
#endif
  free(host_mom);
  free( retGPU ) ;
  finalize();

  return 0;
}
