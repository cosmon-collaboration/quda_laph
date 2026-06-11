#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

#include <quda.h>
#include <timer.h>
#include <blas_lapack.h>
#include <blas_quda.h>
#include <tune_quda.h>
#include <color_spinor_field.h>
#include <contract_quda.h>

#include <cassert>
#include <complex.h>

#include <random>

#include <omp.h>

using namespace LaphEnv ;
using namespace quda ;

// will give zero for these otherwise default to uniform random numbers
//#define PSEUDOCONSTANT
//#define VERBOSE_COMPARISON
//#define GPU_STRESS
#define CPUCROSSCHECK

// cpu color cross
static void
cpuColorCross( void *A , void *B , void *result , const int X[4] )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ;
  std::complex<double> *ptA = (std::complex<double>*)A ;
  std::complex<double> *ptB = (std::complex<double>*)B ;
  std::complex<double> *ptC = (std::complex<double>*)result ;
  //#pragma omp parallel for
  for( size_t i = 0 ; i < (size_t)Nsites ; i++ ) {
    ptC[ 3*i + 0 ] =  ptA[ 3*i + 1 ]*ptB[ 3*i + 2 ] - ptA[ 3*i + 2 ]*ptB[ 3*i + 1 ] ;
    ptC[ 3*i + 1 ] = -ptA[ 3*i + 0 ]*ptB[ 3*i + 2 ] + ptA[ 3*i + 2 ]*ptB[ 3*i + 0 ] ;
    ptC[ 3*i + 2 ] =  ptA[ 3*i + 0 ]*ptB[ 3*i + 1 ] - ptA[ 3*i + 1 ]*ptB[ 3*i + 0 ] ;
  }
}

// cpu color contract
static void
cpuInnerProduct( void *A , void *B , void *result , const int X[4] )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ;
  const double _Complex *ptA = (const double _Complex*)A ;
  const double _Complex *ptB = (const double _Complex*)B ;
  double _Complex *ptC = (double _Complex*)result ;
  for( size_t i = 0 ; i < (size_t)Nsites ; i++ ) {
    #ifdef USE_OPENBLAS
    ptC[i] = cblas_zdotc( 3 , ptA+3*i , 1 , ptB+3*i , 1 ) ;
    #elif (defined USE_GSL_CBLAS)
    cblas_zdotc( 3 , ptA+3*i , 1 , ptB+3*i , 1 , ptC+i ) ;
    #else
    ptC[i]  = conj(ptA[0+3*i])*(ptB[0+3*i]) ;
    ptC[i] += conj(ptA[1+3*i])*(ptB[1+3*i]) ;
    ptC[i] += conj(ptA[2+3*i])*(ptB[2+3*i]) ;
    #endif
  }
}

static inline void
evprodv2( const std::vector<size_t> &ndil ,
	  const std::vector<const double _Complex *> &coeffs ,
	  const void *const *host_evec ,
	  const size_t nEv ,
	  const size_t nsites ,
	  std::vector<double _Complex *>q )
{
#pragma omp for
  for( size_t i = 0 ; i < nsites ; i++ ) {
    //const double _Complex alpha = 1. , beta = 0. ;
    double _Complex pt2[3][nEv] ;
    for( size_t ev = 0 ; ev < nEv ; ev++ ) {
      const double _Complex *pt = (const double _Complex*)host_evec[ev]+3*i ;
      pt2[0][ev] = *(pt+0) ;
      pt2[1][ev] = *(pt+1) ;
      pt2[2][ev] = *(pt+2) ;
    }
    for( size_t nq = 0 ; nq < q.size() ; nq++ ) {
      for( size_t dil = 0 ; dil < ndil[nq] ; dil++ ) {
	const double _Complex *pc = (const double _Complex*)coeffs[nq]+dil*nEv ;
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

// slow loopy version
static void
cpu_code( const std::vector<size_t> &nDil,
	  const std::vector<const double _Complex*> &host_coeffs,
	  const int nMom,
	  const double _Complex *host_mom,
	  const int nEv,
	  void **host_evec, 
	  double _Complex *return_arr,
	  const size_t blockSizeMomProj,
	  const int X[4])
{
  assert( host_coeffs.size() == 4 && ndil.size() == 4 ) ;
  // spatial only
  const size_t nSp   = X[0]*X[1]*X[2] ;
  const size_t nsites = nSp*X[3] ;
  const int n1 = nDil[0] , n2 = nDil[1] , n3 = nDil[2] , n4 = nDil[3] ;
  double _Complex *q1 = (double _Complex*)calloc( n1*nsites*3 , sizeof(double _Complex) ) ;
  double _Complex *q2 = (double _Complex*)calloc( n2*nsites*3 , sizeof(double _Complex) ) ;
  double _Complex *q3 = (double _Complex*)calloc( n3*nsites*3 , sizeof(double _Complex) ) ;
  double _Complex *q4 = (double _Complex*)calloc( n4*nsites*3 , sizeof(double _Complex) ) ;
  std::vector<double _Complex*> q = { q1 , q2 , q3 , q4 } ;
#pragma omp parallel
  {
    evprodv2( nDil , host_coeffs , host_evec , nEv , nsites , q ) ;

    #pragma omp for collapse(2)
    for( size_t aEv = 0 ; aEv < nDil[0] ; aEv++ ) {
      for( size_t bEv = 0 ; bEv < nDil[1] ; bEv++ ) {

	LattField Diq1( FieldSiteType::ColorVector ) , Diq2( FieldSiteType::ColorVector ) ;
	LattField tmp( FieldSiteType::Complex);
	cpuColorCross( q[0]+aEv*nsites*3 , q[1]+bEv*nsites*3 , (void*)Diq1.getDataPtr() , X ) ;

	for( size_t cEv = 0 ; cEv < nDil[2] ; cEv++ ) {
	  for( size_t dEv = 0 ; dEv < nDil[3] ; dEv++ ) {

	    cpuColorCross( q[2]+cEv*nsites*3 , q[3]+dEv*nsites*3 , (void*)Diq2.getDataPtr() , X ) ;

	    // inner product cd are conjugated
	    cpuInnerProduct( Diq1.getDataPtr() , Diq2.getDataPtr() , tmp.getDataPtr() , X ) ;
	    
	    // Matrix mul here mom*d_tmp -> d_ret
	    for( int cb = 0 ; cb < 2 ; cb++ ) {
	      for( int p = 0 ; p < nMom ; p++ ) {
		const double _Complex *p2 = (const double _Complex*)host_mom+nSp*p+cb*nSp/2 ;
		for( int T = 0 ; T < X[3] ; T++ ) {
		  const double _Complex *p1 = (const double _Complex*)tmp.getDataPtr() + nSp*T + cb*nSp/2 ;
		  double _Complex sum = 0.0 ;
                #ifdef USE_OPENBLAS
		  sum = cblas_zdotc( nSp/2 , p1 , 1 , p2 , 1 ) ;
                #elif (defined USE_GSL_CBLAS)
		  cblas_zdotc( nSp/2 , p1 , 1 , p2 , 1 , &sum ) ;
                #else
		  for( size_t i = 0 ; i < nSp/2 ; i++ ) {
		    sum += conj(p2[i])*p1[i] ;
		  }
                #endif
		  const size_t midx = dEv + nDil[3]*( cEv + nDil[2]*( bEv + nDil[1]*aEv ) ) ;
		  return_arr[ T + X[3]*( p + nMom*midx ) ] += sum ;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  free( q1 ) ; free( q2 ) ; free( q3 ) ; free( q4 ) ;
}

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
default_BLAS( const int nMom , const int X[4] , const size_t blockSizeMomProj , const int precision , const bool bDag = false )
{
  const int nSp = X[0]*X[1]*X[2] ;
  const int nSites = nSp*X[3] ;
  QudaBLASParam cublas_param = newQudaBLASParam() ;
  cublas_param.trans_a = QUDA_BLAS_OP_N;
  cublas_param.trans_b = bDag ? QUDA_BLAS_OP_C : QUDA_BLAS_OP_T;
  cublas_param.m = nMom ;
  cublas_param.n = X[3] ;
  cublas_param.k = nSp ;
  cublas_param.lda = nSp ;
  cublas_param.ldb = nSp ;
  cublas_param.ldc = X[3] ;
  cublas_param.a_stride = 0 ;
  cublas_param.b_stride = nSites ;
  cublas_param.c_stride = X[3]*nMom ;
  cublas_param.batch_count = (int)blockSizeMomProj;
  cublas_param.alpha = 1. ; cublas_param.beta = 0. ;
  cublas_param.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param.data_type = ( precision == QUDA_SINGLE_PRECISION ) ? \
    QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z;
  cublas_param.blas_type = QUDA_BLAS_GEMM ;
  return cublas_param ;
}

// compute the DFT
static inline void
doBlasReturn( QudaBLASParam cublas_dft ,
	      void *d_tmp , void *d_ret , void *d_mom , 
	      double _Complex *return_array ,
	      size_t &nInBlock , size_t &blockStart ,
	      const int nMom , const int X[4] , const int precision )
{
  cublas_dft.batch_count = (int)nInBlock;
  //getProfileBLAS().TPSTART(QUDA_PROFILE_COMPUTE);  
  blas_lapack::native::stridedBatchGEMM(d_mom, d_tmp, d_ret,
					cublas_dft,
					QUDA_CUDA_FIELD_LOCATION);
  //getProfileBLAS().TPSTOP(QUDA_PROFILE_COMPUTE);
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_D2H);
  hostreturn( d_ret , return_array + X[3]*nMom*blockStart ,
	      (size_t)nInBlock*X[3]*nMom , precision ) ;
  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_D2H);
  blockStart += nInBlock; nInBlock = 0;
}

static void
allocateMomDevice( void **d_tmp , void **d_ret , void **d_mom ,
		   const size_t nRHS , const size_t blockSizeMomProj ,
		   const int X[4] , const int precision , const int nMom )
{
  const size_t nSp = X[0]*X[1]*X[2] ;
  const size_t nSites = nSp*X[3] ;
  const size_t data_tmp_bytes = std::max( nRHS , blockSizeMomProj )*nSites*2*precision;
  const size_t data_ret_bytes = nMom*blockSizeMomProj*X[3]*2*precision;
  const size_t data_mom_bytes = nMom*nSp*2*precision;
  *d_tmp = pool_device_malloc(data_tmp_bytes);
  *d_ret = pool_device_malloc(data_ret_bytes);
  *d_mom = pool_device_malloc(data_mom_bytes);
  if( getVerbosity() >= QUDA_SUMMARIZE ) {
    const size_t OneGB = 1024*1024*1024;
    const size_t total_bytes = data_tmp_bytes + data_ret_bytes + data_mom_bytes ;
    printfQuda("d_tmp %fGB | d_ret %fGB | d_mom %fGB | total = %fGB\n",
	       (double)data_tmp_bytes/(OneGB), (double)data_ret_bytes/(OneGB),
	       (double)data_mom_bytes/(OneGB), (double)total_bytes/(OneGB)); 
  }
}

// so I tied an onion to my belt, which was the style at the time
void alamode( const std::vector<size_t> &nDil ,
	      const std::vector<const double _Complex*> &host_coeffs ,
	      const int nMom,
	      const double _Complex *host_mom,
	      const int nEv,
	      void **host_evec, 
	      QudaInvertParam inv_param,
	      double _Complex *return_array,
	      const size_t blockSizeMomProj,
	      const int X[4])
{
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_TOTAL);
#if 1
  const size_t n1 = nDil[0] , n2 = nDil[1] , n3 = nDil[2] , n4 = nDil[3] ;
  const size_t nSp    = X[0]*X[1]*X[2];
  const size_t nSites = nSp*X[3] ;
  // these are hard-coded do not change
  const size_t nDiq = 8 ;
  
  // appropriate checks and balances
  //if (sizeof(Complex) != sizeof(double _Complex)) {
  //  errorQuda("Irreconcilable difference between interface and internal complex number conventions");
  //}
  //if( nEvChoose3%blockSizeMomProj != 0 ) {
  //  errorQuda("Block size mom proj needs to divide %zu %d", nEvChoose3 , blockSizeMomProj);
  //}
  if( blockSizeMomProj > n1*n2*n3*n4 ) { errorQuda("blockSizeMomProj too big %zu > %zu" , blockSizeMomProj , n1*n2*n3*n4 ) ; }
  
  if( inv_param.cuda_prec != QUDA_SINGLE_PRECISION &&
      inv_param.cuda_prec != QUDA_DOUBLE_PRECISION ) {
    errorQuda( "Unsupported device precision %d" , inv_param.cuda_prec ) ;
  }
  const int precision = inv_param.cuda_prec ;
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_INIT);

  // Parameter object describing evecs
  const lat_dim_t x = { X[0] , X[1] , X[2] , X[3] } ;
  ColorSpinorParam cpu_evec_param(host_evec, inv_param, x, false, QUDA_CPU_FIELD_LOCATION);
  cpu_evec_param.nSpin = 1;
  std::vector<ColorSpinorField> evec(nEv) ;
  for (int iEv=0; iEv<nEv; ++iEv) {
    cpu_evec_param.v = host_evec[iEv];
    evec[iEv] = ColorSpinorField(cpu_evec_param);
  }
  ColorSpinorParam cuda_evec_param(cpu_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_evec_param.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  // Create q temporaries
  std::vector<ColorSpinorParam> quda_q_param( nDil.size() ) ;
  std::vector<std::vector<std::complex<double>>> coeffs( host_coeffs.size() ) ;
  std::vector<std::vector<ColorSpinorField>> quda_q( host_coeffs.size() ) ;
  for( size_t d = 0 ; d < host_coeffs.size() ; d++ ) {
    quda_q_param[d] = ColorSpinorParam( cuda_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION );
    quda_q_param[d].create = QUDA_ZERO_FIELD_CREATE ;
    quda_q[d].resize( nDil[d] ) ; coeffs[d].resize( nDil[d]*nEv) ;
    for(size_t i=0; i<nDil[d]; i++) {
      quda_q[d][i] = ColorSpinorField( quda_q_param[d] ) ;
      for( int j = 0 ; j < nEv ; j++ ) { coeffs[d][i+j*nDil[d]] = (std::complex<double>)host_coeffs[d][j+i*nEv] ; }
    }
  }
  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_INIT);
  void *d_tmp = NULL , *d_ret = NULL , *d_mom = NULL ; 
  allocateMomDevice( &d_tmp , &d_ret , &d_mom , nDiq , blockSizeMomProj , X , precision , nMom ) ;

  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_H2D);
  device_hostmom( host_mom , d_mom , nMom*nSp , precision ) ;
  //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_H2D);
  
  QudaBLASParam cublas_dft =					\
    default_BLAS( nMom , X , blockSizeMomProj , precision , true ) ;

  apply_noises( evec , cuda_evec_param , quda_q , coeffs ) ;

  // Create device diquark vector
  ColorSpinorParam cuda_diq_param(cpu_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  std::vector< ColorSpinorField > quda_diq1( nDiq ) , quda_diq2( nDiq ) ;
  for( size_t i = 0 ; i < quda_diq1.size() ; i++ ) {
    quda_diq1[i] = ColorSpinorField( cuda_diq_param ) ;
    quda_diq2[i] = ColorSpinorField( cuda_diq_param ) ;
  }

  size_t nInBlock = 0, blockStart = 0;
  for (size_t aEv=0; aEv<n1; aEv++) {
    for(size_t bEv=0; bEv<n2; bEv++ ) {
      const int diqBlk1 = 1 ;
      colorCrossQudaV( quda_q[0][aEv],
		       { quda_q[1].begin() + bEv , quda_q[1].begin() + bEv + diqBlk1 } ,
		       { quda_diq1.begin() , quda_diq1.begin()+diqBlk1 } ) ;      
      for(size_t cEv=0; cEv<n3; cEv++ ) {
	size_t dEv = 0 ;
	while( dEv < n4 ) {
	  // blocking factor for multi RHS colorCross and inner product
	  const int diqBlk2 = std::min( std::min( n4 - dEv , nDiq ) , (blockSizeMomProj-nInBlock ) ) ; ;
	  colorCrossQudaV( quda_q[2][cEv],
			   { quda_q[3].begin() + dEv , quda_q[3].begin() + dEv + diqBlk2 } ,
			   { quda_diq2.begin() , quda_diq2.begin()+diqBlk2 } ) ;      
	  innerProductQudaV( quda_diq1[0], { quda_diq2.begin() , quda_diq2.begin()+diqBlk2 } ,
			     (char*)d_tmp+nSites*nInBlock*2*precision );
	  nInBlock += diqBlk2 ; dEv += diqBlk2 ;
	  if (nInBlock == blockSizeMomProj) {
	    doBlasReturn( cublas_dft , d_tmp , d_ret , d_mom , return_array ,
			  nInBlock , blockStart , nMom , X , precision ) ;
	  }
	}
      }
    }
  }
  // overspill, code is more efficient if you avoid this
  if( nInBlock > 0 ) {
    doBlasReturn( cublas_dft , d_tmp , d_ret , d_mom , return_array ,
		  nInBlock , blockStart , nMom , X , precision ) ;
  }

  // Copy return array back to host
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_TOTAL); 
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_D2H);
  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_D2H);
  
  // Clean up memory allocations
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_FREE);

  pool_device_free(d_tmp);
  pool_device_free(d_mom);
  pool_device_free(d_ret);

  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_FREE);
  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_TOTAL);
#endif
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

#ifdef OPENMP
  const int max_threads = omp_get_max_threads() ;
#else
  const int max_threads = 1 ;
#endif
  std::cout<< "Max threads here" << max_threads << std::endl ;

  int global = 1 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
#endif
  setVerbosityQuda(QUDA_VERBOSE, "#" , stdout ) ;
  
#ifdef GPU_STRESS
  const int Nev = 96 ;
#else
  const int Nev = 32 ;
  const std::vector<size_t> nDil = { 8,8,8,8 } ;
#endif
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);

  std::cout<<"Constant Eigvecs :: "<< Nev <<std::endl ;
  set_constant( laphEigvecs ) ;

  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 64 ;
  const int X[4] = {
    LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int nspat  = X[0]*X[1]*X[2] ;
  std::cout<<"nmom "<<nmom<<" | nEv "<<Nev<<std::endl ;
  
  double _Complex *host_mom = (double _Complex*)calloc( nmom*nspat , sizeof( double _Complex ) ) ;

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

  double _Complex coeffs1[ Nev*nDil[0] ] , coeffs2[ Nev*nDil[1] ] ;
  double _Complex coeffs3[ Nev*nDil[2] ] , coeffs4[ Nev*nDil[3] ] ;
  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt ;
  for( size_t i = 0 ; i < Nev*nDil[0] ; i++ ) coeffs1[i] = unif(mt) + I*unif(mt) ;
  for( size_t i = 0 ; i < Nev*nDil[1] ; i++ ) coeffs2[i] = unif(mt) + I*unif(mt) ;
  for( size_t i = 0 ; i < Nev*nDil[2] ; i++ ) coeffs3[i] = unif(mt) + I*unif(mt) ;
  for( size_t i = 0 ; i < Nev*nDil[3] ; i++ ) coeffs4[i] = unif(mt) + I*unif(mt) ;

  QudaInvertParam inv_param = newQudaInvertParam();
  inv_param.dslash_type = QUDA_WILSON_DSLASH;
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  inv_param.solve_type = QUDA_DIRECT_SOLVE;
  inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
  inv_param.cuda_prec = QUDA_DOUBLE_PRECISION;
  inv_param.dirac_order = QUDA_DIRAC_ORDER;
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

  double _Complex *retGPU = (double _Complex*)calloc( X[3]*nmom*nDil[0]*nDil[1]*nDil[2]*nDil[3] , sizeof( double _Complex  ) );

#ifdef GPU_STRESS
  for( int blockSizeMomProj = 2 ; blockSizeMomProj < 8192 ; blockSizeMomProj *= 2 ) {
    std::cout<<"nmom "<<nmom<<" | block "<<blockSizeMomProj<<std::endl ;
    memset( retGPU , 0.0 , X[3]*nmom*nDil[0]*nDil[1]*nDil[2]*nDil[3]*sizeof( double _Complex )) ;
#else
    const int blockSizeMomProj = 1024 ;
#endif

    //alamode(
    const size_t ndil[4] = { nDil[0] , nDil[1] , nDil[2] , nDil[3] } ;
    const double _Complex *host[4] = { coeffs1, coeffs2 , coeffs3 , coeffs4 } ;

    nKernel(
	    ndil , host,
	    nmom,
	    host_mom ,
	    Nev,
	    evList.data() ,
	    inv_param,
	    retGPU,
	    blockSizeMomProj,
	    X , 4 ) ;

    double GPUtime = 0 ;
    int NP = nmom ;
    //for( int NP = 1 ; NP <= 4096 ; NP*=2 ) {
      StopWatch gpu ;
      gpu.start() ;
      
      nKernel(
	      ndil , host ,
	      nmom,
	      host_mom ,
	      Nev,
	      evList.data() ,
	      inv_param,
	      retGPU,
	      blockSizeMomProj,
	      X ,4 ) ;
    gpu.stop() ;
    GPUtime = gpu.getTimeInSeconds() ;
    printLaph(make_strf("\nGPU modetripletA in = %d %g seconds\n", NP , GPUtime )) ;
    //}
#ifdef GPU_STRESS
  }
#elif (defined CPUCROSSCHECK)
  double _Complex *retCPU = (double _Complex*)calloc( X[3]*nmom*nDil[0]*nDil[1]*nDil[2]*nDil[3] , sizeof( double _Complex  ) );
  StopWatch cpu ;

  printf( "Alloc retCPU %e GB\n" , X[3]*nmom*(size_t)nDil[0]*nDil[1]*nDil[2]*nDil[3]/(1024*1024*1024.) ) ;
  cpu.start() ;
  const std::vector<const double _Complex*> host_coeffs = { coeffs1 , coeffs2 , coeffs3 , coeffs4 } ;
  cpu_code( nDil , host_coeffs, nmom, host_mom , Nev, evList.data(), retCPU, blockSizeMomProj , X ) ;
  cpu.stop() ;
  const double CPUtime = cpu.getTimeInSeconds() ;
  printLaph(make_strf("\nCPU modetripletA in = %g seconds\n", CPUtime ));
  
  printf( "\n*************************************\n" ) ;
  printf( "-----> GPU speedup factor %gx\n" , CPUtime/GPUtime ) ;
  printf( "*************************************\n\n" ) ;
  
  // test outputs
  printf( "CPU == GPU\n" ) ;
  for( size_t p = 0 ; p < nmom ; p++ ) {
    double sum = 0 ;
    for( int T = 0 ; T < X[3] ; T++ ) {
      for( size_t i = 0 ; i < nDil[0]*nDil[1]*nDil[2]*nDil[3] ; i++ ) {
	const size_t idx = T + X[3]*(i + nDil[0]*nDil[1]*nDil[2]*nDil[3]*p) ;
	sum += cabs(( retCPU[idx] - retGPU[idx] )/retCPU[idx] ) ;
#ifdef VERBOSE_COMPARISON
	printf( "%d %zu %zu (%f %f) == (%f %f)\n" , T , i , p ,
		creal(retCPU[idx]) , cimag(retCPU[idx]) ,
		creal(retGPU[idx]) , cimag(retGPU[idx]) ) ;
#endif
      }
    }
    std::cout<<"Summed diff p="<<p<<" "<<sum/(nDil[0]*nDil[1]*nDil[2]*nDil[3]*X[3])<<std::endl ;
  }
  free( retCPU ) ;
#endif
  free( retGPU ) ;
  
  finalize();

  return 0;
}
