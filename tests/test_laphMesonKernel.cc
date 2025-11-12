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

//#define VERBOSE_COMPARISON
//#define GPU_STRESS

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

static void
cpuInner( const double _Complex *host_quark , const double _Complex *host_quark_bar , void *result , const int X[4] )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ;
  const std::complex<double> *ptA = (const std::complex<double>*)host_quark_bar ;
  const std::complex<double> *ptB = (const std::complex<double>*)host_quark ;
  std::complex<double> *ptC = (std::complex<double>*)result ;  
  //#pragma omp parallel for
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

static void cpu_code_v2( const int n1,
			 const int n2,
			 const int nMom,
			 const int blockSizeMomProj,
			 const double _Complex *host_coeffs1, 
			 const double _Complex *host_coeffs2,
			 const int nEv,
			 const void *const *host_evec,
			 const double _Complex *host_mom,
			 QudaInvertParam inv_param,
			 double _Complex *return_array,
			 const int X[4])
{
  const size_t nSp = X[0]*X[1]*X[2] ;
  const size_t nsites = nSp*X[3] ;
  double _Complex *q1 = (double _Complex*)calloc( n1*nsites*3 , sizeof(double _Complex) ) ;
  double _Complex *q2 = (double _Complex*)calloc( n2*nsites*3 , sizeof(double _Complex) ) ;
#pragma omp parallel
  {
    evprod( host_coeffs1 , host_evec , n1 , nEv , nsites , q1 ) ;
    evprod( host_coeffs2 , host_evec , n2 , nEv , nsites , q2 ) ;
  }  
  double _Complex *rt = (double _Complex*)return_array ;
#pragma omp parallel for collapse(2)
  for( int dil1 = 0 ; dil1 < n1 ; dil1++ ) {
    for( int dil2 = 0 ; dil2 < n2 ; dil2++ ) {
      double _Complex *result = (double _Complex*)calloc( nsites , sizeof( double _Complex ) ) ;
      cpuInner( q2 + dil2*3*nsites , q1 + dil1*3*nsites , result , X ) ;
      for( int p = 0 ; p < nMom ; p++ ) {
	double _Complex *pm = (double _Complex*)host_mom + nSp*p ;
	for( int t = 0 ; t < X[3] ; t++ ) {
	  double _Complex *rs = (double _Complex*)result + nSp*t ;
	  double _Complex sum = 0. ;
          #ifdef USE_OPENBLAS
          sum = cblas_zdotu( nSp , pm , 1 , rs , 1 ) ;
          #elif (defined USE_GSL_CBLAS)
          cblas_zdotu( nSp , pm , 1 , rs , 1 , &sum ) ;
          #else
	  for( size_t i = 0 ; i < (size_t)nSp ; i++ ) {
	    sum += pm[i]*rs[i] ;
	  }
	  #endif
	  rt[ t + X[3]*( p + nMom*( dil2 + n2*dil1 )) ] = sum ;
	}
      }
      free( result ) ;
    }
  }
  free( q1 ) ; free( q2 ) ;
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

// new GPU interface with better behaviour
static void alamode( const int n1,
		     const int n2,
		     const int nMom,
		     const int blockSizeMomProj,
		     const double _Complex *host_coeffs1, 
		     const double _Complex *host_coeffs2,
		     const int nEv,
		     void **host_evec,
		     const double _Complex *host_mom,
		     QudaInvertParam inv_param,
		     double _Complex *return_array,
		     const int X[4])
{
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_TOTAL);
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_INIT);  
  // Check we are safe to cast into a Complex (= std::complex<double>)
  //if (sizeof(Complex) != sizeof(double _Complex)) {
  //  errorQuda("Irreconcilable difference between interface and internal complex number conventions");
  //}
  if( blockSizeMomProj > (n1*n2) ) {
    errorQuda("block_size_mom_proj %d > (n1*n2) %d" , blockSizeMomProj, n1*n2 ) ;
  }
  if( inv_param.cuda_prec != QUDA_DOUBLE_PRECISION &&
      inv_param.cuda_prec != QUDA_SINGLE_PRECISION ) {
    errorQuda("Unsupported device precision") ;
  }
  const QudaPrecision precision = inv_param.cuda_prec ;
  // Some common variables
  const size_t nSp = X[0]*X[1]*X[2];
  const size_t nSites = nSp*X[3];
  const lat_dim_t x = { X[0] , X[1] , X[2] , X[3] } ;
  ColorSpinorParam cpu_evec_param(host_evec, inv_param, x, false, QUDA_CPU_FIELD_LOCATION);
  cpu_evec_param.nSpin = 1;
  std::vector<ColorSpinorField> evec(nEv) ;
  for( int i = 0 ; i < nEv ; i++ ) {
    cpu_evec_param.v = host_evec[i] ;
    evec[i] = ColorSpinorField(cpu_evec_param) ;
  }
  ColorSpinorParam cuda_evec_param(cpu_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_evec_param.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  ColorSpinorParam cuda_q1_param( cuda_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION );
  ColorSpinorParam cuda_q2_param( cuda_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION );
  cuda_q1_param.create = cuda_q2_param.create = QUDA_ZERO_FIELD_CREATE;
  std::vector<std::vector<ColorSpinorField>> quda_q(2) ;
  std::vector<std::vector<std::complex<double>>> coeffs(2) ;
  quda_q[0].resize(n1) ; coeffs[0].resize(n1*nEv) ;
  for(int i=0; i<n1; i++) {
    quda_q[0][i] = ColorSpinorField(cuda_q1_param) ;
    for( int j = 0 ; j < nEv ; j++ ) { coeffs[0][j*n1+i] = (std::complex<double>)host_coeffs1[j+i*nEv] ; }
  }
  quda_q[1].resize(n2) ; coeffs[1].resize(n2*nEv) ;
  for(int i=0; i<n2; i++) {
    quda_q[1][i] = ColorSpinorField(cuda_q2_param) ;
    for( int j = 0 ; j < nEv ; j++ ) { coeffs[1][j*n2+i] = (std::complex<double>)host_coeffs2[j+i*nEv] ; }
  }
  // Device array to hold the entire return array
  const size_t data_ret_bytes = nMom*X[3]*n1*n2*2*precision;
  const size_t data_tmp_bytes = nSites*blockSizeMomProj*2*precision;
  const size_t data_mom_bytes = nMom*nSp*2*precision;
  void *d_ret = pool_device_malloc(data_ret_bytes);
  void *d_tmp = pool_device_malloc(data_tmp_bytes);
  void *d_mom = pool_device_malloc(data_mom_bytes);
  // momentum contractions are a batched strided BLAS
  QudaBLASParam cublas_param_mom_sum = newQudaBLASParam();
  cublas_param_mom_sum.trans_a = QUDA_BLAS_OP_N;
  cublas_param_mom_sum.trans_b = QUDA_BLAS_OP_T;
  cublas_param_mom_sum.m = nMom ;
  cublas_param_mom_sum.n = X[3] ;
  cublas_param_mom_sum.k   = nSp ;
  cublas_param_mom_sum.lda = nSp ;
  cublas_param_mom_sum.ldb = nSp ;
  cublas_param_mom_sum.ldc = X[3] ;
  cublas_param_mom_sum.a_stride = 0 ; // mom stays the same
  cublas_param_mom_sum.b_stride = nSp*X[3] ;
  cublas_param_mom_sum.c_stride = X[3]*nMom ;
  cublas_param_mom_sum.batch_count = blockSizeMomProj ;
  cublas_param_mom_sum.alpha = 1.0; cublas_param_mom_sum.beta = 0.0;
  cublas_param_mom_sum.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param_mom_sum.data_type = (inv_param.cuda_prec == QUDA_SINGLE_PRECISION) ? \
    QUDA_BLAS_DATATYPE_C : QUDA_BLAS_DATATYPE_Z ;
  //getProfileCurrentKernel().TPSTOP(QUDA_PROFILE_INIT);
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_H2D);
  device_hostmom( host_mom , d_mom , nMom*nSp , precision ) ;
  //getProfileCurrentKernel().TPSTOP(QUDA_PROFILE_H2D);
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_COMPUTE);
  apply_noises( evec , cuda_evec_param , quda_q , coeffs ) ;
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_COMPUTE);
  // doing too much work here as color contract (di1,dil2) == (dil2,dil1)*
  // so we could do a conjugated GEMM on the lower diagonal for the mom proj - TODO
  int nInBlock = 0 , blockStart = 0 ;
  for (int dil1=0; dil1<n1; dil1++) {
    for (int dil2=0; dil2<n2; dil2++){
      //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_COMPUTE);
      innerProductQuda( quda_q[0][dil1], quda_q[1][dil2], (char*)d_tmp+nSites*nInBlock*2*precision );
      //getProfileCurrentKernel().TPSTOP(QUDA_PROFILE_COMPUTE);
      nInBlock++ ;
      if( nInBlock == blockSizeMomProj ) {
	//getProfileBLAS().TPSTART(QUDA_PROFILE_COMPUTE);
	blas_lapack::native::stridedBatchGEMM( d_mom, d_tmp, (char*)d_ret+blockStart*X[3]*nMom*2*precision,
					       cublas_param_mom_sum, QUDA_CUDA_FIELD_LOCATION);
	blockStart += nInBlock ;
	nInBlock = 0 ;
	//getProfileBLAS().TPSTOP(QUDA_PROFILE_COMPUTE);
      }
    }
  }
  if( nInBlock > 0 ) {
    //getProfileBLAS().TPSTART(QUDA_PROFILE_COMPUTE);
    cublas_param_mom_sum.batch_count = nInBlock ;
    blas_lapack::native::stridedBatchGEMM( d_mom, d_tmp, (char*)d_ret+blockStart*X[3]*nMom*2*precision,
					   cublas_param_mom_sum, QUDA_CUDA_FIELD_LOCATION);
    //getProfileBLAS().TPSTOP(QUDA_PROFILE_COMPUTE);
  }
  // Copy device data back to host
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_D2H);
  hostreturn( d_ret , return_array , nMom*X[3]*n1*n2 , precision ) ;
  //getProfileCurrentKernel().TPSTOP(QUDA_PROFILE_D2H);
  // Clean up memory allocations
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_FREE);
  pool_device_free(d_ret);
  pool_device_free(d_tmp);
  pool_device_free(d_mom);
  //getProfileCurrentKernel().TPSTOP(QUDA_PROFILE_TOTAL);
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
  setVerbosityQuda(QUDA_VERBOSE, "#" , stdout ) ;
  
  // call rephase here
#ifdef GPU_STRESS
  const int Nev = 512 , n1 = 288 , n2 = 288 ;
#else
  const int Nev = 256 , n1 = 64 , n2 = 64 ;
#endif
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);
  set_constant( laphEigvecs ) ;

  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 32 ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int nspat  = X[0]*X[1]*X[2] ;
  std::cout<<"n1,n2 "<< n1 << "," << n2 << " | nmom" << nmom << std::endl ;

  double _Complex coeffs1[ Nev*n1 ] = {} , coeffs2[ Nev*n2 ] = {} ;
  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt ;
  for( size_t i = 0 ; i < Nev*n1 ; i++ ) {
    coeffs1[i] = unif(mt) + I*unif(mt) ;
  }
  for( size_t i = 0 ; i < Nev*n2 ; i++ ) {
    coeffs2[i] = unif(mt) + I*unif(mt) ;
  }
  
  // host_mom should be complex
  double _Complex *host_mom = (double _Complex*)calloc( nmom*nspat , sizeof(double _Complex) ) ;
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

  double _Complex *GPU_ret = (double _Complex*)calloc(n1*n2*nmom*X[3],sizeof(double _Complex)) ;

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

#ifdef GPU_STRESS
  for( int blockSizeMomProj = 2 ; blockSizeMomProj < 8192 ; blockSizeMomProj *= 2 ) {
    std::cout<< "block " << blockSizeMomProj << std::endl ;
    memset( GPU_ret , 0.0 , n1*n2*nmom*X[3]*sizeof(double _Complex));
#else
    const int blockSizeMomProj = 256 ;
#endif

    //alamode(
    laphMesonKernel(
	    n1,n2, nmom, blockSizeMomProj,
	    coeffs1 , coeffs2 ,
	    Nev, evList.data() , 
	    host_mom ,
	    inv_param ,
	    GPU_ret ,
	    X ) ;

    
    // GPU version
    double GPUtime = 0 ;
    //for( int NP = 1 ; NP < nmom ; NP*=2 ) { 
      StopWatch gpu ;
      gpu.start() ;
      //alamode(
      laphMesonKernel(
	      n1,n2, nmom, blockSizeMomProj,
	      coeffs1 , coeffs2 ,
	      Nev, evList.data() , 
	      host_mom ,
	      inv_param ,
	      GPU_ret ,
	      X ) ;
      gpu.stop();
      GPUtime = gpu.getTimeInSeconds();
      printLaph(make_strf("\nGPU (NP%d) current kernel in = %g seconds\n", 1 , GPUtime)) ;
      //}
#ifdef GPU_STRESS
  }
#else
  double _Complex *CPU_ret = (double _Complex*)calloc(n1*n2*nmom*X[3],sizeof(double _Complex)) ;
  // CPU version
  double CPUtime = 0. ;
  //for( int NP = 1 ; NP <= nmom ; NP*=2 ) { 
    StopWatch cpu ;
    cpu.start() ;
    cpu_code_v2(
	      n1,n2, nmom, blockSizeMomProj,
	      coeffs1 , coeffs2 ,
	      Nev, evList.data() , 
	      host_mom ,
	      inv_param ,
	      CPU_ret ,
	      X ) ;
    cpu.stop() ;
    CPUtime = cpu.getTimeInSeconds() ;
    printLaph(make_strf("\nCPU current (NP%d) kernel in = %g seconds\n", 1 , CPUtime));
    //  }
  printf( "\n*************************************\n" ) ;
  printf( "-----> GPU speedup factor %gx\n" , CPUtime/GPUtime ) ;
  printf( "*************************************\n\n" ) ;
  for( int p = 0 ; p < nmom ; p++ ) {
    double sum = 0.0 ;
    for( int dil1 = 0 ; dil1 < n1 ; dil1++ ) {
      for( int dil2 = 0 ; dil2 < n2 ; dil2++ ) {
	#ifdef VERBOSE_COMPARISON
	printf( "(di1,dil2,p) %d,%d,%d\n" , dil1, dil2, p) ;
	#endif
	for( int t = 0 ; t < X[3] ; t++ ) {
	  #ifdef VERBOSE_COMPARISON
	  printf( "(%f %f) == (%f %f)\n" ,
		  creal( CPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] ) ,
		  cimag( CPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] ) ,
		  creal( GPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] ) ,
		  cimag( GPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] )
		  ) ;
	  #endif
	  sum += cabs( (CPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] - GPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )])/CPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1))] ) ; 
	}
      }
    }
    printf( "diff %e\n" , sum/(n1*n2*nmom*X[3]) ) ;
  }
  free( CPU_ret ) ;
#endif
  free( host_mom ) ;
  free( GPU_ret ) ;
  finalize( ) ;
  return 0;
}
