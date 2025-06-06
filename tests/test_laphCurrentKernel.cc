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
  double _Complex *loc_sum = (double _Complex*)calloc( X[3]*nMom*n1*n2 , sizeof( double _Complex) );  
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
  free( loc_sum ) ;
}

// new GPU interface with better behaviour
static void alamode( const int n1,
		     const int n2,
		     const int n_mom,
		     const int block_size_mom_proj,
		     void **host_quark,
		     void **host_quark_bar,
		     const double _Complex *host_mom,
		     QudaInvertParam inv_param,
		     void *ret_arr,
		     const int X[4])
{
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_TOTAL);
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_INIT);
  
  // Check we are safe to cast into a Complex (= std::complex<double>)
  //  if (sizeof(Complex) != sizeof(double _Complex)) {
  //  errorQuda("Irreconcilable difference between interface and internal complex number conventions");
  //}

  if( (n2*n1)%block_size_mom_proj != 0 ) {
    errorQuda("I only support block sizes that are factors of n1*n2") ;
  }

  // Some common variables
  const size_t n_spatial_sites = X[0]*X[1]*X[2];
  const size_t n_sites = n_spatial_sites * X[3];
  const QudaPrecision precision = QUDA_DOUBLE_PRECISION;
  
  const lat_dim_t x = { X[0] , X[1] , X[2] , X[3] } ;

  // Create device vectors for quarks
  ColorSpinorParam cpu_quark_param(host_quark, inv_param, x, false, QUDA_CPU_FIELD_LOCATION);
  cpu_quark_param.nSpin = 1;
  std::vector<ColorSpinorField> quark(n2) ;
  for( int dil2 = 0 ; dil2 < n2 ; dil2++ ) {
    cpu_quark_param.v = host_quark[dil2] ;
    quark[dil2] = ColorSpinorField(cpu_quark_param) ;
  }
  ColorSpinorParam cuda_quark_param(cpu_quark_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_quark_param.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  
  // Create device vectors for quark_bar
  ColorSpinorParam cpu_quark_bar_param(host_quark_bar, inv_param, x, false, QUDA_CPU_FIELD_LOCATION);
  cpu_quark_bar_param.nSpin = 1;
  ColorSpinorParam cuda_quark_bar_param(cpu_quark_bar_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_quark_bar_param.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  
  // Device array to hold the entire return array
  const size_t data_ret_bytes = n_mom * X[3] * n1 * n2 * 2 * precision;
  void *d_ret = pool_device_malloc(data_ret_bytes);

  // Device array to hold the inner product
  const size_t data_tmp_bytes = block_size_mom_proj * n_sites * 2 * precision;
  void *d_tmp = pool_device_malloc(data_tmp_bytes);

  // Device array to hold the momentum
  const size_t data_mom_bytes = n_mom * n_spatial_sites * 2 * precision;
  void *d_mom = pool_device_malloc(data_mom_bytes);
  
  __complex__ double alpha = 1.0 , beta = 0.0;
  QudaBLASParam cublas_param_mom_sum = newQudaBLASParam();
  cublas_param_mom_sum.trans_a = QUDA_BLAS_OP_N;
  cublas_param_mom_sum.trans_b = QUDA_BLAS_OP_T;
  // going to be doing A.B where A is the host_mom array of length
  // nmom x nsites | nsites x T
  // output is an n_mom x X[3] matrix
  cublas_param_mom_sum.m = n_mom ; // # of rows of A -> mom list
  cublas_param_mom_sum.n = X[3] ;

  cublas_param_mom_sum.k   = n_spatial_sites ; // should be lda and ldb
  cublas_param_mom_sum.lda = n_spatial_sites ; // # of cols of A == L^3
  cublas_param_mom_sum.ldb = n_spatial_sites ; // #of rows of B == L^3

  cublas_param_mom_sum.ldc = X[3] ;
  
  cublas_param_mom_sum.a_stride = 0 ; // mom matrix stays the same
  cublas_param_mom_sum.b_stride = n_spatial_sites*X[3] ;
  cublas_param_mom_sum.c_stride = X[3]*n_mom ;

  cublas_param_mom_sum.batch_count = block_size_mom_proj ;
  cublas_param_mom_sum.alpha = (__complex__ double)alpha;  
  cublas_param_mom_sum.beta  = (__complex__ double)beta;
  cublas_param_mom_sum.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param_mom_sum.data_type = QUDA_BLAS_DATATYPE_Z;
  //--------------------------------------------------------------------------------

  // Copy host data to device
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_H2D);

  // pull all of quark here as it is doing more work in the loop qbar is loaded as needed
  std::vector<ColorSpinorField> quda_quark(n2) ;
  for (int dil2=0; dil2<n2; dil2++) {
    quda_quark[dil2] = ColorSpinorField(cuda_quark_param) ;
    quda_quark[dil2] = quark[dil2] ;
  }
  // For the moment, use the chroma_laph defined momenta, then compute on host
  qudaMemcpy(d_mom, host_mom, data_mom_bytes, qudaMemcpyHostToDevice);  
  //getProfileCurrentKernel().TPSTOP(QUDA_PROFILE_H2D);

  // doing too much work here as (di1,dil2) == (dil2,dil1)*
  int n_in_block = 0 , idx_last = 0 ;
  for (int dil1=0; dil1<n1; dil1++) {
    cpu_quark_bar_param.v = host_quark_bar[dil1] ;
    ColorSpinorField quark_bar(cpu_quark_bar_param) ;
    ColorSpinorField quda_quark_bar(cuda_quark_bar_param) ;
    quda_quark_bar = quark_bar ;
    // just block dil2
    for (int dil2=0; dil2<n2; dil2++){
      //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_COMPUTE);
      innerProductQuda( quda_quark_bar, quda_quark[dil2], (std::complex<double>*)d_tmp + n_sites*n_in_block );
      n_in_block++ ;
      //getProfileCurrentKernel().TPSTOP(QUDA_PROFILE_COMPUTE);

      ///getProfileBLAS().TPSTART(QUDA_PROFILE_COMPUTE);
      if( n_in_block == block_size_mom_proj ) {
	blas_lapack::native::stridedBatchGEMM( d_mom, d_tmp,
					       (std::complex<double>*)d_ret+(idx_last)*X[3]*n_mom,
					       cublas_param_mom_sum,
					       QUDA_CUDA_FIELD_LOCATION);
	//getProfileBLAS().TPSTOP(QUDA_PROFILE_COMPUTE);
	idx_last += block_size_mom_proj ;
	n_in_block = 0 ;
      }
    }
  }

  // Copy device data back to host
  qudaMemcpy(ret_arr, d_ret, data_ret_bytes, qudaMemcpyDeviceToHost) ;
  
  // Clean up memory allocations
  //getProfileCurrentKernel().TPSTART(QUDA_PROFILE_FREE);

  pool_device_free(d_ret);
  pool_device_free(d_tmp);
  pool_device_free(d_mom);
  
  //getProfileCurrentKernel().TPSTOP(QUDA_PROFILE_FREE);
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
  assert( global == 1 ) ;
  
  // call rephase here
#ifdef GPU_STRESS
  const int Nev = 256 , n1 = 64 , n2 = 64 ;
#else
  const int Nev = 64 , n1 = 16 , n2 = 16 ;
#endif
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);
  set_constant( laphEigvecs ) ;

  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 40 ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int nspat  = X[0]*X[1]*X[2] ;
  std::cout<<"n1,n2 "<< n1 << "," << n2 << " | nmom" << nmom << std::endl ;
  
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
  for( int blockSizeMomProj = 1 ; blockSizeMomProj < 4096 ; blockSizeMomProj *= 2 ) {
    std::cout<< "block " << blockSizeMomProj << std::endl ;
    memset( GPU_ret , 0.0 , n1*n2*nmom*X[3]*sizeof(double _Complex));
#else
    const int blockSizeMomProj = 32 ;
#endif
    // GPU version
    StopWatch gpu ;
    gpu.start() ;
    laphCurrentKernel( n1, n2,
		       nmom,
		       blockSizeMomProj,
		       evList.data() , 
		       evList.data() ,
		       host_mom ,
		       inv_param ,
		       GPU_ret ,
		       X ) ;
    gpu.stop();
    const double GPUtime = gpu.getTimeInSeconds();
    printLaph(make_strf("\nGPU current kernel in = %g seconds\n", GPUtime)) ;
#ifdef GPU_STRESS
  }
#else
  double _Complex *CPU_ret = (double _Complex*)calloc(n1*n2*nmom*X[3],sizeof(double _Complex)) ;
  // CPU version
  StopWatch cpu ;
  cpu.start() ;
  cpu_code( n1, n2,
	    nmom,
	    blockSizeMomProj,
	    evList.data() , 
	    evList.data() ,
	    host_mom ,
	    CPU_ret ,
	    X ) ;
  cpu.stop() ;
  const double CPUtime = cpu.getTimeInSeconds() ;
  printLaph(make_strf("\nCPU current kernel in = %g seconds\n", CPUtime));
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
	  sum += cabs( CPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] - GPU_ret[t+X[3]*(p+nmom*(dil2+n2*dil1) )] ) ; 
	}
      }
    }
    printf( "diff %e\n" , sum ) ;
  }
  free( CPU_ret ) ;
#endif
  free( host_mom ) ;
  free( GPU_ret ) ;
  finalize( ) ;
  return 0;
}
