/**
   Test the laphBaryonTripletB code in my fork of quda
 */
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

#include <random>
#include <cassert>
#include <complex.h>

using namespace LaphEnv ;
using namespace quda ;

//#define VERY_VERBOSE

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
  const int nSubEv = nEv ;

  // zgemm d_mtb*d_coeffs = d_q3 
  const double _Complex *mtb = host_mode_trip_buf ;

  // is the first product
  double _Complex *q3 = (double _Complex*)calloc( nMom*nSubEv*nEv*n3 , sizeof( double _Complex ) ) ;
  
  // eg this is doing C_ij \sum_k A_ik Bjk strided for nMom*nSubEv
  // C_ij is size Nmom
  // nmom*nSubEv*Ev*Ev | n3*nEv
#pragma omp parallel for collapse(3)
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
  int p ;
#pragma omp parallel for private(p)
  for( p = 0 ; p < nMom ; p++ ) {
    double _Complex *tmp = (double _Complex*)calloc( nSubEv*n2*n3 , sizeof( double _Complex ) ) ;
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
    free( tmp ) ;
  }
  free( q3 ) ;
}

static void
alamode( const int n1,
	 const int n2,
	 const int n3,
	 const int nMom,
	 const int nEv,
	 const double _Complex *host_coeffs1,
	 const double _Complex *host_coeffs2,
	 const double _Complex *host_coeffs3,
	 const double _Complex *host_mode_trip_buf,
	 double _Complex *return_array )
{
  //getProfileBaryonKernelModeTripletsB().TPSTART(QUDA_PROFILE_TOTAL);
  //getProfileBaryonKernelModeTripletsB().TPSTART(QUDA_PROFILE_INIT);
   
  // number of EV indices (in first position) that this rank deals with
  const int nSubEv = nEv ;

  // check we are safe to cast into a Complex (= std::complex<double>)
  //if (sizeof(Complex) != sizeof(double _Complex)) {
  //  errorQuda("Irreconcilable difference between interface and internal complex number conventions");
  //}  
   
  const size_t OneGB = 1024*1024*1024;
  const size_t data_coeffs1_bytes = n1*nEv*2*QUDA_DOUBLE_PRECISION;
  const size_t data_coeffs2_bytes = n2*nEv*2*QUDA_DOUBLE_PRECISION;
  const size_t data_coeffs3_bytes = n3*nEv*2*QUDA_DOUBLE_PRECISION;  
  const size_t data_q3_bytes      = nSubEv*nEv*n3*2*QUDA_DOUBLE_PRECISION;

  // only create one temp
  const size_t data_tmp_bytes = std::max( nSubEv*nEv*nEv , nSubEv*n2*n3 )*2*QUDA_DOUBLE_PRECISION ;
  const size_t data_ret_bytes = nMom*n1*n2*n3*2*QUDA_DOUBLE_PRECISION;

  // Allocate required memory
  const size_t total_bytes = data_tmp_bytes + data_q3_bytes + data_coeffs3_bytes
    +data_coeffs1_bytes+data_coeffs2_bytes+data_tmp_bytes+data_ret_bytes;
  void *d_coeffs1 = pool_device_malloc(data_coeffs1_bytes);
  void *d_coeffs2 = pool_device_malloc(data_coeffs2_bytes);
  void *d_ret     = pool_device_malloc(data_ret_bytes);  
  void *d_tmp     = pool_device_malloc(data_tmp_bytes);
  void *d_q3      = pool_device_malloc(data_q3_bytes);
  void *d_coeffs3 = pool_device_malloc(data_coeffs3_bytes);
  if (getVerbosity() >= QUDA_VERBOSE) {
    printfQuda("coeffs1 %gGB | coeffs2 %gGB | ret %gGB | total %gGB\n",
	       (double)data_coeffs1_bytes/OneGB,
	       (double)data_coeffs2_bytes/OneGB,
	       (double)data_ret_bytes/OneGB,
	       (double)total_bytes/OneGB);
    printfQuda("mtb %gGB | q3 %gGB | coeffs3 %gGB | total %gGB\n",
	       (double)data_tmp_bytes/(OneGB),
	       (double)data_q3_bytes/(OneGB),
	       (double)data_coeffs3_bytes/(OneGB),
	       (double)total_bytes/(OneGB)); 
  }
  
  //getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_INIT);

  //getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_COMPUTE);
     
  // Allocate the rest of the arrays
  //getProfileBaryonKernelModeTripletsB().TPSTART(QUDA_PROFILE_INIT);
  //getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_INIT);

  QudaBLASParam cublas_param_1 = newQudaBLASParam();
  cublas_param_1.trans_a = QUDA_BLAS_OP_N;
  cublas_param_1.trans_b = QUDA_BLAS_OP_T;
  cublas_param_1.m = nEv;
  cublas_param_1.n = n3;
  cublas_param_1.k = nEv;
  cublas_param_1.lda = nEv;
  cublas_param_1.ldb = nEv;
  cublas_param_1.ldc = n3;
  cublas_param_1.a_stride = nEv*nEv ;
  cublas_param_1.b_stride = 0 ;
  cublas_param_1.c_stride = n3*nEv ;
  cublas_param_1.batch_count = nSubEv;
  cublas_param_1.alpha = 1.0 ; cublas_param_1.beta = 0.0 ;
  cublas_param_1.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param_1.data_type = QUDA_BLAS_DATATYPE_Z;

  // Initialise teh final ZGEMMs
  //getProfileBaryonKernelModeTripletsB().TPSTART(QUDA_PROFILE_INIT);
  QudaBLASParam cublas_param_2 = newQudaBLASParam();
  cublas_param_2.trans_a = cublas_param_2.trans_b = QUDA_BLAS_OP_N;
  cublas_param_2.m = n2;
  cublas_param_2.n = n3;
  cublas_param_2.k = nEv;
  cublas_param_2.lda = nEv;
  cublas_param_2.ldb = n3;
  cublas_param_2.ldc = n3;
  cublas_param_2.a_stride = 0 ;
  cublas_param_2.b_stride = n3*nEv ;
  cublas_param_2.c_stride = n2*n3 ;
  cublas_param_2.batch_count = nSubEv ;
  cublas_param_2.alpha = 1.0 ; cublas_param_2.beta = 0.0 ;
  cublas_param_2.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param_2.data_type = QUDA_BLAS_DATATYPE_Z;

  QudaBLASParam cublas_param_3 = newQudaBLASParam();
  cublas_param_3.trans_a = cublas_param_3.trans_b = QUDA_BLAS_OP_N;
  cublas_param_3.m = n1;
  //cublas_param_3.n = n2*n3;
  cublas_param_3.n = n3;
  cublas_param_3.k = nSubEv;
  cublas_param_3.lda = nEv;
  cublas_param_3.ldb = n2*n3;
  cublas_param_3.ldc = n2*n3;

  cublas_param_3.a_stride = 0 ;
  cublas_param_3.b_stride = n3 ;
  cublas_param_3.c_stride = n3 ;
  
  cublas_param_3.batch_count = n2 ;
  cublas_param_3.alpha = 1.0 ; cublas_param_3.beta = 0.0 ;
  cublas_param_3.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param_3.data_type = QUDA_BLAS_DATATYPE_Z;
  //getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_INIT);

  // flush this guy
  qudaMemset( d_tmp , 0 , data_tmp_bytes ) ;

  //getProfileBaryonKernelModeTripletsB().TPSTART(QUDA_PROFILE_H2D);  
  qudaMemcpy(d_coeffs1, host_coeffs1, data_coeffs1_bytes, qudaMemcpyHostToDevice);  
  qudaMemcpy(d_coeffs2, host_coeffs2, data_coeffs2_bytes, qudaMemcpyHostToDevice);  
  qudaMemcpy(d_coeffs3, host_coeffs3, data_coeffs3_bytes, qudaMemcpyHostToDevice);  
  //getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_H2D);

  // Compute ZGEMMs
  //getProfileBaryonKernelModeTripletsB().TPSTART(QUDA_PROFILE_COMPUTE);  
  for(int p=0; p<nMom; p++) {

    qudaMemcpy( d_tmp, (double _Complex*)host_mode_trip_buf+p*nSubEv*nEv*nEv ,
	        data_tmp_bytes, qudaMemcpyHostToDevice ) ;

    blas_lapack::native::stridedBatchGEMM( d_tmp, d_coeffs3, d_q3, cublas_param_1,
					   QUDA_CUDA_FIELD_LOCATION ) ;

    // batched GEMM that actually does this in batches
    blas_lapack::native::stridedBatchGEMM( d_coeffs2, d_q3, d_tmp,
					   cublas_param_2, QUDA_CUDA_FIELD_LOCATION);

    blas_lapack::native::stridedBatchGEMM( (double _Complex*)d_coeffs1,
					   d_tmp,
					   (double _Complex*)d_ret + p*n1*n2*n3 ,
					   cublas_param_3, QUDA_CUDA_FIELD_LOCATION);
  }
  //getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_COMPUTE);
   
  // Copy return array to host
  //getProfileBaryonKernelModeTripletsB().TPSTART(QUDA_PROFILE_D2H);
  qudaMemcpy(return_array, d_ret, data_ret_bytes, qudaMemcpyDeviceToHost);  
  //getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_D2H);
      
  // Clean up all remaining memory allocations
  //getProfileBaryonKernelModeTripletsB().TPSTART(QUDA_PROFILE_FREE);
  pool_device_free(d_coeffs1);
  pool_device_free(d_coeffs2);
  pool_device_free(d_coeffs3);
  pool_device_free(d_tmp);
  pool_device_free(d_q3);
  pool_device_free(d_ret);  
  ///getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_FREE);

  //saveTuneCache();
  //getProfileBaryonKernelModeTripletsB().TPSTOP(QUDA_PROFILE_TOTAL);
}

int main(int argc, char *argv[]) {
  XMLHandler xml_in;

  if( init_quda_laph(argc, argv, xml_in) != 0 ) {
    exit(1) ;
  }
  setVerbosityQuda(QUDA_SUMMARIZE, "#" , stdout ) ;

  int global = 1 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
#endif
  assert( global == 1 ) ; // for now

  const int Nev = 128 ;
  const int NsubEv = Nev ;
  const int nmom = 40 , n1 = 16 , n2 = 16 , n3 = 16 ;

  static std::uniform_real_distribution<double> unif(0.0,1.0) ;
  std::mt19937 mt ;
  
  double _Complex host_coeffs1[ n1*Nev ] = {} ;
  for( int i = 0 ; i < n1*Nev ; i++ ) {
    host_coeffs1[ i ] = unif(mt) + I*unif(mt) ;
  }  
  double _Complex host_coeffs2[ n2*Nev ] = {} ;
  for( int i = 0 ; i < n2*Nev ; i++ ) {
    host_coeffs2[ i ] = unif(mt) + I*unif(mt) ;
  }
  double _Complex host_coeffs3[ n3*Nev ] = {} ;
  for( int i = 0 ; i < n3*Nev ; i++ ) {
    host_coeffs3[ i ] = unif(mt) + I*unif(mt) ;
  }
  double _Complex *host_mode_trip_buf = (double _Complex*)calloc( NsubEv*Nev*Nev*nmom , sizeof( double _Complex));
  for( int i = 0 ; i < NsubEv*Nev*Nev*nmom ; i++ ) {
    host_mode_trip_buf[ i ] = unif(mt) + I*unif(mt) ;
  }

  double _Complex host_ret_arr[ nmom*n1*n2*n3 ] = {} ;
  StopWatch gpu;
  gpu.start() ;
  //alamode
  laphBaryonKernelComputeModeTripletB
      ( n1, n2, n3,
	nmom,
	Nev ,
	host_coeffs1 ,
	host_coeffs2 ,
	host_coeffs3 ,
	host_mode_trip_buf ,
	host_ret_arr ) ;
  gpu.stop() ;
  const double GPUtime = gpu.getTimeInSeconds() ;
  printLaph(make_strf("\nGPU modetripletB in = %g seconds\n", GPUtime )) ;

  StopWatch cpu ;
  cpu.start() ;
  double _Complex cpu_ret_arr[ nmom*n1*n2*n3 ] = {} ;
  cpu_code( n1, n2, n3,
	    nmom,
	    Nev ,
	    host_coeffs1 ,
	    host_coeffs2 ,
	    host_coeffs3 ,
	    host_mode_trip_buf ,
	    cpu_ret_arr ) ;
  cpu.stop() ;
  const double CPUtime = cpu.getTimeInSeconds() ;
  printLaph(make_strf("\nCPU modetripletB in = %g seconds\n", CPUtime )) ;

  printf( "\n*************************************\n" ) ;
  printf( "-----> GPU speedup factor %gx\n" , CPUtime/GPUtime ) ;
  printf( "*************************************\n\n" ) ;
  
  for( int p = 0 ; p < nmom ; p++ ) {
    double Sum = 0. ;
    for( int dil = 0 ; dil < n1*n2*n3 ; dil++ ) {
      Sum += cabs( cpu_ret_arr[ dil + n1*n2*n3*p ] - host_ret_arr[ dil + n1*n2*n3*p ] ) ;
      #ifdef VERY_VERBOSE
      printf( "(%f,%f) == (%f,%f)\n" ,
	      creal( cpu_ret_arr[dil+n1*n2*n3*p] ) , cimag( cpu_ret_arr[dil+n1*n2*n3*p] ) ,
	      creal( host_ret_arr[dil+n1*n2*n3*p] ) , cimag( host_ret_arr[dil+n1*n2*n3*p] ) ) ;
      #endif
    }
    printf( "%d Sub %e\n" , p , Sum ) ;
  }

  free( host_mode_trip_buf ) ;
  finalize();

  return 0;
}
