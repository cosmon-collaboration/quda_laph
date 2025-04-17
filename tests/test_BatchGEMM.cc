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

// slow loopy version
void
cpu_code( const double _Complex *A ,
	  const int M , const int lda ,
	  const double _Complex *B ,
	  const int N , const int ldb ,
	  double _Complex *C ,
	  const int K , const int ldc ,
	  const int batch_count ,
	  const bool is_col_major )
{
  if( is_col_major == true ) {
    const int a_stride = M*K ;
    const int b_stride = K*N ;
    const int c_stride = M*N ;
    for( int p = 0 ; p < batch_count ; p++ ) {
      for( int m = 0 ; m < N ; m++ ) {
	for( int n = 0 ; n < M ; n++ ) {
	  double _Complex c_mnp = 0.0 ;
	  for( int k = 0 ; k < K ; k++ ) {
	    c_mnp += A[m+k*lda+p*a_stride]*B[k+n*ldb+p*b_stride] ;
	  }
	  C[m+n*ldc+p*c_stride] = c_mnp ;
	}
      }
    }
  } else {
    // row-major matrix mul
    const int a_stride = M*K ;
    const int b_stride = K*N ;
    const int c_stride = M*N ;
    for( int p = 0 ; p < batch_count ; p++ ) {
      for( int m = 0 ; m < N ; m++ ) {
	for( int n = 0 ; n < M ; n++ ) {
	  double _Complex c_mnp = 0.0 ;
	  for( int k = 0 ; k < K ; k++ ) {
	    c_mnp += A[k+m*lda+p*a_stride]*B[n+k*ldb+p*b_stride] ;
	  }
	  C[n+m*ldc+p*c_stride] = c_mnp ;
	}
      }
    }   
  }
}

void
gpu_code( const double _Complex *A ,
	  const int m , const int lda ,
	  const double _Complex *B ,
	  const int n , const int ldb ,
	  double _Complex *C ,
	  const int k , const int ldc ,
	  const int batch_count ,
	  const bool is_col_major )
{
  const size_t a_bytes = m*lda*batch_count*sizeof(double _Complex) ;
  const size_t b_bytes = n*ldb*batch_count*sizeof(double _Complex) ;
  const size_t c_bytes = m*n*batch_count*sizeof(double _Complex) ;
  void *d_A = pool_device_malloc( a_bytes ) ;
  void *d_B = pool_device_malloc( b_bytes ) ;
  void *d_C = pool_device_malloc( c_bytes ) ;
  qudaMemcpy( d_A , A , a_bytes , qudaMemcpyHostToDevice ) ;
  qudaMemcpy( d_B , B , b_bytes , qudaMemcpyHostToDevice ) ;

  // do a zgemm
  const __complex__ double alpha = 1.0, beta = 0.0;
  QudaBLASParam cublas_param = newQudaBLASParam();
  cublas_param.trans_a = QUDA_BLAS_OP_N;
  cublas_param.trans_b = QUDA_BLAS_OP_N;
  cublas_param.m = m;
  cublas_param.n = n;
  cublas_param.k = k;
  cublas_param.lda = lda;
  cublas_param.ldb = ldb;
  cublas_param.ldc = ldc;
  cublas_param.a_stride = 1 ; //m*k ;
  cublas_param.b_stride = 1 ; //k*n ;
  cublas_param.c_stride = 1 ; //m*n ;
  cublas_param.batch_count = batch_count;
  cublas_param.alpha = (__complex__ double)alpha;  
  cublas_param.beta  = (__complex__ double)beta;
  if( is_col_major ) {
    cublas_param.data_order = QUDA_BLAS_DATAORDER_COL;
  } else {
    cublas_param.data_order = QUDA_BLAS_DATAORDER_ROW;
  }
  cublas_param.data_type = QUDA_BLAS_DATATYPE_Z;

  blas_lapack::native::stridedBatchGEMM(d_A, d_B, d_C, cublas_param, QUDA_CUDA_FIELD_LOCATION);

  qudaMemcpy( C , d_C , c_bytes , qudaMemcpyDeviceToHost ) ;
  
  pool_device_free( d_A ) ;
  pool_device_free( d_B ) ;
  pool_device_free( d_C ) ;
}

static void
print_matrix( const double _Complex *A , const int M , const int lda , const int num_batch , const bool is_col_order )
{
  if( is_col_order == true ) {
    for( int p = 0 ; p < num_batch ; p++ ) {
      for( int m = 0 ; m < M ; m++ ) {
	for( int n = 0 ; n < lda ; n++ ) {
	  printf( "(%g %g) ", creal( A[m+M*n] ) , cimag( A[m+M*n] ) ) ;
	}
	printf( "\n" ) ;
      }
      A += M*lda ;
      printf( "===============\n" ) ;
    }
  } else {
    for( int p = 0 ; p < num_batch ; p++ ) {
      for( int m = 0 ; m < M ; m++ ) {
	for( int n = 0 ; n < lda ; n++ ) {
	  printf( "(%g %g) ", creal( A[n+M*m] ) , cimag( A[n+M*m] ) ) ;
	}
	printf( "\n" ) ;
      }
      A += M*lda ;
      printf( "===============\n" ) ;
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

  // do a batched blas call example from
  // CUDALibrarySamples/cuBLAS/Level-3/gemmStridedBatched/cublas_gemmStridedBatched_example.cu
  const int m = 2 ;
  const int n = 2 ;
  const int k = 2 ;
  const int lda = 2 ;
  const int ldb = 2 ;
  const int ldc = 2 ;
  const int batch_count = 2 ;

  printf( "Column major odering\n" ) ;
  {
    // should be of size m*lda*batch_count
    const double _Complex A[ 2*4 ] =		\
      { 1.0, 3.0, 2.0, 4.0,
	5.0, 7.0, 6.0, 8.0 } ;
    // should be of size n*ldb*batch_count? Does it need to be???
    const double _Complex B[ 2*4 ] =		\
      { 5.0, 7.0 , 6.0 , 8.0,
	9.0, 11.0, 10.0, 12.0 } ;
    printf( "A - matrix\n" ) ;
    print_matrix( A , m , k , batch_count , true ) ;
    
    printf( "B - matrix\n" ) ;
    print_matrix( B , k , n , batch_count , true ) ;
    
    double _Complex C[ batch_count*m*n ] = {} ;
    gpu_code( A , m , lda ,
	      B , n , ldb ,
	      C , k , ldc ,
	      batch_count , true ) ;
    printf( "\nC - matrix (Quda result)\n" ) ;
    /* Should be according to cublas routines
     *
     *   C = | 19.0 | 22.0 | 111.0 | 122.0 |
     *       | 43.0 | 50.0 | 151.0 | 166.0 |
     */
    print_matrix( C , m , n , batch_count , true ) ;
    
    printf( "\n" ) ;
    memset( C , 0.0 , batch_count*m*n*sizeof(double _Complex));
    cpu_code( A , m , lda ,
	      B , n , ldb ,
	      C , k , ldc ,
	      batch_count , true ) ;
    printf( "C - matrix (cpu)\n" ) ;
    print_matrix( C , m , n , batch_count , true ) ;
  }
  printf( "\nRow major ordering\n" ) ;
  {
    // OK that was all column-major now redo it being row-major in Quda and do a CPU version
    const double _Complex A[ 2*4 ] =		\
      { 1.0, 2.0, 3.0, 4.0,
	5.0, 6.0, 7.0, 8.0 } ;
    // should be of size n*ldb*batch_count? Does it need to be???
    const double _Complex B[ 2*4 ] =		\
      { 5.0, 6.0 , 7.0 , 8.0,
	9.0, 10.0, 11.0, 12.0 } ;
    printf( "A - matrix\n" ) ;
    print_matrix( A , m , k , batch_count , false ) ;
    
    printf( "B - matrix\n" ) ;
    print_matrix( B , k , n , batch_count , false ) ;

    double _Complex C[ batch_count*m*n ] = {} ;
    gpu_code( A , m , lda ,
	      B , n , ldb ,
	      C , k , ldc ,
	      batch_count , false ) ;
    printf( "\nC - matrix (Quda result)\n" ) ;
    /* Should be according to cublas routines
     *
     *   C = | 19.0 | 22.0 | 111.0 | 122.0 |
     *       | 43.0 | 50.0 | 151.0 | 166.0 |
     */
    print_matrix( C , m , n , batch_count , false ) ;

    printf( "\n" ) ;
    memset( C , 0.0 , batch_count*m*n*sizeof(double _Complex));
    cpu_code( A , m , lda ,
	      B , n , ldb ,
	      C , k , ldc ,
	      batch_count , false ) ;
    printf( "C - matrix (cpu)\n" ) ;
    print_matrix( C , m , n , batch_count , false ) ;    
  }
  finalize();

  return 0;
}
