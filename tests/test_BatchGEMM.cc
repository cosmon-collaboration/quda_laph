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

//#define VERBOSE

// slow loopy version
void
cpu_code( const double _Complex *A , const int M , const int lda , const int stride_a ,
	  const double _Complex *B , const int N , const int ldb , const int stride_b ,
	  double _Complex *C       , const int K , const int ldc , const int stride_c ,
	  const int batch_count ,
	  const bool is_col_major )
{
  const int a_stride = stride_a ; //M*K ;
  const int b_stride = stride_b ; //K*N ;
  const int c_stride = stride_c ; //M*N ;
  if( is_col_major == true ) {
    for( int p = 0 ; p < batch_count ; p++ ) {
      for( int m = 0 ; m < M ; m++ ) {
	for( int n = 0 ; n < N ; n++ ) {
	  double _Complex c_mnp = 0.0 ;
	  for( int k = 0 ; k < K ; k++ ) {
	    c_mnp += A[m+k*lda+p*a_stride]*B[k+n*ldb+p*b_stride] ;

          #ifdef VERBOSE
	    printf( "In here A[%d] %f x B[%d] %f\n" ,
		    m+k*lda+p*a_stride , creal(A[m+k*lda+p*a_stride]) ,
		    k+n*ldb+p*b_stride , creal(B[k+n*ldb+p*b_stride]) ) ;
	  #endif	    
	  }
	  C[m+n*ldc+p*c_stride] = c_mnp ;

          #ifdef VERBOSE
	  printf( "C[ %d ] %f\n" , m+n*ldc+p*c_stride , c_mnp ) ;
          #endif

	}
      }
            #ifdef VERBOSE
      printf( "\n" ) ;
      #endif

    }
  } else {
    // row-major matrix mul
    for( int p = 0 ; p < batch_count ; p++ ) {
      for( int m = 0 ; m < M ; m++ ) {
	for( int n = 0 ; n < N ; n++ ) {
	  double _Complex c_mnp = 0.0 ;
	  for( int k = 0 ; k < K ; k++ ) {
	    c_mnp += A[k+m*lda+p*a_stride]*B[n+k*ldb+p*b_stride] ;

          #ifdef VERBOSE
	    printf( "In here A[%d] %f x B[%d] %f\n" ,
		    k+m*lda+p*a_stride , creal(A[k+m*lda+p*a_stride]) ,
		    n+k*ldb+p*b_stride , creal(B[n+k*ldb+p*b_stride]) ) ;
	  #endif
	  }
	  C[n+m*ldc+p*c_stride] = c_mnp ;
	  #ifdef VERBOSE
	  printf( "C[ %d ] %f\n" , n+m*ldc+p*c_stride , c_mnp ) ;
          #endif
	}
      }
      #ifdef VERBOSE
      printf( "\n" ) ;
      #endif
    }   
  }
}

QudaBLASParam
getQBLAS( const int m , const int lda , const int stride_a ,
	  const int n , const int ldb , const int stride_b ,
	  const int k , const int ldc , const int stride_c ,
	  const int batch_count , const bool is_col_major )
{
  // do a zgemm
  const __complex__ double alpha = 1.0, beta = 0.0;
  QudaBLASParam cublas_param = newQudaBLASParam();
  cublas_param.blas_type = QUDA_BLAS_GEMM ;
  cublas_param.inv_mat_size = 0 ;
  cublas_param.trans_a = QUDA_BLAS_OP_N;
  cublas_param.trans_b = QUDA_BLAS_OP_N;
  cublas_param.m = m;
  cublas_param.n = n;
  cublas_param.k = k;
  cublas_param.lda = lda;
  cublas_param.ldb = ldb;
  cublas_param.ldc = ldc;
  cublas_param.a_stride = stride_a ;
  cublas_param.b_stride = stride_b ;
  cublas_param.c_stride = stride_c ;
  cublas_param.batch_count = batch_count;
  cublas_param.alpha = (__complex__ double)alpha;  
  cublas_param.beta  = (__complex__ double)beta;
  if( is_col_major ) {
    cublas_param.data_order = QUDA_BLAS_DATAORDER_COL;
  } else {
    cublas_param.data_order = QUDA_BLAS_DATAORDER_ROW;
  }
  cublas_param.data_type = QUDA_BLAS_DATATYPE_Z;
  return cublas_param ;
}

void
eigen_code( const double _Complex *A , const int m , const int lda , const int stride_a ,
	    const double _Complex *B , const int n , const int ldb , const int stride_b ,
	    double _Complex *C       , const int k , const int ldc , const int stride_c ,
	    const int batch_count ,
	    const bool is_col_major )
{
  QudaBLASParam cublas_param = getQBLAS( m , lda , stride_a ,
					 n , ldb , stride_b ,
					 k , ldc , stride_c ,
					 batch_count , is_col_major ) ;

  blas_lapack::generic::stridedBatchGEMM( (double _Complex*)A, (double _Complex*)B, C, cublas_param, QUDA_CPU_FIELD_LOCATION);
}

	  
void
gpu_code( const double _Complex *A , const int m , const int lda , const int stride_a ,
	  const double _Complex *B , const int n , const int ldb , const int stride_b ,
	  double _Complex *C       , const int k , const int ldc , const int stride_c ,
	  const int batch_count ,
	  const bool is_col_major )
{
  size_t Asize , Bsize , Csize ;
  if( is_col_major ) {
    Asize = (stride_a==0)? lda*k : k*m*batch_count ;
    Bsize = (stride_b==0)? ldb*n : k*n*batch_count ;
    Csize = (stride_c==0)? m*n   : m*n*batch_count ;
  } else {
    Asize = (stride_a==0)? lda*m : k*m*batch_count ;
    Bsize = (stride_b==0)? ldb*k : k*n*batch_count ;
    Csize = (stride_c==0)? m*n   : m*n*batch_count ;
  }
  
  printf( "Asize %zu != %zu\n" , Asize , 2*6 ) ;
  printf( "Bsize %zu != %zu\n" , Bsize , 4*6 ) ;
  printf( "Csize %zu != %zu\n" , Csize , 2*4 ) ;
  
  const size_t a_bytes = Asize*sizeof(double _Complex) ;
  const size_t b_bytes = Bsize*sizeof(double _Complex) ;
  const size_t c_bytes = Csize*sizeof(double _Complex) ;
  void *d_A = pool_device_malloc( a_bytes ) ;
  void *d_B = pool_device_malloc( b_bytes ) ;
  void *d_C = pool_device_malloc( c_bytes ) ;
  qudaMemcpy( d_A , A , a_bytes , qudaMemcpyHostToDevice ) ;
  qudaMemcpy( d_B , B , b_bytes , qudaMemcpyHostToDevice ) ;

  QudaBLASParam cublas_param = getQBLAS( m , lda , stride_a ,
					 n , ldb , stride_b ,
					 k , ldc , stride_c ,
					 batch_count , is_col_major ) ;
 
  blas_lapack::native::stridedBatchGEMM(d_A, d_B, d_C, cublas_param, QUDA_CUDA_FIELD_LOCATION);

  qudaMemcpy( C , d_C , c_bytes , qudaMemcpyDeviceToHost ) ;
  
  pool_device_free( d_A ) ;
  pool_device_free( d_B ) ;
  pool_device_free( d_C ) ;
}

static void
print_matrix( const double _Complex *A , const int M , const int lda , const int stride_a , const int num_batch , const bool is_col_order )
{
  if( is_col_order == true ) {
    for( int p = 0 ; p < num_batch ; p++ ) {
      for( int m = 0 ; m < M ; m++ ) {
	for( int n = 0 ; n < lda ; n++ ) {
	  printf( "(%g %g) ", creal( A[m+M*n+p*stride_a] ) , cimag( A[m+M*n+p*stride_a] ) ) ;
	}
	printf( "\n" ) ;
      }
      printf( "===============\n" ) ;
    }
  } else {
    for( int p = 0 ; p < num_batch ; p++ ) {
      for( int m = 0 ; m < M ; m++ ) {
	for( int n = 0 ; n < lda ; n++ ) {
	  printf( "(%g %g) ", creal( A[n+lda*m+p*stride_a] ) , cimag( A[n+lda*m+p*stride_a] ) ) ;
	}
	printf( "\n" ) ;
      }
      printf( "===============\n" ) ;
    }
  }
}

// batched A*B
static void
run_test( const double _Complex *A , const int m , const int lda , const int stride_a ,
	  const double _Complex *B , const int n , const int ldb , const int stride_b ,
                                     const int k , const int ldc , const int stride_c ,
	  const int batch_count ,
	  const bool is_col_major )
{
  printf( "\n***** A matrix *****\n" ) ;
  print_matrix( A , m , k , stride_a , batch_count , is_col_major ) ;
  printf( "\n***** B matrix *****\n" ) ;
  print_matrix( B , k , n , stride_b , batch_count , is_col_major ) ;
  
  double _Complex Ccpu[ n*m*batch_count ] , Cgpu[ n*m*batch_count ] , Ceig[ n*m*batch_count ] , Cint[ n*m*batch_count ] ;
  cpu_code( A , m , lda , stride_a ,
	    B , n , ldb , stride_b ,
	    Ccpu , k , ldc , stride_c ,
	    batch_count , is_col_major ) ;
  printf( "\n***** C matrix (cpu) *****\n" ) ;
  //print_matrix( Ccpu , m , n , stride_c , batch_count , is_col_major ) ;
  print_matrix( Ccpu , 2 , 4, 1, 1, false ) ;
  
  // GPU ----
  gpu_code( A , m , lda , stride_a ,
	    B , n , ldb , stride_b ,
	    Cgpu , k , ldc , stride_c ,
	    batch_count , is_col_major ) ;
  printf( "\n***** C matrix (gpu) *****\n" ) ;
  print_matrix( Cgpu , 2 , 4, 1, 1, false ) ;
  
  eigen_code( A , m , lda , stride_a ,
	      B , n , ldb , stride_b ,
	      Ceig , k , ldc , stride_c ,
	      batch_count , is_col_major ) ;
  printf( "\n***** C matrix (eig lapack) *****\n" ) ;
  print_matrix( Ceig , 2 , 4, 1, 1, false ) ;
  //print_matrix( Cgpu , m , n , stride_c , batch_count , is_col_major ) ;

  QudaBLASParam cublas_param = getQBLAS( m , lda , stride_a ,
					 n , ldb , stride_b ,
					 k , ldc , stride_c ,
					 batch_count , is_col_major ) ;
  blasGEMMQuda( (void*)A , (void*)B , Cint , QUDA_BOOLEAN_FALSE , cublas_param ) ;
  printf( "\n***** C matrix (interface lapack) *****\n" ) ;
  print_matrix( Cint , 2 , 4, 1, 1, false ) ;
  
  double diff = 0. ;
  for( int i = 0 ; i < n*m*batch_count ; i++ ) {
    diff += cabs( Ccpu[i] - Cgpu[i] ) ;
  }
  printf( "\n->Summed diff GPU-CPU %e<-\n\n" , diff ) ;
  diff = 0. ;
  for( int i = 0 ; i < n*m*batch_count ; i++ ) {
    diff += cabs( Ccpu[i] - Ceig[i] ) ;
  }
  printf( "\n->Summed diff Interface-CPU %e<-\n\n" , diff ) ;
  diff = 0. ;
  for( int i = 0 ; i < n*m*batch_count ; i++ ) {
    diff += cabs( Ccpu[i] - Cint[i] ) ;
  }
  printf( "\n->Summed diff Eigen-CPU %e<-\n\n" , diff ) ;
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

#if 0
  const int m = 2 ;
  const int n = 2 ;
  const int k = 2 ;
  const int lda = 2 ;
  const int ldb = 2 ;
  const int ldc = 2 ;
  const int batch_count = 2 ;

  printf( "Column major odering\n" ) ;
  {
    // do a batched blas call example from
    // CUDALibrarySamples/cuBLAS/Level-3/gemmStridedBatched/cublas_gemmStridedBatched_example.cu

    // should be of size m*lda*batch_count
    const double _Complex A[ 2*2*2 ] =		\
      { 1.0, 3.0, 2.0, 4.0,
	5.0, 7.0, 6.0, 8.0 } ;
    // should be of size n*ldb*batch_count? Does it need to be???
    const double _Complex B[ 2*2*2 ] =		\
      { 5.0, 7.0 , 6.0 , 8.0,
	9.0, 11.0, 10.0, 12.0 } ;
    run_test( A , m , lda , m*k ,
	      B , n , ldb , n*k ,
	      k , ldc , m*n , batch_count , true ) ;
  }
  printf( "\nRow major ordering\n" ) ;
  {
    // OK that was all column-major now redo it being row-major in Quda and do a CPU version
    const double _Complex A[ 2*2*2 ] =		\
      { 1.0, 2.0, 3.0, 4.0,
	5.0, 6.0, 7.0, 8.0 } ;
    // should be of size n*ldb*batch_count? Does it need to be???
    const double _Complex B[ 2*2*2 ] =		\
      { 5.0, 6.0 , 7.0 , 8.0,
	9.0, 10.0, 11.0, 12.0 } ;

    run_test( A , m , lda , m*k ,
	      B , n , ldb , n*k ,
	      k , ldc , m*n , batch_count , false ) ;    
  }
#endif

#if 1
  printf( "\n**************************\n\n" ) ;
  {
    // 2mom 6 sites
    // 0 1 2 3 4 5
    // 6 7 8 9 10 11
    const double _Complex A[ 2*6 ] =		\
      { 0, 6, 1, 7, 2, 8,
	3, 9, 4, 10, 5, 11 } ;
    print_matrix( A , 2 , 6 , 0 , 1 , true ) ;

    // 6 sites 4 dil
    const double _Complex B[ 6*4 ] = \
      { 0,   1,  2,  3,  4,  5,
	6,   7,  8,  9, 10, 11,
	12, 13, 14, 15, 16, 17,
	18, 19, 20, 21, 22, 23 } ;
    print_matrix( B , 4 , 6 , 0 , 1 , true ) ;

    // full ZGEMM
    run_test( A , 2 , 2 , 0 ,
	      B , 4 , 6 , 0 ,
	          6 , 2 , 0 ,
	      1 , true ) ;

    // Batching in A
    run_test( A , 1 , 2 , 1 ,
	      B , 4 , 6 , 0 ,
	          6 , 2 , 1 ,
	      2 , true ) ;

    // Batching in B
    run_test( A , 2 , 2 , 0 ,
	      B , 1 , 6 , 6 ,
	          6 , 2 , 2 ,
	      4 , true ) ;
    
    // Batching in B by factor 2
    run_test( A , 2 , 2 , 0 ,
	      B , 2 , 6 , 12 ,
	          6 , 2 , 4 ,
	      2 , true ) ;
  }
#endif

#if 1
  printf( "\n**************************\n\n" ) ;
  {
    // 2mom 6 sites
    const double _Complex A[ 2*6 ] = \
      { 0 , 1 , 2 , 3 , 4  , 5 ,
	6 , 7 , 8 , 9 , 10 , 11 } ;
    // 6 sites 4 dil
    const double _Complex B[ 6*4 ] = \
      { 0 , 6  , 12 , 18,
	1 , 7  , 13 , 19,
	2 , 8  , 14 , 20,
	3 , 9  , 15 , 21,
	4 , 10 , 16 , 22,
	5 , 11 , 17 , 23 } ;
    // full ZGEMM
    run_test( A , 2 , 6 , 0 ,
	      B , 4 , 4 , 0 ,
	          6 , 4 , 0 ,
	      1 , false ) ; 
    // Batching in A
    run_test( A , 1 , 6 , 6 ,
	      B , 4 , 4 , 0 ,
	          6 , 4 , 4 ,
	      2 , false ) ; 

    // Batching in B
    run_test( A , 2 , 6 , 0 ,
	      B , 1 , 4 , 1 ,
	          6 , 4 , 1 ,
	      4 , false ) ;

    // partial batching in B
    run_test( A , 2 , 6 , 0 ,
	      B , 2 , 4 , 2 ,
	          6 , 4 , 2 ,
	      2 , false ) ; 
  }
#endif
  
  
  finalize();

  return 0;
}
