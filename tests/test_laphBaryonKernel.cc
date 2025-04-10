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

static void
evprod( const double _Complex *coeffs ,
	void **host_evec ,
	const size_t n ,
	const size_t nEv ,
	const size_t nsites ,
	double _Complex *q )
{
#pragma omp parallel for
  for( size_t i = 0 ; i < nsites ; i++ ) {
    for( size_t dil = 0 ; dil < n ; dil++ ) {
      for( size_t c = 0 ; c < 3 ; c++ ) {
	double _Complex sum = 0.0 ;
	for( size_t ev = 0 ; ev < nEv ; ev++ ) {
	  double _Complex *pt = (double _Complex*)host_evec[ev]+c+3*i ; 
	  sum += (*pt) * coeffs[ev+dil*nEv] ;
	}
	q[c+3*(i+nsites*dil)] = sum ;
      }
    }
  }
}

static inline
void diq( const double _Complex *q1 ,
	  const double _Complex *q2 ,
	  double _Complex Diq[3] )
{
  Diq[0] = +q1[1]*q2[2] - q1[2]*q2[1] ;
  Diq[1] = -q1[0]*q2[2] + q1[2]*q2[0] ;
  Diq[2] = +q1[0]*q2[1] - q1[1]*q2[0] ;
}
	
void cpu_code( const int n1, const int n2, const int n3, const int nMom,
	       const double _Complex *host_coeffs1, 
	       const double _Complex *host_coeffs2, 
	       const double _Complex *host_coeffs3,
	       const double _Complex *host_mom, 
	       const int nEv,
	       void **host_evec, 
	       double _Complex *return_arr,
	       const int blockSizeMomProj,
	       const int X[4] )
{
  const size_t nsites = X[0]*X[1]*X[2] ;
  double _Complex q1[ n1*nsites*3 ] = {} , q2[ n2*nsites*3 ] = {} , q3[ n3*nsites*3 ] = {} ;
  evprod( host_coeffs1 , host_evec , n1 , nEv , nsites , q1 ) ;
  evprod( host_coeffs2 , host_evec , n2 , nEv , nsites , q2 ) ;
  evprod( host_coeffs3 , host_evec , n3 , nEv , nsites , q3 ) ;

  // ok and then it is just the baryon triplet contraction with q1,q2,q3
  //#pragma omp parallel for collapse(3) reduction(+:return_arr[:n1*n2*n3*nMom])
  for( size_t i = 0 ; i < nsites ; i++ ) {
    for( int dil1 = 0 ; dil1 < n1 ; dil1++ ) {
      for( int dil2 = 0 ; dil2 < n2 ; dil2++ ) {
	double _Complex Diq[3] ;
	diq( (const double _Complex*)q1 + 3*(i+nsites*dil1) ,
	     (const double _Complex*)q2 + 3*(i+nsites*dil2) ,
	     Diq ) ;
	for( int dil3 = 0 ; dil3 < n3 ; dil3++ ) {
	  const double _Complex Con = \
	    Diq[0]*q3[0+3*(i+nsites*dil3)] +
	    Diq[1]*q3[1+3*(i+nsites*dil3)] +
	    Diq[2]*q3[2+3*(i+nsites*dil3)] ;
	  for( int p = 0 ; p < nMom ; p++ ) {
	    return_arr[dil3+n3*(dil2+n2*(dil1+n1*p))] += host_mom[i+nsites*p]*Con ; 
	  }
	}
      }
    }
  }
}

static void
alamode( const int n1, const int n2, const int n3, const int nMom,
	 const double _Complex *host_coeffs1, 
	 const double _Complex *host_coeffs2, 
	 const double _Complex *host_coeffs3,
	 const double _Complex *host_mom, 
	 const int nEv,
	 void **host_evec,
	 QudaInvertParam inv_param,
	 double _Complex *return_arr,
	 const int blockSizeMomProj,
	 const int X[4] )
{
  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_TOTAL);

  // checks and balances
  //if (sizeof(Complex) != sizeof(double _Complex)) {
  //  errorQuda("Irreconcilable difference between interface and internal complex number conventions");
  // }
  if( (n1*n2*n3)%blockSizeMomProj != 0 ) {
    errorQuda( "Block size mom proj needs to divide %d %d\n" , n1*n2*n3 , blockSizeMomProj ) ;
  }
  // Allocate device memory for evecs. This is done to ensure a contiguous
  // this is a double store of evecs, which is bad  
  const int nSites = X[0]*X[1]*X[2];

  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_INIT);
  lat_dim_t x = { X[0] , X[1] , X[2] , X[3] } ;
  ColorSpinorParam cpu_evec_param(host_evec, inv_param, x, false, QUDA_CPU_FIELD_LOCATION);
  cpu_evec_param.nSpin = 1;
  std::vector<ColorSpinorField> evec(nEv);
  for (int iEv=0; iEv<nEv; ++iEv) {
    cpu_evec_param.v = host_evec[iEv];
    evec[iEv] = ColorSpinorField(cpu_evec_param) ;
  }
  // Create device evecs
  ColorSpinorParam cuda_evec_param(cpu_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_evec_param.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  std::vector<ColorSpinorField*> quda_evec ;
  for (int i=0; i<nEv; i++) {
    quda_evec.push_back( ColorSpinorField::Create(cuda_evec_param) );
    *quda_evec[i] = evec[i] ; // load here because fuck it why not
  }
  // Create q1
  ColorSpinorParam cuda_q1_param(cuda_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_q1_param.create = QUDA_ZERO_FIELD_CREATE;
  std::vector<std::complex<double>> coeffs1(n1*nEv) ;
  std::vector<ColorSpinorField*> quda_q1 ;
  for(int i=0; i<n1; i++) {
    quda_q1.push_back(ColorSpinorField::Create(cuda_q1_param));
    for( int j = 0 ; j < nEv ; j++ ) coeffs1[j*n1+i] = host_coeffs1[j+i*nEv] ;
  }
  // Create q2
  ColorSpinorParam cuda_q2_param(cuda_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_q2_param.create = QUDA_ZERO_FIELD_CREATE;
  std::vector<std::complex<double>> coeffs2(n2*nEv) ;
  std::vector<ColorSpinorField*> quda_q2 ;
  for(int i=0; i<n2; i++) {
    quda_q2.push_back(ColorSpinorField::Create(cuda_q2_param));
    for( int j = 0 ; j < nEv ; j++ ) coeffs2[j*n2+i] = host_coeffs2[j+i*nEv] ;
  }
  // create q3
  ColorSpinorParam cuda_q3_param(cuda_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_q3_param.create = QUDA_ZERO_FIELD_CREATE;
  std::vector<std::complex<double>> coeffs3(n3*nEv) ;
  std::vector<ColorSpinorField*> quda_q3 ;
  for(int i=0; i<n3; i++) {
    quda_q3.push_back(ColorSpinorField::Create(cuda_q3_param));
    for( int j = 0 ; j < nEv ; j++ ) coeffs3[j*n3+i] = host_coeffs3[j+i*nEv] ;
  }

  // device temporaries, momentum and return buffers
  const size_t data_tmp_bytes = blockSizeMomProj*X[0]*X[1]*X[2]*2*quda_q3[0]->Precision();
  void *d_tmp = pool_device_malloc(data_tmp_bytes);

  const size_t data_ret_bytes = nMom*n1*n2*n3*2*quda_q3[0]->Precision();
  void *d_ret = pool_device_malloc(data_ret_bytes);

  const size_t data_mom_bytes = nMom*nSites*2*quda_q3[0]->Precision();
  void *d_mom = pool_device_malloc(data_mom_bytes);
  //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_INIT);  
  
  // Copy host_mom data to device
  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_H2D);
  qudaMemcpy(d_mom, host_mom, data_mom_bytes, qudaMemcpyHostToDevice);  
  //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_H2D);

  // Perfrom the caxpy to compute all q-vectors
  //getProfileAccumulateEvecs().TPSTART(QUDA_PROFILE_COMPUTE);
  quda::blas::legacy::caxpy(coeffs1.data(), quda_evec , quda_q1 ) ;
  quda::blas::legacy::caxpy(coeffs2.data(), quda_evec , quda_q2 ) ;
  quda::blas::legacy::caxpy(coeffs3.data(), quda_evec , quda_q3 ) ;
  //getProfileAccumulateEvecs().TPSTOP(QUDA_PROFILE_COMPUTE);

  // evecs irrelevant and can be purged here I gues if we are really desperate for space

  // Create device diquark vector
  ColorSpinorParam cuda_diq_param( cuda_evec_param , inv_param , QUDA_CUDA_FIELD_LOCATION ) ;
  ColorSpinorField quda_diq( cuda_diq_param ) ;

  // usual momentum contraction
  const __complex__ double alpha = 1.0, beta = 0.0;
  QudaBLASParam cublas_param_mom_sum = newQudaBLASParam();
  cublas_param_mom_sum.trans_a = QUDA_BLAS_OP_N;
  cublas_param_mom_sum.trans_b = QUDA_BLAS_OP_T;
  cublas_param_mom_sum.m = nMom;
  cublas_param_mom_sum.k = nSites;
  cublas_param_mom_sum.lda = nSites;
  cublas_param_mom_sum.ldb = nSites;
  cublas_param_mom_sum.ldc = n1*n2*n3;
  cublas_param_mom_sum.batch_count = 1;
  cublas_param_mom_sum.alpha = (__complex__ double)alpha;  
  cublas_param_mom_sum.beta  = (__complex__ double)beta;
  cublas_param_mom_sum.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param_mom_sum.data_type = QUDA_BLAS_DATATYPE_Z;
  cublas_param_mom_sum.blas_type = QUDA_BLAS_GEMM ;
  
  int nInBlock = 0;
  for( int dil1=0; dil1<n1; dil1++ ) {
    for( int dil2=0; dil2<n2; dil2++ ) {

      //getProfileColorCross().TPSTART(QUDA_PROFILE_COMPUTE);
      colorCrossQuda(*quda_q1[dil1], *quda_q2[dil2], quda_diq);
      //getProfileColorCross().TPSTOP(QUDA_PROFILE_COMPUTE);
      
      for (int dil3=0; dil3<n3; dil3++) {

	//getProfileColorContract().TPSTART(QUDA_PROFILE_COMPUTE);	
	colorContractQuda(quda_diq, *quda_q3[dil3], (std::complex<double>*)d_tmp + nSites*nInBlock);
	//getProfileColorContract().TPSTOP(QUDA_PROFILE_COMPUTE);
	nInBlock++;

	if (nInBlock == blockSizeMomProj ) {
	  // To gauge how to block the calls to remove launch latency.
	  printfQuda("dil1 = %d, dil2 = %d, dil3 = %d, nInBlock = %d\n", dil1, dil2, dil3, nInBlock);
	  cublas_param_mom_sum.n = nInBlock;
	  cublas_param_mom_sum.c_offset = (dil1*n2 + dil2)*n3 + dil3 - nInBlock + 1;
          //cublas_param_mom_sum.b_stride = nInBlock * nSites;
          //cublas_param_mom_sum.c_stride = nInBlock * nSites;
	  //getProfileBLAS().TPSTART(QUDA_PROFILE_COMPUTE);	  
	  blas_lapack::native::stridedBatchGEMM(d_mom, d_tmp, d_ret, cublas_param_mom_sum, QUDA_CUDA_FIELD_LOCATION);
	  //getProfileBLAS().TPSTOP(QUDA_PROFILE_COMPUTE);	  
	  nInBlock = 0;
	}
      }
    }
  }

  // Copy return array back to host
  // getProfileBaryonKernel().TPSTART(QUDA_PROFILE_D2H);
  qudaMemcpy(return_arr, d_ret, data_ret_bytes, qudaMemcpyDeviceToHost);  
  //getProfileBaryonKernel().TPSTOP(QUDA_PROFILE_D2H);
  
  // Clean up memory allocations
  //getProfileBaryonKernel().TPSTART(QUDA_PROFILE_FREE);
  for (int i=0; i<n1; i++ ) delete quda_q1[i];
  for (int i=0; i<n2; i++ ) delete quda_q2[i];
  for (int i=0; i<n3; i++ ) delete quda_q3[i];
  for (int i=0; i<nEv; i++) delete quda_evec[i];
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
  assert( global == 1 ) ;
  
  // call rephase here
  const int Nev = 8 ;
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);

  std::cout<<"Constant Eigvecs"<<std::endl ;
  set_constant( laphEigvecs ) ;
  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 3 , blockSizeMomProj = 4 , n1 = 2 , n2 = 2 , n3 = 2 ;
  const int X[4] = {
    LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int nsites = X[0]*X[1]*X[2] ;
  assert( X[3] == 1 ) ;
  
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
  
  // what is this? exponentiated mom? Which index runs fastest? Fuck knows
  double _Complex host_mom[ nsites*nmom ] = {} ;
    // host_mom should be complex
  for( size_t p = 0 ; p < nmom ; p++ ) {
    const size_t mom[3] = { p+1, p+2, p+3 } ;
    for( size_t i = 0 ; i < (size_t)nsites ; i++ ) {
      const size_t x = i%X[0] ;
      const size_t y = (i/X[0])%X[1] ;
      const size_t z = (i/(X[0]*X[1]))%X[2] ;
      const double arg = 2*M_PI*( mom[0]*x/(double)X[0] +
				  mom[1]*y/(double)X[1] +
				  mom[2]*z/(double)X[2] ) ;
      double s = 1 ,c = 1 ;
      sincos( arg , &s , &c ) ;
      host_mom[ i + nsites*p ] = c+I*s ;
    }
  }
  
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
  
  double _Complex retGPU[ n1*n2*n3*nmom ] = {} ;
  /*
  laphBaryonKernel( n1, n2, n3,
		    nmom,
		    coeffs1 ,
		    coeffs2 ,
		    coeffs3 ,
		    host_mom ,
		    Nev ,
		    evList.data(),
		    retGPU,
		    blockSizeMomProj,
		    X ) ;
  */
  StopWatch GPU ;
  GPU.start() ;
  memset( retGPU , 0.0 , n1*n2*n3*sizeof(double _Complex)) ;
  alamode( n1, n2, n3,
	   nmom,
	   coeffs1 ,
	   coeffs2 ,
	   coeffs3 ,
	   host_mom ,
	   Nev ,
	   evList.data(),
	   inv_param ,
	   retGPU,
	   blockSizeMomProj,
	   X ) ;
  GPU.stop() ;
  printLaph(make_strf("\nGPU baryonkernel in = %g seconds\n", GPU.getTimeInSeconds()));

  StopWatch CPU ;
  CPU.start() ;
  double _Complex retCPU[ n1*n2*n3*nmom ] = {} ;
  cpu_code( n1 , n2 , n3 ,
	    nmom,
	    coeffs1, 
	    coeffs2, 
	    coeffs3,
	    host_mom, 
	    Nev,
	    evList.data(),
	    retCPU,
	    blockSizeMomProj,
	    X ) ;
  CPU.stop() ;
  printLaph(make_strf("\nCPU baryonkernel in = %g seconds\n", CPU.getTimeInSeconds()));
  
  // test outputs
  printf( "CPU == GPU\n" ) ;
  for( size_t p = 0 ; p < nmom ; p++ ) {
    double sum = 0 ;
    for( size_t i = 0 ; i < n1*n2*n3 ; i++ ) {
      const size_t idx = i + (n1*n2*n3)*p ;
      sum += cabs( retCPU[idx] - retGPU[idx] ) ;
      #ifdef VERBOSE_COMPARISON
      printf( " (%f %f) == (%f %f)\n" ,
	      creal(retCPU[idx]) , cimag(retCPU[idx]) ,
	      creal(retGPU[idx]) , cimag(retGPU[idx]) ) ;
      #endif
    }
    std::cout<<"Summed diff p="<<p<<" "<<sum<<std::endl ;
  }
  
  finalize();

  return 0;
}
