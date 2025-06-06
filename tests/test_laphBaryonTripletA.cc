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
#define GPU_STRESS

// cpu color cross
static void
cpuColorCross( void *A , void *B , void *result , const int X[4] )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ;
  std::complex<double> *ptA = (std::complex<double>*)A ;
  std::complex<double> *ptB = (std::complex<double>*)B ;
  std::complex<double> *ptC = (std::complex<double>*)result ;
#pragma omp parallel for
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
  #pragma omp parallel for
  for( size_t i = 0 ; i < (size_t)Nsites ; i++ ) {
    // simple inner product over color A_{c}B_{c}
    ptC[i]  = ptA[0+3*i]*(ptB[0+3*i]) ;
    ptC[i] += ptA[1+3*i]*(ptB[1+3*i]) ;
    ptC[i] += ptA[2+3*i]*(ptB[2+3*i]) ;
  }
}

// slow loopy version
static void
cpu_code( const int nMom,
	  const int nEv,
	  const int blockSizeMomProj,
	  void **host_evec, 
	  const double _Complex *host_mom,
	  double _Complex *return_arr,
	  const int X[4])
{  
  // spatial only
  const size_t nSp   = X[0]*X[1]*X[2] ;
  const size_t nEvC3 = nEv*(nEv-1)*(nEv-2)/6 ;

  LattField Diq( FieldSiteType::ColorVector);
  LattField tmp( FieldSiteType::Complex);

  size_t idx = 0 ;
  for( int aEv = 0 ; aEv < nEv ; aEv++ ) {
    for( int bEv = aEv+1 ; bEv < nEv ; bEv++ ) {
      cpuColorCross( host_evec[aEv] , host_evec[bEv] , (void*)Diq.getDataPtr() , X ) ;
      for( int cEv = bEv+1 ; cEv < nEv ; cEv++ ) {
	// color contract
	cpuColorContract( (void*)Diq.getDataPtr() , host_evec[cEv] , (void*)tmp.getDataPtr() , X ) ;

	// Matrix mul here mom*d_tmp -> d_ret
	for( int p = 0 ; p < nMom ; p++ ) {
	  const double _Complex *p2 = (const double _Complex*)host_mom+nSp*p ;
	  for( int T = 0 ; T < X[3] ; T++ ) {
	    const double _Complex *p1 = (const double _Complex*)tmp.getDataPtr() + nSp*T ;
	    double _Complex sum = 0.0 ;
            #pragma omp parallel for reduction(+:sum)
	    for( size_t i = 0 ; i < nSp ; i++ ) {
	      sum += p2[i]*p1[i] ;
	    }
	    return_arr[ T+X[3]*(idx + nEvC3*p) ] = sum ;
	  }
	}
	idx++ ;
      }
    }
  }
}

// slow loopy version
static void
cpu_codev2( const int nMom,
	    const int nEv,
	    const int blockSizeMomProj,
	    void **host_evec, 
	    const double _Complex *host_mom,
	    double _Complex *return_arr,
	    const int X[4])
{  
  // spatial only
  const size_t nSp    = X[0]*X[1]*X[2] ;  
  const size_t nEvC3 = nEv*(nEv-1)*(nEv-2)/6 ;
  memset( return_arr , 0.0 , X[3]*nEvC3*nMom*sizeof( double _Complex ) ) ;
  LattField Diq( FieldSiteType::ColorVector);
  size_t idx = 0 ;
  for( int aEv = 0 ; aEv < nEv ; aEv++ ) {
    for( int bEv = aEv+1 ; bEv < nEv ; bEv++ ) {
      cpuColorCross( host_evec[aEv] , host_evec[bEv] , (void*)Diq.getDataPtr() , X ) ;
      for( int cEv = bEv+1 ; cEv < nEv ; cEv++ ) {

	// inline contraction and momproj
	double _Complex sum[nMom*X[3]] ;
	memset( sum , 0 , nMom*X[3]*sizeof(double _Complex) ) ;
	const double _Complex *ptA = (const double _Complex*)Diq.getDataPtr() ;
	const double _Complex *ptB = (const double _Complex*)host_evec[cEv] ;

	#pragma omp parallel
	{
	  const double _Complex *p2 = (const double _Complex*)host_mom  ;
          #pragma omp for reduction(+:sum[:nMom*X[3]]) collapse(2)
	  for( int T = 0 ; T < X[3] ; T++ ) {
	    for( size_t i = 0 ; i < (size_t)nSp ; i++ ) {
	      const size_t idx = i + T*nSp ;
	      // is the color contraction
	      const double _Complex Con = \
		ptA[0+3*idx]*(ptB[0+3*idx]) +\
		ptA[1+3*idx]*(ptB[1+3*idx]) +\
		ptA[2+3*idx]*(ptB[2+3*idx]) ;
	      for( int p = 0 ; p < nMom ; p++ ) {
		sum[T+X[3]*p] += p2[i+p*nSp]*Con ;
	      }
	    }
	  }
	}
	// copy back
	for( int p = 0 ; p < nMom ; p++ ) {
	  for( int T = 0 ; T < X[3] ; T++ ) {
	    return_arr[ T + X[3]*(idx + nEvC3*p) ] = sum[T + X[3]*p] ;
	  }
	}
	idx++ ;
      }
    }
  }
}

// slow loopy version
static void
cpu_codev3( const int nMom,
	    const int nEv,
	    const int blockSizeMomProj,
	    void **host_evec, 
	    const double _Complex *host_mom,
	    double _Complex *return_arr,
	    const int X[4])
{
  // spatial only
  const size_t nSp   = X[0]*X[1]*X[2] ;  
  const size_t nEvC3 = nEv*(nEv-1)*(nEv-2)/6 ;
  memset( return_arr , 0.0 , X[3]*nEvC3*nMom*sizeof( double _Complex ) ) ;
  // loop over sites is far more efficient as we are only opening a parallel region once
#pragma omp parallel for reduction(+:return_arr[:X[3]*nEvC3*nMom]) collapse(2)
  for( int T = 0 ; T < X[3] ; T++ ) {
    for( size_t i = 0 ; i < nSp ; i++ ) {
      // precache these 
      const double _Complex *pt[ nEv ] ;
      for( int ev = 0 ; ev < nEv ; ev++ ) {
	pt[ ev ] = (const double _Complex*)host_evec[ev] + 3*(i+T*nSp) ;
      }
      size_t idx = 0 ;
      for( int aEv = 0 ; aEv < nEv ; aEv++ ) {
	for( int bEv = aEv+1 ; bEv < nEv ; bEv++ ) {
	  const double _Complex Diq[3] = {
	    +pt[aEv][1]*pt[bEv][2] - pt[aEv][2]*pt[bEv][1] ,
	    -pt[aEv][0]*pt[bEv][2] + pt[aEv][2]*pt[bEv][0] ,	  
	    +pt[aEv][0]*pt[bEv][1] - pt[aEv][1]*pt[bEv][0] } ;
	  for( int cEv = bEv+1 ; cEv < nEv ; cEv++ ) {
	    const double _Complex *p2 = (const double _Complex*)host_mom+i ;
	    // color contraction
	    const double _Complex Con = Diq[0]*(pt[cEv][0])+Diq[1]*(pt[cEv][1])+Diq[2]*(pt[cEv][2]) ;
	    for( int p = 0 ; p < nMom ; p++ ) {
	      return_arr[T+X[3]*(idx+nEvC3*p)] += p2[p*nSp]*Con ;
	    }
	    idx++ ;
	  }
	}
      }
    }
  }
}

// so I tied an onion to my belt, which was the style at the time
void alamode( const int nMom,
	      const int nEv,
	      const int blockSizeMomProj,
	      void **host_evec, 
	      const double _Complex *host_mom,
	      QudaInvertParam inv_param,
	      double _Complex *return_arr,
	      const int X[4])
{
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_TOTAL);
  
  const size_t nSp    = X[0]*X[1]*X[2];
  const size_t nSites = nSp*X[3] ;
  const size_t nEvChoose3 = nEv*(nEv-1)/2*(nEv-2)/3;
  
  // appropriate checks and balances
  //if (sizeof(Complex) != sizeof(double _Complex)) {
  //  errorQuda("Irreconcilable difference between interface and internal complex number conventions");
  //}
  //if( nEvChoose3%blockSizeMomProj != 0 ) {
  //  errorQuda("Block size mom proj needs to divide %zu %d", nEvChoose3 , blockSizeMomProj);
  //}
  if( (size_t)blockSizeMomProj > nEvChoose3 ) { errorQuda("blockSizeMomProj too big %d > %zu" , blockSizeMomProj , nEvChoose3 ) ; }
  
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
  // chuck all the evecs on the GPU
  ColorSpinorParam cuda_evec_param(cpu_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  cuda_evec_param.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);
  std::vector<ColorSpinorField> quda_evec(nEv);
  for (int i=0; i<nEv; i++) {
    quda_evec[i] = ColorSpinorField(cuda_evec_param) ;
    quda_evec[i] = evec[i] ; // CPU -> GPU
  }
  
  // Create device diquark vector
  ColorSpinorParam cuda_diq_param(cpu_evec_param,inv_param,QUDA_CUDA_FIELD_LOCATION);
  ColorSpinorField quda_diq(cuda_diq_param) ;
  
  // Device side temp array (complBuf in chroma_laph)
  const size_t data_tmp_bytes = blockSizeMomProj*nSites*2*quda_evec[0].Precision();
  const size_t data_ret_bytes = nEvChoose3*nMom*X[3]*2*quda_evec[0].Precision();
  const size_t data_mom_bytes = nMom*nSp*2*quda_evec[0].Precision();
  void *d_tmp = pool_device_malloc(data_tmp_bytes);
  void *d_ret = pool_device_malloc(data_ret_bytes);
  void *d_mom = pool_device_malloc(data_mom_bytes);
  if( getVerbosity() >= QUDA_SUMMARIZE ) {
    const size_t OneGB = 1024*1024*1024;
    const size_t total_bytes = data_tmp_bytes + data_ret_bytes + data_mom_bytes ;
    printfQuda("d_tmp %fGB | d_ret %fGB | d_mom %fGB | total = %fGB\n",
	       (double)data_tmp_bytes/(OneGB), (double)data_ret_bytes/(OneGB),
	       (double)data_mom_bytes/(OneGB), (double)total_bytes/(OneGB)); 
  }
  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_INIT);

  // Copy host data to device
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_H2D);
  qudaMemcpy(d_mom, host_mom, data_mom_bytes, qudaMemcpyHostToDevice);  
  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_H2D);

  // idea here like always is to do several ev-blocks at once in a zgemm
  QudaBLASParam cublas_param_mom_sum = newQudaBLASParam();
  cublas_param_mom_sum.trans_a = QUDA_BLAS_OP_N;
  cublas_param_mom_sum.trans_b = QUDA_BLAS_OP_T;
  cublas_param_mom_sum.m = nMom ;
  cublas_param_mom_sum.k = nSp ;
  cublas_param_mom_sum.n = X[3] ;
  cublas_param_mom_sum.lda = nSp ;
  cublas_param_mom_sum.ldb = nSp ;
  cublas_param_mom_sum.ldc = X[3]*nEvChoose3;
  cublas_param_mom_sum.a_stride = 0 ;
  cublas_param_mom_sum.b_stride = nSites ;
  cublas_param_mom_sum.c_stride = X[3] ;
  cublas_param_mom_sum.batch_count = blockSizeMomProj;
  cublas_param_mom_sum.alpha = 1. ; cublas_param_mom_sum.beta = 0. ;
  cublas_param_mom_sum.data_order = QUDA_BLAS_DATAORDER_ROW;
  cublas_param_mom_sum.data_type = QUDA_BLAS_DATATYPE_Z;

  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_TOTAL);

  int nInBlock = 0, blockStart = 0;
  for (int aEv=0; aEv<nEv; aEv++) {
    for (int bEv=aEv+1; bEv<nEv; bEv++) {
      //getProfileColorCross().TPSTART(QUDA_PROFILE_COMPUTE);
      colorCrossQuda(quda_evec[aEv], quda_evec[bEv], quda_diq);
      //getProfileColorCross().TPSTOP(QUDA_PROFILE_COMPUTE);
      for (int cEv=bEv+1; cEv<nEv; cEv++) {
	//getProfileColorContract().TPSTART(QUDA_PROFILE_COMPUTE);
	colorContractQuda(quda_diq, quda_evec[cEv],(std::complex<double>*)d_tmp + nSites*nInBlock);
	//getProfileColorContract().TPSTOP(QUDA_PROFILE_COMPUTE);
	nInBlock++;
	// eh todo this can be cleaned up
	if (nInBlock == blockSizeMomProj) {
	  //getProfileBLAS().TPSTART(QUDA_PROFILE_COMPUTE);  
	  blas_lapack::native::stridedBatchGEMM(d_mom, d_tmp,
						(std::complex<double>*)d_ret+X[3]*blockStart,
						cublas_param_mom_sum,
						QUDA_CUDA_FIELD_LOCATION);
	  //getProfileBLAS().TPSTOP(QUDA_PROFILE_COMPUTE);
	  blockStart += nInBlock;
	  nInBlock = 0;
	}
      }
    }
  }
  // overspill, code is more efficient if you avoid this
  if( nInBlock > 0 ) {
    cublas_param_mom_sum.batch_count = nInBlock;
    blas_lapack::native::stridedBatchGEMM(d_mom, d_tmp,
					  (std::complex<double>*)d_ret+X[3]*blockStart,
					  cublas_param_mom_sum,
					  QUDA_CUDA_FIELD_LOCATION);
  }

  // Copy return array back to host
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_TOTAL); 
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_D2H);
  qudaMemcpy(return_arr, d_ret, data_ret_bytes, qudaMemcpyDeviceToHost);  
  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_D2H);
  
  // Clean up memory allocations
  //getProfileBaryonKernelModeTripletsA().TPSTART(QUDA_PROFILE_FREE);

  pool_device_free(d_tmp);
  pool_device_free(d_mom);
  pool_device_free(d_ret);

  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_FREE);
  //getProfileBaryonKernelModeTripletsA().TPSTOP(QUDA_PROFILE_TOTAL);
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

  const int max_threads = omp_get_max_threads() ;
  std::cout<< "Max threads here" << max_threads << std::endl ;

  int global = 1 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
#endif
  setVerbosityQuda(QUDA_VERBOSE, "#" , stdout ) ;
  
#ifdef GPU_STRESS
  const int Nev = 96 ;
#else
  const int Nev = 64 ;
#endif
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);

  std::cout<<"Constant Eigvecs :: "<< Nev <<std::endl ;
  set_constant( laphEigvecs ) ;

  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 40 ;
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

  const size_t nEvChoose3 = Nev*(Nev-1)*(Nev-2)/6;
  double _Complex *retGPU = (double _Complex*)calloc( X[3]*nmom*nEvChoose3 , sizeof( double _Complex  ) );

#ifdef GPU_STRESS
  for( int blockSizeMomProj = 2 ; blockSizeMomProj < 8192 ; blockSizeMomProj *= 2 ) {
    std::cout<<"nmom "<<nmom<<" | block "<<blockSizeMomProj<<std::endl ;
    memset( retGPU , 0.0 , X[3]*nmom*nEvChoose3*sizeof( double _Complex )) ;
#else
    const int blockSizeMomProj = 64 ;
#endif
    StopWatch gpu ;
    gpu.start() ;
    //alamode(
    laphBaryonKernelComputeModeTripletA(
					nmom,
					Nev,
					blockSizeMomProj,
					evList.data() ,
					host_mom ,
					inv_param,
					retGPU,
					X ) ;
    gpu.stop() ;
    const double GPUtime = gpu.getTimeInSeconds() ;
    printLaph(make_strf("\nGPU modetripletA in = %g seconds\n", GPUtime )) ;
#ifdef GPU_STRESS
  }
#else
  double _Complex *retCPU = (double _Complex*)calloc( X[3]*nmom*nEvChoose3 , sizeof( double _Complex  ) );
  StopWatch cpu ;
  cpu.start() ;
  //cpu_code( nmom, Nev, blockSizeMomProj, evList.data() , host_mom, retCPU, X ) ;
  cpu_codev2( nmom, Nev, blockSizeMomProj, evList.data() , host_mom, retCPU, X ) ;
  //cpu_codev3( nmom, Nev, blockSizeMomProj, evList.data() , host_mom, retCPU, X ) ;
  cpu.stop() ;
  const double CPUtime = cpu.getTimeInSeconds() ;
  printLaph(make_strf("\nCPUv3 modetripletA in = %g seconds\n", CPUtime ));
  
  printf( "\n*************************************\n" ) ;
  printf( "-----> GPU speedup factor %gx\n" , CPUtime/GPUtime ) ;
  printf( "*************************************\n\n" ) ;
  
  // test outputs
  printf( "CPU == GPU\n" ) ;
  for( size_t p = 0 ; p < nmom ; p++ ) {
    double sum = 0 ;
    for( int T = 0 ; T < X[3] ; T++ ) {
      for( size_t i = 0 ; i < nEvChoose3 ; i++ ) {
	const size_t idx = T + X[3]*(i + nEvChoose3*p) ;
	sum += cabs( retCPU[idx] - retGPU[idx] ) ;
#ifdef VERBOSE_COMPARISON
	printf( "%d %zu %d (%f %f) == (%f %f)\n" , T , i , p ,
		creal(retCPU[idx]) , cimag(retCPU[idx]) ,
		creal(retGPU[idx]) , cimag(retGPU[idx]) ) ;
#endif
      }
    }
    std::cout<<"Summed diff p="<<p<<" "<<sum<<std::endl ;
  }
  free( retCPU ) ;
#endif
  free( retGPU ) ;
  
  finalize();

  return 0;
}
