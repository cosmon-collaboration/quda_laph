/**
   Test that color cross product is sensible and then do a stress test
 **/
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

using namespace LaphEnv ;
using namespace quda ;

static void
cpuColorCross( void **host_evec , void *result , const int nEv , const int X[4] , const int A , const int B )
{
  std::complex<double> *ptA = (std::complex<double>*)host_evec[A] ;
  std::complex<double> *ptB = (std::complex<double>*)host_evec[B] ;
  std::complex<double> *ptC = (std::complex<double>*)result ;  
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < (size_t)LayoutInfo::getRankLatticeNumSites() ; i++ ) {
    // first row determinant
    ptC[ 3*i + 0 ] =  ptA[ 3*i + 1 ]*ptB[ 3*i + 2 ] - ptA[ 3*i + 2 ]*ptB[ 3*i + 1 ] ;
    ptC[ 3*i + 1 ] = -ptA[ 3*i + 0 ]*ptB[ 3*i + 2 ] + ptA[ 3*i + 2 ]*ptB[ 3*i + 0 ] ;
    ptC[ 3*i + 2 ] =  ptA[ 3*i + 0 ]*ptB[ 3*i + 1 ] - ptA[ 3*i + 1 ]*ptB[ 3*i + 0 ] ; 
  }
}

static void
cpuColorCrossStress( void **host_evec , void *result , const int nEv , const int X[4] )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ; 
#pragma omp parallel for collapse(3)
  for( int j = 0 ; j < nEv ; j++ ) {
    for( int k = 0 ; k < nEv ; k++ ) {
      for( int i = 0 ; i < Nsites ; i++ ) {
	std::complex<double> *ptA = (std::complex<double>*)host_evec[j] ;
	std::complex<double> *ptB = (std::complex<double>*)host_evec[k] ;
	std::complex<double> *ptC = (std::complex<double>*)result ;  
	// first row determinant
	ptC[ 3*i + 0 ] =  ptA[ 3*i + 1 ]*ptB[ 3*i + 2 ] - ptA[ 3*i + 2 ]*ptB[ 3*i + 1 ] ;
	ptC[ 3*i + 1 ] = -ptA[ 3*i + 0 ]*ptB[ 3*i + 2 ] + ptA[ 3*i + 2 ]*ptB[ 3*i + 0 ] ;
	ptC[ 3*i + 2 ] =  ptA[ 3*i + 0 ]*ptB[ 3*i + 1 ] - ptA[ 3*i + 1 ]*ptB[ 3*i + 0 ] ; 
      }
    }
  }
}

static void
gpuColorCross( void **host_evec , void *result , const int nEv , const int X[4] , const int A , const int B )
{  
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

  const lat_dim_t x = { X[0] , X[1] , X[2] , X[3] } ;
  ColorSpinorParam cpu_evec_param( host_evec , inv_param , x , false , QUDA_CPU_FIELD_LOCATION ) ;
  cpu_evec_param.gammaBasis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  cpu_evec_param.nSpin = 1 ;

  // cpu copies of all evecs
  std::vector<ColorSpinorField> evec(nEv) ;
  for (int iEv=0; iEv<nEv; ++iEv) {
    cpu_evec_param.v = host_evec[iEv];
    evec[iEv] = ColorSpinorField(cpu_evec_param) ;
  }
  
  // GPU -side params
  ColorSpinorParam cuda_evec_param(cpu_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION);
  cuda_evec_param.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);

  std::vector<ColorSpinorField> quda_evec(2) ;
  quda_evec[0] = ColorSpinorField(cuda_evec_param);
  quda_evec[0] = evec[A] ;
  quda_evec[1] = ColorSpinorField(cuda_evec_param);
  quda_evec[1] = evec[B] ;
  
  // create output params
  ColorSpinorParam cpu_res_param( result , inv_param , x , false , QUDA_CPU_FIELD_LOCATION ) ;
  cpu_res_param.gammaBasis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  cpu_res_param.nSpin = 1 ;
  cpu_res_param.v = result ; // point to the result
  ColorSpinorField res(cpu_res_param) ;
  ColorSpinorParam cuda_res_param(cpu_res_param, inv_param, QUDA_CUDA_FIELD_LOCATION);
  ColorSpinorField Diq(cuda_res_param);
 
  // just cross 0 with 1 for now
  colorCrossQuda( quda_evec[0] , quda_evec[1] , Diq ) ;

  // download the evecs
  res = Diq ;
}

static void
gpuColorCrossStress( void **host_evec , void *result , const int nEv , const int X[4] )
{  
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

  const lat_dim_t x = { X[0] , X[1] , X[2] , X[3] } ;
  ColorSpinorParam cpu_evec_param( host_evec , inv_param , x , false , QUDA_CPU_FIELD_LOCATION ) ;
  cpu_evec_param.gammaBasis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  cpu_evec_param.nSpin = 1 ;

  // cpu copies of all evecs
  std::vector<ColorSpinorField> evec(nEv) ;
  for (int iEv=0; iEv<nEv; ++iEv) {
    cpu_evec_param.v = host_evec[iEv];
    evec[iEv] = ColorSpinorField(cpu_evec_param) ;
  }
  
  // GPU -side params
  ColorSpinorParam cuda_evec_param(cpu_evec_param, inv_param, QUDA_CUDA_FIELD_LOCATION);
  cuda_evec_param.setPrecision(inv_param.cuda_prec, inv_param.cuda_prec, true);

  std::vector<ColorSpinorField> quda_evec(nEv) ;
  for( int i = 0 ; i < nEv ; i++ ) {
    quda_evec[i] = ColorSpinorField(cuda_evec_param);
    quda_evec[i] = evec[i] ;
  }
  
  // create output params
  ColorSpinorParam cpu_res_param( result , inv_param , x , false , QUDA_CPU_FIELD_LOCATION ) ;
  cpu_res_param.gammaBasis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  cpu_res_param.nSpin = 1 ;
  cpu_res_param.v = result ; // point to the result
  ColorSpinorField res(cpu_res_param) ;
  ColorSpinorParam cuda_res_param(cpu_res_param, inv_param, QUDA_CUDA_FIELD_LOCATION);
  ColorSpinorField Diq(cuda_res_param);
 
  // cross color products over NeV
  for( int i = 0 ; i < nEv ; i++ ) {
    for( int j = 0 ; j < nEv ; j++ ) {
      colorCrossQuda( quda_evec[i] , quda_evec[j] , Diq ) ;
    }
  }
  
  // download the evecs
  res = Diq ;
}

static void
set_constant( std::vector<LattField> &laphEigvecs )
{
  int myrank = 0 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_rank( MPI_COMM_WORLD , &myrank ) ;
#endif
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
	const std::complex<double> z( n+laphEigvecs.size()*(c+FieldNcolor*(i+V*myrank))+1 ,
				      n+laphEigvecs.size()*(c+FieldNcolor*(i+V*myrank))+1 ) ;
	*ptr = z ;
	ptr++ ;
      }
    }
    #endif
  }
}

static void
is_zero( void *result , const int X[4] )
{
  const size_t nsites = X[0]*X[1]*X[2]*X[3] ;
  const std::complex<double>*pt = (const std::complex<double>*)result ;
  std::complex<double> Sum = 0. ;
  for( size_t i = 0 ; i < 3*nsites ; i++ ) {
    Sum += abs(pt[i]) ;
  }
  std::cout<<"Should be zero "<<Sum<<std::endl ;
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
  
  const int Nev = 16 ;
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);
  set_constant( laphEigvecs ) ;
  LattField resultCPU( FieldSiteType::ColorVector ) ;
  LattField resultGPU( FieldSiteType::ColorVector ) ;
  
  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;

  StopWatch timerGPU , timerCPU ;

  // test this is zero CPU side
  cpuColorCross( evList.data() , (void*)resultCPU.getDataPtr() , Nev , X , 0 , 0 ) ;
  is_zero( (void*)resultCPU.getDataPtr() , X ) ;

  // test this is zero GPU side
  gpuColorCross( evList.data() , (void*)resultGPU.getDataPtr() , Nev , X , 0 , 0 ) ;
  is_zero( (void*)resultGPU.getDataPtr() , X ) ;
  
  // GPU code
  timerGPU.start() ;
  gpuColorCross( evList.data() , (void*)resultGPU.getDataPtr() , Nev , X , 0 , 1 ) ;
  timerGPU.stop() ;
  printLaph(make_strf("\nGPU color cross in = %g seconds\n", timerGPU.getTimeInSeconds()));
  
  // CPU code
  timerCPU.start() ;
  cpuColorCross( evList.data() , (void*)resultCPU.getDataPtr() , Nev , X , 0 , 1 ) ;
  timerCPU.stop() ;
  printLaph(make_strf("\nCPU color cross in = %g seconds\n", timerCPU.getTimeInSeconds()));

  // compare them
  std::complex<double> sum = 0 ;
  std::complex<double>*cpu = (std::complex<double>*)resultCPU.getDataPtr() ;
  std::complex<double>*gpu = (std::complex<double>*)resultGPU.getDataPtr() ;
  for( size_t i = 0 ; i < (size_t)LayoutInfo::getRankLatticeNumSites() ; i++ ) {
    for( int c = 0 ; c < 3 ; c++ ) {
      sum += fabs( cpu[c+3*i] - gpu[c+3*i] ) ;
    }
  }
  std::cout<<"CPU - GPU :: "<<sum<<std::endl ;

  timerGPU.reset() ;
  timerGPU.start() ;
  gpuColorCrossStress( evList.data() , (void*)resultGPU.getDataPtr() , Nev , X ) ;  
  timerGPU.stop() ;
  printLaph(make_strf("\nGPU stress NeV*Nev in = %g seconds\n", timerGPU.getTimeInSeconds()));

  timerCPU.reset() ;
  timerCPU.start() ;
  cpuColorCrossStress( evList.data() , (void*)resultCPU.getDataPtr() , Nev , X ) ;  
  timerCPU.stop() ;
  printLaph(make_strf("\nCPU stress NeV*Nev in = %g seconds\n", timerCPU.getTimeInSeconds()));
  
  finalize();

  return 0;
}
