/**
   Test that color color contraction is sensible and then do a stress test
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

// ok so Dean's code wanted #2 to be daggered already instead of doing the inner product, weird and perhaps wrong ...
static void
cpuColorContract( void **host_evec , void *result , const int nEv , const int X[4] , const int A , const int B )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ;
  std::complex<double> *ptA = (std::complex<double>*)host_evec[A] ;
  std::complex<double> *ptB = (std::complex<double>*)host_evec[B] ;
  std::complex<double> *ptC = (std::complex<double>*)result ;  
#pragma omp parallel for
  for( size_t i = 0 ; i < (size_t)Nsites ; i++ ) {
    // simple inner product over color A_{c}B*_{c}
    ptC[i]  = ptA[0+3*i]*(ptB[0+3*i]) ;
    ptC[i] += ptA[1+3*i]*(ptB[1+3*i]) ;
    ptC[i] += ptA[2+3*i]*(ptB[2+3*i]) ;  
  }
}

static void
gpuColorContract( void **host_evec , void *result , const int nEv , const int X[4] , const int A , const int B )
{
  const int Nsites = X[0]*X[1]*X[2]*X[3] ;

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
  const size_t dbytes = Nsites*sizeof(std::complex<double>) ;
  void *d_tmp = pool_device_malloc(dbytes);
  
  // just cross 0 with 1 for now
  colorContractQuda( quda_evec[0] , quda_evec[1] , d_tmp ) ;

  qudaMemcpy(result,d_tmp,dbytes,qudaMemcpyDeviceToHost) ;

  // pool's closed
  pool_device_free(d_tmp) ;
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

int main(int argc, char *argv[]) {
  XMLHandler xml_in;

  if( init_quda_laph(argc, argv, xml_in) != 0 ) {
    exit(1) ;
  }

  int global = 1 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
#endif

  const int Nev = 16 ;
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);
  set_constant( laphEigvecs ) ;
  LattField resultCPU( FieldSiteType::Complex ) ;
  LattField resultGPU( FieldSiteType::Complex ) ;
  
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
  cpuColorContract( evList.data() , (void*)resultCPU.getDataPtr() , Nev , X , 0 , 1 ) ;

  gpuColorContract( evList.data() , (void*)resultGPU.getDataPtr() , Nev , X , 0 , 1 ) ;

  // compare them
  std::complex<double> sum = 0 ;
  std::complex<double>*cpu = (std::complex<double>*)resultCPU.getDataPtr() ;
  std::complex<double>*gpu = (std::complex<double>*)resultGPU.getDataPtr() ;
  for( size_t i = 0 ; i < (size_t)LayoutInfo::getRankLatticeNumSites() ; i++ ) {
    sum += fabs( cpu[i] - gpu[i] ) ;
  }
  std::cout<<"CPU - GPU :: "<<sum<<std::endl ;
  
  finalize();

  return 0;
}
