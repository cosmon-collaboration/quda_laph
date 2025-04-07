#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

#include <cassert>
#include <complex.h>

#include <random>

using namespace LaphEnv ;

// will give zero for these otherwise default to uniform random numbers
//#define PSEUDOCONSTANT

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
  std::complex<double> *ptA = (std::complex<double>*)A ;
  std::complex<double> *ptB = (std::complex<double>*)B ;
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
  assert( X[3] == 1 ) ;
  
  // spatial only
  const size_t nSites = X[0]*X[1]*X[2] ;
  
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
	  double _Complex sum = 0.0 ;
	  const double _Complex *p1 = (const double _Complex*)tmp.getDataPtr() ;
	  const double _Complex *p2 = (const double _Complex*)host_mom+nSites*p ;
          #pragma omp parallel for reduction(+:sum)
	  for( size_t i = 0 ; i < nSites ; i++ ) {
	    sum += p2[i]*p1[i] ;
	  }
	  return_arr[ idx + nEvC3*p ] = sum ;
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
  assert( X[3] == 1 ) ;
  
  // spatial only
  const size_t nSites = X[0]*X[1]*X[2] ;
  
  const size_t nEvC3 = nEv*(nEv-1)*(nEv-2)/6 ;

  LattField Diq( FieldSiteType::ColorVector);

  size_t idx = 0 ;
  for( int aEv = 0 ; aEv < nEv ; aEv++ ) {
    for( int bEv = aEv+1 ; bEv < nEv ; bEv++ ) {
      cpuColorCross( host_evec[aEv] , host_evec[bEv] , (void*)Diq.getDataPtr() , X ) ;
      for( int cEv = bEv+1 ; cEv < nEv ; cEv++ ) {

	// inline contraction and momproj
	double _Complex sum[nMom] ;
	memset( sum , 0 , nMom*sizeof(double _Complex) ) ;
	const double _Complex *ptA = (const double _Complex*)Diq.getDataPtr() ;
	const double _Complex *ptB = (const double _Complex*)host_evec[cEv] ;

	#pragma omp parallel
	{
	  const double _Complex *p2 = (const double _Complex*)host_mom  ;
          #pragma omp for reduction(+:sum[:nMom])
	  for( size_t i = 0 ; i < (size_t)nSites ; i++ ) {
	    // is the color contraction
	    double _Complex Con = ptA[0+3*i]*(ptB[0+3*i]) ;
 	    Con += ptA[1+3*i]*(ptB[1+3*i]) ;
	    Con += ptA[2+3*i]*(ptB[2+3*i]) ;
            #pragma unroll(nMom)
	    for( int p = 0 ; p < nMom ; p++ ) {
	      sum[p] += p2[i+p*nSites]*Con ;
	    }
	  }
	}
	// copy back
	for( int p = 0 ; p < nMom ; p++ ) {
	  return_arr[ idx + nEvC3*p ] = sum[p] ;
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
  assert( X[3] == 1 ) ;
  // spatial only
  const size_t nSites = X[0]*X[1]*X[2] ;  
  const size_t nEvC3 = nEv*(nEv-1)*(nEv-2)/6 ;
  memset( return_arr , 0.0 , nEvC3*nMom*sizeof( double _Complex ) ) ;
  // loop over sites is far more efficient as we are only opening a parallel region once
#pragma omp parallel for reduction(+:return_arr[:nEvC3*nMom])
  for( size_t i = 0 ; i < nSites ; i++ ) {
    // precache these 
    const double _Complex *pt[ nEv ] ;
    for( int ev = 0 ; ev < nEv ; ev++ ) {
      pt[ ev ] = (const double _Complex*)host_evec[ev] + 3*i ;
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
	  #pragma unroll
	  for( int p = 0 ; p < nMom ; p++ ) {
	    return_arr[idx+nEvC3*p] += p2[p*nSites]*Con ;
	  }
	  idx++ ;
	}
      }
    }
  }
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
  assert( global == 1 ) ;
  
  // call rephase here
  const int Nev = 16 ;
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);

  std::cout<<"Constant Eigvecs"<<std::endl ;
  set_constant( laphEigvecs ) ;

  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 4 ;
  const int blockSizeMomProj = nmom ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int nspat  = X[0]*X[1]*X[2] ;
  
  double _Complex host_mom[ nmom*nspat ] = {} ;

  // host_mom should be complex
  for( size_t p = 0 ; p < nmom ; p++ ) {
    const int mom[3] = { p+1, p+2, p+3 } ;
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

  const size_t nEvChoose3 = Nev*(Nev-1)*(Nev-2)/6;
  double _Complex retGPU[ nmom*nEvChoose3 ] = {} ;

  laphBaryonKernelComputeModeTripletA( nmom, Nev,
				       blockSizeMomProj,
				       evList.data() ,
				       host_mom ,
				       retGPU, X ) ;
  memset( retGPU , 0.0 , nmom*nEvChoose3*sizeof( double _Complex )) ;

  StopWatch gpu ;
  gpu.start() ;
  laphBaryonKernelComputeModeTripletA( nmom,
				       Nev,
				       blockSizeMomProj,
				       evList.data() ,
				       host_mom ,
				       retGPU,
				       X ) ;
  gpu.stop() ;
  printLaph(make_strf("\nGPU modetripletA in = %g seconds\n", gpu.getTimeInSeconds()));

  {
    double _Complex retCPU[ nmom*nEvChoose3 ] = {} ;
    StopWatch cpu ;
    cpu.start() ;
    cpu_code( nmom, Nev, blockSizeMomProj, evList.data() , host_mom, retCPU, X ) ;
    cpu.stop() ;
    printLaph(make_strf("\nCPUv1 modetripletA in = %g seconds\n", cpu.getTimeInSeconds()));
  }
  {
    double _Complex retCPU[ nmom*nEvChoose3 ] = {} ;
    StopWatch cpu ;
    cpu.start() ;
    cpu_codev2( nmom, Nev, blockSizeMomProj, evList.data() , host_mom, retCPU, X ) ;
    cpu.stop() ;
    printLaph(make_strf("\nCPUv2 modetripletA in = %g seconds\n", cpu.getTimeInSeconds()));
  }

  double _Complex retCPU[ nmom*nEvChoose3 ] = {} ;
  StopWatch cpu ;
  cpu.start() ;
  cpu_codev3( nmom, Nev, blockSizeMomProj, evList.data() , host_mom, retCPU, X ) ;
  cpu.stop() ;
  printLaph(make_strf("\nCPUv3 modetripletA in = %g seconds\n", cpu.getTimeInSeconds()));
  
  // test outputs
  for( size_t p = 0 ; p < nmom ; p++ ) {
    double sum = 0 ;
    for( size_t i = 0 ; i < nEvChoose3 ; i++ ) {
      const size_t idx = i + nEvChoose3*p ;
      sum += cabs( retCPU[idx] - retGPU[idx] ) ;
      #if 0
      printf( " (%f %f) == (%f %f)\n" ,
	      creal(retCPU[idx1]) , cimag(retCPU[idx1]) ,
	      creal(retGPU[idx2]) , cimag(retGPU[idx2]) ) ;
      #endif
    }
    std::cout<<"Summed diff p="<<p<<" "<<sum<<std::endl ;
  }
  
  finalize();

  return 0;
}
