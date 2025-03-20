#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

using namespace LaphEnv ;

void
cpu_code( std::vector<LattField> &laph_evecs )
{
}

static void
printevecs( std::vector<LattField> laphEigvecs , const int rank )
{
#ifdef PRINT_EVECS
  int myrank = rank ;
#ifdef ARCH_PARALLEL
  MPI_Comm_rank(  MPI_COMM_WORLD  , &myrank ) ;
#endif
  for( size_t n = 0 ; n < laphEigvecs.size() ; n++ ) {
    if( myrank == rank ) {
      std::cout<<"Ev "<<n<<" | rank "<<rank<<std::endl ;
      std::complex<double> *pt = (std::complex<double>*)laphEigvecs[n].getDataPtr() ;
      for( size_t i = 0 ; i < (size_t)LayoutInfo::getRankLatticeNumSites() ; i++ ) {
	for( size_t c = 0 ; c < (size_t)FieldNcolor ; c++ ) {
	  std::cout<<i<<"/"<<c<<" "<<pt[c+FieldNcolor*i]<<std::endl ;
	}
      }
      std::cout<<std::endl ;
      fflush( stdout ) ;
    }
  }
#endif
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
  
  // call rephase here
  const int Nev = 8 ;
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);

  std::cout<<"Constant Eigvecs"<<std::endl ;
  set_constant( laphEigvecs ) ;
  for( int i = 0 ; i < global ; i++ ) {
    printevecs( laphEigvecs , i ) ;
  }

  std::vector<void*> evList( Nev ) ;
  for( int i = 0 ; i < Nev ; i++ ) {
    evList[i] = (void*)laphEigvecs[i].getDataPtr() ;
  }

  const int nmom = 1 , blockSizeMomProj = 1 ;
  const int X[4] = { LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3]	} ;
  double _Complex host_mom[512] = {} ;
  double _Complex return_array[512] = {} ;
  
  laphBaryonKernelComputeModeTripletA( nmom , Nev , blockSizeMomProj ,
				       evList.data() ,
				       host_mom ,
				       return_array ,
				       X ) ;
  
  finalize();

  return 0;
}
