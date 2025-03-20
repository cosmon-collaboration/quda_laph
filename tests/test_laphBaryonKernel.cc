#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

using namespace LaphEnv ;

void
old_code( std::vector<LattField> &laph_evecs )
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

  const int nmom = 1 , blockSizeMomProj = 1 , n1 = 1 , n2 = 1 , n3 = 1 ;
  const int X[4] = {
    LayoutInfo::getRankLattExtents()[0],
    LayoutInfo::getRankLattExtents()[1],
    LayoutInfo::getRankLattExtents()[2],
    LayoutInfo::getRankLattExtents()[3] } ;
  const int nsites = X[0]*X[1]*X[2] ;
    
  double _Complex coeffs1[ Nev*n1*2 ] = {} ;
  double _Complex coeffs2[ Nev*n2*2 ] = {} ;
  double _Complex coeffs3[ Nev*n3*2 ] = {} ;

  std::cout<<"Memsets?"<<std::endl ;

  memset( coeffs1 , 1.0 , Nev*n1*sizeof( double _Complex) ) ;
  memset( coeffs2 , 1.0 , Nev*n2*sizeof( double _Complex) ) ;
  memset( coeffs3 , 1.0 , Nev*n3*sizeof( double _Complex) ) ;
  
  // what is this? exponentiated mom? Which index runs fastest? Fuck knows
  double _Complex host_mom[ nsites*nmom*2 ] = {} ;
  for( int p = 0 ; p < nmom ; p++ ) {
    for( int z = 0 ; z < X[2] ; z++ ) {
      for( int y = 0 ; y < X[1] ; y++ ) {
	for( int x = 0 ; x < X[0] ; x++ ) {
	  // lexi site ignoring evenodd because this is just a nonzero bullshit entry
	  const int i = x + X[0]*( y + X[1]*z ) ;
	  host_mom[ i + nsites*p ] = cos( i*2*3.14159/X[0] ) ;
	}
      }
    }
  }

  std::cout<<"LaphBaryonKernel"<<std::endl ;
  
  double _Complex return_array[ n1*n2*n3*nmom*2 ] = {} ;
  void **ev_list = (void**)evList.data() ;
  
  laphBaryonKernel( n1, n2, n3,
		    nmom,
		    coeffs1 ,
		    coeffs2 ,
		    coeffs3 ,
		    host_mom ,
		    1 , //Nev ,
		    ev_list,
		    return_array ,
		    blockSizeMomProj,
		    X ) ;
  
  finalize();

  return 0;
}
