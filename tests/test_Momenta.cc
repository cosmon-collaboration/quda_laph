#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "momenta.h"

using namespace LaphEnv ;

int main(int argc, char *argv[]) {
  XMLHandler xml_in;

  if( init_quda_laph(argc, argv, xml_in) != 0 ) {
    exit(1) ;
  }

  int global = 1 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
#endif

  const std::vector<Momentum> Mom = { Momentum( 1 , 2 , 3 ) ,
				      Momentum( 2 , 1 , 3 ) ,
				      Momentum( 2 , 5 , 3 ) } ;
  // what Quda contractor wants
  Momenta<double> p( Mom ) ;
  const std::complex<double> *host_mom = p.phases.data() ;
  for( size_t i = 0 ; i < p.phases.size() ; i++ ) {
    std::cout<<host_mom[i]<<std::endl ;
  }
  
  finalize();

  return 0;
}
