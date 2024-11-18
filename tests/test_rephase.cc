#include "init_quda_laph.h"

int main( int argc , char *argv[])
{
  XMLHandler xml_in ;

  const int ntasks = init_quda_laph( argc , argv , xml_in ) ;

  // call rephase here

  
  finalize() ;
  
  return 0 ;
}
