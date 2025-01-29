#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"

#include "laph_noise.h"
#include <iomanip>
#include <cfloat>

std::vector<std::complex<double>>
old_code( const int zn )
{
  std::vector<std::complex<double>>values(zn) ;
  const double one = 1.0, zero = 0.0;
  if (zn == 4) {
    values[0] = std::complex<double>(one, zero);
    values[1] = std::complex<double>(zero, one);
    values[2] = std::complex<double>(-one, zero);
    values[3] = std::complex<double>(zero, -one);
  } else if (zn == 8) {
    const double sqrthalf = 0.70710678118654752440;
    values[0] = std::complex<double>(one, zero);
    values[1] = std::complex<double>(sqrthalf, sqrthalf);
    values[2] = std::complex<double>(zero, one);
    values[3] = std::complex<double>(-sqrthalf, sqrthalf);
    values[4] = std::complex<double>(-one, zero);
    values[5] = std::complex<double>(-sqrthalf, -sqrthalf);
    values[6] = std::complex<double>(zero, -one);
    values[7] = std::complex<double>(sqrthalf, -sqrthalf);
  } else if (zn == 32) {
    const double c0 = 0.98078528040323044912; // cos(Pi/16)
    const double c1 = 0.92387953251128675613; // cos(Pi/8)
    const double c2 = 0.83146961230254523707; // cos(3*Pi/16)
    const double c3 = 0.70710678118654752440; // sqrt(1/2);
    const double c4 = 0.55557023301960222475; // cos(5*Pi/16)
    const double c5 = 0.38268343236508977173; // cos(3*Pi/8)
    const double c6 = 0.19509032201612826785; // cos(7*Pi/16)
    values[0] = std::complex<double>(one, zero);
    values[1] = std::complex<double>(c0, c6);
    values[2] = std::complex<double>(c1, c5);
    values[3] = std::complex<double>(c2, c4);
    values[4] = std::complex<double>(c3, c3);
    values[5] = std::complex<double>(c4, c2);
    values[6] = std::complex<double>(c5, c1);
    values[7] = std::complex<double>(c6, c0);
    values[8] = std::complex<double>(zero, one);
    values[9] = std::complex<double>(-c6, c0);
    values[10] = std::complex<double>(-c5, c1);
    values[11] = std::complex<double>(-c4, c2);
    values[12] = std::complex<double>(-c3, c3);
    values[13] = std::complex<double>(-c2, c4);
    values[14] = std::complex<double>(-c1, c5);
    values[15] = std::complex<double>(-c0, c6);
    values[16] = std::complex<double>(-one, zero);
    values[17] = std::complex<double>(-c0, -c6);
    values[18] = std::complex<double>(-c1, -c5);
    values[19] = std::complex<double>(-c2, -c4);
    values[20] = std::complex<double>(-c3, -c3);
    values[21] = std::complex<double>(-c4, -c2);
    values[22] = std::complex<double>(-c5, -c1);
    values[23] = std::complex<double>(-c6, -c0);
    values[24] = std::complex<double>(zero, -one);
    values[25] = std::complex<double>(c6, -c0);
    values[26] = std::complex<double>(c5, -c1);
    values[27] = std::complex<double>(c4, -c2);
    values[28] = std::complex<double>(c3, -c3);
    values[29] = std::complex<double>(c2, -c4);
    values[30] = std::complex<double>(c1, -c5);
    values[31] = std::complex<double>(c0, -c6);
  } else if (zn == 1) {
    values[0] = std::complex<double>(one, zero);
  } else {
    std::cout<<zn<<" not supported"<<std::endl ;
    exit(1) ;
  }
  return values ;
}

void
do_test( const int zn )
{
  // Z4
  LaphEnv::LaphZnNoise noise4( zn , 1234 ) ;
  std::cout<<"Elements of Z" << zn <<std::endl ;
  std::vector<std::complex<double>> res = old_code( zn ) ;
  bool fail = false ;
  for( int i = 0 ; i < zn ; i++ ) {
    if( abs(noise4.values[i] - res[i] ) > 0.0 ) {
      std::cout << std::setprecision(19) << "Noise " << noise4.values[i] << " " << res[i] << " " << abs(noise4.values[i] - res[i])<<std::endl ;
      fail = true ;
    }
  }
  if( fail == true ) {
    std::cout<<"ZN test failed"<<std::endl ;
    exit(1) ;
  }
  return ;
}

int main(int argc, char *argv[]) {
  XMLHandler xml_in;

  if( init_quda_laph(argc, argv, xml_in) != 0 ) {
    exit(1) ;
  }

  do_test( 4 ) ;
  do_test( 8 ) ;
  do_test( 32 ) ;
  
  finalize();
  
  return 0;
}
