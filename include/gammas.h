#include <array>
#include <iostream>
#include <cassert>
#include <cstdint>

typedef enum {
  Gamma_X , Gamma_Y , Gamma_Z , Gamma_T , Identity , Gamma_5 
} GAMMAS ;
  
typedef enum {
  p1 , pi , m1 , mi
} Z4 ;

// complex conjugate
static inline Z4 gconj( const Z4 &g ) {
  switch( g ) {
  case p1 : return p1 ;
  case pi : return mi ;
  case m1 : return m1 ;
  case mi : return pi ;
  }
  assert( 1 ) ; return p1 ;
}
// multiply by a set of z4
static inline Z4 z4mul( const Z4 &g , const Z4 &z4 ) {
  switch( (g + z4 )&3 ){
  case 0 : return p1 ;
  case 1 : return pi ;
  case 2 : return m1 ;
  case 3 : return mi ;
  }
  assert( 1 ) ; return p1 ;
}
  
// individual gamma matrix
class NRgamma
{
public :
  // value is +1,i,-1,-i mapped to 0,1,2,3
  std::array<Z4,4> g ;
  // index of each element
  std::array<uint8_t,4> ig ;
  std::string name ;
  
  // constructor
  NRgamma( const GAMMAS idx = Identity ) {
    assert( idx >= 0 && idx < 6 ) ;
    switch( idx ) {
    case Gamma_X : // gamma_x
      ig[0] = 3 ; g[0] = mi ;
      ig[1] = 2 ; g[1] = mi ;
      ig[2] = 1 ; g[2] = pi ;
      ig[3] = 0 ; g[3] = pi ;
      name = "gx" ;
      break ;
    case Gamma_Y : // gamma_y
      ig[0] = 3 ; g[0] = m1 ;
      ig[1] = 2 ; g[1] = p1 ;
      ig[2] = 1 ; g[2] = p1 ;
      ig[3] = 0 ; g[3] = m1 ;
      name = "gy" ;
      break ;
    case Gamma_Z : // gamma_z
      ig[0] = 2 ; g[0] = mi ;
      ig[1] = 3 ; g[1] = pi ;
      ig[2] = 0 ; g[2] = pi ;
      ig[3] = 1 ; g[3] = mi ;
      name = "gz" ;
      break ;
    case Gamma_T : // gamma_t
      ig[0] = 0 ; g[0] = p1 ;
      ig[1] = 1 ; g[1] = p1 ;
      ig[2] = 2 ; g[2] = m1 ;
      ig[3] = 3 ; g[3] = m1 ;
      name = "gt" ;
      break ;
    case Identity : // identity
      ig[0] = 0 ; g[0] = p1 ;
      ig[1] = 1 ; g[1] = p1 ;
      ig[2] = 2 ; g[2] = p1 ;
      ig[3] = 3 ; g[3] = p1 ;
      name = "I" ;
      break ;
    case Gamma_5 : // gamma_5
      ig[0] = 2 ; g[0] = m1 ;
      ig[1] = 3 ; g[1] = m1 ;
      ig[2] = 0 ; g[2] = m1 ;
      ig[3] = 1 ; g[3] = m1 ;
      name = "g5" ;
      break ;
    }
  }

  NRgamma( const std::array<Z4,4> &_g ,
	   const std::array<uint8_t,4> &_ig ,
	   const std::string &_name ) :
    g(_g), ig(_ig) , name(_name) { }

  friend NRgamma operator*( NRgamma b , const NRgamma &c ) {
    NRgamma a ;
    a.ig[0] = c.ig[ b.ig[0] ] ; a.g[0] = z4mul( b.g[0] , c.g[ b.ig[0] ] ) ;
    a.ig[1] = c.ig[ b.ig[1] ] ; a.g[1] = z4mul( b.g[1] , c.g[ b.ig[1] ] ) ;
    a.ig[2] = c.ig[ b.ig[2] ] ; a.g[2] = z4mul( b.g[2] , c.g[ b.ig[2] ] ) ;
    a.ig[3] = c.ig[ b.ig[3] ] ; a.g[3] = z4mul( b.g[3] , c.g[ b.ig[3] ] ) ;
    a.name = b.name + "." + c.name ;
    return a ;
  }
  
  // this is not a dagger, just a complex conjugate of the elements
  NRgamma conjugate()
  {
    NRgamma a ;
    a.g[0] = gconj( g[0] ) ; a.ig[0] = ig[0] ;
    a.g[1] = gconj( g[1] ) ; a.ig[1] = ig[1] ;
    a.g[2] = gconj( g[2] ) ; a.ig[2] = ig[2] ;
    a.g[3] = gconj( g[3] ) ; a.ig[3] = ig[3] ;
    a.name = "(" + name + ")^*" ;
    return a ;
  }

  NRgamma mul( const Z4 z4 )
  {
    NRgamma a(g,ig,name) ;
    for( int mu = 0 ; mu < 4 ; mu++ ) {
      a.g[mu] = z4mul( g[mu] , z4 ) ;
    }
    switch( z4 ) {
    case p1 : a.name = name ; break ;
    case pi : a.name = "+i" + name ; break ;
    case m1 : a.name = "-1" + name ; break ;
    case mi : a.name = "-i" + name ; break ;
    }
    return a ;
  }

  NRgamma dagger()
  {
    NRgamma a ;
    for( int mu = 0 ; mu < 4 ; mu++ ) {
      a.g[mu] = gconj( g[ ig[mu] ] ) ; a.ig[ ig[mu] ] = mu ;
    }
    a.name = "(" + name + ")^dag" ;
    return a ;
  }

  NRgamma transpose()
  {
    NRgamma a ;
    for( int mu = 0 ; mu < 4 ; mu++ ) {
      a.g[mu] = g[ ig[mu] ] ; a.ig[ ig[mu] ] = mu ;
    }
    a.name = "(" + name + ")^T" ;
    return a ;
  }

  NRgamma gt_Gdagger_gt()
  {
    NRgamma a(g,ig,name) , gt( Gamma_T ) ;
    return gt*a.dagger()*gt ;
  }

  // i.gy.gt.Gmu
  NRgamma cGmu()
  {
    NRgamma a(g,ig,name) , gy( Gamma_Y ) , gt( Gamma_T ) ;
    return gy.mul(pi)*gt*a ;
  }
  
  bool operator==(const NRgamma &b ) const 
  {
    for( int i = 0 ; i < 4 ; i++ ) {
      if( g[i] != b.g[i] || ig[i] != b.ig[i] ) return false ;
    }
    return true ;
  }
 
  void print( )
  {
    std::cout<<std::endl<<name<<std::endl ;
    for( int i = 0 ; i < 4 ; i++ ) {
      for( int j = 0 ; j < 4 ; j++ ) {
	if( ig[ i ] == j ) {
	  switch( g[i] ) {
	  case p1 : printf( " +1 " ) ; break ;
	  case pi : printf( " +i " ) ; break ;
	  case m1 : printf( " -1 " ) ; break ;
	  case mi : printf( " -i " ) ; break ;
	  }
	} else {
	  printf( " 0 " ) ;
	}
      }
      printf( "\n" ) ;
    }
  }
} ;

class Gbasis {
private :
  static void
  unit( const bool passed , const std::string &message )
  {
    std::cout<< message ;
    if( passed == true ) {
      std::cout<< " passed" << std::endl ;
      return ;
    }
    std::cout<< " failed" << std::endl ;
    assert( passed ) ;
  }
public :
  std::array<NRgamma,16> G ;
  Gbasis( const bool verbose = false ) {
    G[ Gamma_X  ] = NRgamma( Gamma_X ) ;
    G[ Gamma_Y  ] = NRgamma( Gamma_Y ) ;
    G[ Gamma_Z  ] = NRgamma( Gamma_Z ) ;
    G[ Gamma_T  ] = NRgamma( Gamma_T ) ;
    G[ Identity ] = NRgamma( Identity ) ;
    G[ Gamma_5  ] = NRgamma( Gamma_5 ) ;
    // all the rest are axials 
    G[6] = G[ Gamma_X ]*G[ Gamma_5 ] ; // gxg5
    G[7] = G[ Gamma_Y ]*G[ Gamma_5 ] ; // gyg5
    G[8] = G[ Gamma_Z ]*G[ Gamma_5 ] ; // gzg5
    G[9] = G[ Gamma_T ]*G[ Gamma_5 ] ; // gtg5
    // tensors
    G[10] = G[ Gamma_Y ]*G[ Gamma_Z ] ; // gygz \propto gx
    G[11] = G[ Gamma_Z ]*G[ Gamma_X ] ; // gzgx \propto gy 
    G[12] = G[ Gamma_X ]*G[ Gamma_Y ] ; // gygx \propto gz
    // sigma
    G[13] = G[ Gamma_X ]*G[ Gamma_T ] ;
    G[14] = G[ Gamma_Y ]*G[ Gamma_T ] ;
    G[15] = G[ Gamma_Z ]*G[ Gamma_T ] ;
    if( verbose == true ) {	
      for( int i = 0 ; i < 16 ; i++ ) {
	G[i].print() ;
      }
    }
  }

  // unit tests of this guy
  void
  unit_tests( void )
  {
    std::cout<<std::endl<<"Unit-testing gamma identities"<<std::endl<<std::endl ;
    unit( G[ Gamma_X ]*G[ Gamma_Y ]*G[ Gamma_Z ]*G[ Gamma_T ] == G[ Gamma_5 ] ,
	"g5 == gx.gy.gz.gt" ) ;

    for( int mu = 0 ; mu < 6 ; mu++ ) {
      NRgamma tmp = G[mu] * G[mu] ;
      unit( tmp == G[Identity] , tmp.name+" == " + G[Identity].name ) ; 
    }
    
    // test conjugate identity for euclidean matrices
    for( int mu = 0 ; mu < 6 ; mu++ ) {
      NRgamma tmp = G[mu].dagger() ;
      unit( tmp == G[mu] , tmp.name+" == "+G[mu].name ) ;
    }
    
    // test C.C^\dagger == I
    NRgamma tmp = G[ Identity ].cGmu( ) ;
    unit( tmp*tmp.dagger() == G[Identity] ,
	  (tmp*tmp.dagger()).name + " == " + G[Identity].name ) ;
    unit( tmp == tmp.transpose().mul( m1 ) ,
	  tmp.name + " == " + tmp.transpose().mul( m1 ).name ) ;
    
    for( int mu = 0 ; mu < 4 ; mu++ ) {
      NRgamma tmp1 = (G[mu].cGmu())*(G[Identity].cGmu()) ;
      NRgamma tmp2 = G[mu].transpose().mul(m1) ;
      unit( tmp1 == tmp2 , tmp1.name + " == " + tmp2.name ) ;
    }
    
    // test C^2=1
    {
      NRgamma tmp1 = G[Identity].cGmu() ;
      tmp1 = tmp1*tmp1 ;
      unit( tmp1 == G[Identity] , tmp1.name + " == " + G[Identity].name ) ;
    }
    
    // test gt.I^\dagger.gt = I
    {
      NRgamma tmp1 = G[Identity].gt_Gdagger_gt() ;
      unit( tmp1 == G[Identity] , tmp1.name + " == " + G[Identity].name ) ;
    }
  }
} ;
