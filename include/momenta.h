#ifndef MOMENTA_H
#define MOMENTA_H

#include <set>
#include <vector>
#include <complex>
#include "multi_compare.h"
#include "layout_info.h"

namespace LaphEnv {
  struct Momentum {
    int x,y,z;
    Momentum(){}
    Momentum(const int px, const int py, const int pz) : x(px), y(py), z(pz) {}
    Momentum(const Momentum& rhs) : x(rhs.x), y(rhs.y), z(rhs.z) {}
    Momentum& operator=(const Momentum& rhs){x=rhs.x; y=rhs.y; z=rhs.z; return *this;}
    bool operator<(const Momentum& rhs) const {return multiLessThan(x,rhs.x,  y,rhs.y,  z,rhs.z);}
    bool operator==(const Momentum& rhs) const {return multiEqual(x,rhs.x,  y,rhs.y,  z,rhs.z);}
    bool operator!=(const Momentum& rhs) const {return multiNotEqual(x,rhs.x,  y,rhs.y,  z,rhs.z);}
    std::vector<int> getMomentumVec() const {
      return std::vector<int>{ x , y , z } ;
    }
  };

  // want to get a lattice-wide set of phases for contraction for GPU offload
  template <typename T>
  class Momenta
  {
  public :
    // great big flat array of exponentiated phases
    std::vector<std::complex<T>> phases ;
    // constructor
    Momenta( const std::vector<struct Momentum> &mom ) {
      const int locNX = LayoutInfo::getRankLattExtents()[0];
      const int locNY = LayoutInfo::getRankLattExtents()[1];
      const int locNZ = LayoutInfo::getRankLattExtents()[2];
      const int locNT = LayoutInfo::getRankLattExtents()[3];
      const size_t Lvol = (size_t)(locNX * locNY)*(size_t)(locNZ * locNT) ;
      // global values
      const int Lx = LayoutInfo::getLattExtents()[0] ;
      const int Ly = LayoutInfo::getLattExtents()[1] ;
      const int Lz = LayoutInfo::getLattExtents()[2] ;
      const int gx = LayoutInfo::getMyCommCoords()[0] ;
      const int gy = LayoutInfo::getMyCommCoords()[1] ;
      const int gz = LayoutInfo::getMyCommCoords()[2] ;
      const double TWOPI = 6.2831853071795864770 ;
      phases.resize( mom.size()*Lvol ) ;
      // parallel loop over phases, could flatten and do Lz*phases if we wished
      #pragma omp parallel for
      for( size_t i = 0 ; i < mom.size() ; i++ ) {
	// get momenta and precompute twiddles
	const std::vector<int> p = mom[i].getMomentumVec() ;
	const std::vector<double> twids = { p[0]*TWOPI/Lx , p[1]*TWOPI/Ly , p[2]*TWOPI/Lz } ;
	std::vector<std::complex<T>> lexi( Lvol ) ;
	// loop local xyzt
	for( int z = 0 ; z < locNZ ; z++ ) {
	  for( int y = 0 ; y < locNY ; y++ ) {
	    for( int x = 0 ; x < locNX ; x++ ) {
	      const double pdotx = (x+gx)*twids[0] + (y+gy)*twids[1] + (z+gz)*twids[2] ;	    
	      double sn , cn ;
	      sincos( pdotx , &sn , &cn ) ;
	      const std::complex<T> phase( (T)cn , -(T)sn ) ;
	      for( int t = 0 ; t < locNT ; t++ ) {
		lexi[ x + locNX*( y + locNY*( z + locNZ*t ) ) ] = phase ;
	      }
	    }
	  }
	}
	// set to evenodd
	LayoutInfo::lexico_to_evenodd( (char*)phases.data() + i*Lvol*sizeof( std::complex<T> ) ,
				       (char*)lexi.data() , sizeof( std::complex<T> ) ) ;
      }
    } ;
    ~Momenta() {
      phases.resize(0) ;
    }
  } ;
}
#endif
