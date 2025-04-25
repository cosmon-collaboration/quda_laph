#ifndef MOMENTA_H
#define MOMENTA_H

#include "multi_compare.h"
#include <vector>


namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Defines a class "Momentum" which is just a convenient struct  *
// *   for storing a three-momentum defined as three integers (since *
// *   momentum is quantized on a toroid).                           *
// *                                                                 *
// *******************************************************************


struct Momentum
{
 int x,y,z;

 Momentum(){}
 Momentum(int px, int py, int pz) : x(px), y(py), z(pz) {}
 Momentum(const Momentum& rhs) : x(rhs.x), y(rhs.y), z(rhs.z) {}
 Momentum& operator=(const Momentum& rhs)
  {x=rhs.x; y=rhs.y; z=rhs.z; return *this;}
 bool operator<(const Momentum& rhs) const
  {return multiLessThan(x,rhs.x,  y,rhs.y,  z,rhs.z);}
 bool operator==(const Momentum& rhs) const
  {return multiEqual(x,rhs.x,  y,rhs.y,  z,rhs.z);}
 bool operator!=(const Momentum& rhs) const
  {return multiNotEqual(x,rhs.x,  y,rhs.y,  z,rhs.z);}

 std::vector<int> getMomentumVec() const;
};


   //   Returns true if (px,py,pz) is an allowed momentum ray
   //   and returns the momentum ray string in "ray".  To determine
   //   the ray, the momentum is rescaled so that (2,2,2) becomes (1,1,1).
   //   In other words, the momentum is divided by the GCD of the
   //   three components.  An allowed ray is one in which the
   //   sum of the squares is less than or equal to maxmomsq.

//bool getMomentumRay(int px, int py, int pz, std::string& ray, uint maxmomsq=6);

// **********************************************************
}
#endif  
