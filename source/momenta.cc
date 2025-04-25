#include "momenta.h"

using namespace std;

namespace LaphEnv {

vector<int> Momentum::getMomentumVec() const
{
  vector<int> res(3);
  res[0] = x;
  res[1] = y;
  res[2] = z;
  return res;
}

   //  compute the GCD of two positive integers
/*   
uint eval_gcd(uint a, uint b) 
{
 int aa=a;
 int bb=b;
 int temp;
 while (bb>0){
    temp=bb;
    bb=aa % bb;
    aa=temp;}
 return aa;
}


    //  compute the GCD of three positive integers
    
uint eval_gcd(uint a, uint b, uint c)
{
 return eval_gcd(eval_gcd(a,b),c);
}


   //   Returns true if (px,py,pz) is an allowed momentum ray
   //   and returns the momentum ray string in "ray".  To determine
   //   the ray, the momentum is rescaled so that (2,2,2) becomes (1,1,1).
   //   In other words, the momentum is divided by the GCD of the
   //   three components.  An allowed ray is one in which the
   //   sum of the squares is less than or equal to maxmomsq.
 
 
bool getMomentumRay(int px, int py, int pz, string& ray, uint maxmomsq)
{
 ray.clear();
 int gcd=eval_gcd((px<0)?-px:px,(py<0)?-py:py,(pz<0)?-pz:pz);
 int ppx=px, ppy=py, ppz=pz;
 if (gcd>1){
    ppx/=gcd; ppy/=gcd; ppz/=gcd;}
 uint momsq=ppx*ppx+ppy*ppy+ppz*ppz;
 if (momsq>maxmomsq) return false;
 if (momsq>8) return false;     // prevent ppx,ppy,ppz from being 3
 char dir[5]={'=','-','0','+','#'};
 ray="mom_ray_"; ray+=dir[ppx+2]; ray+=dir[ppy+2]; ray+=dir[ppz+2];
 return true;
}

*/
  // *****************************************************************
}
