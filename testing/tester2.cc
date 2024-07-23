#include <iostream>
#include <complex>
#include <functional>


using namespace std;

  struct GaugeFieldParam  {
    int nColor = 3;
    int nFace = 0;
    double anisotropy = 1.0;
    double tadpole = 1.0;
    double i_mu = 0.0;
 };

int main(int argc, char** argv)
{

 std::_Mem_fn<double GaugeFieldParam::*> g=std::mem_fn(&GaugeFieldParam::anisotropy);
 
 GaugeFieldParam gfp;
 gfp.tadpole=5.0;
 
 g(&gfp)=2.0;
 
 cout << gfp.nColor<<endl;
 cout << gfp.nFace<<endl;
 cout << gfp.anisotropy<<endl;
 cout << gfp.tadpole<<endl;
 cout << gfp.i_mu<<endl;
 
 return 0;
}
