#include <iostream>
#include <vector>
using namespace std;

//   (x,y,z,t) column-major

int linearSiteIndex_lexico(const vector<int>& loc_coords, const vector<int>& loc_sizes)
{
 return loc_coords[0]+loc_sizes[0]*(loc_coords[1]+loc_sizes[1]*(loc_coords[2]
        +loc_sizes[2]*loc_coords[3]));
}

int linearSiteIndex_evenodd(const vector<int>& loc_coords, const vector<int>& loc_sizes,
                            int loc_nsites, int start_parity)
{
 return (loc_coords[0]+loc_sizes[0]*(loc_coords[1]+loc_sizes[1]*(loc_coords[2]
        +loc_sizes[2]*loc_coords[3])))/2 + loc_nsites/2*((start_parity+loc_coords[0]
        +loc_coords[1]+loc_coords[2]+loc_coords[3])%2);
}



int main(int argc, char** argv)
{
 vector<vector<int>> tests;
   // all even
 tests.push_back(vector<int>{12,14,22,24});
   // one odd
 tests.push_back(vector<int>{13,14,22,24});
 tests.push_back(vector<int>{12,15,22,24});
 tests.push_back(vector<int>{12,14,21,24});
 tests.push_back(vector<int>{12,14,22,25});
   // two odd
 tests.push_back(vector<int>{13,15,22,24});
 tests.push_back(vector<int>{13,14,21,24});
 tests.push_back(vector<int>{13,14,22,23});
 tests.push_back(vector<int>{12,15,21,24});
 tests.push_back(vector<int>{12,15,22,23});
 tests.push_back(vector<int>{12,14,21,23});
   // three odd
 tests.push_back(vector<int>{12,15,21,23});
 tests.push_back(vector<int>{13,14,21,23});
 tests.push_back(vector<int>{13,15,22,23});
 tests.push_back(vector<int>{13,15,21,24});

 for (int kk=0;kk<tests.size();++kk)
 for (int start_parity=0;start_parity<2;++start_parity){

 vector<int> locSizes=tests[kk];

 int loc_nsites=locSizes[0]*locSizes[1]*locSizes[2]*locSizes[3];
 vector<int> loc_coords(4);
 
 cout << "locSizes = ("<<locSizes[0]<<","<<locSizes[1]<<","<<locSizes[2]<<","<<locSizes[3]<<")"<<endl<<endl;
 cout << "start parity = "<<start_parity<<endl<<endl;
    // lexicographic check
 cout << endl<<"Lexicographic checks"<<endl<<endl;
 int count=0;
 bool success=true;
 for (loc_coords[3]=0;loc_coords[3]<locSizes[3];++loc_coords[3])
 for (loc_coords[2]=0;loc_coords[2]<locSizes[2];++loc_coords[2])
 for (loc_coords[1]=0;loc_coords[1]<locSizes[1];++loc_coords[1])
 for (loc_coords[0]=0;loc_coords[0]<locSizes[0];++loc_coords[0],++count){
    int calc=linearSiteIndex_lexico(loc_coords,locSizes);
    cout << "local site ("<<loc_coords[0]<<","<<loc_coords[1]<<","<<loc_coords[2]<<","<<loc_coords[3]
         <<")  count = "<<count<<"  calc = "<<calc;
    if (calc!=count){ cout << " MISMATCH"<<endl; success=false;}
    else {cout<<endl;}}
 if (success){ cout << "Was successful!"<<endl;}
 else{ cout << "Some comparisons failed"<<endl;}
     // even-odd checks
 cout << endl<<"Even-odd checkerboard checks"<<endl<<endl;
 count=0;
 success=true;
 for (int par=0;par<2;++par)
 for (loc_coords[3]=0;loc_coords[3]<locSizes[3];++loc_coords[3])
 for (loc_coords[2]=0;loc_coords[2]<locSizes[2];++loc_coords[2])
 for (loc_coords[1]=0;loc_coords[1]<locSizes[1];++loc_coords[1])
 for (loc_coords[0]=0;loc_coords[0]<locSizes[0];++loc_coords[0])
 if (((loc_coords[0]+loc_coords[1]+loc_coords[2]+loc_coords[3]+start_parity)%2)==par){
    int calc=linearSiteIndex_evenodd(loc_coords,locSizes,loc_nsites,start_parity);
    cout << "local site ("<<loc_coords[0]<<","<<loc_coords[1]<<","<<loc_coords[2]<<","<<loc_coords[3]
         <<")  count = "<<count<<"  calc = "<<calc;
    if (calc!=count){ cout << " MISMATCH"<<endl; success=false;}
    else {cout<<endl;}
    ++count;}
 if (success){ cout << "Was successful!"<<endl;}
 else{ cout << "Some comparisons failed"<<endl;}
 }

 return 0;
}
