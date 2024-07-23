#include "task_tests.h"
#include "array.h"
#include <cstdio> 
#include <cstdlib>
#include <ctime>
#include <complex>
#include "laph_stdio.h"
#include "multi_looper.h"
#include "utils.h"

//#include "stop_watch.h"
//#include <unistd.h>

using namespace std;
using namespace LaphEnv;

typedef std::complex<double> cmplx;
typedef std::complex<float> fcmplx;

namespace QLTestEnv {

// ************************************************

unsigned int inc(unsigned int& k, unsigned int maximum)
{
 unsigned int v=k;
 k++;
 if (k==maximum) k=0;
 return v;
}

bool neq(const uint& i1, const uint& i2)
{return (i1!=i2);}

bool neq(const int& i1, const int& i2)
{return (i1!=i2);}

bool neq(const float& x1, const float& x2)
{const double eps=1e-6;
 return std::abs(x1-x2)>std::max(eps,0.5*eps*(std::abs(x1)+std::abs(x2)));}

bool neq(const double& x1, const double& x2)
{const double eps=1e-12;
 return std::abs(x1-x2)>std::max(eps,0.5*eps*(std::abs(x1)+std::abs(x2)));}

bool neq(const fcmplx& x1, const fcmplx& x2)
{const double eps=1e-6;
 return std::abs(x1-x2)>std::max(eps,0.5*eps*(std::abs(x1)+std::abs(x2)));}

bool neq(const cmplx& x1, const cmplx& x2)
{const double eps=1e-12;
 return std::abs(x1-x2)>std::max(eps,0.5*eps*(std::abs(x1)+std::abs(x2)));}
 
bool eq(const uint& i1, const uint& i2)
{return (i1==i2);}

bool eq(const int& i1, const int& i2)
{return (i1==i2);}

bool eq(const float& x1, const float& x2)
{const double eps=1e-6;
 return std::abs(x1-x2)<std::max(eps,0.5*eps*(std::abs(x1)+std::abs(x2)));}

bool eq(const double& x1, const double& x2)
{const double eps=1e-12;
 return std::abs(x1-x2)<std::max(eps,0.5*eps*(std::abs(x1)+std::abs(x2)));}

bool eq(const fcmplx& x1, const fcmplx& x2)
{const double eps=1e-6;
 return std::abs(x1-x2)<std::max(eps,0.5*eps*(std::abs(x1)+std::abs(x2)));}

bool eq(const cmplx& x1, const cmplx& x2)
{const double eps=1e-12;
 return std::abs(x1-x2)<std::max(eps,0.5*eps*(std::abs(x1)+std::abs(x2)));}

  //   WARNING:  make sure SAFETY_FLAG is defined in array.h
  
  // Performs tests of Array<T> for dimension "ndim".  Integer values in "idata" are used
  // for sizing the test arrays, and "datacmp" is used for assigning values.

template <typename T>
bool dotestArray(unsigned int ndim, const vector<T>& datacmp, const vector<unsigned int>& idata)
{
 unsigned int successcount=0;
 unsigned int failcount=0;
 unsigned int datacount=0;
 unsigned int icount=0;
 unsigned int dmax=datacmp.size();
 unsigned int imax=idata.size();

    // test number of dimensions and sizes 
 Array<T> *Aptr=0;
 vector<unsigned int> csizes;
 if (ndim<1){
    Aptr=new Array<T>;}
 else if (ndim==1){
    Aptr=new Array<T>(2);
    csizes.resize(ndim); csizes[0]=2;}
 else if (ndim==2){
    Aptr=new Array<T>(2,3);
    csizes.resize(ndim); csizes[0]=2; csizes[1]=3;}
 else if (ndim==3){
    Aptr=new Array<T>(2,3,4);
    csizes.resize(ndim); csizes[0]=2; csizes[1]=3; csizes[2]=4;}
 else if (ndim==4){
    Aptr=new Array<T>(2,3,4,2);
    csizes.resize(ndim); csizes[0]=2; csizes[1]=3; csizes[2]=4;
    csizes[3]=2;}
 else if (ndim==5){
    Aptr=new Array<T>(2,3,4,2,5);
    csizes.resize(ndim); csizes[0]=2; csizes[1]=3; csizes[2]=4;
    csizes[3]=2; csizes[4]=5;}
 else if (ndim==6){
    Aptr=new Array<T>(2,3,4,2,5,3);
    csizes.resize(ndim); csizes[0]=2; csizes[1]=3; csizes[2]=4;
    csizes[3]=2; csizes[4]=5; csizes[5]=3;}
 else if (ndim==7){
    Aptr=new Array<T>(2,3,4,2,5,3,2);
    csizes.resize(ndim); csizes[0]=2; csizes[1]=3; csizes[2]=4;
    csizes[3]=2; csizes[4]=5; csizes[5]=3; csizes[6]=2;}
 else if (ndim==8){
    Aptr=new Array<T>(2,3,4,2,5,3,2,6);
    csizes.resize(ndim); csizes[0]=2; csizes[1]=3; csizes[2]=4;
    csizes[3]=2; csizes[4]=5; csizes[5]=3; csizes[6]=2; csizes[7]=6;}
 if (ndim<=8){
    unsigned int dimcheck=Aptr->numDimensions();
    if (dimcheck==ndim) successcount++;
    else{ failcount++; printf("number of dimensions test failed\n");}
    unsigned int sizecheck=Aptr->size();
    unsigned int ss=1; for (uint k=0;k<ndim;++k) ss*=csizes[k];
    if (ndim==0) ss=0;
    if (ss==sizecheck) successcount++;
    else{ failcount++; printf("size test failed\n");}
    const std::vector<uint>& scheck=Aptr->sizes();
    if (scheck==csizes) successcount++;
    else{ failcount++; printf("sizes test failed\n");}
    for (uint i=0;i<10;i++){
       try{
          unsigned int isize=Aptr->size(i);
          if ((i<ndim) && (isize==csizes[i])) successcount++;
          else{ failcount++; printf("size(i) test failed\n");}}
       catch(...){
          if (i>=ndim) successcount++;
          else { failcount++; printf("size(i) test failed\n");}}}
    delete Aptr; Aptr=0;}

 if (ndim==0){
    printLaph(make_strf("\n Number of successful tests = %d",successcount));
    printLaph(make_strf("     Number of FAILED tests = %d",failcount));
    return (failcount==0);}

   // use "idata" to determine sizes of arrays to test now
 csizes.resize(ndim);
 for (uint i=0;i<ndim;++i)
    csizes[i]=idata[inc(icount,imax)];
 Array<T> A0(csizes);

 for (int jj=0;jj<3;jj++){
 
  if (jj!=0){
     for (uint i=0;i<ndim;++i)
        csizes[i]=idata[inc(icount,imax)];
     A0.resize(csizes);}

    // more size tests
 unsigned int dimcheck=A0.numDimensions();
 if (dimcheck==ndim) successcount++;
 else{ failcount++; printf("number of dimensions test failed\n");}
 unsigned int sizecheck=A0.size();
 unsigned int ss=1; for (uint k=0;k<ndim;++k) ss*=csizes[k];
 if (ndim==0) ss=0;
 if (ss==sizecheck) successcount++;
 else{ failcount++; printf("size test failed\n");}
 const std::vector<uint>& scheck=A0.sizes();
 if (scheck==csizes) successcount++;
 else{ failcount++; printf("sizes test failed\n");}
 for (uint i=0;i<10;i++){
    try{
       unsigned int isize=A0.size(i);
       if ((i<ndim) && (isize==csizes[i])) successcount++;
       else{ failcount++; printf("size(i) test failed\n");}}
    catch(...){
       if (i>=ndim) successcount++;
       else { failcount++; printf("size(i) test failed\n");}}}

    // tests of content of arrays (value assignments)
 datacount=0;
 unsigned int dmax=datacmp.size();
 if (ndim==1){
    for (uint k0=0;k0<csizes[0];++k0){
       A0(k0)=datacmp[inc(datacount,dmax)]; 
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==2){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1){
       A0(k0,k1)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==3){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2){
       A0(k0,k1,k2)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==4){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3){
       A0(k0,k1,k2,k3)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==5){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3)
    for (uint k4=0;k4<csizes[4];++k4){
       A0(k0,k1,k2,k3,k4)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==6){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3)
    for (uint k4=0;k4<csizes[4];++k4)
    for (uint k5=0;k5<csizes[5];++k5){
       A0(k0,k1,k2,k3,k4,k5)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==7){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3)
    for (uint k4=0;k4<csizes[4];++k4)
    for (uint k5=0;k5<csizes[5];++k5)
    for (uint k6=0;k6<csizes[6];++k6){
       A0(k0,k1,k2,k3,k4,k5,k6)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==8){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3)
    for (uint k4=0;k4<csizes[4];++k4)
    for (uint k5=0;k5<csizes[5];++k5)
    for (uint k6=0;k6<csizes[6];++k6)
    for (uint k7=0;k7<csizes[7];++k7){
       A0(k0,k1,k2,k3,k4,k5,k6,k7)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}

 datacount=0;
 if (ndim==1){
    for (uint k0=0;k0<csizes[0];++k0){
       if (neq(A0(k0),datacmp[inc(datacount,dmax)])) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==2){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1){
       if (neq(A0(k0,k1),datacmp[inc(datacount,dmax)])) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==3){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2){
       if (neq(A0(k0,k1,k2),datacmp[inc(datacount,dmax)])) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==4){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3){
       if (neq(A0(k0,k1,k2,k3),datacmp[inc(datacount,dmax)])) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==5){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3)
    for (uint k4=0;k4<csizes[4];++k4){
       if (neq(A0(k0,k1,k2,k3,k4),datacmp[inc(datacount,dmax)])) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==6){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3)
    for (uint k4=0;k4<csizes[4];++k4)
    for (uint k5=0;k5<csizes[5];++k5){
       if (neq(A0(k0,k1,k2,k3,k4,k5),datacmp[inc(datacount,dmax)])) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==7){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3)
    for (uint k4=0;k4<csizes[4];++k4)
    for (uint k5=0;k5<csizes[5];++k5)
    for (uint k6=0;k6<csizes[6];++k6){
       if (neq(A0(k0,k1,k2,k3,k4,k5,k6),datacmp[inc(datacount,dmax)])) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==8){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3)
    for (uint k4=0;k4<csizes[4];++k4)
    for (uint k5=0;k5<csizes[5];++k5)
    for (uint k6=0;k6<csizes[6];++k6)
    for (uint k7=0;k7<csizes[7];++k7){
       if (neq(A0(k0,k1,k2,k3,k4,k5,k6,k7),datacmp[inc(datacount,dmax)])) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 if (failcount>0) printf("Failure in value assignment\n");
 }

 vector<unsigned int> isizes(ndim);
 for (uint k=0;k<ndim;++k) isizes[k]=idata[inc(icount,imax)];
 T cf=datacmp[inc(datacount,dmax)];
 Array<T> B0(isizes,cf);
 MultiIntLooper ind(isizes);
 for (ind.start();ind.notdone();++ind){
    if (eq(B0(ind.get_current()),cf)) successcount++; else failcount++;}
 if (failcount>0) printf("Value failed at point  2\n");

 cf=datacmp[inc(datacount,dmax)];
 cf*=datacmp[inc(datacount,dmax)];
 cf-=datacmp[inc(datacount,dmax)];
 B0=cf;
 for (ind.start();ind.notdone();++ind){
    if (eq(B0(ind.get_current()),cf)) successcount++; else failcount++;}
 if (failcount>0) printf("Value failed at point  3\n");

 T cf2=datacmp[inc(datacount,dmax)];
 cf+=cf2;
 B0+=cf2;
 for (ind.start();ind.notdone();++ind){
    if (eq(B0(ind.get_current()),cf)) successcount++; else failcount++;}
 if (failcount>0) printf("Value failed at point  4\n");

 cf2=datacmp[inc(datacount,dmax)];
 cf*=cf2;
 B0*=cf2;
 for (ind.start();ind.notdone();++ind){
    if (eq(B0(ind.get_current()),cf)) successcount++; else failcount++;}
 if (failcount>0) printf("Value failed at point  5\n");

 cf2=datacmp[inc(datacount,dmax)];
 cf-=cf2;
 B0-=cf2;
 for (ind.start();ind.notdone();++ind){
    if (eq(B0(ind.get_current()),cf)) successcount++; else failcount++;}
 if (failcount>0) printf("Value failed at point  6\n");

 cf2=datacmp[inc(datacount,dmax)];
 cf/=cf2;
 B0/=cf2;
 for (ind.start();ind.notdone();++ind){
    if (eq(B0(ind.get_current()),cf)) successcount++; else failcount++;}
 if (failcount>0) printf("Value failed at point  7\n");

 Array<T> W(B0);
 for (ind.start();ind.notdone();++ind){
    if (eq(W(ind.get_current()),B0(ind.get_current()))) successcount++; else failcount++;}
 if (failcount>0) printf("Value failed at point  8\n");

 for (uint k=0;k<ndim;++k) isizes[k]=idata[inc(icount,imax)];
 B0.resize(isizes);
 MultiIntLooper dil(isizes);
 datacount=0;
 for (dil.start();dil.notdone();++dil)
    B0(dil.get_current())=datacmp[inc(datacount,dmax)];
 datacount=dmax/2;
 Array<T> B1(isizes);
 for (dil.start();dil.notdone();++dil)
    B1(dil.get_current())=datacmp[inc(datacount,dmax)];

     //  test +=  and *=
 Array<T> B2(B0);
 B2+=B1;
 datacount=0;
 unsigned int dcount2=dmax/2;
 for (dil.start();dil.notdone();++dil){
    if (eq(B2(dil.get_current()),(datacmp[inc(datacount,dmax)]+datacmp[inc(dcount2,dmax)]))) successcount++; else failcount++;}
 if (failcount>0) printf("Failure in += \n");

 Array<T> B3(B0);
 B3*=B1;
 datacount=0;
 dcount2=dmax/2;
 for (dil.start();dil.notdone();++dil){
    if (eq(B3(dil.get_current()),(datacmp[inc(datacount,dmax)]*datacmp[inc(dcount2,dmax)]))) successcount++; else failcount++;}
 if (failcount>0) printf("Failure in *= \n");

#if ARCH_PARALLEL
 comm_barrier();
 successcount=globalMin(successcount);
 failcount=globalMax(failcount);
#endif
 printLaph(make_strf("\n Number of successful tests = %d",successcount));
 printLaph(make_strf("     Number of FAILED tests = %d",failcount));
 return (failcount==0);   
}   
    


bool dotestArrayConj(unsigned int ndim, const vector<complex<double> >& datacmp, 
                     const vector<unsigned int>& idata)
{
 unsigned int successcount=0;
 unsigned int failcount=0;
 unsigned int datacount=0;
 unsigned int icount=0;
 unsigned int dmax=datacmp.size();
 unsigned int imax=idata.size();

 vector<unsigned int> csizes(ndim);
 for (uint i=0;i<ndim;++i)
    csizes[i]=idata[inc(icount,imax)];
 Array<complex<double>> A0(csizes);

 if (ndim==1){
    for (uint k0=0;k0<csizes[0];++k0){
       A0(k0)=datacmp[inc(datacount,dmax)]; 
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==2){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1){
       A0(k0,k1)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==3){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2){
       A0(k0,k1,k2)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==4){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3){
       A0(k0,k1,k2,k3)=datacmp[inc(datacount,dmax)];
       if (datacount>=dmax) datacount=0;}}

 A0.conj();  // take the complex conjugate

 datacount=0;
 if (ndim==1){
    for (uint k0=0;k0<csizes[0];++k0){
       if (neq(A0(k0),conj(datacmp[inc(datacount,dmax)]))) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==2){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1){
       if (neq(A0(k0,k1),conj(datacmp[inc(datacount,dmax)]))) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==3){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2){
       if (neq(A0(k0,k1,k2),conj(datacmp[inc(datacount,dmax)]))) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}
 else if (ndim==4){
    for (uint k0=0;k0<csizes[0];++k0)
    for (uint k1=0;k1<csizes[1];++k1)
    for (uint k2=0;k2<csizes[2];++k2)
    for (uint k3=0;k3<csizes[3];++k3){
       if (neq(A0(k0,k1,k2,k3),conj(datacmp[inc(datacount,dmax)]))) failcount++; else successcount++;
       if (datacount>=dmax) datacount=0;}}

#if ARCH_PARALLEL
 comm_barrier();
 successcount=globalMin(successcount);
 failcount=globalMax(failcount);
#endif
 printLaph(make_strf("\n Number of successful tests = %d",successcount));
 printLaph(make_strf("     Number of FAILED tests = %d",failcount)); 
 return (failcount==0);
}  
    



void testArray(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestArray")==0)
 return;

/*
StopWatch bulova; bulova.start();
sleep(5);
bulova.stop();
cout << "bulova time 1 is "<<bulova.getTimeInSeconds()<<endl;
cout << "bulova last time interval 1 is "<<bulova.getLastIntervalInSeconds()<<endl;
sleep(4);
bulova.start();
sleep(3);
bulova.start();
bulova.stop();
bulova.stop();
cout << "bulova time 2 is "<<bulova.getTimeInSeconds()<<endl;
cout << "bulova last time interval 2 is "<<bulova.getLastIntervalInSeconds()<<endl;
bulova.reset();
sleep(4);
bulova.start();
sleep(7);
bulova.reset();
bulova.stop();
cout << "bulova time 3 is "<<bulova.getTimeInSeconds()<<endl;
cout << "bulova last time interval 3 is "<<bulova.getLastIntervalInSeconds()<<endl;
*/

 //unsigned int successcount=0;
 printLaph("Starting TestArray\n");
 bool success=true;
 //srand (time(NULL));
 unsigned int nintegers=2056;
 vector<unsigned int> idata(nintegers);
 for (unsigned int k=0;k<nintegers;++k)
    idata[k]=rand() % 8 + 1;

 unsigned int ndata=5048; 
 vector<int> intdata(ndata);
 for (unsigned int k=0;k<ndata;++k)
    intdata[k]= (rand() % 4096 + 1)-2047;

 for (int ndim=0;ndim<=9;++ndim){
    printLaph(make_strf("Tests: integer type, ndim = %d",ndim));
    success=success && dotestArray(ndim,intdata,idata);
    printLaph("\n");}

 ndata=5048;
 vector<unsigned int> uintdata(ndata);
 for (unsigned int k=0;k<ndata;++k)
    uintdata[k]= (rand() % 4096 + 1);

 for (int ndim=0;ndim<=9;++ndim){
    printLaph(make_strf("\n\nTests: unsigned integer type, ndim = %d",ndim));
    success=success && dotestArray(ndim,uintdata,idata);
    printLaph("\n");}

 ndata=5048;
 vector<float> fdata(ndata);
 for (unsigned int k=0;k<ndata;++k)
    fdata[k]= 1.434*((rand() % 4096 + 1)-2047)- 0.54511*((rand() % 128)-63);

 for (int ndim=0;ndim<=9;++ndim){
    printLaph(make_strf("\nTests: float type, %d",ndim));
    success=success && dotestArray(ndim,fdata,idata);
    printLaph("\n");}

 ndata=5048;
 vector<double> ddata(ndata);
 for (unsigned int k=0;k<ndata;++k)
    ddata[k]= -2.3418*((rand() % 4096 + 1)-2047)- 0.22485*((rand() % 128)-63);

 for (int ndim=0;ndim<=9;++ndim){
    printLaph(make_strf("\nTests: double type, %d",ndim));
    success=success && dotestArray(ndim,ddata,idata);
    printLaph("\n");}

 ndata=5048;
 vector<complex<float> > zfdata(ndata);
 for (unsigned int k=0;k<ndata;++k)
    zfdata[k]=complex<float>(1.434*((rand() % 4096 + 1)-2047),- 0.54511*((rand() % 128)-63));

 for (int ndim=0;ndim<=9;++ndim){
    printLaph(make_strf("\nTests: complex<float> type, %d",ndim));
    success=success && dotestArray(ndim,zfdata,idata);
    printLaph("\n");}

 ndata=5048;
 vector<complex<double> > zddata(ndata);
 for (unsigned int k=0;k<ndata;++k)
    zddata[k]=complex<double>(-2.3418*((rand() % 4096 + 1)-2047), 0.22485*((rand() % 128)-63));

 for (int ndim=0;ndim<=9;++ndim){
    printLaph(make_strf("\nTests: complex<double> type, %d",ndim));
    success=success && dotestArray(ndim,zddata,idata);
    printLaph("\n");}
 
 ndata=5048;
 vector<complex<double> > zdcdata(ndata);
 for (unsigned int k=0;k<ndata;++k)
    zdcdata[k]=complex<double>(-2.3418*((rand() % 4096 + 1)-2047), 0.22485*((rand() % 128)-63));

 for (int ndim=1;ndim<=4;++ndim){
    printLaph(make_strf("\nTests: conj complex<double> type, %d",ndim));
    success=success && dotestArrayConj(ndim,zdcdata,idata);
    printLaph("\n");}
 
// cout << "Testing Array outer product\n");
/* uint d1=4, d2=6;
 Array<double> W1(d1), W2(d2);
 for (uint k=0;k<d1;k++) W1(k)=10.0*((rand() % 4096 + 1)-2047);
 for (uint k=0;k<d2;k++) W2(k)=-3.2*((rand() % 4096 + 1)-2047);
 Array<double> W12; outerProduct(W1,W2,W12);
 for (uint k1=0;k1<d1;k1++)
 for (uint k2=0;k2<d2;k2++){
    if (W12(k1,k2)!=(W1(k1)*W2(k2))) success=false;} 

 d1=3; d2=5; uint d3=4, d4=5, d5=2; 
 Array<double> WW1(d1,d2,d3), WW2(d4,d5);
 for (uint k1=0;k1<d1;k1++) 
 for (uint k2=0;k2<d2;k2++) 
 for (uint k3=0;k3<d3;k3++) 
    WW1(k1,k2,k3)=1.04*((rand() % 4096 + 1)-2047);
 for (uint k4=0;k4<d4;k4++) 
 for (uint k5=0;k5<d5;k5++) 
    WW2(k4,k5)=1.04*((rand() % 4096 + 1)-2047);
 Array<double> WW12; outerProduct(WW1,WW2,WW12);
 for (uint k1=0;k1<d1;k1++)
 for (uint k2=0;k2<d2;k2++)
 for (uint k3=0;k3<d3;k3++)
 for (uint k4=0;k4<d4;k4++)
 for (uint k5=0;k5<d5;k5++){
    if (WW12(k1,k2,k3,k4,k5)!=(WW1(k1,k2,k3)*WW2(k4,k5))) success=false;}

 cout << "testing insertLeft 1\n");
 d1=3; d2=5; d3=4, d4=5, d5=2;
 Array<complex<double> > Z(d1,d2,d3,d4),Zcheck(d1,d2,d3,d4);
 Array<complex<double> > Zsubs(d1,d2,d3);
 for (uint k4=0;k4<d4;k4++){
 for (uint k3=0;k3<d3;k3++)
 for (uint k2=0;k2<d2;k2++)
 for (uint k1=0;k1<d1;k1++){
    double zr=(((k4*100.0+k3)*100.0+k2)*100.0+k1)/1321712.0;
    double zi=(((k1*98.0+k2)*98.0+k3)*98.0+k4)/2516231.0;
    Zsubs(k1,k2,k3)=complex<double>(zr,zi);
    Zcheck(k1,k2,k3,k4)=complex<double>(zr,zi);}
    Z.insertLeft(Zsubs,k4);}
 for (uint k4=0;k4<d4;k4++)
 for (uint k3=0;k3<d3;k3++)
 for (uint k2=0;k2<d2;k2++)
 for (uint k1=0;k1<d1;k1++){
    if (std::abs(Zcheck(k1,k2,k3,k4)-Z(k1,k2,k3,k4))>1e-13){
       cout << "insertLeft failed\n"); success=false;}}

 cout << "testing insertLeft 2\n");
 d1=3; d2=5; d3=4, d4=5, d5=2;
 Array<complex<double> > ZZ(d1,d2,d3,d4,d5),ZZcheck(d1,d2,d3,d4,d5);
 Array<complex<double> > ZZsubs(d1,d2,d3);
 for (uint k5=0;k5<d5;k5++)
 for (uint k4=0;k4<d4;k4++){
 for (uint k3=0;k3<d3;k3++)
 for (uint k2=0;k2<d2;k2++)
 for (uint k1=0;k1<d1;k1++){
    double zr=(((k4*100.0+k3)*100.0+k2)*100.0+k1)/1321712.0;
    double zi=(((k1*98.0+k2)*98.0+k3)*98.0+k4)/2516231.0;
    ZZsubs(k1,k2,k3)=complex<double>(zr,zi);
    ZZcheck(k1,k2,k3,k4,k5)=complex<double>(zr,zi);}
    ZZ.insertLeft(ZZsubs,k4,k5);}
 for (uint k5=0;k5<d5;k5++)
 for (uint k4=0;k4<d4;k4++)
 for (uint k3=0;k3<d3;k3++)
 for (uint k2=0;k2<d2;k2++)
 for (uint k1=0;k1<d1;k1++){
    if (std::abs(ZZcheck(k1,k2,k3,k4,k5)-ZZ(k1,k2,k3,k4,k5))>1e-13){
       cout << "insertLeft failed\n"); success=false;}}
*/
 success=globalAnd(success);
 if (success) printLaph("ALL TESTS PASSED!!");
 else printLaph("Some tests FAILED");
 
 printLaph("Testing map_insert_new\n");
 map<int,double> amap;
 int key=0;
 double val=0.0;
 double *ptr;
 success=true;
 try{
    ptr=map_insert_new(amap,key,val);
    if (*ptr!=val){ printLaph("TEST FAILED"); success=false;}}
 catch(...){printLaph("insertion failed"); success=false;}
 try{
    key=1; val=1.0;
    ptr=map_insert_new(amap,key,val);
    if (*ptr!=val){ printLaph("TEST FAILED"); success=false;}}
 catch(...){printLaph("insertion failed"); success=false;}
 try{
    key=5; val=5.0;
    ptr=map_insert_new(amap,key,val);
    if (*ptr!=val){ printLaph("TEST FAILED"); success=false;}}
 catch(...){printLaph("insertion failed"); success=false;}
 try{
    key=3; val=3.0;
    ptr=map_insert_new(amap,key,val);
    if (*ptr!=val){ printLaph("TEST FAILED"); success=false;}}
 catch(...){printLaph("insertion failed"); success=false;}
 try{
    key=3; val=3.0;
    ptr=map_insert_new(amap,key,val);
    printLaph("INSERTION SHOULD HAVE FAILED"); success=false;}
 catch(...){}  // This insertion SHOULD fail
 success=globalAnd(success);
 if (success) printLaph("All tests passed for map_insert_new");
 else printLaph("Some tests FAILED for map_insert_new");
}


void rand_array(uint ndim, uint ndil, Array<cmplx>& arr)
{
 vector<uint> msizes(ndim,ndil);
 arr.resize(msizes);
 MultiIntLooper v(msizes);
 for (v.start();v.notdone();++v){
    arr(v.get_current())=complex<double>(-2.3418*((rand() % 4096 + 1)-2047), 
                0.22485*((rand() % 128)-63));}
}

/*
void testTensorDataBuilder(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestTensorDataBuilder")==0)
 return;

 bool success=true;
 int ndil=24;

 cout << "Testing TensorDataBuilder"<<endl;

 vector<cmplx> coefs(3);
 coefs[0]=cmplx(1.0,-2.0);
 coefs[1]=cmplx(0.2,0.7);
 coefs[2]=cmplx(-2.2,0.44);

 vector<uint> arr_dims(2); arr_dims[0]=2; arr_dims[1]=3;
 uint nmeas=4;


 Array<Array<cmplx> > A0(3,nmeas),A1(3,nmeas);
 for (uint k=0;k<3;++k)
 for (uint j=0;j<nmeas;++j){
    rand_array(2,ndil,A0(k,j));
    rand_array(3,ndil,A1(k,j));}

 TensorDataBuilder TDB(coefs,arr_dims,nmeas);

 cout << "number of TDB terms = "<<TDB.getNumberOfTerms()<<endl;
 cout << "number of TDB meas = "<<TDB.getNumberOfMeasurements()<<endl;
 cout << "number of TDB hadrons = "<<TDB.getNumberOfHadrons()<<endl;
 const vector<cmplx>& cf=TDB.getCoefs();
 for (uint k=0;k<cf.size();++k) cout << "TDB coef["<<k<<"] = "<<cf[k]<<endl;

 for (uint k=0;k<3;++k)
 for (uint j=0;j<nmeas;++j)
    TDB.insert(A0(k,j),A1(k,j),k,j);

 for (uint k=0;k<3;++k)
 for (uint j=0;j<nmeas;++j){
    const Array<cmplx>& tref=TDB.getArray(0,k,j);
    if (tref(3,7)!=A0(k,j)(3,7)) success=false;
    const Array<cmplx>& t2ref=TDB.getArray(1,k,j);
    if (t2ref(3,7,5)!=A1(k,j)(3,7,5)) success=false;}


 cout << endl<<endl;
 if (success) cout << "ALL TENSORDATABUILDER TESTS PASSED!!"<<endl<<endl;
 else cout << "Some TensorDataBuilder tests FAILED"<<endl<<endl;
} */

// ***********************************************
}
