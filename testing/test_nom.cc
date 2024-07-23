#include "task_tests.h"
#include "latt_field.h"
#include <complex>
#include "named_obj_map.h"
#include "laph_stdio.h"

using namespace std;
using namespace LaphEnv;

typedef std::complex<double> cmplx;
typedef std::complex<float> fcmplx;

namespace QLTestEnv {

// ************************************************

void testNamedObjectMap(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestNamedObjectMap")==0)
 return;

 XMLHandler dxml("RootValue","Double");
 double value1=4.5;
 NamedObjMap::insert("key1",value1,dxml);

// XMLHandler ixml("RootValue","Integer");
// NamedObject<int> iobj(7,ixml);
 
 XMLHandler vcdxml("RootValue","VectorOfCmplxDouble");
 std::vector<std::complex<double>> vcdvals;
 vcdvals.push_back(cmplx(3.265,-8.1232));
 vcdvals.push_back(cmplx(0.535,3.194));
 vcdvals.push_back(cmplx(-7.75,6.22));
 NamedObjMap::insert("key2",vcdvals,vcdxml);

 double res1=NamedObjMap::getData<double>("key1");
 printLaph(make_strf("res1 = %f",res1));
 
 NamedObjMap::getData<double>("key1")=8.3;
 printLaph(make_strf("res1 = %f",NamedObjMap::getData<double>("key1")));
 XMLHandler xml1; NamedObjMap::getXMLInfo("key1",xml1);
 printLaph(make_strf("dobj xml info is %s",xml1.str()));
 
 XMLHandler xmld("Dummy");
 xmld.put_child("Pet","Rabbit");
 xmld.put_child("Job","Butcher");
 NamedObjMap::setXMLInfo("key1",xmld); NamedObjMap::getXMLInfo("key1",xml1);
 printLaph(make_strf("dobj xml info is %s",xml1.str()));
 
 const vector<cmplx>& res2=NamedObjMap::getData<vector<cmplx>>("key2");
 for (uint k=0;k<res2.size();++k){
    printLaph(make_strf("res2[%d] = (%f,%f)",k,res2[k].real(),res2[k].imag()));}

 vector<cmplx>& res2b=NamedObjMap::getData<vector<cmplx>>("key2");
 res2b[1]=cmplx(12.5,56.0);
 res2b.push_back(cmplx(99.9,88.8));
 for (uint k=0;k<res2b.size();++k){
    printLaph(make_strf("res2b[%d] = (%f,%f)",k,res2b[k].real(),res2b[k].imag()));}
 
 NamedObjMap::getXMLInfo("key2",xml1);
 printLaph(make_strf("vcdobj xml info is %s",xml1.str()));

 
 vector<int>& ivals=NamedObjMap::insert<vector<int>>("key3");
 ivals.push_back(6);
 ivals.push_back(3);
 ivals.push_back(9);

 //cmplx& z=NamedObjMap::insert("key4",cmplx(4.5,-6.2));
 
 {cmplx z1(88.5,-44.2);
 NamedObjMap::insert("key5",z1);}

 printLaph(make_strf("key5 has value (%f,%f)",NamedObjMap::getData<cmplx>("key5").real(),
      NamedObjMap::getData<cmplx>("key5").imag()));

 double& dv=NamedObjMap::getData<double>("key1");
 printLaph(make_strf("dv has value %f",dv));
 dv=7;
 
 double& dx=NamedObjMap::getData<double>("key1");
 printLaph(make_strf("dx has value %f",dx));

 double *dptr=0;
 NamedObjMap::getData<double>("key1",dptr);
 
 vector<int> *ivptr=0;
 NamedObjMap::getData("key3",ivptr);
 (*ivptr)[1]=18;
 ivptr->push_back(33);
 for (uint k=0;k<ivptr->size();++k)
    printLaph(make_strf("iv[%d]=%d",k,(*ivptr)[k]));
 XMLHandler xmlget0;
 printLaph(make_strf("empty xml handler has output %s",xmlget0.output()));
 NamedObjMap::getXMLInfo("key3",xmlget0);
 printLaph(make_strf("key3 initially has xmlinfo %s",xmlget0.output()));
 XMLHandler xmlnew("IVector","Donuts");
 NamedObjMap::setXMLInfo("key3",xmlnew);
 XMLHandler xmlget;
 NamedObjMap::getXMLInfo("key3",xmlget);
 printLaph(make_strf("key3 now has xmlinfo %s",xmlget.output()));


 const vector<int> *icvptr=0;
 NamedObjMap::getData("key3",icvptr);
 for (uint k=0;k<icvptr->size();++k)
    printLaph(make_strf("icv[%d]=%d",k,(*icvptr)[k]));

 vector<int> *ivptr2=0;
 NamedObjMap::getDataAndXMLInfo("key3",ivptr2,xml1);
 (*ivptr2)[1]=24;
 ivptr2->push_back(11);
 for (uint k=0;k<ivptr2->size();++k)
    printLaph(make_strf("iv2[%d]=%d",k,(*ivptr2)[k]));
 printLaph(make_strf(" xmlinfo is %s",xml1.output()));

 const vector<int> *icvptr2=0;
 XMLHandler xml1c;
 NamedObjMap::getDataAndXMLInfo("key3",icvptr2,xml1c);
 for (uint k=0;k<icvptr2->size();++k)
    printLaph(make_strf("icv2[%d]=%d",k,(*icvptr2)[k]));
 printLaph(make_strf(" xmlinfo is %s",xml1c.output()));

 NamedObjMap::erase("key3");
 
 try{
    const vector<int> *icvptr3=0;
    NamedObjMap::getDataAndXMLInfo("key3",icvptr3,xml1c);}
 catch(const std::exception& x){
    printLaph("could not find key3, which is correct!");}
 
 printLaph(make_strf(" query for key1 is %u",NamedObjMap::query("key1")));
 printLaph(make_strf(" query for key2 is %u",NamedObjMap::query("key2")));
 printLaph(make_strf(" query for key3 is %u",NamedObjMap::query("key3")));
 printLaph(make_strf(" query for key4 is %u",NamedObjMap::query("key4")));
 printLaph(make_strf(" query for key5 is %u",NamedObjMap::query("key5")));
 printLaph(make_strf(" query for key6 is %u",NamedObjMap::query("key6")));

 set<string> idkeys=NamedObjMap::getIds();
 printLaph("The keys in the NamedObjmap as now as follows:");
 for (set<string>::const_iterator it=idkeys.begin();it!=idkeys.end();++it){
    printLaph(make_strf("%s",*it));}

 LatticeAssigner Q(0);
 {LattField tmp;
 tmp.reset(FieldSiteType::ColorVector);
 Q.assign_field(tmp,"First ColorVector");
 NamedObjMap::insert("LattColVector",tmp);
 Q.reSeed(673215432);
 Q.assign_field(tmp,"Second ColorVector");
 XMLHandler xmllat("Lattice","SU3ColorVector");
 NamedObjMap::insert("LattColVector2",tmp,xmllat);}

 const LattField *latptr=0;
 NamedObjMap::getData("LattColVector2",latptr);
 printLaph(make_strf("LattColVector2 has elements per site of %d",latptr->elemsPerSite()));
 LatticeChecker LC(Q);
 LC.check_field(*latptr,"Second ColorVector",true);

 NamedObjMap::getData("LattColVector",latptr);
 Q.reSeed(0);
 LC.check_field(*latptr,"First ColorVector",true);

}

// ***********************************************
}

