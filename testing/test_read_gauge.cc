#include "task_tests.h"
#include <fstream>
#include <vector>
#include <complex>
#include "latt_field.h"
#include "utils.h"
#include "read_gauge_field.h"
#include "layout_info.h"
#include "laph_stdio.h"
#include "quda_info.h"

using namespace std; 
using namespace LaphEnv;


namespace QLTestEnv {

void convert_sp_to_dp(const vector<LattField>& GFsp, vector<LattField>& GFdp)
{
 GFdp.resize(LayoutInfo::Ndim);
 int nreal=2*GFsp[0].elemsPerSite()*LayoutInfo::getRankLatticeNumSites();
 for (int dir=0;dir<LayoutInfo::Ndim;++dir){
    GFdp[dir].reset_by_precision(FieldSiteType::ColorMatrix,'D');
    const float* src=reinterpret_cast<const float*>(GFsp[dir].getDataConstPtr());
    double* dest=reinterpret_cast<double*>(GFdp[dir].getDataPtr());
    for (int k=0;k<nreal;++k,++src,++dest){
       *dest=*src;}}
}

// ************************************************

    // This is a slow function for testing purposes only; do not use
    // too large of a lattice.
 
void writeCERNconfig(const vector<LattField>& u, const string& fname)
{
 if (u.size()!=4){
    throw(std::invalid_argument("Can only writeCERNconfig if LattField is 4 dimensional"));}
 for (int dir=0;dir<4;++dir){
    if (u[dir].bytesPerSite()!=9*sizeof(complex<double>)){
       throw(std::invalid_argument("Can only writeCERNconfig if LattField is double precision"));}}
 printLaph(make_strf("Dummy Lattice to be written to file %s\n",fname));
 int NT = LayoutInfo::getLattSizes()[3];
 int NZ = LayoutInfo::getLattSizes()[2];
 int NY = LayoutInfo::getLattSizes()[1];
 int NX = LayoutInfo::getLattSizes()[0];
 int su3bytes=2*FieldNcolor*FieldNcolor*sizeof(double);
 ofstream fout(fname,ios::binary);
 if (isPrimaryRank()){
    fout.write((char*)(&NT),sizeof(int));
    fout.write((char*)(&NX),sizeof(int));
    fout.write((char*)(&NY),sizeof(int));
    fout.write((char*)(&NZ),sizeof(int));}
 double avgplaq=0.9852;
 if (isPrimaryRank()){
    fout.write((char*)(&avgplaq),sizeof(double));}
#ifdef ARCH_PARALLEL
 comm_barrier();
#endif
 for (int it=0;it<NT;it++)
   for (int ix=0;ix<NX;ix++)
     for (int iy=0;iy<NY;iy++)
       for (int iz=0;iz<NZ;iz++)
         if ((ix+iy+iz+it)%2){
            vector<int> site(LayoutInfo::Ndim);
            site[0]=ix; site[1]=iy; site[2]=iz; site[3]=it;
            int dn=ix-1; if (dn<0) dn=NX-1;
            vector<int> xdn(LayoutInfo::Ndim);
            xdn[0]=dn; xdn[1]=iy; xdn[2]=iz; xdn[3]=it;
            dn=iy-1; if (dn<0) dn=NY-1;
            vector<int> ydn(LayoutInfo::Ndim);
            ydn[0]=ix; ydn[1]=dn; ydn[2]=iz; ydn[3]=it;
            dn=iz-1; if (dn<0) dn=NZ-1;
            vector<int> zdn(LayoutInfo::Ndim);
            zdn[0]=ix; zdn[1]=iy; zdn[2]=dn; zdn[3]=it;
            dn=it-1; if (dn<0) dn=NT-1;
            vector<int> tdn(LayoutInfo::Ndim);
            tdn[0]=ix; tdn[1]=iy; tdn[2]=iz; tdn[3]=dn;
            vector<char> temp1(u[3].getSiteData(site));
            vector<char> temp2(u[3].getSiteData(tdn));
            vector<char> temp3(u[0].getSiteData(site));
            vector<char> temp4(u[0].getSiteData(xdn));
            vector<char> temp5(u[1].getSiteData(site));
            vector<char> temp6(u[1].getSiteData(ydn));
            vector<char> temp7(u[2].getSiteData(site));
            vector<char> temp8(u[2].getSiteData(zdn));
#ifdef ARCH_PARALLEL
            comm_barrier();
#endif
            if (isPrimaryRank()){
               fout.write((char*)(temp1.data()),su3bytes);
               fout.write((char*)(temp2.data()),su3bytes);
               fout.write((char*)(temp3.data()),su3bytes);
               fout.write((char*)(temp4.data()),su3bytes);
               fout.write((char*)(temp5.data()),su3bytes);
               fout.write((char*)(temp6.data()),su3bytes);
               fout.write((char*)(temp7.data()),su3bytes);
               fout.write((char*)(temp8.data()),su3bytes);}
#ifdef ARCH_PARALLEL
            comm_barrier();
#endif
            }
 fout.close();
 printLaph(make_strf("DONE: Dummy Lattice written to file %s\n",fname));
}


template <typename T>
void compare_configs(const vector<LattField>& u1, const vector<LattField>& u2,
                     const T& tolerance)
{  
 printLaph("Comparing lattice gauge fields");
 int nsites=LayoutInfo::getRankLatticeNumSites();
 int ndble=2*nsites*u1[0].elemsPerSite();
 T maxdiff=0.0;
 for (int dir=0;dir<LayoutInfo::Ndim;++dir){
    const T* p1=reinterpret_cast<const T*>(u1[dir].getDataConstPtr());
    const T* p2=reinterpret_cast<const T*>(u2[dir].getDataConstPtr());
    for (int k=0;k<ndble;++k,++p1,++p2){
       T df=std::abs((*p1)-(*p2));
       if (df>maxdiff){maxdiff=df; printf("MISMATCH dir = %d k = %d %g %g\n",dir,k,*p1,*p2);}
       //else {printf("dir = %d k = %d %g %g\n",dir,k,*p1,*p2);}
       
       }}
 T glmaxdiff=globalMax(maxdiff);
 if (glmaxdiff>tolerance){
    printLaph("The two configurations DISAGREE\n");
    printLaph(make_strf("Maximum difference is %g\n",glmaxdiff));}
 else{
    printLaph("The two configurations agree\n");}
}


void compare_configs(const vector<LattField>& u1, const vector<LattField>& u2)
{
 if (QudaInfo::get_cpu_prec()==QUDA_DOUBLE_PRECISION){
    double tolerance=1e-11;
    compare_configs(u1,u2,tolerance);}
 else{
    float tolerance=1e-6;
    compare_configs(u1,u2,tolerance);}
}


void testReadGaugeConfig(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestReadGaugeConfig")==0)
 return;

 printLaph("Starting testReadGaugeConfig\n");
 uint MTseed1=0, MTseed2=0, MTseed3=0, MTseed4=0;
 XMLHandler xmlr(xml_in,"TestReadGaugeConfig");
 xmlreadif(xmlr,"MTSeed1",MTseed1,"TestReadGaugeConfig");
 xmlreadif(xmlr,"MTSeed2",MTseed2,"TestReadGaugeConfig");
 xmlreadif(xmlr,"MTSeed3",MTseed3,"TestReadGaugeConfig");
 xmlreadif(xmlr,"MTSeed4",MTseed4,"TestReadGaugeConfig");
 string cfg_file_CERNformat;
 xmlreadif(xmlr,"CERNConfigFileName",cfg_file_CERNformat,"TestReadGaugeConfig");

 printLaph(make_strf("MTSeed1 is %d\n",MTseed1));
 printLaph(make_strf("MTSeed2 is %d\n",MTseed2));
 printLaph(make_strf("MTSeed3 is %d\n",MTseed3));
 printLaph(make_strf("MTSeed4 is %d\n",MTseed4));

 GaugeFieldAssigner GC(MTseed1,MTseed2,MTseed3,MTseed4);
 vector<LattField> GF1;
 GC.assign_gauge_field(GF1,"GaugeField1");

     // Test CERN read:  first, write out a config in CERN format

 if (GF1[0].bytesPerSite()==18*sizeof(double)){
    writeCERNconfig(GF1,cfg_file_CERNformat);}
 else{
    printf("Single precision mode: convert to double precision to write CERN config\n");
    vector<LattField> GF1dp;
    convert_sp_to_dp(GF1,GF1dp);
    writeCERNconfig(GF1dp,cfg_file_CERNformat);}
#ifdef ARCH_PARALLEL
 comm_barrier();
#endif

     //  Now read the config
 GaugeCERNConfigReader GR; 
 vector<LattField> Uread;
 GR.read(Uread,cfg_file_CERNformat);

     // Check the correctness of the read
 compare_configs(Uread,GF1);  
}


// ***********************************************
}
