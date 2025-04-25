#include "field_ops.h"
#include "laph_stdio.h"
#include "utils.h"
#include "array.h"
#include <cstring>

#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif

#if defined(USE_GSL_CBLAS)
#include "gsl_cblas.h"
#elif defined(USE_OPENBLAS)
#include "cblas.h"
#endif

using namespace std;

namespace LaphEnv {


// **************************************************************************
// *                                                                        *
// *   This file contains routines for performing various lattice field     *
// *   operations on the CPU.  Typically, these routines are used for       *
// *   performing checks on solutions obtained by QUDA.  These routines     *
// *   are not particularly slow, but they are not particularly fast.       *
// *   They should only by used for occasional checks.                      *
// *                                                                        *
// **************************************************************************

    //  This routine applies an inner product 
    //         conj(leftfield).rightfield     if "left_conj" is true,
    //              leftfield .rightfield     if "left_conj" is false,
    //  where leftfield and rightfield are fields of the same site type, but the
    //  inner product is taken over the time slices.  The ntime inner
    //  products are returned.  

std::vector<complex<double>> getTimeSlicedInnerProducts(const LattField& leftfield, 
                                                        const LattField& rightfield, 
                                                        bool left_conj)
{
 if (leftfield.getFieldSiteType()!=rightfield.getFieldSiteType()){
    errorLaph("getTimeSlicedInnerProducts only supports fields of the same site type");}
 if (leftfield.bytesPerWord()!=rightfield.bytesPerWord()){
    errorLaph("getTimeSlicedInnerProducts requires both fields have same precision");}
 bool dp=(leftfield.bytesPerWord()==sizeof(std::complex<double>));
 int loc_nsites=LayoutInfo::getRankLatticeNumSites();
 int loc_npsites=loc_nsites/2;
 int start_parity=LayoutInfo::getMyStartParity();
 int nloctime=LayoutInfo::getRankLattExtents()[3];
 int tstride=LayoutInfo::getRankLattExtents()[0]*LayoutInfo::getRankLattExtents()[1]
            *LayoutInfo::getRankLattExtents()[2];
 int cbytes=(dp)?sizeof(std::complex<double>):sizeof(std::complex<float>);
 int bps=leftfield.bytesPerSite();
 int nelems_per_site=leftfield.elemsPerSite();

 vector<char> iprods1(nloctime*cbytes);
 vector<char> iprods2(nloctime*cbytes);
 const char* lp=reinterpret_cast<const char*>(leftfield.getDataConstPtr());
 const char* rp=reinterpret_cast<const char*>(rightfield.getDataConstPtr());
 char* ip1=reinterpret_cast<char*>(iprods1.data());
 char* ip2=reinterpret_cast<char*>(iprods2.data());
 for (int tloc=0;tloc<nloctime;++tloc,ip1+=cbytes,ip2+=cbytes){
    int parshift=loc_npsites*((start_parity+tloc)%2);
    int start1=((tstride*tloc)/2) + parshift;
    int stop1=((1+tstride*(tloc+1))/2) + parshift;
    int n1=nelems_per_site*(stop1-start1);
    parshift=loc_npsites*((start_parity+1+tloc)%2);
    int start2=((1+tstride*tloc)/2) + parshift;
    int stop2=((tstride*(tloc+1))/2) + parshift;
    int n2=nelems_per_site*(stop2-start2);
    const char* x1=lp+bps*start1;
    const char* x2=lp+bps*start2;
    const char* y1=rp+bps*start1;
    const char* y2=rp+bps*start2;
    if (dp){
       if (left_conj){
          cblas_zdotc_sub(n1,x1,1,y1,1,ip1);
          cblas_zdotc_sub(n2,x2,1,y2,1,ip2);}
       else{
          cblas_zdotu_sub(n1,x1,1,y1,1,ip1);
          cblas_zdotu_sub(n2,x2,1,y2,1,ip2);}}
    else{
       if (left_conj){
          cblas_cdotc_sub(n1,x1,1,y1,1,ip1);
          cblas_cdotc_sub(n2,x2,1,y2,1,ip2);}
       else{
          cblas_cdotu_sub(n1,x1,1,y1,1,ip1);
          cblas_cdotu_sub(n2,x2,1,y2,1,ip2);}}}

 int ntime=LayoutInfo::getLattExtents()[3];
 vector<complex<double>> iprods(ntime);
 for (int t=0;t<ntime;++t){
    iprods[t]=complex<double>(0.0,0.0);}
 int mytmin=LayoutInfo::getMyCommCoords()[3]*LayoutInfo::getRankLattExtents()[3];
 if (dp){
    complex<double>* z1=reinterpret_cast<complex<double>*>(iprods1.data());
    complex<double>* z2=reinterpret_cast<complex<double>*>(iprods2.data());
    for (int tloc=0;tloc<nloctime;++tloc,++z1,++z2){
       iprods[tloc+mytmin]+=(*z1)+(*z2);}}
 else{
    float* f1=reinterpret_cast<float*>(iprods1.data());
    float* f2=reinterpret_cast<float*>(iprods2.data());
    for (int tloc=0;tloc<nloctime;++tloc,++f1,++f2){
       double zr=(*f1)+(*f2); ++f1; ++f2;
       double zi=(*f1)+(*f2);
       iprods[tloc+mytmin]+=complex<double>(zr,zi);}}

#ifdef ARCH_PARALLEL
 vector<complex<double>> results(ntime);
 int status=MPI_Allreduce(iprods.data(),results.data(),2*ntime,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
 if (status!=MPI_SUCCESS){
    errorLaph("Problem occurred in getTimeSlicedInnerProducts");}
 return results;
#else
 return iprods;
#endif
}



void setConstantField(LattField& field, const std::complex<double>& zconst)
{
 bool dp=(field.bytesPerWord()==sizeof(complex<double>));
 int n=field.elemsPerSite()*LayoutInfo::getRankLatticeNumSites();
 if (dp){
    complex<double>* z=reinterpret_cast<complex<double>*>(field.getDataPtr());
    for (int k=0;k<n;++k,++z){
       *z=zconst;}}
 else{
    complex<float> zfconst(float(real(zconst)),float(imag(zconst)));
    complex<float>* z=reinterpret_cast<complex<float>*>(field.getDataPtr());
    for (int k=0;k<n;++k,++z){
       *z=zfconst;}}
}


void setUnitField(LattField& field)
{
 setConstantField(field,complex<double>(1.0,0.0));
}


void setZeroField(LattField& field)
{
 setConstantField(field,complex<double>(0.0,0.0));
}

    //  Returns the field   exp( -I * pvec.x )  in  "phases"

void makeMomentumPhaseField(LattField& phases, const Momentum& pvec)
{
 phases.reset(FieldSiteType::Complex);
 int locNt=LayoutInfo::getRankLattExtents()[3];
 int locNz=LayoutInfo::getRankLattExtents()[2];
 int locNy=LayoutInfo::getRankLattExtents()[1];
 int locNx=LayoutInfo::getRankLattExtents()[0];
 vector<complex<double>> xphases(locNx);
 vector<complex<double>> yphases(locNy);
 vector<complex<double>> zphases(locNz);
 double twopi=6.2831853071795864770;
 double momquantum=twopi*pvec.x / double(LayoutInfo::getLattExtents()[0]);
 int gx=LayoutInfo::getMyCommCoords()[0]*locNx;
 for (int x=0;x<locNx;++x,++gx){
    double xpx=momquantum*double(gx);
    xphases[x]=complex<double>(std::cos(xpx), -std::sin(xpx));}
 momquantum=twopi*pvec.y / double(LayoutInfo::getLattExtents()[1]);
 int gy=LayoutInfo::getMyCommCoords()[1]*locNy;
 for (int y=0;y<locNy;++y,++gy){
    double ypy=momquantum*double(gy);
    yphases[y]=complex<double>(std::cos(ypy), -std::sin(ypy));}
 momquantum=twopi*pvec.z / double(LayoutInfo::getLattExtents()[2]);
 int gz=LayoutInfo::getMyCommCoords()[2]*locNz;
 for (int z=0;z<locNz;++z,++gz){
    double zpz=momquantum*double(gz);
    zphases[z]=complex<double>(std::cos(zpz), -std::sin(zpz));}
   //  make phases in lexicographic order
   //  make phases first for time = 0
 LattField lexico(FieldSiteType::Complex);
 if (lexico.bytesPerWord()==sizeof(complex<double>)){
    complex<double>* rp=reinterpret_cast<complex<double>*>(lexico.getDataPtr());
    for (int z=0;z<locNz;++z){
       for (int y=0;y<locNy;++y){
          complex<double> zyphase=zphases[z]*yphases[y];
          for (int x=0;x<locNx;++x,++rp){
             *rp=xphases[x]*zyphase;}}}}
 else{
    complex<float>* rp=reinterpret_cast<complex<float>*>(lexico.getDataPtr());
    for (int z=0;z<locNz;++z){
       for (int y=0;y<locNy;++y){
          complex<double> zyphase=zphases[z]*yphases[y];
          for (int x=0;x<locNx;++x,++rp){
             *rp=xphases[x]*zyphase;}}}}
   //  copy phases from time = 0 to all other times
 const char* src=lexico.getDataConstPtr();
 char* dest=lexico.getDataPtr();
 size_t tslicebytes=locNx*locNy*locNz*lexico.bytesPerSite();
 for (int t=1;t<locNt;++t){
    dest+=tslicebytes;
    std::memcpy(dest,src,tslicebytes);}
   //  now rearrange into even-odd checkerboard
 LayoutInfo::lexico_to_evenodd(phases.getDataPtr(),lexico.getDataConstPtr(),
                               lexico.bytesPerWord());
}

           //  Evaluates sum_a conj( qbar[a](x,t) ) q[a](x,t)   a = color index
           //  (This could easily be threaded, but the fast version will run on the device)

void doColorContract(LattField& result, const LattField& qbar, const LattField& q)
{
 if ((qbar.getFieldSiteType()!=FieldSiteType::ColorVector)
     ||(q.getFieldSiteType()!=FieldSiteType::ColorVector)){
    errorLaph("doColorContract only works with ColorVector lattice fields");}
 result.reset(FieldSiteType::Complex);

 int s1=qbar.get_cpu_prec_bytes();
 int s2=q.get_cpu_prec_bytes();
 if (s1!=s2){
    errorLaph("doColorContract requires lattice fields have same precision");}
 int nsites=LayoutInfo::getRankLatticeNumSites();
 if (s1==sizeof(complex<double>)){
    const complex<double>* z1=reinterpret_cast<const complex<double>*>(qbar.getDataConstPtr());
    const complex<double>* z2=reinterpret_cast<const complex<double>*>(q.getDataConstPtr());
    complex<double>* rp=reinterpret_cast<complex<double>*>(result.getDataPtr());
    for (int k=0;k<nsites;++k){
       complex<double> res=std::conj(*z1)*(*z2); ++z1; ++z2;
       res+=std::conj(*z1)*(*z2); ++z1; ++z2;
       res+=std::conj(*z1)*(*z2); ++z1; ++z2;
       *rp=res; ++rp;}}
 else{
    const complex<float>* z1=reinterpret_cast<const complex<float>*>(qbar.getDataConstPtr());
    const complex<float>* z2=reinterpret_cast<const complex<float>*>(q.getDataConstPtr());
    complex<float>* rp=reinterpret_cast<complex<float>*>(result.getDataPtr());
    for (int k=0;k<nsites;++k){
       complex<float> res=std::conj(*z1)*(*z2); ++z1; ++z2;
       res+=std::conj(*z1)*(*z2); ++z1; ++z2;
       res+=std::conj(*z1)*(*z2); ++z1; ++z2;
       *rp=res; ++rp;}}
}


           //  Evaluates sum_a  q1[a](x,t)  q2[a](x,t)   a = color index
           //  (This could easily be threaded, but the fast version will run on the device)

void doColorVectorContract(LattField& result, const LattField& q1, const LattField& q2)
{
 if  ((q1.getFieldSiteType()!=FieldSiteType::ColorVector)
    ||(q2.getFieldSiteType()!=FieldSiteType::ColorVector)){
    errorLaph("doColorContract only works with ColorVector lattice fields");}
 result.reset(FieldSiteType::Complex);

 int s1=q1.get_cpu_prec_bytes();
 int s2=q2.get_cpu_prec_bytes();
 if (s1!=s2){
    errorLaph("doColorContract requires lattice fields have same precision");}
 int nsites=LayoutInfo::getRankLatticeNumSites();
 if (s1==sizeof(complex<double>)){
    const complex<double>* z1=reinterpret_cast<const complex<double>*>(q1.getDataConstPtr());
    const complex<double>* z2=reinterpret_cast<const complex<double>*>(q2.getDataConstPtr());
    complex<double>* rp=reinterpret_cast<complex<double>*>(result.getDataPtr());
    for (int k=0;k<nsites;++k){
       complex<double> res=(*z1)*(*z2); ++z1; ++z2;
       res+=(*z1)*(*z2); ++z1; ++z2;
       res+=(*z1)*(*z2); ++z1; ++z2;
       *rp=res; ++rp;}}
 else{
    const complex<float>* z1=reinterpret_cast<const complex<float>*>(q1.getDataConstPtr());
    const complex<float>* z2=reinterpret_cast<const complex<float>*>(q2.getDataConstPtr());
    complex<float>* rp=reinterpret_cast<complex<float>*>(result.getDataPtr());
    for (int k=0;k<nsites;++k){
       complex<float> res=(*z1)*(*z2); ++z1; ++z2;
       res+=(*z1)*(*z2); ++z1; ++z2;
       res+=(*z1)*(*z2); ++z1; ++z2;
       *rp=res; ++rp;}}
}

           // Evaluates   sum_a,b,c epsilson(a,b,c)  q1[a](x,t) q2[b](x,t) q3[c](x,t) 
           //  (This could easily be threaded, but the fast version will run on the device)

void doColorContract(LattField& result, const LattField& q1, const LattField& q2, 
                     const LattField& q3)
{
 if  ((q1.getFieldSiteType()!=FieldSiteType::ColorVector)
    ||(q2.getFieldSiteType()!=FieldSiteType::ColorVector)
    ||(q3.getFieldSiteType()!=FieldSiteType::ColorVector)){
    errorLaph("doColorContract only works with ColorVector lattice fields");}
 result.reset(FieldSiteType::Complex);
 int s1=q1.get_cpu_prec_bytes();
 int s2=q2.get_cpu_prec_bytes();
 int s3=q3.get_cpu_prec_bytes();
 if ((s1!=s2)||(s2!=s3)){
    errorLaph("doColorContract requires lattice fields have same precision");}
 int nsites=LayoutInfo::getRankLatticeNumSites();
 if (s1==sizeof(complex<double>)){
    const complex<double>* z1a=reinterpret_cast<const complex<double>*>(q1.getDataConstPtr());
    const complex<double>* z2a=reinterpret_cast<const complex<double>*>(q2.getDataConstPtr());
    const complex<double>* z3a=reinterpret_cast<const complex<double>*>(q3.getDataConstPtr());
    const complex<double>* z1b=z1a+1;  const complex<double>* z1c=z1b+1;
    const complex<double>* z2b=z2a+1;  const complex<double>* z2c=z2b+1;
    const complex<double>* z3b=z3a+1;  const complex<double>* z3c=z3b+1;
    complex<double>* rp=reinterpret_cast<complex<double>*>(result.getDataPtr());
    for (int k=0;k<nsites;++k){
       *rp=(*z1a)*((*z2b)*(*z3c)-(*z2c)*(*z3b))
          +(*z1b)*((*z2c)*(*z3a)-(*z2a)*(*z3c))
          +(*z1c)*((*z2a)*(*z3b)-(*z2b)*(*z3a));
       z1a+=3; z2a+=3; z3a+=3; z1b+=3; z2b+=3; z3b+=3; 
       z1c+=3; z2c+=3; z3c+=3; ++rp;}}
 else{
    const complex<float>* z1a=reinterpret_cast<const complex<float>*>(q1.getDataConstPtr());
    const complex<float>* z2a=reinterpret_cast<const complex<float>*>(q2.getDataConstPtr());
    const complex<float>* z3a=reinterpret_cast<const complex<float>*>(q3.getDataConstPtr());
    const complex<float>* z1b=z1a+1;  const complex<float>* z1c=z1b+1;
    const complex<float>* z2b=z2a+1;  const complex<float>* z2c=z2b+1;
    const complex<float>* z3b=z3a+1;  const complex<float>* z3c=z3b+1;
    complex<float>* rp=reinterpret_cast<complex<float>*>(result.getDataPtr());
    for (int k=0;k<nsites;++k){
       *rp=(*z1a)*((*z2b)*(*z3c)-(*z2c)*(*z3b))
          +(*z1b)*((*z2c)*(*z3a)-(*z2a)*(*z3c))
          +(*z1c)*((*z2a)*(*z3b)-(*z2b)*(*z3a));
       z1a+=3; z2a+=3; z3a+=3; z1b+=3; z2b+=3; z3b+=3; 
       z1c+=3; z2c+=3; z3c+=3; ++rp;}}
}


           // Evaluates   sum_a,b epsilson(a,b,c)  q1[a](x,t) q2[b](x,t)
           //  (This could easily be threaded, but the fast version will run on the device)

void doColorCrossProduct(LattField& result, const LattField& q1, const LattField& q2)
{
 if  ((q1.getFieldSiteType()!=FieldSiteType::ColorVector)
    ||(q2.getFieldSiteType()!=FieldSiteType::ColorVector)){
    errorLaph("doColorContract only works with ColorVector lattice fields");}
 result.reset(FieldSiteType::ColorVector);
 int s1=q1.get_cpu_prec_bytes();
 int s2=q2.get_cpu_prec_bytes();
 if (s1!=s2){
    errorLaph("doColorContract requires lattice fields have same precision");}
 int nsites=LayoutInfo::getRankLatticeNumSites();
 if (s1==sizeof(complex<double>)){
    const complex<double>* z1a=reinterpret_cast<const complex<double>*>(q1.getDataConstPtr());
    const complex<double>* z2a=reinterpret_cast<const complex<double>*>(q2.getDataConstPtr());
    const complex<double>* z1b=z1a+1;  const complex<double>* z1c=z1b+1;
    const complex<double>* z2b=z2a+1;  const complex<double>* z2c=z2b+1;
    complex<double>* rp=reinterpret_cast<complex<double>*>(result.getDataPtr());
    for (int k=0;k<nsites;++k){
       *rp=(*z1b)*(*z2c)-(*z1c)*(*z2b); ++rp;
       *rp=(*z1c)*(*z2a)-(*z1a)*(*z2c); ++rp;
       *rp=(*z1a)*(*z2b)-(*z1b)*(*z2a); ++rp;
       z1a+=3; z2a+=3; z1b+=3; z2b+=3; z1c+=3; z2c+=3;}}
 else{
    const complex<float>* z1a=reinterpret_cast<const complex<float>*>(q1.getDataConstPtr());
    const complex<float>* z2a=reinterpret_cast<const complex<float>*>(q2.getDataConstPtr());
    const complex<float>* z1b=z1a+1;  const complex<float>* z1c=z1b+1;
    const complex<float>* z2b=z2a+1;  const complex<float>* z2c=z2b+1;
    complex<float>* rp=reinterpret_cast<complex<float>*>(result.getDataPtr());
    for (int k=0;k<nsites;++k){
       *rp=(*z1b)*(*z2c)-(*z1c)*(*z2b); ++rp;
       *rp=(*z1c)*(*z2a)-(*z1a)*(*z2c); ++rp;
       *rp=(*z1a)*(*z2b)-(*z1b)*(*z2a); ++rp;
       z1a+=3; z2a+=3; z1b+=3; z2b+=3; z1c+=3; z2c+=3;}}
}



   // Sets the zeroth element at a site to exp(I*sum(site[k]*momfactors[k],k))
   // and each subsequent element at the site gets multiplied by another
   // factor of exp(I*local_phase)

void setVariablePhaseField(LattField& field, const std::vector<double>& momfactors,
                           const double& local_phase)
{
 bool dp=(field.bytesPerWord()==sizeof(complex<double>));
 vector<char> sitedata(field.bytesPerSite());
 vector<int> start(LayoutInfo::Ndim);
 vector<int> stop(LayoutInfo::Ndim);
 for (int d=0;d<LayoutInfo::Ndim;++d){
    start[d]=LayoutInfo::getMyCommCoords()[d]*LayoutInfo::getRankLattExtents()[d];
    stop[d]=start[d]+LayoutInfo::getRankLattExtents()[d];}
 complex<double> localfactor(std::cos(local_phase),std::sin(local_phase));
 vector<int> site(LayoutInfo::Ndim);
 for (site[3]=start[3];site[3]<stop[3];++site[3]){
    double tphase=momfactors[3]*double(site[3]);
    complex<double> tfactor(std::cos(tphase),std::sin(tphase));
    for (site[2]=start[2];site[2]<stop[2];++site[2]){
       double zphase=momfactors[2]*double(site[2]);
       complex<double> zfactor(std::cos(zphase),std::sin(zphase));
       zfactor*=tfactor;
       for (site[1]=start[1];site[1]<stop[1];++site[1]){
          double yphase=momfactors[1]*double(site[1]);
          complex<double> yfactor(std::cos(yphase),std::sin(yphase));
          yfactor*=zfactor;
          for (site[0]=start[0];site[0]<stop[0];++site[0]){
             double xphase=momfactors[0]*double(site[0]);
             complex<double> factor(std::cos(xphase),std::sin(xphase));
             factor*=yfactor;
             double* dptr=reinterpret_cast<double*>(sitedata.data());
             float* sptr=reinterpret_cast<float*>(sitedata.data());
             for (int j=0;j<int(field.elemsPerSite());++j,++dptr,++sptr){
                if (dp){ *dptr=factor.real(); ++dptr; *dptr=factor.imag();}
                else { *sptr=float(factor.real()); ++sptr; *sptr=float(factor.imag());}
                factor*=localfactor;}
             field.putSiteData(site,sitedata);}}}}
}


void compare_latt_fields(const LattField& src1, const LattField& src2)
{
 bool flag=true;
 int n1=src1.elemsPerSite();
 int n2=src2.elemsPerSite();
 double df=0.0;
 if (n1!=n2){ flag=false;}
 else{
    int s1=src1.get_cpu_prec_bytes();
    int s2=src2.get_cpu_prec_bytes();
    if (s1!=s2){flag=false;}
    else{
       if (s1==sizeof(complex<double>)){
          const complex<double>* z1=reinterpret_cast<const complex<double>*>(src1.getDataConstPtr());
          const complex<double>* z2=reinterpret_cast<const complex<double>*>(src2.getDataConstPtr());
          int ncomp=n1*LayoutInfo::getRankLatticeNumSites();
          for (int k=0;k<ncomp;++k,++z1,++z2){
             double dff=std::abs((*z1)-(*z2));
             if (dff>df){df=dff;}}}
       else{
          const complex<float>* z1=reinterpret_cast<const complex<float>*>(src1.getDataConstPtr());
          const complex<float>* z2=reinterpret_cast<const complex<float>*>(src2.getDataConstPtr());
          int ncomp=n1*LayoutInfo::getRankLatticeNumSites();
          for (int k=0;k<ncomp;++k,++z1,++z2){
             double dff=double(std::abs((*z1)-(*z2)));
             if (dff>df){df=dff;}}}}}
 bool gflag=globalAnd(flag);
 double gdf=globalMax(df);
 if (gflag){ printLaph(make_strf("Fields AGREE within %g",gdf));}
 else printLaph("Fields DIFFER in precision or site type");
}



template <typename T>
void su3color_mult(complex<T>* prod, const complex<T>* cdata1, const complex<T>* cdata2,
                   int krange, int kstride, int istride, int jstride)
{
 for (int k=0;k<krange;++k)
 for (int i=0;i<FieldNcolor;++i){
    complex<T> z(0.0,0.0);
    for (int j=0;j<FieldNcolor;++j){
       z+=cdata1[FieldNcolor*i+j]*cdata2[jstride*j+kstride*k];}
    prod[istride*i+kstride*k]=z;}
}

template <typename T>
void su3color_adjmult(complex<T>* prod, const complex<T>* cdata1, const complex<T>* cdata2,
                      int krange, int kstride, int istride, int jstride)
{
 for (int k=0;k<krange;++k)
 for (int i=0;i<FieldNcolor;++i){
    complex<T> z(0.0,0.0);
    for (int j=0;j<FieldNcolor;++j){
       z+=std::conj(cdata1[FieldNcolor*j+i])*cdata2[jstride*j+kstride*k];}
    prod[istride*i+kstride*k]=z;}
}

    //  Lattice-site-wise color-matrix multiplies outfield = fieldL * fieldR.
    //  fieldL must be of type color-matrix. Any spin indices go along untouched.

void su3color_multiplier(LattField& outfield, const LattField& fieldL, const LattField& fieldR,
                         char Lmat)
{
 if (fieldL.getFieldSiteType()!=FieldSiteType::ColorMatrix){
    errorLaph("left field must be ColorMatrix in su3color_mult");}
 if (fieldR.getFieldSiteType()==FieldSiteType::Complex){
    errorLaph("right field in su3color_mult cannot be complex (non-color) field");}
 outfield.reset(fieldR.getFieldSiteType());
 if (fieldR.bytesPerWord()!=fieldL.bytesPerWord()){
    errorLaph("su3color_mult requires same precision between left and right fields");}
 int oinc=outfield.elemsPerSite();
 int rinc=fieldR.elemsPerSite();
 int linc=fieldL.elemsPerSite();
 int nsites=LayoutInfo::getRankLatticeNumSites();
 int kextent=0, kstride=0, jstride=0, istride=0;
 if (fieldR.getFieldSiteType()==FieldSiteType::ColorMatrix){
    kstride=1; kextent=3; istride=3; jstride=3;}
 else if (fieldR.getFieldSiteType()==FieldSiteType::ColorVector){
    kstride=1; kextent=1; istride=1; jstride=1;}
 else if (fieldR.getFieldSiteType()==FieldSiteType::ColorSpinVector){
    kstride=3; kextent=4; istride=1; jstride=1;}
 if (fieldR.bytesPerWord()==sizeof(complex<double>)){
    complex<double>* op=reinterpret_cast<complex<double>*>(outfield.getDataPtr());
    const complex<double>* rp=reinterpret_cast<const complex<double>*>(fieldR.getDataConstPtr());
    const complex<double>* lp=reinterpret_cast<const complex<double>*>(fieldL.getDataConstPtr());
    void (*multfunc)(complex<double>*,const complex<double>*,const complex<double>*,int,int,int,int)
         = (Lmat=='m')? &su3color_mult<double> : &su3color_adjmult<double>;
    for (int ind=0;ind<nsites;++ind,op+=oinc,rp+=rinc,lp+=linc){
       multfunc(op,lp,rp,kextent,kstride,istride,jstride);}}
 else{
    complex<float>* op=reinterpret_cast<complex<float>*>(outfield.getDataPtr());
    const complex<float>* rp=reinterpret_cast<const complex<float>*>(fieldR.getDataConstPtr());
    const complex<float>* lp=reinterpret_cast<const complex<float>*>(fieldL.getDataConstPtr());
    void (*multfunc)(complex<float>*,const complex<float>*,const complex<float>*,int,int,int,int)
         = (Lmat=='m')? &su3color_mult<float> : &su3color_adjmult<float>;
    for (int ind=0;ind<nsites;++ind,op+=oinc,rp+=rinc,lp+=linc){
       multfunc(op,lp,rp,kextent,kstride,istride,jstride);}}
}


    //  Lattice-site-wise color-matrix multiplies outfield = fieldL * fieldR.
    //  fieldL must be of type color-matrix. Any spin indices go along untouched.

void su3color_mult(LattField& outfield, const LattField& fieldL, const LattField& fieldR)
{
 su3color_multiplier(outfield,fieldL,fieldR,'m');
}

    //  Lattice-site-wise color-matrix multiplies outfield = colorAdj(fieldL) * fieldR.
    //  fieldL must be of type color-matrix.  Any spin indices go along untouched.

void su3color_adjmult(LattField& outfield, const LattField& fieldL, const LattField& fieldR)
{
 su3color_multiplier(outfield,fieldL,fieldR,'a');
}



   // assigns y[d] from x[d]; if d==dir, y[d]=(x[d]+1) % N[d]

void index_increaser(int dir, int d, vector<int>& y, const vector<int>& x, 
                     const vector<int>& N)
{
 if (dir==d){ y[d]=(x[d]==(N[d]-1))?0:x[d]+1;}
 else{ y[d]=x[d];}
}

   // assigns y[d] from x[d]; if d==dir, y[d]=(x[d]-1+N[d]) % N[d]

void index_decreaser(int dir, int d, vector<int>& y, const vector<int>& x, 
                     const vector<int>& N)
{
 if (dir==d){ y[d]=(x[d]==0)?N[d]-1:x[d]-1;}
 else{ y[d]=x[d];}
}


void flipsign(char* sitedata, int bps, bool dp)
{
 if (dp){
    int nelem=bps/sizeof(complex<double>);
    complex<double>* op=reinterpret_cast<complex<double>*>(sitedata);
    for (int k=0;k<nelem;++k,++op){
       *op=-(*op);}}
 else{
    int nelem=bps/sizeof(complex<float>);
    complex<float>* op=reinterpret_cast<complex<float>*>(sitedata);
    for (int k=0;k<nelem;++k,++op){
       *op=-(*op);}}
}


#ifdef ARCH_SERIAL


      //  This routine is slow: meant only for debugging, testing
      //   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
      //   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

      // linear index for site (x,y,z,t) is
      //      (x+Nx*(y+Ny*(z+Nz*t)))/2 + (nsites/2)*((x+y+z+t)%2);

void latt_shifter(LattField& outfield, const LattField& infield, int dir, 
                  void (*indexfunc)(int,int,vector<int>&,const vector<int>&, 
                                    const vector<int>&), char fwd_or_bwd,
                  bool apply_antiperiodic_fermbc=false)
{
 if ((dir<0)||(dir>3)){
    errorLaph("Invalid direction in lattice field shift");}
 outfield.reset(infield.getFieldSiteType());
 const vector<int>& N=LayoutInfo::getLattExtents();
 int nsites=LayoutInfo::getLatticeNumSites();
 int npsites=nsites/2;
 int bps=infield.bytesPerSite();
 char* oute=outfield.getDataPtr();
 char* outo=oute+npsites*bps;
 const char* ine=infield.getDataConstPtr();
 const char* ino=ine+npsites*bps;
 int fliptime=-1;
 if ((infield.getFieldSiteType()==FieldSiteType::ColorSpinVector)
     && (apply_antiperiodic_fermbc) && (dir==3)){
    fliptime=(fwd_or_bwd=='B')?0:N[3]-1;}
 bool dp=(infield.bytesPerWord()==sizeof(complex<double>)); 
 vector<int> x(LayoutInfo::Ndim);
 vector<int> y(LayoutInfo::Ndim);
 for (x[3]=0;x[3]<N[3];++x[3]){
    indexfunc(dir,3,y,x,N);
    bool signflip=(x[3]==fliptime);
    for (x[2]=0;x[2]<N[2];++x[2]){
       indexfunc(dir,2,y,x,N);
       for (x[1]=0;x[1]<N[1];++x[1]){
          indexfunc(dir,1,y,x,N);
          for (x[0]=0;x[0]<N[0];++x[0]){
             indexfunc(dir,0,y,x,N);
             int shift=bps*((y[0]+N[0]*(y[1]+N[1]*(y[2]+N[2]*y[3])))/2);
             if ((x[0]+x[1]+x[2]+x[3])%2){
                const char* inp=ine+shift;
                std::memcpy(outo,inp,bps);
                if (signflip){ flipsign(outo,bps,dp);}
                outo+=bps;}
             else{
                const char* inp=ino+shift;
                std::memcpy(oute,inp,bps);
                if (signflip){ flipsign(oute,bps,dp);}
                oute+=bps;}}}}}
}

      //  This routine is slow: meant only for debugging, testing
      //   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
      //   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

void lattice_shift(LattField& outfield, const LattField& infield, int dir, char fwd_or_bwd,
                   bool apply_antiperiodic_fermbc=false)
{
 if (fwd_or_bwd=='F'){
    latt_shifter(outfield,infield,dir,&index_increaser,fwd_or_bwd,apply_antiperiodic_fermbc);}
 else if (fwd_or_bwd=='B'){
    latt_shifter(outfield,infield,dir,&index_decreaser,fwd_or_bwd,apply_antiperiodic_fermbc);}
 else{
    errorLaph("lattice shift needs F or B for forward/backward");}
}


#else                // parallel version now


bool fwd_stay_local(int dir, const vector<int>& x, const vector<int>& N, bool nocomm)
{
 if (x[dir]<(N[dir]-1)) return true;
 return nocomm;
}

bool bwd_stay_local(int dir, const vector<int>& x, const vector<int>& N, bool nocomm)
{
 if (x[dir]>0) return true;
 return nocomm;
}


void get_fwd_comm_to_from(int dir, int& send_to, int& recv_from)
{
 vector<int> rank_coord(LayoutInfo::getMyCommCoords());  // rank coords of this node
       // get mpi rank where to send
 if (rank_coord[dir]>0) rank_coord[dir]--;
 else rank_coord[dir]=LayoutInfo::getCommNumPartitions()[dir]-1;
 send_to=LayoutInfo::getRankFromCommCoords(rank_coord);
       // get mpi rank where to receive from
 rank_coord=LayoutInfo::getMyCommCoords();
 rank_coord[dir]++;
 if (rank_coord[dir]==LayoutInfo::getCommNumPartitions()[dir]) rank_coord[dir]=0;
 recv_from=LayoutInfo::getRankFromCommCoords(rank_coord);
}

void get_bwd_comm_to_from(int dir, int& send_to, int& recv_from)
{
 vector<int> rank_coord(LayoutInfo::getMyCommCoords());  // rank coords of this node
       // get mpi rank where to receive from
 if (rank_coord[dir]>0) rank_coord[dir]--;
 else rank_coord[dir]=LayoutInfo::getCommNumPartitions()[dir]-1;
 recv_from=LayoutInfo::getRankFromCommCoords(rank_coord);
       // get mpi rank where to send
 rank_coord=LayoutInfo::getMyCommCoords();
 rank_coord[dir]++;
 if (rank_coord[dir]==LayoutInfo::getCommNumPartitions()[dir]) rank_coord[dir]=0;
 send_to=LayoutInfo::getRankFromCommCoords(rank_coord);
}


      //  This routine is slow: meant only for debugging, testing
      //   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
      //   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

      // linear index for site (x,y,z,t) is
      //      (x+Nx*(y+Ny*(z+Nz*t)))/2 + (nsites/2)*((x+y+z+t)%2);


void latt_shifter(LattField& outfield, const LattField& infield, int dir, 
                  void (*indexfunc)(int,int,vector<int>&,const vector<int>&, 
                                    const vector<int>&),
                  bool (*staylocal)(int,const vector<int>&,const vector<int>&,bool),
                  void (*neighbors)(int,int&,int&), char fwd_or_bwd,
                  bool apply_antiperiodic_fermbc=false)
{
 int start_parity=LayoutInfo::getMyStartParity();
 if ((dir<0)||(dir>3)){
    errorLaph("Invalid direction in lattice field shift");}
 outfield.reset(infield.getFieldSiteType());
 const vector<int>& N=LayoutInfo::getRankLattExtents();
 int nsites=LayoutInfo::getRankLatticeNumSites();
 int npsites=nsites/2;
 int bps=infield.bytesPerSite();
 char* oute=outfield.getDataPtr();
 char* outo=oute+npsites*bps;
 const char* ine=infield.getDataConstPtr();
 const char* ino=ine+npsites*bps;
 vector<int> x(LayoutInfo::Ndim);
 vector<int> y(LayoutInfo::Ndim);
 bool needcomm=(LayoutInfo::getCommNumPartitions()[dir]>1);
 bool nocm=!needcomm;
 vector<char> sendbuffer,recvbuffer;
 char *snde=0;
 char *sndo=0;
 int nbuf=0;
 bool nospflip=((N[dir]%2)==0);
 if (needcomm){
    nbuf=nsites/N[dir]; nbuf+=nbuf%2;
    sendbuffer.resize(nbuf*bps);
    recvbuffer.resize(nbuf*bps);
    snde=sendbuffer.data();
    sndo=snde+(nbuf/2)*bps;
    if (!nospflip){
       sndo=sendbuffer.data();
       snde=sndo+(nbuf/2)*bps;}}
 int fliptime=-1;
 if ((infield.getFieldSiteType()==FieldSiteType::ColorSpinVector)
     && (apply_antiperiodic_fermbc) && (dir==3)){
    if ((fwd_or_bwd=='B')&&(LayoutInfo::getMyCommCoords()[3]==(LayoutInfo::getCommNumPartitions()[3]-1))){
       fliptime=0;}
    else if ((fwd_or_bwd=='F')&&(LayoutInfo::getMyCommCoords()[3]==0)){
       fliptime=N[3]-1;}}
 bool dp=(infield.bytesPerWord()==sizeof(complex<double>)); 
        // make local changes and put data in send buffer
 for (x[3]=0;x[3]<N[3];++x[3]){
    indexfunc(dir,3,y,x,N);
    bool signflip=(x[3]==fliptime);
    for (x[2]=0;x[2]<N[2];++x[2]){
       indexfunc(dir,2,y,x,N);
       for (x[1]=0;x[1]<N[1];++x[1]){
          indexfunc(dir,1,y,x,N);
          for (x[0]=0;x[0]<N[0];++x[0]){
             indexfunc(dir,0,y,x,N);
             int site_parity=(start_parity+x[0]+x[1]+x[2]+x[3])%2;
             int shift=bps*((y[0]+N[0]*(y[1]+N[1]*(y[2]+N[2]*y[3])))/2);
             if (site_parity){
                if (staylocal(dir,x,N,nocm)){
                   const char* inp=ine+shift; 
                   std::memcpy(outo,inp,bps);
                   if (signflip){ flipsign(outo,bps,dp);}}
                else{
                   const char* inp=(nospflip?ine:ino)+shift;
                   std::memcpy(sndo,inp,bps);
                   if (signflip){ flipsign(sndo,bps,dp);}
                   sndo+=bps;}
                outo+=bps;}
             else{
                if (staylocal(dir,x,N,nocm)){
                   const char* inp=ino+shift; 
                   std::memcpy(oute,inp,bps);
                   if (signflip){ flipsign(oute,bps,dp);}}
                else{
                   const char* inp=(nospflip?ino:ine)+shift;
                   std::memcpy(snde,inp,bps);
                   if (signflip){ flipsign(snde,bps,dp);}
                   snde+=bps;}
                oute+=bps;}}}}}
 if (!needcomm){
    return;}

      // do the communication (send and receive) 
 int send_to=0;
 int recv_from=0;
 neighbors(dir,send_to,recv_from);
 MPI_Status mpistatus;
 int status=MPI_Sendrecv(sendbuffer.data(), nbuf*bps, MPI_BYTE, send_to, LayoutInfo::getMyRank(), 
                         recvbuffer.data(), nbuf*bps, MPI_BYTE, recv_from, recv_from, 
                         MPI_COMM_WORLD, &mpistatus);
 if (status!=MPI_SUCCESS){
    errorLaph("Error in communication while doing shift");}

     // put data received into local memory appropriately
 sendbuffer.clear();
 oute=outfield.getDataPtr();
 outo=oute+npsites*bps;
 const char* rcve=recvbuffer.data();
 const char* rcvo=rcve+(nbuf/2)*bps;
 for (x[3]=0;x[3]<N[3];++x[3]){
    for (x[2]=0;x[2]<N[2];++x[2]){
       for (x[1]=0;x[1]<N[1];++x[1]){
          for (x[0]=0;x[0]<N[0];++x[0]){
             if ((start_parity+x[0]+x[1]+x[2]+x[3])%2){
                if (!staylocal(dir,x,N,nocm)){ std::memcpy(outo,rcvo,bps); rcvo+=bps;}
                outo+=bps;}
             else{
                if (!staylocal(dir,x,N,nocm)){ std::memcpy(oute,rcve,bps); rcve+=bps;}
                oute+=bps;}}}}}
}


     //  This routine is slow: meant only for debugging, testing
      //   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
      //   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

void lattice_shift(LattField& outfield, const LattField& infield, int dir, char fwd_or_bwd,
                   bool apply_antiperiodic_fermbc=false)
{
 if (fwd_or_bwd=='F'){
    latt_shifter(outfield,infield,dir,&index_increaser,&fwd_stay_local,&get_fwd_comm_to_from,
                 fwd_or_bwd,apply_antiperiodic_fermbc);}
 else if (fwd_or_bwd=='B'){
    latt_shifter(outfield,infield,dir,&index_decreaser,&bwd_stay_local,&get_bwd_comm_to_from,
                 fwd_or_bwd,apply_antiperiodic_fermbc);}
 else{
    errorLaph("lattice shift needs F or B for forward/backward");}
 MPI_Barrier(MPI_COMM_WORLD);
}


#endif 
    //   outfield(x) <=  U_dir(x) infield(x+mu)            if fwd_or_bwd=='F'
    //   outfield(x) <=  U^dag_dir(x-mu) infield(x-mu)     if fwd_or_bwd=='B'

    //   for 'F', su3mult( U[dir], shift(infield, mu, 'F') )
    //   for 'B', shift(  su3mult( adj(U[dir]), infield ), mu, 'B')

void lattice_cov_shift(LattField& outfield, const LattField& infield, 
                       const vector<LattField>& gauge_field,
                       int dir, char fwd_or_bwd,
                       bool apply_antiperiodic_fermbc=false)
{
 if (int(gauge_field.size())!=LayoutInfo::Ndim){
    errorLaph("invalid gauge field in lattice_cov_shift");}
 for (uint dir=0;dir<gauge_field.size();++dir){
    if (gauge_field[dir].getFieldSiteType()!=FieldSiteType::ColorMatrix){
       errorLaph("invalid gauge field in lattice_cov_shift");}}
 if (infield.getFieldSiteType()==FieldSiteType::Complex){
    errorLaph("cannot covariantly shift a complex (non-color) field");}
 LattField tmp;
 if (fwd_or_bwd=='F'){
    lattice_shift(tmp,infield,dir,'F',apply_antiperiodic_fermbc);
    su3color_mult(outfield,gauge_field[dir],tmp);}
 else if (fwd_or_bwd=='B'){
    su3color_adjmult(tmp,gauge_field[dir],infield);
    lattice_shift(outfield,tmp,dir,'B',apply_antiperiodic_fermbc);}
 else{
    errorLaph("lattice shift needs F or B for forward/backward");}
}


     //  covariant shift by a path of directions: each direction is an integer
     //  having values    1,2,3,4 (forward directions) and -1,-2,-3,-4 (backwards)   
     //  where 1 means x, 2 means y, 3 means z, and 4 means t

void lattice_cov_shift(LattField& outfield, const LattField& infield, 
                       const vector<LattField>& gauge_field,
                       vector<int>& path)
{
 LattField tmp;
 LattField *p1=0,*p2=0, *sw=0;
 const LattField *pc;
 if (path.size()%2){ p2=&outfield; pc=&infield;}
 else{ p2=&tmp; pc=&infield;}
 for (int seg=int(path.size())-1;seg>=0;--seg){
    lattice_cov_shift(*p2,*pc,gauge_field,abs(path[seg])-1,path[seg]>0?'F':'B');
    if (seg==int(path.size())-1){
       if (path.size()%2){ p1=&tmp;}
       else{ p1=&outfield;}}
    sw=p1; p1=p2; p2=sw; pc=p1;}
}


template <typename T>
void su3color_adjcopy(complex<T>* out, const complex<T>* in)
{
 for (int i=0;i<FieldNcolor;++i)
 for (int j=0;j<FieldNcolor;++j){
    out[FieldNcolor*i+j]=std::conj(in[FieldNcolor*j+i]);}
}


void su3color_adjcopy(LattField& outfield, const LattField& infield)
{
 if (infield.getFieldSiteType()!=FieldSiteType::ColorMatrix){
    errorLaph("field must be ColorMatrix in su3color_adjcopy");}
 outfield.reset(FieldSiteType::ColorMatrix);
 int inc=infield.elemsPerSite();
 int nsites=LayoutInfo::getRankLatticeNumSites();
 if (infield.bytesPerWord()==sizeof(complex<double>)){
    complex<double>* op=reinterpret_cast<complex<double>*>(outfield.getDataPtr());
    const complex<double>* ip=reinterpret_cast<const complex<double>*>(infield.getDataConstPtr());
    for (int ind=0;ind<nsites;++ind,op+=inc,ip+=inc){
       su3color_adjcopy(op,ip);}}
 else{
    complex<float>* op=reinterpret_cast<complex<float>*>(outfield.getDataPtr());
    const complex<float>* ip=reinterpret_cast<const complex<float>*>(infield.getDataConstPtr());
    for (int ind=0;ind<nsites;++ind,op+=inc,ip+=inc){
       su3color_adjcopy(op,ip);}}
}



void lattice_link_path(LattField& outfield, const vector<LattField>& gaugefield,
                       vector<int>& dir_path)
{
 uint nshifts=dir_path.size();
 if (nshifts==0){
    errorLaph("lattice_link_path requires a non-trivial path");}
 LattField temp;
 vector<int>::const_reverse_iterator seg=dir_path.rbegin();
 bool doprod=false;
 for (uint k=0;k<dir_path.size();++k){
    int dir=*seg;
    if (dir>0){
       if (doprod){
          lattice_shift(temp,outfield,dir-1,'F');
          su3color_mult(outfield,gaugefield[dir-1],temp);}
       else{
          outfield=gaugefield[dir-1];
          doprod=true;}}
    else{
       if (doprod){
          su3color_adjmult(temp,gaugefield[-dir-1],outfield);}
       else{
          su3color_adjcopy(temp,gaugefield[-dir-1]);
          doprod=true;}
       lattice_shift(outfield,temp,-dir-1,'B');}
    ++seg;}
}


void lattice_addto(LattField& outfield, const LattField& infield, 
                   const complex<double>& zcoef=complex<double>(1.0,0.0))
{
 if (outfield.getFieldSiteType()!=infield.getFieldSiteType()){
    errorLaph("fields must be of same type to add");}
 if (outfield.bytesPerWord()!=infield.bytesPerWord()){
    errorLaph("addition requires same precision of fields");}
 int nelem=LayoutInfo::getRankLatticeNumSites()*infield.elemsPerSite();
 if (infield.bytesPerWord()==sizeof(complex<double>)){
    complex<double>* op=reinterpret_cast<complex<double>*>(outfield.getDataPtr());
    const complex<double>* ip=reinterpret_cast<const complex<double>*>(infield.getDataConstPtr());
    for (int k=0;k<nelem;++k,++op,++ip){
       *op+=zcoef*(*ip);}}
  else{
    complex<float>* op=reinterpret_cast<complex<float>*>(outfield.getDataPtr());
    const complex<float>* ip=reinterpret_cast<const complex<float>*>(infield.getDataConstPtr());
    complex<float> zf(float(zcoef.real()),float(zcoef.imag()));
    for (int k=0;k<nelem;++k,++op,++ip){
       *op+=zf*(*ip);}}
}

   //  *op = Gamma(spin_matrix_index) * (*ip)  site-wise
   //
   //        spin_matrix_index     matrix  (DeGrand-Rossi basis)
   //               1              gamma[1]
   //               2              gamma[2]     sigma[mu,nu] = i/2 [gamma[mu], gamma[nu]]
   //               3              gamma[3]
   //               4              gamma[4]
   //               5              gamma[5] = gamma[4]*gamma[1]*gamma[2]*gamma[3]
   //               6              sigma[1,2]
   //               7              sigma[1,3]
   //               8              sigma[1,4]
   //               9              sigma[2,3]
   //              10              sigma[2,4]
   //              11              sigma[3,4]

template <typename T>
void spin_mult(complex<T>* op, const complex<T>* ip, const vector<pair<int,complex<T>>>& spinmat)
{
 int spinstride=FieldNcolor;
 complex<T> *opp=op;
 const complex<T> *ipp=ip;
 for (int color=0;color<FieldNcolor;++color,++opp,++ipp){
    for (int outspin=0;outspin<FieldNspin;++outspin){
       *(opp+outspin*spinstride)=*(ipp+spinmat[outspin].first*spinstride)*spinmat[outspin].second;}}
}


   //  Spin matrices in the DeGrand-Rossi basis
template <typename T>
void assign_spin_matrix(int spin_matrix_index, vector<pair<int,complex<T>>>& spin_mat)
{
 spin_mat.resize(4);
 complex<T> I(0.0,1.0);
 complex<T> one(1.0,0.0);
 if (spin_matrix_index==1){
    spin_mat[0]=pair<int,complex<T>>(3,I); spin_mat[1]=pair<int,complex<T>>(2,I);
    spin_mat[2]=pair<int,complex<T>>(1,-I); spin_mat[3]=pair<int,complex<T>>(0,-I);}
 else if (spin_matrix_index==2){
    spin_mat[0]=pair<int,complex<T>>(3,-one); spin_mat[1]=pair<int,complex<T>>(2,one);
    spin_mat[2]=pair<int,complex<T>>(1,one); spin_mat[3]=pair<int,complex<T>>(0,-one);}
 else if (spin_matrix_index==3){
    spin_mat[0]=pair<int,complex<T>>(2,I); spin_mat[1]=pair<int,complex<T>>(3,-I);
    spin_mat[2]=pair<int,complex<T>>(0,-I); spin_mat[3]=pair<int,complex<T>>(1,I);}
 else if (spin_matrix_index==4){
    spin_mat[0]=pair<int,complex<T>>(2,one); spin_mat[1]=pair<int,complex<T>>(3,one);
    spin_mat[2]=pair<int,complex<T>>(0,one); spin_mat[3]=pair<int,complex<T>>(1,one);}
 else if (spin_matrix_index==5){
    spin_mat[0]=pair<int,complex<T>>(0,-one); spin_mat[1]=pair<int,complex<T>>(1,-one);
    spin_mat[2]=pair<int,complex<T>>(2,one); spin_mat[3]=pair<int,complex<T>>(3,one);}
 else if (spin_matrix_index==6){
    spin_mat[0]=pair<int,complex<T>>(0,one); spin_mat[1]=pair<int,complex<T>>(1,-one);
    spin_mat[2]=pair<int,complex<T>>(2,one); spin_mat[3]=pair<int,complex<T>>(3,-one);}
 else if (spin_matrix_index==7){
    spin_mat[0]=pair<int,complex<T>>(1,-I); spin_mat[1]=pair<int,complex<T>>(0,I);
    spin_mat[2]=pair<int,complex<T>>(3,-I); spin_mat[3]=pair<int,complex<T>>(2,I);}
 else if (spin_matrix_index==8){
    spin_mat[0]=pair<int,complex<T>>(1,-one); spin_mat[1]=pair<int,complex<T>>(0,-one);
    spin_mat[2]=pair<int,complex<T>>(3,one); spin_mat[3]=pair<int,complex<T>>(2,one);}
 else if (spin_matrix_index==9){
    spin_mat[0]=pair<int,complex<T>>(1,one); spin_mat[1]=pair<int,complex<T>>(0,one);
    spin_mat[2]=pair<int,complex<T>>(3,one); spin_mat[3]=pair<int,complex<T>>(2,one);}
 else if (spin_matrix_index==10){
    spin_mat[0]=pair<int,complex<T>>(1,-I); spin_mat[1]=pair<int,complex<T>>(0,I);
    spin_mat[2]=pair<int,complex<T>>(3,I); spin_mat[3]=pair<int,complex<T>>(2,-I);}
 else if (spin_matrix_index==11){
    spin_mat[0]=pair<int,complex<T>>(0,-one); spin_mat[1]=pair<int,complex<T>>(1,one);
    spin_mat[2]=pair<int,complex<T>>(2,one); spin_mat[3]=pair<int,complex<T>>(3,-one);}
 else{
    errorLaph("Unsupported Dirac gamma spin index");}
}


   //  Multiply infield by a spin matrices in the DeGrand-Rossi basis, returning
   //  result in outfield.
void lattice_spin_multiply(LattField& outfield, const LattField& infield, int spin_matrix_index)
{
 if (infield.getFieldSiteType()!=FieldSiteType::ColorSpinVector){
    errorLaph("spin_multiply can only be done on a ColorSpinVector");}
 outfield.reset(FieldSiteType::ColorSpinVector);
 int nsites=LayoutInfo::getRankLatticeNumSites();
 int nelem=infield.elemsPerSite();
 if (infield.bytesPerWord()==sizeof(complex<double>)){
    vector<pair<int,complex<double>>> spin_mat(4);
    assign_spin_matrix<double>(spin_matrix_index,spin_mat);
    complex<double>* op=reinterpret_cast<complex<double>*>(outfield.getDataPtr());
    const complex<double>* ip=reinterpret_cast<const complex<double>*>(infield.getDataConstPtr());
    for (int k=0;k<nsites;++k,op+=nelem,ip+=nelem){
       spin_mult<double>(op,ip,spin_mat);}}
  else{
    vector<pair<int,complex<float>>> spin_mat(4);
    assign_spin_matrix<float>(spin_matrix_index,spin_mat);
    complex<float>* op=reinterpret_cast<complex<float>*>(outfield.getDataPtr());
    const complex<float>* ip=reinterpret_cast<const complex<float>*>(infield.getDataConstPtr());
    for (int k=0;k<nsites;++k,op+=nelem,ip+=nelem){
       spin_mult<float>(op,ip,spin_mat);}}
}


void get_spin_transform_matrix(Array<double>& transmat, const std::string& from_to)
{
 double vcp=1.0/sqrt(2.0);
 double vcm=-1.0/sqrt(2.0);
 if (from_to=="GR_to_DP"){
    transmat(0,0)=0.0; transmat(0,1)=vcm; transmat(0,2)=0.0; transmat(0,3)=vcm;
    transmat(1,0)=vcp; transmat(1,1)=0.0; transmat(1,2)=vcp; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=vcp; transmat(2,2)=0.0; transmat(2,3)=vcm;
    transmat(3,0)=vcm; transmat(3,1)=0.0; transmat(3,2)=vcp; transmat(3,3)=0.0;}
 else if (from_to=="DP_to_GR"){
    transmat(0,0)=0.0; transmat(0,1)=vcp; transmat(0,2)=0.0; transmat(0,3)=vcm;
    transmat(1,0)=vcm; transmat(1,1)=0.0; transmat(1,2)=vcp; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=vcp; transmat(2,2)=0.0; transmat(2,3)=vcp;
    transmat(3,0)=vcm; transmat(3,1)=0.0; transmat(3,2)=vcm; transmat(3,3)=0.0;}
 else if (from_to=="UK_to_GR"){
    transmat(0,0)=0.0; transmat(0,1)=vcm; transmat(0,2)=0.0; transmat(0,3)=vcm;
    transmat(1,0)=vcp; transmat(1,1)=0.0; transmat(1,2)=vcp; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=vcm; transmat(2,2)=0.0; transmat(2,3)=vcp;
    transmat(3,0)=vcp; transmat(3,1)=0.0; transmat(3,2)=vcm; transmat(3,3)=0.0;}
 else if (from_to=="GR_to_UK"){
    transmat(0,0)=0.0; transmat(0,1)=vcp; transmat(0,2)=0.0; transmat(0,3)=vcp;
    transmat(1,0)=vcm; transmat(1,1)=0.0; transmat(1,2)=vcm; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=vcp; transmat(2,2)=0.0; transmat(2,3)=vcm;
    transmat(3,0)=vcm; transmat(3,1)=0.0; transmat(3,2)=vcp; transmat(3,3)=0.0;}
 else if ((from_to=="CH_to_GR")||(from_to=="GR_to_CH")){
    transmat(0,0)=0.0; transmat(0,1)=0.0; transmat(0,2)=0.0; transmat(0,3)=-1.0;
    transmat(1,0)=0.0; transmat(1,1)=0.0; transmat(1,2)=1.0; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=1.0; transmat(2,2)=0.0; transmat(2,3)=0.0;
    transmat(3,0)=-1.0;transmat(3,1)=0.0; transmat(3,2)=0.0; transmat(3,3)=0.0;}
 else if (from_to=="UK_to_CH"){
    transmat(0,0)=vcm; transmat(0,1)=0.0; transmat(0,2)=vcp; transmat(0,3)=0.0;
    transmat(1,0)=0.0; transmat(1,1)=vcm; transmat(1,2)=0.0; transmat(1,3)=vcp;
    transmat(2,0)=vcp; transmat(2,1)=0.0; transmat(2,2)=vcp; transmat(2,3)=0.0;
    transmat(3,0)=0.0; transmat(3,1)=vcp; transmat(3,2)=0.0; transmat(3,3)=vcp;}
 else if (from_to=="CH_to_UK"){
    transmat(0,0)=vcm; transmat(0,1)=0.0; transmat(0,2)=vcp; transmat(0,3)=0.0;
    transmat(1,0)=0.0; transmat(1,1)=vcm; transmat(1,2)=0.0; transmat(1,3)=vcp;
    transmat(2,0)=vcp; transmat(2,1)=0.0; transmat(2,2)=vcp; transmat(2,3)=0.0;
    transmat(3,0)=0.0; transmat(3,1)=vcp; transmat(3,2)=0.0; transmat(3,3)=vcp;}
 else if ((from_to=="DP_to_UK")||(from_to=="UK_to_DP")){
    transmat(0,0)=-1.0;transmat(0,1)=0.0; transmat(0,2)=0.0; transmat(0,3)=0.0;
    transmat(1,0)=0.0; transmat(1,1)=-1.0;transmat(1,2)=0.0; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=0.0; transmat(2,2)=1.0; transmat(2,3)=0.0;
    transmat(3,0)=0.0; transmat(3,1)=0.0; transmat(3,2)=0.0; transmat(3,3)=1.0;}
}

template <typename T>
void convertSpinBasis(complex<T>* zGR, const complex<T>* zDP,
                      const Array<double>& transmat)
{
 complex<T>* pL=zGR;
 const complex<T>* pR=zDP;
 for (int c=0;c<3;c++,++pL,++pR){
    complex<T>* ppL=pL;
    for (int s=0;s<4;++s,ppL+=3){
       complex<T> res=complex<T>(0.0,0.0);
       const complex<T>* ppR=pR;
       for (int ss=0;ss<4;++ss,ppR+=3){
          res+=T(transmat(s,ss))*(*ppR);}
       *ppL=res;}}
}

void convertSpinBasis(std::vector<char>& siteGR, const std::vector<char>& siteDP, bool dp,
                      const Array<double>& transmat)
{
 if (dp){
    complex<double>* zGR=reinterpret_cast<complex<double>*>(siteGR.data());
    const complex<double>* zDP=reinterpret_cast<const complex<double>*>(siteDP.data());
    convertSpinBasis<double>(zGR,zDP,transmat);}
 else{
    complex<float>* zGR=reinterpret_cast<complex<float>*>(siteGR.data());
    const complex<float>* zDP=reinterpret_cast<const complex<float>*>(siteDP.data());
    convertSpinBasis<float>(zGR,zDP,transmat);}
}


void convertSpinBasis(LattField& outferm, const LattField& inferm, const std::string& from_to)
{
 bool dp=(inferm.bytesPerWord()==sizeof(complex<double>)); 
 //int nelem_site=inferm.elemsPerSite();
 Array<double> transmat(4,4);
 get_spin_transform_matrix(transmat,from_to);
 std::vector<char> siteDataOut(inferm.bytesPerSite());
 std::vector<int> coord(LayoutInfo::Ndim);
 for (coord[3]=0;coord[3]<LayoutInfo::getLattExtents()[3];++coord[3])
 for (coord[2]=0;coord[2]<LayoutInfo::getLattExtents()[2];++coord[2])
 for (coord[1]=0;coord[1]<LayoutInfo::getLattExtents()[1];++coord[1])
 for (coord[0]=0;coord[0]<LayoutInfo::getLattExtents()[0];++coord[0]){
    std::vector<char> siteDataIn(inferm.getSiteData(coord));
    convertSpinBasis(siteDataOut,siteDataIn,dp,transmat);
    outferm.putSiteData(coord,siteDataOut);}
}


void calcCloverLeaves(LattField& cloverleaf, const vector<LattField>& gaugefield,
                      int dir1, int dir2)
{
 LattField Utmp;
 vector<int> path(4);
      // upper right leaf
 path[0]=dir1; path[1]=dir2; path[2]=-dir1; path[3]=-dir2;
 lattice_link_path(cloverleaf,gaugefield,path);
      // upper left leaf
 path[0]=dir2; path[1]=-dir1; path[2]=-dir2; path[3]=dir1;
 lattice_link_path(Utmp,gaugefield,path);
 lattice_addto(cloverleaf,Utmp);
      // lower left leaf
 path[0]=-dir1; path[1]=-dir2; path[2]=dir1; path[3]=dir2;
 lattice_link_path(Utmp,gaugefield,path);
 lattice_addto(cloverleaf,Utmp);
      // lower right leaf
 path[0]=-dir2; path[1]=dir1; path[2]=dir2; path[3]=-dir1;
 lattice_link_path(Utmp,gaugefield,path);
 lattice_addto(cloverleaf,Utmp);
}


   //   Applies the clover Dirac operation to "infield", returning
   //   the result in "outfield".  This operation is
   //
   //    outfield =  [ 1/(2*kappa) - (1/2) Dterm  + CFterm ] infield
   //
   //   where
   //
   //        Dterm = sum_mu [ (1-gamma_mu) U  + (1+gamma_mu) * U^dag ]
   //
   //        CFterm = csw (i/4)  sigma[mu,nu] F[mu,nu]
   //
   //            sigma[mu,nu] (i/2) [gamma_mu, gamma_nu]
   //
   //            F[mu,nu] = (1/8) ( Q[mu,nu]-Q[nu,mu] )
   //
   //            Q[mu,nu] = U(mu,nu,-mu,-nu) + U(nu,-mu,-nu,mu)
   //                     + U(-mu,-nu,mu,nu) + U(-nu,mu,nu,-mu)
   //
   //   Lattice shifts must take the fermion temporal boundary
   //   conditions into account.
   //
   //   The routine below assumes the DeGrand-Rossi Dirac matrix convention.
   

void applyCloverDiracDGR(LattField& outfield, const LattField& infield,
                         const vector<LattField>& gauge_field,
                         const GaugeConfigurationInfo& gaction,
                         const QuarkActionInfo& qaction)
{
 if (int(gauge_field.size())!=LayoutInfo::Ndim){
    errorLaph("invalid gauge field in applyCloverDirac");}
 for (uint dir=0;dir<gauge_field.size();++dir){
    if (gauge_field[dir].getFieldSiteType()!=FieldSiteType::ColorMatrix){
       errorLaph("invalid gauge field in applyCloverDirac");}}
 if (infield.getFieldSiteType()!=FieldSiteType::ColorSpinVector){
    errorLaph("can apply clover Dirac operator only to a color-spin vector field");}
 bool tbc1=qaction.isFermionTimeBCAntiPeriodic();
 bool tbc2=gaction.isFermionTimeBCAntiPeriodic();
 if (tbc1!=tbc2){
    errorLaph("Inconsistent fermion time boundary conditions in QuarkActionInfo and GaugeConfigurationInfo",true);}
 outfield.reset(FieldSiteType::ColorSpinVector);

 if (qaction.getName()!="WILSON_CLOVER"){
    errorLaph("Can only applyCloverDirac if action name is WILSON_CLOVER");}
 //bool timebc_antiperiodic=qaction.isFermionTimeBCAntiPeriodic();
 bool timebc_antiperiodic=false;       // minus sign is stored in gauge field both on host and device
 double kappa=qaction.getRValues()[2];
 double csw=qaction.getRValues()[4];
 double anisotropy=qaction.getRValues()[3];
 if (std::abs(anisotropy-1.0)>1e-12){
    errorLaph("Current applyCloverDirac only applies for isotropic actions");}

    //  Now for the clover Dirac operator acting on a color-spin field
 LattField cloverleaf(FieldSiteType::ColorMatrix);
 LattField cloverleaf2(FieldSiteType::ColorMatrix);
 LattField phi(FieldSiteType::ColorSpinVector);
 LattField phi2(FieldSiteType::ColorSpinVector);
 LattField CFterm(FieldSiteType::ColorSpinVector);
 LattField Dterm(FieldSiteType::ColorSpinVector);
 setZeroField(CFterm);
 setZeroField(Dterm);
 int spin_index=6;
 complex<double> cfclover(0.0,csw/16.0);
 for (int dir1=1;dir1<=LayoutInfo::Ndim;++dir1)
 for (int dir2=dir1+1;dir2<=LayoutInfo::Ndim;++dir2,++spin_index){
    calcCloverLeaves(cloverleaf,gauge_field,dir1,dir2);
    calcCloverLeaves(cloverleaf2,gauge_field,dir2,dir1);
    lattice_addto(cloverleaf,cloverleaf2,complex<double>(-1.0,0.0));
    su3color_mult(phi,cloverleaf,infield);
    lattice_spin_multiply(phi2,phi,spin_index);
    lattice_addto(CFterm,phi2,cfclover);}

    //  Now the leading derivative and the Wilson term
 for (int dir=1;dir<=4;++dir){
    lattice_cov_shift(phi,infield,gauge_field,dir-1,'F',timebc_antiperiodic);
    lattice_addto(Dterm,phi);
    lattice_spin_multiply(phi2,phi,dir);
    lattice_addto(Dterm,phi2,complex<double>(-1.0,0.0));
    lattice_cov_shift(phi,infield,gauge_field,dir-1,'B',timebc_antiperiodic);
    lattice_addto(Dterm,phi);
    lattice_spin_multiply(phi2,phi,dir);
    lattice_addto(Dterm,phi2);}

 setZeroField(outfield);
 lattice_addto(outfield,infield,complex<double>(1.0/(2.0*kappa),0.0));
 lattice_addto(outfield,Dterm,complex<double>(-0.5,0.0));
 lattice_addto(outfield,CFterm);
}


   //   This routine assumes the Dirac-Pauli convention for the Dirac gamma-matrices.
   
void applyCloverDirac(LattField& outfield, const LattField& infield,
                      const vector<LattField>& gauge_field,
                      const GaugeConfigurationInfo& gaction,
                      const QuarkActionInfo& qaction)
{
 LattField tmp(FieldSiteType::ColorSpinVector);
 convertSpinBasis(outfield,infield,"DP_to_GR");
 applyCloverDiracDGR(tmp,outfield,gauge_field,gaction,qaction);
 convertSpinBasis(outfield,tmp,"GR_to_DP");
}
 

    //  Applies the 3d spatial Laplacian with a smeared gauge field
    //  onto "infield", returning result in "outfield".  These fields
    //  must be color vectors

void applyMinusSpatialLaplacian(LattField& outfield, const LattField& infield,
                                const vector<LattField>& smeared_gauge_field)
{
 if (int(smeared_gauge_field.size())!=LayoutInfo::Ndim){
    errorLaph("invalid number of components in smeared gauge field in applySpatialLaplacian");}
 for (uint dir=0;dir<smeared_gauge_field.size();++dir){
    if (smeared_gauge_field[dir].getFieldSiteType()!=FieldSiteType::ColorMatrix){
       errorLaph(make_strf("invalid field type in smeared gauge field[%d] in applySpatialLaplacian",dir));}}
 if (infield.getFieldSiteType()!=FieldSiteType::ColorVector){
    errorLaph("can apply spatial laplacian operator only to a color vector field");}
 outfield.reset(FieldSiteType::ColorVector);

 setZeroField(outfield);
 lattice_addto(outfield,infield,complex<double>(6.0,0.0));
 LattField phi(FieldSiteType::ColorVector);

    //  Apply the spatial covariant shifts
 for (int dir=1;dir<=3;++dir){
    lattice_cov_shift(phi,infield,smeared_gauge_field,dir-1,'F');
    lattice_addto(outfield,phi,complex<double>(-1.0,0.0));
    lattice_cov_shift(phi,infield,smeared_gauge_field,dir-1,'B');
    lattice_addto(outfield,phi,complex<double>(-1.0,0.0));}
}


// **************************************************************************
}

