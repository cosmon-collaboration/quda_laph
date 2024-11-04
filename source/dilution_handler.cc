#include "dilution_handler.h"
#include "laph_stdio.h"
#include "layout_info.h"

using namespace std;

namespace LaphEnv {


// *************************************************************

DilutionHandler::DilutionHandler() : dilPtr(0), Textent(0), //minTime(-1),
                 /*maxTime(-1),*/ nEigvecs(0), nSpinProjectors(0),
                 nEigvecProjectors(0), nTimeProjectors(0), Nspin(4) {}


DilutionHandler::DilutionHandler(const DilutionSchemeInfo& dilScheme,
                                 const QuarkSmearingInfo& qSmear,
                                 bool UpperSpinComponentsOnly)
{
 set_info(dilScheme,qSmear,UpperSpinComponentsOnly);
}


void DilutionHandler::setInfo(const DilutionSchemeInfo& dilScheme,
                              const QuarkSmearingInfo& qSmear,
                              bool UpperSpinComponentsOnly)
{
 clear();
 set_info(dilScheme,qSmear,UpperSpinComponentsOnly);
}


DilutionHandler::~DilutionHandler()
{
 clear();
}


void DilutionHandler::clear()
{
 try{
    delete dilPtr;}
 catch(const std::exception& xp){
    errorLaph("abort");}
 dilPtr=0;
 Textent=0;
// minTime=-1;
// maxTime=-1;
 nEigvecs=0;
 nSpinProjectors=0;
 nEigvecProjectors=0;
 nTimeProjectors=0;
 spinProjs.clear();
 eigvecProjs.clear();
 timeProjs.clear();
 spin_proj_indices.clear();
 eigvec_proj_indices.clear();
 which_time_proj.clear();
 which_spin_proj.clear();
 which_eigvec_proj.clear();
}


void DilutionHandler::set_info(const DilutionSchemeInfo& dilScheme,
                               const QuarkSmearingInfo& qSmear,
                               bool UpperSpinComponentsOnly)
{
 try{
    dilPtr=new DilutionSchemeInfo(dilScheme);}
 catch(const std::exception& xp){
    errorLaph("could not allocate memory for DilutionSchemeInfo");}
 Textent=LayoutInfo::getLattExtents()[3];
// minTime=gaugeinfo.getMinTime();
// maxTime=gaugeinfo.getMaxTime();
 nEigvecs=qSmear.getNumberOfLaplacianEigenvectors(); 
 Nspin=UpperSpinComponentsOnly ? 2 : 4;

 setProjectorMasks(spinProjs, dilPtr->spinDilutionType, Nspin, 
                   which_spin_proj, nSpinProjectors);
 setProjectorMasks(eigvecProjs, dilPtr->eigvecDilutionType, nEigvecs, 
                   which_eigvec_proj, nEigvecProjectors);
 setProjectorMasks(timeProjs, dilPtr->timeDilutionType, Textent, 
                   which_time_proj, nTimeProjectors);

 // Throw out times outside [minTime, maxTime]
 // leaves us with empty dilution projectors; get-routines throw if one
 // of these is actually requested
/* for (vector<list<int>>::iterator tprojIt = timeProjs.begin();
	tprojIt != timeProjs.end(); ++tprojIt) {
   for (list<int>::iterator tIt = tprojIt->begin(); tIt != tprojIt->end();) {
     if (*tIt < minTime || *tIt > maxTime) tIt = tprojIt->erase(tIt);
     else tIt++;
   }
 }*/

 spin_proj_indices.resize(nSpinProjectors*nEigvecProjectors);
 eigvec_proj_indices.resize(nSpinProjectors*nEigvecProjectors);
 int count=0;
 for (int v=0;v<nEigvecProjectors;v++)
    for (int s=0;s<nSpinProjectors;s++){
       spin_proj_indices[count]=s;
       eigvec_proj_indices[count]=v;
       count++;}
}



bool DilutionHandler::isInfoSet() const
{
 return ((dilPtr!=0)&&(Textent>0)&&(nEigvecs>0)  /*&&(minTime>-1)&&(maxTime>-1)*/ );
}


const DilutionSchemeInfo& DilutionHandler::getDilutionSchemeInfo() const
{
 check_info_set("getDilutionSchemeInfo");
 return *dilPtr;
}

int DilutionHandler::getNumberOfSpinEigvecProjectors() const
{
 check_info_set("getNumberOfSpinEigvecProjectors");
 return nSpinProjectors*nEigvecProjectors;
}

int DilutionHandler::getNumberOfSpinProjectors() const
{
 check_info_set("getNumberOfSpinProjectors");
 return nSpinProjectors;
}

int DilutionHandler::getNumberOfEigvecProjectors() const
{
 check_info_set("getNumberOfEigvecProjectors");
 return nEigvecProjectors;
}

int DilutionHandler::getNumberOfTimeProjectors() const
{
 check_info_set("getNumberOfTimeProjectors");
 return nTimeProjectors;
}

int DilutionHandler::getSpinProjectorIndex(int spineigvec_index) const
{
 check_info_set("getSpinProjectorIndex");
 check_valid_spineig(spineigvec_index);
 return spin_proj_indices[spineigvec_index];
}

int DilutionHandler::getEigvecProjectorIndex(int spineigvec_index) const
{
 check_info_set("getEigvecProjectorIndex");
 check_valid_spineig(spineigvec_index);
 return eigvec_proj_indices[spineigvec_index];
}

int DilutionHandler::getTimeProjectorIndex(int time_val) const
{
 check_info_set("getTimeProjectorIndex");
 if ((time_val<0)||(time_val>=Textent)){
 //if ((time_val<minTime)||(time_val>maxTime)){
    throw(std::invalid_argument("bad time value"));}
 return which_time_proj[time_val];
}

const list<int>& DilutionHandler::getOnSpinIndices(int spineigvec_index) const
{
 check_info_set("getOnSpinIndices");
 check_valid_spineig(spineigvec_index);
 return spinProjs[spin_proj_indices[spineigvec_index]];
}

const list<int>& DilutionHandler::getOnEigvecIndices(int spineigvec_index) const
{
 check_info_set("getOnEigvecIndices");
 check_valid_spineig(spineigvec_index);
 return eigvecProjs[eigvec_proj_indices[spineigvec_index]];
}

const list<int>& DilutionHandler::getSourceOnEigvecIndices(int eigvec_index) const
{
 check_info_set("getSourceOnEigvecIndices");
 if ((eigvec_index<0)||(eigvec_index>=nEigvecProjectors))
    throw(std::invalid_argument("bad eigvec_index"));
 return eigvecProjs[eigvec_index];
}

const list<int>& DilutionHandler::getOnTimes(int time_proj_index) const
{
 check_info_set("getOnTimes");
 check_valid_timeproj(time_proj_index);
 return timeProjs[time_proj_index];
}


bool DilutionHandler::isOnSpin(int spineigvec_index, int spin_val) const
{
 check_info_set("isOnSpin");
 check_valid_spineig(spineigvec_index);
 if ((spin_val<0)||(spin_val>=int(Nspin))){
    throw(std::invalid_argument("invalid spin index"));}
 return (which_spin_proj[spin_val]==spin_proj_indices[spineigvec_index]);
}

bool DilutionHandler::isOnEigvec(int spineigvec_index, int eigvec_index) const
{
 check_info_set("isOnEigvec");
 check_valid_spineig(spineigvec_index);
 if ((eigvec_index<0)||(eigvec_index>=nEigvecs)){
    throw(std::invalid_argument("invalid eigenvector index"));}
 return (which_eigvec_proj[eigvec_index]==eigvec_proj_indices[spineigvec_index]);
}

bool DilutionHandler::isOnTime(int time_proj_index, int time_val) const
{
 check_info_set("isOnTime");
 check_valid_timeproj(time_proj_index);
 if ((time_val<0)||(time_val>=Textent)){
 //if ((time_val<minTime)||(time_val>maxTime)){
    throw(std::invalid_argument("invalid time index"));}
 return (which_time_proj[time_val]==time_proj_index);
}


void DilutionHandler::setProjectorMasks(vector<list<int> >& projs, 
                                        int dil_type, int nBasis,
                                        vector<int>& projind, int& nproj)
{
 projs.clear();
 projind.resize(nBasis);
 if (dil_type==0){         // no dilution
    nproj=1;
    for (int k=0;k<nBasis;k++) projind[k]=0;
    }
 else if (dil_type==1){    // full dilution
    nproj=nBasis;
    for (int k=0;k<nBasis;k++) projind[k]=k;
    }
 else if (dil_type>1){     // block dilution
    nproj=dil_type;
    if (nproj>(nBasis/2)){
       throw(std::invalid_argument("number of blocked dilution projectors too large"));}
    int bksize=nBasis/nproj;
    vector<int> blocksize(nproj,bksize);
    int tb=nBasis-bksize*nproj;
    for (int sb=0;sb<tb;sb++)
       blocksize[sb]++;
    int jb=0, bc=0;
    for (int k=0;k<nBasis;k++){
       projind[k]=jb; bc++;
       if (bc==blocksize[jb]){ jb++; bc=0;}}
    }
 else if (dil_type<-1){    // interlace dilution
    nproj=-dil_type;
    if (nproj>(nBasis/2)){
       throw(std::invalid_argument("number of interlaced dilution projectors too large"));}
    int jp=0;
    for (int k=0;k<nBasis;k++){
       projind[k]=jp++;
       if (jp==nproj) jp=0;}
    }

 projs.resize(nproj);
 for (int k=0;k<nBasis;k++)
    projs[projind[k]].push_back(k);
}
 
void DilutionHandler::check_info_set(const string& name) const
{
 if (!isInfoSet()){
    errorLaph(make_strf("error in DilutionHandler: must setInfo before calling %s",name));}
}

void DilutionHandler::check_valid_spineig(int spineigvec_index) const
{
 if ((spineigvec_index<0)||(spineigvec_index>=(nSpinProjectors*nEigvecProjectors)))
    throw(std::invalid_argument("bad spineigvec_index"));
}

void DilutionHandler::check_valid_timeproj(int timeproj_index) const
{
 if ((timeproj_index<0)||(timeproj_index>=nTimeProjectors)||(timeProjs[timeproj_index].empty()))
    throw(std::invalid_argument("bad time projector index"));
}

bool DilutionHandler::isValidTimeProjectorIndex(int timeproj_index) const
{
 if ((timeproj_index<0)||(timeproj_index>=nTimeProjectors)||(timeProjs[timeproj_index].empty())) return false;
 else return true;
}

bool DilutionHandler::isValidSpinEigvecProjectorIndex(int spineigvec_index) const
{
 if ((spineigvec_index<0)||(spineigvec_index>=(nSpinProjectors*nEigvecProjectors)))
    return false;
 else return true;
}

bool DilutionHandler::isFullTimeDilution() const
{
 return (dilPtr->timeDilutionType==1);
}

// *************************************************************
}
