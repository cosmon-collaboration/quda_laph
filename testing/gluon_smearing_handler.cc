#include "gluon_smearing_handler.h"
#include "stop_watch.h"

using namespace std;
using namespace quda;

namespace LaphEnv {


 // **************************************************************
 // *                                                            *
 // *                                                            *
 // *           GluonSmearingHandler implementation              *
 // *                                                            *
 // *                                                            *
 // **************************************************************

   // constructors

GluonSmearingHandler::GluonSmearingHandler() 
     : smearPtr(0), uPtr(0), m_read_mode(true), dh_ptr(0)  {}



GluonSmearingHandler::GluonSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                                           const GaugeConfigurationInfo& gauge,
                                           const string& smearedFieldFileName,
                                           bool readmode)
{
 dh_ptr = 0;
 set_info(gluon_smearing,gauge,smearedFieldFileName,readmode);
}


void GluonSmearingHandler::setInfo(const GluonSmearingInfo& gluon_smearing,
                                   const GaugeConfigurationInfo& gauge,
                                   const string& smearedFieldFileName,
                                   bool readmode)
{
 clear();
 set_info(gluon_smearing,gauge,smearedFieldFileName,readmode);
}


void GluonSmearingHandler::set_info(const GluonSmearingInfo& gluon_smearing,
                                    const GaugeConfigurationInfo& gauge,
                                    const string& smearedFieldFileName,
                                    bool readmode)
{
 m_read_mode=readmode;
 string ftype("Laph--SmearedGaugeField4D");
 h_filename=tidyString(smearedFieldFileName);
 if (emptyFileName(h_filename)){
    errorLaph("empty file name in GluonSmearingHandler");}
 try{
    uPtr=new GaugeConfigurationInfo(gauge);
    smearPtr=new GluonSmearingInfo(gluon_smearing);}
 catch(const std::exception& xp){
    printLaph("problem in setInfo in GluonSmearingHandler");
    throw(std::runtime_error("No smeared gauge file found"));}

    // check that the file exists and contains the appropriate data
 if (m_read_mode){
    printLaph(make_strf("Checking file: %s",h_filename));
    if (!fileExists(h_filename)){
         string mesg("file "); mesg+=h_filename+" does not exist\n";
         mesg+=" reading the smeared gauge time slices\n";
         mesg+="   will not be possible\n";
         mesg+="No smeared gauge file found\n";
         clear(); 
         throw(std::runtime_error(mesg));}}

      // Open for reading (reads header info and checks)

 if (m_read_mode){
    dh_ptr = new DataGetHandlerSFO<GluonSmearingHandler,RecordKey,
                                   vector<LattField> >(*this,h_filename,ftype);
    printLaph(make_strf("File %s successfully opened and header matches\n",h_filename));}
}

    // destructor

GluonSmearingHandler::~GluonSmearingHandler()
{
 clear();
}


void GluonSmearingHandler::clear()
{
 clearData();
 if (dh_ptr) delete dh_ptr; 
 dh_ptr=0;
 try{
    delete smearPtr;
    delete uPtr;} 
 catch(const std::exception& xp){
    errorLaph("abort");}
 h_filename.clear();
 smearPtr=0;
 uPtr=0;
}



 // **********************************************************

bool GluonSmearingHandler::isInfoSet() const
{
 return ((smearPtr!=0)&&(uPtr!=0));
}

void GluonSmearingHandler::check_info_set(const string& name) const
{
 if (!isInfoSet()){
    string mesg("error in GluonSmearingHandler:");
    mesg+="  must setInfo before calling "+name;
    errorLaph(mesg);}
}

const GluonSmearingInfo& GluonSmearingHandler::getGluonSmearingInfo() const
{
 check_info_set("getGluonSmearingInfo");
 return *smearPtr;
}

const GaugeConfigurationInfo& GluonSmearingHandler::getGaugeConfigurationInfo() const
{
 check_info_set("getGaugeConfigurationInfo");
 return *uPtr;
}

void GluonSmearingHandler::filefail(const string& message)
{
 clear();
 errorLaph(message);
}


 // **********************************************************

           // Compute the smeared gauge field and put into the file
           // whose name is stored in "h_filename".  Use file name
           // starting with "NOM_" to place in the NamedObjMap.

void GluonSmearingHandler::computeSmearedGaugeField()
{
 check_info_set("computeSmearedGaugeField");
 if (m_read_mode)
    throw(std::runtime_error("Cannot computeSmearedGaugeField in read mode"));

    // check if output file exists; if so, we do not want to overwrite
    // so quit

 if (fileExists(h_filename)){
    printLaph(make_strf("file %s already exists",h_filename));
    printLaph(" will not overwrite; no output written\n");
    return;}

 printLaph("Computation of smeared gauge field commencing in GluonSmearingHandler");
 printLaph(smearPtr->output());

 QudaGaugeSmearParam gauge_smear_param = newQudaGaugeSmearParam();
 smearPtr->setQudaGaugeSmearParam(gauge_smear_param);

 int nsmear=gauge_smear_param.n_steps/gauge_smear_param.meas_interval+1;
 vector<QudaGaugeObservableParam> gauge_obs_param(nsmear);
 for (int i=0;i<nsmear;++i){
    gauge_obs_param[i] = newQudaGaugeObservableParam();
    gauge_obs_param[i].compute_plaquette = QUDA_BOOLEAN_FALSE;   // QUDA_BOOLEAN_TRUE;
    gauge_obs_param[i].compute_qcharge = QUDA_BOOLEAN_FALSE;     // QUDA_BOOLEAN_TRUE;
    gauge_obs_param[i].su_project = QUDA_BOOLEAN_TRUE;}

 StopWatch bulova; bulova.start();
 performGaugeSmearQuda(&gauge_smear_param, gauge_obs_param.data());

// printLaph(make_strf(" Computed mean plaquette = %15.12f",gauge_obs_param[nsmear-1].plaquette[0]));
// printLaph(make_strf("  spatial mean plaquette = %15.12f",gauge_obs_param[nsmear-1].plaquette[1]));
// printLaph(make_strf(" temporal mean plaquette = %15.12f",gauge_obs_param[nsmear-1].plaquette[2]));

 bulova.stop();
 printLaph(make_strf("Smearing computation done: time was %g seconds",bulova.getTimeInSeconds()));
 QudaInfo::smeared_gauge_on_device=true;

    // Now save to host
 bulova.reset(); bulova.start();
 QudaGaugeParam smeared_gauge_info = newQudaGaugeParam();
 uPtr->setQudaGaugeParam(smeared_gauge_info);
 smeared_gauge_info.type = QUDA_SMEARED_LINKS;

 std::vector<LattField> stoutlinks(LayoutInfo::Ndim,FieldSiteType::ColorMatrix);
 std::vector<char*> gaugeptrs(LayoutInfo::Ndim);
 for (int k=0;k<LayoutInfo::Ndim;++k) gaugeptrs[k]=stoutlinks[k].getDataPtr();
 saveGaugeQuda(gaugeptrs.data(), &smeared_gauge_info);
 printLaph("Smeared links saved to host");

     //  Now output to file or the NamedObjMap
 DataPutHandlerSFO<GluonSmearingHandler,UIntKey,vector<LattField>> 
           dhput(*this,h_filename,"Laph--SmearedGaugeField4D");
 UIntKey dummykey(0);
 dhput.putData(dummykey,stoutlinks);
 bulova.stop();
 printLaph("Write to file complete");
 printLaph(make_strf("Time to write to host and then to file was %g seconds",bulova.getTimeInSeconds()));

 printLaph("\n\nSmeared gauge field computation done");
 printLaph(make_strf("Output file: %s\n",h_filename));
}


const vector<LattField>& GluonSmearingHandler::getSmearedGaugeField() const
{
 if (!m_read_mode)
    throw(std::runtime_error("Cannot getSmearedGaugeField in write mode"));
 check_info_set("setSmearedGaugeFieldTimeSlice");
 RecordKey dummykey(0);
 return dh_ptr->getData(dummykey);
}


void GluonSmearingHandler::clearData()
{
 if (dh_ptr) dh_ptr->clearData();
}


void GluonSmearingHandler::writeHeader(XMLHandler& xmlout)
{
 xmlout.set_root("SmearedGaugeField");
 XMLHandler xmlh;
 smearPtr->output(xmlh);
 xmlout.put_child(xmlh);
 uPtr->output(xmlh);
 xmlout.put_child(xmlh);
}

bool GluonSmearingHandler::checkHeader(XMLHandler& xmlr)
{
 GaugeConfigurationInfo gauge_check(xmlr);
 GluonSmearingInfo gsmear_check(xmlr);
 try { uPtr->checkEqual(gauge_check); 
       smearPtr->checkEqual(gsmear_check);}
 catch(const std::exception& xp){ 
    return false;}
 return true;
}

void GluonSmearingHandler::copyDataToDevice()
{
 check_info_set("copyDataToDevice");
 if (QudaInfo::smeared_gauge_on_device) return;

      // Create the array of pointers, and call the load function
 const vector<LattField>& stoutlinks=getSmearedGaugeField();
 vector<const char*> gauge_ptrs(LayoutInfo::Ndim);
 for (int dir=0;dir<LayoutInfo::Ndim;++dir){
    gauge_ptrs[dir]=stoutlinks[dir].getDataConstPtr();}

      // create the quda gauge params
 QudaGaugeParam smeared_gauge_info = newQudaGaugeParam();
 uPtr->setQudaGaugeParam(smeared_gauge_info);
 smeared_gauge_info.type = QUDA_SMEARED_LINKS;

 smeared_gauge_info.type = QUDA_WILSON_LINKS;     //   JUST FOR TESTING PURPOSES: CHANGE CHANGE


      // copy to the gpu device ( gaugeSmeared is pointer to it)
 loadGaugeQuda((void *)gauge_ptrs.data(), &smeared_gauge_info);
 QudaInfo::smeared_gauge_on_device=true;
}


void GluonSmearingHandler::eraseDataOnDevice()
{
 if (QudaInfo::smeared_gauge_on_device){
    freeGaugeSmearedQuda();
    QudaInfo::smeared_gauge_on_device=false;}
}

// ******************************************************************
}
