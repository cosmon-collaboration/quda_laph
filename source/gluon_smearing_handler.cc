#include "gluon_smearing_handler.h"
#include "host_global.h"
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
     : smearPtr(0), uPtr(0), m_read_mode(true) {}



GluonSmearingHandler::GluonSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                                           const GaugeConfigurationInfo& gauge,
                                           const string& smearedFieldFileName,
                                           bool readmode)
{
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
 h_filename=tidyNameOfFile(smearedFieldFileName);
 try{
    uPtr=new GaugeConfigurationInfo(gauge);
    smearPtr=new GluonSmearingInfo(gluon_smearing);}
 catch(const std::exception& xp){
    printLaph("problem in setInfo in GluonSmearingHandler");
    throw(std::runtime_error("No smeared gauge file found"));}
}

    // destructor

GluonSmearingHandler::~GluonSmearingHandler()
{
 clear();
}


void GluonSmearingHandler::clear()
{
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

           // Compute the smeared gauge field into HostGlobal and put into
           // the file whose name is stored in "h_filename" unless the
           // file name is empty

void GluonSmearingHandler::computeSmearedGaugeField()
{
 check_info_set("computeSmearedGaugeField");
 if (m_read_mode)
    throw(std::runtime_error("Cannot computeSmearedGaugeField in read mode"));

 if (GB::theSmearedGaugeConfigIsSet()){
    return;}
 
    // check if output file exists; if so, we do not want to overwrite
    // so return

 if ((!h_filename.empty())&&(fileExists(h_filename))){
    printLaph(make_strf("file %s already exists",h_filename));
    printLaph(" will not overwrite; no output written\n");
    return;}

 printLaph("Computation of smeared gauge field commencing in GluonSmearingHandler");
 printLaph(smearPtr->output());

    // get the gauge configuration and put onto the device (config will remain
    // in HostGlobal and on the device after this handler is destroyed)
 GaugeConfigurationHandler GH(*uPtr);
 GH.setData();
 XMLHandler gauge_xmlinfo;
 GH.getXMLInfo(gauge_xmlinfo);
 printLaph("XML info for the gauge configuration:");
 printLaph(make_strf("%s\n",gauge_xmlinfo.output()));
 GH.copyDataToDevice();

 QudaGaugeSmearParam gauge_smear_param = newQudaGaugeSmearParam();
 smearPtr->setQudaGaugeSmearParam(gauge_smear_param);

 int calcobs=gauge_smear_param.n_steps/gauge_smear_param.meas_interval+1;
 vector<QudaGaugeObservableParam> gauge_obs_param(calcobs);
 for (int i=0;i<calcobs;++i){
    gauge_obs_param[i] = newQudaGaugeObservableParam();
    gauge_obs_param[i].compute_plaquette = QUDA_BOOLEAN_TRUE;
    gauge_obs_param[i].compute_qcharge = QUDA_BOOLEAN_TRUE;
    gauge_obs_param[i].su_project = QUDA_BOOLEAN_TRUE;
    gauge_obs_param[i].struct_size = sizeof(gauge_obs_param[i]);}

 StopWatch bulova; bulova.start();
 performGaugeSmearQuda(&gauge_smear_param, gauge_obs_param.data());

 bulova.stop();
 printLaph(make_strf("Smearing computation done: time was %g seconds",bulova.getTimeInSeconds()));
 QudaInfo::smeared_gauge_on_device=true;

    // Now save smeared gauge to host
 bulova.reset(); bulova.start();
 QudaGaugeParam smeared_gauge_info = newQudaGaugeParam();
 uPtr->setQudaGaugeParam(smeared_gauge_info);
 smeared_gauge_info.type = QUDA_SMEARED_LINKS;

 GB::theSmearedGaugeConfig.resize(LayoutInfo::Ndim,FieldSiteType::ColorMatrix);
 std::vector<char*> gaugeptrs(LayoutInfo::Ndim);
 for (int k=0;k<LayoutInfo::Ndim;++k) gaugeptrs[k]=GB::theSmearedGaugeConfig[k].getDataPtr();
 saveGaugeQuda(gaugeptrs.data(), &smeared_gauge_info);

     //  Now output to file if requested
 if (!h_filename.empty()){
    DataPutHandlerSF<GluonSmearingHandler,UIntKey,vector<LattField>> 
               dhput(*this,h_filename,"Laph--SmearedGaugeField4D");
    UIntKey dummykey(0);
    dhput.putData(dummykey,GB::theSmearedGaugeConfig);
    printLaph("Write to file complete");}

 bulova.stop();
 printLaph(make_strf("Time to write to host and then to file was %g seconds",bulova.getTimeInSeconds()));
 printLaph("\n\nSmeared gauge field computation done");
 printLaph(make_strf("Output file: %s\n",h_filename));
 QudaInfo::clearDeviceSmearedGaugeConfiguration();   // remove from device
}


const vector<LattField>& GluonSmearingHandler::getSmearedGaugeField()
{
 if (!GB::theSmearedGaugeConfigIsSet()){
    if (!loadSmearedGaugeField()){
       errorLaph("Was unable to get the smeared gauge field");}}
 return GB::theSmearedGaugeConfig;
}

bool GluonSmearingHandler::querySmearedGaugeField()
{
 if (!GB::theSmearedGaugeConfigIsSet()){
    return loadSmearedGaugeField();}
 return true;
}

   // load from file into HostGlobal

bool GluonSmearingHandler::loadSmearedGaugeField()
{
 check_info_set("loadSmearedGaugeField");
 if (!m_read_mode) return false;
 if (GB::theSmearedGaugeConfigIsSet()) return true;
 if (h_filename.empty()) return false;

    // check that the file exists and contains the appropriate data
 printLaph(make_strf("Checking file: %s",h_filename));
 if (!fileExists(h_filename)) return false;

      // Open for reading (reads header info and checks)
 try{
    string ftype("Laph--SmearedGaugeField4D");
    DataGetHandlerSF<GluonSmearingHandler,RecordKey,vector<LattField>> DH(*this,h_filename,ftype);
    printLaph(make_strf("File %s successfully opened and header matches\n",h_filename));
    RecordKey dummykey(0);
    GB::theSmearedGaugeConfig=DH.getData(dummykey);
    return true;}
 catch(const std::exception& xp){
   printLaph("unable to read smeared gauge from file(s)");
   return false;}
 return true;
}


void GluonSmearingHandler::clearData()
{
 GB::theSmearedGaugeConfig.clear();
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
 if (QudaInfo::smeared_gauge_on_device){ return;}

      // Create the array of pointers, and call the load function
 const vector<LattField>& stoutlinks=getSmearedGaugeField();
 vector<const char*> gauge_ptrs(LayoutInfo::Ndim);
 for (int dir=0;dir<LayoutInfo::Ndim;++dir){
    gauge_ptrs[dir]=stoutlinks[dir].getDataConstPtr();}

      // create the quda gauge params
 QudaGaugeParam smeared_gauge_info = newQudaGaugeParam();
 uPtr->setQudaGaugeParam(smeared_gauge_info);
 smeared_gauge_info.type = QUDA_SMEARED_LINKS;

      // copy to the gpu device ( gaugeSmeared is pointer to it)
 loadGaugeQuda((void *)gauge_ptrs.data(), &smeared_gauge_info);
 QudaInfo::smeared_gauge_on_device=true;
}


void GluonSmearingHandler::eraseDataOnDevice()
{
 QudaInfo::clearDeviceSmearedGaugeConfiguration();   // remove from device
}

// ******************************************************************
}
