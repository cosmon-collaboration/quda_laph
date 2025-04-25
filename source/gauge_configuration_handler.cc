#include "gauge_configuration_handler.h"
#include "host_global.h"
#include "laph_stdio.h"
#include "read_gauge_field.h"

using namespace std;

namespace LaphEnv {

 // **********************************************************

GaugeConfigurationHandler::GaugeConfigurationHandler()
     : m_gauge_info(0), m_cfg(0) {}


GaugeConfigurationHandler::GaugeConfigurationHandler(
                      const GaugeConfigurationInfo& gaugeinfo)
     : m_gauge_info(0), m_cfg(0)
{
 set_info(gaugeinfo);
}


void GaugeConfigurationHandler::setInfo(const GaugeConfigurationInfo& gaugeinfo)
{
 clear();
 set_info(gaugeinfo);
}

void GaugeConfigurationHandler::set_info(const GaugeConfigurationInfo& gaugeinfo)
{
 try{
    m_gauge_info = new GaugeConfigurationInfo(gaugeinfo);}
 catch(...){
    errorLaph("problem allocating GaugeConfigurationInfo");}    
}

void GaugeConfigurationHandler::setData()
{
 check_info_set("setData");
      // check first to see if configuration has already been read
 if (GB::theGaugeConfigIsSet(*m_gauge_info)){
    m_cfg = &GB::getGaugeConfig(*m_gauge_info); }
 else{      // if not already in the HostGlobal or data in HostGlobal
    initialize_config();}          //  has different info, then must read it now
}


GaugeConfigurationHandler::~GaugeConfigurationHandler()
{
 clear();
}

     // clears this handler, but leaves config in HostGlobal

void GaugeConfigurationHandler::clear()
{
 try {
    delete m_gauge_info;} 
 catch(const std::exception& xp) {
    errorLaph("Fatal error in GaugeConfigurationHandler");}
 m_gauge_info=0;
 m_cfg=0;
}

     // removes config from HostGlobal
     
void GaugeConfigurationHandler::eraseDataOnHost()
{
 GB::clearGaugeConfig();
}
 

const std::vector<LattField>& GaugeConfigurationHandler::getData() 
{
 check_info_set("getData");
 if (!isDataSet()) setData();
 return *m_cfg;
}


const GaugeConfigurationInfo& GaugeConfigurationHandler::getGaugeConfigurationInfo() const
{
 check_info_set("getGaugeConfigurationInfo");
 return *m_gauge_info;
}


void GaugeConfigurationHandler::check_info_set(const string& name) const
{
 if (!isInfoSet()){
    printLaph("error in GaugeConfigurationHandler:");
    errorLaph(make_strf("  must setInfo before calling %s",name));}
}


void GaugeConfigurationHandler::getXMLInfo(XMLHandler& gauge_xmlinfo) const
{
 check_info_set("getXMLInfo");
 m_gauge_info->output(gauge_xmlinfo);
}


void GaugeConfigurationHandler::initialize_config()
{
 try{
    std::vector<LattField>& the_config(GB::setGaugeConfig(*m_gauge_info));
    GaugeConfigReader GR;
    XMLHandler gauge_xmlinfo;
    GR.read(the_config,gauge_xmlinfo,*m_gauge_info);
    m_cfg=&the_config;
        // if the fermion field is antiperiodic in time, then multiply the temporal
        // links on the last time slice by -1.  This does not change the plaquette
        // nor any spatial smearing or hadron operators defined on single time slices.
    if (m_gauge_info->isFermionTimeBCAntiPeriodic()){
       applyFermionTemporalAntiPeriodic();}}
 catch(const std::exception& xp){
    errorLaph(make_strf("Gauge configuration initialization failed: %s",xp.what()));}
}

      // multiplies host temporal links on time Nt-1 by -1 to help with fermion 
      // antiperiodic temporal boundary condition

void GaugeConfigurationHandler::applyFermionTemporalAntiPeriodic()
{
 const int Tdir=LayoutInfo::Ndim-1;
 check_info_set("getData");
 if (!isDataSet()){
    errorLaph("Gauge configuration must already be in the HostGlobal to apply fermion temporal antiperiodic b.c.");}
 try{
    std::vector<LattField>& the_config(GB::setGaugeConfig(*m_gauge_info));
    the_config[Tdir].applyFermionTemporalAntiPeriodic();}
 catch(const std::exception& xp){
    errorLaph(make_strf("applyFermionTemporalAntiPeriodic failed: %s",xp.what()));}
}
 

    // load the configuration onto the gpu device

void GaugeConfigurationHandler::copyDataToDevice(bool removeFromHost)
{
 check_info_set("copyDataToDevice");
 if (QudaInfo::gauge_config_on_device) return;

      // Create the array of pointers, and call the load function
 vector<const char*> gauge_ptrs(LayoutInfo::Ndim);
 for (int dir=0;dir<LayoutInfo::Ndim;++dir){
    gauge_ptrs[dir]=(*m_cfg)[dir].getDataConstPtr();}

      // create the quda gauge params
 
 QudaGaugeParam gauge_param = newQudaGaugeParam();
 m_gauge_info->setQudaGaugeParam(gauge_param);
 gauge_param.location = QUDA_CPU_FIELD_LOCATION;

      // load onto gpu device (stored in gaugePrecise pointer)
 loadGaugeQuda((void *)gauge_ptrs.data(), &gauge_param); 
 QudaInfo::gauge_config_on_device=true;

 QudaGaugeObservableParam param = newQudaGaugeObservableParam();
 param.compute_plaquette = QUDA_BOOLEAN_TRUE;
 gaugeObservablesQuda(&param);
 printLaph("Gauge configuration copied to the device");
 printLaph(make_strf(" Computed mean plaquette = %15.12f",param.plaquette[0]));
 printLaph(make_strf("  spatial mean plaquette = %15.12f",param.plaquette[1]));
 printLaph(make_strf(" temporal mean plaquette = %15.12f",param.plaquette[2]));

      // remove out of the HostGlobal if requested
 if (removeFromHost){ eraseDataOnHost();}
}



void GaugeConfigurationHandler::eraseDataOnDevice()
{
 if (QudaInfo::gauge_config_on_device){
    freeGaugeQuda();}
 QudaInfo::gauge_config_on_device=false;
}



 // **********************************************************************
}

