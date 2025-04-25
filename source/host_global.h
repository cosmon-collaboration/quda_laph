#ifndef HOST_GLOBAL_H
#define HOST_GLOBAL_H

#include <memory>
#include "latt_field.h"
#include "gauge_configuration_info.h"
#include "field_smearing_info.h"


namespace LaphEnv {

// *********************************************************************************
// *                                                                               *
// *   Class "HostGlobal" is a singleton which stores important data which needs   *
// *   to be persistent across several tasks. In particular, the gauge             *
// *   configuration, the smeared gauge config, the LapH eigenvectors are stored.  *
// *                                                                               *
// *   Global data (persistent across tasks) (declared in quda_laph.cc)            *
// *                                                                               *
// *     std::vector<LattField> HostGlobal::theGaugeConfig;                        *
// *     std::vector<LattField> HostGlobal::theSmearedGaugeConfig;                 *
// *     std::vector<LattField> HostGlobal::theLaphEigvecs;                        *
// *     std::unique_ptr<GaugeConfigurationInfo> theGaugeConfigInfoPtr;            *
// *     std::unique_ptr<GluonSmearingInfo> theGluonSmearingInfoPtr;               *
// *     std::unique_ptr<QuarkSmearingInfo> theQuarkSmearingInfoPtr;               *
// *                                                                               *
// *   The "set" routines return a reference to the data for the user to then      *
// *   use for assignment purposes; if the info matches, nothing is done to any    *
// *   data already in the HostGlobal; if the info does not match, the data is     *
// *   first cleared. Setting the smeared gauge configuration requires             *
// *   agreement with the info of the current gauge configuration; setting         *
// *   the Laph eigenvectors requires info agreement with the gauge config and     *
// *   smeared gauge config.  If you set the gauge config and the info does not    *
// *   match the current info, the smeared gauge and the Laph eigenvectors are     *
// *   cleared.  If you set the smeared config, then the Laph eigenvectors are     *
// *   also cleared.                                                               *
// *                                                                               *
// *     vector<LattField>& ref(GB::setGaugeConfig(gaugeinfo));                    *
// *     vector<LattField>& ref(GB::setSmearedGaugeConfig(gluonsmearinfo,          *
// *                                                      gaugeinfo));             *
// *     vector<LattField>& ref(GB::setLaphEigvecs(quarksmearinfo,gluonsmearinfo,  *
// *                                               gaugeinfo));                    *
// *                                                                               *
// *   The "get" routines return a const reference to the data and the             *
// *   info must match or an exception is thrown.                                  *
// *                                                                               *
// *     const vector<LattField>& ref(GB::getGaugeConfig(gaugeinfo));              *
// *     const vector<LattField>& ref(GB::getSmearedGaugeConfig(gluonsmearinfo));  *
// *     const vector<LattField>& ref(GB::getLaphEigvecs(quarksmearinfo));         *
// *                                                                               *
// *   Data and info can be erased from HostGlobal using the "clear" members:      *
// *                                                                               *
// *     GB::clear();                                                              *
// *     GB::clearGaugeConfig();                                                   *
// *     GB::clearSmearedGaugeConfig();                                            *
// *     GB::clearLaphEigvecs();                                                   *
// *                                                                               *
// *   Data existence can be queried with the following members:                   *
// *                                                                               *
// *     bool flag=GB::theGaugeConfigIsSet(gaugeinfo);                             *
// *     bool flag=GB::theSmearedGaugeConfigIsSet(gluonsmearinfo,gaugeinfo);       *
// *     bool flag=GB::theLaphEigvecsAreSet(quarksmearinfo,gluonsmearinfo,         *
// *                                        gaugeinfo);                            *
// *                                                                               *
// *   Data but NOT info can be erased from HostGlobal using the "clearData"       *
// *   members:                                                                    *
// *                                                                               *
// *     GB::clearData();                                                          *
// *     GB::clearGaugeConfigData();                                               *
// *     GB::clearSmearedGaugeConfigData();                                        *
// *     GB::clearLaphEigvecsData();                                               *
// *                                                                               *
// *********************************************************************************


class HostGlobal
{

   HostGlobal() {}

   ~HostGlobal() {}

   static std::vector<LattField> theGaugeConfig;
   
   static std::unique_ptr<GaugeConfigurationInfo> theGaugeConfigInfoPtr;

   static std::vector<LattField> theSmearedGaugeConfig;
   
   static std::unique_ptr<GluonSmearingInfo> theGluonSmearingInfoPtr;

   static std::vector<LattField> theLaphEigvecs;
   
   static std::unique_ptr<QuarkSmearingInfo> theQuarkSmearingInfoPtr;


   HostGlobal(const HostGlobal &in) = delete;           // no copy constructor

   HostGlobal(HostGlobal &in) = delete;                 // no copy constructor

   HostGlobal& operator=(const HostGlobal &in) = delete; // not assignable

   HostGlobal& operator=(HostGlobal &in) = delete;       // not assignable


   static bool gauge_config_info_matches(const GaugeConfigurationInfo& gaugeinfo)
   {
    return theGaugeConfigInfoPtr && (*theGaugeConfigInfoPtr==gaugeinfo);
   }

   static bool gluon_smearing_info_matches(const GluonSmearingInfo& gluonsmearinfo)
   {
    return theGluonSmearingInfoPtr && (*theGluonSmearingInfoPtr==gluonsmearinfo);
   }

   static bool quark_smearing_info_matches(const QuarkSmearingInfo& quarksmearinfo)
   {
    return theQuarkSmearingInfoPtr && (*theQuarkSmearingInfoPtr==quarksmearinfo);
   }


 public:

   static std::vector<LattField>& setGaugeConfig(const GaugeConfigurationInfo& gaugeinfo)
   {
    if (!gauge_config_info_matches(gaugeinfo)){
       clear();   // clear gauge, smeared gauge, Laph eigenvecs if putting in new configuration
       theGaugeConfigInfoPtr=std::make_unique<GaugeConfigurationInfo>(gaugeinfo);}
    return theGaugeConfig;
   }
   
   static const std::vector<LattField>& getGaugeConfig(const GaugeConfigurationInfo& gaugeinfo)
   {
    if (!gauge_config_info_matches(gaugeinfo)){
       throw(std::invalid_argument("Incorrect GaugeConfigurationInfo in getGaugeConfig from HostGlobal"));}
    return theGaugeConfig;
   }

   static std::vector<LattField>& setSmearedGaugeConfig(const GluonSmearingInfo& gluonsmearinfo, 
                                                        const GaugeConfigurationInfo& gaugeinfo)
   {
    if (!theGaugeConfigInfoPtr){
       theGaugeConfigInfoPtr=std::make_unique<GaugeConfigurationInfo>(gaugeinfo);}
    else if (!gauge_config_info_matches(gaugeinfo)){
       throw(std::invalid_argument("Incorrect GaugeConfigurationInfo in setSmearedGaugeConfig from HostGlobal"));}
    if (!gluon_smearing_info_matches(gluonsmearinfo)){
       clearSmearedGaugeConfig(); clearLaphEigvecs();
       theGluonSmearingInfoPtr=std::make_unique<GluonSmearingInfo>(gluonsmearinfo);}
    return theSmearedGaugeConfig;
   }
   
   static const std::vector<LattField>& getSmearedGaugeConfig(const GluonSmearingInfo& gluonsmearinfo,
                                                              const GaugeConfigurationInfo& gaugeinfo)
   {
    if (!gauge_config_info_matches(gaugeinfo)){
       throw(std::invalid_argument("Incorrect GaugeConfigurationInfo in getSmearedGaugeConfig from HostGlobal"));}
    if (!gluon_smearing_info_matches(gluonsmearinfo)){
       throw(std::invalid_argument("Incorrect GluonSmearingInfo in getSmearedGaugeConfig from HostGlobal"));}
    return theSmearedGaugeConfig;
   }

   static std::vector<LattField>& setLaphEigvecs(const QuarkSmearingInfo& quarksmearinfo,
                                                 const GluonSmearingInfo& gluonsmearinfo, 
                                                 const GaugeConfigurationInfo& gaugeinfo)
   {
    if (!theGaugeConfigInfoPtr){
       theGaugeConfigInfoPtr=std::make_unique<GaugeConfigurationInfo>(gaugeinfo);}
    else if (!gauge_config_info_matches(gaugeinfo)){
       throw(std::invalid_argument("Incorrect GaugeConfigurationInfo in setLaphEigvecs from HostGlobal"));}
    if (!theGluonSmearingInfoPtr){
       theGluonSmearingInfoPtr=std::make_unique<GluonSmearingInfo>(gluonsmearinfo);}
    else if (!gluon_smearing_info_matches(gluonsmearinfo)){
       throw(std::invalid_argument("Incorrect GluonSmearingInfo in setLaphEigvecs from HostGlobal"));}
    if (!quark_smearing_info_matches(quarksmearinfo)){
       clearLaphEigvecs();
       theQuarkSmearingInfoPtr=std::make_unique<QuarkSmearingInfo>(quarksmearinfo);}
    return theLaphEigvecs;
   }
   
   static const std::vector<LattField>& getLaphEigvecs(const QuarkSmearingInfo& quarksmearinfo,
                                                       const GluonSmearingInfo& gluonsmearinfo, 
                                                       const GaugeConfigurationInfo& gaugeinfo)
   {
    if (!gauge_config_info_matches(gaugeinfo)){
       throw(std::invalid_argument("Incorrect GaugeConfigurationInfo in getLaphEigvecs from HostGlobal"));}
    if (!gluon_smearing_info_matches(gluonsmearinfo)){
       throw(std::invalid_argument("Incorrect GluonSmearingInfo in getLaphEigvecs from HostGlobal"));}
    if (!quark_smearing_info_matches(quarksmearinfo)){
       throw(std::invalid_argument("Incorrect QuarkSmearingInfo in getLaphEigvecs from HostGlobal"));}
    return theLaphEigvecs;
   }


   static void clear()
   {
    theGaugeConfig.clear(); 
    theGaugeConfigInfoPtr.reset();
    theSmearedGaugeConfig.clear(); 
    theGluonSmearingInfoPtr.reset();
    theLaphEigvecs.clear(); 
    theQuarkSmearingInfoPtr.reset();
   }

   static void clearGaugeConfig()
   {
    theGaugeConfig.clear();
    theGaugeConfigInfoPtr.reset();
   }

   static void clearSmearedGaugeConfig()
   {
    theSmearedGaugeConfig.clear();
    theGluonSmearingInfoPtr.reset();
   }

   static void clearLaphEigvecs()
   {
    theLaphEigvecs.clear(); 
    theQuarkSmearingInfoPtr.reset();
   }


   static void clearData()   // does not clear infos
   {
    theGaugeConfig.clear(); 
    theSmearedGaugeConfig.clear(); 
    theLaphEigvecs.clear(); 
   }

   static void clearGaugeConfigData()
   {
    theGaugeConfig.clear();
   }

   static void clearSmearedGaugeConfigData()
   {
    theSmearedGaugeConfig.clear();
   }

   static void clearLaphEigvecsData()
   {
    theLaphEigvecs.clear(); 
   }


   static bool theGaugeConfigIsSet(const GaugeConfigurationInfo& gaugeinfo)
   { return gauge_config_info_matches(gaugeinfo) && (!theGaugeConfig.empty());}

   static bool theSmearedGaugeConfigIsSet(const GluonSmearingInfo& gluonsmearinfo,
                                          const GaugeConfigurationInfo& gaugeinfo)
   { return gauge_config_info_matches(gaugeinfo) && gluon_smearing_info_matches(gluonsmearinfo) 
            && (!theSmearedGaugeConfig.empty());}

   static bool theLaphEigvecsAreSet(const QuarkSmearingInfo& quarksmearinfo,
                                    const GluonSmearingInfo& gluonsmearinfo,
                                    const GaugeConfigurationInfo& gaugeinfo)
   { return gauge_config_info_matches(gaugeinfo) && gauge_config_info_matches(gaugeinfo) 
            && quark_smearing_info_matches(quarksmearinfo) && (!theLaphEigvecs.empty());}


};

typedef HostGlobal  GB;


// *************************************************************
}
#endif
