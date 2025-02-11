#ifndef LAPH_GAUGE_CONFIG_HANDLER_H
#define LAPH_GAUGE_CONFIG_HANDLER_H

#include "gauge_configuration_info.h"
#include "latt_field.h"
#include "quda_info.h"

namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *  Class "GaugeConfigurationHandler" manages information about    *
// *  and access to the gauge configuration in the Laph environment. *
// *  It contains a "GaugeConfigurationInfo" that holds the gauge    *
// *  header, as well as a reference to the actual cfg.              *
// *                                                                 *
// *    Basic usage:                                                 *
// *                                                                 *
// *      (a) declare a GaugeConfigurationHandler                    *
// *      (b) set the info  -- via constructor or setInfo(..)        *
// *      (c) set the data  -- via setData(..)                       *
// *      (d) use getData(..) for access to gauge configuration      *
// *                                                                 *
// *     -- getInfo(..) provides access to configuration info        *
// *     -- can be cleared and info/data reset                       *
// *                                                                 *
// *    Example:                                                     *
// *       GaugeConfigurationInfo Gin;                               *
// *       GaugeConfigurationHandler uHandler;                       *
// *       uHandler.setInfo(Gin);                                    *
// *       uHandler.setData();                                       *
// *       const multi1d<LatticeColorMatrix>& U=uHandler.getData();  *
// *    --> access to config through U                               *
// *                                                                 *
// *******************************************************************


class GaugeConfigurationHandler
{

      // pointer to the info about the gauge config (internally
      // managed by this handler)
      
    const GaugeConfigurationInfo* m_gauge_info;

      // pointer to the gauge field (external: in HostGlobal)
 
    const std::vector<LattField>* m_cfg;

      // prevent copying
    GaugeConfigurationHandler(const GaugeConfigurationHandler& gch);
    GaugeConfigurationHandler& operator=(const GaugeConfigurationHandler& gch);


  public:

    GaugeConfigurationHandler();
    
    GaugeConfigurationHandler(const GaugeConfigurationInfo& gaugeinfo);

    void setInfo(const GaugeConfigurationInfo& gaugeinfo);

    ~GaugeConfigurationHandler();   // destructor, but leaves config in HostGlobal
                                    // if one device, leaves it on the device too
    void clear();   // clears this handler, but leaves config in HostGlobal
                    // if one device, leaves it on the device too
    void eraseDataOnHost();   // removes config from HostGlobal on host
  
    void setData();
    
    
      // access to the gauge configuration and its info

    const std::vector<LattField>& getData();

    const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

    bool isInfoSet() const { return (m_gauge_info!=0);}

    bool isDataSet() const { return (m_cfg!=0);}

    void getXMLInfo(XMLHandler& gauge_xmlinfo) const;

    void applyFermionTemporalAntiPeriodic();  // multiplies host temporal links on time Nt-1
                                              // by -1 as expected by Quda to handle fermion 
                                              // antiperiodic temporal boundary condition

       //  putting to and erasing from the GPU device

    void copyDataToDevice(bool removeFromHost=false);
    
    void eraseDataOnDevice();
    
    bool isDataOnDevice() const {return QudaInfo::gauge_config_on_device;}


  private:

    void set_info(const GaugeConfigurationInfo& gaugeinfo);
  
    void check_info_set(const std::string& name) const;

    void initialize_config();
 
};

// *************************************************************************
}
#endif
