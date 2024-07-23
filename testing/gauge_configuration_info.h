#ifndef LAPH_GAUGE_CONFIGURATION_INFO_H
#define LAPH_GAUGE_CONFIGURATION_INFO_H

#include "xml_handler.h"
#include "quda.h"
#include "layout_info.h"

namespace LaphEnv {

// ********************************************************************
// *                                                                  *
// *  Class "GaugeConfigurationInfo" holds information about a        *
// *  4-dimensional gauge configuration.  The constructor requires    *
// *  an XMLHandler.  The XML input must have the following form:     *
// *                                                                  *
// *      <GaugeConfigurationInfo>                                    *
// *         <EnsembleName>...</EnsembleName>                         *
// *         <FileFormat> .... </FileFormat>                          *
// *         <ConfigType> .... </ConfigType> (optional)               *
// *         <FileName> .... </FileName>                              *
// *         <ConfigNumber> .... </ConfigNumber>                      *
// *         <MarkovChainNumber>...</MarkovChainNumber> (optional)    *
// *         <NOMId> ...</NOMId> (optional)                           *
// *         <GluonAnisotropy> ... </GluonAnisotropy> (optional)      *
// *         <FermionTimeBC> ... </FermionTimeBC> (optional)          *
// *      </GaugeConfigurationInfo>                                   *
// *                                                                  *
// *  Information about the individual tags now follows:              *
// *                                                                  *
// *    <EnsembleName>                                                *
// *       - simple descriptive name, such as CLS_E250                *
// *       - subsequent tasks will have to use the same name          *
// *    <FileFormat>                                                  *
// *       - supported options:                                       *
// *          -  CERN  or  CLS                                        *
// *          -  SZIN_SP or USQCD_SP  (DP for double precision)       *
// *    <ConfigType>                                                  *
// *       - supported options currently:                             *
// *          -  Wilson                                               *
// *          -  WilsonImproved (default)                             *
// *    <FileName>                                                    *
// *       - name of file containing the configuration data           *
// *    <ConfigNumber>                                                *
// *       - integer number of the config in the ensemble             *
// *    <MarkovChainNumber>                                           *
// *       - integer number of the Markov chain (0 default)           *
// *    <NOMId>                                                       *
// *       - Id string to use in the NamedObjMap                      *
// *       - if absent, "default_gauge_id" is used                    *
// *    <GluonAnisotropy>                                             *
// *       - real number: spatial spacing/ temporal spacing           *
// *       - default 1.0                                              *
// *    <FermionTimeBC>                                               *
// *       - temporal boundary conditions for the fermion fields      *
// *       - supported options:                                       *
// *          - antiperiodic (default), periodic                      *
// *                                                                  *
// *                                                                  *
// ********************************************************************


class GaugeConfigurationInfo
{

  std::string ensemble_name;
  std::string file_format;
  std::string config_type;
  std::string file_name;
  int config_num;
  int mc_chain_num;
  std::string nom_id;
  double gluon_anisotropy;
  char fermion_time_bc;


 public:

  GaugeConfigurationInfo(const XMLHandler& xml_rdr);

  GaugeConfigurationInfo(const GaugeConfigurationInfo& rhs);
                 
  GaugeConfigurationInfo& operator=(const GaugeConfigurationInfo& rhs);

  ~GaugeConfigurationInfo(){}
  
  void checkEqual(const GaugeConfigurationInfo& rhs) const;

  bool operator==(const GaugeConfigurationInfo& rhs) const;

  std::string getEnsembleName() const { return ensemble_name; }

  std::string getConfigType() const { return config_type; }

  std::string getFileName() const { return file_name; }

  int getConfigNumber() const { return config_num; }

  int getMarkovChainNumber() const { return mc_chain_num; }
  
  std::string getNOMId() const { return nom_id; }

  int getTimeExtent() const { return LayoutInfo::getLattSizes()[3]; }
  
  bool isFermionTimeBCAntiPeriodic() const {return (fermion_time_bc=='A');}

  std::string output(int indent = 0) const;

  void output(XMLHandler& xmlout) const;

  void setQudaGaugeParam(QudaGaugeParam &gauge_param) const;


 private:

  void set_info(XMLHandler& xmlg);

  void check_valid(const std::string& sdata, const std::string& tag, 
                   const std::vector<std::string>& allowed);

  friend class GaugeConfigReader;

};


// **********************************************************
}
#endif
