#ifndef GLUON_SMEARING_HANDLER_H
#define GLUON_SMEARING_HANDLER_H

#include "gauge_configuration_handler.h"
#include "field_smearing_info.h"
#include "data_io_handler.h"

namespace LaphEnv {


// *****************************************************************
// *                                                               *
// *  "GluonSmearingHandler" handles access to the smeared gauge   *
// *  field.  See declarations below for available member          *
// *  functions.                                                   *
// *                                                               *
// *  GluonSmearingHandler:                                        *
// *                                                               *
// *     -- computes smeared gauge field and places in the         *
// *        HostGlobal, and optionally writes to a file            *
// *                                                               *
// *                                                               *
// *****************************************************************



class GluonSmearingHandler
{

   typedef UIntKey  RecordKey;   // dummy key for IOMap
   
       // pointers to internal infos (managed by this handler
       // with new and delete)

   const GluonSmearingInfo *smearPtr;
   const GaugeConfigurationInfo *uPtr;
   std::string h_filename;
   bool m_read_mode;

       // prevent copying ... handler might contain large
       // amounts of data

   GluonSmearingHandler(const GluonSmearingHandler&) = delete;
   GluonSmearingHandler& operator=(const GluonSmearingHandler&) = delete;


 public:


   GluonSmearingHandler();

   GluonSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                        const GaugeConfigurationInfo& gauge,
                        const std::string& smearedFieldFileName,
                        bool readmode=true); 

   void setInfo(const GluonSmearingInfo& gluon_smearing,
                const GaugeConfigurationInfo& gauge,
                const std::string& smearedFieldFileName,
                bool readmode=true);

   ~GluonSmearingHandler();

   void clear();  // clears everything in the handler


           // access to the info

   bool isInfoSet() const;

   const GluonSmearingInfo& getGluonSmearingInfo() const;

   const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

   std::string getFileName() const {return h_filename;}

           // Computes the smeared gauge field into HostGlobal and 
           // optionally puts it into the file whose name is stored in 
           // "h_filename".

   void computeSmearedGaugeField();

           // Returns a reference to the smeared gauge field.
           // Checks first if in HostGlobal; if not, reads from
           // file and stores it in HostGlobal

   const std::vector<LattField>& getSmearedGaugeField();

   bool querySmearedGaugeField();
   
           // clear the smeared config from HostGlobal
   void clearData();

   void copyDataToDevice();
    
   void eraseDataOnDevice();
    
   bool isDataOnDevice() const {return QudaInfo::smeared_gauge_on_device;}


 private:

   void set_info(const GluonSmearingInfo& gluon_smearing,
                 const GaugeConfigurationInfo& gauge,
                 const std::string& smearedFieldFileName,
                 bool readmode);
   
   bool loadSmearedGaugeField();  // load from file, put into HostGlobal

   void filefail(const std::string& message);

   void check_info_set(const std::string& name) const;
 
   bool checkHeader(XMLHandler& xmlr);

   void writeHeader(XMLHandler& xmlw);

   friend class DataGetHandlerSF<GluonSmearingHandler,UIntKey,
                                 std::vector<LattField>>;
   friend class DataPutHandlerSF<GluonSmearingHandler,UIntKey,
                                 std::vector<LattField>>;

};



// *******************************************************************
}
#endif  
