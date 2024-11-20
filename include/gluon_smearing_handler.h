#ifndef GLUON_SMEARING_HANDLER_H
#define GLUON_SMEARING_HANDLER_H

#include "data_io_handler.h"
#include "field_smearing_info.h"
#include "gauge_configuration_handler.h"

namespace LaphEnv {

// *****************************************************************
// *                                                               *
// *  "GluonSmearingHandler" handles access to the smeared gauge   *
// *  field.  See declarations below for available member          *
// *  functions.                                                   *
// *                                                               *
// *  GluonSmearingHandler:                                        *
// *                                                               *
// *     -- computes smeared gauge field and either writes to      *
// *        a file or inserts into the NamedObjMap                 *
// *                                                               *
// *                                                               *
// *****************************************************************

class GluonSmearingHandler {

  typedef UIntKey RecordKey; // dummy key for IOMap

  // pointers to internal infos (managed by this handler
  // with new and delete)

  const GluonSmearingInfo *smearPtr;
  const GaugeConfigurationInfo *uPtr;
  std::string h_filename;
  bool m_read_mode;

  // storage and/or references to internal data

  DataGetHandlerSFO<GluonSmearingHandler, RecordKey, std::vector<LattField>>
      *dh_ptr;

  // prevent copying ... handler might contain large
  // amounts of data

  GluonSmearingHandler(const GluonSmearingHandler &) = delete;
  GluonSmearingHandler &operator=(const GluonSmearingHandler &) = delete;

public:
  GluonSmearingHandler();

  GluonSmearingHandler(const GluonSmearingInfo &gluon_smearing,
                       const GaugeConfigurationInfo &gauge,
                       const std::string &smearedFieldFileName,
                       bool readmode = true);

  void setInfo(const GluonSmearingInfo &gluon_smearing,
               const GaugeConfigurationInfo &gauge,
               const std::string &smearedFieldFileName, bool readmode = true);

  ~GluonSmearingHandler();

  void clear(); // clears everything in the handler

  // access to the info

  bool isInfoSet() const;

  const GluonSmearingInfo &getGluonSmearingInfo() const;

  const GaugeConfigurationInfo &getGaugeConfigurationInfo() const;

  std::string getFileName() const { return h_filename; }

  // Computes the smeared gauge field and puts it into the file
  // whose name is stored in "h_filename".  If "h_filename"
  // starts with "NOM_", results are actually put into the
  // "NamedObjMap" with identifier "Laph--SmearedGaugeField".

  void computeSmearedGaugeField();

  const std::vector<LattField> &getSmearedGaugeField() const;

  // remove one time slice, or clear all of the smeared gauge
  // fields time slices from memory

  void clearData();

  void copyDataToDevice();

  void eraseDataOnDevice();

  bool isDataOnDevice() const { return QudaInfo::smeared_gauge_on_device; }

private:
  void set_info(const GluonSmearingInfo &gluon_smearing,
                const GaugeConfigurationInfo &gauge,
                const std::string &smearedFieldFileName, bool readmode);

  void filefail(const std::string &message);

  void check_info_set(const std::string &name) const;

  bool checkHeader(XMLHandler &xmlr);

  void writeHeader(XMLHandler &xmlw);

  friend class DataGetHandlerSF<GluonSmearingHandler, UIntKey,
                                std::vector<LattField>>;
  friend class DataPutHandlerSF<GluonSmearingHandler, UIntKey,
                                std::vector<LattField>>;
  friend class DataGetHandlerSFNOM<GluonSmearingHandler, UIntKey,
                                   std::vector<LattField>>;
  friend class DataPutHandlerSFNOM<GluonSmearingHandler, UIntKey,
                                   std::vector<LattField>>;
};
} // namespace LaphEnv
#endif
