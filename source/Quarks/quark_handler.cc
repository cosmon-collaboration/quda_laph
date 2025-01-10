#include "QudaLaphIncludes.h"

#include "quark_handler.h"
#include "array.h"
#include "field_ops.h"
#include "multi_compare.h"

// STILL TO DO:  single asynchronous thread for output so can continue next
// computations
//                   while output occurs

typedef std::complex<double> dcmplx;
typedef std::complex<float> fcmplx;

namespace LaphEnv {

void QuarkHandler::RecordKey::output(XMLHandler &xmlw) const {
  xmlw.set_root("RecordKey");
  xmlw.put_child("Spin", make_string(getSpin()));
  xmlw.put_child("Time", make_string(getTime()));
  xmlw.put_child("SpinLaphEigvecIndex", make_string(getSpinLaphEigvecIndex()));
}

QuarkHandler::FileKey::FileKey(const LaphNoiseInfo &in_noise, const int tprojind)
    : noise(in_noise), time_proj_index(tprojind) {}

QuarkHandler::FileKey::FileKey(XMLHandler &xmlr) {
  try {
    XMLHandler xmlf(xmlr, "FileKey");
    noise = LaphNoiseInfo(xmlf);
    xmlread(xmlf, "TimeProjectorIndex", time_proj_index,
            "QuarkHandler::FileKey");
  } catch (const std::exception &xp) {
    errorLaph("Could not read QuarkHandler::FileKey");
  }
}

QuarkHandler::FileKey::FileKey(const FileKey &rhs)
    : noise(rhs.noise), time_proj_index(rhs.time_proj_index) {}

QuarkHandler::FileKey &QuarkHandler::FileKey::operator=(const FileKey &rhs) {
  noise = rhs.noise;
  time_proj_index = rhs.time_proj_index;
  return *this;
}

void QuarkHandler::FileKey::output(XMLHandler &xmlw) const {
  xmlw.set_root("FileKey");
  XMLHandler xmln;
  noise.output(xmln);
  xmlw.put_child(xmln);
  xmlw.put_child("TimeProjectorIndex", make_string(time_proj_index));
}

bool QuarkHandler::FileKey::operator<(const FileKey &rhs) const {
  return multiLessThan(time_proj_index, rhs.time_proj_index, noise, rhs.noise);
}

bool QuarkHandler::FileKey::operator==(const FileKey &rhs) const {
  return multiEqual(time_proj_index, rhs.time_proj_index, noise, rhs.noise);
}

bool QuarkHandler::FileKey::operator!=(const FileKey &rhs) const {
  return multiNotEqual(time_proj_index, rhs.time_proj_index, noise, rhs.noise);
}

QuarkHandler::StorageKey::StorageKey(const FileKey &in_fkey,
                                     const RecordKey &in_stkey,
                                     const DirPath &in_disp, int in_disp_length)
    : fkey(in_fkey), rkey(in_stkey), disp(in_disp),
      disp_length(in_disp_length) {
  if (disp.Length() == 0)
    disp_length = 0;
  if ((disp.Length() > 0) && (disp_length <= 0)) {
    errorLaph("Invalid displacement length");
  }
}

QuarkHandler::StorageKey::StorageKey(const FileKey &in_fkey,
                                     const RecordKey &in_stkey)
    : fkey(in_fkey), rkey(in_stkey), disp_length(0) {}

QuarkHandler::StorageKey::StorageKey(const QuarkHandler::StorageKey &in)
    : fkey(in.fkey), rkey(in.rkey), disp(in.disp), disp_length(in.disp_length) {
}

QuarkHandler::StorageKey &
QuarkHandler::StorageKey::operator=(const QuarkHandler::StorageKey &in) {
  fkey = in.fkey;
  rkey = in.rkey;
  disp = in.disp;
  disp_length = in.disp_length;
  return *this;
}

bool QuarkHandler::StorageKey::operator<(
    const QuarkHandler::StorageKey &rhs) const {
  return multiLessThan(rkey, rhs.rkey, disp, rhs.disp, disp_length,
                       rhs.disp_length, fkey, rhs.fkey);
}

bool QuarkHandler::StorageKey::operator==(
    const QuarkHandler::StorageKey &rhs) const {
  return multiEqual(rkey, rhs.rkey, disp, rhs.disp, disp_length,
                    rhs.disp_length, fkey, rhs.fkey);
}

bool QuarkHandler::StorageKey::operator!=(
    const QuarkHandler::StorageKey &rhs) const {
  return multiNotEqual(rkey, rhs.rkey, disp, rhs.disp, disp_length,
                       rhs.disp_length, fkey, rhs.fkey);
}

QuarkHandler::QuarkHandler()
    : uPtr(0), gSmearPtr(0), qSmearPtr(0), dilPtr(0), qactionPtr(0), fPtr(0),
      invertPtr(0), compute_mode(false), preconditioner(0), dilHandler(0),
      DHputPtr(0), DHgetPtr(0), normal_mode(true) {}

QuarkHandler::QuarkHandler(
    const GaugeConfigurationInfo &gaugeinfo,
    const GluonSmearingInfo &gluonsmear,
    const QuarkSmearingInfo &quarksmear,
    const DilutionSchemeInfo &dil,
    const QuarkActionInfo &quark,
    const FileListInfo &flist,
    const std::string &smeared_quark_filestub,
    const std::string &smeared_gauge_filename,
    const bool setComputeMode)
    : invertPtr(0), compute_mode(setComputeMode), preconditioner(0),
      dilHandler(0), DHputPtr(0), DHgetPtr(0), normal_mode(true) {
  set_info(gaugeinfo, gluonsmear, quarksmear, dil, quark, flist,
           smeared_quark_filestub, smeared_gauge_filename);
}

void QuarkHandler::setInfo(
    const GaugeConfigurationInfo &gaugeinfo,
    const GluonSmearingInfo &gluonsmear,
    const QuarkSmearingInfo &quarksmear,
    const DilutionSchemeInfo &dil,
    const QuarkActionInfo &quark,
    const FileListInfo &flist,
    const std::string &smeared_quark_filestub,
    const std::string &smeared_gauge_filename,
    const bool setComputeMode) {
  clear();
  compute_mode = setComputeMode;
  set_info(gaugeinfo, gluonsmear, quarksmear, dil, quark, flist,
           smeared_quark_filestub, smeared_gauge_filename);
}

void QuarkHandler::set_info(const GaugeConfigurationInfo &gaugeinfo,
                            const GluonSmearingInfo &gluonsmear,
                            const QuarkSmearingInfo &quarksmear,
                            const DilutionSchemeInfo &dil,
                            const QuarkActionInfo &quark,
                            const FileListInfo &flist,
                            const std::string &smeared_quark_filestub,
                            const std::string &smeared_gauge_filename) {
  try {
    uPtr = new GaugeConfigurationInfo(gaugeinfo);
    gSmearPtr = new GluonSmearingInfo(gluonsmear);
    qSmearPtr = new QuarkSmearingInfo(quarksmear);
    dilPtr = new DilutionSchemeInfo(dil);
    qactionPtr = new QuarkActionInfo(quark);
    fPtr = new FileListInfo(flist);
    if (compute_mode) {
      DHputPtr =
          new DataPutHandlerMF<QuarkHandler, FileKey, RecordKey, DataType>(
              *this, *fPtr, "Laph--QuarkSink", "QuarkHandlerDataFile");
    } else {
      DHgetPtr =
          new DataGetHandlerMF<QuarkHandler, FileKey, RecordKey, DataType>(
              *this, *fPtr, "Laph--QuarkSink", "QuarkHandlerDataFile");
    }
  } catch (const std::exception &xp) {
    errorLaph(make_strf("allocation problem in QuarkHandler: %s", xp.what()));
  }

  connectGaugeConfigurationHandler();
  connectGluonSmearingHandler(smeared_gauge_filename);
  connectQuarkSmearingHandler(smeared_quark_filestub);
  connectDilutionHandler();
}

void QuarkHandler::setInverter(const InverterInfo &invinfo) {
  if (!compute_mode) {
    errorLaph("Cannot setInverter in QuarkHandler unless in compute mode");
  }
  if (!isInfoSet()) {
    errorLaph("Cannot setInverter in QuarkHandler unless info is set");
  }
  try {
    delete invertPtr;
    invertPtr = new InverterInfo(invinfo);
  } catch (const std::exception &xp) {
    errorLaph("allocation error in QuarkHandler::setInverter");
  }

  // set up the inverter params for quda
  invertPtr->setQudaInvertParam(quda_inv_param, *qactionPtr);

  // check that Dirac-Pauli basis is used
  if (quda_inv_param.gamma_basis != QUDA_DIRAC_PAULI_GAMMA_BASIS) {
    errorQuda("The Dirac-Pauli basis must be used");
  }
}

const InverterInfo &QuarkHandler::getInverterInfo() const {
  if (invertPtr != 0) {
    errorLaph("error in QuarkHandler: must setInverter before calling "
              "getInverterInfo");
  }
  return *invertPtr;
}

void QuarkHandler::clearSinkComputations() { sinkComps.computations.clear(); }

void QuarkHandler::setSinkComputations(const XMLHandler &xmlin) {
  if (!compute_mode) {
    errorLaph(
        "Cannot setSinkComputations in QuarkHandler unless in compute mode");
  }
  if (!isInfoSet()) {
    errorLaph("Cannot setSinkComputations in QuarkHandler unless info is set");
  }
  XMLHandler xmlrdr(xmlin);
  if (!sinkComps.computations.empty())
    sinkComps.computations.clear();

  if (xml_tag_count(xmlrdr, "SinkComputations") == 1) {

    XMLHandler xmlrd(xmlrdr, "SinkComputations");

    uint nSinkLaphBatch, nSinkQudaBatch, nEigQudaBatch;
    xmlread(xmlrd, "NumSinksBeforeProject", nSinkLaphBatch,
            "LAPH_QUARK_LINE_ENDS");
    xmlread(xmlrd, "NumSinksInProjectBatch", nSinkQudaBatch,
            "LAPH_QUARK_LINE_ENDS");
    xmlread(xmlrd, "NumEigsInProjectBatch", nEigQudaBatch,
            "LAPH_QUARK_LINE_ENDS");

    int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
    int ndil = dilHandler->getNumberOfSpinEigvecProjectors();
    if ((int(nSinkLaphBatch) > ndil) || (nSinkLaphBatch == 0)) {
      errorLaph(make_strf("Invalid value %d for <NumSinksBeforeProject>: must "
                          "be between 1 and ndil=%d",
                          nSinkLaphBatch, ndil));
    }
    if ((int(nEigQudaBatch) > nEigs) || (nEigQudaBatch == 0)) {
      errorLaph(make_strf("Invalid value %d for <NumEigsInProjectBatch>: must "
                          "be between 1 and nEigs=%d",
                          nEigQudaBatch, nEigs));
    }
    if ((nSinkQudaBatch > nSinkLaphBatch) || (nSinkQudaBatch == 0)) {
      errorLaph(make_strf("Invalid value %d for <NumSinksInProjectBatch>: must "
                          "be between 1 and %d",
                          nSinkQudaBatch, nSinkLaphBatch));
    }
    sinkComps.nSinkLaphBatch = nSinkLaphBatch;
    sinkComps.nSinkQudaBatch = nSinkQudaBatch;
    sinkComps.nEigQudaBatch = nEigQudaBatch;

    if (xml_tag_count(xmlrd, "NoiseList_TimeProjIndexList") == 1) {
      XMLHandler xmlr(xmlrd, "NoiseList_TimeProjIndexList");
      std::vector<int> time_proj_inds;
      if (xml_tag_count(xmlr, "TimeProjIndexList") == 1) {
        XMLHandler xmltpi(xmlr, "TimeProjIndexList");
        if (xml_tag_count(xmltpi, "All") == 1) {
          int tpnum = dilHandler->getNumberOfTimeProjectors();
          time_proj_inds.resize(tpnum);
          for (int t = 0; t < tpnum; t++)
            time_proj_inds[t] = t;
        } else {
          xmlread(xmltpi, "Values", time_proj_inds, "LAPH_QUARK_LINE_ENDS");
        }
      }
      XMLHandler xmln(xmlr, "LaphNoiseList");
      std::list<XMLHandler> xmlns(xmln.find("LaphNoiseInfo"));
      for (std::list<XMLHandler>::iterator it = xmlns.begin(); it != xmlns.end();
           ++it) {
        LaphNoiseInfo aNoise(*it);
        for (int t = 0; t < int(time_proj_inds.size()); t++) {
          sinkComps.computations.push_back(
              SinkComputation(aNoise, time_proj_inds[t]));
        }
      }
    }

    if (xml_tag_count(xmlrd, "ComputationList") == 1) {
      XMLHandler xmlr(xmlrd, "ComputationList");
      std::list<XMLHandler> xmlcs(xmlr.find("Computation"));
      for (std::list<XMLHandler>::iterator it = xmlcs.begin(); it != xmlcs.end();
           ++it) {
        LaphNoiseInfo aNoise(*it);
        int time_proj_index;
        xmlread(*it, "TimeProjIndex", time_proj_index, "LAPH_QUARK_LINE_ENDS");
        sinkComps.computations.push_back(
            SinkComputation(aNoise, time_proj_index));
      }
    }
  }

  printLaph("\nLAPH_QUARK_LINE_ENDS sink computations:\n");
  printLaph(
      make_strf("  NumSinksBeforeProject = %d", sinkComps.nSinkLaphBatch));
  printLaph(
      make_strf(" NumSinksInProjectBatch = %d", sinkComps.nSinkQudaBatch));
  printLaph(
      make_strf("  NumEigsInProjectBatch = %d\n", sinkComps.nEigQudaBatch));
  printLaph(make_strf(" Number of sink computations = %d",
                      sinkComps.computations.size()));
  int count = 0;
  for (std::list<SinkComputation>::const_iterator it =
	 sinkComps.computations.begin();
       it != sinkComps.computations.end(); count++, it++) {
    XMLHandler xmlout("SinkComputation");
    XMLHandler xmln;
    it->Noise.output(xmln);
    xmlout.put_child(xmln);
    xmlout.put_child("TimeProjIndex", make_string(it->TimeProjIndex));
    printLaph(make_strf("\nSinkComputation %d:", count));
    printLaph(make_strf("%s", xmlout.output()));
  }
  printLaph("\n");
}

QuarkHandler::~QuarkHandler() { clear(); }

void QuarkHandler::clear() {
  try {
    delete uPtr;
    delete gSmearPtr;
    delete qSmearPtr;
    delete dilPtr;
    delete qactionPtr;
    delete fPtr;
    delete invertPtr;
    if (preconditioner) {
      destroyMultigridQuda(preconditioner);
    }
    clearSinkComputations();
  } catch (const std::exception &xp) {
    errorLaph("abort");
  }
  uPtr = 0;
  gSmearPtr = 0;
  qSmearPtr = 0;
  dilPtr = 0;
  qactionPtr = 0;
  fPtr = 0;
  invertPtr = 0;
  normal_mode = true;
  preconditioner = 0;

  disconnectGaugeConfigurationHandler();
  disconnectGluonSmearingHandler();
  disconnectQuarkSmearingHandler();
  disconnectDilutionHandler();
  clearData();
  delete DHgetPtr;
  DHgetPtr = 0;
  delete DHputPtr;
  DHputPtr = 0;
}

void QuarkHandler::connectGaugeConfigurationHandler() {
  if ((gaugeCounter == 0) && (gaugeHandler.get() == 0)) {
    try {
      gaugeHandler.reset(new GaugeConfigurationHandler(*uPtr));
      gaugeCounter = 1;
    } catch (const std::exception &xp) {
      errorLaph("allocation problem in "
                "QuarkHandler::connectGaugeConfigurationHandler");
    }
  } else {
    try {
      if (gaugeHandler.get() == 0)
        throw(std::invalid_argument("error"));
      uPtr->checkEqual(gaugeHandler->getGaugeConfigurationInfo());
      gaugeCounter++;
    } catch (const std::exception &xp) {
      errorLaph("inconsistent QuarkHandler::connectGaugeConfigurationHandler");
    }
  }
}

void QuarkHandler::disconnectGaugeConfigurationHandler() {
  gaugeCounter--;
  if (gaugeCounter == 0) {
    try {
      gaugeHandler.reset();
    } catch (const std::exception &xp) {
      errorLaph("delete problem in "
                "QuarkHandler::disconnectGaugeConfigurationHandler");
    }
  }
}

void QuarkHandler::connectGluonSmearingHandler(
    const std::string &smeared_gauge_filename) {
  if ((gSmearCounter == 0) && (gSmearHandler.get() == 0)) {
    try {
      gSmearHandler.reset(
          new GluonSmearingHandler(*gSmearPtr, *uPtr, smeared_gauge_filename));
      gSmearCounter = 1;
    } catch (const std::exception &xp) {
      clear();
      errorLaph(
          "allocation problem in QuarkHandler::connectGluonSmearingHandler");
    }
  } else {
    try {
      if (gSmearHandler.get() == 0)
        throw(std::invalid_argument("error"));
      uPtr->checkEqual(gSmearHandler->getGaugeConfigurationInfo());
      gSmearPtr->checkEqual(gSmearHandler->getGluonSmearingInfo());
      gSmearCounter++;
    } catch (const std::exception &xp) {
      errorLaph("inconsistent QuarkHandler::connectGluonSmearingHandler");
    }
  }
}

void QuarkHandler::disconnectGluonSmearingHandler() {
  gSmearCounter--;
  if (gSmearCounter == 0) {
    try {
      gSmearHandler.reset();
    } catch (const std::exception &xp) {
      errorLaph(
          "delete problem in QuarkHandler::disconnectGluonSmearingHandler");
    }
  }
}

void QuarkHandler::connectQuarkSmearingHandler(
    const std::string &smeared_quark_filestub) {
  if ((qSmearCounter == 0) && (qSmearHandler.get() == 0)) {
    try {
      qSmearHandler.reset(new QuarkSmearingHandler(
          *gSmearPtr, *uPtr, *qSmearPtr, smeared_quark_filestub));
      qSmearCounter = 1;
    } catch (const std::exception &xp) {
      errorLaph(
          "allocation problem in QuarkHandler::connectQuarkSmearingHandler");
    }
  } else {
    try {
      if (qSmearHandler.get() == 0)
        throw(std::invalid_argument("error"));
      uPtr->checkEqual(qSmearHandler->getGaugeConfigurationInfo());
      qSmearHandler->updateSmearing(*qSmearPtr); // increase eigvecs if needed
      qSmearCounter++;
    } catch (const std::exception &xp) {
      errorLaph("inconsistent QuarkHandler::connectQuarkSmearingHandler");
    }
  }
  qSmearHandler->checkAllLevelFilesExist();
}

void QuarkHandler::disconnectQuarkSmearingHandler() {
  qSmearCounter--;
  if (qSmearCounter == 0) {
    try {
      qSmearHandler.reset();
    } catch (const std::exception &xp) {
      errorLaph(
          "delete problem in QuarkHandler::disconnectQuarkSmearingHandler");
    }
  }
}

void QuarkHandler::connectDilutionHandler() const {
  if (dilHandler != 0) {
    errorLaph("QuarkHandler::connectDilutionHandler already connected");
  }
  try {
    dilHandler = new DilutionHandler(*dilPtr, *qSmearPtr);
  } catch (const std::exception &xp) {
    errorLaph("allocation problem in QuarkHandler::connectDilutionHandler");
  }
}

void QuarkHandler::disconnectDilutionHandler() const {
  try {
    delete dilHandler;
  } catch (const std::exception &xp) {
    errorLaph("delete problem in QuarkHandler::disconnectDilutionHandler");
  }
  dilHandler = 0;
}

void QuarkHandler::setNormalMode() {
  if (normal_mode)
    return;
  clearData();
  normal_mode = true;
}

void QuarkHandler::setGamma5HermiticityMode() {
  if (!normal_mode)
    return;
  clearData();
  normal_mode = false;
}

std::map<int, QuarkHandler::FileKey> QuarkHandler::getSuffixMap() const {
  check_info_set("getSuffixMap");
  if (compute_mode) {
    return DHputPtr->getSuffixMap();
  } else {
    return DHgetPtr->getSuffixMap();
  }
}

void QuarkHandler::outputSuffixMap() {
  check_info_set("getSuffixMap");
  std::map<int, QuarkHandler::FileKey> suffixmap = getSuffixMap();
  printLaph("\nSuffix map:");
  for (std::map<int, QuarkHandler::FileKey>::const_iterator it = suffixmap.begin();
       it != suffixmap.end(); ++it) {
    printLaph(make_strf(
        "suffix  %d:  LaphNoiseInfo seed = %d  time proj index = %d", it->first,
        it->second.noise.getSeed(), it->second.time_proj_index));
  }
  printLaph("\n");
}

bool QuarkHandler::isInfoSet() const {
  return ((uPtr != 0) && (gSmearPtr != 0) && (qSmearPtr != 0) && (fPtr != 0) &&
          (dilPtr != 0) && (qactionPtr != 0));
}

void QuarkHandler::check_info_set(const std::string &name) const {
  if (!isInfoSet()) {
    errorLaph(make_strf("error in QuarkHandler: must setInfo before calling %s",
                        name));
  }
}

bool QuarkHandler::isComputeReady() const {
  return ((isInfoSet()) && (compute_mode) && (invertPtr != 0));
}

void QuarkHandler::check_compute_ready(const std::string &name) const {
  if (!isComputeReady()) {
    errorLaph(make_strf("error in QuarkHandler: must setInfo before calling %s",
                        name));
  }
}

const GaugeConfigurationInfo &QuarkHandler::getGaugeConfigurationInfo() const {
  check_info_set("getGaugeConfigurationInfo");
  return *uPtr;
}

const GluonSmearingInfo &QuarkHandler::getGluonSmearingInfo() const {
  check_info_set("getGluonSmearingInfo");
  return *gSmearPtr;
}

const QuarkSmearingInfo &QuarkHandler::getQuarkSmearingInfo() const {
  check_info_set("getQuarkSmearingInfo");
  return *qSmearPtr;
}

const DilutionSchemeInfo &QuarkHandler::getDilutionSchemeInfo() const {
  check_info_set("getDilutionSchemeInfo");
  return *dilPtr;
}

const QuarkActionInfo &QuarkHandler::getQuarkActionInfo() const {
  check_info_set("getQuarkActionInfo");
  return *qactionPtr;
}

const FileListInfo &QuarkHandler::getFileListInfo() const {
  check_info_set("getFileListInfo");
  return *fPtr;
}

int QuarkHandler::getNumberOfSpinEigvecDilutionProjectors() const {
  check_info_set("getNumberOfSpinEigvecDilutionProjectors");
  return dilHandler->getNumberOfSpinEigvecProjectors();
}

int QuarkHandler::getNumberOfTimeDilutionProjectors() const {
  check_info_set("getNumberOfTimeDilutionProjectors");
  return dilHandler->getNumberOfTimeProjectors();
}

int QuarkHandler::getTimeDilutionProjectorIndex(const int time_val) const {
  check_info_set("getTimeDilutionProjectorIndex");
  return dilHandler->getTimeProjectorIndex(time_val);
}

const std::list<int> &QuarkHandler::getOnTimes(const int time_proj_index) const {
  check_info_set("getOnTimes");
  return dilHandler->getOnTimes(time_proj_index);
}

bool QuarkHandler::isFullTimeDilution() const {
  check_info_set("isFullTimeDilution");
  return dilHandler->isFullTimeDilution();
}

int QuarkHandler::getTimeExtent() const {
  check_info_set("getTimeExtent");
  return uPtr->getTimeExtent();
}

bool QuarkHandler::checkHeader(XMLHandler &xmlin, const int suffix) {
  XMLHandler xml_in(xmlin);
  if (xml_tag_count(xml_in, "QuarkHandlerDataFile") != 1)
    return false;
  XMLHandler xmlr(xml_in, "QuarkHandlerDataFile");
  GaugeConfigurationInfo gauge_check(xmlr);
  GluonSmearingInfo gsmear_check(xmlr);
  QuarkSmearingInfo qsmear_check(xmlr);
  DilutionSchemeInfo dil_check(xmlr);
  QuarkActionInfo qaction_check(xmlr);
  try {
    uPtr->checkEqual(gauge_check);
    gSmearPtr->checkEqual(gsmear_check);
    qSmearPtr->checkEqual(qsmear_check);
    dilPtr->checkEqual(dil_check);
    qactionPtr->checkEqual(qaction_check);
  } catch (const std::exception &xp) {
    return false;
  }
  return true;
}

void QuarkHandler::writeHeader(XMLHandler &xmlout,
                               const QuarkHandler::FileKey &fkey,
                               const int suffix) {
  xmlout.set_root("QuarkHandlerDataFile");
  XMLHandler xmltmp;
  uPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  gSmearPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  qSmearPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  dilPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  qactionPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  fkey.output(xmltmp);
  xmlout.put_child(xmltmp);
}

void QuarkHandler::setUpPreconditioning(QudaInvertParam &invParam) {
  if (invertPtr->getName() == "GCR_MULTIGRID") {
    if (invertPtr->QudaMGInfoPtr) {
      if (preconditioner) {
        destroyMultigridQuda(preconditioner);
      }
      printLaph("Initializing the fabulous Multigrid preconditioner!");
      preconditioner = newMultigridQuda(&invertPtr->QudaMGInfoPtr->mg_param);
      invParam.preconditioner = preconditioner;
    } else {
      throw(std::invalid_argument(
          "setUpPreconditioning failed since invParam not set"));
    }
  }
}

void QuarkHandler::computeSinks(const bool verbose,
                                const bool extra_soln_check) {
  check_compute_ready("computeSink");
  StopWatch totaltime;
  totaltime.start();

  printLaph("Task computeSinks beginnning\n");
  printLaph("Quda Projection batch parameters:");
  printLaph(make_strf("      Number of inversions before projecting is %d",
                      sinkComps.nSinkLaphBatch));
  printLaph(make_strf("  Number of solutions projected in one batch is %d",
                      sinkComps.nSinkQudaBatch));
  printLaph(make_strf("    Number of eigvecs projected in one batch is %d",
                      sinkComps.nEigQudaBatch));

  // check temporal boundary conditions in InverterInfo and
  // GaugeConfigurationInfo are the same
  const bool tbc1 = qactionPtr->isFermionTimeBCAntiPeriodic();
  const bool tbc2 = uPtr->isFermionTimeBCAntiPeriodic();
  if (tbc1 != tbc2) {
    errorLaph("Inconsistent fermion time boundary conditions in "
              "QuarkActionInfo and GaugeConfigurationInfo",
              true);
  }

  // load the gauge configuration onto host and device
  StopWatch rolex;
  rolex.start();
  gaugeHandler->setData();
  // quda hack to handle antifermionic temporal boundary conditions
  const bool fermbchack = qactionPtr->isFermionTimeBCAntiPeriodic();
  if (fermbchack) {
    gaugeHandler->eraseDataOnDevice(); // erase on device to apply b.c. below to
                                       // host, then copy to device
    gaugeHandler->applyFermionTemporalAntiPeriodic();
  }
  gaugeHandler->copyDataToDevice();
  // undo the hack on the host
  if (fermbchack) {
    gaugeHandler->applyFermionTemporalAntiPeriodic();
  }
  rolex.stop();
  double grtime = rolex.getTimeInSeconds();
  printLaph("...gauge configuration loaded on host and device");
  printLaph(
      make_str(" time to load gauge configuration was ", grtime, " seconds"));

  // load the clover term if needed
  double clovertime = 0.0;
  if (quda_inv_param.dslash_type == QUDA_CLOVER_WILSON_DSLASH) {
    if (!QudaInfo::clover_on_device) {
      rolex.reset();
      rolex.start();
      loadCloverQuda(
          NULL, NULL,
          &quda_inv_param); // allocates space for and compute clover field
      printLaph("...clover term set up and loaded on device");
      printLaph("...clover term should be positive-definite, then clover "
                "inverse calculated using Cholesky");
      rolex.stop();
      clovertime = rolex.getTimeInSeconds();
      printLaph(
          make_str(" Time to set up clover term was ", clovertime, " seconds"));
      QudaInfo::clover_on_device = true;
    }
  }

  // set up any preconditioning in the inverter
  setUpPreconditioning(quda_inv_param);

  // load the LapH eigenvectors onto host, store pointers
  rolex.reset();
  rolex.start();
  const int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
  std::vector<void *> evList(nEigs);
  for (int n = 0; n < nEigs; n++) {
    evList[n] =
        (void *)(qSmearHandler->getLaphEigenvector(n).getDataConstPtr());
    printLaph(make_strf("read LapH eigvec level %d", n));
  }
  qSmearHandler->closeLaphLevelFiles();
  rolex.stop();
  double evreadtime = rolex.getTimeInSeconds();
  printLaph(make_str("All Laph eigvecs read in ", evreadtime, " seconds\n"));

  const int ncomp = sinkComps.computations.size();
  double srctime = 0.0, invtime = 0.0, evprojtime = 0.0, writetime = 0.0;
  int count = 0;

  for (std::list<SinkComputation>::const_iterator it =
           sinkComps.computations.begin();
       it != sinkComps.computations.end(); count++, it++) {
    printLaph(
        "\n\n *************************************************************");
    printLaph(make_strf(
        " *\n *  Now starting sink computation %d (with %d as last):", count,
        ncomp - 1));
    printLaph(" *");
    printLaph(" *************************************************************");

    computeSinks(it->Noise, it->TimeProjIndex, evList, verbose,
                 extra_soln_check, srctime, invtime, evprojtime, writetime);
  }

  totaltime.stop();
  printLaph("\n\n");
  printLaph("computeQuarkSink: ran successfully");
  printLaph(make_str("                       Total time = ",
                     totaltime.getTimeInSeconds(), " seconds"));
  printLaph(
      make_str(" Time to load gauge configuration = ", grtime, " seconds"));
  printLaph(
      make_str("            LapH eigvec read time = ", evreadtime, " seconds"));
  if (quda_inv_param.dslash_type == QUDA_CLOVER_WILSON_DSLASH) {
    printLaph(make_str("       Time to set up clover term = ", clovertime,
                       " seconds"));
  }
  printLaph(
      make_str("         Total source set up time = ", srctime, " seconds"));
  printLaph(
      make_str("           Total time in inverter = ", invtime, " seconds"));
  printLaph(
      make_str("Projection onto LapH eigvecs time = ", evprojtime, " seconds"));
  printLaph(make_str("          Total file writing time = ", writetime,
                     " seconds\n\n"));
}
//  do inversions for all spin-laph-eigenvector dilution indices
//  but for one noise, one time projector index

void QuarkHandler::computeSinks(const LaphNoiseInfo &noise,
                                const int time_proj_index,
                                const std::vector<void *> &evList,
                                const bool verbose, const bool extra_soln_check,
                                double &src_time, double &inv_time,
                                double &evproj_time, double &write_time) {
  if (!dilHandler->isValidTimeProjectorIndex(time_proj_index)) {
    printLaph(
        make_strf("invalid time projector index %d for compute in QuarkHandler",
                  time_proj_index));
    printLaph(" ...skipping this computation");
    return;
  }

  StopWatch bulova;
  bulova.start();
  printLaph("\nQuark sink computation for all dilutions, one noise,");
  printLaph(" one time dilution projector beginning");
  printLaph(make_strf(" Time dilution projector index = %d", time_proj_index));

  const int Textent = uPtr->getTimeExtent();
  const int Nspin = FieldNspin;
  const int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
  const int ndil = dilHandler->getNumberOfSpinEigvecProjectors();
  const int minTime = 0;
  const int maxTime = Textent - 1;
  const uint nSinkLaphBatch = sinkComps.nSinkLaphBatch;
  const uint nSinkQudaBatch = sinkComps.nSinkQudaBatch;
  const uint nEigQudaBatch = sinkComps.nEigQudaBatch;

  // generate the source noise field
  LaphZnNoise rho(noise.getZNGroup(), noise.getSeed(*uPtr));
  Array<cmplx> laph_noise =
      rho.generateLapHQuarkSourceForSink(Textent, Nspin, nEigs);
  // time dilution masks
  const std::list<int> &on_times = dilHandler->getOnTimes(time_proj_index);
  const double soln_rescale = qactionPtr->getSolutionRescaleFactor();

  // allocate space for batched solutions and make pointers suitable for quda
  int iSinkBatch = 0;
  std::vector<int> sinkBatchInds(nSinkLaphBatch);
  std::vector<LattField> sinkBatchData(nSinkLaphBatch,
                                  FieldSiteType::ColorSpinVector);
  std::vector<void *> sinkList(nSinkLaphBatch);
  for (int iSink = 0; iSink < int(nSinkLaphBatch); iSink++) {
    sinkList[iSink] = (void *)(sinkBatchData[iSink].getDataPtr());
  }
  void **sinks_ptr = (void **)sinkList.data();
  void **evs_ptr = (void **)evList.data();

  // get the file key and open file for writing

  FileKey fkey(noise, time_proj_index);
  DHputPtr->open(fkey);
  bulova.stop();
  double srctime = bulova.getTimeInSeconds();
  double invtime = 0.0, writetime = 0.0, evprojtime = 0.0;

  // loop over dilutions
  for (int dil = 0; dil < ndil; dil++) {

    printLaph(
        make_strf("\nStarting dilution %d (with %d as last)", dil, ndil - 1));
    bool doneflag = true;
    for (int t = minTime; t <= maxTime; t++)
      for (int s = 1; s <= Nspin; s++) {
        if (!DHputPtr->queryData(RecordKey(s, t, dil))) {
          doneflag = false;
          s = Nspin + 1;
          t = Textent;
        }
      }

    if ((!fPtr->isModeOverwrite()) && (doneflag)) {
      printLaph("warning: these quark sinks already computed...");
      printLaph("  skip re-computing since fileMode not overwrite");
    } else {

      // get the lists of which spins and which eigenvectors are
      // "on" for this dilution projector

      bulova.reset();
      bulova.start();
      const std::list<int> &on_spins = dilHandler->getOnSpinIndices(dil);
      const std::list<int> &on_eigs = dilHandler->getOnEigvecIndices(dil);

      //  initialize source (include gamma_4) in Dirac-Pauli basis, allocate
      //  sink
      LattField ferm_src;
      make_source(ferm_src, laph_noise, evList, on_times, on_spins, on_eigs);
      bulova.stop();
      double addsrctime = bulova.getTimeInSeconds();
      printLaph(make_str(" fermion source set up in Dirac-Pauli basis in ",
                         addsrctime, " seconds"));
      srctime += addsrctime;

      // now do the inversion (source and sink should be on host in Dirac-Pauli
      // basis)
      bulova.reset();
      bulova.start();
      void *spinor_src = (void *)(ferm_src.getDataPtr());
      void *spinor_snk = (void *)(sinkBatchData[iSinkBatch].getDataPtr());

      invertQuda(spinor_snk, spinor_src, &quda_inv_param);

      bulova.stop();
      double addinvtime = bulova.getTimeInSeconds();
      invtime += addinvtime;

      printLaph(make_strf("Inversion done:  number of iterations = %d",
                          quda_inv_param.iter));
      printLaph(make_strf(" Residual = %g", quda_inv_param.true_res[0]));
      printLaph(make_str("  Inversion time was ", addinvtime, " seconds"));
      printLaph(make_str("  Inverter used was ", invertPtr->getName()));

      if (quda_inv_param.true_res[0] > quda_inv_param.tol) {
        printLaph("Inverter failed to reach required tolerance");
      }

      // additional checks of the solution

      if (extra_soln_check) {
        printLaph(
            "Performing additional check on solution using quda::MatQuda");
        LattField sol_check(FieldSiteType::ColorSpinVector);
        void *sol_checkptr = sol_check.getDataPtr();
        MatQuda(sol_checkptr, spinor_snk, &quda_inv_param);
        compare_latt_fields(sol_check, ferm_src);
      }

      if ((quda_inv_param.iter > 0) &&
          (quda_inv_param.iter < int(invertPtr->getMaxIterations()))) {
        sinkBatchInds[iSinkBatch] = dil;
        iSinkBatch++;
      } else {
        printLaph(
            "\n\nInversion FAILED to converge before max iteration reached");
        printLaph("Solution NOT WRITTEN to file\n\n");
      }
    }

    // carry out projections
    if ((iSinkBatch == int(nSinkLaphBatch)) ||
        ((dil == (ndil - 1)) && (iSinkBatch > 0))) {

      bulova.reset();
      bulova.start();
      const int nSinks = iSinkBatch;
      int nSinksBatch = std::min(nSinks, int(nSinkQudaBatch));

      Array<dcmplx> qudaRes(Nspin, Textent, nEigs,
                            nSinks); // quda laph reversing major order
      __complex__ double *qudaResPtr =
          (__complex__ double *)(&qudaRes(0, 0, 0, 0));

      // do the projections
      printLaph("projecting batch of solutions onto LapH eigenvectors");
      laphSinkProject(qudaResPtr, sinks_ptr, nSinks, nSinksBatch, evs_ptr,
                      nEigs, nEigQudaBatch, &quda_inv_param,
                      LayoutInfo::getRankLattExtents().data());

      bulova.stop();
      double addevprojtime = bulova.getTimeInSeconds();
      printLaph(make_str(" this batch of projections onto Laph evs took ",
                         addevprojtime, " seconds"));
      evprojtime += addevprojtime;

      // rearrange data then output to file
      bulova.reset();
      bulova.start();
      for (int iSink = 0; iSink < nSinks; ++iSink) {
        for (int t = minTime; t <= maxTime; t++) {
          for (int iSpin = 0; iSpin < Nspin; ++iSpin) {
	    std::vector<dcmplx> quark_sink(nEigs);
            for (int iEv = 0; iEv < nEigs; ++iEv) {
              quark_sink[iEv] = soln_rescale * qudaRes(iSpin, t, iEv, iSink);
            }
            DHputPtr->putData(RecordKey(iSpin + 1, t, sinkBatchInds[iSink]),
                              quark_sink);
            if (verbose) {
              printLaph(make_strf("dil = %d, spin = %d, time = %d",
                                  sinkBatchInds[iSink], iSpin + 1, t));
              for (int n = 0; n < nEigs; n++) {
                printLaph(make_strf("coef for eigenlevel %d = (%14.8f, %14.8f)",
                                    n, real(quark_sink[n]),
                                    imag(quark_sink[n])));
              }
            }
          }
        }
      }
      bulova.stop();
      double otime = bulova.getTimeInSeconds();
      printLaph(
          make_str(" Output of this batch to file took ", otime, " seconds"));
      writetime += otime;

      iSinkBatch = 0;
    } // batch end
    DHputPtr->flush();
  } // dilution loop end

  inv_time += invtime;
  src_time += srctime;
  evproj_time += evprojtime;
  write_time += writetime;
}

// Creates the source (multiplied by gamma_4) in Dirac-Pauli basis

void QuarkHandler::make_source(LattField &ferm_src,
                               const Array<cmplx> &laph_noise,
                               const std::vector<void *> &evList,
                               const std::list<int> &on_times,
                               const std::list<int> &on_spins,
                               const std::list<int> &on_eigs) {
  printLaph(" Making source for this inversion...");
  if (evList.size() == 0) {
    errorLaph(
        "Cannot make fermion source since no Laph eigenvectors available");
  }
  ferm_src.reset(FieldSiteType::ColorSpinVector);
  const bool dp = (ferm_src.bytesPerWord() == sizeof(std::complex<double>));
  // initialize source field to zero
  const int loc_nsites = LayoutInfo::getRankLatticeNumSites();
  const int ncmplx_per_site = ferm_src.elemsPerSite();
  const int ncmplx = ncmplx_per_site * loc_nsites;
  int cbytes;
  dcmplx zrhodp;
  fcmplx zrhosp;
  char *zrho;
  if (dp) {
    double *z0 = reinterpret_cast<double *>(ferm_src.getDataPtr());
    std::fill(z0, z0 + 2 * ncmplx, 0.0);
    cbytes = sizeof(std::complex<double>);
    zrho = reinterpret_cast<char *>(&zrhodp);
  } else {
    float *z0 = reinterpret_cast<float *>(ferm_src.getDataPtr());
    std::fill(z0, z0 + 2 * ncmplx, 0.0);
    cbytes = sizeof(std::complex<float>);
    zrho = reinterpret_cast<char *>(&zrhosp);
  }

  const int loc_npsites = loc_nsites / 2;
  const int start_parity = LayoutInfo::getMyStartParity();
  const int mytmin =
      LayoutInfo::getMyCommCoords()[3] * LayoutInfo::getRankLattExtents()[3];
  const int mytmax = mytmin + LayoutInfo::getRankLattExtents()[3] - 1;
  const int tstride = LayoutInfo::getRankLattExtents()[0] *
                      LayoutInfo::getRankLattExtents()[1] *
                      LayoutInfo::getRankLattExtents()[2];
  const int incx = FieldNcolor;
  const int incy = FieldNcolor * FieldNspin;

  for (std::list<int>::const_iterator it0 = on_times.begin(); it0 != on_times.end();
       it0++) {
    if (((*it0) >= mytmin) && ((*it0) <= mytmax)) {

      const int tloc = ((*it0) - mytmin);
      int parshift = loc_npsites * ((start_parity + tloc) % 2);
      const int start1 = ((tstride * tloc) / 2) + parshift;
      const int stop1 = ((1 + tstride * (tloc + 1)) / 2) + parshift;
      const int n1 = stop1 - start1;

      parshift = loc_npsites * ((start_parity + 1 + tloc) % 2);
      const int start2 = ((1 + tstride * tloc) / 2) + parshift;
      const int stop2 = ((tstride * (tloc + 1)) / 2) + parshift;
      const int n2 = stop2 - start2;
      const int xstart1 = start1 * incx * cbytes;
      const int xstart2 = start2 * incx * cbytes;
      char *ystart1 = ferm_src.getDataPtr() + start1 * incy * cbytes;
      char *ystart2 = ferm_src.getDataPtr() + start2 * incy * cbytes;

      for (std::list<int>::const_iterator vmask = on_eigs.begin();
           vmask != on_eigs.end(); vmask++) {
        char *x0 = reinterpret_cast<char *>(evList[*vmask]);
        for (std::list<int>::const_iterator smask = on_spins.begin();
             smask != on_spins.end(); smask++) {
          zrhodp = laph_noise(*it0, *smask, *vmask);
          if (*smask > 1) {
            zrhodp = -zrhodp;
          } // multiply by gamma_4
          if (!dp) {
            zrhosp = std::complex<float>(real(zrhodp), imag(zrhodp));
          }
          char *x1 = x0 + xstart1;
          char *x2 = x0 + xstart2;
          char *y1 = ystart1 + (*smask) * incx * cbytes;
          char *y2 = ystart2 + (*smask) * incx * cbytes;
          for (int c = 0; c < FieldNcolor; ++c) {
            // branching in the hot part of the code is inneficient
            if (dp) {
              cblas_zaxpy(n1, (dcmplx *)(zrho), (dcmplx *)(x1), incx,
                          (dcmplx *)(y1), incy);
              cblas_zaxpy(n2, (dcmplx *)(zrho), (dcmplx *)(x2), incx,
                          (dcmplx *)(y2), incy);
            } else {
              cblas_caxpy(n1, (fcmplx *)(zrho), (fcmplx *)(x1), incx,
                          (fcmplx *)(y1), incy);
              cblas_caxpy(n2, (fcmplx *)(zrho), (fcmplx *)(x2), incx,
                          (fcmplx *)(y2), incy);
            }
            x1 += cbytes;
            y1 += cbytes;
            x2 += cbytes;
            y2 += cbytes;
          }
        }
      }
    }
  }
  printLaph("Source for this inversion created");
}

void QuarkHandler::clearData() {
  if (compute_mode)
    return;
  for (std::map<StorageKey, LattField *>::iterator it = store_sources.begin();
       it != store_sources.end(); ++it) {
    if (it->second != 0)
      delete it->second;
  }
  store_sources.clear();
  for (std::map<StorageKey, LattField *>::iterator it = store_sinks.begin();
       it != store_sinks.end(); ++it) {
    delete it->second;
  }
  store_sinks.clear();
  DHgetPtr->clearData();
}

//  static pointers (set to null in default constructor)
std::unique_ptr<GluonSmearingHandler> QuarkHandler::gSmearHandler;
std::unique_ptr<QuarkSmearingHandler> QuarkHandler::qSmearHandler;
std::unique_ptr<GaugeConfigurationHandler> QuarkHandler::gaugeHandler;

int QuarkHandler::gSmearCounter = 0;
int QuarkHandler::qSmearCounter = 0;
int QuarkHandler::gaugeCounter = 0;
} // namespace LaphEnv
