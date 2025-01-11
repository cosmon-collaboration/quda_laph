#include "QudaLaphIncludes.h"

#include "perambulator_handler.h"
#include "field_ops.h"
#include "multi_compare.h"

#include <cassert>

namespace LaphEnv {

void PerambulatorHandler::RecordKey::output(XMLHandler &xmlw) const {
  xmlw.set_root("RecordKey");
  xmlw.put_child("SinkSpin", make_string(getSinkSpin()));
  xmlw.put_child("SinkTime", make_string(getSinkTime()));
  xmlw.put_child("SourceLaphEigvecIndex",
                 make_string(getSourceLaphEigvecIndex()));
}

PerambulatorHandler::FileKey::FileKey(const int in_srctime,
                                      const int in_srcspin)
    : src_time(in_srctime), src_spin(in_srcspin) {}

PerambulatorHandler::FileKey::FileKey(XMLHandler &xmlr) {
  try {
    XMLHandler xmlf(xmlr, "FileKey");
    xmlread(xmlf, "SourceTime", src_time, "PerambulatorHandler::FileKey");
    xmlread(xmlf, "SourceSpin", src_spin, "PerambulatorHandler::FileKey");
  } catch (...) {
    errorLaph("Could not read PerambulatorHandler::FileKey");
  }
}

PerambulatorHandler::FileKey::FileKey(const FileKey &rhs)
    : src_time(rhs.src_time), src_spin(rhs.src_spin) {}

PerambulatorHandler::FileKey &
PerambulatorHandler::FileKey::operator=(const FileKey &rhs) {
  src_time = rhs.src_time;
  src_spin = rhs.src_spin;
  return *this;
}

void PerambulatorHandler::FileKey::output(XMLHandler &xmlw) const {
  xmlw.set_root("FileKey");
  xmlw.put_child("SourceTime", make_string(src_time));
  xmlw.put_child("SourceSpin", make_string(src_spin));
}

bool PerambulatorHandler::FileKey::operator<(const FileKey &rhs) const {
  return multiLessThan(src_time, rhs.src_time, src_spin, rhs.src_spin);
}

bool PerambulatorHandler::FileKey::operator==(const FileKey &rhs) const {
  return multiEqual(src_time, rhs.src_time, src_spin, rhs.src_spin);
}

bool PerambulatorHandler::FileKey::operator!=(const FileKey &rhs) const {
  return multiNotEqual(src_time, rhs.src_time, src_spin, rhs.src_spin);
}

PerambulatorHandler::PerambulatorHandler()
    : uPtr(0), gSmearPtr(0), qSmearPtr(0), qactionPtr(0), fPtr(0), invertPtr(0),
      Nspin(4), mode(ReadOnly), preconditioner(0), DHputPtr(0), DHgetPtr(0) {}

PerambulatorHandler::PerambulatorHandler(
    const GaugeConfigurationInfo &gaugeinfo,
    const GluonSmearingInfo &gluonsmear, const QuarkSmearingInfo &quarksmear,
    const QuarkActionInfo &quark, const FileListInfo &flist,
    const std::string &smeared_quark_filestub, const bool upper_spin_components_only,
    const Mode in_mode, const std::string &gauge_str)
    : invertPtr(0), preconditioner(0), DHputPtr(0), DHgetPtr(0) {
  set_info(gaugeinfo, gluonsmear, quarksmear, quark, flist,
           smeared_quark_filestub, upper_spin_components_only, gauge_str,
           in_mode);
}

void PerambulatorHandler::setInfo(const GaugeConfigurationInfo &gaugeinfo,
                                  const GluonSmearingInfo &gluonsmear,
                                  const QuarkSmearingInfo &quarksmear,
                                  const QuarkActionInfo &quark,
                                  const FileListInfo &flist,
                                  const std::string &smeared_quark_filestub,
                                  const bool upper_spin_components_only,
                                  const Mode in_mode, const std::string &gauge_str) {
  clear();
  set_info(gaugeinfo, gluonsmear, quarksmear, quark, flist,
           smeared_quark_filestub, upper_spin_components_only, gauge_str,
           in_mode);
}

void PerambulatorHandler::set_info(
    const GaugeConfigurationInfo &gaugeinfo,
    const GluonSmearingInfo &gluonsmear,
    const QuarkSmearingInfo &quarksmear,
    const QuarkActionInfo &quark,
    const FileListInfo &flist,
    const std::string &smeared_quark_filestub,
    const bool upper_spin_components_only,
    const std::string &gauge_str,
    const Mode in_mode) {
  try {
    uPtr = new GaugeConfigurationInfo(gaugeinfo);
    gSmearPtr = new GluonSmearingInfo(gluonsmear);
    qSmearPtr = new QuarkSmearingInfo(quarksmear);
    qactionPtr = new QuarkActionInfo(quark);
    fPtr = new FileListInfo(flist);
    Nspin = (upper_spin_components_only) ? 2 : 4;
    mode = in_mode;
    if (mode != ReadOnly) {
      DHputPtr = new DataPutHandlerMF<PerambulatorHandler, FileKey, RecordKey,
                                      DataType>(
          *this, *fPtr, "Laph--QuarkPeramb", "PerambulatorHandlerDataFile");
    }
    if (mode != Compute) {
      DHgetPtr = new DataGetHandlerMF<PerambulatorHandler, FileKey, RecordKey,
                                      DataType>(
          *this, *fPtr, "Laph--QuarkPeramb", "PerambulatorHandlerDataFile");
    }
  } catch (const std::exception &xp) {
    errorLaph(
        make_strf("allocation problem in PerambulatorHandler: %s", xp.what()));
  }

  if (mode == Compute) {
    connectGaugeConfigurationHandler();
  }
  connectQuarkSmearingHandler(smeared_quark_filestub);
}

PerambulatorHandler::~PerambulatorHandler() { clear(); }

void PerambulatorHandler::clear() {
  try {
    delete uPtr;
    delete gSmearPtr;
    delete qSmearPtr;
    delete qactionPtr;
    delete fPtr;
    delete invertPtr;
    if (preconditioner) {
      destroyMultigridQuda(preconditioner);
    }
    clearComputationSet();
  } catch (const std::exception &xp) {
    errorLaph("abort");
  }
  uPtr = 0;
  gSmearPtr = 0;
  qSmearPtr = 0;
  qactionPtr = 0;
  fPtr = 0;
  invertPtr = 0;
  mode = ReadOnly;
  preconditioner = 0;
  disconnectGaugeConfigurationHandler();
  disconnectQuarkSmearingHandler();
  delete DHgetPtr;
  DHgetPtr = 0;
  delete DHputPtr;
  DHputPtr = 0;
}

void PerambulatorHandler::connectGaugeConfigurationHandler() {
  if ((gaugeCounter == 0) && (gaugeHandler.get() == 0)) {
    try {
      gaugeHandler.reset(new GaugeConfigurationHandler(*uPtr));
      gaugeCounter = 1;
    } catch (const std::exception &xp) {
      errorLaph("allocation problem in "
                "PerambulatorHandler::connectGaugeConfigurationHandler");
    }
  } else {
    try {
      if (gaugeHandler.get() == 0)
        throw(std::invalid_argument("error"));
      uPtr->checkEqual(gaugeHandler->getGaugeConfigurationInfo());
      gaugeCounter++;
    } catch (const std::exception &xp) {
      errorLaph(
          "inconsistent PerambulatorHandler::connectGaugeConfigurationHandler");
    }
  }
}

void PerambulatorHandler::disconnectGaugeConfigurationHandler() {
  gaugeCounter--;
  if (gaugeCounter == 0) {
    try {
      gaugeHandler.reset();
    } catch (const std::exception &xp) {
      errorLaph("delete problem in "
                "PerambulatorHandler::disconnectGluonSmearingHandler");
    }
  }
}

void PerambulatorHandler::connectQuarkSmearingHandler(
    const std::string &smeared_quark_filestub) {
  if ((qSmearCounter == 0) && (qSmearHandler.get() == 0)) {
    try {
      qSmearHandler.reset(new QuarkSmearingHandler(
          *gSmearPtr, *uPtr, *qSmearPtr, smeared_quark_filestub));
      qSmearCounter = 1;
    } catch (const std::exception &xp) {
      errorLaph("allocation problem in "
                "PerambulatorHandler::connectQuarkSmearingHandler");
    }
  } else {
    try {
      if (qSmearHandler.get() == 0)
        throw(std::invalid_argument("error"));
      uPtr->checkEqual(qSmearHandler->getGaugeConfigurationInfo());
      qSmearHandler->updateSmearing(*qSmearPtr); // increase eigvecs if needed
      qSmearCounter++;
    } catch (const std::exception &xp) {
      errorLaph(
          "inconsistent PerambulatorHandler::connectQuarkSmearingHandler");
    }
  }
  qSmearHandler->checkAllLevelFilesExist();
}

void PerambulatorHandler::disconnectQuarkSmearingHandler() {
  qSmearCounter--;
  if (qSmearCounter == 0) {
    try {
      qSmearHandler.reset();
    } catch (const std::exception &xp) {
      errorLaph("delete problem in "
                "PerambulatorHandler::disconnectQuarkSmearingHandler");
    }
  }
}

void PerambulatorHandler::setInverter(const InverterInfo &invinfo) {
  if (mode != Compute) {
    errorLaph(
        "Cannot setInverter in PerambulatorHandler unless in compute mode");
  }
  if (!isInfoSet()) {
    errorLaph("Cannot setInverter in PerambulatorHandler unless info is set");
  }
  try {
    delete invertPtr;
    invertPtr = new InverterInfo(invinfo);
  } catch (const std::exception &xp) {
    errorLaph("allocation error in PerambulatorHandler::setInverter");
  }

  // set up the inverter params for quda
  invertPtr->setQudaInvertParam(quda_inv_param, *qactionPtr);

  // check that Dirac-Pauli basis is used
  if (quda_inv_param.gamma_basis != QUDA_DIRAC_PAULI_GAMMA_BASIS) {
    errorQuda("The Dirac-Pauli basis must be used");
  }
}

const InverterInfo &PerambulatorHandler::getInverterInfo() const {
  if (invertPtr != 0) {
    errorLaph("error in PerambulatorHandler: must setInverter before calling "
              "getInverterInfo");
  }
  return *invertPtr;
}

void PerambulatorHandler::clearComputationSet() {
  perambComps.computations.clear();
}

void PerambulatorHandler::setComputationSet(const XMLHandler &xmlin) {
  if (mode != Compute) {
    errorLaph("Cannot setComputationSet in PerambulatorHandler unless in "
              "compute mode");
  }
  if (!isInfoSet()) {
    errorLaph(
        "Cannot setComputationSet in PerambulatorHandler unless info is set");
  }
  XMLHandler xmlrdr(xmlin);
  if (!perambComps.computations.empty())
    perambComps.computations.clear();
  const uint Textent = uPtr->getTimeExtent();
  const uint nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();

  if (xml_tag_count(xmlrdr, "ComputationSet") == 1) {

    XMLHandler xmlrd(xmlrdr, "ComputationSet");

    uint nSinkLaphBatch, nSinkQudaBatch, nEigQudaBatch;
    xmlread(xmlrd, "NumSinksBeforeProject", nSinkLaphBatch,
            "LAPH_PERAMBULATORS");
    xmlread(xmlrd, "NumSinksInProjectBatch", nSinkQudaBatch,
            "LAPH_PERAMBULATORS");
    xmlread(xmlrd, "NumEigsInProjectBatch", nEigQudaBatch,
            "LAPH_PERAMBULATORS");

    if ((nSinkLaphBatch > nEigs * Nspin) || (nSinkLaphBatch == 0)) {
      errorLaph(make_strf("Invalid value %d for <NumSinksBeforeProject>: must "
                          "be between 1 and nEigs*Nspin=%d",
                          nSinkLaphBatch, nEigs * Nspin));
    }
    if ((nEigQudaBatch > nEigs) || (nEigQudaBatch == 0)) {
      errorLaph(make_strf("Invalid value %d for <NumEigsInProjectBatch>: must "
                          "be between 1 and nEigs=%d",
                          nEigQudaBatch, nEigs));
    }
    if ((nSinkQudaBatch > nSinkLaphBatch) || (nSinkQudaBatch == 0)) {
      errorLaph(make_strf("Invalid value %d for <NumSinksInProjectBatch>: must "
                          "be between 1 and %d",
                          nSinkQudaBatch, nSinkLaphBatch));
    }
    perambComps.nSinkLaphBatch = nSinkLaphBatch;
    perambComps.nSinkQudaBatch = nSinkQudaBatch;
    perambComps.nEigQudaBatch = nEigQudaBatch;

    std::list<XMLHandler> xmlcs(xmlrd.find("Computation"));
    for (std::list<XMLHandler>::iterator it = xmlcs.begin(); it != xmlcs.end();
         ++it) {
      int source_time;
      xmlread(*it, "SourceTime", source_time, "LAPH_PERAMBULATORS");
      if ((source_time < 0) || (source_time >= int(Textent))) {
        errorLaph(make_strf("Invalid source time %d", source_time));
      }
      std::set<int> srcev_indices;
      if (xml_tag_count(*it, "SourceLaphEigvecIndices") == 1) {
	std::vector<int> srcev_inds;
        xmlread(*it, "SourceLaphEigvecIndices", srcev_inds,
                "LAPH_PERAMBULATORS");
        for (int k = 0; k < int(srcev_inds.size()); ++k) {
          if ((srcev_inds[k] < 0) || (srcev_inds[k] >= int(nEigs))) {
            errorLaph(
                make_strf("Invalid src laph eigvec index %d", srcev_inds[k]));
          }
          srcev_indices.insert(srcev_inds[k]);
        }
      } else if ((xml_tag_count(*it, "SourceLaphEigvecIndexMin") == 1) &&
                 (xml_tag_count(*it, "SourceLaphEigvecIndexMax") == 1)) {
        int sevmin = -1, sevmax = -1;
        xmlread(*it, "SourceLaphEigvecIndexMin", sevmin,
                "LAPH_PERAMBULATORS");
        xmlread(*it, "SourceLaphEigvecIndexMax", sevmax,
                "LAPH_PERAMBULATORS");
        if (sevmin < 0)
          sevmin = 0;
        if (sevmax >= int(nEigs))
          sevmax = nEigs - 1;
        if (sevmax < sevmin)
          sevmax = sevmin;
        for (int ev = sevmin; ev <= sevmax; ++ev) {
          srcev_indices.insert(ev);
        }
      } else {
        for (int ev = 0; ev < int(nEigs); ++ev) {
          srcev_indices.insert(ev);
        }
      }
      if (srcev_indices.empty()) {
        errorLaph("Empty src laph eigvec indices set");
      }
      perambComps.computations.push_back(
          PerambComputation(source_time, srcev_indices));
    }
  }

  printLaph("\nLAPH_PERAMBULATORS sink computations:\n");
  printLaph(
      make_strf("  NumSinksBeforeProject = %d", perambComps.nSinkLaphBatch));
  printLaph(
      make_strf(" NumSinksInProjectBatch = %d", perambComps.nSinkQudaBatch));
  printLaph(
      make_strf("  NumEigsInProjectBatch = %d\n", perambComps.nEigQudaBatch));
  printLaph(make_strf(" Number of computations = %d",
                      perambComps.computations.size()));
  int count = 0;
  for (std::list<PerambComputation>::const_iterator it =
           perambComps.computations.begin();
       it != perambComps.computations.end(); count++, it++) {
    XMLHandler xmlout("Computation");
    xmlout.put_child("SourceTime", make_string(it->src_time));
    std::string srcevindstr;
    const std::set<int> &indices(it->src_lapheigvec_indices);
    std::set<int>::const_iterator vt = indices.begin();
    int rangestart = *vt, rangestop = *vt + 1;
    ++vt;
    bool keep_going = true;
    while (keep_going) {
      if ((vt != indices.end()) &&
          (rangestop == (*vt))) { // in a started range, increase rangestop
        rangestop++;
        ++vt;
      } else {                               // end of a range
        if (rangestop == (rangestart + 1)) { // only one in range
          srcevindstr += make_string(rangestart);
        } else { // more than one in range
          srcevindstr +=
              make_string(rangestart) + "-" + make_string(rangestop - 1);
        }
        if (vt != indices.end()) {
          rangestart = *vt;
          rangestop = rangestart + 1;
          srcevindstr += " ";
          ++vt;
        } else {
          keep_going = false;
        }
      }
    }
    xmlout.put_child("SourceLaphEigvecIndices", srcevindstr);
    printLaph(make_strf("\nComputation %d:", count));
    printLaph(make_strf("%s", xmlout.output()));
  }
  printLaph("\n");
}

void PerambulatorHandler::getFileMap(XMLHandler &xmlout) const {
}

std::map<int, PerambulatorHandler::FileKey>
PerambulatorHandler::getSuffixMap() const {
  check_info_set("getSuffixMap");
  if (mode != ReadOnly) {
    return DHputPtr->getSuffixMap();
  } else {
    return DHgetPtr->getSuffixMap();
  }
}

void PerambulatorHandler::outputSuffixMap() {
  check_info_set("getSuffixMap");
  std::map<int, PerambulatorHandler::FileKey> suffixmap = getSuffixMap();
  printLaph("\nSuffix map:");
  for (std::map<int, PerambulatorHandler::FileKey>::const_iterator it =
           suffixmap.begin();
       it != suffixmap.end(); ++it) {
    printLaph(make_strf("suffix  %d:  source time = %d  source spin = %d",
                        it->first, it->second.src_time, it->second.src_spin));
  }
  printLaph("\n");
}

bool PerambulatorHandler::isInfoSet() const {
  return ((uPtr != 0) && (gSmearPtr != 0) && (qSmearPtr != 0) && (fPtr != 0) &&
          (qactionPtr != 0));
}

void PerambulatorHandler::check_info_set(const std::string &name) const {
  if (!isInfoSet()) {
    errorLaph(make_strf(
        "error in PerambulatorHandler: must setInfo before calling %s", name));
  }
}

const GaugeConfigurationInfo &
PerambulatorHandler::getGaugeConfigurationInfo() const {
  check_info_set("getGaugeConfigurationInfo");
  return *uPtr;
}

const GluonSmearingInfo &PerambulatorHandler::getGluonSmearingInfo() const {
  check_info_set("getGluonSmearingInfo");
  return *gSmearPtr;
}

const QuarkSmearingInfo &PerambulatorHandler::getQuarkSmearingInfo() const {
  check_info_set("getQuarkSmearingInfo");
  return *qSmearPtr;
}

const QuarkActionInfo &PerambulatorHandler::getQuarkActionInfo() const {
  check_info_set("getQuarkActionInfo");
  return *qactionPtr;
}

const FileListInfo &PerambulatorHandler::getFileListInfo() const {
  check_info_set("getFileListInfo");
  return *fPtr;
}

uint PerambulatorHandler::getNumberOfLaplacianEigenvectors() const {
  return qSmearPtr->getNumberOfLaplacianEigenvectors();
}

int PerambulatorHandler::getTimeExtent() const {
  check_info_set("getTimeExtent");
  return uPtr->getTimeExtent();
}

bool PerambulatorHandler::checkHeader(XMLHandler &xmlin, const int suffix) {
  XMLHandler xml_in(xmlin);
  if (xml_tag_count(xml_in, "PerambulatorHandlerDataFile") != 1)
    return false;
  XMLHandler xmlr(xml_in, "PerambulatorHandlerDataFile");
  GaugeConfigurationInfo gauge_check(xmlr);
  GluonSmearingInfo gsmear_check(xmlr);
  QuarkSmearingInfo qsmear_check(xmlr);
  uint numspin;
  xmlread(xmlr, "NumSpinComponents", numspin, "PerambulatorHandler");
  QuarkActionInfo qaction_check(xmlr);
  try {
    uPtr->checkEqual(gauge_check);
    gSmearPtr->checkEqual(gsmear_check);
    qSmearPtr->checkEqual(qsmear_check);
    if (numspin != Nspin) {
      throw(std::invalid_argument(
          "Perambulator checkEqual failed...NumSpinComponents mismatch"));
    }
    qactionPtr->checkEqual(qaction_check);
  } catch (const std::exception &xp) {
    return false;
  }
  return true;
}

void PerambulatorHandler::writeHeader(XMLHandler &xmlout,
                                      const PerambulatorHandler::FileKey &fkey,
                                      const int suffix) {
  xmlout.set_root("PerambulatorHandlerDataFile");
  XMLHandler xmltmp;
  uPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  gSmearPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  qSmearPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  qactionPtr->output(xmltmp);
  xmlout.put_child(xmltmp);
  fkey.output(xmltmp);
  xmlout.put_child(xmltmp);
  xmlout.put_child("NumSpinComponents", make_string(Nspin));
}

void PerambulatorHandler::setUpPreconditioning(QudaInvertParam &invParam) {
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

//  do "perambulator" inversions for all spin indices, one time source and a set
//  of source laph eigenvector numbers.  All sink spin indices, sink times, and
//  laph_eigvec indices at the sink are computed.

void PerambulatorHandler::computePerambulators(const bool verbose,
                                               const bool extra_soln_check) {
  if ((!isInfoSet()) || (invertPtr == 0) || (mode != Compute)) {
    errorLaph(make_str(
        "cannot computePerambulators in PerambulatorHandler until info ",
        "and inverter set and in compute mode"));
  }
  StopWatch totaltime;
  totaltime.start();

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

  double srctime = 0.0, invtime = 0.0, evprojtime = 0.0, writetime = 0.0;
  int count = 0, ncomp = perambComps.computations.size();

  for (std::list<PerambComputation>::const_iterator it =
           perambComps.computations.begin();
       it != perambComps.computations.end(); count++, it++) {
    printLaph(
        "\n\n *************************************************************");
    printLaph(make_strf(
        " *\n *  Now starting computation %d (with %d as last):", count,
        ncomp - 1));
    printLaph(" *");
    printLaph(" *************************************************************");

    computePerambulators(it->src_time, it->src_lapheigvec_indices, evList,
                         verbose, extra_soln_check, srctime, invtime,
                         evprojtime, writetime);
  }

  totaltime.stop();
  printLaph("\n\n");
  printLaph("computePerambulators: ran successfully");
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

//  do inversions for one source time, one requested selection
//  of source laph-eigenvector dilution indices (all source spins)

void PerambulatorHandler::computePerambulators(const int src_time,
					       const std::set<int> &src_evindices,
					       const std::vector<void *> &evList,
					       const bool verbose,
					       const bool extra_soln_check,
					       double &makesrc_time,
					       double &inv_time,
					       double &evproj_time,
					       double &write_time) {
  StopWatch bulova;
  bulova.start();
  printLaph("\nQuark perambulator computation for one source time,");
  printLaph(" one set of source eigvec indices beginning");
  printLaph(make_strf(" Source time = %d", src_time));

  const int Textent = uPtr->getTimeExtent();
  const int nEigs = qSmearPtr->getNumberOfLaplacianEigenvectors();
  int minTime = 0;
  const int maxTime = Textent - 1;
  const uint nSinkLaphBatch = perambComps.nSinkLaphBatch;
  const uint nSinkQudaBatch = perambComps.nSinkQudaBatch;
  const uint nEigQudaBatch = perambComps.nEigQudaBatch;
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

  bulova.stop();
  double srctime = bulova.getTimeInSeconds();
  double invtime = 0.0;
  double writetime = 0.0;
  double evprojtime = 0.0;
  const int ninv = src_evindices.size();

  // loop over source laph eigvec indices and source spin
  for (int srcspin = 1; srcspin <= int(Nspin); ++srcspin) {

    // get the file key and open file for writing
    FileKey fkey(src_time, srcspin);
    DHputPtr->open(fkey);
    int invcount = 0;

    for (std::set<int>::const_iterator vt = src_evindices.begin();
         vt != src_evindices.end(); ++vt, ++invcount) {

      int srcev_ind = *vt;
      printLaph(make_strf("\nStarting inversion number %d (with %d as last) "
                          "with source spin %d",
                          invcount, ninv - 1, srcspin));
      bool doneflag = true;
      for (int sink_time = 0; sink_time < Textent; ++sink_time) {
        for (int sinkspin = 1; sinkspin <= int(Nspin); ++sinkspin) {
          if (!DHputPtr->queryData(RecordKey(sinkspin, sink_time, srcev_ind))) {
            doneflag = false;
            sinkspin = Nspin + 1;
            sink_time = Textent;
          }
        }
      }

      if ((!fPtr->isModeOverwrite()) && (doneflag)) {
        printLaph("warning: these quark sinks already computed...");
        printLaph("  skip re-computing since fileMode not overwrite");
      } else {

        //  initialize source (include gamma_4) in Dirac-Pauli basis, allocate
        //  sink
        bulova.reset();
        bulova.start();
        LattField ferm_src;
        make_source(ferm_src, evList[srcev_ind], src_time, srcspin);
        bulova.stop();
        double addsrctime = bulova.getTimeInSeconds();
        printLaph(make_str(" fermion source set up in Dirac-Pauli basis in ",
                           addsrctime, " seconds"));
        srctime += addsrctime;

        // now do the inversion (source and sink should be on host in
        // Dirac-Pauli basis)
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
          sinkBatchInds[iSinkBatch] = srcev_ind;
          iSinkBatch++;
        } else {
          printLaph(
              "\n\nInversion FAILED to converge before max iteration reached");
          printLaph("Solution NOT WRITTEN to file\n\n");
        }
      }

      // carry out projections
      if ((iSinkBatch == int(nSinkLaphBatch)) ||
          ((invcount == (ninv - 1)) && (iSinkBatch > 0))) {

        bulova.reset();
        bulova.start();
        int nSinks = iSinkBatch;
        int nSinksBatch = std::min(nSinks, int(nSinkQudaBatch));

        Array<std::complex<double>> qudaRes(Nspin, Textent, nEigs,
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
            for (int iSpin = 0; iSpin < int(Nspin); ++iSpin) {
	      std::vector<std::complex<double>> quark_sink(nEigs);
              for (int iEv = 0; iEv < nEigs; ++iEv) {
                quark_sink[iEv] = soln_rescale * qudaRes(iSpin, t, iEv, iSink);
              }
              DHputPtr->putData(RecordKey(iSpin + 1, t, sinkBatchInds[iSink]),
                                quark_sink);
              if (verbose) {
                printLaph(make_strf("srcev_index = %d, spin = %d, time = %d",
                                    sinkBatchInds[iSink], iSpin + 1, t));
                for (int n = 0; n < nEigs; n++) {
                  printLaph(
                      make_strf("coef for eigenlevel %d = (%14.8f, %14.8f)", n,
                                real(quark_sink[n]), imag(quark_sink[n])));
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
    }
  } // src_ind, src_spin loop end

  inv_time += invtime;
  makesrc_time += srctime;
  evproj_time += evprojtime;
  write_time += writetime;
}

// T is expected to be double or float
template <class T>
static void new_make_source(LattField &ferm_src,
			    const void *ev_src_ptr,
			    const int src_time,
			    const int src_spin,
			    void (*blasfunc)( const int N, const void *alpha, const void *x, const int incX,
					      void* y, const int incy) ) {
  printLaph(" Making source for this inversion...");
  assert(src_spin>0);
  ferm_src.reset(FieldSiteType::ColorSpinVector);

  // initialize source field to zero
  const int loc_nsites = LayoutInfo::getRankLatticeNumSites();
  std::complex<T> zrho = std::complex<T>( src_spin>2?-1:1 , 0 ) ;

  // fill it with zeros
  T *z0 = reinterpret_cast<T *>(ferm_src.getDataPtr());
  std::fill(z0, z0 + 2*loc_nsites*ferm_src.elemsPerSite(), 0.0);
  
  const int loc_npsites = loc_nsites/2;
  const int start_parity = LayoutInfo::getMyStartParity();
  const int mytmin = LayoutInfo::getMyCommCoords()[3] * LayoutInfo::getRankLattExtents()[3];
  const int mytmax = mytmin + LayoutInfo::getRankLattExtents()[3] - 1;
  const int tstride = LayoutInfo::getRankLattExtents()[0] *
                      LayoutInfo::getRankLattExtents()[1] *
                      LayoutInfo::getRankLattExtents()[2];
  
  // could be more efficient
  if ((src_time >= mytmin) && (src_time <= mytmax)) {
    const int tloc = src_time - mytmin;
    int parshift = loc_npsites * ((start_parity + tloc) % 2);
    const int start1 = ((tstride * tloc) / 2) + parshift;
    const int stop1 = ((1 + tstride * (tloc + 1)) / 2) + parshift;
    parshift = loc_npsites * ((start_parity + 1 + tloc) % 2);
    const int start2 = ((1 + tstride * tloc) / 2) + parshift;
    const int stop2 = ((tstride * (tloc + 1)) / 2) + parshift;

    const std::complex<T> *x0 = reinterpret_cast<const std::complex<T>*>(ev_src_ptr);
    const std::complex<T>*x1 = x0 + start1*FieldNcolor ;
    const std::complex<T>*x2 = x0 + start2*FieldNcolor ;

    std::complex<T> *y0 = reinterpret_cast<std::complex<T>*>(ferm_src.getDataPtr());
    std::complex<T>* y1 = y0 + FieldNcolor*((src_spin-1) + start1*FieldNspin) ;
    std::complex<T>* y2 = y0 + FieldNcolor*((src_spin-1) + start2*FieldNspin) ;
    
    for (int c = 0; c < FieldNcolor; ++c) {
      blasfunc( (stop1-start1), &zrho, x1+c, FieldNcolor, y1+c, FieldNcolor*FieldNspin ) ;
      blasfunc( (stop2-start2), &zrho, x2+c, FieldNcolor, y2+c, FieldNcolor*FieldNspin ) ;
    }
  }
  printLaph("Source for this inversion created");
}
  
// Creates the source (multiplied by gamma_4) in Dirac-Pauli basis
//   src_spin here is 1,2,3,4

void PerambulatorHandler::make_source(LattField &ferm_src,
                                      const void *ev_src_ptr,
                                      const int src_time,
				      const int src_spin) {
  const bool dp = (ferm_src.bytesPerWord() == sizeof(std::complex<double>));
  if( dp ) {
    new_make_source<double>(ferm_src, ev_src_ptr, src_time, src_spin, cblas_zaxpy ) ;
  } else {
    new_make_source<float>( ferm_src, ev_src_ptr, src_time, src_spin, cblas_caxpy ) ;
  }
}

//  static pointers (set to null in default constructor)

std::unique_ptr<QuarkSmearingHandler> PerambulatorHandler::qSmearHandler;
std::unique_ptr<GaugeConfigurationHandler> PerambulatorHandler::gaugeHandler;

int PerambulatorHandler::qSmearCounter = 0;
int PerambulatorHandler::gaugeCounter = 0;

} // namespace LaphEnv
