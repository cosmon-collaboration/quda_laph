#include "QudaLaphIncludes.h"

#include "quark_smearing_handler.h"
#include "field_ops.h"
#include "gluon_smearing_handler.h"

// STILL TO DO:  single asynchronous thread for output so can continue next
// computations while output occurs
using namespace std;

namespace LaphEnv {

// constructors
QuarkSmearingHandler::QuarkSmearingHandler()
    : uPtr(0), gSmearPtr(0), qSmearPtr(0), m_read_mode(true), dh_ptr(0) {}

QuarkSmearingHandler::QuarkSmearingHandler(
    const GluonSmearingInfo &gluon_smearing,
    const GaugeConfigurationInfo &gauge,
    const QuarkSmearingInfo &quark_smearing,
    const string &smeared_quark_file_stub, const bool read_mode) {
  set_info(gluon_smearing, gauge, quark_smearing, smeared_quark_file_stub,
           read_mode);
}

void QuarkSmearingHandler::setInfo(const GluonSmearingInfo &gluon_smearing,
                                   const GaugeConfigurationInfo &gauge,
                                   const QuarkSmearingInfo &quark_smearing,
                                   const string &smeared_quark_file_stub,
                                   const bool read_mode) {
  clear();
  set_info(gluon_smearing, gauge, quark_smearing, smeared_quark_file_stub,
           read_mode);
}

void QuarkSmearingHandler::set_info(const GluonSmearingInfo &gluon_smearing,
                                    const GaugeConfigurationInfo &gauge,
                                    const QuarkSmearingInfo &quark_smearing,
                                    const string &smeared_quark_file_stub,
                                    const bool read_mode) {
  smearedQuarkFileStub = tidyString(smeared_quark_file_stub);
  if (smearedQuarkFileStub.empty()) {
    errorLaph("empty file name or stub in QuarkSmearingHandler");
  }

  try {
    uPtr = new GaugeConfigurationInfo(gauge);
    gSmearPtr = new GluonSmearingInfo(gluon_smearing);
    qSmearPtr = new QuarkSmearingInfo(quark_smearing);
  } catch (const std::exception &xp) {
    errorLaph("problem in setInfo in QuarkSmearingHandler");
  }

  m_read_mode = read_mode;

  if (m_read_mode) {
    // read from files (or the NamedObjMap)
    FileListInfo files(smearedQuarkFileStub + "_level", 0,
                       qSmearPtr->getNumberOfLaplacianEigenvectors() - 1,
                       false);
    try {
      dh_ptr = new DataGetHandlerMFO<QuarkSmearingHandler, LevelKey, LevelKey,
                                     LattField>(
          *this, files, "Laph--SmearedQuarkLevelFile", "LaphEigenvectors");
      if (dh_ptr->queryFile(LevelKey(0)))
        printLaph("QuarkSmearingHandler file stub appears to be okay");
      else {
        errorLaph(make_strf("problem with QuarkSmearingHandler file stub: %s",
                            smeared_quark_file_stub));
      }
    } catch (const std::exception &xp) {
      errorLaph("unable to allocate data handler in QuarkSmearingHandler");
    }
  } else {
    dh_ptr = 0;
  }
}

// use this to increase the number of Laplacian eigenvectors allowed;
// need inside QuarkHandler since smearing sub-handler is static.
void QuarkSmearingHandler::updateSmearing(
    const QuarkSmearingInfo &quark_smearing) {
  check_info_set("updateQuarkSmearingInfo");
  qSmearPtr->increaseUpdate(quark_smearing);
}

QuarkSmearingHandler::~QuarkSmearingHandler() { clear(); }

void QuarkSmearingHandler::clear() {
  try {
    delete gSmearPtr;
    delete qSmearPtr;
    delete uPtr;
    delete dh_ptr;
  } catch (const std::exception &xp) {
    errorLaph("abort");
  }
  smearedQuarkFileStub.clear();
  gSmearPtr = 0;
  qSmearPtr = 0;
  uPtr = 0;
  Eigenvalues.clear();
}

// access to the info
bool QuarkSmearingHandler::isInfoSet() const {
  return ((gSmearPtr != 0) && (uPtr != 0) && (qSmearPtr != 0) &&
          (!smearedQuarkFileStub.empty()));
}

// check_mode = 0 means no check, 1 means check for read mode,
// 2 means check for write mode
void QuarkSmearingHandler::check_info_set(const string &name,
                                          const int check_mode) const {
  if (!isInfoSet()) {
    errorLaph(make_strf(
        "error in QuarkSmearingHandler: must setInfo before calling %s", name));
  }
  if (check_mode == 1) {
    if (!m_read_mode) {
      errorLaph(make_strf(
          "error in QuarkSmearingHandler: must be in read mode when calling %s",
          name));
    }
  } else if (check_mode == 2) {
    if (m_read_mode) {
      errorLaph(make_strf("error in QuarkSmearingHandler: must not be in read "
                          "mode when calling %s",
                          name));
    }
  }
}

const GluonSmearingInfo &QuarkSmearingHandler::getGluonSmearingInfo() const {
  check_info_set("getGluonSmearingInfo");
  return *gSmearPtr;
}

const QuarkSmearingInfo &QuarkSmearingHandler::getQuarkSmearingInfo() const {
  check_info_set("getQuarkSmearingInfo");
  return *qSmearPtr;
}

const GaugeConfigurationInfo &
QuarkSmearingHandler::getGaugeConfigurationInfo() const {
  check_info_set("getGaugeConfigurationInfo");
  return *uPtr;
}

const string &QuarkSmearingHandler::getSmearedQuarkFieldFileStub() const {
  check_info_set("getSmearedQuarkFieldFileStub");
  return smearedQuarkFileStub;
}

void QuarkSmearingHandler::failure(const string &message) {
  errorLaph(make_strf("%s in QuarkSmearingHandler", message));
}

// Compute the Laph Eigenvectors simultaneously on all
// time slices (from mintime to maxtime in the *u_ptr).
// Put results into TheNamedObjMap or to level files.
void QuarkSmearingHandler::computeLaphEigenvectors(
    const LaphEigenSolverInfo &solver_info, const string &smeared_gauge_file) {
  check_info_set("computeLaphEigenvectors", 2);

  StopWatch rolex, bulova;
  rolex.start();
  double iotime = 0.0;

  const int nTime = uPtr->getTimeExtent();
  const int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();

  // check to see if output files exist already
  for (int v = 0; v < nEigvecs; ++v) {
    ostringstream oss;
    oss << smearedQuarkFileStub << "_level." << v;
    if (fileExists(oss.str())) {
      printLaph(make_str("file ", oss.str(), " already exists:\n",
                         " ... will not overwrite so computation aborted\n\n"));
      return;
    }
  }

  // fire up the GluonSmearingHandler
  GluonSmearingHandler gHandler(*gSmearPtr, *uPtr, smeared_gauge_file);

  printLaph(make_str("Computation of Laplacian eigenvectors commencing",
                     " in FieldSmearingHandler"));
  printLaph(make_str("  number of requested eigenvectors = ", nEigvecs));
  printLaph(make_str("            time extent of lattice = ", nTime));

  // fire up the data handler
  FileListInfo writeFiles(smearedQuarkFileStub + "_level", 0, nEigvecs - 1);
  DataPutHandlerMFO<QuarkSmearingHandler, LevelKey, LevelKey, LattField> DHput(
      *this, writeFiles, "Laph--SmearedQuarkLevelFile", "LaphEigenvectors",
      true, false);
  // read smeared gauge field
  StopWatch iotimer;
  iotimer.start();
  gHandler.copyDataToDevice();
  printLaph(make_str("\n\nSmeared gauge field obtained from ",
                     smeared_gauge_file, "\n"));

  iotimer.stop();
  iotime += iotimer.getTimeInSeconds();

  // use the Krylov-spectral restarted Lanczos method (blocked version) to
  // compute the eigenvalues/eigenvectors of minus the 3d Laplacian on all time
  // slices
  bulova.start();
  QudaInvertParam eig_inv_param;
  QudaEigParam eig_param;
  solver_info.setQudaParam(eig_inv_param, eig_param, *qSmearPtr);

  // allocate space for eigenvectors and eigenvalues on the host
  vector<LattField> laphEigvecs(nEigvecs, FieldSiteType::ColorVector);
  vector<complex<double>> laphEigvals(nEigvecs * nTime);
  __complex__ double *h_evals = (__complex__ double *)(laphEigvals.data());
  vector<void *> h_evecs(nEigvecs);
  for (int v = 0; v < nEigvecs; ++v) {
    h_evecs[v] = laphEigvecs[v].getDataPtr();
  }

  // set starting vector to be a constant one or set to zero for random
  if (solver_info.getStartingVectorType() == "equal_components") {
    complex<double> z(
        1.0 / double(FieldNcolor * LayoutInfo::getLatticeNumSites()), 0.0);
    setConstantField(laphEigvecs[0], z);
  } else {
    setZeroField(laphEigvecs[0]);
  }

  // now determine the eigenvectors!!
  eigensolveQuda(h_evecs.data(), h_evals, &eig_param);

  bulova.stop();
  printLaph(make_strf("\n\n SUCCESS: calculation time = %g seconds\n",
                      bulova.getTimeInSeconds()));
  if (solver_info.getChebyshevOrder() >= 2) {
    printLaph(make_str("  Chebyshev acceleration was used"));
  }

  //  multiply each eigenvector by a
  //  phase so that the zero-th color component at site (0,0,0)
  //  is real and positive; this uniquely specifies the overall
  //  phase for each eigenvector, which is important since a
  //  change in these phases can change the effective Laph noise.
  //  Using this convention makes the results robust against
  //  eigenvector deletion and subsequent reconstruction.
  bulova.reset();
  bulova.start();
  applyLaphPhaseConvention(laphEigvecs);
  bulova.stop();
  printLaph(make_strf(" Phase shift convention applied in %g seconds\n",
                      bulova.getTimeInSeconds()));

  if (solver_info.getOutputVerbosity() > 0) {
    printLaph("Eigenvalues of -Laplacian:");
    for (int t = 0; t < nTime; ++t)
      for (int v = 0; v < nEigvecs; ++v) {
        printLaph(make_strf(" eigenvalue[%3d](t=%3d) = %20.15f", v, t,
                            6.0 * real(laphEigvals[nEigvecs * t + v])));
      }
  }

  //  Do some checks of the results if requested
  if (solver_info.doCheckSolution()) {
    printLaph("Performing extra slower checks");
    checkLaphEigvecComputation(laphEigvecs, gHandler.getSmearedGaugeField());
  }

  iotimer.reset();
  iotimer.start();
  printLaph("Outputting results to file or named obj map");
  for (int k = 0; k < nEigvecs; k++) {
    DHput.open(LevelKey(k)); // open file, write header
    DHput.putData(LevelKey(k), laphEigvecs[k]);
    DHput.close(); // finalize current file
    printLaph(make_str("Level ", k, " with max ", nEigvecs - 1, " done"));
    laphEigvecs[k].clear();
  }

  gHandler.clearData(); // clear smeared gauge time slices
  iotimer.stop();
  printLaph(make_strf("Time to write laph eigenvectors = %g seconds",
                      iotimer.getTimeInSeconds()));
  iotime += iotimer.getTimeInSeconds();
  printLaph("Computation done \n");

  if (!QudaInfo::gauge_config_on_device) {
    freeGaugeQuda();
  } // undoes the hack above in copyDataToDevice

  rolex.stop();
  printLaph(make_strf("computeLaphEigenvectors: total time = %g seconds",
                      rolex.getTimeInSeconds()));
  printLaph(make_strf("Total file I/O time = %g seconds", iotime));
  printLaph("ran successfully\n");
}

// broadcast these rephase factors to the ranks that need them
template <class T> static void communicate_phase(vector<T> &rephase) {
#ifdef ARCH_PARALLEL
  vector<int> comm_coords = {
      0, 0, 0, LayoutInfo::getMyCommCoords()[LayoutInfo::Ndim - 1]};
  const int orig_sender_rank = LayoutInfo::getRankFromCommCoords(comm_coords);
  int sender_rank = 0, count;
  vector<int> broadcast_ranks;
  for (comm_coords[0] = 0;
       comm_coords[0] < LayoutInfo::getCommNumPartitions()[0];
       ++comm_coords[0]) {
    for (comm_coords[1] = 0;
         comm_coords[1] < LayoutInfo::getCommNumPartitions()[1];
         ++comm_coords[1]) {
      for (comm_coords[2] = 0;
           comm_coords[2] < LayoutInfo::getCommNumPartitions()[2];
           ++comm_coords[2], ++count) {
        const int next_rank = LayoutInfo::getRankFromCommCoords(comm_coords);
        if (next_rank == orig_sender_rank) {
          sender_rank = count;
        }
        broadcast_ranks.push_back(next_rank);
      }
    }
  }
  if (broadcast_ranks.size() > 1) {
    if (int(broadcast_ranks.size()) == LayoutInfo::getNumRanks()) {
      comm_broadcast(rephase.data(), rephase.size() * sizeof(rephase[0]),
                     orig_sender_rank);
    } else {
      MPI_Group world_group;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);
      MPI_Group time_group;
      MPI_Group_incl(world_group, broadcast_ranks.size(),
                     broadcast_ranks.data(), &time_group);
      MPI_Comm time_comm;
      MPI_Comm_create_group(MPI_COMM_WORLD, time_group, sender_rank,
                            &time_comm);
      if (time_comm != MPI_COMM_NULL) {
        const int status =
            MPI_Bcast(rephase.data(), rephase.size() * sizeof(rephase[0]),
                      MPI_BYTE, sender_rank, time_comm);
        if (status != MPI_SUCCESS) {
          errorLaph("Broadcast of re-phase factor in applyLaphPhaseConvention "
                    "failed");
        }
      }
      MPI_Group_free(&time_group);
      MPI_Group_free(&world_group);
      MPI_Comm_free(&time_comm);
    }
  }
#endif
}

template <class T>
static void getPhase(vector<T> &rephase, vector<LattField> &laph_evecs) {
  const int nev = laph_evecs.size();
  const int nloctime = LayoutInfo::getRankLattExtents()[3];
  const int loc_nsites = LayoutInfo::getRankLatticeNumSites();
  const int loc_npsites = loc_nsites / 2;
  const int start_parity = LayoutInfo::getMyStartParity();
  const int tstride = LayoutInfo::getRankLattExtents()[0] *
                      LayoutInfo::getRankLattExtents()[1] *
                      LayoutInfo::getRankLattExtents()[2];

  // rephase appears to be normalised V^{i}_{c=0}(0,t)
  if ((LayoutInfo::getMyCommCoords()[0] == 0) &&
      (LayoutInfo::getMyCommCoords()[1] == 0) &&
      (LayoutInfo::getMyCommCoords()[2] == 0)) {
    int v;
#pragma omp parallel for private(v) collapse(2)
    for (v = 0; v < nev; ++v) {
      for (int tloc = 0; tloc < nloctime; ++tloc) {
        const T *fp = reinterpret_cast<const T *>(laph_evecs[v].getDataPtr());
        const int parshift = loc_npsites * ((start_parity + tloc) % 2);
        const int start1 = (tstride * tloc) / 2 + parshift;
        const T x1 = fp[FieldNcolor * start1];
        rephase[tloc + v * nloctime] = conj(x1) / abs(x1);
      }
    }
  }
}

template <class T>
static void applyPhase(vector<LattField> &laph_evecs, const vector<T> rephase,
                       void (*blasfunc)(const int N, const void *alpha, void *x,
                                        const int incX)) {
  const int nev = laph_evecs.size();
  const int nloctime = LayoutInfo::getRankLattExtents()[3];
  const int loc_nsites = LayoutInfo::getRankLatticeNumSites();
  const int loc_npsites = loc_nsites / 2;
  const int start_parity = LayoutInfo::getMyStartParity();
  const int tstride = LayoutInfo::getRankLattExtents()[0] *
                      LayoutInfo::getRankLattExtents()[1] *
                      LayoutInfo::getRankLattExtents()[2];

  // now apply the rephase factors
  int v;
#pragma omp parallel for private(v) collapse(2)
  for (v = 0; v < nev; ++v) {
    for (int tloc = 0; tloc < nloctime; ++tloc) {
      T *fp = reinterpret_cast<T *>(laph_evecs[v].getDataPtr());
      for (int cb = 0; cb < 2; cb++) {
        const int parshift = loc_npsites * ((start_parity + cb + tloc) % 2);
        const int start = (((cb ^ 1) + tstride * tloc) / 2) + parshift;
        const int stop = (((cb ^ 0) + tstride * (tloc + 1)) / 2) + parshift;
        T *x = fp + FieldNcolor * start;
        blasfunc(FieldNcolor * (stop - start), &rephase[tloc + v * nloctime], x,
                 1);
      }
    }
  }
}

//  This is temporary; should be done on the device
void QuarkSmearingHandler::applyLaphPhaseConvention(
    vector<LattField> &laph_evecs) {
  const int nev = laph_evecs.size();
  if (nev == 0) {
    return;
  }
  for (int v = 0; v < nev; ++v) {
    if (laph_evecs[v].getFieldSiteType() != FieldSiteType::ColorVector) {
      errorLaph(
          "Applying phase convention can only be done to ColorVector fields");
    }
  }

  const bool dp =
      (laph_evecs[0].bytesPerWord() == sizeof(std::complex<double>));
  const int nloctime = LayoutInfo::getRankLattExtents()[3];

  if (dp) {
    vector<complex<double>> rephasedp(nev * nloctime);
    getPhase<complex<double>>(rephasedp, laph_evecs);

    communicate_phase<complex<double>>(rephasedp);

    applyPhase<complex<double>>(laph_evecs, rephasedp, cblas_zscal);

  } else {
    vector<complex<float>> rephasesp(nev * nloctime);
    getPhase<complex<float>>(rephasesp, laph_evecs);

    // (possibly) send it to the other nodes
    communicate_phase<complex<float>>(rephasesp);

    applyPhase<complex<float>>(laph_evecs, rephasesp, cblas_cscal);
  }

#ifdef ARCH_PARALLEL
  comm_barrier();
#endif
}

void QuarkSmearingHandler::checkLaphEigvecComputation(
    const vector<LattField> &laphEigvecs,
    const vector<LattField> &smeared_gauge_field) {
  int nTime = uPtr->getTimeExtent();
  int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();
  Array<double> block_eigvals(nEigvecs, nTime);
  vector<double> offmaxmag(nTime, 0.0);
  LattField minus_delta_on_ev(FieldSiteType::ColorVector);
  for (int v = 0; v < nEigvecs; ++v) {
    // act with the spatial Laplacian on this eigenvector  M=-Delta  M x =
    // lambda x
    applyMinusSpatialLaplacian(minus_delta_on_ev, laphEigvecs[v],
                               smeared_gauge_field);
    vector<complex<double>> eigvals(
        getTimeSlicedInnerProducts(laphEigvecs[v], minus_delta_on_ev));
    vector<complex<double>> norm(
        getTimeSlicedInnerProducts(laphEigvecs[v], laphEigvecs[v]));
    for (int t = 0; t < nTime; ++t) {
      block_eigvals(v, t) = real(eigvals[t]) / real(norm[t]);
    }
    for (int vv = 0; vv < nEigvecs; ++vv) {
      if (vv != v) {
        vector<complex<double>> offdiag(
            getTimeSlicedInnerProducts(laphEigvecs[vv], minus_delta_on_ev));
        for (int t = 0; t < nTime; ++t) {
          double offdiagmag = std::abs(offdiag[t]);
          if (offdiagmag > offmaxmag[t]) {
            offmaxmag[t] = offdiagmag;
          }
        }
      }
    }
  }
  for (int t = 0; t < nTime; ++t) {
    printLaph(make_str("\nCheck of Laph eigenvectors for time = ", t));
    printLaph(make_strf("  Maximum magnitude of offdiagonal <eigvec[v] | "
                        "-Delta | eigvec[vv]> = %g",
                        offmaxmag[t]));
    for (int v = 0; v < nEigvecs; ++v) {
      const bool orderok =
          (v == 0) || (block_eigvals(v, t) >= block_eigvals(v - 1, t));
      const string ok = (orderok) ? "" : " Wrong order";
      printLaph(make_strf("   Eigenvalue[%3d] = %20.15f %s", v,
                          block_eigvals(v, t), ok));
    }
  }
}

const LattField &
QuarkSmearingHandler::getLaphEigenvector(const int eigpair_num) {
  check_info_set("getLaphEigenvector", 1);
  LevelKey key(eigpair_num);
  return dh_ptr->getData(key, key);
}

bool QuarkSmearingHandler::queryLaphEigenvector(const int eigpair_num) {
  check_info_set("queryLaphEigenvector", 1);
  LevelKey key(eigpair_num);
  return dh_ptr->queryData(key, key);
}

void QuarkSmearingHandler::removeLaphEigenvector(int eigpair_num) {
  if (!m_read_mode)
    return;
  LevelKey key(eigpair_num);
  return dh_ptr->removeData(key, key);
}

void QuarkSmearingHandler::clearLaphEigenvectors() {
  if (!m_read_mode)
    return;
  dh_ptr->clearData();
  Eigenvalues.resize(0);
}

void QuarkSmearingHandler::closeLaphLevelFiles() {
  if (!m_read_mode)
    return;
  dh_ptr->closeAll();
}

bool QuarkSmearingHandler::checkAllLevelFilesExist() {
  check_info_set("checkAllLevelFilesExist");

  for (int level = 0; level < qSmearPtr->getNumberOfLaplacianEigenvectors();
       level++) {
    ostringstream oss;
    oss << smearedQuarkFileStub << "_level." << level;
    if (!fileExists(oss.str())) {
      printLaph(make_strf("needed file %s does not exists", oss.str()));
      return false;
    }
  }
  return true;
}

void QuarkSmearingHandler::writeHeader(XMLHandler &xmlw, const LevelKey &fkey,
                                       const int suffix) {
  if (suffix != fkey.value) {
    errorLaph("level suffix does not match file key");
  }
  xmlw.set_root("LaphEigenvectors");
  XMLHandler xmltmp;
  gSmearPtr->output(xmltmp);
  xmlw.put_child(xmltmp);
  uPtr->output(xmltmp);
  xmlw.put_child(xmltmp);
  qSmearPtr->output(xmltmp);
  xmlw.put_child(xmltmp);
  fkey.output(xmltmp);
  xmlw.put_child(xmltmp);
}

bool QuarkSmearingHandler::checkHeader(XMLHandler &xml_in, const int suffix) {
  if (xml_tag_count(xml_in, "LaphEigenvectors") != 1)
    return false;
  XMLHandler xmlr(xml_in, "LaphEigenvectors");
  try {
    // GaugeConfigurationInfo gauge_check(xmlr);   //  TO MATCH USQCD CHANGE
    // LATER!! uPtr->checkEqual(gauge_check);
    GluonSmearingInfo gsmear_check(xmlr);
    gSmearPtr->checkEqual(gsmear_check);
    QuarkSmearingInfo qsmear_check(xmlr);
    if (m_read_mode) {
      qSmearPtr->checkEqual(qsmear_check);
      LevelKey level_check(xmlr);
      if (level_check.value != suffix) {
        throw(std::invalid_argument("error"));
      }
    }
  } catch (const std::exception &xp) {
    return false;
  }
  return true;
}
} // namespace LaphEnv
