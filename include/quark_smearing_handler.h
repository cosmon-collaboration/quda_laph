#ifndef QUARK_SMEARING_HANDLER_H
#define QUARK_SMEARING_HANDLER_H

#include "array.h"
#include "data_io_handler.h"
#include "field_smearing_info.h"
#include "gauge_configuration_info.h"
#include "laph_eigen_info.h"
#include "xml_handler.h"

namespace LaphEnv {

// *****************************************************************
// *                                                               *
// *  QuarkSmearingHandler:                                        *
// *                                                               *
// *     -- must be constructed in read mode or write mode         *
// *     -- in 3-d write mode, computes the Laplacian eigenvectors,*
// *          saving results to "*_time.nn" files, each file       *
// *          containing all eigenvectors for one time slice OR    *
// *          inserting into TheNamedObjMap                        *
// *     -- in 3-d read mode, reads the eigenvectors from the      *
// *          "*_time.nn" files and makes them available for the   *
// *          hadron handler to do quark displacements OR gets     *
// *          from TheNamedObjMap                                  *
// *     -- in 4-d write mode, either reads the "*_time.nn" files  *
// *          and combines each level on all time slices into one  *
// *          4-d file, writing out "_level.nn" files, OR it       *
// *          computes them on all requested time slices and       *
// *          inserts results into TheNamedObjMap                  *
// *     -- in 4-d read mode, reads the "*_level.nn" files and     *
// *          makes them available for the quark handler to        *
// *          compute quark sinks, or reads from TheNamedObjMap    *
// *          (if "level" files are not available, it will read    *
// *           the "time" files)                                   *
// *                                                               *
// *  A note concerning the phases multiplying each eigenvector:   *
// *                                                               *
// *   Due to the way that noise is introduced, changing the       *
// *   overall phase of any given Laph eigenvector changes the     *
// *   value of the quark line for one particular noise.  The      *
// *   effect of the phase change is to change the noise (the      *
// *   noise is effectively U(1)).  This is not a problem, but     *
// *   erroneous results can occur if the original eigenvector     *
// *   files used to determine the quark sinks get deleted and     *
// *   the eigenvectors have to be reconstructed for making the    *
// *   hadrons.  With different run parameters, the eigensolver    *
// *   could produce a different phase.  The introduction of a     *
// *   phase convention eliminates this potential problem.         *
// *   The phase convention adopted here is as follows:            *
// *                                                               *
// *   *** each eigenvector is multiplied by a phase so that       *
// *       the 0-th color component at site (0,0,0)                *
// *       is real and positive (abort if component very small)    *
// *                                                               *
// *                                                               *
// *  All Laph Handlers follow the member naming convention:       *
// *                                                               *
// *    compute....()  to do original computation                  *
// *    set...()       to internally set from file or NamedObjMap  *
// *                                                               *
// *    get...()       provides access to results                  *
// *                                                               *
// *****************************************************************

class QuarkSmearingHandler {

  struct TimeKey {
    int value;

    TimeKey() : value(0) {}
    TimeKey(int in_val) : value(in_val) {}
    TimeKey(XMLHandler &xmlr) {
      xmlread(xmlr, "TimeSlice", value, "QuarkSmearingHandler::TimeKey");
    }
    TimeKey(const TimeKey &in) : value(in.value) {}
    ~TimeKey() {}

    bool operator<(const TimeKey &rhs) const { return (value < rhs.value); }

    bool operator==(const TimeKey &rhs) const { return (value == rhs.value); }

    void output(XMLHandler &xmlw) const {
      xmlw.set_root("TimeSlice", make_string(value));
    }

    explicit TimeKey(const unsigned int *buf) : value(*buf) {}
    static int numints() { return 1; }
    size_t numbytes() const { return sizeof(unsigned int); }
    void copyTo(unsigned int *buf) const { *buf = value; }
  };

  struct LevelKey {
    int value;

    LevelKey() : value(0) {}
    LevelKey(int in_val) : value(in_val) {}
    LevelKey(XMLHandler &xmlr) {
      xmlread(xmlr, "Level", value, "QuarkSmearingHandler::LevelKey");
    }
    LevelKey(const LevelKey &in) : value(in.value) {}
    ~LevelKey() {}

    bool operator<(const LevelKey &rhs) const { return (value < rhs.value); }

    bool operator==(const LevelKey &rhs) const { return (value == rhs.value); }

    void output(XMLHandler &xmlw) const {
      xmlw.set_root("Level", make_string(value));
    }

    explicit LevelKey(const unsigned int *buf) : value(*buf) {}
    static int numints() { return 1; }
    size_t numbytes() const { return sizeof(unsigned int); }
    void copyTo(unsigned int *buf) const { *buf = value; }
  };

  const GaugeConfigurationInfo *uPtr;
  const GluonSmearingInfo *gSmearPtr;
  QuarkSmearingInfo *qSmearPtr;
  std::string smearedQuarkFileStub;
  bool m_read_mode;

  // storage and/or references to internal data, and other handlers
  Array<double> Eigenvalues;
  DataGetHandlerMFO<QuarkSmearingHandler, LevelKey, LevelKey, LattField>
      *dh_ptr;

  // prevent copying ... handler might contain large
  // amounts of data
  QuarkSmearingHandler(const QuarkSmearingHandler &);
  QuarkSmearingHandler &operator=(const QuarkSmearingHandler &);

public:
  QuarkSmearingHandler();

  QuarkSmearingHandler(const GluonSmearingInfo &gluon_smearing,
                       const GaugeConfigurationInfo &gauge,
                       const QuarkSmearingInfo &quark_smearing,
                       const std::string &smeared_quark_file_stub,
                       bool read_mode = true);

  void setInfo(const GluonSmearingInfo &gluon_smearing,
               const GaugeConfigurationInfo &gauge,
               const QuarkSmearingInfo &quark_smearing,
               const std::string &smeared_quark_file_stub,
               bool read_mode = true);

  // update if number of eigenvectors needs to be increased
  void updateSmearing(const QuarkSmearingInfo &quark_smearing);

  ~QuarkSmearingHandler();

  void clear(); // clears everything in the handler

  // access to the info

  bool isInfoSet() const;

  const GluonSmearingInfo &getGluonSmearingInfo() const;

  const QuarkSmearingInfo &getQuarkSmearingInfo() const;

  const GaugeConfigurationInfo &getGaugeConfigurationInfo() const;

  const std::string &getSmearedQuarkFieldFileStub() const;

  uint getNumberOfLaplacianEigenvectors() const {
    return qSmearPtr->getNumberOfLaplacianEigenvectors();
  }

  void computeLaphEigenvectors(const LaphEigenSolverInfo &solver_info,
                               const std::string &smeared_gauge_file);

  // get and query data when in read mode
  const LattField &getLaphEigenvector(int eigpair_num);

  bool queryLaphEigenvector(int eigpair_num);

  void removeLaphEigenvector(int eigpair_num);

  void clearLaphEigenvectors();

  void closeLaphLevelFiles();

  // checks to see if all _level files exist.
  bool checkAllLevelFilesExist();

  // public so I can test it
  void applyLaphPhaseConvention(std::vector<LattField> &laph_eigvecs);

private:
  void set_info(const GluonSmearingInfo &gluon_smearing,
                const GaugeConfigurationInfo &gauge,
                const QuarkSmearingInfo &quark_smearing,
                const std::string &smearedQuarkFileStub, bool read_mode);

  void check_info_set(const std::string &name, int check_mode = 0) const;

  void failure(const std::string &message);

  bool checkHeader(XMLHandler &xmlr, int suffix);

  void
  checkLaphEigvecComputation(const std::vector<LattField> &laphEigvecs,
                             const std::vector<LattField> &smeared_gauge_field);

  void writeHeader(XMLHandler &xmlw, const LevelKey &fkey, int suffix);

  friend class DataGetHandlerMF<QuarkSmearingHandler, LevelKey, LevelKey,
                                LattField>;

  friend class DataPutHandlerMF<QuarkSmearingHandler, LevelKey, LevelKey,
                                LattField>;

  friend class DataGetHandlerMFNOM<QuarkSmearingHandler, LevelKey, LevelKey,
                                   LattField>;

  friend class DataPutHandlerMFNOM<QuarkSmearingHandler, LevelKey, LevelKey,
                                   LattField>;
};
} // namespace LaphEnv
#endif
