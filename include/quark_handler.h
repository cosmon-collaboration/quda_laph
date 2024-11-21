#ifndef QUARK_HANDLER_H
#define QUARK_HANDLER_H

#include "data_io_handler.h"
#include "dilution_handler.h"
#include "dilution_scheme_info.h"
#include "dir_path.h"
#include "field_smearing_info.h"
#include "filelist_info.h"
#include "gauge_configuration_handler.h"
#include "gauge_configuration_info.h"
#include "gluon_smearing_handler.h"
#include "inverter_info.h"
#include "laph_noise.h"
#include "laph_noise_info.h"
#include "laph_stdio.h"
#include "quark_action_info.h"
#include "quark_smearing_handler.h"
#include "quda.h"
#include "xml_handler.h"
#include <list>

namespace LaphEnv {

// *****************************************************************
// *                                                               *
// *  "QuarkHandler" handles computation of and subsequent access  *
// *  to the quark source and sink functions.                      *
// *                                                               *
// *  One of these handlers deals with quark sources/sinks         *
// *  for **one** set of info parameters, given by                 *
// *                                                               *
// *         GaugeConfigurationInfo                                *
// *         GluonSmearingInfo                                     *
// *         QuarkSmearingInfo                                     *
// *         DilutionSchemeInfo                                    *
// *         QuarkActionInfo                                       *
// *         FileListInfo                                          *
// *                                                               *
// *  but one handler deals with different                         *
// *                                                               *
// *         time dilution projectors                              *
// *         laph noises                                           *
// *         all spin-eigvec dilution indices                      *
// *                                                               *
// *  File structure and contents:                                 *
// *                                                               *
// *   - Results manipulated by one handler are contained in       *
// *     several files.  Each file has the same stub, but          *
// *     different positive integer suffices.                      *
// *             stub.0                                            *
// *             stub.1                                            *
// *             ...                                               *
// *             stub.N                                            *
// *     The files included are specified in a FileListInfo.       *
// *                                                               *
// *   - The header info in each file has a common part and a      *
// *     part that is specific to that one file:                   *
// *        <QuarkHandlerDataFile>                                 *
// *           - common part                                       *
// *           - specific part                                     *
// *        </QuarkHandlerDataFile>                                *
// *                                                               *
// *   - The common header info includes                           *
// *         GaugeConfigurationInfo                                *
// *         GluonSmearingInfo                                     *
// *         QuarkSmearingInfo                                     *
// *         DilutionSchemeInfo                                    *
// *         QuarkActionInfo                                       *
// *                                                               *
// *   - The specific header info includes (FileKey)               *
// *         LaphNoiseInfo noise;                                  *
// *         int time_projector_index;                             *
// *                                                               *
// *   - Each file contains several records whose access is given  *
// *     by a RecordKey.  The RecordKey contains spin, time, and   *
// *     the spin-eigvec dilution index.  The data in each record  *
// *     (DataType) is a multi1d<Complex> containing "nev"         *
// *     complex numbers, where "nev" is the number of Laph        *
// *     eigenvectors.  Hence, each quark sink is stored by its    *
// *     superposition coefficients in terms of the Laph           *
// *     eigenvectors.   Only sinks are stored.  Sources are made  *
// *     on the fly.  The Dirac-Pauli spin convention is used.     *
// *     Note that these coefficients are **gauge invariant.**     *
// *                                                               *
// *   - In 3d, the "get" functions must re-construct the full     *
// *     LatticeColorVector by evaluating the superposition of     *
// *     Laph eigenvectors.  Hence, the Laph eigenvectors are      *
// *     needed at this stage too. Quark displacements are done    *
// *     on the fly as well, so the smeared gauge field is needed. *
// *                                                               *
// *  All Laph Handlers follow the member naming convention:       *
// *                                                               *
// *    compute....()  to do original computation                  *
// *    set...()       to internally set from file or NamedObjMap  *
// *                                                               *
// *    get...()       provides access to results                  *
// *                                                               *
// *  When using "set" and "get", individual spin components and   *
// *  time slices are returned, and this handler computes the      *
// *  requested covariant displacements on the fly.                *
// *                                                               *
// *  The usual use of a QuarkHandler is as follows:               *
// *                                                               *
// *   - to compute the quark sources/sinks:                       *
// *                                                               *
// *       QuarkHandler Q;  // declare                             *
// *       Q.setInfo(...);            // input common info         *
// *       Q.setInverter(...);        // input inverter info       *
// *       Q.computeSink(...);                                     *
// *                                                               *
// *   - to subsequently use results to compute                    *
// *     hadron source/sinks:                                      *
// *                                                               *
// *       QuarkHandler Q;                                         *
// *       Q.setInfo(...);                                         *
// *       Q.getSink(...);      // stores in internal map          *
// *       Q.getSources(...);   // also does displacements         *
// *                                                               *
// *       Q.querySources(...);  // query if available             *
// *       Q.querySink(...);                                       *
// *                                                               *
// *       Q.removeSources(...);  // remove from internal map      *
// *       Q.removeSink(...);                                      *
// *       Q.clearOnlyDisplacedData();                             *
// *       Q.clearData();                                          *
// *                                                               *
// *  Gamma-5 Hermiticity can be applied by calling                *
// *                                                               *
// *       Q.setGamma5HermiticityMode();                           *
// *                                                               *
// *  This essentially multiplies sinks by gamma_5*gamma_4 and     *
// *  sources by -gamma_5*gamma_4.                                 *
// *                                                               *
// *****************************************************************

class QuarkHandler {

public:
  // Spin is 1,2,3,4  (but stored as 0,1,2,3 in right-most two bits)
  // time (next 15 bits)
  // spinlev_index (next 15 bits)

  class RecordKey {
    unsigned long code; // at least 32 bit integer

  public:
    RecordKey() : code(0) {}
    RecordKey(int spin, int time, int splev_ind) {
      encode(spin, time, splev_ind);
    }
    RecordKey(const RecordKey &in) : code(in.code) {}
    RecordKey &operator=(const RecordKey &in) {
      code = in.code;
      return *this;
    }
    RecordKey &set(int spin, int time, int splev_ind) {
      encode(spin, time, splev_ind);
      return *this;
    }
    ~RecordKey() {}

    bool operator<(const RecordKey &rhs) const { return (code < rhs.code); }
    bool operator==(const RecordKey &rhs) const { return (code == rhs.code); }
    bool operator!=(const RecordKey &rhs) const { return (code != rhs.code); }

    unsigned int getSpin() const {
      const unsigned long spinmask = 0x3ul;
      return (code & spinmask) + 1;
    }
    unsigned int getTime() const {
      const unsigned long mask15bit = 0x7FFFul;
      return (code >> 2) & mask15bit;
    }
    unsigned int getSpinLaphEigvecIndex() const { return (code >> 17); }

    void output(XMLHandler &xmlw) const;

    explicit RecordKey(const unsigned int *buf) { code = *buf; }
    static int numints() { return 1; }
    size_t numbytes() const { return sizeof(unsigned int); }
    void copyTo(unsigned int *buf) const { *buf = code; }

    void applyGamma5Herm(int &g5sign) {
      const unsigned long spinmask = 0x3ul;
      g5sign = ((code & spinmask) < 2) ? -1 : 1;
      code ^= 0x2ul;
    } // flip the bit second from right

  private:
    void encode(int spin, int time, int splev_ind) {
      const int imax = 32768;
      if ((spin < 1) || (spin > 4) || (time < 0) || (splev_ind < 0)) {
        errorLaph("invalid indices in QuarkHandler::RecordKey");
      }
      if ((time >= imax) || (splev_ind >= imax)) {
        errorLaph("indices in QuarkHandler::RecordKey exceed maximum");
      }
      code = splev_ind;
      code <<= 15;
      code |= time;
      code <<= 2;
      code |= spin - 1;
    }
  };

  struct FileKey {
    LaphNoiseInfo noise;
    int time_proj_index;

    FileKey() : time_proj_index(0) {}
    FileKey(const LaphNoiseInfo &in_noise, int tprojind);
    FileKey(XMLHandler &xmlr);
    FileKey(const FileKey &rhs);
    FileKey &operator=(const FileKey &rhs);
    ~FileKey() {}
    void output(XMLHandler &xmlw) const;
    bool operator<(const FileKey &rhs) const;
    bool operator==(const FileKey &rhs) const;
    bool operator!=(const FileKey &rhs) const;
  };

  typedef std::vector<std::complex<double>>
      DataType; // store in double precision even if single precision

private:
  struct StorageKey // used for mapping in set,get routines
  {
    FileKey fkey;
    RecordKey rkey;
    DirPath disp;
    int disp_length;

    StorageKey(const StorageKey &in);
    StorageKey(const FileKey &in_fkey, const RecordKey &in_rkey,
               const DirPath &in_disp, int in_disp_length);
    StorageKey(const FileKey &in_fkey, const RecordKey &in_rkey);
    StorageKey &operator=(const StorageKey &in);
    ~StorageKey() {}

    bool operator<(const StorageKey &rhs) const;
    bool operator==(const StorageKey &rhs) const;
    bool operator!=(const StorageKey &rhs) const;
  };

  struct SinkComputation {

    LaphNoiseInfo Noise;
    int TimeProjIndex;

    SinkComputation(const LaphNoiseInfo &in_noise, int in_time_proj_index)
        : Noise(in_noise), TimeProjIndex(in_time_proj_index) {}
  };

  struct SinkComputations {

    std::list<SinkComputation> computations;
    uint nSinkLaphBatch; // number of inversions before projecting on Laph evs
    uint nSinkQudaBatch; // number of quark sinks to project at a time as one
                         // batch
    uint nEigQudaBatch;  // number of Laph evs to project at a time as one batch
  };

  // pointers to internal infos (managed by this handler
  // with new and delete)

  const GaugeConfigurationInfo *uPtr;
  const GluonSmearingInfo *gSmearPtr;
  const QuarkSmearingInfo *qSmearPtr;
  const DilutionSchemeInfo *dilPtr;
  const QuarkActionInfo *qactionPtr;
  const FileListInfo *fPtr;
  const InverterInfo *invertPtr;
  bool compute_mode;
  void *preconditioner;

  // sub-handler pointers

  static std::unique_ptr<GluonSmearingHandler> gSmearHandler;
  static std::unique_ptr<QuarkSmearingHandler> qSmearHandler;
  static std::unique_ptr<GaugeConfigurationHandler> gaugeHandler;
  mutable DilutionHandler *dilHandler;

  static int gSmearCounter;
  static int qSmearCounter;
  static int gaugeCounter;
  // static bool keepInMemory;

  // structure containing the computations to perform
  SinkComputations sinkComps;

  // necessary quda data
  QudaInvertParam quda_inv_param;

  // Prevent copying ... handler might contain large
  // amounts of data

  QuarkHandler(const QuarkHandler &);
  QuarkHandler &operator=(const QuarkHandler &);

  // data I/O handler pointers

  DataPutHandlerMF<QuarkHandler, FileKey, RecordKey, DataType> *DHputPtr;
  DataGetHandlerMF<QuarkHandler, FileKey, RecordKey, DataType> *DHgetPtr;

  //  internal storage for sources/sinks (used in
  //  hadron construction)
  //  (NOT declared static to make clearing data easier when
  //  using multiple handler; but take care with memory footprint)

  std::map<StorageKey, LattField *> store_sinks;
  std::map<StorageKey, LattField *> store_sources;
  bool normal_mode;

public:
  QuarkHandler();

  QuarkHandler(const GaugeConfigurationInfo &gaugeinfo,
               const GluonSmearingInfo &gluonsmear,
               const QuarkSmearingInfo &quarksmear,
               const DilutionSchemeInfo &dil, const QuarkActionInfo &quark,
               const FileListInfo &flist,
               const std::string &smeared_quark_filestub,
               const std::string &smeared_gauge_filename,
               bool setComputeMode = true);

  void setInfo(const GaugeConfigurationInfo &gaugeinfo,
               const GluonSmearingInfo &gluonsmear,
               const QuarkSmearingInfo &quarksmear,
               const DilutionSchemeInfo &dil, const QuarkActionInfo &quark,
               const FileListInfo &flist,
               const std::string &smeared_quark_filestub,
               const std::string &smeared_gauge_filename,
               bool setComputeMode = true);

  ~QuarkHandler();

  void clear();

  // other set info members

  void setInverter(const InverterInfo &invinfo);

  // static void setSubHandlersKeepInMemory() {keepInMemory=true;}

  // static void unsetSubHandlersKeepInMemory() {keepInMemory=false;}

  void clearSinkComputations();

  void setSinkComputations(const XMLHandler &xmlcmp);

  void setNormalMode();

  void setGamma5HermiticityMode();

  bool isInfoSet() const;

  bool isComputeReady() const;

  const GaugeConfigurationInfo &getGaugeConfigurationInfo() const;

  const GluonSmearingInfo &getGluonSmearingInfo() const;

  const QuarkSmearingInfo &getQuarkSmearingInfo() const;

  const DilutionSchemeInfo &getDilutionSchemeInfo() const;

  const QuarkActionInfo &getQuarkActionInfo() const;

  const FileListInfo &getFileListInfo() const;

  int getNumberOfSpinEigvecDilutionProjectors() const;

  int getNumberOfTimeDilutionProjectors() const;

  int getTimeDilutionProjectorIndex(int time_val) const;

  const std::list<int> &getOnTimes(int time_proj_index) const;

  bool isFullTimeDilution() const;

  int getTimeExtent() const;

  //   void getHeader(XMLHandler& xmlout) const;

  //   void getFileMap(XMLHandler& xmlout) const;

  void outputSuffixMap();

  // void outputSuffixMap(TextFileWriter& fout);

  const InverterInfo &getInverterInfo() const;

  void setUpPreconditioning(QudaInvertParam &invParam);

  // if "verbose" is set to true, the quark sinks
  // will be output to standard output too

  void computeSinks(bool verbose = false, bool extra_soln_check = false);

  bool isNormalMode() const { return normal_mode; }

  bool isGamma5HermiticityMode() const { return !normal_mode; }

  //   void outputKeys(XMLHandler& xmlout);

  //   std::set<FileKey> getNoisesAndTimeProjectors() const;

  //   void getNoisesAndTimeProjectors(std::set<LaphNoiseInfo>& noises,
  //                                   std::set<int>& time_proj_indices) const;

  //   void getNoisesAndSourceTimes(std::set<LaphNoiseInfo>& noises,
  //                                std::set<int>& source_times) const;

  std::map<int, FileKey> getSuffixMap() const;

  // Access to the displaced sources/sinks.  If a source is zero,
  // a null pointer is returned.
  /*
     const LatticeColorVector* getData(bool source, const LaphNoiseInfo& noise,
                                       int time_proj_index,
                                       int spinlev_dilution_index, int spin,
                                       const DirPath& displace, int disp_length,
                                       int time);

     multi1d<Complex> getCoefficients(bool source, const LaphNoiseInfo& noise,
                                      int time_proj_index,
                                      int spinlev_dilution_index, int spin,
                                      int time);

     bool queryData(bool source, const LaphNoiseInfo& noise, int
     time_proj_index, int spinlev_dilution_index, int spin=1, int time=-1);


     void removeData(bool source, const LaphNoiseInfo& noise,
                     int time_proj_index, int spinlev_dilution_index, int spin,
                     const DirPath& displace, int disp_length, int time);


     void clearOnlyDisplacedData();
  */
  void clearData();
  /*
     void clearGaugeData();

     DataType getLaphEigenvectorSinkCoefficients(const LaphNoiseInfo& noise,
                                                 int time_proj_index,
                                                 int spinlev_dilution_index,
                                                 int spin, int time);

     QuarkSmearingHandler& getQuarkSmearingHandler()
      {return *qSmearHandler;}

     GluonSmearingHandler& getGluonSmearingHandler()
      {return *gSmearHandler;}
  */

private:
  void set_info(const GaugeConfigurationInfo &gaugeinfo,
                const GluonSmearingInfo &gluonsmear,
                const QuarkSmearingInfo &quarksmear,
                const DilutionSchemeInfo &dil, const QuarkActionInfo &quark,
                const FileListInfo &flist,
                const std::string &smeared_quark_filestub,
                const std::string &smeared_gauge_filename);

  bool checkHeader(XMLHandler &xmlr, int suffix);
  void writeHeader(XMLHandler &xmlout, const FileKey &fkey, int suffix);

  void check_info_set(const std::string &name) const;

  void check_compute_ready(const std::string &name) const;

  //  sub-handler connections

  void connectGluonSmearingHandler(const std::string &smeared_gauge_filename);
  void connectGaugeConfigurationHandler();
  void connectQuarkSmearingHandler(const std::string &smeared_quark_filestub);
  void connectDilutionHandler() const;

  void disconnectGluonSmearingHandler();
  void disconnectGaugeConfigurationHandler();
  void disconnectQuarkSmearingHandler();
  void disconnectDilutionHandler() const;

  void computeSinks(const LaphNoiseInfo &noise, int time_proj_index,
                    const std::vector<void *> &evList, bool verbose,
                    bool extra_soln_check, double &srctime, double &invtime,
                    double &evprojtime, double &writetime);

  // Makes the source in the Dirac-Pauli basis.  The source is
  // gamma_4 times the usual source since we use the chi = psi-bar gamma_4
  // field operator.
  void make_source(LattField &ferm_src, const Array<cmplx> &laph_noise,
                   const std::vector<void *> &evList,
                   const std::list<int> &on_times,
                   const std::list<int> &on_spins,
                   const std::list<int> &on_eigs);

  friend class DataPutHandlerMF<QuarkHandler, FileKey, RecordKey, DataType>;
  friend class DataGetHandlerMF<QuarkHandler, FileKey, RecordKey, DataType>;
};

typedef QuarkHandler::FileKey NoiseAndTimeProjector;

// ***************************************************************
} // namespace LaphEnv
#endif
