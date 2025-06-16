#ifndef PERAMBULATOR_HANDLER_H
#define PERAMBULATOR_HANDLER_H

#include "gauge_configuration_handler.h"
#include "quark_smearing_handler.h"
#include "inverter_info.h"
#include "sparse_grid.h"
#include "array.h"


namespace LaphEnv {


// *****************************************************************
// *                                                               *
// *  "PerambulatorHandler" handles computation of and subsequent  *
// *  access to the quark perambulators.                           *
// *                                                               *
// *  This handler operators in one of four modes:                 *
// *      enum Mode { ReadOnly, Merge, Compute, Check };           *
// *                                                               *
// *  Computations are done in "Compute" mode.  If computations    *
// *  for all source eigvec indices cannot be done in a single     *
// *  run and are sent to separate files, then these results       *
// *  can be merged when in "Merge" mode.  In "Check" mode,        *
// *  the results in the files are checked that they are all       *
// *  available and are not NaNs.  "ReadOnly" mode is used for     *
// *  accessing the results when computing correlators, for        *
// *  example.                                                     *
// *                                                               *
// *  One of these handlers deals with quark perambulators         *
// *  for **one** set of info parameters, given by                 *
// *                                                               *
// *         GaugeConfigurationInfo                                *
// *         GluonSmearingInfo                                     *
// *         QuarkSmearingInfo                                     *
// *         Nspin  (4 or 2)                                       *
// *         QuarkActionInfo                                       *
// *         FileListInfo                                          *
// *                                                               *
// *  but one handler deals with different                         *
// *                                                               *
// *         src_time, snk_time, spin and Laph eigvec indices      *
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
// *        <PerambulatorHandlerDataFile>                          *
// *           - common part                                       *
// *           - specific part                                     *
// *        </PerambulatorHandlerDataFile>                         *
// *                                                               *
// *   - The common header info includes                           *
// *         GaugeConfigurationInfo                                *
// *         GluonSmearingInfo                                     *
// *         QuarkSmearingInfo                                     *
// *         QuarkActionInfo                                       *
// *         Nspin  (4 or 2)                                       *
// *                                                               *
// *   - The specific header info includes (FileKey)               *
// *         int src_time, src_spin                                *
// *                                                               *
// *   - Each file contains several records whose access is given  *
// *     by a RecordKey.  The RecordKey contains an integer        *
// *     that specifies sink spin, sink time, and source eigenvec  *
// *     index.  The data in each record (DataType) is a           *
// *     vector<complex<double>> containing "nev" complex numbers, *
// *     where "nev" is the number of Laph eigenvectors.   The     *
// *     Dirac-Pauli spin convention is used.                      *
// *                                                               *
// *  All Laph Handlers follow the member naming convention:       *
// *                                                               *
// *    compute....()  to do original computation                  *
// *                                                               *
// *    get...()       provides access to results                  *
// *                                                               *
// *  When using "set" and "get", individual spin components and   *
// *  time slices are returned.  At this time, no quark            *
// *  covariant displacements are accommodated.                    *
// *                                                               *
// *  The usual use of a PerambulatorHandler is as follows:        *
// *                                                               *
// *   - to compute the quark sources/sinks:                       *
// *                                                               *
// *       PerambulatorHandler Q;  // declare                      *
// *       Q.setInfo(...);            // input common info         *
// *       Q.setInverter(...);        // input inverter info       *
// *       Q.computePerambulators(...);                            *
// *                                                               *
// *   - to subsequently use results to compute                    *
// *     hadron source/sinks:                                      *
// *                                                               *
// *       PerambulatorHandler Q;                                  *
// *       Q.setInfo(...);                                         *
// *       Q.getData(...);                                         *
// *                                                               *
// *  Usually all of the perambulators for all source ev indices   *
// *  cannot be done in a single run, and so for safety, the       *
// *  data for different sets of source ev indices will get saved  *
// *  into separate files.  These then need to be combined into    *
// *  a single file for subsequent reading to compute the pion     *
// *  and nucleon correlator.  The method                          *
// *                                                               *
// *     mergeData(...)                                            *
// *                                                               *
// *  is provided to accomplish this.                              *
// *                                                               *
// *****************************************************************
class SparseGridHandler { 
	const PerambulatorHandler &pHand;
	RandomSparseGrid &grid; 

	SparseGridHandler(const PerambulatorHandler& _pHand, 
			RandomSparseGrid& _grid) : pHand(_pHand), grid(_grid) {} 

	bool checkHeader(XMLHandler& xmlr, int suffix);
	void writeHeader(XMLHandler& xmlout, const FileKey& fkey,
                    int suffix);
}; 

class PerambulatorHandler {

 public:

   enum Mode { ReadOnly, Merge, Compute, Check };


       // sink spin is 1,2,3,4  (but stored as 0,1,2,3 in right-most two bits)
       // sink time (next 15 bits)
       // source_laphev_index (next 15 bits)

   class RecordKey {

      unsigned long code;   // at least 32 bit integer

    public:

      RecordKey() : code(0) {}
      RecordKey(int snkspin, int snktime, int srcev_index)
       {encode(snkspin,snktime,srcev_index);}
      RecordKey(const RecordKey& in) : code(in.code) {}
      RecordKey& operator=(const RecordKey& in) {code=in.code; return *this;}
      RecordKey& set(int snkspin, int snktime, int srcev_index)
       {encode(snkspin,snktime,srcev_index); return *this;}
      ~RecordKey() {}

      bool operator<(const RecordKey& rhs) const {return (code<rhs.code);}
      bool operator==(const RecordKey& rhs) const {return (code==rhs.code);}
      bool operator!=(const RecordKey& rhs) const {return (code!=rhs.code);}

      unsigned int getSinkSpin() const
       {const unsigned long spinmask=0x3ul; return (code&spinmask)+1;}
      unsigned int getSinkTime() const 
       {const unsigned long mask15bit=0x7FFFul; return (code>>2)&mask15bit;}
      unsigned int getSourceLaphEigvecIndex() const {return (code>>17);}

      void output(XMLHandler& xmlw) const;
      std::string output(int indent=0) const;
      
      explicit RecordKey(const unsigned int* buf) {code=*buf;}
      static int numints() {return 1;}
      size_t numbytes() const {return sizeof(unsigned int);}
      void copyTo(unsigned int* buf) const {*buf=code;}

    //  void applyGamma5Herm(int& g5sign)
    //  {const unsigned long spinmask=0x3ul;
    //   g5sign=((code&spinmask)<2)?-1:1;
    //   code^=0x2ul;}   // flip the bit second from right

    private:

      void encode(int snkspin, int snktime, int srcev_index)
      {
       const int imax=32768;
       if  ((snkspin<1)||(snkspin>4)||(snktime<0)||(srcev_index<0)){
          errorLaph("invalid indices in QuarkHandler::RecordKey");}
       if  ((snktime>=imax)||(srcev_index>=imax)){
          errorLaph("indices in QuarkHandler::RecordKey exceed maximum");}
       code=srcev_index; code<<=15; 
       code|=snktime; code<<=2; code|=snkspin-1;
      }

   };

   struct FileKey
   {
      int src_time;
      int src_spin;

      FileKey() : src_time(0), src_spin(1)  {}
      FileKey(int in_srctime, int in_srcspin);
      FileKey(XMLHandler& xmlr);
      FileKey(const FileKey& rhs);
      FileKey& operator=(const FileKey& rhs);
      ~FileKey(){}
      virtual void output(XMLHandler& xmlw) const;
      virtual std::string output(int indent=0) const;
      bool operator<(const FileKey& rhs) const;
      bool operator==(const FileKey& rhs) const;
      bool operator!=(const FileKey& rhs) const;
   };

	 struct PerambComputation {
      int src_time;
      std::set<int> src_lapheigvec_indices;
      PerambComputation(int in_src_time, const std::set<int>& src_evinds)
               : src_time(in_src_time), src_lapheigvec_indices(src_evinds) {}
   };

   struct PerambComputations {
      std::list<PerambComputation> computations;
      uint nSinkLaphBatch;      // number of inversions before projecting on Laph evs
      uint nSinkQudaBatch;      // number of quark sinks to project at a time as one batch
      uint nEigQudaBatch;       // number of Laph evs to project at a time as one batch
      bool useMultiSrcInverter;  
   };

   typedef std::vector<std::complex<double>> DataType;   // store in double precision even if single precision

 private:

	 // pointers to internal infos (managed by this handler
	 // with new and delete)

	 const GaugeConfigurationInfo *uPtr;
	 const GluonSmearingInfo *gSmearPtr;
	 const QuarkSmearingInfo *qSmearPtr;
	 const QuarkActionInfo *qactionPtr;
	 const FileListInfo *fPtr;
	 const FileListInfo *fPtrSparseGrid;
	 const InverterInfo *invertPtr;
	 uint Nspin;
	 Mode mode;
	 void* preconditioner;

	 // structure containing the computations to perform
	 PerambComputations perambComps;

	 // necessary quda data
	 QudaInvertParam quda_inv_param;

	 // sub-handler pointers
   static std::unique_ptr<SparseGridHandler> sgHandler; 	 
	 static std::unique_ptr<QuarkSmearingHandler> qSmearHandler;
	 static std::unique_ptr<GaugeConfigurationHandler> gaugeHandler;

	 // Prevent copying ... handler might contain large
	 // amounts of data

	 PerambulatorHandler(const PerambulatorHandler&);
	 PerambulatorHandler& operator=(const PerambulatorHandler&);


	 // data I/O handler pointers

	 DataPutHandlerMF<PerambulatorHandler,FileKey,RecordKey,DataType> *DHputPtr;
	 DataPutHandlerMF<SparseGridHandler,FileKey,RecordKey,DataType> *DHputPtrSpGrid;

	 DataGetHandlerMF<PerambulatorHandler,FileKey,RecordKey,DataType> *DHgetPtr;


 public:


	 PerambulatorHandler();

	 PerambulatorHandler(const GaugeConfigurationInfo& gaugeinfo,
			 const GluonSmearingInfo& gluonsmear,
			 const QuarkSmearingInfo& quarksmear,
			 const QuarkActionInfo& quark,
			 const RandomSparseGrid& rsgrid, 
			 const FileListInfo& flist,
			 const FileListInfo& flist_sparse_grid,
			 const std::string& smeared_quark_filestub,
			 bool upper_spin_components_only=false,
			 Mode in_mode=ReadOnly,
			 const std::string& gauge_str="default_gauge_field");

	 void setInfo(const GaugeConfigurationInfo& gaugeinfo,
			 const GluonSmearingInfo& gluonsmear,
			 const QuarkSmearingInfo& quarksmear,
			 const QuarkActionInfo& quark,
			 const RandomSparseGrid& rsgrid, 
			 const FileListInfo& flist,
			 const FileListInfo& flist_sparse_grid,
			 const std::string& smeared_quark_filestub,
			 bool upper_spin_components_only=false,
			 Mode in_mode=ReadOnly,
			 const std::string& gauge_str="default_gauge_field");

	 ~PerambulatorHandler();

	 void clear();

	 bool isInfoSet() const;

	 const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

	 const GluonSmearingInfo& getGluonSmearingInfo() const;

	 const QuarkSmearingInfo& getQuarkSmearingInfo() const;

	 const QuarkActionInfo& getQuarkActionInfo() const;

	 const FileListInfo& getFileListInfo() const;

	 const FileListInfo& getFileListInfoSparseGrid() const;

	 uint getNumberOfLaplacianEigenvectors() const;

	 int getTimeExtent() const;

	 void getHeader(XMLHandler& xmlout) const;

	 void getFileMap(XMLHandler& xmlout) const;

	 void outputSuffixMap();

	 std::map<int,FileKey> getSuffixMap() const;

	 void outputSuffixMap(XMLHandler& xmlout);

   void setInverter(const InverterInfo& invinfo);

   const InverterInfo& getInverterInfo() const;

   bool setUpPreconditioning(QudaInvertParam& invParam);

        // compute quark perambulators

   void clearComputationSet();
      
   void setComputationSet(const XMLHandler& xmlcmp);

   void computePerambulators(bool extra_soln_check=false,
                             bool print_coeffs=false,
                             bool report_gflops=false);

        // get routines

   const DataType& getData(int snk_time, int snk_spin, int src_time, int src_spin,
                           int src_eigvec_index) const;

   Array<std::complex<double>> getFullData(int snk_time, int snk_spin, 
                                           int src_time, int src_spin, int nEigsUse) const;

   bool queryData(int snk_time, int snk_spin, int src_time, int src_spin,
                  int src_eigvec_index) const;

   bool queryFullData(int snk_time, int snk_spin, int src_time, int src_spin, int nEigsUse) const;

        // merge data

   void mergeData(const FileListInfo& input_files);
   
        // check data
   
   void setChecks(const XMLHandler& xmlin);
   
   void doChecks(const std::string& logfilestub, bool verbose_output=true);
   

 private:


   void set_info(const GaugeConfigurationInfo& gaugeinfo,
                 const GluonSmearingInfo& gluonsmear,
                 const QuarkSmearingInfo& quarksmear,
                 const QuarkActionInfo& quark,
								 const RandomSparseGrid& rsgrid, 
                 const FileListInfo& flist,
                 const std::string& smeared_quark_filestub,
                 bool upper_spin_components_only,
                 const std::string& gauge_str, Mode in_mode);


   bool checkHeader(XMLHandler& xmlr, int suffix);
   void writeHeader(XMLHandler& xmlout, const FileKey& fkey,
                    int suffix);

   void check_info_set(const std::string& name) const;

         //  sub-handler connections

   void connectGaugeConfigurationHandler();
   void connectQuarkSmearingHandler(const std::string& smeared_quark_filestub);

   void disconnectGaugeConfigurationHandler();
   void disconnectQuarkSmearingHandler();

   void computePerambulators(int src_time, const std::set<int>& src_evindices,
                             const std::vector<void*>& evList, bool print_coeffs,
                             bool extra_soln_check, bool report_gflops, 
                             double& makesrc_time, double& inv_time,
                             double& evproj_time, double& write_time,
                             bool use_multisrc_inverter);

   void computePerambulatorsMS(int src_time, const std::set<int>& src_evindices,
                               const std::vector<void*>& evList, bool print_coeffs,
                               bool extra_soln_check, bool report_gflops, 
                               double& makesrc_time, double& inv_time,
                               double& evproj_time, double& write_time);

   void computePerambulatorsSS(int src_time, const std::set<int>& src_evindices,
                               const std::vector<void*>& evList, bool print_coeffs,
                               bool extra_soln_check, bool report_gflops, 
                               double& makesrc_time, double& inv_time,
                               double& evproj_time, double& write_time);

          // Makes the source in the Dirac-Pauli basis.  The source is
	  // gamma_4 times the usual source since we use the chi = psi-bar gamma_4
	  // field operator.   "src_spin" is 1,2,3,4

   void make_source(LattField& ferm_src, const void* ev_src_ptr, 
                    int src_time, int src_spin);

   void doACheck(const PerambComputation& pcomp, int src_spin,
                 const std::string& logfile, bool verbose_output);

   friend class DataPutHandlerMF<PerambulatorHandler,FileKey,RecordKey,DataType>;
   friend class DataGetHandlerMF<PerambulatorHandler,FileKey,RecordKey,DataType>;

};



// **************************************************************************************
}
#endif  
