#ifndef QUARK_SMEARING_HANDLER_H
#define QUARK_SMEARING_HANDLER_H

#include "gauge_configuration_info.h"
#include "data_io_handler.h"
#include "laph_eigen_info.h"

namespace LaphEnv {


// *****************************************************************
// *                                                               *
// *  QuarkSmearingHandler:                                        *
// *                                                               *
// *     -- computes the Laplacian eigenvectors, saving to the     *
// *          HostGlobal or writing to files                       *
// *     -- also used to read from files                           *
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
// *    set...()       to internally set from file or HostGlobal   *
// *                                                               *
// *    get...()       provides access to results                  *
// *                                                               *
// *****************************************************************


class QuarkSmearingHandler
{

   struct TimeKey
   {
    int value;

    TimeKey() : value(0) {}
    TimeKey(int in_val) : value(in_val) {}
    TimeKey(XMLHandler& xmlr)
    {xmlread(xmlr,"TimeSlice",value,"QuarkSmearingHandler::TimeKey");}
    TimeKey(const TimeKey& in) : value(in.value) {}
    ~TimeKey() {}

    bool operator<(const TimeKey& rhs) const
    {return (value<rhs.value);}

    bool operator==(const TimeKey& rhs) const
    {return (value==rhs.value);}

    void output(XMLHandler& xmlw) const 
    {xmlw.set_root("TimeSlice",make_string(value));}

    explicit TimeKey(const unsigned int* buf) : value(*buf) {}
    static int numints() {return 1;} 
    size_t numbytes() const {return sizeof(unsigned int);}
    void copyTo(unsigned int* buf) const { *buf=value;}

   };

   struct LevelKey
   {
    int value;

    LevelKey() : value(0) {}
    LevelKey(int in_val) : value(in_val) {}
    LevelKey(XMLHandler& xmlr)
    {xmlread(xmlr,"Level",value,"QuarkSmearingHandler::LevelKey");}
    LevelKey(const LevelKey& in) : value(in.value) {}
    ~LevelKey() {}

    bool operator<(const LevelKey& rhs) const
    {return (value<rhs.value);}

    bool operator==(const LevelKey& rhs) const
    {return (value==rhs.value);}

    void output(XMLHandler& xmlw) const 
    {xmlw.set_root("Level",make_string(value));}

    explicit LevelKey(const unsigned int* buf) : value(*buf) {}
    static int numints() {return 1;} 
    size_t numbytes() const {return sizeof(unsigned int);}
    void copyTo(unsigned int* buf) const { *buf=value;}

   };

      // Key that combines the eigenvector level number with the
      // details of a covariant displacement.  Must fit in a 32-bit 
      // unsigned integer.
/*
   struct LevelDispKey
   {
    uint value;

    LevelDispKey() : value(0) {}
    LevelDispKey(uint in_val, const DirPath& displace);
    LevelDispKey(XMLHandler& xmlr);
    LevelDispKey(const LevelDispKey& in) : value(in.value) {}
    ~LevelDispKey() {}

    bool operator<(const LevelDispKey& rhs) const
    {return (value<rhs.value);}

    bool operator==(const LevelDispKey& rhs) const
    {return (value==rhs.value);}

    void output(XMLHandler& xmlw) const;
    void encode(uint in_val, const DirPath& displace);
    
    uint getLevel() const;
    DirPath getDisplacement() const;

    explicit LevelDispKey(const unsigned int* buf) : value(*buf) {}
    static int numints() {return 1;} 
    size_t numbytes() const {return sizeof(unsigned int);}
    void copyTo(unsigned int* buf) const { *buf=value;}

   };
*/
       // pointers to internal infos (managed by this handler
       // with new and delete)

   const GaugeConfigurationInfo *uPtr;
   const GluonSmearingInfo *gSmearPtr;
   QuarkSmearingInfo *qSmearPtr;
   std::string smearedQuarkFileStub; 
   bool m_read_mode;

       // storage and/or references to internal data, and other handlers
/*
#if (QDP_ND == 3)
   multi1d<double> Eigenvalues;
   DataGetHandlerMFO<QuarkSmearingHandler,TimeKey,LevelKey,
                     LattField> *dh_ptr;
   DataGetHandlerMFO<QuarkSmearingHandler,TimeKey,LevelDispKey,
                     LattField> *dph_ptr;
   std::string smearedDispQuarkFileStub; 
#elif (QDP_ND == 4) */

//   Array<double> Eigenvalues;
//   DataGetHandlerMF<QuarkSmearingHandler,LevelKey,LevelKey,
//                    LattField> *dh_ptr;


       // prevent copying ... handler might contain large
       // amounts of data

   QuarkSmearingHandler(const QuarkSmearingHandler&);
   QuarkSmearingHandler& operator=(const QuarkSmearingHandler&);


 public:


   QuarkSmearingHandler();

   QuarkSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                        const GaugeConfigurationInfo& gauge,
                        const QuarkSmearingInfo& quark_smearing,
                        const std::string& smeared_quark_file_stub,
                        bool read_mode=true); 

   void setInfo(const GluonSmearingInfo& gluon_smearing,
                const GaugeConfigurationInfo& gauge,
                const QuarkSmearingInfo& quark_smearing,
                const std::string& smeared_quark_file_stub,
                bool read_mode=true);

          // update if number of eigenvectors needs to be increased
//   void updateSmearing(const QuarkSmearingInfo& quark_smearing);

   ~QuarkSmearingHandler();

   void clear();      // clears everything in the handler


           // access to the info

   bool isInfoSet() const;

   const GluonSmearingInfo& getGluonSmearingInfo() const;

   const QuarkSmearingInfo& getQuarkSmearingInfo() const;

   const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

   const std::string& getSmearedQuarkFieldFileStub() const;

   uint getNumberOfLaplacianEigenvectors() const
    {return qSmearPtr->getNumberOfLaplacianEigenvectors();}

//   void getFileMap(XMLHandler& xmlout) const;

//   void outputKeys(XMLHandler& xmlout);


           // Compute the Laph Eigenvectors, store in HostGlobal, 
           // optionally output to file

   void computeLaphEigenvectors(const LaphEigenSolverInfo& solver_info,
                                const std::string& smeared_gauge_file,
                                bool print_eigvals);
/*
           // Compute covariantly displaced Laph Eigenvectors

   void displaceLaphEigenvectors(const std::set<DirPath>& displacements, int disp_length,
                                 uint timeval, const std::string& smeared_gauge_file,
                                 const std::string& disp_lapheigvec_file_stub);

   void displaceLaphEigenvectors(const std::set<DirPath>& displacements, int disp_length,
                                 uint timeval, GluonSmearingHandler& gHandler,
                                 const std::string& disp_lapheigvec_file_stub);

           // get and query data when in read mode

   const LattField& getLaphEigenvector(int time, int eigpair_num);

   bool queryLaphEigenvector(int time, int eigpair_num);

   void removeLaphEigenvector(int time, int eigpair_num);*/

   void clearLaphEigenvectors();

/*   const LattField& getDisplacedLaphEigenvector(int time, int eigpair_num, 
                                                         const DirPath& dirpath);

   bool queryDisplacedLaphEigenvector(int time, int eigpair_num, const DirPath& dirpath);

   void removeDisplacedLaphEigenvector(int time, int eigpair_num, const DirPath& dirpath);

   void clearDisplacedLaphEigenvectors();
*/
/*
#elif (QDP_ND == 4)


           // Compute the Laph Eigenvectors and store in HostGlobal

   void computeLaphEigenvectors(const LaphEigenSolverInfo& solver_info,
                                const std::string& smeared_gauge_file,
                                int striping_factor=1, int striping_unit=0,
                                int mintime=0, int maxtime=-1);   // -1 means as large as possible
*/

           // if eigenvectors are in the HostGlobal, returns a reference to them;
           // if not, reads the eigenvectors from file(s) and places them in
           // the HostGlobal

   const std::vector<LattField>& getLaphEigenvectors();

   bool queryLaphEigenvectors();
/*
   const LattField& getLaphEigenvector(int eigpair_num);

   bool queryLaphEigenvector(int eigpair_num);

   void removeLaphEigenvector(int eigpair_num);

   void clearLaphEigenvectors();
   
   void closeLaphLevelFiles();

           // checks to see if all _level files exist.

   bool checkAllLevelFilesExist();
*/
/*

#endif
*/
 private:

   void set_info(const GluonSmearingInfo& gluon_smearing,
                 const GaugeConfigurationInfo& gauge,
                 const QuarkSmearingInfo& quark_smearing,
                 const std::string& smearedQuarkFileStub,
                 bool read_mode);

   void check_info_set(const std::string& name,
                       int check_mode=0) const;

   bool loadLaphEigenvectors();

   void failure(const std::string& message);

   bool checkHeader(XMLHandler& xmlr, int suffix);
   
   void applyLaphPhaseConvention(std::vector<LattField>& laph_eigvecs);

   void checkLaphEigvecComputation(const std::vector<LattField>& laphEigvecs,
                          const std::vector<LattField>& smeared_gauge_field);

   void writeHeader(XMLHandler& xmlw, const LevelKey& fkey, int suffix);

   friend class DataGetHandlerMF<QuarkSmearingHandler,LevelKey,LevelKey,LattField>;
   friend class DataPutHandlerMF<QuarkSmearingHandler,LevelKey,LevelKey,LattField>;
};

// ***************************************************************
}
#endif  
