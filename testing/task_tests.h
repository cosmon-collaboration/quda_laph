#ifndef QUDA_LAPH_TESTS_H
#define QUDA_LAPH_TESTS_H
#ifdef TESTING

#include "xml_handler.h"
#include "latt_field.h"
#include "laph_noise.h"
#include "gauge_configuration_info.h"
#include "quark_action_info.h"

using namespace LaphEnv;


namespace QLTestEnv {

// ************************************************

void doQudaLaphTests(XMLHandler& xml_rdr);


// ************************************************

void testXMLHandler(XMLHandler& xmlr);

void testLayoutInfo(XMLHandler& xmlr);

void testIOHandlerFM(XMLHandler& xmlr);

void testLatticeField(XMLHandler& xmlr);

void testNamedObjectMap(XMLHandler& xmlr);

void testReadGaugeConfig(XMLHandler& xmlr);

void testIOMap(XMLHandler& xmlr);

void testArray(XMLHandler& xmlr);

void testFileListInfo(XMLHandler& xmlr);

void testDataIOHandler(XMLHandler& xmlr);

void testMakeQuarkSource(XMLHandler& xmlr);

void testShift(XMLHandler& xmlr);

void testCloverDiracMatrix(XMLHandler& xmlr);

void testSpinConversions(XMLHandler& xmlr);

void testMinusSpatialLaplacian(XMLHandler& xmlr);

/*
void testmulticompare(XMLHandler& xmlr);

void testmultilooper(XMLHandler& xmlr);

void testMomenta(XMLHandler& xmlr);

void testGaugeConfigurationInfo(XMLHandler& xmlr);

void testDilutionSchemeInfo(XMLHandler& xmlr);

void testMCEnsembleInfo(XMLHandler& xmlr);

void testLaphNoiseInfo(XMLHandler& xmlr);

void testBaryonOperatorInfo(XMLHandler& xmlr);

void testTetraquarkOperatorInfo(XMLHandler& xmlr);

void testMesonOperatorInfo(XMLHandler& xmlr);

void testGlueballOperatorInfo(XMLHandler& xmlr);

void testOperatorInfo(XMLHandler& xmlr);

void testCorrelatorInfo(XMLHandler& xmlr);

void testQuarkLineEndInfo(XMLHandler& xmlr);

void testQuarkInfo(XMLHandler& xmlr);

void testQuarkInfoMap(XMLHandler& xmlr);

void testQuarkActionInfo(XMLHandler& xmlr);

void testDilutionHandler(XMLHandler& xmlr);

void testOperatorHandler(XMLHandler& xmlr);

void testHadronHandler(XMLHandler& xmlr);

void testHadronLineEndInfo(XMLHandler& xmlr);


void testTensorDataBuilder(XMLHandler& xmlr);

void testContractions(XMLHandler& xmlr);

void testBLASContractions(XMLHandler& xmlr);


void testXMLHandler(XMLHandler& xmlr);

void testTimesInfo(XMLHandler& xmlr);

void testIOHandler(XMLHandler& xmlr);

void testLaphEigenSolverInfo(XMLHandler& xml_in);

void testInverterInfo(XMLHandler& xml_in); 

void testGluonSmearingInfo(XMLHandler& xml_in);

void testQuarkSmearingInfo(XMLHandler& xml_in);

void convertQuarkSinksToNewFormat(XMLHandler& xml_in);

void checkEqualQuarkSinks(XMLHandler& xml_in);

void testLayout(XMLHandler& xml_in);

void testPIOHandler(XMLHandler& xml_in);


#if (QDP_ND == 4)

void makeGaugeTransformation(XMLHandler& xml_in);

void applyGaugeTransformation(XMLHandler& xml_in);

void printSmearedLinks(const multi1d<LatticeColorMatrix>& usmear,
                       const GaugeConfigurationInfo& gauge_info,
                       const std::string& outfile);

void printAllSmearedLinks(const multi1d<LatticeColorMatrix>& usmear,
                          const GaugeConfigurationInfo& gauge_info,
                          const std::string& outfile);

void printLaphEigenvectors(const GluonSmearingInfo& gSmear,
                           const GaugeConfigurationInfo& gaugeinfo,
                           const QuarkSmearingInfo& qSmear,
                           const std::string& smeared_quark_filestub,
                           const std::string& outfile);

void printQuarkSinks(const LaphNoiseInfo& noise, 
                     int time_proj_index, int dilution_index, 
                     const multi1d<LatticeColorVector>& phi,
                     const GaugeConfigurationInfo& gaugeinfo);

void testPion(XMLHandler& xml_in);

void computePerambulators(XMLHandler& xml_in);



#elif (QDP_ND == 3)

void applyGaugeTransformation(XMLHandler& xml_in);

void printSmearedLinks(const GaugeConfigurationInfo& gauge_info,
                       const GluonSmearingInfo& gmear,
                       const std::string& smeared_gauge_filename,
                       const std::string& outfile);

void printLaphEigenvectors(LaphVectorHandler& laphEigvecs, int t,
                           const GaugeConfigurationInfo& gaugeinfo,
                           int nEigvecs, TextFileWriter& fout);

void printQuarkSinks(QuarkHandler& Q,
                     const LaphNoiseInfo& noise,
                     int time_proj_index, int ndilutions,
                     const GaugeConfigurationInfo& gaugeinfo);

void printQuarkSinkCoefficients(QuarkHandler& Q,
                     const LaphNoiseInfo& noise,
                     int time_proj_index, int ndilutions,
                     const GaugeConfigurationInfo& gaugeinfo);

void testQuarkSinkRead(XMLHandler& xml_in);

void testBaryonRead(XMLHandler& xml_in);

void testMomentumSums(XMLHandler& xml_in);

void testMesonIntLoops(XMLHandler& xml_in);

void testLaphPion(XMLHandler& xml_in);

void autoLaphPion(XMLHandler& xml_in);

void autoLaphMesons(XMLHandler& xml_in);

void testLaphNucleon(XMLHandler& xml_in);

void autoLaphNucleon(XMLHandler& xml_in);

void testLaphBaryon123uds(XMLHandler& xml_in);

void compressQuarkSinks(XMLHandler& xml_in);

void testSmearedPion(XMLHandler& xml_in);

void testSmearedNucleon(XMLHandler& xml_in);

void testSmearedBaryon123uds(XMLHandler& xml_in);

void testLaphPion2(XMLHandler& xml_in);

void testBaryonLineCalc(XMLHandler& xmlr);

void checkEqualMesonSinks(XMLHandler& xml_in);

void checkEqualBaryonSinks(XMLHandler& xml_in);

void convertSmearedGaugeToNewFormat(XMLHandler& xml_in);

#endif
*/


class LatticeAssigner
{
 
    bool random_set;

    uint mtseed;
    
    UniformDeviate32 Q;

    double minval,maxval,range;

    static bool dprec;
    
    double gcoef;

 public:
 
    LatticeAssigner(uint in_mtseed=0, double in_minval=-10.0, double in_maxval=10.0, 
                    bool in_random=false);

    void reSeed(uint in_mtseed);

    static bool isDblePrec() {return dprec;}

    void assign_field(LattField& latfield, const std::string& fieldname);

    std::vector<int> getRandomSite();

 private:

    double generate();

    void setup_gen();

    uint seedval(uint inseed, bool inrandom);

    void assign_site(const std::vector<int>& coord, std::vector<char>& result, int cmplx_per_site);

    void site_basic_assigner(const std::vector<int>& coord, std::vector<std::complex<double>>& result);

    void site_random_assigner(const std::vector<int>& coord, std::vector<std::complex<double>>& result);

    void reunitarize(std::vector<std::complex<double>>& result);

    void cross_prod(const std::complex<double>& A1, const std::complex<double>& A2, const std::complex<double>& A3,
                    const std::complex<double>& B1, const std::complex<double>& B2, const std::complex<double>& B3,
                    std::complex<double>& C1, std::complex<double>& C2, std::complex<double>& C3);

    friend class LatticeChecker;
};


class LatticeChecker
{

    LatticeAssigner& latassigner;

 public:

    LatticeChecker(LatticeAssigner& inassigner) : latassigner(inassigner) {}

    void reSeed(uint in_mtseed){ latassigner.reSeed(in_mtseed);}

         // tol<=0.0 means set tolerance from whether double or single precision
    bool check_field(const LattField& latfield, const std::string& fieldname, bool full,
                     double tol=-0.1);

 private:

    template <typename T>
    bool do_check_field(const LattField& latfield, const std::string& fieldname, bool full, const T& tolerance);

    template <typename T>
    bool check_field_full(const LattField& latfield, const std::string& fieldname, const T& tolerance);
 
    template <typename T>
    bool check_field_random(const LattField& latfield, const std::string& fieldname, const T& tolerance);
 
    template <typename T>
    bool do_site_checker(const LattField& latfield, const std::vector<int>& coord, const T& tolerance);

};


class GaugeFieldAssigner
{
    LatticeAssigner LA;

    std::vector<uint> mtseeds;

 public:

    GaugeFieldAssigner(uint inseed0, uint inseed1, uint inseed2, uint inseed3,
                       double in_minval=-10.0, double in_maxval=10.0)
        : LA(0,in_minval,in_maxval), mtseeds(4)
          {mtseeds[0]=inseed0; mtseeds[1]=inseed1; mtseeds[2]=inseed2; mtseeds[3]=inseed3;}
    
    void reSeed(uint inseed0, uint inseed1, uint inseed2, uint inseed3)
     {mtseeds[0]=inseed0; mtseeds[1]=inseed1; mtseeds[2]=inseed2; mtseeds[3]=inseed3;}

    void assign_gauge_field(std::vector<LattField>& latfield, const std::string& fieldname);

    friend class GaugeFieldChecker;
};


class GaugeFieldChecker
{

    GaugeFieldAssigner& GF_assigner;

    LatticeChecker LC;

 public:

    GaugeFieldChecker(GaugeFieldAssigner& inassigner) : GF_assigner(inassigner), LC(GF_assigner.LA) {}

    void reSeed(uint in_mtseed1, uint in_mtseed2, uint in_mtseed3, uint in_mtseed4)
    {GF_assigner.reSeed(in_mtseed1,in_mtseed2,in_mtseed3,in_mtseed4);}

         // tol<=0.0 means set tolerance from whether double or single precision
    bool check_gauge_field(const std::vector<LattField>& gffield, const std::string& fieldname, 
                           bool full, double tol=-0.1);

};


  // This is useful for debugging only; very slow.  "nsites" is the number of
  // sites to print in order from (0,0,0,0).  An input value of -1 means print all.

void printField(const LattField& field, const std::string& fieldname, int nsites=-1);

void printFieldToFile(const LattField& field, const std::string& fieldname, const std::string& filename);

void printField(const LattField& field, const std::string& fieldname, 
                const std::vector<std::vector<int>>& sites);

void compare_fields(const LattField& src1, const LattField& src2);

void convertSpinBasis(LattField& outferm, const LattField& inferm, const std::string& from_to);

void setUnitField(LattField& field);

void setZeroField(LattField& field);

void setConstantField(LattField& field, const std::complex<double>& zconst);

      //  This routine is slow: meant only for debugging, testing
      //   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
      //   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

void shift(LattField& outfield, const LattField& infield, int dir, char fwd_or_bwd);

   //   Applies the clover Dirac operation to "infield", returning
   //   the result in "outfield".  This operation is
   //
   //    outfield =  [ 1/(2*kappa) - (1/2) Dterm  + CFterm ] infield
   //
   //   where
   //
   //        Dterm = sum_mu [ (1-gamma_mu) U  + (1+gamma_mu) * U^dag ]
   //
   //        CFterm = csw (i/4)  sigma[mu,nu] F[mu,nu]
   //
   //            sigma[mu,nu] (i/2) [gamma_mu, gamma_nu]
   //
   //            F[mu,nu] = (1/8) ( Q[mu,nu]-Q[nu,mu] )
   //
   //            Q[mu,nu] = U(mu,nu,-mu,-nu) + U(nu,-mu,-nu,mu)
   //                     + U(-mu,-nu,mu,nu) + U(-nu,mu,nu,-mu)
   //
   //   Lattice shifts must take the fermion temporal boundary
   //   conditions into account.
   

void applyCloverDirac(LattField& outfield, const LattField& infield,
                      const std::vector<LattField>& gauge_field,
                      const GaugeConfigurationInfo& gaction,
                      const QuarkActionInfo& qaction);


    //  Applies the 3d spatial Laplacian with a smeared gauge field
    //  onto "infield", returning result in "outfield".  These fields
    //  must be color vectors

void applyMinusSpatialLaplacian(LattField& outfield, const LattField& infield,
                                const std::vector<LattField>& smeared_gauge_field);


    //  This routine applies an inner product conj(leftfield).rightfield,
    //  where leftfield and rightfield are color-vector fields, but the
    //  inner product is taken over the time slices.  The ntime inner
    //  products are returned.

std::vector<std::complex<double>> getTimeSlicedInnerProducts(const LattField& leftfield, 
                                                             const LattField& rightfield);


// ***********************************************
}
#endif
#endif
