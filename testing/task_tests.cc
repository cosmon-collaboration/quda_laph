#include "task_tests.h"
#ifdef TESTING

namespace QLTestEnv {

// *********************************************************************

void doQudaLaphTests(XMLHandler& xml_rdr) 
{

 testXMLHandler(xml_rdr);
 testLayoutInfo(xml_rdr);
 testIOHandlerFM(xml_rdr);
 testLatticeField(xml_rdr);
 testNamedObjectMap(xml_rdr);
 testReadGaugeConfig(xml_rdr);
 testArray(xml_rdr);
 testIOMap(xml_rdr);
 testFileListInfo(xml_rdr);
 testDataIOHandler(xml_rdr);
 testMakeQuarkSource(xml_rdr);
 testShift(xml_rdr);
 testCloverDiracMatrix(xml_rdr);
 testSpinConversions(xml_rdr);
 testMinusSpatialLaplacian(xml_rdr);

// testmulticompare(xml_rdr);
// testmultilooper(xml_rdr);
// testMomenta(xml_rdr);
// testGaugeConfigurationInfo(xml_rdr);
// testLaphNoiseInfo(xml_rdr);
// testBaryonOperatorInfo(xml_rdr);
// testMesonOperatorInfo(xml_rdr);
// testGlueballOperatorInfo(xml_rdr);
// testTetraquarkOperatorInfo(xml_rdr);
// testOperatorInfo(xml_rdr);
// testCorrelatorInfo(xml_rdr);
// testDilutionSchemeInfo(xml_rdr);
// testQuarkLineEndInfo(xml_rdr);
// testQuarkActionInfo(xml_rdr);
// testQuarkInfo(xml_rdr);
// testQuarkInfoMap(xml_rdr);
// testTimesInfo(xml_rdr);
// testIOHandler(xml_rdr);
// testLaphEigenSolverInfo(xml_rdr);
// testInverterInfo(xml_rdr); */
// testGluonSmearingInfo(xml_rdr);
// testQuarkSmearingInfo(xml_rdr);
// testTensorDataBuilder(xml_rdr);
// testContractions(xml_rdr);
// testBLASContractions(xml_rdr);
// testMCEnsembleInfo(xml_rdr);
// testDilutionHandler(xml_rdr);
// testHadronHandler(xml_rdr);
// testHadronLineEndInfo(xml_rdr);
// testOperatorHandler(xml_rdr);
// testPIOHandler(xml_rdr);

/* testQuarkSinkRead(xml_rdr);
 testBaryonRead(xml_rdr);
 testMomentumSums(xml_rdr);
 testMesonIntLoops(xml_rdr);
 testLaphPion(xml_rdr);
 testLaphNucleon(xml_rdr);
 testLaphBaryon123uds(xml_rdr);
// testLaphEigenvectors(xml_rdr);
 compressQuarkSinks(xml_rdr);
 testSmearedPion(xml_rdr);
 testSmearedNucleon(xml_rdr);
 testSmearedBaryon123uds(xml_rdr);
 testLaphPion2(xml_rdr);
 autoLaphPion(xml_rdr);
 autoLaphMesons(xml_rdr);
 autoLaphNucleon(xml_rdr);
 testBaryonLineCalc(xml_rdr);
 checkEqualMesonSinks(xml_rdr);
 checkEqualBaryonSinks(xml_rdr);
 convertSmearedGaugeToNewFormat(xml_rdr); */

}

// ******************************************************************
}
#endif
