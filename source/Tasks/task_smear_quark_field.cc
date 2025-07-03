#include "quark_smearing_handler.h"
#include "stop_watch.h"

namespace LaphEnv {

// ************************************************************************************
// * *
// *   Task to smear the quark field with LapH smearing and then either write to
// *
// *   file or to the NamedObjMap. The input XML is expected to have the
// following    *
// *   form: *
// *
// *    <Task> *
// *       <Name>SMEAR_QUARK_FIELD</Name> *
// *       <GaugeConfigurationInfo> *
// *          <EnsembleName>CLS_A653</EnsembleName> *
// *          <FileFormat>CERN</FileFormat> *
// *          <ConfigType>WilsonImproved</ConfigType> *
// *          <FileName>/path/cfgs/A653r000n1</FileName> *
// *          <ConfigNumber>1</ConfigNumber> *
// *       </GaugeConfigurationInfo> *
// *       <GluonStoutSmearingInfo> *
// *          <LinkIterations>20</LinkIterations> *
// *          <LinkStapleWeight>0.1</LinkStapleWeight> *
// *       </GluonStoutSmearingInfo> *
// * <SmearedGaugeFileName>/path/smeared_gauge_field</SmearedGaugeFileName> *
// *       <QuarkLaphSmearingInfo> *
// *          <LaphSigmaCutoff>0.32</LaphSigmaCutoff> *
// *          <NumberLaphEigvecs> 164 </NumberLaphEigvecs> *
// *       </QuarkLaphSmearingInfo> *
// * <SmearedQuarkFileStub>/path/smeared_quark_field</SmearedQuarkFileStub> *
// *       <LaphEigenSolverInfo> *
// *          <ResidualTolerance>1.0e-11</ResidualTolerance> *
// *          <MaxIterations>200</MaxIterations> *
// *          <KrylovDimension>200</KrylovDimension> *
// *          <ChebyshevOrder>6</ChebyshevOrder> *
// *          <MaxEigenvalue>15.0</MaxEigenvalue> *
// *          <CutoffEigenvalue>0.5</CutoffEigenvalue> *
// *          <StartingVectorType>equal_components</StartingVectorType> *
// *          <CheckSolution>true</CheckSolution> *
// *          <OutputVerbosity>2</OutputVerbosity> *
// *       </LaphEigenSolverInfo> *
// *    </Task> *
// * *
// *   The LapH eigenvectors are removed from the gpu device memory since it is
// *
// *   assumed quark propagators will subsequently be computed, so gpu memory *
// *   will be needed for those computations. *
// * *
// ************************************************************************************

void doSmearQuarkField(XMLHandler &xmltask) {
  GaugeConfigurationInfo gaugeinfo(xmltask);
  GluonSmearingInfo gsmear(xmltask);
  std::string smeared_gauge_filename;
  xmlread(xmltask, "SmearedGaugeFileName", smeared_gauge_filename,
          "SMEAR_QUARK_FIELD");
  QuarkSmearingInfo qsmear(xmltask);
  std::string smeared_quark_filestub;
  xmlread(xmltask, "SmearedQuarkFileStub", smeared_quark_filestub,
          "SMEAR_QUARK_FIELD");

  bool laph_solve = true;
  LaphEigenSolverInfo *eigsolveinfo = 0;
  if (xml_tag_count(xmltask, "LaphEigenSolverInfo") == 1)
    eigsolveinfo = new LaphEigenSolverInfo(xmltask);
  else
    laph_solve = false;

  printLaph("\n");
  printLaph(" *************************************************************");
  printLaph(" *                                                           *");
  printLaph(" *   Laph Task: Smear the quark field by computing the       *");
  printLaph(" *              eigenvectors of the smeared-covariant        *");
  printLaph(" *              Laplacian, write to file or the NamedObjMap  *");
  printLaph(" *                                                           *");
  printLaph(
      " *************************************************************\n\n");
  printLaph(make_strf("%s", gaugeinfo.output()));
  printLaph(make_strf("\n%s", gsmear.output()));
  printLaph(
      make_strf("Smeared gauge field file name = %s", smeared_gauge_filename));
  printLaph(make_strf("\n%s", qsmear.output()));
  printLaph(
      make_strf("Smeared quark field file stub = %s", smeared_quark_filestub));
  if (laph_solve) {
    printLaph(make_strf("\n%s", eigsolveinfo->output()));
  } else {
    printLaph("\nEstimate largest eigenvalue of -Laplacian only");
  } // Not implemented yet

  // create the handler (in write mode)

  QuarkSmearingHandler Q(gsmear, gaugeinfo, qsmear, smeared_quark_filestub,
                         false);

  StopWatch outer;
  outer.start();

  if (laph_solve)
    Q.computeLaphEigenvectors(*eigsolveinfo, smeared_gauge_filename);
  /* else{
      double
     lambda_max=Q.estimateLargestLaplacianEigenvalue(smeared_gauge_filename);
      QDPIO::cout << endl<<"Estimate of largest eigenvalue of -Laplacian is "
                  << lambda_max <<endl<<endl;}
  */
  delete eigsolveinfo;
  outer.stop();
  printLaph(make_strf("SMEAR_QUARK_FIELD task: total time = %g seconds",
                      outer.getTimeInSeconds()));
  printLaph("SMEARED_QUARK_FIELD task: ran successfully\n");
}

// ******************************************************************
} // namespace LaphEnv
