#include "dilution_scheme_info.h"
#include "filelist_info.h"
#include "gauge_configuration_info.h"
#include "inverter_info.h"
#include "quark_action_info.h"
#include "quark_handler.h"
#include "quark_smearing_handler.h"
#include "stop_watch.h"
#include "tasks.h"

using namespace std;
using namespace quda;

namespace LaphEnv {

// *****************************************************************************************
// * *
// *   Task to compute stochastic LapH quark line sinks and write to file or *
// *   the NamedObjMap.   The input XML is expected to have the following  form:
// *
// * *
// *      <Task> *
// *         <Name>LAPH_QUARK_LINE_ENDS</Name> *
// *         <QuarkLineEndInfo> *
// *            <GaugeConfigurationInfo> *
// *               <EnsembleName>CLS_A653</EnsembleName> *
// *               <FileFormat>CERN</FileFormat> *
// *               <ConfigType>WilsonImproved</ConfigType> *
// *               <FileName>/path/A653r000n1</FileName> *
// *               <ConfigNumber>1</ConfigNumber> *
// *               <FermionTimeBC>antiperiodic</FermionTimeBC> *
// *            </GaugeConfigurationInfo> *
// *            <GluonStoutSmearingInfo> *
// *               <LinkIterations>20</LinkIterations> *
// *               <LinkStapleWeight>0.1</LinkStapleWeight> *
// *            </GluonStoutSmearingInfo> *
// * <SmearedGaugeFileName>/path/smeared_gauge_field_1</SmearedGaugeFileName> *
// *            <QuarkLaphSmearingInfo> *
// *               <LaphSigmaCutoff>0.32</LaphSigmaCutoff> *
// *               <NumberLaphEigvecs> 32 </NumberLaphEigvecs> *
// *            </QuarkLaphSmearingInfo> *
// * <SmearedQuarkFileStub>/path/smeared_quark_field_1</SmearedQuarkFileStub> *
// *            <QuarkActionInfo> *
// *               <Name> WILSON_CLOVER </Name> *
// *               <Flavor> ud </Flavor> *
// *               <Kappa>0.1365716</Kappa> *
// *               <clovCoeff>2.0668577</clovCoeff> *
// *               <TimeBC> antiperiodic </TimeBC> *
// *             </QuarkActionInfo> *
// *            <FileListInfo> *
// *               <FileNameStub>/path/quark_sinks_1</FileNameStub> *
// *               <MinFileNumber>0</MinFileNumber> *
// *               <MaxFileNumber>100</MaxFileNumber> *
// *            </FileListInfo> *
// *            <InverterInfo> *
// *               <Name>BICGSTAB</Name> *
// *               <Tolerance>1.0e-11</Tolerance> *
// *               <MaxIterations>35000</MaxIterations> *
// *            </InverterInfo> *
// *            <LaphDilutionScheme> *
// *               <TimeDilution> *
// *                  <DilutionType>full</DilutionType> *
// *               </TimeDilution> *
// *               <EigvecDilution> *
// *                  <DilutionType> interlace </DilutionType> *
// *                  <NumberProjectors> 4 </NumberProjectors> *
// *               </EigvecDilution> *
// *               <SpinDilution> *
// *                  <DilutionType>full</DilutionType> *
// *               </SpinDilution> *
// *            </LaphDilutionScheme> *
// *            <Verbosity>full</Verbosity> *
// *            <ExtraSolutionChecks/> (optional) *
// *         </QuarkLineEndInfo> *
// *         <SinkComputations> *
// *            <NumSinksBeforeProject>  4 </NumSinksBeforeProject> *
// *            <NumSinksInProjectBatch> 4 </NumSinksInProjectBatch> *
// *            <NumEigsInProjectBatch>  4 </NumEigsInProjectBatch> *
// *            <NoiseList_TimeProjIndexList> *
// *               <LaphNoiseList> *
// *                  <LaphNoiseInfo> *
// *                     <Seed>21117</Seed> *
// *                     <ZNGroup>4</ZNGroup> *
// *                  </LaphNoiseInfo> *
// *                  <LaphNoiseInfo> *
// *                     <ZNGroup>4</ZNGroup> *
// *                     <Seed>9494</Seed> *
// *                  </LaphNoiseInfo> *
// *                  <LaphNoiseInfo> *
// *                     <ZNGroup>4</ZNGroup> *
// *                     <Seed>20019</Seed> *
// *                  </LaphNoiseInfo> *
// *               </LaphNoiseList> *
// *               <TimeProjIndexList> *
// *                  <Values>23</Values> *
// *               </TimeProjIndexList> *
// *            </NoiseList_TimeProjIndexList> *
// *         </SinkComputations> *
// *      </Task> *
// * *
// *   If the tag "Verbosity" is included with value "full", then the quark sink
// *
// *   solutions will be echoed to standard output as coefficients of the LapH *
// *   eigenvectors. *
// * *
// *****************************************************************************************

void doLaphQuarkLineEnds(XMLHandler &xmltask) {
  if (xml_tag_count(xmltask, "QuarkLineEndInfo") != 1) {
    errorLaph("Must have one <QuarkLineEndInfo> tag");
  }
  XMLHandler xmlr(xmltask, "QuarkLineEndInfo");
  GaugeConfigurationInfo gaugeinfo(xmlr);
  GluonSmearingInfo gSmear(xmlr);
  QuarkSmearingInfo qSmear(xmlr);
  string smeared_quark_filestub;
  xmlread(xmlr, "SmearedQuarkFileStub", smeared_quark_filestub,
          "LAPH_QUARK_LINE_ENDS");
  string smeared_gauge_filename;
  xmlread(xmlr, "SmearedGaugeFileName", smeared_gauge_filename,
          "LAPH_QUARK_LINE_ENDS");
  DilutionSchemeInfo dil(xmlr);
  QuarkActionInfo quark(xmlr);
  FileListInfo files(xmlr);
  InverterInfo invinfo(xmlr);
  bool verbose = false; // output final results
  string verbosity;
  xmlreadif(xmlr, "Verbosity", verbosity, "LAPH_QUARK_LINE_ENDS");
  if (tidyString(verbosity) == "full") {
    verbose = true;
    setVerbosity(QUDA_VERBOSE);
  } else if (verbosity == "low") {
    setVerbosity(QUDA_SUMMARIZE);
  } else if (verbosity == "none") {
    setVerbosity(QUDA_SILENT);
  }
  bool extra_soln_check = false;
  if (xml_tag_count(xmlr, "ExtraSolutionChecks") >= 1) {
    extra_soln_check = true;
  }

  printLaph("\n");
  printLaph(" ***********************************************************");
  printLaph(" *                                                         *");
  printLaph(" *   Laph Task: Compute the quark line ends                *");
  printLaph(" *              and write to file as time slices           *");
  printLaph(" *                                                         *");
  printLaph(" ***********************************************************\n");
  printLaph(make_strf("\n%s\n", gaugeinfo.output()));
  printLaph(make_strf("\n\nGluon Smearing:\n%s\n", gSmear.output()));
  printLaph(make_strf("\n\nQuark Smearing:\n%s\n", qSmear.output()));
  printLaph(make_strf("SmearedQuarkFileStub: %s", smeared_quark_filestub));
  printLaph(make_strf("SmearedGaugeFileName: %s", smeared_gauge_filename));
  printLaph(make_strf("\nDilution Scheme Info:\n%s\n", dil.output()));
  printLaph(make_str("\nQuarkAction:\n", quark.output()));
  printLaph(make_str("\nInverter Info:\n", invinfo.output()));
  if (extra_soln_check) {
    printLaph("Extra solution checks will be performed");
  }
  bool compute_mode = true;

  // create quark handler
  QuarkHandler Q(gaugeinfo, gSmear, qSmear, dil, quark, files,
                 smeared_quark_filestub, smeared_gauge_filename, compute_mode);

  // read the list of computations (noises, time sources, file indices)
  // as well as the batching parameters
  XMLHandler xmlcmp(xmltask, "SinkComputations");
  Q.setSinkComputations(xmlcmp);

  // set the inverter info

  Q.setInverter(invinfo);
  printLaph("Inverter initialized in QuarkHandler");
  Q.outputSuffixMap();

  // now do the computations!
  StopWatch outer;
  outer.start();
  Q.computeSinks(verbose, extra_soln_check);
  outer.stop();
  printLaph(make_strf("LAPH_QUARK_LINE_ENDS: total time = %g secs",
                      outer.getTimeInSeconds()));
  printLaph("LAPH_QUARK_LINE_ENDS: ran successfully\n");
}

// ******************************************************************
} // namespace LaphEnv
