#include "tasks.h"
#include "perambulator_handler.h"
#include "stop_watch.h"

using namespace std;

namespace LaphEnv {

	
// *****************************************************************************************
// *                                                                                       *
// *   Task to compute the LapH quark perambulators and write to file(s).                  *
// *   The input XML is expected to have the following  form:                              *
// *                                                                                       *
// *      <Task>                                                                           *
// *         <Name>LAPH_QUARK_PERAMBULATORS</Name>                                         *
// *         <QuarkPerambulatorInfo>                                                       *
// *            <GaugeConfigurationInfo>                                                   *
// *               <EnsembleName>CLS_A653</EnsembleName>                                   *
// *               <FileFormat>CERN</FileFormat>                                           *
// *               <ConfigType>WilsonImproved</ConfigType>                                 *
// *               <FileName>/path/A653r000n1</FileName>                                   *
// *               <ConfigNumber>1</ConfigNumber>                                          *
// *               <FermionTimeBC>antiperiodic</FermionTimeBC>                             *
// *            </GaugeConfigurationInfo>                                                  *
// *            <GluonStoutSmearingInfo>                                                   *
// *               <LinkIterations>20</LinkIterations>                                     *
// *               <LinkStapleWeight>0.1</LinkStapleWeight>                                *
// *            </GluonStoutSmearingInfo>                                                  *
// *            <QuarkLaphSmearingInfo>                                                    *
// *               <LaphSigmaCutoff>0.32</LaphSigmaCutoff>                                 *
// *               <NumberLaphEigvecs> 32 </NumberLaphEigvecs>                             *
// *            </QuarkLaphSmearingInfo>                                                   *
// *            <SmearedQuarkFileStub>/path/sq_field</SmearedQuarkFileStub> (see below)    *
// *            <QuarkActionInfo>                                                          *
// *               <Name> WILSON_CLOVER </Name>                                            *
// *               <Flavor> ud </Flavor>                                                   *
// *               <Kappa>0.1365716</Kappa>                                                *
// *               <clovCoeff>2.0668577</clovCoeff>                                        *
// *               <TimeBC> antiperiodic </TimeBC>                                         *
// *             </QuarkActionInfo>                                                        *
// *            <FileListInfo>/path/quark_perambs_1</FileNameStub>                         *
// *               <MinFileNumber>0</MinFileNumber>                                        *
// *               <MaxFileNumber>200</MaxFileNumber>                                      *
// *            </FileListInfo>                                                            *
// *            <InverterInfo>                                                             *
// *               <Name>BICGSTAB</Name>                                                   *
// *               <Tolerance>1.0e-11</Tolerance>                                          *
// *               <MaxIterations>35000</MaxIterations>                                    *
// *            </InverterInfo>                                                            *
// *            <ExtraSolutionChecks/> (optional)                                          *
// *         </QuarkPerambulatorInfo>                                                      *
// *         <ComputationSet>                                                              *
// *            <NumSinksBeforeProject>  4 </NumSinksBeforeProject>                        *
// *            <NumSinksInProjectBatch> 4 </NumSinksInProjectBatch>                       *
// *            <NumEigsInProjectBatch>  4 </NumEigsInProjectBatch>                        *
// *            <Computation>                                                              *
// *               <SourceTime>8</SourceTime>                                              *
// *               <SourceLaphEigvecIndices>12 17 23 24 25</SourceLaphEigvecIndices>       *
// *            </Computation>                                                             *
// *            <Computation>                                                              *
// *               <SourceTime>12</SourceTime>                                             *
// *               <SourceLaphEigvecIndexMin>24</SourceLaphEigvecIndexMin>                 *
// *               <SourceLaphEigvecIndexMax>32</SourceLaphEigvecIndexMax>                 *
// *            </Computation>                                                             *
// *         </ComputationSet>                                                             *
// *         <Verbosity>full</Verbosity>  (optional: override default)                     *
// *         <PrintCoefficients/>  (optional)                                              *
// *      </Task>                                                                          *
// *                                                                                       *
// *   If the tag <SmearedQuarkFileStub> is set, then the LapH eigenvectors will be read   *
// *   from file.  Otherwise, this task expects that the LapH eigenvectors already         *
// *   reside in HostGlobal; otherwise, an error results.                                  *
// *                                                                                       *
// *   The default verbosity from the main program is used, unless overridden by           *
// *   the <Verbosity> tag above.                                                          *
// *                                                                                       *
// *   If the <PrintCoefficients/> is present, then the coefficients are printed           *
// *   to standard output.                                                                 *
// *                                                                                       *
// *   Often, all of the perambulators for all source Laph eigenvector indices cannot be   *
// *   done in a single run, and so for safety, the data for different sets of source ev   *
// *   indices will get saved into separate files.  These then need to be combined into    *
// *   a single file for subsequent reading to compute correlators.  The <Merge> task      *
// *   is needed to accomplish this.                                                       *
// *                                                                                       *
// *****************************************************************************************


void doLaphQuarkPerambulators(XMLHandler& xmltask)
{
 if (xml_tag_count(xmltask,"QuarkPerambulatorInfo")!=1){
    errorLaph("Must have one <QuarkPerambulatorInfo> tag");}
 XMLHandler xmlr(xmltask,"QuarkPerambulatorInfo");
 GaugeConfigurationInfo gaugeinfo(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 string smeared_quark_filestub;
 xmlreadif(xmlr,"SmearedQuarkFileStub",smeared_quark_filestub,"LAPH_QUARK_PERAMBULATORS");
 QuarkActionInfo quark(xmlr);
 FileListInfo files(xmlr);
 InverterInfo invinfo(xmlr);
 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>0){
    upper_spin_only=true;}
 bool extra_soln_check=false;
 if (xml_tag_count(xmlr,"ExtraSolutionChecks")>0){
    extra_soln_check=true;}
    // change verbosity from default if requested
 Verbosity task_verbosity(getVerbosity());
 if (xml_read_if(xmltask,task_verbosity)){
    setVerbosity(task_verbosity.getQudaValue());}
 bool print_coeffs=false;
 if (xml_tag_count(xmltask,"PrintCoefficients")>0){
    print_coeffs=true;}
 PerambulatorHandler::Mode mode=PerambulatorHandler::Compute;

 printLaph("\n");
 printLaph(" ***********************************************************");
 printLaph(" *                                                         *");
 printLaph(" *   Laph Task: Compute the quark perambulators            *");
 printLaph(" *              and write to file                          *");
 printLaph(" *                                                         *");
 printLaph(" ***********************************************************\n");
 printLaph(make_strf("\n%s\n",gaugeinfo.output()));
 printLaph(make_strf("\n\nGluon Smearing:\n%s\n",gSmear.output()));
 printLaph(make_strf("\n\nQuark Smearing:\n%s\n",qSmear.output()));
 if (!smeared_quark_filestub.empty()){
    printLaph(make_strf("SmearedQuarkFileStub: %s",smeared_quark_filestub));}
 printLaph(make_str("\nQuarkAction:\n",quark.output()));
 printLaph(make_str("\nInverter Info:\n",invinfo.output()));
 if (upper_spin_only){
    printLaph("Only upper spin components used");}
 else{
    printLaph("All spin components used");}
 if (extra_soln_check){
    printLaph("Extra solution checks will be performed");}

     // create handler
 PerambulatorHandler Q(gaugeinfo,gSmear,qSmear,quark,files,
                       smeared_quark_filestub,upper_spin_only,mode);

    // read the set of computations (time sources, src eigvec indices)
    // as well as the batching parameters
 XMLHandler xmlcmp(xmltask,"ComputationSet");
 Q.setComputationSet(xmlcmp);

     // set the inverter info

 Q.setInverter(invinfo);
 printLaph("Inverter initialized in PerambulatorHandler");
 Q.outputSuffixMap();

     // now do the computations!
 StopWatch outer; outer.start();
 Q.computePerambulators(extra_soln_check,print_coeffs);
 outer.stop();
 printLaph(make_strf("LAPH_PERAMBULATORS: total time = %g secs",
           outer.getTimeInSeconds()));
 printLaph("LAPH_PERAMBULATORS: ran successfully\n"); 
} 


// *****************************************************************************
}
