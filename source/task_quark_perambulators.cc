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
// *            <UseMultiSrcInverter> yes </UseMultiSrcInverter> (optional: yes default)   *
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
// *         <ReportGflops/>  (optional)                                                   *
// *      </Task>                                                                          *
// *                                                                                       *
// *   If the tag <SmearedQuarkFileStub> is set, then the LapH eigenvectors will be read   *
// *   from file.  Otherwise, this task expects that the LapH eigenvectors already         *
// *   reside in HostGlobal; otherwise, an error results.                                  *
// *                                                                                       *
// *   The default verbosity from the main program is used, unless overridden by           *
// *   the <Verbosity> tag above.                                                          *
// *                                                                                       *
// *   The following tags relate to control of the computations:                           *
// *      <UseMultiSrcInverter>  -- yes or no (yes can speed up, uses more memory)         *
// *      <NumSinksBeforeProject>  -- number of inversions to do before projecting         *
// *                                  onto Laph evs (the number to do in multi-src)        *
// *      <NumSinksInProjectBatch> -- while projecting onto Laph evs, the number of        *
// *                                  quark sinks to project at a time as one batch        *
// *                                  (generally should be same as value for previous tag) *
// *      <NumEigsInProjectBatch> -- number of Laph evs to project at a time as one batch  *
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
 bool report_gflops=false;
 if (xml_tag_count(xmltask,"ReportGflops")>0){
    report_gflops=true;}
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
 Q.computePerambulators(extra_soln_check,print_coeffs,report_gflops);
 outer.stop();
 printLaph(make_strf("LAPH_PERAMBULATORS: total time = %g secs",
           outer.getTimeInSeconds()));
 printLaph("LAPH_PERAMBULATORS: ran successfully\n"); 
} 


// *****************************************************************************************
// *                                                                                       *
// *   Task to check the LapH quark perambulators previously computed and written to       *
// *   file.  The input XML is expected to have the following  form:                       *
// *                                                                                       *
// *      <Task>                                                                           *
// *         <Name>LAPH_CHECK_PERAMBULATORS</Name>                                         *
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
// *            <SmearedQuarkFileStub>/path/sq_field</SmearedQuarkFileStub>                *
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
// *         </QuarkPerambulatorInfo>                                                      *
// *         <CheckSet>                                                                    *
// *            <Check>                                                                    *
// *               <SourceTime>8</SourceTime>                                              *
// *               <SourceLaphEigvecIndices>12 17 23 24 25</SourceLaphEigvecIndices>       *
// *            </Check>                                                                   *
// *            <Check>                                                                    *
// *               <SourceTime>12</SourceTime>                                             *
// *               <SourceLaphEigvecIndexMin>24</SourceLaphEigvecIndexMin>                 *
// *               <SourceLaphEigvecIndexMax>32</SourceLaphEigvecIndexMax>                 *
// *            </Check>                                                                   *
// *         </CheckSet>                                                                   *
// *         <LogFileStub>/path/logfilestub</LogFileStub>                                  *
// *         <VerboseOutput/> (optional)                                                   *
// *      </Task>                                                                          *
// *                                                                                       *
// *****************************************************************************************

void doLaphCheckPerambulators(XMLHandler& xmltask)
{
 if (xml_tag_count(xmltask,"QuarkPerambulatorInfo")!=1){
    errorLaph("Must have one <QuarkPerambulatorInfo> tag");}
 XMLHandler xmlr(xmltask,"QuarkPerambulatorInfo");
 GaugeConfigurationInfo gaugeinfo(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 QuarkActionInfo quark(xmlr);
 FileListInfo files(xmlr);
 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>0){
    upper_spin_only=true;}
 bool verbose_output=false;
 if (xml_tag_count(xmltask,"VerboseOutput")>0){
    verbose_output=true;}
 string logfilestub("perambs_check");
 xmlreadif(xmltask,"LogFileStub",logfilestub,"LAPH_CHECK_PERAMBULATORS");
 PerambulatorHandler::Mode mode=PerambulatorHandler::Check;

 printLaph("\n");
 printLaph(" ***********************************************************");
 printLaph(" *                                                         *");
 printLaph(" *   Laph Task: Check the quark perambulators and          *");
 printLaph(" *              output results to log files                *");
 printLaph(" *                                                         *");
 printLaph(" ***********************************************************\n");
 printLaph(make_strf("\n%s\n",gaugeinfo.output()));
 printLaph(make_strf("\n\nGluon Smearing:\n%s\n",gSmear.output()));
 printLaph(make_strf("\n\nQuark Smearing:\n%s\n",qSmear.output()));
 printLaph(make_str("\nQuarkAction:\n",quark.output()));
 if (upper_spin_only){
    printLaph("Only upper spin components used");}
 else{
    printLaph("All spin components used");}
 printLaph(make_strf("LogFileStub: %s",logfilestub));

     // create handler
 PerambulatorHandler Q(gaugeinfo,gSmear,qSmear,quark,files,
                       "",upper_spin_only,mode);

    // read the set of computations (time sources, src eigvec indices)
    // as well as the batching parameters
 XMLHandler xmlchk(xmltask,"CheckSet");
 Q.setChecks(xmlchk);

     // now do the checks!
 StopWatch outer; outer.start();
 Q.doChecks(logfilestub,verbose_output);
 outer.stop();
 printLaph(make_strf("LAPH_CHECK_PERAMBULATORS: total time = %g secs",
           outer.getTimeInSeconds()));
 printLaph("LAPH_CHECK_PERAMBULATORS: ran successfully\n"); 
} 


// *****************************************************************************************
// *                                                                                       *
// *   Task to merge the LapH quark perambulators previously computed and written to       *
// *   file.  The input XML is expected to have the following  form:                       *
// *                                                                                       *
// *      <Task>                                                                           *
// *         <Name>LAPH_MERGE_PERAMBULATORS</Name>                                         *
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
// *            <SmearedQuarkFileStub>/path/sq_field</SmearedQuarkFileStub>                *
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
// *         </QuarkPerambulatorInfo>                                                      *
// *         <InputFileListInfos>                                                          *
// *            <FileListInfo>/path/quark_perambs_1</FileNameStub>                         *
// *               <MinFileNumber>0</MinFileNumber>                                        *
// *               <MaxFileNumber>200</MaxFileNumber>                                      *
// *            </FileListInfo>                                                            *
// *            <FileListInfo>/path/quark_perambs_2</FileNameStub>                         *
// *               <MinFileNumber>0</MinFileNumber>                                        *
// *               <MaxFileNumber>200</MaxFileNumber>                                      *
// *            </FileListInfo>                                                            *
// *                    ...                                                                *
// *         </InputFileListInfos>                                                         *
// *      </Task>                                                                          *
// *                                                                                       *
// *****************************************************************************************

void doLaphMergePerambulators(XMLHandler& xmltask)
{
 if (xml_tag_count(xmltask,"QuarkPerambulatorInfo")!=1){
    errorLaph("Must have one <QuarkPerambulatorInfo> tag");}
 XMLHandler xmlr(xmltask,"QuarkPerambulatorInfo");
 GaugeConfigurationInfo gaugeinfo(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 QuarkActionInfo quark(xmlr);
 FileListInfo files(xmlr);
 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>0){
    upper_spin_only=true;}
 PerambulatorHandler::Mode mode=PerambulatorHandler::Merge;

 printLaph("\n");
 printLaph(" ***********************************************************");
 printLaph(" *                                                         *");
 printLaph(" *   Laph Task: Merge the quark perambulators              *");
 printLaph(" *                                                         *");
 printLaph(" ***********************************************************\n");
 printLaph(make_strf("\n%s\n",gaugeinfo.output()));
 printLaph(make_strf("\n\nGluon Smearing:\n%s\n",gSmear.output()));
 printLaph(make_strf("\n\nQuark Smearing:\n%s\n",qSmear.output()));
 printLaph(make_str("\nQuarkAction:\n",quark.output()));
 if (upper_spin_only){
    printLaph("Only upper spin components used");}
 else{
    printLaph("All spin components used");}

 if (xml_tag_count(xmltask,"InputFileListInfos")!=1){
    errorLaph("Must have one <InputFileListInfos> tag");}
 XMLHandler xmlinf(xmltask,"InputFileListInfos");
 list<XMLHandler> infxmls=xmlinf.find_among_children("FileListInfo");
 list<FileListInfo> inflos;
 for (list<XMLHandler>::iterator it=infxmls.begin();it!=infxmls.end();++it){
    inflos.push_back(FileListInfo(*it));}

     // create handler
 PerambulatorHandler Q(gaugeinfo,gSmear,qSmear,quark,files,
                       "",upper_spin_only,mode);

    // read the set of computations (time sources, src eigvec indices)
    // as well as the batching parameters
 XMLHandler xmlchk(xmltask,"CheckSet");
 Q.setChecks(xmlchk);

     // now do the checks!
 StopWatch outer; outer.start();
 for (list<FileListInfo>::iterator it=inflos.begin();it!=inflos.end();++it){
    Q.mergeData(*it);}
 outer.stop();
 printLaph(make_strf("LAPH_MERGE_PERAMBULATORS: total time = %g secs",
           outer.getTimeInSeconds()));
 printLaph("LAPH_MERGE_PERAMBULATORS: ran successfully\n"); 
} 


// *****************************************************************************
}
