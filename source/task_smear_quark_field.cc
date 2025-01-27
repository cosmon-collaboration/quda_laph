#include "tasks.h"
#include "quark_smearing_handler.h"
#include "stop_watch.h"

using namespace std;

namespace LaphEnv {


// ************************************************************************************
// *                                                                                  *
// *   Task to smear the quark field with LapH smearing and optionally write to       *
// *   file.  The eigenvectors will be stored in HostGlobal.  The input XML is        *
// *   expected to have the following form:                                           *
// *                                                                                  *
// *    <Task>                                                                        *
// *       <Name>SMEAR_QUARK_FIELD</Name>                                             *
// *       <GaugeConfigurationInfo>                                                   *
// *          <EnsembleName>CLS_A653</EnsembleName>                                   *
// *          <FileFormat>CERN</FileFormat>                                           *
// *          <ConfigType>WilsonImproved</ConfigType>                                 *
// *          <FileName>/path/cfgs/A653r000n1</FileName>                              *
// *          <ConfigNumber>1</ConfigNumber>                                          *
// *       </GaugeConfigurationInfo>                                                  *
// *       <GluonStoutSmearingInfo>                                                   *
// *          <LinkIterations>20</LinkIterations>                                     *
// *          <LinkStapleWeight>0.1</LinkStapleWeight>                                *
// *       </GluonStoutSmearingInfo>                                                  *
// *       <SmearedGaugeFileName>/path/smg_field</SmearedGaugeFileName> (see below)   *
// *       <QuarkLaphSmearingInfo>                                                    *
// *          <LaphSigmaCutoff>0.32</LaphSigmaCutoff>                                 *
// *          <NumberLaphEigvecs> 164 </NumberLaphEigvecs>                            *
// *       </QuarkLaphSmearingInfo>                                                   *
// *       <SmearedQuarkFileStub>/path/smq_field</SmearedQuarkFileStub> (see below)   *
// *       <LaphEigenSolverInfo>                                                      *
// *          <ResidualTolerance>1.0e-11</ResidualTolerance>                          *
// *          <MaxIterations>200</MaxIterations>                                      *
// *          <KrylovDimension>200</KrylovDimension>                                  *
// *          <ChebyshevOrder>6</ChebyshevOrder>                                      *
// *          <MaxEigenvalue>15.0</MaxEigenvalue>                                     *
// *          <CutoffEigenvalue>0.5</CutoffEigenvalue>                                *
// *          <StartingVectorType>equal_components</StartingVectorType>               *
// *          <CheckSolution>true</CheckSolution>                                     *
// *          <Verbosity>low</Verbosity> (optional)                                   *
// *       </LaphEigenSolverInfo>                                                     *
// *    </Task>                                                                       *
// *                                                                                  *
// *   If the tag <SmearedGaugeFileName> is set, then the smeared gauge field is      *
// *   read from file.  Otherwise, this routine expects that the smeared gauge        *
// *   field resides in HostGlobal; if not there, an error results.                   *
// *                                                                                  *
// *   If the tag <SmearedQuarkFileStub> is set, then the LapH eigenvectors are       *
// *   written out to file(s).  In any case, the eigenvectors are returned in         *
// *   HostGlobal.                                                                    *
// *                                                                                  *
// *   The LapH eigenvectors are removed from the gpu device memory since it is       *
// *   assumed quark propagators will subsequently be computed, so gpu memory         *
// *   will be needed for those computations.                                         *
// *                                                                                  *
// *   The default verbosity from the main program is used, unless overridden by      *
// *   the <Verbosity> tag above.                                                     *
// *                                                                                  *
// ************************************************************************************
	

void doSmearQuarkField(XMLHandler& xmltask)
{
 GaugeConfigurationInfo gaugeinfo(xmltask);
 GluonSmearingInfo gsmear(xmltask);
 string smeared_gauge_filename;
 xmlreadif(xmltask,"SmearedGaugeFileName",smeared_gauge_filename,"SMEAR_QUARK_FIELD");
 QuarkSmearingInfo qsmear(xmltask);
 string smeared_quark_filestub;
 xmlreadif(xmltask,"SmearedQuarkFileStub",smeared_quark_filestub,"SMEAR_QUARK_FIELD");
 LaphEigenSolverInfo eigsolveinfo(xmltask);

 printLaph("\n");
 printLaph(" *************************************************************");
 printLaph(" *                                                           *");
 printLaph(" *   Laph Task: Smear the quark field by computing the       *");
 printLaph(" *              eigenvectors of the smeared-covariant        *");
 printLaph(" *              Laplacian, write to file or the Host Global  *");
 printLaph(" *                                                           *");
 printLaph(" *************************************************************\n\n");
 printLaph(make_strf("%s",gaugeinfo.output()));
 printLaph(make_strf("\n%s",gsmear.output()));
 if (!smeared_gauge_filename.empty()){
    printLaph(make_strf("Smeared gauge field file name = %s",smeared_gauge_filename));}
 printLaph(make_strf("\n%s",qsmear.output()));
 if (!smeared_quark_filestub.empty()){
    printLaph(make_strf("Smeared quark field file stub = %s",smeared_quark_filestub));}
 printLaph(make_strf("\n%s",eigsolveinfo.output()));

    // create the handler (in write mode)

 QuarkSmearingHandler Q(gsmear,gaugeinfo,qsmear,smeared_quark_filestub,false);
 StopWatch outer; outer.start();
 Q.computeLaphEigenvectors(eigsolveinfo,smeared_gauge_filename);
 outer.stop();
 printLaph(make_strf("SMEAR_QUARK_FIELD task: total time = %g seconds",
                     outer.getTimeInSeconds()));
 printLaph("SMEARED_QUARK_FIELD task: ran successfully\n");
} 

// ******************************************************************
}
