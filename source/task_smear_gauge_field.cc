#include "tasks.h"
#include "gluon_smearing_handler.h"
#include "stop_watch.h"

using namespace std;

namespace LaphEnv {


// *********************************************************************************
// *                                                                               *
// *   Task to compute the smeared gauge field with stout link smearing and        *
// *   optionially write to file. The gauge field will be stored in HostGlobal.    *
// *   The input XML is expected to have the following form:                       *
// *                                                                               *
// *   <Task>                                                                      *
// *      <Name>SMEAR_GAUGE_FIELD</Name>                                           *
// *      <GaugeConfigurationInfo>                                                 *
// *         <EnsembleName>CLS_A653</EnsembleName>                                 *
// *         <FileFormat>CERN</FileFormat>                                         *
// *         <ConfigType>WilsonImproved</ConfigType>                               *
// *         <FileName>/path/cls21_A653_r000/cfgs/A653r000n1</FileName>            *
// *         <ConfigNumber>1</ConfigNumber>                                        *
// *      </GaugeConfigurationInfo>                                                *
// *      <GluonStoutSmearingInfo>                                                 *
// *         <LinkIterations>20</LinkIterations>                                   *
// *         <LinkStapleWeight>0.1</LinkStapleWeight>                              *
// *      </GluonStoutSmearingInfo>                                                *
// *      <SmearedGaugeFileName>smg_field</SmearedGaugeFileName> (see below)       *
// *   </Task>                                                                     *
// *                                                                               *
// *   If the tag <SmearedGaugeFileName> is set, then the smeared gauge field is   *
// *   written out to file.  In any case, the smeared gauge field is returned in   *
// *   HostGlobal.                                                                 *
// *                                                                               *
// *   The smeared gauge remains in the gpu device memory after this task.         * 
// *                                                                               *
// *********************************************************************************
	

void doSmearGaugeField(XMLHandler& xmltask)
{
 GaugeConfigurationInfo gaugeinfo(xmltask);
 GluonSmearingInfo smear(xmltask);
 string smeared_filename;
 xmlreadif(xmltask,"SmearedGaugeFileName",smeared_filename,"SMEAR_GAUGE_FIELD");

 printLaph("\n");
 printLaph(" ***********************************************************");
 printLaph(" *                                                         *");
 printLaph(" *   Laph Task:   Stout smear the gauge field and write    *");
 printLaph(" *                to file or the HostGlobal                *");
 printLaph(" *                                                         *");
 printLaph(" ***********************************************************\n\n");
 printLaph(make_strf("%s",gaugeinfo.output()));
 printLaph(make_strf("\n%s\n",smear.output()));

    // create the handler in write mode

 StopWatch outer; outer.start();
 GluonSmearingHandler G(smear,gaugeinfo,smeared_filename,false);

    //  Compute the stout-smeared gauge field, place in the HostGlobal,
    //  optionally write to file. Leaves the smeared gauge 
    //  on the gpu device. 

 G.computeSmearedGaugeField();

 outer.stop();
 printLaph(make_strf("SMEAR_GAUGE_FIELD task: total time = %g seconds",
           outer.getTimeInSeconds()));
 printLaph("SMEAR_GAUGE_FIELD task: ran successfully\n");
} 

// ******************************************************************
}
