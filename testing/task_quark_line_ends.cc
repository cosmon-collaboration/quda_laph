#include "tasks.h"
#include "gauge_configuration_info.h"
#include "quark_smearing_handler.h"
#include "dilution_scheme_info.h"
#include "quark_action_info.h"
#include "inverter_info.h"
#include "filelist_info.h"
#include "quark_handler.h"
#include "stop_watch.h"

using namespace std;
using namespace quda;

namespace LaphEnv {


// *********************************************************************
	
     // Subroutine to compute the LapH quark line ends


void doLaphQuarkLineEnds(XMLHandler& xmltask)
{
 if (xml_tag_count(xmltask,"QuarkLineEndInfo")!=1){
    errorLaph("Must have one <QuarkLineEndInfo> tag");}
 XMLHandler xmlr(xmltask,"QuarkLineEndInfo");
 GaugeConfigurationInfo gaugeinfo(xmlr);
 GluonSmearingInfo gSmear(xmlr);
 QuarkSmearingInfo qSmear(xmlr);
 string smeared_quark_filestub;
 xmlread(xmlr,"SmearedQuarkFileStub",smeared_quark_filestub,"LAPH_QUARK_LINE_ENDS");
 string smeared_gauge_filename;
 xmlread(xmlr,"SmearedGaugeFileName",smeared_gauge_filename,"LAPH_QUARK_LINE_ENDS");
 DilutionSchemeInfo dil(xmlr);
 QuarkActionInfo quark(xmlr);
 FileListInfo files(xmlr);
 InverterInfo invinfo(xmlr);
 bool verbose=false;
 string verbosity;
 xmlreadif(xmlr,"Verbosity",verbosity,"LAPH_QUARK_LINE_ENDS");
 if (tidyString(verbosity)=="full") verbose=true;

 printLaph("\n");
 printLaph(" ***********************************************************");
 printLaph(" *                                                         *");
 printLaph(" *   Laph Task: Compute the quark line ends                *");
 printLaph(" *              and write to file as time slices           *");
 printLaph(" *                                                         *");
 printLaph(" ***********************************************************\n");
 printLaph(make_strf("%",gaugeinfo.output()));
 printLaph(make_strf("\n\nGluon Smearing:\n%s\n",gSmear.output()));
 printLaph(make_strf("\n\nQuark Smearing:\n%s\n",qSmear.output()));
 printLaph(make_strf("SmearedQuarkFileStub: %s",smeared_quark_filestub));
 printLaph(make_strf("SmearedGaugeFileName: %s",smeared_gauge_filename));
 printLaph(make_strf("\nDilution Scheme Info:\n%s\n",dil.output()));
 printLaph(make_str("\nQuarkAction:\n",quark.output()));
 printLaph(make_str("\nInverter Info:\n",invinfo.output()));
 bool compute_mode=true;

     // create quark handler
 QuarkHandler Q(gaugeinfo,gSmear,qSmear,dil,quark,files,
                smeared_quark_filestub,smeared_gauge_filename,
                compute_mode);

    // read the list of computations (noises, time sources, file indices)
    // as well as the batching parameters
 XMLHandler xmlcmp(xmltask,"SinkComputations");
 Q.setSinkComputations(xmlcmp);

     // set the inverter info

 Q.setInverter(invinfo);
 printLaph("Info initialized in QuarkHandler");
 Q.outputSuffixMap();

     // now do the computations!
 StopWatch outer; outer.start();
 Q.computeSinks(verbose);
 outer.stop();
 printLaph(make_strf("LAPH_QUARK_LINE_ENDS: total time = %g secs",
           outer.getTimeInSeconds()));
 printLaph("LAPH_QUARK_LINE_ENDS: ran successfully\n"); 
} 

// ******************************************************************
}
