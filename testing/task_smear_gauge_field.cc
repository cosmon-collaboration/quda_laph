#include "tasks.h"
#include "gauge_configuration_info.h"
#include "gauge_configuration_handler.h"
#include "gluon_smearing_handler.h"
#include "layout_info.h"
#include "latt_field.h"
#include "field_smearing_info.h"
#include "stop_watch.h"

using namespace std;
using namespace quda;

namespace LaphEnv {


// *********************************************************************
	
     // Subroutine to smear the gauge field with stout link smearing


void doSmearGaugeField(XMLHandler& xmltask)
{
 GaugeConfigurationInfo gaugeinfo(xmltask);
 GluonSmearingInfo smear(xmltask);
 string smeared_filename;
 xmlread(xmltask,"SmearedGaugeFileName",smeared_filename,"SMEAR_GAUGE_FIELD");

 printLaph("\n");
 printLaph(" ***********************************************************");
 printLaph(" *                                                         *");
 printLaph(" *   Laph Task:   Stout smear the gauge field and write    *");
 printLaph(" *                to file or the NamedObjMap               *");
 printLaph(" *                                                         *");
 printLaph(" ***********************************************************\n\n");
 printLaph(make_strf("%s",gaugeinfo.output()));
 printLaph(make_strf("\n%s\n",smear.output()));

    // create the handler in write mode

 GaugeConfigurationHandler GH(gaugeinfo);
 GH.setData();
 XMLHandler gauge_xmlinfo;
 GH.getXMLInfo(gauge_xmlinfo);
 printLaph("XML info for the gauge configuration:");
 printLaph(make_strf("%s\n",gauge_xmlinfo.output()));

 GH.copyDataToDevice();
 GluonSmearingHandler G(smear,gaugeinfo,smeared_filename,false);

 StopWatch outer; outer.start();

    //  Compute the stout-smeared gauge field, write to file or NamedObjMap.
    //  Leaves the smeared gauge on the gpu device. 

 G.computeSmearedGaugeField();
 
 outer.stop();
 printLaph(make_strf("SMEAR_GAUGE_FIELD task: total time = %g seconds",
           outer.getTimeInSeconds())); 
 printLaph("SMEAR_GAUGE_FIELD task: ran successfully\n");
} 

// ******************************************************************
}
