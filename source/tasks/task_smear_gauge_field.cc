#include "field_smearing_info.h"
#include "gauge_configuration_handler.h"
#include "gauge_configuration_info.h"
#include "gluon_smearing_handler.h"
#include "latt_field.h"
#include "layout_info.h"
#include "stop_watch.h"
#include "tasks.h"

using namespace std;
using namespace quda;

namespace LaphEnv {

// *********************************************************************************
// * *
// *   Task to compute the smeared gauge field with stout link smearing and *
// *   either write to file or to the NamedObjMap.  The input XML is expected to
// *
// *   have the following form: *
// * *
// *   <Task> *
// *      <Name>SMEAR_GAUGE_FIELD</Name> *
// *      <GaugeConfigurationInfo> *
// *         <EnsembleName>CLS_A653</EnsembleName> *
// *         <FileFormat>CERN</FileFormat> *
// *         <ConfigType>WilsonImproved</ConfigType> *
// *         <FileName>/path/cls21_A653_r000/cfgs/A653r000n1</FileName> *
// *         <ConfigNumber>1</ConfigNumber> *
// *      </GaugeConfigurationInfo> *
// *      <GluonStoutSmearingInfo> *
// *         <LinkIterations>20</LinkIterations> *
// *         <LinkStapleWeight>0.1</LinkStapleWeight> *
// *      </GluonStoutSmearingInfo> *
// *      <SmearedGaugeFileName>NOM_smeared_gauge_field</SmearedGaugeFileName> *
// *   </Task> *
// * *
// *   The smeared gauge remains in the gpu device memory after this task. *
// * *
// *********************************************************************************

void doSmearGaugeField(XMLHandler &xmltask) {
  GaugeConfigurationInfo gaugeinfo(xmltask);
  GluonSmearingInfo smear(xmltask);
  string smeared_filename;
  xmlread(xmltask, "SmearedGaugeFileName", smeared_filename,
          "SMEAR_GAUGE_FIELD");

  printLaph("\n");
  printLaph(" ***********************************************************");
  printLaph(" *                                                         *");
  printLaph(" *   Laph Task:   Stout smear the gauge field and write    *");
  printLaph(" *                to file or the NamedObjMap               *");
  printLaph(" *                                                         *");
  printLaph(" ***********************************************************\n\n");
  printLaph(make_strf("%s", gaugeinfo.output()));
  printLaph(make_strf("\n%s\n", smear.output()));

  // create the handler in write mode

  GaugeConfigurationHandler GH(gaugeinfo);
  GH.setData();
  XMLHandler gauge_xmlinfo;
  GH.getXMLInfo(gauge_xmlinfo);
  printLaph("XML info for the gauge configuration:");
  printLaph(make_strf("%s\n", gauge_xmlinfo.output()));

  GH.copyDataToDevice();
  GluonSmearingHandler G(smear, gaugeinfo, smeared_filename, false);

  StopWatch outer;
  outer.start();

  //  Compute the stout-smeared gauge field, write to file or NamedObjMap.
  //  Leaves the smeared gauge on the gpu device.

  G.computeSmearedGaugeField();

  outer.stop();
  printLaph(make_strf("SMEAR_GAUGE_FIELD task: total time = %g seconds",
                      outer.getTimeInSeconds()));
  printLaph("SMEAR_GAUGE_FIELD task: ran successfully\n");
}
} // namespace LaphEnv
