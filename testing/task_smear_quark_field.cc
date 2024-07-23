#include "tasks.h"
#include "gauge_configuration_info.h"
#include "gauge_configuration_handler.h"
#include "quark_smearing_handler.h"
#include "layout_info.h"
#include "latt_field.h"
#include "field_smearing_info.h"
#include "stop_watch.h"
#include "laph_eigen_info.h"

using namespace std;
using namespace quda;

namespace LaphEnv {


// *********************************************************************
	
     // Subroutine to smear the quark field with LapH smearing


void doSmearQuarkField(XMLHandler& xmltask)
{
 GaugeConfigurationInfo gaugeinfo(xmltask);
 GluonSmearingInfo gsmear(xmltask);
 string smeared_gauge_filename;
 xmlread(xmltask,"SmearedGaugeFileName",smeared_gauge_filename,"SMEAR_QUARK_FIELD");
 QuarkSmearingInfo qsmear(xmltask);
 string smeared_quark_filestub;
 xmlread(xmltask,"SmearedQuarkFileStub",smeared_quark_filestub,"SMEAR_QUARK_FIELD");

 bool laph_solve=true;
 LaphEigenSolverInfo *eigsolveinfo=0;
 if (xml_tag_count(xmltask,"LaphEigenSolverInfo")==1)
    eigsolveinfo = new LaphEigenSolverInfo(xmltask);
 else 
    laph_solve=false;

// int striping_factor=1,striping_unit=0;
// xmlreadif(xml_rdr,"StripingFactor",striping_factor,"SMEAR_QUARK_FIELD_TIMESLICES");
// xmlreadif(xml_rdr,"StripingUnit",striping_unit,"SMEAR_QUARK_FIELD_TIMESLICES");
/*
 int mintime=0, maxtime=-1;    // optional... -1 maxtime means go to Textent-1
 xmlreadif(xml_rdr,"CalcMinTime",mintime,"SMEAR_QUARK_FIELD_TIMESLICES");
 xmlreadif(xml_rdr,"CalcMaxTime",maxtime,"SMEAR_QUARK_FIELD_TIMESLICES");
*/
 printLaph("\n");
 printLaph(" *************************************************************");
 printLaph(" *                                                           *");
 printLaph(" *   Laph Task: Smear the quark field by computing the       *");
 printLaph(" *              eigenvectors of the smeared-covariant        *");
 printLaph(" *              Laplacian, write to file or the NamedObjMap  *");
 printLaph(" *                                                           *");
 printLaph(" *************************************************************\n\n");
 printLaph(make_strf("%s",gaugeinfo.output()));
 printLaph(make_strf("\n%s",gsmear.output()));
 printLaph(make_strf("Smeared gauge field file name = %s",smeared_gauge_filename));
 printLaph(make_strf("\n%s",qsmear.output()));
 printLaph(make_strf("Smeared quark field file stub = %s",smeared_quark_filestub));
// QDPIO::cout << "Striping factor = "<<striping_factor<<endl;
// QDPIO::cout << "Striping unit = "<<striping_unit<<endl;
// QDPIO::cout << "Min time = "<<mintime<<endl;
// if (maxtime>=0) QDPIO::cout << "Max time = "<<maxtime<<endl;
// else QDPIO::cout << "Max time = NT-1"<<endl;
 if (laph_solve){
    printLaph(make_strf("\n%s",eigsolveinfo->output()));}
 else{
    printLaph("\nEstimate largest eigenvalue of -Laplacian only");}

    // create the handler (in write mode)

 QuarkSmearingHandler Q(gsmear,gaugeinfo,qsmear,smeared_quark_filestub,false);

 StopWatch outer; outer.start();

 if (laph_solve) Q.computeLaphEigenvectors(*eigsolveinfo,smeared_gauge_filename);
/* else{
    double lambda_max=Q.estimateLargestLaplacianEigenvalue(smeared_gauge_filename);
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
}
