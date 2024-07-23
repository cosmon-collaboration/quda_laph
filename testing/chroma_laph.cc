#include "chroma.h"
#include "inline_erase_named_obj.h"
#include "inline_smear_gauge_field.h"
#include "inline_smear_quark_field.h"
#include "inline_quark_line_ends.h"
#include "inline_quark_line_check.h"
#include "inline_quark_perambulators.h"
#include "inline_quark_line_undilute.h"
#include "inline_baryon_line_ends.h"
#include "inline_meson_line_ends.h"
#include "inline_current_line_ends.h"
#include "inline_meson_int_loops.h"
#include "inline_glueball_line_ends.h"
#include "inline_static_light.h"
#include "inline_tetraquark_line_ends.h"
#include "inline_tetraquark_one_int_loop.h"
#include "inline_tetraquark_two_int_loop.h"
#include "inline_db_output.h"
#include "inline_hdf5_gauge_field.h"
#include "inline_scale.h"
#include "inline_laph_smear_tune.h"
#include "inline_stout_smear_tune.h"
#include "inline_laph_smear_spatial_dist.h"
#include "inline_laph_dilution_tune.h"
#include "xml_help.h"
#include "meas/smear/hyp_link_smearing.h"

using namespace std;
using namespace Chroma;

struct Params_t
{
  multi1d<int>    nrow;
  std::string     inline_measurement_xml;
};

struct Inline_input_t
{
  Params_t        param;
  GroupXML_t      cfg;
  QDP::Seed       rng_seed;
  bool            config_read;
};

void read(XMLReader& xml, const std::string& path, Params_t& p) 
{
  XMLReader paramtop(xml, path);
  read(paramtop, "nrow", p.nrow);

  XMLReader measurements_xml(paramtop, "InlineMeasurements");
  std::ostringstream inline_os;
  measurements_xml.print(inline_os);
  p.inline_measurement_xml = inline_os.str();
  QDPIO::cout << "InlineMeasurements are: " << endl;
  QDPIO::cout << p.inline_measurement_xml << endl;
}

void read(XMLReader& xml, const std::string& path, Inline_input_t& p) 
{
 try{
    XMLReader paramtop(xml, path);      
    read(paramtop, "Param", p.param);
#if (QDP_ND == 4)
    if (LaphEnv::xml_tag_count(paramtop,"Cfg")>0){
       p.config_read=true; 
       p.cfg = readXMLGroup(paramtop, "Cfg", "cfg_type");}
    else{
       p.config_read=false;}
#endif
    //read(paramtop, "RNG", p.rng_seed);
 }
  catch( const std::exception& e ){
    QDPIO::cerr << "Error reading XML : " << e.what() << endl;
    QDP_abort(1);}
}

bool linkageHack(void)
{
  bool success = true;

  // Inline Measurements
  success &= InlineAggregateEnv::registerAll();
  success &= GaugeInitEnv::registerAll();
#if (QDP_ND == 4)
  success &= InlineSmearGaugeFieldEnv::registerAll();
  success &= InlineStochLaphQuarkEnv::registerAll();
  success &= InlineStochLaphStaticLightEnv::registerAll();
  success &= InlineQuarkPerambulatorsEnv::registerAll();
  success &= InlineStoutSmearTuneEnv::registerAll();
#ifdef LAPH_USE_HDF5
  success &= InlineHDF5GaugeFieldConvertEnv::registerAll();
  success &= InlineHDF5GaugeFieldReadEnv::registerAll();
#endif
#elif (QDP_ND == 3)
  success &= InlineStochLaphBaryonEnv::registerAll();
  success &= InlineStochLaphMesonEnv::registerAll();
  success &= InlineStochLaphMesonIntLoopEnv::registerAll();
  success &= InlineStochCheckQuarkEnv::registerAll();
  success &= InlineStochUndiluteQuarkEnv::registerAll();
  success &= InlineStochLaphDBOutputEnv::registerAll();
  success &= InlineScaleLaphEnv::registerAll();
  success &= InlineLaphSmearTuneEnv::registerAll();
  success &= InlineStochLaphGlueballEnv::registerAll();
  success &= InlineStochLaphTetraquarkEnv::registerAll();
  success &= InlineStochLaphTetraquarkOneIntLoopEnv::registerAll();
  success &= InlineStochLaphTetraquarkTwoIntLoopEnv::registerAll();
  success &= InlineLaphSmearSpatialEnv::registerAll();
  success &= InlineLaphDilutionTuneEnv::registerAll();
#endif
  success &= InlineStochCheckQuarkEnv::registerAll();
  success &= InlineSmearQuarkFieldEnv::registerAll();
  success &= InlineEraseObjectEnv::registerAll();
  success &= HypLinkSmearingEnv::registerAll();
#if (QDP_ND == 4)
#if QDP_USE_LEXICO_LAYOUT == 1
  success &= InlineStochLaphCurrentEnv::registerAll();
#endif
#endif


  return success;
}



// *************************************************
// *                                               *
// *   Main program to run all measurement codes   *
// *                                               *
// *************************************************


int main(int argc, char *argv[]) 
{
  // Chroma Init stuff
  Chroma::initialize(&argc, &argv);
//  QDP::RNG::lattice_ran_mult=0;
  START_CODE();

#if (QDP_ND == 4)
  QDPIO::cout << endl << "Starting ChromaLaph in 4 space-time dimensions"<<endl<<endl;
#elif (QDP_ND == 3)
  QDPIO::cout << endl << "Starting ChromaLaph in 3 space-time dimensions"<<endl<<endl;
#endif
  linkageHack();

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  XMLReader xml_in;

  // Input parameter structure
  Inline_input_t  input;
  try{
    xml_in.open(Chroma::getXMLInputFileName());
    read(xml_in, "/chroma", input);}
  catch(const std::exception& e){
    QDPIO::cerr << "CHROMA_LAPH: Caught Exception reading XML: " << e.what() << endl;
    QDP_abort(1);}
  catch(std::exception& e){
    QDPIO::cerr << "CHROMA_LAPH: Caught standard library exception: " << e.what() << endl;
    QDP_abort(1);}
  catch(const std::string& e){
    QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
    QDP_abort(1);}
  catch(...){
    QDPIO::cerr << "CHROMA: caught generic exception reading xml" << std::endl;
    QDPIO::cerr << "Rethrowing" << std::endl;
    throw;
    QDP_abort(1);}

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "chroma");

  Layout::setLattSize(input.param.nrow);
  Layout::create();

  proginfo(xml_out);    // Print out basic program info
#if (ARCH_PARSCALAR == 1)
  QDPIO::cout << "Number of MPI ranks = "<<Layout::numNodes()<<endl<<endl;
#endif

  // Initialise the RNG
  //QDP::RNG::setrn(input.rng_seed);
  //write(xml_out,"RNG", input.rng_seed);
  StopWatch swatch;

  // Write out the input
  write(xml_out, "Input", xml_in);

#if (QDP_ND == 4)

  if (input.config_read){

  // Start up the config
  swatch.reset();
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Start up the gauge field
  QDPIO::cout << "Attempt to read gauge field" << endl;
  swatch.start();
  try{
    std::istringstream  xml_c(input.cfg.xml);
    XMLReader  cfgtop(xml_c);
    QDPIO::cout << "Gauge initialization: cfg_type = " << input.cfg.id << endl;
    Handle< GaugeInit >
      gaugeInit(TheGaugeInitFactory::Instance().createObject(
                input.cfg.id,cfgtop,input.cfg.path));
    (*gaugeInit)(gauge_file_xml, gauge_xml, u);}
  catch(std::bad_cast){
    QDPIO::cerr << "CHROMA_LAPH: caught cast error" << endl;
    QDP_abort(1);}
  catch(std::bad_alloc){ 
    // This might happen on any node, so report it
    cerr << "CHROMA_LAPH: caught bad memory allocation" << endl;
    QDP_abort(1);}
  catch(const std::exception& e){
    QDPIO::cerr << "CHROMA_LAPH: Caught Exception: " << e.what() << endl;
    QDP_abort(1);}
  catch(...){
    // This might happen on any node, so report it
    cerr << "CHROMA_LAPH: caught generic exception during gaugeInit" << endl;
    QDPIO::cerr << "Rethrowing" << std::endl;
    throw;
    QDP_abort(1);}
  swatch.stop();

  QDPIO::cout << "Gauge field successfully read: time= " 
              << swatch.getTimeInSeconds()<< " secs" << endl;

  XMLBufferWriter config_xml;
  config_xml << gauge_xml;

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

    // Reset and set the default gauge field
  InlineDefaultGaugeField::reset();
  InlineDefaultGaugeField::set(u, config_xml);

  }

#endif
  
  // Get the measurements
  try 
  {
    std::istringstream Measurements_is(input.param.inline_measurement_xml);
    XMLReader MeasXML(Measurements_is);
    multi1d < Handle< AbsInlineMeasurement > > the_measurements;
    read(MeasXML, "/InlineMeasurements", the_measurements);

    QDPIO::cout << "There are " << the_measurements.size() << " measurements " << endl;

    // Measure inline observables 
    push(xml_out, "InlineObservables");
    xml_out.flush();

    QDPIO::cout << "Doing " << the_measurements.size() 
                <<" measurements" << endl;
    swatch.start();
    unsigned long cur_update = 0;
    for (int m=0; m < the_measurements.size(); m++){
       AbsInlineMeasurement& the_meas = *(the_measurements[m]);
       push(xml_out, "elem");
       the_meas(cur_update, xml_out);
       pop(xml_out); 
       xml_out.flush();}
    swatch.stop();

    QDPIO::cout << "CHROMA_LAPH measurements: time= " << swatch.getTimeInSeconds() << " secs" << endl;
    pop(xml_out); // pop("InlineObservables");

    // Reset the default gauge field
    InlineDefaultGaugeField::reset();}
  catch(std::bad_cast){
    QDPIO::cerr << "CHROMA_LAPH: caught cast error" << endl;
    QDP_abort(1);}
  catch(std::bad_alloc){ 
    // This might happen on any node, so report it
    cerr << "CHROMA_LAPH: caught bad memory allocation" << endl;
    QDP_abort(1);}
  catch(const std::exception& e){
    QDPIO::cerr << "CHROMA_LAPH: Caught Exception: " << e.what() << endl;
    QDP_abort(1);}
  catch(...){
    // This might happen on any node, so report it
    cerr << "CHROMA_LAPH: caught generic exception during measurement" << endl;
    QDPIO::cerr << "Rethrowing" << std::endl;
    throw;
    QDP_abort(1);}
  pop(xml_out);
  xml_out.flush();

  snoop.stop();
  QDPIO::cout << "CHROMA_LAPH: total time = " << snoop.getTimeInSeconds() << " secs" << endl;
  QDPIO::cout << "CHROMA_LAPH: ran successfully" << endl;

#if (QDP_ND == 4)

  // Trigger clean up of quark sink and openQCD handler explicitly
  if (TheNamedObjMap::Instance().check("OpenQCDHandler"))
    TheNamedObjMap::Instance().erase("OpenQCDHandler");

  if (TheNamedObjMap::Instance().check(LaphEnv::CurrentHandler::quarkSinkMapID))
    TheNamedObjMap::Instance().erase(LaphEnv::CurrentHandler::quarkSinkMapID);

#endif


  END_CODE();
  delete QDP::RNG::lattice_ran_mult;  
  Chroma::finalize();
  xmlCleanupParser();   // clean up libxml2 globals
 
}

