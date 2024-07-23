#include "xml_handler.h"
#include <cstdio>
#include <ctime>
#include <map>
#include "quda.h"
#include "utils.h"
#include "laph_stdio.h"
#include "stop_watch.h"
#include "layout_info.h"
#include "named_obj_map.h"
#include "quda_info.h"
#if defined(QUDA_OPENMP)
#include <omp.h>
#endif
#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif
#include "tasks.h"
#ifdef TESTING
#include "task_tests.h"
using namespace QLTestEnv;
#endif
using namespace std;
using namespace LaphEnv;


// **********************************************************************
// *                                                                    *
// *         QudaLaph:  main driver program to run all tasks            *
// *                                                                    *
// *   QudaLaph does stochastic LapH tasks using the QUDA library for   *
// *   lattice field theory.  QUDA should be compiled with MPI support  *
// *   when running in parallel mode.  Quda should also be compiled     *
// *   with the QDP interface turned ON.  Lattice quantities will be in *
// *   QDP format on the hosts (cpus).  The program requires the        *
// *   following command line arguments:                                *
// *                                                                    *
// *      (mpirun ...) quda_laph  -i input.xml -npartitions 2 2 2 3     *
// *                                                                    *
// *   The file "input.xml" must contain the appropriate input          *
// *   information (discussed below) in XML format.  The npartitions    *
// *   argument indicates how to split up the 4-dim space-time lattice  *
// *   into identical sublattices on each MPI rank.  The npartitions    *
// *   argument is ignored in serial mode.  If the npartitions argument *
// *   is absent, 1 1 1 1 is assumed (serial mode).  NOTE: currently    *
// *   it is required that there is one and only one device (gpu) for   *
// *   each MPI rank. Standard output is used.  Use redirection for     *
// *   output to a log file.                                            *
// *                                                                    *
// *   The environment variable QUDA_RESOURCE_PATH should be set to     *
// *   the directory (which should exist) where Quda will store its     *
// *   kernel tuning information.                                       *
// *                                                                    *
// *   The input file must contain a single XML document with root tag  *
// *   named <QudaLaph>. Inside the root tag should be one              *
// *   <LatticeLayoutInfo> tag, and one or more <Task> tags, and inside *
// *   each <Task> element should be a <Name> tag whose content is the  *
// *   name of the task.  The name must be one of the allowed names     *
// *   specified in the "do_task" subroutine.  If a tag <EchoXML/>      *
// *   is present as a child of the root tag, then the XML input        *
// *   is echoed to standard output. A <Verbosity> tag can only be      *
// *   included.                                                        *
// *                                                                    *
// *   Sample input XML:                                                *
// *                                                                    *
// *    <QudaLaph>                                                      *
// *       <LatticeLayoutInfo>                                          *
// *          ....                                                      *
// *       </LatticeLayoutInfo>                                         *
// *       <QudaInfo>    (optional)                                     *
// *          ....                                                      *
// *       </QudaInfo>                                                  *
// *       <EchoXML/>  (optional)                                       *
// *       <Verbosity>...</Verbosity>  (none, low, medium, high)        *
// *       <Task>                                                       *
// *          <Name>Task 1</Name>                                       *
// *              ....                                                  *
// *       </Task>                                                      *
// *       <Task>                                                       *
// *          <Name>Task 2</Name>                                       *
// *       </Task>                                                      *
// *    </QudaLaph>                                                     *
// *                                                                    *
// *                                                                    *
// *   The <LatticeLayoutInfo> tag specifies the lattice size.          *
// *                                                                    *
// *       <LatticeLayoutInfo>                                          *
// *         <XYZTSizes>24 24 24 96</XYZTSizes>                         *
// *       </LatticeLayoutInfo>                                         *
// *                                                                    *
// *   On the hosts (cpus), QDP format of sites, color, and spin is     *
// *   always assumed, and even-odd checkerboard with site ordering     *
// *   xyzt with x varying fastest is always assumed.  QUDA converts    *
// *   the CPU ordering to its preferred ordering on the devices        *
// *   (gpus) which is done very quickly, according to Kate.            *
// *                                                                    *
// *   The default precisions are set in the <QudaInfo> tag:            *
// *                                                                    *
// *   <QudaInfo>                                                       *
// *      <CPUPrecision>double</CPUPrecision> (or single)               *
// *      <CUDAPrecision>double</CUDAPrecision> (or single)             *
// *      <CUDASloppyPrecision>double</CUDASloppyPrecision> (or single) *
// *   </QudaInfo>                                                      *
// *                                                                    *
// **********************************************************************

  // Parse the command line options:
  // One mandatory:  -i <input.xml>   -> return in "inputxmlfile"
  // Communications partitioning  -npartitions xn yn zn tn  -> return in "npartitions"

void parse_args(int* argc, char ***argv, std::vector<int>& npartitions, std::string& inputxmlfile)
{
 npartitions.clear();
 inputxmlfile.clear();
 int nargs=*argc;
 int Nd=4;      // four dimensions required

 for (int i=0; i<nargs; ++i){
    string argv_i( (*argv)[i] );
      //  get the input xml file name
    if ((argv_i==std::string("-i"))&&((i+1)<nargs)){
       inputxmlfile=std::string( (*argv)[i+1] );
       ++i;}
    else if ((argv_i=="-npartitions")&&((i+Nd)<nargs)){
       npartitions.resize(Nd); 
       for (int j=0;j<Nd;j++){
          int uu;
          sscanf((*argv)[++i], "%d", &uu);
          if (uu<=0){
             npartitions.clear();
             return;}
          npartitions[j] = uu;}}}
}


void output_datetime()
{
 time_t rawtime;
 struct tm *timeinfo;
 time(&rawtime);
 timeinfo=localtime(&rawtime);
 printLaph(make_strf("  Current date/time: %s\n",asctime(timeinfo)));
}


  //  Constructor sets up the known tasks.  Call
  //  member function "do_task" to perform the task.
  //  "TaskMap" is a map that associates a string
  //  containing a task name to a function pointer.

class Tasker
{
   typedef void (*task_ptr)(XMLHandler&);
   std::map<string, task_ptr> TaskMap;

 public:
    
   Tasker();
   ~Tasker(){}
   void do_task(XMLHandler& xml_in, bool echo);

};

     // set up the known tasks

Tasker::Tasker()
{
 TaskMap["SMEAR_GAUGE_FIELD"]=&doSmearGaugeField;
 TaskMap["SMEAR_QUARK_FIELD"]=&doSmearQuarkField;
 TaskMap["LAPH_QUARK_LINE_ENDS"]=&doLaphQuarkLineEnds;
#ifdef TESTING
 TaskMap["QUDA_LAPH_TEST"]=&QLTestEnv::doQudaLaphTests;
#endif
};


void Tasker::do_task(XMLHandler& xml_task, bool echo)
{
#ifdef ARCH_PARALLEL
 comm_barrier();
#endif
 if (echo){
    printLaph("Input XML for this task:");
    printLaph(xml_task.output());}
 XMLHandler xmlt(xml_task);
 xml_child_assert(xmlt,"Name","do_task");
 if (!xmlt.is_simple_element()){
    throw(std::runtime_error("<Name> tag absent or is not simple XML element"));}
 string task_name=xmlt.get_text_content();
 printLaph(make_strf("  Task name = %s\n",task_name));

 map<string,task_ptr >::iterator taskit=TaskMap.find(task_name);
 if (taskit!=TaskMap.end()){
    (*(taskit->second))(xml_task);}  // do the task!!
 else{
    throw(std::runtime_error("Unknown task name"));}   // unknown task?
}


void initRand()
{
  int rank = 0;
#if defined(ARCH_PARALLEL)
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  srand(17 * rank + 137);
}


void readVerbosity(XMLHandler& xmlin, QudaVerbosity_s& verbosity)
{
 string tagvalue;
 xmlread(xmlin,"Verbosity",tagvalue,"QudaLaph");
 if (tagvalue=="none"){
    verbosity=QUDA_SILENT;
    printLaph("Verbosity set to none: SILENT\n");}
 else if (tagvalue=="low"){
    verbosity=QUDA_SUMMARIZE;
    printLaph("Verbosity set to medium: SUMMARIZE\n");}
 else if ((tagvalue=="medium")||(tagvalue=="med")){
    verbosity=QUDA_VERBOSE;
    printLaph("Verbosity set to medium: VERBOSE\n");}
 else if (tagvalue=="high"){
    verbosity=QUDA_DEBUG_VERBOSE;
    printLaph("Verbosity set to high: DEBUG_VERBOSE\n");}
}


std::map<std::string, NamedObjBase*> NamedObjMap::the_map;



// ****************************************************************
// *                                                              *
// *                           Main                               *
// *                                                              *
// ****************************************************************


int main(int argc, char** argv)
{
   // parse the command line options -i and -npartitions 
 vector<int> npartitions;
 string inputxmlfile;
 parse_args(&argc,&argv,npartitions,inputxmlfile);

#ifdef ARCH_PARALLEL
 MPI_Init(&argc,&argv);
 if ((int(npartitions.size())!=LayoutInfo::Ndim)||(inputxmlfile.empty())){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank==0){
       std::cout <<endl<<endl<< "Invalid command line options: should be "
                 << "-i <inputxmlfile> -npartitions <nx> <ny> <nz> <nt>"<<endl<<endl;}
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);}
 initCommsGridQuda(LayoutInfo::Ndim,npartitions.data(),nullptr, nullptr);
 comm_barrier();
#else
 if (inputxmlfile.empty()){
    std::cout <<endl<<endl<< "Invalid command line options: should be "
              << "-i <inputxmlfile> -npartitions <nx> <ny> <nz> <nt>"<<endl<<endl;
    std::exit(1);}
#endif

 fflush(stdout);
 initRand();  
 printLaph("\nStarting QUDA_LAPH");
 output_datetime();
  
 StopWatch rolex;
 rolex.start();
 printLaph(make_strf("  Name of input XML file: %s\n\n",inputxmlfile));

 XMLHandler xml_in;
 xml_in.set_from_file(inputxmlfile);
 if (xml_in.fail()){
    errorLaph("  Unable to read/parse XML content in input file ... exiting...");}
 if (xml_in.get_node_name()!="QudaLaph"){
    errorLaph("  Root tag of input XML must be named \"QudaLaph\" ... exiting...");}

 char *qrpath = getenv("QUDA_RESOURCE_PATH");
 if (!qrpath) {
    printLaph("Environment variable QUDA_RESOURCE_PATH is not set.");
    printLaph("Ensure it is set to allow caching of tuned Quda parameters.");}
 else{
    printLaph(make_strf("QUDA_RESOURCE_PATH is set to %s\n",qrpath));}

 bool echo=false;
 bool layoutinfo=false;
 bool qudainfo=false;
 QudaVerbosity_s verbosity=QUDA_SILENT;
 XMLHandler xml_layoutinfo;
 XMLHandler xml_qudainfo;
 int ntasks = 0;
 try{
    xml_in.seek_first_child();
    while (xml_in.good()){
       if (xml_in.get_node_name()=="Task") ntasks++;
       else if (xml_in.get_node_name()=="EchoXML") echo=true;
       else if (xml_in.get_node_name()=="LatticeLayoutInfo"){
          layoutinfo=true;
          xml_layoutinfo.set(xml_in);}
       else if (xml_in.get_node_name()=="QudaInfo"){
          qudainfo=true;
          xml_qudainfo.set(xml_in);}
       else if (xml_in.get_node_name()=="Verbosity"){
          readVerbosity(xml_in,verbosity);}
       else throw(std::runtime_error("Invalid XML input"));
       xml_in.seek_next_sibling();}}
 catch(const std::exception& err){
    errorLaph(make_strf("  Error: %s ... exiting...\n",err.what()));}

 if (!layoutinfo){
    errorLaph("No <LatticeLayoutInfo> tag was found as child of root tag");}
 try{
    LayoutInfo::init(xml_layoutinfo,npartitions,echo);}
 catch(const std::exception& err){
    errorLaph("Problem reading contents of <LatticeLayoutInfo>");}
      // if non-default quda parameters requested, set them
 if (qudainfo){
    QudaInfo::init(xml_qudainfo,echo);}
     // set the verbosity level for Quda output
 setVerbosity(verbosity); 

   // Quda initializations
 //initComms(argc, argv, gridsize_from_cmdline);
 initQuda(QudaInfo::getDeviceOrdinal());

 printLaph(make_strf("\n\n  Number of tasks is %d",ntasks));
#if defined(QUDA_OPENMP)
 printLaph(make_strf("  Maximum number of threads is %d",omp_get_max_threads));
#endif
#ifdef ARCH_PARALLEL
 printLaph(make_strf("  Number of MPI ranks = %d",comm_size()));
#endif

 xml_in.seek_root();
 xml_in.seek_first_child();
 Tasker T;

 for (int task=1;task<=ntasks;++task){
    StopWatch swatch;
    swatch.start();
    try{
       printLaph(make_strf("\n\nStarting Task %d\n",task));
       if (xml_in.fail()) throw(std::runtime_error("XML input error"));
       while (xml_in.get_node_name()!="Task"){
          xml_in.seek_next_sibling();}
       XMLHandler xml_task(xml_in);
           // the main task
       T.do_task(xml_task,echo);
       }
    catch(const std::exception& err){
       errorLaph(make_strf("Error on Task %d: %s\n",task,err.what()));}
    //catch(...){
    //   errorLaph(make_strf("Error on Task %d: generic\n",task));}
    swatch.stop();
    printLaph(make_strf("Task %d done using time = %g secs\n",
              task,swatch.getTimeInSeconds()));
    xml_in.seek_next_sibling();
    }

#ifdef ARCH_PARALLEL
 comm_barrier();
#endif
 QudaInfo::clearDevice();
 rolex.stop();
 printLaph(make_strf("\n\nQUDA_LAPH: total time = %g secs",
           rolex.getTimeInSeconds()));
 printLaph("QUDA_LAPH: completion");
 output_datetime();

 endQuda();
#ifdef ARCH_PARALLEL
 quda::comm_finalize();
 MPI_Finalize();
#endif

 return 0;
}

