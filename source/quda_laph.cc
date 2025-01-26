#include "tasks.h"
#include "utils.h"
#include "laph_stdio.h"
#include "stop_watch.h"
#include "quda_info.h"
#include "verbosity_info.h"
#include "host_global.h"
#include "layout_info.h"
#if defined(QUDA_OPENMP)
#include <omp.h>
#endif
#ifdef ARCH_PARALLEL
#include <mpi.h>
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
// *   is echoed to standard output. A <Verbosity> tag can also be      *
// *   included: this will be the default verbosity in any task that    *
// *   does not have its own verbosity set.                             *
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
// *         <XYZTExtents>24 24 24 96</XYZTExtents>                     *
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
 TaskMap["LAPH_QUARK_PERAMBULATORS"]=&doLaphQuarkPerambulators;
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
  srand(time(NULL) * rank + 1073741827);
}


 //  Global data (persistent across tasks)

std::vector<LattField> HostGlobal::theGaugeConfig;
std::vector<LattField> HostGlobal::theSmearedGaugeConfig;
std::vector<LattField> HostGlobal::theLaphEigvecs;


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
 int rank;
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 if ((int(npartitions.size())!=LayoutInfo::Ndim)||(inputxmlfile.empty())){
    if (rank==0){
       std::cout <<endl<<endl<< "Invalid command line options: should be "
                 << "-i <inputxmlfile> -npartitions <nx> <ny> <nz> <nt>"<<endl<<endl;}
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);}

      // We want the MPI ranks on each node to correspond to sublattices separated 
      // in space x,y,z, not in time.  
 char comm_dist=LayoutInfo::getMPIRankDistributionOnNodes();
 QudaCommsMap commfunc=nullptr;
 void *commfuncdata=nullptr;
 if (comm_dist=='B'){
    commfunc=LayoutInfo::lex_rank_from_commcoord_x_fastest;
    commfuncdata=npartitions.data();}
 else if (comm_dist=='C'){
    commfunc=LayoutInfo::lex_rank_from_commcoord_t_fastest;
    commfuncdata=npartitions.data();}
 MPI_Barrier(MPI_COMM_WORLD);

 initCommsGridQuda(LayoutInfo::Ndim,npartitions.data(),commfunc,commfuncdata);
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

#ifdef ARCH_PARALLEL
 if (comm_dist=='B'){
    printLaph("\nMPI rank distribution is block (consecutive ranks on same node)");
    printLaph("Sublattices laid out with x comm coordinate varying fastest\n");}
 else if (comm_dist=='C'){
    printLaph("\nMPI rank distribution is cyclic (consecutive ranks on consecutive nodes)");
    printLaph("Sublattices laid out with time comm coordinate varying fastest\n");}
 else{
    printLaph("\nMPI rank distribution is not block nor cyclic");
    printLaph("Sublattices laid out with QUDA default\n");}
#endif

 bool echo=false;
 bool layoutinfo=false;
 bool qudainfo=false;
 Verbosity verbosity(0);
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
          verbosity=Verbosity(xml_in);}
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
 setVerbosity(verbosity.getQudaValue()); 

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
       {const int lsize=62;
       string liner0(lsize+2,'#');
       string liner1(lsize,' ');
       printLaph(make_strf("\n\n  %s",liner0.c_str()));
       printLaph(make_strf("  %c%s%c",'#',liner1.c_str(),'#'));
       string liner2b("Starting Task "); liner2b+=make_string(task);
       string liner2a((lsize-liner2b.size())/2,' ');
       string liner2c(lsize-liner2a.size()-liner2b.size(),' ');
       printLaph(make_str("  #",liner2a,liner2b,liner2c,"#"));
       printLaph(make_strf("  %c%s%c",'#',liner1.c_str(),'#'));
       printLaph(make_strf("  %s\n\n",liner0.c_str()));}
       if (xml_in.fail()) throw(std::runtime_error("XML input error"));
       while (xml_in.get_node_name()!="Task"){
          xml_in.seek_next_sibling();}
       XMLHandler xml_task(xml_in);
           // the main task
       T.do_task(xml_task,echo);
       setVerbosity(verbosity.getQudaValue());   // set verbosity back to default after task may have changed it
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

