#include "laph_stdio.h"
#include "layout_info.h"
#include "named_obj_map.h"
#include "quda.h"
#include "quda_info.h"
#include "stop_watch.h"
#include "tasks.h"
#include "utils.h"
#include "xml_handler.h"
#include <cstdio>
#include <ctime>
#include <map>
#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif

using namespace std;
using namespace LaphEnv;

map<string, NamedObjBase *> NamedObjMap::the_map;

static StopWatch rolex;

static bool echo = false;
static bool layoutinfo = false;
static bool qudainfo = false;
static QudaVerbosity_s verbosity = QUDA_SILENT;

static void parse_args(int *argc, char ***argv, std::vector<int> &npartitions,
                       std::string &inputxmlfile) {
  npartitions.clear();
  inputxmlfile.clear();
  int nargs = *argc;
  const int Nd = 4; // four dimensions required

  for (int i = 0; i < nargs; ++i) {
    string argv_i((*argv)[i]);
    //  get the input xml file name
    if ((argv_i == std::string("-i")) && ((i + 1) < nargs)) {
      inputxmlfile = std::string((*argv)[i + 1]);
      ++i;
    } else if ((argv_i == "-npartitions") && ((i + Nd) < nargs)) {
      npartitions.resize(Nd);
      for (int j = 0; j < Nd; j++) {
        int uu;
        sscanf((*argv)[++i], "%d", &uu);
        if (uu <= 0) {
          npartitions.clear();
          return;
        }
        npartitions[j] = uu;
      }
    }
  }
}

static void output_datetime() {
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  printLaph(make_strf("  Current date/time: %s\n", asctime(timeinfo)));
}

//  Constructor sets up the known tasks.  Call
//  member function "do_task" to perform the task.
//  "TaskMap" is a map that associates a string
//  containing a task name to a function pointer.

class Tasker {
  typedef void (*task_ptr)(XMLHandler &);
  std::map<string, task_ptr> TaskMap;

public:
  Tasker();
  ~Tasker() {}
  void do_task(XMLHandler &xml_in, bool echo);
};

// set up the known tasks
Tasker::Tasker() {
  TaskMap["SMEAR_GAUGE_FIELD"] = &doSmearGaugeField;
  TaskMap["SMEAR_QUARK_FIELD"] = &doSmearQuarkField;
  TaskMap["LAPH_QUARK_LINE_ENDS"] = &doLaphQuarkLineEnds;
  TaskMap["LAPH_QUARK_PERAMBULATORS"] = &doLaphQuarkPerambulators;
};

void Tasker::do_task(XMLHandler &xml_task, bool echo) {
#ifdef ARCH_PARALLEL
  comm_barrier();
#endif
  if (echo) {
    printLaph("Input XML for this task:");
    printLaph(xml_task.output());
  }
  XMLHandler xmlt(xml_task);
  xml_child_assert(xmlt, "Name", "do_task");
  if (!xmlt.is_simple_element()) {
    throw(std::runtime_error("<Name> tag absent or is not simple XML element"));
  }
  string task_name = xmlt.get_text_content();
  printLaph(make_strf("  Task name = %s\n", task_name));

  map<string, task_ptr>::iterator taskit = TaskMap.find(task_name);
  if (taskit != TaskMap.end()) {
    (*(taskit->second))(xml_task);
  } // do the task!!
  else {
    throw(std::runtime_error("Unknown task name"));
  } // unknown task?
}

static void initRand() {
  int rank = 0;
#if defined(ARCH_PARALLEL)
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  srand(17 * rank + 137);
}

static void readVerbosity(XMLHandler &xmlin, QudaVerbosity_s &verbosity) {
  string tagvalue;
  xmlread(xmlin, "Verbosity", tagvalue, "QudaLaph");
  if (tagvalue == "none") {
    verbosity = QUDA_SILENT;
    printLaph("Verbosity set to none: SILENT\n");
  } else if (tagvalue == "low") {
    verbosity = QUDA_SUMMARIZE;
    printLaph("Verbosity set to medium: SUMMARIZE\n");
  } else if ((tagvalue == "medium") || (tagvalue == "med")) {
    verbosity = QUDA_VERBOSE;
    printLaph("Verbosity set to medium: VERBOSE\n");
  } else if (tagvalue == "high") {
    verbosity = QUDA_DEBUG_VERBOSE;
    printLaph("Verbosity set to high: DEBUG_VERBOSE\n");
  }
}

// inits MPI and Quda, returns an XMLHandler object
int init_quda_laph(int argc, char *argv[], XMLHandler &xml_in) {
  // parse the command line options -i and -npartitions
  vector<int> npartitions;
  string inputxmlfile;
  parse_args(&argc, &argv, npartitions, inputxmlfile);

  // init MPI
#ifdef ARCH_PARALLEL
  MPI_Init(&argc, &argv);
  if ((int(npartitions.size()) != LayoutInfo::Ndim) || (inputxmlfile.empty())) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      std::cout << endl
                << endl
                << "Invalid command line options: should be "
                << "-i <inputxmlfile> -npartitions <nx> <ny> <nz> <nt>" << endl
                << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  initCommsGridQuda(LayoutInfo::Ndim, npartitions.data(), nullptr, nullptr);
  comm_barrier();
#else
  if (inputxmlfile.empty()) {
    std::cout << endl
              << endl
              << "Invalid command line options: should be "
              << "-i <inputxmlfile> -npartitions <nx> <ny> <nz> <nt>" << endl
              << endl;
    std::exit(1);
  }
#endif

  fflush(stdout);
  initRand();
  printLaph("\nStarting QUDA_LAPH");
  output_datetime();

  rolex.start();
  printLaph(make_strf("  Name of input XML file: %s\n\n", inputxmlfile));

  xml_in.set_from_file(inputxmlfile);
  if (xml_in.fail()) {
    errorLaph(
        "  Unable to read/parse XML content in input file ... exiting...");
  }
  if (xml_in.get_node_name() != "QudaLaph") {
    errorLaph(
        "  Root tag of input XML must be named \"QudaLaph\" ... exiting...");
  }

  char *qrpath = getenv("QUDA_RESOURCE_PATH");
  if (!qrpath) {
    printLaph("Environment variable QUDA_RESOURCE_PATH is not set.");
    printLaph("Ensure it is set to allow caching of tuned Quda parameters.");
  } else {
    printLaph(make_strf("QUDA_RESOURCE_PATH is set to %s\n", qrpath));
  }

  XMLHandler xml_layoutinfo, xml_qudainfo;
  int ntasks = 0;
  try {
    xml_in.seek_first_child();
    while (xml_in.good()) {
      if (xml_in.get_node_name() == "Task")
        ntasks++;
      else if (xml_in.get_node_name() == "EchoXML")
        echo = true;
      else if (xml_in.get_node_name() == "LatticeLayoutInfo") {
        layoutinfo = true;
        xml_layoutinfo.set(xml_in);
      } else if (xml_in.get_node_name() == "QudaInfo") {
        qudainfo = true;
        xml_qudainfo.set(xml_in);
      } else if (xml_in.get_node_name() == "Verbosity") {
        readVerbosity(xml_in, verbosity);
      } else
        throw(std::runtime_error("Invalid XML input"));
      xml_in.seek_next_sibling();
    }
  } catch (const std::exception &err) {
    errorLaph(make_strf("  Error: %s ... exiting...\n", err.what()));
  }

  if (!layoutinfo) {
    errorLaph("No <LatticeLayoutInfo> tag was found as child of root tag");
  }
  try {
    LayoutInfo::init(xml_layoutinfo, npartitions, echo);
  } catch (const std::exception &err) {
    errorLaph("Problem reading contents of <LatticeLayoutInfo>");
  }
  // if non-default quda parameters requested, set them
  if (qudainfo) {
    QudaInfo::init(xml_qudainfo, echo);
  }
  // set the verbosity level for Quda output
  setVerbosity(verbosity);

  // Quda initializations
  initQuda(QudaInfo::getDeviceOrdinal());

  printLaph(make_strf("\n\n  Number of tasks is %d", ntasks));
#ifdef OPENMP
  printLaph(
      make_strf("  Maximum number of threads is %d", omp_get_max_threads));
#endif
#ifdef ARCH_PARALLEL
  printLaph(make_strf("  Number of MPI ranks = %d", comm_size()));
#endif

  return ntasks;
}

void run_tasks(XMLHandler &xml_in, const int ntasks) {
  xml_in.seek_root();
  xml_in.seek_first_child();
  Tasker T;

  for (int task = 1; task <= ntasks; ++task) {
    StopWatch swatch;
    swatch.start();
    try {
      printLaph(make_strf("\n\nStarting Task %d\n", task));
      if (xml_in.fail())
        throw(std::runtime_error("XML input error"));
      while (xml_in.get_node_name() != "Task") {
        xml_in.seek_next_sibling();
      }
      XMLHandler xml_task(xml_in);
      // the main task
      T.do_task(xml_task, echo);
      setVerbosity(verbosity); // set verbosity back to default after task may
                               // have changed it
    } catch (const std::exception &err) {
      errorLaph(make_strf("Error on Task %d: %s\n", task, err.what()));
    }
    swatch.stop();
    printLaph(make_strf("Task %d done using time = %g secs\n", task,
                        swatch.getTimeInSeconds()));
    xml_in.seek_next_sibling();
  }
  return;
}

void finalize(void) {
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
  return;
}
