#include "utils.h"
#include "laph_stdio.h"
#include "named_obj_map.h"
#include "quda.h"
#include "quda_info.h"
#include "util_quda.h"
#include "xml_handler.h"
#include <unistd.h>

#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif

using namespace quda;
using namespace std;

namespace LaphEnv {

//  Tests if file having name "file_name" exists on disk: works
//  in global or local mode.

bool doesFileExist(const std::string &file_name, bool global_mode) {
  bool result;
  if ((isPrimaryRank()) || (!global_mode)) {
    result = (access(file_name.c_str(), F_OK) == 0) ? true : false;
  }
#ifdef ARCH_PARALLEL
  if (global_mode)
    comm_broadcast(&result, sizeof(bool), 0);
#endif
  return result;
}

//  Tests if file having name "file_name" exists: works
//  in global or local mode. If starts with "NOM_", checks to
//  see if exists in the NamedObjMap.

bool fileExists(const std::string &file_name, bool global_mode) {
  bool result = false;
  if (file_name.find("NOM_") != 0) {
    if ((isPrimaryRank()) || (!global_mode)) {
      result = (access(file_name.c_str(), F_OK) == 0) ? true : false;
    }
#ifdef ARCH_PARALLEL
    if (global_mode) {
      comm_broadcast(&result, sizeof(bool), 0);
    }
#endif
  } else {
    string objname(tidyString(file_name));
    objname.erase(0, 4);
    result = NamedObjMap::query(objname);
  }
  return result;
}

bool emptyFileName(const std::string &file_name) {
  string ftidy(tidyString(file_name));
  return ((ftidy.empty()) || (ftidy == "NOM_"));
}

bool isFileOnDisk(const std::string &file_name) {
  return (tidyString(file_name).find("NOM_") != 0);
}

// Returns true if this rank is the primary one (rank 0), false otherwise;
// always returns true if serial (non-parallel)

bool isPrimaryRank() {
#ifdef ARCH_PARALLEL
  return (comm_rank() == 0);
#else
  return true;
#endif
}

// Sends a string from rank "send_rank" to all other ranks

void broadcastString(std::string &stringData, const int send_rank) {
#ifdef ARCH_PARALLEL
  char *cbuf;
  int nchar;

  // primary rank has the string
  if (isPrimaryRank()) {
    nchar = stringData.length();
  }

  // broadcast size of string
  comm_broadcast(&nchar, sizeof(int), send_rank);

  // each rank allocates space for string
  cbuf = new (std::nothrow) char[nchar];
  if (cbuf == 0x0) {
    errorLaph("Unable to allocate buffer in broadcastString\n");
  }

  if (isPrimaryRank()) {
    memcpy(cbuf, stringData.c_str(), nchar);
  }

  // broadcast char array out to all ranks
  comm_broadcast((void *)cbuf, nchar, send_rank);

  // all non-primary nodes can now make a string from the char buffer
  if (!isPrimaryRank()) {
    stringData.assign(cbuf, nchar);
  }

  // clean-up temporary buffer
  delete[] cbuf;
#endif
}

//  This does a logical "and" of all of the boolean variables "flag"
//  on each rank, and broadcasts the results to all ranks.

bool globalAnd(const bool &flag) {
#ifdef ARCH_SERIAL
  return flag;
#elif defined(ARCH_PARALLEL)
  bool result = flag;
  int status =
      MPI_Allreduce(&flag, &result, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in globalAnd");
  }
  return result;
#endif
}

//  This returns the maximum value of "ival" on all ranks,
//  and returns the result to the primary rank.

uint globalMax(const uint &uival) {
#ifdef ARCH_SERIAL
  return uival;
#elif defined(ARCH_PARALLEL)
  uint result = uival;
  int status =
      MPI_Allreduce(&uival, &result, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in globalMax");
  }
  return result;
#endif
}

float globalMax(const float &rval) {
#ifdef ARCH_SERIAL
  return rval;
#elif defined(ARCH_PARALLEL)
  float result = rval;
  int status =
      MPI_Allreduce(&rval, &result, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in globalMax");
  }
  return result;
#endif
}

double globalMax(const double &rval) {
#ifdef ARCH_SERIAL
  return rval;
#elif defined(ARCH_PARALLEL)
  double result = rval;
  int status =
      MPI_Allreduce(&rval, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in globalMax");
  }
  return result;
#endif
}

//  This returns the minimum value of "ival" on all ranks,
//  and returns the result to the primary rank.

uint globalMin(const uint &uival) {
#ifdef ARCH_SERIAL
  return uival;
#elif defined(ARCH_PARALLEL)
  uint result = uival;
  int status =
      MPI_Allreduce(&uival, &result, 1, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in globalMax");
  }
  return result;
#endif
}

// Returns true if integer "ival" has the same value on all ranks; false
// otherwise.

bool isIntegerSameAllRanks(const int &ival) {
#ifdef ARCH_SERIAL
  return true;
#elif defined(ARCH_PARALLEL)
  int p[2];
  p[0] = -ival;
  p[1] = ival;
  int status =
      MPI_Allreduce(MPI_IN_PLACE, p, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in isIntegerSameAllRanks");
  }
  return (p[0] == -p[1]);
#endif
}

// Returns true if boolean "bval" has the same value on all ranks; false
// otherwise.

bool isBooleanSameAllRanks(bool bval) {
#ifdef ARCH_SERIAL
  return true;
#elif defined(ARCH_PARALLEL)
  bool check = bval;
  bool result;
  int status =
      MPI_Allreduce(&check, &result, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in globalAnd");
  }
  if (result)
    return true;
  check = !bval;
  status =
      MPI_Allreduce(&check, &result, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in globalAnd");
  }
  return result;
#endif
}

// convert a vector of int to a vector of unsigned int

std::vector<uint> unsign(const std::vector<int> &ivector) {
  return vector<uint>(ivector.begin(), ivector.end());
}

#ifdef ARCH_PARALLEL

// Send from one mpi rank to another (blocking)

void commSendRecv(void *buffer, const size_t nbytes, const int send_rank,
                  const int recv_rank) {
  if (send_rank == recv_rank)
    return;
  const int nranks = comm_size();
  const int myrank = comm_rank();
  if ((send_rank >= nranks) || (recv_rank >= nranks)) {
    throw(std::runtime_error("Invalid MPI ranks given in commSendRecv"));
  }
  if (myrank == send_rank) {
    MPI_Send(buffer, nbytes, MPI_CHAR, recv_rank, 0, MPI_COMM_WORLD);
  } else if (myrank == recv_rank) {
    MPI_Status statusinfo;
    int status = MPI_Recv(buffer, nbytes, MPI_CHAR, send_rank, 0,
                          MPI_COMM_WORLD, &statusinfo);
    if (status != MPI_SUCCESS) {
      errorLaph("Problem occurred in commSendRecv");
    }
  }
  comm_barrier();
}

size_t comm_size() // number of MPI ranks
{
  return quda::comm_size();
}

void comm_broadcast(void *data, const size_t nbytes, const int root) {
  quda::comm_broadcast(data, nbytes, root);
}

void comm_barrier() { quda::comm_barrier(); }

int comm_rank() { return quda::comm_rank(); }

int commCoords(const int rank) { return quda::commCoords(rank); }

#endif

void laph_abort() {
  QudaInfo::clearDevice();
  errorQuda("aborting");
}

//  Extracts the integer value from the child tag named "tagname"
//  in the XMLhandler "xmlin", and saves it in "ivalues[index]",
//  incrementing "index" afterwards. If "optional" is true,
//  the "default_value" is stored if the child tag cannot be
//  found.  If not optional and the child tag cannot be found,
//  an exception is thrown.

void xmlsetQLInt(XMLHandler &xmlin, const std::string &tagname,
                 vector<int> &ivalues, int &index, const bool optional,
                 const int default_value) {
  if (index >= int(ivalues.size())) {
    throw(std::runtime_error("ivalues vector not large enough for QL int set"));
  }
  XMLHandler xmlh(xmlin);
  try {
    bool state = xmlh.get_exceptions();
    xmlh.set_exceptions_on();
    xmlh.seek_child(tagname);
    std::string content = xmlh.get_text_content();
    extract_from_string(content, ivalues[index]);
    xmlh.set_exceptions(state);
  } catch (const std::exception &err_msg) {
    if (optional) {
      ivalues[index] = default_value;
    } else {
      throw(std::invalid_argument(string("Invalid QL integer read: ") +
                                  err_msg.what()));
    }
  }
  ++index;
}

//  Same as the above routine, except the value is of type bool.
//  Will be stored as int though.

void xmlsetQLBool(XMLHandler &xmlin, const std::string &tagname,
                  vector<int> &ivalues, int &index, const bool optional,
                  const bool default_value) {
  if (index >= int(ivalues.size())) {
    throw(
        std::runtime_error("ivalues vector not large enough for QL bool set"));
  }
  XMLHandler xmlh(xmlin);
  try {
    bool state = xmlh.get_exceptions();
    xmlh.set_exceptions_on();
    xmlh.seek_child(tagname);
    std::string content = xmlh.get_text_content();
    if (content == "true") {
      ivalues[index] = 1;
    } else if (content == "false") {
      ivalues[index] = 0;
    } else {
      throw(std::invalid_argument("Invalid bool"));
    }
    xmlh.set_exceptions(state);
  } catch (const std::exception &err_msg) {
    if (optional) {
      ivalues[index] = int(default_value);
    } else {
      throw(std::invalid_argument(string("Invalid QL bool read: ") +
                                  err_msg.what()));
    }
  }
  ++index;
}

//  Same as the above routine, except the value is of type double.

void xmlsetQLReal(XMLHandler &xmlin, const std::string &tagname,
                  vector<double> &rvalues, int &index, const bool optional,
                  const double default_value) {
  if (index >= int(rvalues.size())) {
    throw(
        std::runtime_error("rvalues vector not large enough for QL real set"));
  }
  XMLHandler xmlh(xmlin);
  try {
    bool state = xmlh.get_exceptions();
    xmlh.set_exceptions_on();
    xmlh.seek_child(tagname);
    std::string content = xmlh.get_text_content();
    extract_from_string(content, rvalues[index]);
    xmlh.set_exceptions(state);
  } catch (const std::exception &err_msg) {
    if (optional) {
      rvalues[index] = default_value;
    } else {
      throw(std::invalid_argument(string("Invalid QL real read: ") +
                                  err_msg.what()));
    }
  }
  ++index;
}

void xmlsetQLIntVector(XMLHandler &xmlin, const std::string &tagname,
                       vector<int> &ivalues, int &index, const int nvalues,
                       const bool optional, const vector<int> &default_value) {
  if (index >= (int(ivalues.size()) - nvalues + 1)) {
    throw(std::runtime_error(
        "ivalues vector not large enough for QL int vector set"));
  }
  if ((optional) && (int(default_value.size()) != nvalues)) {
    throw(std::runtime_error("default value wrong size for QL int vector set"));
  }
  XMLHandler xmlh(xmlin);
  try {
    bool state = xmlh.get_exceptions();
    xmlh.set_exceptions_on();
    xmlh.seek_child(tagname);
    std::string content = xmlh.get_text_content();
    vector<int> values_read;
    extract_from_string(content, values_read);
    if (nvalues != int(values_read.size())) {
      throw(std::invalid_argument(
          "Number of values is incorrect for int vector read"));
    }
    for (int j = 0; j < nvalues; ++j) {
      ivalues[index] = values_read[j];
      ++index;
    }
    xmlh.set_exceptions(state);
  } catch (const std::exception &err_msg) {
    if (optional) {
      for (int j = 0; j < nvalues; ++j) {
        ivalues[index] = default_value[j];
        ++index;
      }
    } else {
      throw(std::invalid_argument(string("Invalid QL int vector read: ") +
                                  err_msg.what()));
    }
  }
}

//  Similar to the above routine, except the value is an enum.
//  The actual value to be read must be a string, and the string
//  must match one of the strings in "xmlchoices".  If no match
//  found and not optional, an exception is thrown.  Note that
//  the **index** in xmlchoices is saved, not the string itself.
//  If optional and no tag found, the default index is saved.

void xmlsetQLEnum(XMLHandler &xmlin, const std::string &tagname,
                  const vector<string> &xmlchoices, vector<int> &ivalues,
                  int &index, const bool optional, const int default_index) {
  if (index >= int(ivalues.size())) {
    throw(
        std::runtime_error("ivalues vector not large enough for QL enum set"));
  }
  XMLHandler xmlh(xmlin);
  bool invalid = false;
  try {
    bool state = xmlh.get_exceptions();
    xmlh.set_exceptions_on();
    xmlh.seek_child(tagname);
    std::string content = xmlh.get_text_content();
    invalid = true;
    for (uint j = 0; j < xmlchoices.size(); ++j) {
      if (content == xmlchoices[j]) {
        ivalues[index] = j;
        invalid = false;
        break;
      }
    }
    if (invalid) {
      throw(std::invalid_argument("Not one of the allowed selections"));
    }
    xmlh.set_exceptions(state);
  } catch (const std::exception &err_msg) {
    if ((optional) && (!invalid)) {
      ivalues[index] = default_index;
    } else {
      throw(std::invalid_argument(string("Invalid QL enum read: ") +
                                  err_msg.what()));
    }
  }
  ++index;
}

//  Outputs the integer value into "xmlout" using value stored in
//  "ivalues[index]".  "index" is incremented afterwards.

XMLHandler xmloutputQLInt(const std::string &tagname,
                          const vector<int> &ivalues, int &index) {
  if (index >= int(ivalues.size())) {
    throw(std::runtime_error("ivalues vector not large enough for QL output"));
  }
  string val(make_string(ivalues[index]));
  ++index;
  return XMLHandler(tagname, val);
}

//  Same as the above routine, except the value is of type bool.

XMLHandler xmloutputQLBool(const std::string &tagname,
                           const vector<int> &ivalues, int &index) {
  if (index >= int(ivalues.size())) {
    throw(std::runtime_error("ivalues vector not large enough for QL output"));
  }
  string val((ivalues[index] == 0) ? "false" : "true");
  ++index;
  return XMLHandler(tagname, val);
}

//  Same as the above routine, except the value is of type double.

XMLHandler xmloutputQLReal(const std::string &tagname,
                           const vector<double> &rvalues, int &index) {
  if (index >= int(rvalues.size())) {
    throw(std::runtime_error("rvalues vector not large enough for QL output"));
  }
  string val(make_string(rvalues[index]));
  ++index;
  return XMLHandler(tagname, val);
}

XMLHandler xmloutputQLIntVector(const std::string &tagname,
                                const vector<int> &ivalues, int &index,
                                const int nvalues) {
  if (index >= (int(ivalues.size()) - nvalues + 1)) {
    throw(std::runtime_error("ivalues vector not large enough for QL output"));
  }
  vector<int> res(nvalues);
  for (int j = 0; j < nvalues; ++j) {
    res[j] = ivalues[index];
    ++index;
  }
  string val(make_string(res));
  return XMLHandler(tagname, val);
}

//  Same as above, except the enum string is output.

XMLHandler xmloutputQLEnum(const std::string &tagname,
                           const vector<string> &xmlchoices,
                           const vector<int> &ivalues, int &index) {
  if (index >= int(ivalues.size())) {
    throw(std::runtime_error("ivalues vector not large enough for QL output"));
  }
  string val(xmlchoices[ivalues[index]]);
  ++index;
  return XMLHandler(tagname, val);
}

//  Outputs the integer value into "xmlout" using value stored in
//  "ivalues[index]".  "index" is incremented afterwards.
//  The "tagname" is not used but is informational.

int xmlputQLInt(const std::string &tagname, const vector<int> &ivalues,
                int &index) {
  if (index >= int(ivalues.size())) {
    throw(std::runtime_error("ivalues vector not large enough for QL put"));
  }
  int result = ivalues[index];
  ++index;
  return result;
}

//  Same as the above routine, except the value is of type bool.

bool xmlputQLBool(const std::string &tagname, const vector<int> &ivalues,
                  int &index) {
  if (index >= int(ivalues.size())) {
    throw(std::runtime_error("ivalues vector not large enough for QL put"));
  }
  bool result = (ivalues[index] == 0) ? false : true;
  ++index;
  return result;
}

//  Same as the above routine, except the value is of type double.

double xmlputQLReal(const std::string &tagname, const vector<double> &rvalues,
                    int &index) {
  if (index >= int(rvalues.size())) {
    throw(std::runtime_error("rvalues vector not large enough for QL put"));
  }
  double result = rvalues[index];
  ++index;
  return result;
}

vector<int> xmlputQLIntVector(const std::string &tagname,
                              const vector<int> &ivalues, int &index,
                              const int nvalues) {
  if (index >= (int(ivalues.size()) - nvalues + 1)) {
    throw(std::runtime_error("ivalues vector not large enough for QL put"));
  }
  vector<int> res(nvalues);
  for (int j = 0; j < nvalues; ++j) {
    res[j] = ivalues[index];
    ++index;
  }
  return res;
}

//  Same as above, except the enum integer is output.

int xmlputQLEnum(const std::string &tagname, const std::vector<int> &quda_enums,
                 const vector<int> &ivalues, int &index) {
  if (index >= int(ivalues.size())) {
    throw(std::runtime_error("ivalues vector not large enough for QL put"));
  }
  int result = quda_enums[ivalues[index]];
  ++index;
  return result;
}

//  This routine looks for the tag "tagname" in the current children elements
//  and returns the XML content under it, returning an XMLHandler with "tagname"
//  as the root; if "tagname" is not found, then an XMLHandler is returned
//  with just the root tag of "tagname".

XMLHandler getXML_nofail(XMLHandler xmlin, const std::string &tagname) {
  XMLHandler xmlh(xmlin);
  try {
    xmlh.set_exceptions_on();
    xmlh.seek_child(tagname);
    return XMLHandler(xmlh, tagname);
  } catch (const std::exception &xp) {
    xmlh.set(tagname.c_str());
  }
  return xmlh;
}
} // namespace LaphEnv
