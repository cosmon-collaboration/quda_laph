#ifndef LAPH_UTILS_H
#define LAPH_UTILS_H

#include <string>
#include <cstdarg>
#include "byte_handler.h"
#include "xml_handler.h"

namespace LaphEnv {

// *********************************************************

    //  Tests if file having name "file_name" exists on disk: works
    //  in global or local mode. 
    
bool doesFileExist(const std::string& file_name, bool global_mode=true);

    //  Tests if file having name "file_name" exists: works
    //  in global or local mode. If starts with "NOM_", checks to
    //  see if exists in the NamedObjMap.

bool fileExists(const std::string& file_name, bool global_mode=true);

    //  Returns true if name of file is empty or "NOM_" only
    
bool emptyFileName(const std::string& file_name);

    //  Returns true if name of disk file, false is name
    //  starts with "NOM_" indicating the NamedObjMap
    
bool isFileOnDisk(const std::string& file_name);


// *********************************************************

   // Returns true if this rank is the primary one (rank 0), false otherwise;
   // always returns true if serial (non-parallel)

bool isPrimaryRank();


// *********************************************************

   // Sends a string from from rank "send_rank" to all other ranks

void broadcastString(std::string& stringData, int send_rank=0);


// *********************************************************

    //  This does a logical "and" of all of the boolean variables "flag"
    //  on each rank, and returns the result to the primary rank.
    
bool globalAnd(const bool& flag);

    //  This returns the maximum value of "ival" on all ranks,
    //  and returns the result to the primary rank.
    
uint globalMax(const uint& uival);
float globalMax(const float& rval);
double globalMax(const double& rval);

    //  This returns the minimum value of "ival" on all ranks,
    //  and returns the result to the primary rank.
    
uint globalMin(const uint& uival);

     // Returns true if "ival" has the same value on all ranks; false otherwise.

bool isIntegerSameAllRanks(const int& ival);

     // Returns true if boolean "bval" has the same value on all ranks; false otherwise.

bool isBooleanSameAllRanks(bool bval);


// *********************************************************

#ifdef ARCH_PARALLEL

    // Send from one mpi rank to another (blocking)

void commSendRecv(void *buffer, size_t nbytes, int send_rank, int recv_rank);

size_t comm_size();  // number of MPI ranks

void comm_broadcast(void *data, size_t nbytes, int root = 0);

void comm_barrier();

int comm_rank();

int commCoords(int);

#endif

void laph_abort();


// **********************************************************************

   //  Extracts the integer value from the child tag named "tagname"
   //  in the XMLhandler "xmlin", and saves it in "ivalues[index]",
   //  incrementing "index" afterwards. If "optional" is true,
   //  the "default_value" is stored if the child tag cannot be 
   //  found.  If not optional and the child tag cannot be found,
   //  an exception is thrown.

void xmlsetQLInt(XMLHandler& xmlin, const std::string& tagname, 
                 std::vector<int>& ivalues, int& index,
                 bool optional=false, int default_value=0);

void xmlsetQLBool(XMLHandler& xmlin, const std::string& tagname, 
                  std::vector<int>& ivalues, int& index,
                  bool optional=false, bool default_value=false);

void xmlsetQLReal(XMLHandler& xmlin, const std::string& tagname, 
                  std::vector<double>& rvalues, int& index,
                  bool optional=false, double default_value=0.0);

void xmlsetQLIntVector(XMLHandler& xmlin, const std::string& tagname, 
                       std::vector<int>& ivalues, int& index, int nvalues,
                       bool optional=false, const std::vector<int>& default_value={});

   //  Similar to the above routine, except the value is an enum.
   //  The actual value to be read must be a string, and the string
   //  must match one of the strings in "xmlchoices".  If no match
   //  found and not optional, an exception is thrown.  Note that
   //  the **index** in xmlchoices is saved, not the string itself.
   //  If optional and no tag found, the default index is saved.

void xmlsetQLEnum(XMLHandler& xmlin, const std::string& tagname,
                  const std::vector<std::string>& xmlchoices,
                  std::vector<int>& ivalues, int& index,
                  bool optional=false, int default_index=0);


   //  Outputs the integer value into "xmlout" using value stored in
   //  "ivalues[index]".  "index" is incremented afterwards.

XMLHandler xmloutputQLInt(const std::string& tagname, const std::vector<int>& ivalues, 
                          int& index);

XMLHandler xmloutputQLBool(const std::string& tagname, const std::vector<int>& ivalues, 
                           int& index);

XMLHandler xmloutputQLReal(const std::string& tagname, const std::vector<double>& rvalues, 
                           int& index);

XMLHandler xmloutputQLIntVector(const std::string& tagname, const std::vector<int>& ivalues, 
                                int& index, int nvalues);

XMLHandler xmloutputQLEnum(const std::string& tagname, const std::vector<std::string>& xmlchoices,
                           const std::vector<int>& ivalues, int& index);


   //  Outputs the integer value into "xmlout" using value stored in
   //  "ivalues[index]".  "index" is incremented afterwards.
   //  The "tagname" is not used but is informational.

int xmlputQLInt(const std::string& tagname, const std::vector<int>& ivalues, int& index);

bool xmlputQLBool(const std::string& tagname, const std::vector<int>& ivalues, int& index);

double xmlputQLReal(const std::string& tagname, const std::vector<double>& rvalues, int& index);

std::vector<int> xmlputQLIntVector(const std::string& tagname, const std::vector<int>& ivalues, 
                                   int& index, int nvalues);

int xmlputQLEnum(const std::string& tagname, const std::vector<int>& quda_enums, 
                 const std::vector<int>& ivalues, int& index);

   //  This routine looks for the tag "tagname" in the current children elements
   //  and returns the XML content descending under it, returning an XMLHandler with 
   //  "tagname" as the root; if "tagname" is not found, then an XMLHandler is returned
   //  with just the root tag of "tagname". 

XMLHandler getXML_nofail(XMLHandler xmlin, const std::string& tagname);

// *************************************************************************************
}
#endif
