#ifndef VERBOSITY_INFO_H
#define VERBOSITY_INFO_H

#include "xml_handler.h"
#include "quda.h"

namespace LaphEnv {


// ****************************************************************
// *                                                              *
// *  "Verbosity" stores information about the requested level    *
// *  of output. The different levels are "none", "low",          *
// *  "med" or "medium", and "high" or "full".                    *
// *                                                              *
// *  Required XML input for setting the info:                    *
// *                                                              *
// *   <Verbosity>low</Verbosity>                                 *
// *                                                              *
// ****************************************************************


class Verbosity
{

   QudaVerbosity_s level;

 public:

   Verbosity(XMLHandler& xml_in)
   {set_info(xml_in);}
   
   Verbosity(const std::string& tagstring)
   {set_info(tagstring);}
   
   Verbosity(int ivalue)
   {set_info(ivalue);}
   
   Verbosity() {level=QUDA_SILENT;}

   Verbosity(const Verbosity& vin) : level(vin.level) {}

   Verbosity& operator=(const Verbosity& vin)
   {level=vin.level; return *this;}
   
   Verbosity& setToNone()
   {level=QUDA_SILENT; return *this;}

   Verbosity& setToLow()
   {level=QUDA_SUMMARIZE; return *this;}

   Verbosity& setToMedium()
   {level=QUDA_VERBOSE; return *this;}

   Verbosity& setToFull()
   {level=QUDA_DEBUG_VERBOSE; return *this;}
   
   Verbosity& setToHigh()
   {level=QUDA_DEBUG_VERBOSE; return *this;}
   
       // set using "none", "low", "med", "medium", "high", "full"
   Verbosity& setByString(const std::string& invalue)
   {set_info(invalue); return *this;}

       // set using 0 (none), 1 (low), 2 (medium), 3 (high)
   Verbosity& setByInt(int ivalue)
   {set_info(ivalue); return *this;}

   ~Verbosity(){}


   QudaVerbosity_s getQudaValue() const {return level;}
   
   QudaVerbosity_s operator()() const {return level;}
   
   bool isNone() const {return (level==QUDA_SILENT);}

   bool isLow() const {return (level==QUDA_SUMMARIZE);}

   bool isMedium() const {return (level==QUDA_VERBOSE);}

   bool isHigh() const {return (level==QUDA_DEBUG_VERBOSE);}

   bool isLowOrHigher() const {return (level>=QUDA_SUMMARIZE);}

   bool isMediumOrHigher() const {return (level>=QUDA_VERBOSE);}

   bool operator==(const Verbosity& rhs) const
   {return level==rhs.level;}

   void output(XMLHandler& xmlout) const;

   std::string str() const;
   
   std::string message() const;
   

 private:

   void set_info(XMLHandler& xml_in);

   void set_info(const std::string& instr);

   void set_info(int ivalue);

};


void xml_read_if(XMLHandler& xmlin, Verbosity& v);


// *************************************************************
}
#endif
