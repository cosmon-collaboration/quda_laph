#ifndef LAPH_MODE_DOUBLET_INFO_H
#define LAPH_MODE_DOUBLET_INFO_H

#include "quda.h"
#include <set>
#include <array>
#include "momenta.h"
#include "xml_handler.h"

namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Objects of class "ModeDoubletInfo" store identifying info     *
// *   about one particular mode doublet.                            *
// *   The XML input into the Chroma tasks must inform Chroma        *
// *   which mode doublets to compute, so a pre-defined input        *
// *   format must be followed.  The needed input format and the     *
// *   required storage format is described below.                   *
// *                                                                 *
// *   Chroma input format for each mode doublet:                    *
// *                                                                 *
// *       <ModeDoublet>                                             *
// *          <Momentum>  0 0 0  </Momentum>                         *
// *       at the moment, only SS and SD operators are supported     *
//        the quark is displaced in an SD mode doublet. 
//            <DispLength> 3 </DispLength>                           *
// *          <DispDir> -1 </DispDir> 
// *      0 for no displacement, +/-1,2,3 for the direction of the 
//        displaced quark. 
// *       </ModeDoublet>                                            *
// *                                                                 *
// *   READING multiple mode doublets:                               *
// *                                                                 *
// *      <ModeDoublets>                                             *
// *           <ModeDoublet>                                         *
// *            ...                                                  *
// *           </ModeDoublet>                                        *
// *           <ModeDoublet>                                         *
// *            ...                                                  *
// *           </ModeDoublet>                                        *
// *         ...                                                     *
// *      </ModeDoublets>                                            *
// *                                                                 *
// *   Given the above XML input, one just needs to issue the        *
// *   following:                                                    *
// *                                                                 *
// *      XmlReader xml_in;  // assigned somehow                     *
// *      set<ModeDoubletInfo> mt_set;                               *
// *      readModeDoublets(xml_in,mt_set);                           *
// *                                                                 *
// *   The above code initializes the class, then reads all of       *
// *   mode doublets inside the <ModeDoublets> tags in               *
// *   the input XML file, storing the results in a set.             *
// *                                                                 *
// *******************************************************************

class ModeDoubletInfo {

   unsigned icode; 

 public:

  ModeDoubletInfo(XMLHandler& xml_in);

  ModeDoubletInfo() : icode{0} {}

    // output functions

  Momentum getMomentum() const; 
  int getXMomentum() const;
  int getYMomentum() const;
  int getZMomentum() const;

  int getDisplacementLength() const;
  int getDisplacementDir() const;

  std::string output(int indent) const;  // XML output of  tags

  void output(XMLHandler& xmlout) const;  // XML output of baryon tags

  bool operator==(const ModeDoubletInfo& rhs) const;
  bool operator!=(const ModeDoubletInfo& rhs) const;
  bool operator<(const ModeDoubletInfo& rhs) const;

  friend class ModeDoubletHandler;

};



// **************************************************

   //  useful routine for reading multiple ModeDoublets

void readModeDoublets(XMLHandler& xml_in, 
                         std::set<ModeDoubletInfo>& mdoubs);


// **************************************************
  }
#endif  
