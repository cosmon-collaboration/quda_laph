#ifndef LAPH_MODE_TRIPLET_INFO_H
#define LAPH_MODE_TRIPLET_INFO_H

#include "quda.h"
#include <set>
#include <array>
#include "momenta.h"
#include "xml_handler.h"

namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Objects of class "ModeTripletInfo" store identifying info     *
// *   about one particular mode triplet.                            *
// *   The XML input into the Chroma tasks must inform Chroma        *
// *   which mode triplets to compute, so a pre-defined input        *
// *   format must be followed.  The needed input format and the     *
// *   required storage format is described below.                   *
// *                                                                 *
// *   input format for each mode triplet:                           *
// *                                                                 *
// *       <ModeTriplet>                                             *
// *          <Momentum>  0 0 0  </Momentum>                         *
// *       at the moment, only SS and SD operators are supported     *
// *          <DisplacedQuark> <0,1,2,3> </DisplacedQuark>
// *      0 = SS (no displacements), 1,2, or 3 is which of the quarks
//        is displaced in an SD mode triplet. 
// *          
//            <DispLength> 3 </DispLength>                           *
// *          <DispDir> -1 </DispDir> 
// *      0 for no displacement, +/-1,2,3 for the direction of the 
//        displaced quark. 
// *       </ModeTriplet>                                            *
// *                                                                 *
// *   READING multiple mode triplets:                               *
// *                                                                 *
// *   input format:                                                 *
// *                                                                 *
// *      <ModeTriplets>                                             *
// *           <ModeTriplet>                                         *
// *            ...                                                  *
// *           </ModeTriplet>                                        *
// *           <ModeTriplet>                                         *
// *            ...                                                  *
// *           </ModeTriplet>                                        *
// *         ...                                                     *
// *      </ModeTriplets>                                            *
// *                                                                 *
// *   Given the above XML input, one just needs to issue the        *
// *   following:                                                    *
// *                                                                 *
// *      XmlReader xml_in;  // assigned somehow                     *
// *      set<ModeTripletInfo> mt_set;                               *
// *      readModeTriplets(xml_in,mt_set);                           *
// *                                                                 *
// *   The above code initializes the class, then reads all of       *
// *   mode triplets inside the <ModeTriplets> tags in               *
// *   the input XML file, storing the results in a set.             *
// *                                                                 *
// *******************************************************************

class ModeTripletInfo {

   unsigned icode; 

 public:

  ModeTripletInfo(XMLHandler& xml_in);

  ModeTripletInfo() : icode{0} {}

    // output functions

  Momentum getMomentum() const; 
  int getXMomentum() const;
  int getYMomentum() const;
  int getZMomentum() const;

  int getDisplacementLength() const;
  int getDisplacementDir() const;
	int getDisplacedQuark() const;

  std::string output(int indent) const;  // XML output of  tags

  void output(XMLHandler& xmlout) const;  // XML output of baryon tags

  bool operator==(const ModeTripletInfo& rhs) const;
  bool operator!=(const ModeTripletInfo& rhs) const;
  bool operator<(const ModeTripletInfo& rhs) const;

  friend class ModeTripletHandler;

};



// **************************************************

   //  useful routine for reading multiple ModeTriplets

void readModeTriplets(XMLHandler& xml_in, 
                         std::set<ModeTripletInfo>& mtrips);


// **************************************************
}
#endif  
