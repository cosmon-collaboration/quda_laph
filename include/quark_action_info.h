#ifndef LAPH_QUARK_ACTION_INFO_H
#define LAPH_QUARK_ACTION_INFO_H

#include "quda.h"
#include "xml_handler.h"

namespace LaphEnv {

//  Class "QuarkActionInfo" holds information about the valence     *
//  quark action.  String, integer, and real data values are stored *
//  in "svalues", "ivalues", and "rvalues", respectively.  These    *
//  are each a std::vector of string, integer, and doubles.         *
//  The identifying name of the action is stored in svalues[0].     *
//  There should be an object of this class for the light "ud"      *
//  quarks, and one for the strange "s" quarks. Supported actions   *
//  and their required XML input are shown below.                   *
//                                                                  *
//    <QuarkActionInfo>                                             *
//       <Name> WILSON_CLOVER </Name>  (mandatory)                  *
//       <Flavor> ud </Flavor>  (or "s")                            *
//       <Mass> 0.321 </Mass> (or <Kappa>..</Kappa>)                *
//       <Anisotropy> 3.5 </Anisotropy> (default 1)                 *
//       <clovCoeffSS>2.0668577</clovCoeffSS>                       *
//       <clovCoeffST>2.0668577</clovCoeffST>                       *
//         or just <clovCoeff>..</clovCoeff> is isotropic           *
//       <Tadpole>0.556</Tadpole> (default 1)                       *
//       <TimeBC> antiperiodic </TimeBC> (or periodic)              *
//    </QuarkActionInfo>                                            *

class QuarkActionInfo {

  std::vector<std::string> svalues;

  std::vector<int> ivalues;

  std::vector<double> rvalues;

public:
  QuarkActionInfo(const XMLHandler &xmlin);

  QuarkActionInfo(const QuarkActionInfo &rhs);

  QuarkActionInfo &operator=(const QuarkActionInfo &rhs);

  ~QuarkActionInfo() {}

  std::string getName() const { return svalues[0]; }

  bool isFermionTimeBCAntiPeriodic() const { return ivalues[1] == 0; }

  double getSolutionRescaleFactor() const { return rvalues[0]; }

  void checkEqual(const QuarkActionInfo &rhs) const;

  bool operator==(const QuarkActionInfo &rhs) const;

  std::string output(int indent = 0) const;

  void output(XMLHandler &xmlout) const;

  void setQudaInvertParam(QudaInvertParam &invParam) const;

  const std::vector<double> &getRValues() const { return rvalues; }

private:
  void set_info_wilson_clover(XMLHandler &xmlr);

  void output_wilson_clover(XMLHandler &xmlout) const;

  void setQudaInvertParam_wilson_clover(QudaInvertParam &invParam) const;
};
} // namespace LaphEnv
#endif
