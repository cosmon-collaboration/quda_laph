#ifndef LAPH_NOISE_INFO_H
#define LAPH_NOISE_INFO_H

#include "gauge_configuration_info.h"
#include <cstdint>

namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   An object of class "LaphNoiseInfo" stores identifying info    *
// *   about a quark source noise in the Laph subspace.  The         *
// *   notion of a noise identity extends to an entire ensemble.     *
// *   The XML input must have the format                            *
// *                                                                 *
// *         <LaphNoiseInfo>                                         *
// *            <ZNGroup> 4 </ZNGroup>                               *
// *            <Seed> 315 </Seed>                                   *
// *         </LaphNoiseInfo>                                        *
// *                                                                 *
// *   The group Z(N) is used, so the value of "N" must be specified *
// *   in the tag named "ZNGroup".  Only the values 4,8,32 are       *
// *   currently supported (1 is allowed for debugging purposes).    *
// *                                                                 *
// *   "Seed" must be an integer with value in the range 0 to 65535  *
// *   (which is 2^16-1).                                            *
// *                                                                 *
// *   A 32-bit Mersenne Twister is used to generate the Laph noise. *
// *   It was found that the linear congruential generator in        *
// *   chroma/qdp++ is not adequate for generating the noise and     *
// *   leads to seriously incorrect results in some instances.       *
// *   This class assumes that the HMC trajectory number lies        *
// *   between 0 and 2^16-1 = 65535 (otherwise an error occurs).     *
//                                                                   *
// *   The 32-bit Mersenne Twister needs an unsigned 32-bit integer  *
// *   as a seed.  For a given HMC trajectory number "k", this       *
// *   class constructs an MT seed as follows:                       *
// *                                                                 *
// *    (8 least sig bits of "Seed") (8 most sig bits of "k")        *
// *       (8 most sig bits of "Seed") (8 least sig bits "k")        *
// *                                                                 *
// *                                                                 *
// *   Example usage:                                                *
// *                                                                 *
// *     XMLHandler xml_in(...);                                      *
// *     LaphNoiseInfo rho(xml_in);                                  *
// *                                                                 *
// *     LaphNoiseInfo rho2(....);                                   *
// *     rho.checkEqual(rho2);   // throws exception if rho2 != rho  *
// *     if (rho==rho2) ...      // returns boolean                  *
// *                                                                 *
// *     GaugeConfigurationInfo G(....)                              *
// *     Seed s = rho.getSeed(G);                                    *
// *                                                                 *
// *     string str = rho.output();   // xml string                  *
// *                                                                 *
// *******************************************************************

   // internal representation: right 6 bits are ZN group, other bits
   // are the seed

class LaphNoiseInfo
{

  uint32_t store;

 public:  

  LaphNoiseInfo(); 

  LaphNoiseInfo(const XMLHandler& xml_in);

  LaphNoiseInfo(int znGroup, int seed);

  LaphNoiseInfo(const LaphNoiseInfo& in);

  LaphNoiseInfo& operator=(const LaphNoiseInfo& in);

  ~LaphNoiseInfo(){}

  void checkEqual(const LaphNoiseInfo& in) const;

  bool operator==(const LaphNoiseInfo& in) const;

  bool operator!=(const LaphNoiseInfo& in) const;

  bool operator<(const LaphNoiseInfo& in) const;


    // output functions

  uint32_t getSeed(const GaugeConfigurationInfo& G) const;

  unsigned int getZNGroup() const 
   { const unsigned int GP=0x3Fu; return (store & GP); }

  uint32_t getSeed() const { return (store>>6); }

  std::string output(int indent = 0) const;

  std::string str() const;

  void output(XMLHandler& xmlout) const;

  std::string getHeader() const { return str();}

  void getHeader(XMLHandler& xmlout) const { output(xmlout);}

 private:

  void extract_info_from_reader(XMLHandler& xml_in);

  void encode(int znGroup, int seed);

  void check_assignment(int znGroup, int seed);

  friend class MesonHandler;
  friend class CurrentHandler;
  friend class BaryonHandler;
  friend class MesonInternalLoopHandler;
  friend class TetraquarkHandler;
  friend class TetraquarkOneInternalLoopHandler;
  friend class TetraquarkTwoInternalLoopHandler;

};



// **************************************************
}
#endif
