#ifndef HOST_GLOBAL_H
#define HOST_GLOBAL_H

#include "latt_field.h"

namespace LaphEnv {

// ***********************************************************************
// *                                                                     *
// *   Class "HostGlobal" is a singleton which stores important data     *
// *   which needs to be persistent across several tasks. In             *
// *   particular, the gauge configuration, the smeared gauge config,    *
// *   the LapH eigenvectors are stored.  Only friends of this class     *
// *   can access its static members.                                    *
// *                                                                     *
// *   Global data (persistent across tasks) (declared in quda_laph.cc)  *
// *                                                                     *
// *   std::vector<LattField> HostGlobal::theGaugeConfig;                *
// *   std::vector<LattField> HostGlobal::theSmearedGaugeConfig;         *
// *   std::vector<LattField> HostGlobal::theLaphEigvecs;                *
// *                                                                     *
// ***********************************************************************


class HostGlobal
{

   HostGlobal() {}

   ~HostGlobal() {}

   static std::vector<LattField> theGaugeConfig;

   static std::vector<LattField> theSmearedGaugeConfig;

   static std::vector<LattField> theLaphEigvecs;

   HostGlobal(const HostGlobal &in) = delete;           // no copy constructor

   HostGlobal(HostGlobal &in) = delete;                 // no copy constructor

   HostGlobal& operator=(const HostGlobal &in) = delete; // not assignable

   HostGlobal& operator=(HostGlobal &in) = delete;       // not assignable
    

   static void clear()
   {
    theGaugeConfig.clear();
    theSmearedGaugeConfig.clear();
    theLaphEigvecs.clear();}

   static void clearGaugeConfig()
   {
    theGaugeConfig.clear();}

   static void clearSmearedGaugeConfig()
   {
    theSmearedGaugeConfig.clear();}

   static bool theGaugeConfigIsSet()
   { return !theGaugeConfig.empty();}

   static bool theSmearedGaugeConfigIsSet()
   { return !theSmearedGaugeConfig.empty();}

   static bool theLaphEigvecsAreSet()
   { return !theLaphEigvecs.empty();}

   friend class GaugeConfigurationHandler;
   friend class GluonSmearingHandler;
   friend class QuarkSmearingHandler;
   friend class QuarkHandler;
   friend class PerambulatorHandler;

};

typedef HostGlobal  GB;


// *************************************************************
}
#endif
