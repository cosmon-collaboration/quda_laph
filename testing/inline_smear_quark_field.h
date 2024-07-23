#ifndef __INLINE_SMEAR_QUARK_FIELD_H__
#define __INLINE_SMEAR_QUARK_FIELD_H__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "field_smearing_handler.h"

// ******************************************************************
// *                                                                *
// * Driver inline measurement that handles Laph smearing the quark *
// * field.                                                         *
// *                                                                *
// *    Workflow 1: This should be run first in 3D chroma_laph      *
// * to compute the needed low-lying eigenvectors of the            *
// * gauge-smeared covariant Laplacian, producing a set of          *
// * "xxx.time_nn" files. Then this should be re-run in 4D          *
// * chroma_laph to combine the ".time_nn" files into ".level_nn"   *
// * files suitable for the inline task that computes the quark     *
// * sources and sinks. However, it is now possible to compute      *
// * the quark sinks using the "time" files only!!                  *
// *                                                                *
// *    Workflow 2: The Laph eigenvectors on all requested time     *
// * slices can be computed in 4D chroma_laph (see below), then     *
// * written to level files or inserted in TheNamedObjMap.  If      *
// * the file name begins with "NOM_", then TheNamedObjMap is used. *
// *                                                                *
// * XML input for the 3D version must have the form:               *
// *                                                                *
// *  <chroma>                                                      *
// *   <RNG><Seed> ... </Seed></RNG>                                *
// *   <Param>                                                      *
// *    <nrow>12 12 12</nrow>   # lattice size Nx Ny Nz             *
// *    <InlineMeasurements>                                        *
// *    <elem>                                                      *
// *                                                                *
// *     <Name>SMEAR_QUARK_FIELD_TIMESLICES</Name>                  *
// *       <GaugeConfigurationInfo> ... </GaugeConfigurationInfo>   *
// *       <GluonStoutSmearingInfo> ... </GluonStoutSmearingInfo>   *
// *       <SmearedGaugeFileName> ... </SmearedGaugeFileName>       *
// *       <QuarkLaphSmearingInfo> ... </QuarkLaphSmearingInfo>     *
// *       <SmearedQuarkFileStub> ... </SmearedQuarkFileStub>       *
// *       <LaphEigenSolverInfo> ... </LaphEigenSolverInfo>         *
// *       <CalcMinTime>...</CalcMinTime>  (optional, for timings)  *
// *       <CalcMaxTime>...</CalcMaxTime>  (optional)               *
// *                                                                *
// *    </elem>                                                     *
// *    </InlineMeasurements>                                       *
// *   </Param>                                                     *
// *  </chroma>                                                     *
// *                                                                *
// * If <LaphEigenSolverInfo> is omitted above, then the task       *
// * computes an estimate of largest eigenvalue of -Laplacian only. *
// *                                                                *
// * If <CalcMinTime>,<CalcMaxTime> are omitted, then the task      *
// * computes eigenvectors on all time slices from min to max in    *
// * the GaugeConfigurationInfo.                                    *
// *                                                                *
// * XML input for the 4D version must have the form:               *
// *                                                                *
// *  <chroma>                                                      *
// *   <RNG><Seed> ... </Seed></RNG>                                *
// *   <Param>                                                      *
// *    <nrow>12 12 12 96</nrow>   # lattice size Nx Ny Nz Nt       *
// *    <InlineMeasurements>                                        *
// *    <elem>                                                      *
// *                                                                *
// *     <Name>SMEAR_QUARK_FIELD_TIMESLICES</Name>                  *
// *       <GaugeConfigurationInfo> ... </GaugeConfigurationInfo>   *
// *       <GluonStoutSmearingInfo> ... </GluonStoutSmearingInfo>   *
// *       <QuarkLaphSmearingInfo> ... </QuarkLaphSmearingInfo>     *
// *       <SmearedQuarkFileStub> ... </SmearedQuarkFileStub>       *
// *       <LaphEigenSolverInfo>...</LaphEigenSolverInfo>(optional) *
// *       <CalcMinTime>...</CalcMinTime>  (optional, for debugging)*
// *       <CalcMaxTime>...</CalcMaxTime>  (optional)               *
// *                                                                *
// *    </elem>                                                     *
// *    </InlineMeasurements>                                       *
// *   </Param>                                                     *
// *  </chroma>                                                     *
// *                                                                *
// * If <LaphEigenSolverInfo> is omitted, 4D chroma_laph combines   *
// * the ".time_nn" files into ".level_nn" files suitable for the   *
// * inline task that computes the quark sources and sinks.         *
// * If <LaphEigenSolverInfo> is included, then Laph eigenvectors   *
// * on all requested time slices are simultaneously computed.      *
// *                                                                *
// * If <CalcMinTime>,<CalcMaxTime> are omitted, then the task      *
// * computes eigenvectors on all time slices from min to max in    *
// * the GaugeConfigurationInfo.                                    *
// *                                                                *
// ******************************************************************



    //  The name of this inline measurement.   This is the name 
    //  with which the createMeasurement function is associated in the 
    //  Object Factory. You must include this name in the XML input 
    //  to Chroma through
    //     <InlineMeasurements>
    //        <elem>
    //            <Name> SMEAR_GAUGE_FIELD_TIMESLICES </Name>
    //             ...
    //        </elem>
    //    </InlineMeasurements>


namespace Chroma { 

  namespace InlineSmearQuarkFieldEnv {

 // **************************************************************

extern const std::string name;
bool registerAll();

    /*! \ingroup inlinehadron */

class SmearQuarkFieldInlineMeas : public AbsInlineMeasurement 
{

   XMLReader xml_rd;   // holds the XML input for this inline
                        // measurement, for use by the operator()
                        // member below

 public:

   SmearQuarkFieldInlineMeas(XMLReader& xml_in, const std::string& path) 
                              : xml_rd(xml_in, path) {}

   ~SmearQuarkFieldInlineMeas() {}
      
      //! Do the measurement
   void operator()(const unsigned long update_no, XMLWriter& xmlout); 

   unsigned long getFrequency() const {return 1;}
   
};
	

// ***********************************************************
  }
}
#endif
