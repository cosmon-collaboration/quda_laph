#ifndef __INLINE_QUARK_LINE_ENDS_H__
#define __INLINE_QUARK_LINE_ENDS_H__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "quark_handler.h"
#include "xml_help.h"

// ******************************************************************
// *                                                                *
// * Driver inline measurement that constructs the quark line ends  *
// * (sinks only).  Assumes that the eigenvectors of the            *
// * gauge-smeared covariant Laplacian are available.  Must be run  *
// * in 4D chroma_laph.  XML input must have the form:              *
// *                                                                *
// *  <chroma>                                                      *
// *   <RNG><Seed> ... </Seed></RNG>                                *
// *   <Cfg> ... </Cfg>                                             *
// *   <Param>                                                      *
// *    <nrow>12 12 12 96</nrow>   # lattice size Nx Ny Nz Nt       *
// *    <InlineMeasurements>                                        *
// *    <elem>                                                      *
// *                                                                *
// *     <Name> LAPH_QUARK_LINE_ENDS </Name>                        *
// *     <QuarkLineEndInfo>                                         *
// *       <GaugeConfigurationInfo> ... </GaugeConfigurationInfo>   *
// *       <GluonStoutSmearingInfo> ... </GluonStoutSmearingInfo>   *
// *       <QuarkLaphSmearingInfo> ... </QuarkLaphSmearingInfo>     *
// *       <SmearedQuarkFileStub> ... </SmearedQuarkFileStub>       *
// *       <LaphDilutionScheme> ... </LaphDilutionScheme>           *
// *       <QuarkActionInfo> ... </QuarkActionInfo>                 *
// *       <FileListInfo> ... </FileListInfo>                       *
// *       <InvertParam> ... </InvertParam>                         *
// *       <Verbosity> ... </Verbosity>  (optional)                 *
// *     </QuarkLineEndInfo>                                        *
// *                                                                *
// *     <SinkComputations> ... </SinkComputations>                 *
// *                                                                *
// *    </elem>                                                     *
// *    </InlineMeasurements>                                       *
// *   </Param>                                                     *
// *  </chroma>                                                     *
// *                                                                *
// * The sink computations to be carried out must be specified as   *
// * follows:                                                       *
// *                                                                *
// *   <NoiseList_TimeProjIndexList>                                *
// *     <TimeProjIndexList>                                        *
// *        <Values>3 5 9</Values>  or  </All>                      *
// *     </TimeProjIndexList>                                       *
// *     <LaphNoiseList>                                            *
// *        <LaphNoiseInfo> ... </LaphNoiseInfo>                    *
// *       ... (other LaphNoiseInfo tags)                           *
// *     </LaphNoiseList>                                           *
// *   </NoiseList_TimeProjIndexList>                               *
// *                                                                *
// *   <ComputationList>                                            *
// *     <Computation>                                              *
// *       <LaphNoiseInfo>...</LaphNoiseInfo>                       *
// *       <TimeProjIndex>...</TimeProjIndex>                       *
// *     </Computation>                                             *
// *      ...  (other Computation tags)                             *
// *   </ComputationList>                                           *
// *                                                                *
// * If the tag "Verbosity" is included with value "full", then the *
// * quark sink solutions will be echoed to standard output.        * 
// *                                                                *
// ******************************************************************


namespace Chroma {
                 
#if (QDP_ND == 4) 

  namespace InlineStochLaphQuarkEnv {

 // **************************************************************


extern const std::string name;
bool registerAll();

    /*! \ingroup inlinehadron */

class StochLaphQuarkInlineMeas : public AbsInlineMeasurement 
{

   XMLReader xml_rd;   // holds the XML input for this inline
                       // measurement, for use by the operator()
                       // member below

   struct SinkComputation {
      LaphEnv::LaphNoiseInfo Noise;
      int TimeProjIndex;
      bool postSink; 

    SinkComputation(const LaphEnv::LaphNoiseInfo& in_noise, int in_time_proj_index, bool in_post_sink)
       : Noise(in_noise), TimeProjIndex(in_time_proj_index), postSink(in_post_sink) {}
   };

   std::list<SinkComputation> sinkComputations;

 public:

   StochLaphQuarkInlineMeas(XMLReader& xml_in, const std::string& path) 
                              : xml_rd(xml_in, path) {}

   ~StochLaphQuarkInlineMeas() {}
      
   void setSinkComputations(int NumTimeDilutionProjs);
   void clearSinkComputations();

      //! Do the measurement
   void operator()(const unsigned long update_no, XMLWriter& xmlout); 

   unsigned long getFrequency() const {return 1;}
   
};
	

// ***********************************************************
  }
#endif
}

#endif
