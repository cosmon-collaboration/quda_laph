#ifndef TASKS_H
#define TASKS_H

#include "xml_handler.h"

// ******************************************************************
// *                                                                *
// * Driver inline measurement that stout smears the gauge field    *
// * in 4D chroma_laph, and either outputs the results in 3D or     *
// * 4D format.  3D format can only go to a file, whereas 4D        *
// * format can either go to a file or TheNamedObjMap (as usual     *
// * by starting the file name with "NOM_").  XML input must        *
// * have the form:                                                 *
// *                                                                *
// *  <chroma>                                                      *
// *   <RNG><Seed> ... </Seed></RNG>                                *
// *   <Cfg> ... </Cfg>                                             *
// *   <Param>                                                      *
// *    <nrow>12 12 12 96</nrow>   # lattice size Nx Ny Nz Nt       *
// *    <InlineMeasurements>                                        *
// *    <elem>                                                      *
// *                                                                *
// *     <Name>SMEAR_GAUGE_FIELD</Name>                             *
// *       <GaugeConfigurationInfo> ... </GaugeConfigurationInfo>   *
// *       <GluonStoutSmearingInfo> ... </GluonStoutSmearingInfo>   *
// *       <SmearedGaugeFileName> ... </SmearedGaugeFileName>       *
// *       <OutputFormat>3D</OutputFormat>  (or 4D) (default 3D)    *
// *                                                                *
// *    </elem>                                                     *
// *    </InlineMeasurements>                                       *
// *   </Param>                                                     *
// *  </chroma>                                                     *
// *                                                                *
// ******************************************************************


namespace LaphEnv {

// ************************************************


void doSmearGaugeField(XMLHandler& xmltask);

void doSmearQuarkField(XMLHandler& xmltask);

void doLaphQuarkLineEnds(XMLHandler& xmltask);

void doLaphQuarkPerambulators(XMLHandler& xmltask);



// ***********************************************************
}
#endif
