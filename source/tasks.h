#ifndef TASKS_H
#define TASKS_H

#include "xml_handler.h"

// ******************************************************************
// *                                                                *
// *   The various tasks that quda_laph can do are declared here.   *
// *   Consult the individual cc files for each task for more       *
// *   information on each task and the required XML input format   *
// *                                                                *
// ******************************************************************


namespace LaphEnv {

// ************************************************


void doSmearGaugeField(XMLHandler& xmltask);

void doSmearQuarkField(XMLHandler& xmltask);

void doLaphQuarkLineEnds(XMLHandler& xmltask);

void doLaphQuarkPerambulators(XMLHandler& xmltask);

void doLaphCheckPerambulators(XMLHandler& xmltask);

void doLaphMergePerambulators(XMLHandler& xmltask);

//void doLaphSmearTune(XMLHandler& xmltask);


// ***********************************************************
}
#endif
