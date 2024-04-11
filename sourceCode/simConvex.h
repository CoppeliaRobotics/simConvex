#pragma once

#include <simLib/simTypes.h>
#include <simLib/simExp.h>

SIM_DLLEXPORT int simInit(SSimInit*);
SIM_DLLEXPORT void simCleanup();
SIM_DLLEXPORT void simMsg(SSimMsg*);

SIM_DLLEXPORT void simHACD(void*);
SIM_DLLEXPORT void simVHACD(void*);
SIM_DLLEXPORT void simQhull(void*);
