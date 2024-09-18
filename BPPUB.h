#ifndef BPPUB_H
#define BPPUB_H

#include "helper_functions.h"
#include "AFCSP.h"
#include "RM2GIFFRE.h"

Solution solveBPPUB(const Instance& inst, bool stopEarly = true, double timeLimit = 90);

#endif