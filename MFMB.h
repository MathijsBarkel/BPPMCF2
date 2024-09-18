#ifndef MFMB_H
#define MFMB_H

#include "helper_functions.h"
#include "AFCSP.h"

Solution solveMFMB(const Instance& inst, const Solution& solUB, bool warmStart = false, double timeLimit = 1800, bool reflect = true, bool constrL2 = true, bool constrSym = true);

#endif