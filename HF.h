#ifndef HF_H
#define HF_H

#include "helper_functions.h"
#include "AFCSP.h"
#include<memory>

Solution solveHF(const Instance& inst, bool warmStart = false, double timeLimit = 1800, bool reflect = false, bool constrL2 = true);

#endif
