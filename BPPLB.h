#ifndef BPPLB_H
#define BPPLB_H

#include "helper_functions.h"
#include "AFCSP.h"

Solution solveBPPLB(const Instance& inst, vector<int>& LBs, bool concentration = true, bool reflect = true, bool lexicographic = false);

#endif
