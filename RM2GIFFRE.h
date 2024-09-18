#ifndef RM2GIFFRE_H
#define RM2GIFFRE_H

#include "helper_functions.h"
#include "LF.h"

struct GraphRM2GIFF {
	int Nvert, Narcs;		   // number of vertices and arcs
	vector<vector<vector<int>>> V_N2I; // set of vertices: [V_N2I]_jdg gives the 1D index for node (j,d,g)
	vector<vector<int>> V_I2N; // set of vertices: [V_I2N]_idx gives the 3D index for the idx'th node
	vector<vector<int>> A;	   // set of arcs: (1D index of tail, 1D index of head, type)
	void print();
};

GraphRM2GIFF graphConstructionRM2GIFF(const Instance& inst, bool terminal=false);
GraphRM2GIFF graphConstructionRM2GIFFUnitLoss(const Instance& inst);
Solution solveRM2GIFFRE(const Instance& inst, bool warmStart=false, double timeLimit=1800, bool terminal = false, bool unitLoss = false);
Solution solveRM2GIFFSetupMehrani(const Instance& inst, bool warmStart=false, double timeLimit=1800);
#endif