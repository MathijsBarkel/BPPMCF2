#ifndef LF_H
#define LF_H

#include "helper_functions.h"
#include "RM2GIFFRE.h"

struct Graph2D {
	int Nvert, Narcs;		   // number of vertices and arcs
	vector<vector<int>> V_N2I; // set of vertices: [V_N2I]_vc gives the 1D index for node v_c
	vector<vector<int>> V_I2N; // set of vertices: [V_I2N]_idx gives the 2D index for the idx'th node
	vector<vector<int>> A;	   // set of arcs: (1D index of tail, 1D index of head, type)
	void print();
};

vector<int> determineColorOrder(const Instance& inst);
Graph2D graphConstructionLF(const Instance& inst, int transitionType);
Solution solveLF(const Instance& inst, bool warmStart = false, double timeLimit = 1800, int transitionType = 0, bool constrL2 = true, bool constrSym = true);
Solution solveModelAfterOrdering(Instance inst, bool warmStart = false, double timeLimit = 1800, string model = "LF", int transitionType = 0, bool constrL2 = true, bool constrSym = true);
#endif
