#ifndef AFCSP_H
#define AFCSP_H

#include "helper_functions.h"

struct Graph {
	vector<bool> V;			// set of vertices
	vector<vector<int>> A;	// set of arcs
	void print();
};

Graph graphConstructionCSP(const vector<vector<int>>& wq, int W, bool terminalZero = false, bool useLoss = true);
Solution solveAFCSP(const Instance& inst, const vector<vector<int>>& wq, int L2 = 0, int c = 0, int early = 0);
Solution solveAFCSPConcentration(const Instance& inst, const vector<vector<int>>& wq, int L2 = 0, int c = 0, bool startWithReflect= true, int searchStrategy = 0, bool lexicographic = false);
Graph graphConstructionCSPReflect(const vector<vector<int>>& wq, int W);
Solution solveReflectCSP(const Instance& inst, const vector<vector<int>>& wq, int L2 = 0, int c = 0, int early = 0);
#endif