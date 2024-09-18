#ifndef TS_H
#define TS_H

#include "helper_functions.h"
#include "BPPLB.h"

Solution solveTS(const Instance& inst, Solution solBPPUB, int tabuTime = 50, int focusIter = 40, int maxNIterNoImprovement = 5000, bool diversification = true, int maxCFIncrease = 1000000);
void printSolution(const vector<vector<vector<int>>>& curPacking, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor, int iter);
double findSuperSwapScore(int cX, int bX, int cY, int bY, bool lossConcentration, int W, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor);
double findSwapScore(int iX, int iY, bool lossConcentration, int W, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor);
double findSuperTransferScore(int cX, int bX, int bY, bool lossConcentration, int W, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor);
double findTransferScore(int iX, int bY, bool lossConcentration, int W, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor);
#endif