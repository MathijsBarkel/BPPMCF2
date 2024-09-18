#include "BPPUB.h"

Solution solveBPPUB(const Instance& inst, bool stopEarly, double timeLimit) {
	// This function applies the algorithm BPP-UB to the given instance

	double start = getCPUTime();   // starting time
	Solution sol;				   // initialize Solution struct
	sol.method = "BPP-UB";		   // method declaration
	if (not stopEarly) { sol.method += "(NoEarlyStop)"; }

	// set a random seed
	srand(53);                   // option 1: every time the same
	// srand((unsigned)time(NULL));	// option 2: every time different

	// set the number of tries that the second phase is attempted
	int nTries = 100;

	vector<vector<int>>	wq = inst.colorlessItems; // combine all items (ignoring color)

	// Create graph
	Graph G = graphConstructionCSP(wq, inst.W);
	// G.print(); // print the graph

	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();							// remove Gurobi message
	env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);         // create a new model

	// declaration of the variables for the model
	vector<GRBVar> f(G.A.size());

	// declaration of the variables for the model
	// for a in A: fa = # arc a is used
	int idx;
	for (int a = 0; a < G.A.size(); a++) {
		idx = G.A[a][2];
		if (idx != -1) {
			f[a] = model.addVar(0, wq[idx][1], 0, GRB_INTEGER);
		}
		else {
			f[a] = model.addVar(0, inst.B, 0, GRB_INTEGER);
		}
	}
	model.update();

	// declare linear expressions
	vector<GRBLinExpr> fIn(inst.W + 1, 0);    		// the amount of flow entering each vertex
	vector<GRBLinExpr> fOut(inst.W + 1, 0);   		// the amount of flow leaving each vertex
	vector<GRBLinExpr> typeUsed(wq.size(), 0); 	// the amount of arcs used of each item type 

	// calculate the linear expressions	
	for (int a = 0; a < G.A.size(); a++) {   				// loop over all arcs
		fIn[G.A[a][1]] += f[a];        						// inflow
		fOut[G.A[a][0]] += f[a];       						// outflow
		if (G.A[a][2] >= 0) typeUsed[G.A[a][2]] += f[a];	// number of items used of certain type
	}
	model.update();

	// set the objective: minimize the number of flow going out of 0
	model.setObjective(fOut[0], GRB_MINIMIZE);
	// constraints 1: flow conservation
	for (int v = 1; v < inst.W; v++) {       	// loop over all vertices
		if (G.V[v])
			model.addConstr(fIn[v] == fOut[v]); // inflow = outflow
	}
	model.addConstr(fOut[0] == fIn[inst.W]);	// total flow
	// constraints 2: item type quantities
	for (int j = 0; j < wq.size(); j++) {			// loop over all item types
		model.addConstr(typeUsed[j] == wq[j][1]);	// demand met
	}

	model.update();
	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	// (optionally:) stop early if the required number of bins is reached
	if (stopEarly) { model.getEnv().set(GRB_DoubleParam_BestObjStop, inst.B); }

	// find the optimal solution
	model.optimize();

	// store the results in a Solution object
	sol.Nvar = model.get(GRB_IntAttr_NumVars);       // number of variables
	sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); // number of constraints
	sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      // number of non-zero coefficients	
	sol.feas = -1;
	sol.opt = false;
	sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);

	// if the instance is infeasible
	if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		sol.feas = 0;
		sol.opt = false;
	}
	// if a solution has been found within the time limit
	else if (model.get(GRB_IntAttr_SolCount) >= 1) {
		int nBinsUsed = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		// if the LB is more than the amount of bins --> infeasible instance
		if (sol.LB > inst.B) {
			sol.feas = 0;
			sol.opt = false;
		}
		// if time ran out and the LB is at most the amount of bins, but the UB is more, we don't know if the instance is feasible
		else if (nBinsUsed > inst.B) {
			sol.feas = -1;
			sol.opt = false;
		}
		// if the solution is optimal for the CSP instance and the amount of bin used is at most the amount of bins available
		else if (nBinsUsed <= inst.B) {
			// decompose the solution in a smart way to find the amount of color fragmentation

			sol.feas = 1;
			// first store the used arcs by tail
			int fval;
			vector<vector<vector<int>>> AByTailOriginal(inst.W + 1);
			for (int a = 0; a < G.A.size(); a++) {
				fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
				if (fval > 0) {
					AByTailOriginal[G.A[a][0]].push_back({ G.A[a][1], G.A[a][2], fval }); // add arc to set of arcs with same tail
				}
			}

			// for every type of arc we want the available colors and their quantities
			// (as an intermediate step: every arc type belongs to an item size, which belongs to some colors)
			// colorsPerSize_wc gives the number of items of size w and color c
			vector<vector<int>> colorsPerSizeOriginal(inst.W + 1, vector<int>(inst.C, 0));
			// inst.countedItems contains for every item type the color, size and quantity (c,w,q)
			for (int j = 0; j < inst.countedItems.size(); j++) {
				colorsPerSizeOriginal[inst.countedItems[j][1]][inst.countedItems[j][0]] = inst.countedItems[j][2];
			}
			// colorsPerType_jc gives the number of items of index j and color c
			vector<vector<int>> colorsPerTypeOriginal(wq.size());
			for (int j = 0; j < wq.size(); j++) {
				colorsPerTypeOriginal[j] = colorsPerSizeOriginal[wq[j][0]];
			}

			// declare variables
			vector<vector<vector<int>>> AByTail;
			vector<vector<int>> colorsPerType;
			vector<vector<vector<vector<int>>>> currentPacking, bestPacking;
			vector<vector<vector<int>>> currentBin;
			vector<vector<int>> goodCombinations, badCombinations;
			vector<int> a;
			int tail, head, idx, loss, chosenComb, chosenColor, chosenArc;
			vector<bool> usedColors;
			vector<vector<int>> itemsInBinPerColor;
			int curVal, bestVal = int('inf');

			// try decomposing the flow nTries times
			for (int t = 0; t < nTries; t++) {
				// reset the solution
				currentPacking.clear();
				curVal = 0;
				AByTail = AByTailOriginal;
				colorsPerType = colorsPerTypeOriginal;

				// decompose the flow one-by-one
				for (int b = 0; b < nBinsUsed; b++) { // loop over the to-be-filled bins
					usedColors.clear(); usedColors.resize(inst.C, false);
					itemsInBinPerColor.clear(); itemsInBinPerColor.resize(inst.C);
					currentBin.clear(); currentBin.resize(inst.C);
					tail = 0;				 // start a new path at 0
					while (tail != inst.W) { // while the terminanal vertex hasn't been reached
						loss = -1;			 // contains index of loss arc from current tail, or -1 if no such arc exists
						goodCombinations.clear(); badCombinations.clear();
						for (int a_idx = 0; a_idx < AByTail[tail].size(); a_idx++) { // check each possible outgoing arc
							idx = AByTail[tail][a_idx][1];			 // find the index of the arc
							if (idx != -1) {						 // if the arc is not a loss arc
								for (int c = 0; c < inst.C; c++) {   // check the possible colors
									if (colorsPerType[idx][c] > 0) { // if the color is available
										if (usedColors[c]) {		 // if the color was already used
											goodCombinations.push_back({ a_idx, c }); // add to good combinations
										}
										else {						 // if the color wasn't used yet 
											badCombinations.push_back({ a_idx, c }); // add to bad combinations
										}
									}
								}
							}
							else {				// if the arc is a loss arc
								loss = a_idx;	// save the index of the arc
							}
						}

						// if a good move is possible (i.e., using an arc of an already used color)
						if (goodCombinations.size() > 0) {
							chosenComb = rand() % goodCombinations.size();
							chosenArc = goodCombinations[chosenComb][0];
							chosenColor = goodCombinations[chosenComb][1];
						}
						// otherwise, use a loss arc if possible
						else if (loss != -1) {
							chosenArc = loss;
							chosenColor = -1;
						}
						// otherwise, make a bad move (i.e., using an arc of a new color)
						else {
							chosenComb = rand() % badCombinations.size();
							chosenArc = badCombinations[chosenComb][0];
							chosenColor = badCombinations[chosenComb][1];
						}

						// take the chosen arc + remove afterwards
						head = AByTail[tail][chosenArc][0];
						idx = AByTail[tail][chosenArc][1];
						AByTail[tail][chosenArc][2]--;
						if (AByTail[tail][chosenArc][2] == 0) {
							AByTail[tail].erase(AByTail[tail].begin() + chosenArc);
						}
						tail = head;	// move to the new point

						if (chosenColor != -1) {
							usedColors[chosenColor] = true;    // update list containing colors used in path
							colorsPerType[idx][chosenColor]--; // remove one arc of that size and color
							itemsInBinPerColor[chosenColor].push_back(wq[idx][0]); // add the item
						}

					}

					for (int c = 0; c < inst.C; c++) {			// for every color
						if (itemsInBinPerColor[c].size() > 0) { // if the color is contained in the bin
							curVal++;							// add 1 to the objective value
							currentBin[c] = sortAndCount(itemsInBinPerColor[c]); // group and count the items of the current color in the bin
						}
					}
					currentPacking.push_back(currentBin);		// add the last bin to the set of all bins
				}

				// update the best solution
				if (curVal < bestVal) {
					bestVal = curVal;
					bestPacking = currentPacking;
					// terminate if the solution is proven optimal
					if (bestVal == inst.LBtot) { break; }
				}
			}
			// save the best solution found
			sol.binPacking = bestPacking;
			sol.UB = bestVal;
			sol.LB = inst.LBtot;				// the lower bound is simply set to the L2 bound
			if (sol.UB == sol.LB) { sol.opt = true; }
		}
	}

	sol.timeT = getCPUTime() - start;	// save the total time
	return sol;							// return the solution
}