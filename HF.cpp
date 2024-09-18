#include "HF.h"

Solution solveHF(const Instance& inst, bool warmStart, double timeLimit, bool reflect, bool constrL2) {
	// This function solves an instance using the HierarchyFlow method
	// This function also supports using Reflect for the secondary flow (but this is turned off by default)

	double start = getCPUTime();     // starting time
	Solution sol;                    // initialize Solution struct
	
	// method declaration
	sol.method = "HF";
	if (warmStart) sol.method += "(W)";
	if (reflect) sol.method += "(PartialReflect)";
	if (!constrL2) sol.method += "(no L2 constr)";

	// create primary graphs
	vector<Graph> G(inst.C + 1);
	vector<vector<int>> itemsSubset;
	int totSizePrimary = 0;
	for (int c = 0; c < inst.C; c++) {
		itemsSubset = inst.groupedItems[c]; // select only items of current color
		G[c] = graphConstructionCSP(itemsSubset, inst.W, true); // create subgraph
		totSizePrimary += G[c].A.size();
	}

	// find the terminal vertices/'secondary items' and their max demands
	vector<int> qPerW(inst.W + 1);
	int a, head, tail, idx;
	for (int c = 0; c < inst.C; c++) {				// for every color
		for (int v = 1; v <= inst.W; v++) {			// for every vertex in primary graph G
			if (G[c].V[v]) {						// if the vertex exists
				qPerW[v]++;							// increment the corresponding secondary demand by 1
			}
		}
	}

	vector<vector<int>> wq;
	for (int w = inst.W; w > 0; w--) {
		if (qPerW[w] > 0) wq.push_back({ w,qPerW[w] });
	}

	// print2DVector(wq, "wq", "item");

	// construct the secondary graph
	if (reflect) { G[inst.C] = graphConstructionCSPReflect(wq, inst.W); }
	else { G[inst.C] = graphConstructionCSP(wq, inst.W, false, false); }

	// G[inst.C].print();


	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();                       	// remove Gurobi message
	// env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);         // create a new model

	// declaration of the variables for the model
	vector<vector<GRBVar>> f(inst.C); // primary flow variables for each color
	// for c in C, for a in A_c: f^c_a = # arc a is used
	for (int c = 0; c < inst.C; c++) {
		for (int a = 0; a < G[c].A.size(); a++) {
			idx = G[c].A[a][2];
			if (idx >= 0) { // for item arcs
				f[c].push_back(model.addVar(0, inst.groupedItems[c][idx][1], 0, GRB_INTEGER));
			}
			else { // for feedback arcs to 0
				f[c].push_back(model.addVar(0, inst.B * inst.C, 0, GRB_INTEGER));
			}
		}
	}

	vector<GRBVar> g(G[inst.C].A.size()); // secondary flow variables
	if (not reflect) {
		for (int a = 0; a < G[inst.C].A.size(); a++) {
			idx = G[inst.C].A[a][2];
			if (idx >= 0) { // for item arcs
				g[a] = model.addVar(0, inst.B, 0, GRB_INTEGER);
			}
			//else { // for loss arcs
			//	g[a] = model.addVar(0, inst.B, 0, GRB_INTEGER);
			//}
		}
	}
	else {
		int w;
		for (int a = 0; a < G[inst.C].A.size() - 1; a++) {
			idx = G[inst.C].A[a][2];
			if (idx > 0) { // for item arcs
				g[a] = model.addVar(0, wq[idx - 1][1], 0, GRB_INTEGER);
			}
			else if (idx < 0) { // for reflected arcs
				w = wq[-idx - 1][0];
				g[a] = model.addVar(0, qPerW[w], 0, GRB_INTEGER);
			}
			else { // for loss arcs
				g[a] = model.addVar(0, 2 * inst.B, 0, GRB_INTEGER);
			}
		}
		g[G[inst.C].A.size() - 1] = model.addVar(0, inst.B, 0, GRB_INTEGER);
	}

	model.update();

	// declare and calculate linear expressions
	// primary flows
	vector<vector<GRBLinExpr>> fIn(inst.C, vector<GRBLinExpr>(inst.W + 1, 0));  // the amount of primary flow entering each vertex
	vector<vector<GRBLinExpr>> fOut(inst.C, vector<GRBLinExpr>(inst.W + 1, 0)); // the amount of primary flow leaving each vertex
	vector<vector<GRBLinExpr>> fTypeUsed(inst.C); 	// the number of primary arcs used of each item type 
	vector<GRBLinExpr> fPathsEnded(inst.W + 1);		// the number of primary paths ending at each vertex
	GRBLinExpr nbPrimaryFlows = 0;					// the total outflow from all primary 0 vertices
	for (int c = 0; c < inst.C; c++) {				// loop over all colors
		fTypeUsed[c].resize(inst.groupedItems[c].size(), 0);
		for (int a = 0; a < G[c].A.size(); a++) {   // loop over all arcs
			tail = G[c].A[a][0], head = G[c].A[a][1], idx = G[c].A[a][2]; // decompose the arc
			fIn[c][head] += f[c][a];				// inflow
			fOut[c][tail] += f[c][a];				// outflow
			if (idx >= 0) {							// for item arcs
				fTypeUsed[c][idx] += f[c][a];		// number of items used of certain type
			}
			else {									// for loss arcs
				fPathsEnded[tail] += f[c][a];		// number of paths ending at certain tail
			}
		}
		nbPrimaryFlows += fOut[c][0];
	}

	// secondary flows
	int R;
	GRBLinExpr nbBinsUsed = 0;
	vector<GRBLinExpr> gIn, gOut, gTypeUsed, gInS, gInR;

	if (not reflect) {
		gIn.resize(inst.W + 1, 0);		// the amount of secondary flow entering each vertex
		gOut.resize(inst.W + 1, 0);		// the amount of secondary flow leaving each vertex
		gTypeUsed.resize(wq.size(), 0);	// the number of secondary arcs used of each type
		for (int a = 0; a < G[inst.C].A.size(); a++) {   // loop over all arcs
			tail = G[inst.C].A[a][0], head = G[inst.C].A[a][1], idx = G[inst.C].A[a][2]; // decompose the arc
			gIn[head] += g[a];        // inflow
			gOut[tail] += g[a];       // outflow
			if (idx >= 0) {
				gTypeUsed[idx] += g[a]; // number of items used of certain type
			}
		}
		nbBinsUsed = gOut[0];
	}
	else {
		R = ceil(inst.W / double(2));   // find halfway point
		gInS.resize(R + 1, 0);    		// the amount of secondary flow entering each vertex using standard arcs (item + loss)
		gInR.resize(R + 1, 0);    		// the amount of secondary flow entering each vertex using reflected arcs
		gOut.resize(R + 1, 0);   		// the amount of secondary flow leaving each vertex
		gTypeUsed.resize(wq.size(), 0);	// the number of secondary arcs used of each item type 
		for (int a = 0; a < G[inst.C].A.size() - 1; a++) {   				 // loop over all arcs (treat final arc as special case)
			tail = G[inst.C].A[a][0], head = G[inst.C].A[a][1], idx = G[inst.C].A[a][2]; // decompose the arc
			if (idx > 0) {					// for item arcs
				gInS[head] += g[a];        	// inflow using standard arcs (item + loss)
				gTypeUsed[idx - 1] += g[a];	// number of items used of certain type
			}
			else if (idx == 0) {			// for loss arcs
				gInS[head] += g[a];        	// inflow using standard arcs (item + loss)
			}
			else {							// for reflected arcs
				gInR[head] += g[a];        	// inflow using reflected arcs
				gTypeUsed[-idx - 1] += g[a];// number of items used of certain type
				nbBinsUsed += g[a];   // number of reflected arcs
			}
			gOut[tail] += g[a];       		// outflow
		}
		// the final loss arc
		int aSpecial = G[inst.C].A.size() - 1;
		gInR[R] += g[aSpecial], gOut[R] += g[aSpecial], nbBinsUsed += g[aSpecial];
	}
	model.update();

	// set the objective: minimize the amount of color fragmentation
	model.setObjective(nbPrimaryFlows, GRB_MINIMIZE);

	// constraints 1: primary flow conservation
	for (int c = 0; c < inst.C; c++) {						// loop over all colors
		for (int v = 0; v <= inst.W; v++) {					// loop over all vertices
			if (G[c].V[v]) {
				model.addConstr(fIn[c][v] == fOut[c][v]);	// inflow = outflow
			}
		}
	}

	// constraints 2: item type quantities
	for (int c = 0; c < inst.C; c++) { // loop over all colors
		for (int j = 0; j < inst.groupedItems[c].size(); j++) {				// loop over all item types
			model.addConstr(fTypeUsed[c][j] == inst.groupedItems[c][j][1]);	// demand met
		}
	}

	// constraints 3: limited number of bins
	model.addConstr(nbBinsUsed <= inst.B);

	// constraints 4: secondary flow conservation
	if (not reflect) {
		for (int v = 1; v < inst.W; v++) {     // loop over all vertices
			if (G[inst.C].V[v]) {
				if (gOut[v].size() > 0) { model.addConstr(gIn[v] >= gOut[v]); }
			}
		}
	}
	else {
		for (int v = 1; v <= R; v++) {     // loop over all vertices
			if (G[inst.C].V[v]) {
				model.addConstr(gInS[v] == gInR[v] + gOut[v]); // standard inflow = reflected inflow + outflow
			}
		}
		model.addConstr(gOut[0] == 2 * nbBinsUsed);	// each bin consists of 2 flows, one of which uses a relect arc
	}

	// constraints 5: linking constraints
	for (int j = 0; j < wq.size(); j++) {
		model.addConstr(gTypeUsed[j] == fPathsEnded[wq[j][0]]);
	}

	// (optional) constraints 6: lower bound on number of sub-paths per color
	if (constrL2) {
		for (int c = 0; c < inst.C; c++) {			    // loop over all colors
			model.addConstr(fOut[c][0] >= inst.LBs[c]);	// mimimum number of subpaths
		}
	}

	model.update();

	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	// start with provided feasible solution
	if (warmStart && inst.feasSol.size() > 0) {
		// Find a mapping from the arcs (i.e., the two endpoints) to their indices
		vector<vector<vector<int>>> AN2I(inst.C + 1, vector<vector<int>>(inst.W + 1, vector<int>(inst.W + 1, -1)));
		vector<int> AL2I; 
		if (not reflect) { AL2I.resize(inst.W + 1); }
		else { AL2I.resize(R + 1); }
		vector<vector<int>> AR2I; if (reflect) { AR2I.resize(R + 1, vector<int>(R + 1, -1)); }
		int zeroLossI;

		for (int c = 0; c < inst.C; c++) {				  // loop over all colors
			for (int a = 0; a < G[c].A.size(); a++) {     // loop over all primary arcs in Gc
				tail = G[c].A[a][0], head = G[c].A[a][1]; //idx = G[c].A[a][2]; // find tail and head
				AN2I[c][tail][head] = a;
			}
		}
		for (int a = 0; a < G[inst.C].A.size(); a++) {     // loop over all secondary arcs
			tail = G[inst.C].A[a][0], head = G[inst.C].A[a][1], idx = G[inst.C].A[a][2]; // decompose the arc
			if (not reflect) { // if reflect is not used
				if (idx != -1) AN2I[inst.C][tail][head] = a;
				else AL2I[tail] = a;
			}
			else { // if reflect is used
				if (idx > 0) { AN2I[inst.C][tail][head] = a; }
				else if (idx < 0) { AR2I[tail][head] = a; }
				else if (head == 0) { zeroLossI = a; }
				else { AL2I[tail] = a; }
			}
		}

		int cur_v, vCur, new_v, vNew, w, q, a;
		vector<vector<int>> fInit(inst.C);
		vector<int> subPathLengths;
		vector<vector<int>> orderedSubPathLengths;

		for (int c = 0; c < inst.C; c++) {
			fInit[c].resize(G[c].A.size(), 0);
		}
		vector<int> gInit(G[inst.C].A.size(), 0);
		vector<int> I1, I2;
		int loadI1, loadI2;
		int collidePoint;
		bool reflectUsed, reflectForI1;

		for (int b = 0; b < inst.feasSol.size(); b++) {		// for every bin
			subPathLengths.clear();							// start a new path
			for (int c = 0; c < inst.C; c++) {				// for every color
				if (inst.feasSol[b][c].size() > 0) {		// if the color is used in the bin
					cur_v = 0;
					// find the primary subpath
					for (int j = 0; j < inst.feasSol[b][c].size(); j++) { // for every item type
						w = inst.feasSol[b][c][j][0], q = inst.feasSol[b][c][j][1];
						for (int k = 0; k < q; k++) {	// for every copy of the item type
							new_v = cur_v + w;			// find new value
							a = AN2I[c][cur_v][new_v];	// find 1D index of arc
							fInit[c][a]++;				// add 1 to the flow on that arc
							cur_v = new_v;				// move to the new point
						}
					}
					// use a primary loss arc back to 0
					a = AN2I[c][cur_v][0];
					fInit[c][a]++;  // add 1 to the flow on that arc

					// save the lengths of the subpaths
					subPathLengths.push_back(cur_v);
				}
			}

			// Add the secondary arcs (one for each subpath/color used)
			orderedSubPathLengths = sortAndCount(subPathLengths);

			if (not reflect) { // if reflect is not used
				vCur = 0;
				for (int i = 0; i < orderedSubPathLengths.size(); i++) {	// for every subpath length
					for (int k = 0; k < orderedSubPathLengths[i][1]; k++) { // for all copies of the same length
						vNew = vCur + orderedSubPathLengths[i][0];			// find next point
						a = AN2I[inst.C][vCur][vNew];						// find 1D index of arc
						gInit[a]++;											// add 1 to the flow on that arc
						vCur = vNew;										// move to the new point
					}
				}
				//if (vCur != inst.W) {	// add a loss arc to W if W wasn't reached yet
				//	a = AL2I[vCur];		// find 1D index of the loss arc
				//	gInit[a]++;			// add 1 to the flow on that arc
				//}
			}
			else if (reflect) {			 // if reflect is used
				// divide the items over two half bins
				I1.clear(); I2.clear(); loadI1 = 0; loadI2 = 0;
				for (int i = 0; i < orderedSubPathLengths.size(); i++) {	// for every subpath length
					w = orderedSubPathLengths[i][0];
					for (int k = 0; k < orderedSubPathLengths[i][1]; k++) { // for all copies of the same length
						if (loadI1 >= R) {
							I2.push_back(w); loadI2 += w;
						}
						else if (loadI2 >= R) {
							I1.push_back(w); loadI1 += w;
						}
						else if (loadI1 + w < R || 2 * loadI1 + w <= inst.W) {
							I1.push_back(w); loadI1 += w;
						}
						else {
							I2.push_back(w); loadI2 += w;
						}
					}
				}
				if (loadI1 >= R) {
					reflectForI1 = true;
					collidePoint = inst.W - loadI1;
				}
				else {
					reflectForI1 = false;
					if (loadI2 >= R) { collidePoint = inst.W - loadI2; }
					else { collidePoint = R; }
				}

				// turn the two half bins into flows
				// Halfbin I1
				vCur = 0;
				reflectUsed = false;
				for (int i = 0; i < I1.size(); i++) {
					if (vCur + I1[i] <= R) {
						vNew = vCur + I1[i];
						a = AN2I[inst.C][vCur][vNew];
					}
					else {
						vNew = inst.W - loadI1;
						a = AR2I[vCur][vNew];
						reflectUsed = true;
					}
					gInit[a]++; // add 1 to the flow on that arc
					vCur = vNew; // move to the new point
				}

				if (not reflectForI1 && collidePoint == 0) {
					a = zeroLossI;
					gInit[a]++;
				}

				while (vCur != collidePoint) { // take loss arcs until reaching the collide point
					a = AL2I[vCur]; // take a loss arc
					gInit[a]++; // add 1 to the flow on that arc
					vCur = G[inst.C].A[a][1]; // move to the head of the loss arc
				}
				if (reflectForI1 && not reflectUsed) { // if reflect must still be used for I1
					a = AL2I[R]; // take the final reflected loss arc
					gInit[a]++; // add 1 to the flow on that arc
				}

				// Halfbin I2
				vCur = 0;
				for (int i = 0; i < I2.size(); i++) {
					if (vCur + I2[i] <= R) {
						vNew = vCur + I2[i];
						a = AN2I[inst.C][vCur][vNew];
					}
					else {
						vNew = inst.W - loadI2;
						a = AR2I[vCur][vNew];
						reflectUsed = true;
					}
					gInit[a]++; // add 1 to the flow on that arc
					vCur = vNew; // move to the new point
				}

				if (reflectForI1 && collidePoint == 0) {
					a = zeroLossI;
					gInit[a]++;
				}

				while (vCur != collidePoint) { // take loss arcs until reaching the collide point
					a = AL2I[vCur]; // take a loss arc
					gInit[a]++; // add 1 to the flow on that arc
					vCur = G[inst.C].A[a][1]; // move to the head of the loss arc
				}
				if (not reflectForI1 && not reflectUsed) { // if reflect must still be used for I2
					a = AL2I[R]; // take the final reflected loss arc
					gInit[a]++; // add 1 to the flow on that arc
				}
			}
		}

		// set the initial values
		for (int c = 0; c < inst.C; c++) {
			for (int a = 0; a < G[c].A.size(); a++) {
				f[c][a].set(GRB_DoubleAttr_Start, fInit[c][a]);
			}
		}
		for (int a = 0; a < G[inst.C].A.size(); a++) {
			g[a].set(GRB_DoubleAttr_Start, gInit[a]);
		}

		model.update();
	}

	// find the optimal solution
	model.optimize();

	// store the results in a Solution object
	sol.Nvar = model.get(GRB_IntAttr_NumVars);       // number of variables
	sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); // number of constraints
	sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      // number of non-zero coefficients	
	sol.feas = -1;
	sol.opt = 0;
	sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
	sol.UB = 1000000;

	// if the instance is infeasible
	if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		sol.feas = 0;
	}

	// if a solution has been found within the time limit
	else if (model.get(GRB_IntAttr_SolCount) >= 1) {
		sol.feas = 1;
		sol.UB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		// if the solution is optimal
		if (sol.LB == sol.UB) {
			sol.opt = 1;
		}

		vector<vector<vector<vector<int>>>> subpathsByLength(inst.W + 1);
		vector<vector<int>> colorsOfSubpaths(inst.W + 1);
		vector<int> sizesInBin;
		int fval, end, nFlows;
		vector<vector<vector<int>>> AByTail(inst.W + 1);
		vector<int> a;

		// for each color, decompose all primary subpaths
		for (int c = 0; c < inst.C; c++) {
			// first store the used arcs by tail
			for (int a = 0; a < G[c].A.size(); a++) {
				fval = ceil(f[c][a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
				for (int k = 0; k < fval; k++) {
					AByTail[G[c].A[a][0]].push_back(G[c].A[a]); // add arc to set of arcs with same tail
				}
			}

			// create a subpath one by one
			nFlows = AByTail[0].size();
			for (int b = 0; b < nFlows; b++) {
				tail = 0;						// start at s
				idx = 0;
				while (true) {				// end after using a terminal arc
					a = AByTail[tail].back();	// select an arc
					head = a[1];				// find the head of the arc
					idx = a[2];					// find the index of the arc
					AByTail[tail].pop_back();	// remove the used arc
					if (idx >= 0) {				// for item arcs: 
						sizesInBin.push_back(inst.groupedItems[c][idx][0]); // add item to the bin
						tail = head;										// move to the new point
					}
					else { // after a loss arc
						end = tail; // save the end of the path
						break;
					}
				}

				// add the bin to the set of subpaths with the same end
				subpathsByLength[end].push_back(sortAndCount(sizesInBin));
				colorsOfSubpaths[end].push_back(c);

				// reset the bin
				sizesInBin.clear();
			}

			// reset the tail list
			AByTail.clear();
			AByTail.resize(inst.W + 1);
		}

		// Now combine subpaths using the secondary flow

		// first store the used secondary arcs by tail
		int gval;
		AByTail.clear();
		if (not reflect) { AByTail.resize(inst.W + 1); }
		else { AByTail.resize(R + 1); }
		for (int a = 0; a < G[inst.C].A.size(); a++) {
			gval = ceil(g[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
			for (int k = 0; k < gval; k++) {
				AByTail[G[inst.C].A[a][0]].push_back(G[inst.C].A[a]); // add arc to set of arcs with same tail
			}
		}

		if (not reflect) {
			// store the items assigned to each bin(sorted by color)
			vector<vector<vector<int>>> currentBin(inst.C);
			// combine the subpaths to construct full 0-W paths 
			vector<vector<int>> subpath;
			nFlows = AByTail[0].size();
			int subpathCol, len;
			for (int b = 0; b < nFlows; b++) { // loop over the to-be-filled bins
				tail = 0;						// start at 0
				while (AByTail[tail].size() > 0) {	// while we can continue the path
					a = AByTail[tail].back();	// select an arc
					head = a[1];				// find the head of the arc
					idx = a[2];					// find the index of the arc
					if (idx >= 0) {				// for non-loss arcs: add the subpath to the bin
						len = wq[idx][0];
						subpath = subpathsByLength[len].back();
						subpathCol = colorsOfSubpaths[len].back();
						currentBin[subpathCol] = subpath;
						subpathsByLength[len].pop_back();
						colorsOfSubpaths[len].pop_back();
					}
					AByTail[tail].pop_back();	// remove the used arc
					tail = head;				// move to the new point
				}
				sol.binPacking.push_back(currentBin);		// add the last bin to the set of all bins
				currentBin.clear();							// start a new empty bin
				currentBin.resize(inst.C);
			}
		}
		else { // if reflect is used
			// then fill the half bins one by one by finding a path every iteration
			vector<int> halfBin; // The items packed in the current half bin 
			vector<vector<vector<int>>> standardHalfBinsPerEnd(R + 1); // The half bins not containing a reflected arc per end vertex 
			vector<vector<vector<int>>> reflectedHalfBinsPerEnd(R + 1); // The half bins containing a reflected arc per end vertex
			int tail, head, idx;
			int nBinsUsed = ceil(nbBinsUsed.getValue() - EPSILON); // find the number of bins used;
			bool reflectUsed; // indicates if a reflected arc has been used on the path or not
			for (int b = 0; b < 2 * nBinsUsed; b++) { // loop over the to-be-filled bins
				tail = 0;
				reflectUsed = false;
				while (AByTail[tail].size() > 0 && !reflectUsed) {
					a = AByTail[tail].back();	// select an arc
					AByTail[tail].pop_back();	// remove the used arc
					head = a[1];				// find the head of the arc
					idx = a[2];					// find the index of the arc
					if (idx > 0) {				// for item arcs
						halfBin.push_back(wq[idx - 1][0]);
					}
					else if (idx < 0) { // for reflected arcs
						halfBin.push_back(wq[-idx - 1][0]);
						reflectUsed = true;
					}
					else if (tail == R) { // for the loss/reflected arc from R to R
						reflectUsed = true;
					}
					tail = head;				// move to the new point
				}
				if (reflectUsed) reflectedHalfBinsPerEnd[head].push_back(halfBin);
				else standardHalfBinsPerEnd[head].push_back(halfBin);
				halfBin.clear();
			}
			vector<int> reflectedHalf, standardHalf;
			vector<vector<vector<int>>> currentBin;
			for (int end = 0; end <= R; end++) {
				while (reflectedHalfBinsPerEnd[end].size() > 0) {
					currentBin.clear();							// start a new empty bin
					currentBin.resize(inst.C);
					currentBin.resize(inst.C);
					reflectedHalf = reflectedHalfBinsPerEnd[end].back();
					reflectedHalfBinsPerEnd[end].pop_back();
					standardHalf = standardHalfBinsPerEnd[end].back();
					standardHalfBinsPerEnd[end].pop_back();
					for (const int& p : reflectedHalf) {
						currentBin[colorsOfSubpaths[p].back()] = subpathsByLength[p].back();
						subpathsByLength[p].pop_back();
						colorsOfSubpaths[p].pop_back();
					}
					for (const int& p : standardHalf) {
						currentBin[colorsOfSubpaths[p].back()] = subpathsByLength[p].back();
						subpathsByLength[p].pop_back();
						colorsOfSubpaths[p].pop_back();
					}
					sol.binPacking.push_back(currentBin);
				}
			}
		}

	}

	sol.timeT = getCPUTime() - start;	// save the total time

	// Seperately solve the LP relaxation
	if (not warmStart) {
		GRBModel modelRelaxed = model.relax();
		model.reset(1);
		modelRelaxed.optimize();
		if (modelRelaxed.get(GRB_IntAttr_Status) == GRB_OPTIMAL) { sol.LPrel = modelRelaxed.get(GRB_DoubleAttr_ObjVal); }
		sol.timeLP = getCPUTime() - start - sol.timeT;
	}

	return sol;										 // return the solution
}