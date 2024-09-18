#include "MFMB.h"

Solution solveMFMB(const Instance& inst, const Solution& solUB, bool warmStart, double timeLimit, bool reflect, bool constrL2, bool constrSym) {
	// This function solves an instance using the MonoFlow-MultiBin method, possibly using Reflect
	double start = getCPUTime();     // starting time
	Solution sol;                    // initialize Solution struct

	// method declaration
	sol.method = "MFMB";			     
	if (warmStart) sol.method += "(W)";
	if (!reflect) sol.method += "(NoReflect)";
	
	// find a valid upper bound on the required number of multichromatic bins
	int Bplus = min(inst.B, max(solUB.UB - inst.B,0));

	// create a graph for each color
	vector<Graph> G(inst.C);
	vector<vector<int>> itemsSubset;
	for (int c = 0; c < inst.C; c++) { // loop over all colors
		if (reflect) { G[c] = graphConstructionCSPReflect(inst.groupedItems[c], inst.W); }
		else { G[c] = graphConstructionCSP(inst.groupedItems[c], inst.W); }
	}

	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();                       	// remove Gurobi message
	// env.set(GRB_IntParam_LogToConsole, 0);	// turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);			// create a new model

	// declaration of the variables for the model

	// flow variables for mono-chromatic bins
	vector<vector<GRBVar>> f(inst.C); // flow variables for each color
	int idx;
	// for c in C, for a in A_c: f^c_a = # arc a is used
	for (int c = 0; c < inst.C; c++) {				// loop over all colors
		for (int a = 0; a < G[c].A.size(); a++) {	// loop over all arcs in G[c]
			idx = G[c].A[a][2];						// find the index of the arc
			if (reflect) {							// if reflect is used
				if (a == G[c].A.size() - 1) {		// the special final arc
					f[c].push_back(model.addVar(0, solUB.UB, 0, GRB_INTEGER));
				}
				else if (idx > 0) {					// for regular item arcs
					f[c].push_back(model.addVar(0, inst.groupedItems[c][idx - 1][1], 0, GRB_INTEGER));
				}
				else if (idx < 0) {					// for reflected arcs
					f[c].push_back(model.addVar(0, inst.groupedItems[c][-idx - 1][1], 0, GRB_INTEGER));
				}
				else {								// for loss arcs
					f[c].push_back(model.addVar(0, 2 * solUB.UB, 0, GRB_INTEGER));
				}
			}
			else {									// if reflect is not used
				if (idx >= 0) {						// for item arcs
					f[c].push_back(model.addVar(0, inst.groupedItems[c][idx][1], 0, GRB_INTEGER));
				}
				else {								// for loss arcs
					f[c].push_back(model.addVar(0, solUB.UB, 0, GRB_INTEGER));
				}
			}
		}
	}

	// traditional variables for multi-chromatic bins
	vector<vector<GRBVar>> x(inst.J, vector<GRBVar>(Bplus + 1));
	vector<vector<GRBVar>> y(inst.C, vector<GRBVar>(Bplus));
	vector<GRBVar> z(Bplus);
	for (int b = 0; b <= Bplus; b++) {
		for (int j = 0; j < inst.J; j++) {
			x[j][b] = model.addVar(0, inst.countedItems[j][2], 0, GRB_INTEGER);
		}
		if (b != Bplus) {
			for (int c = 0; c < inst.C; c++) {
				y[c][b] = model.addVar(0, 1, 0, GRB_BINARY);
			}
			z[b] = model.addVar(0, 1, 0, GRB_BINARY);
		}
	}

	model.update();

	// declare and calculate linear expressions
	vector<vector<GRBLinExpr>> fIn, fInS, fInR;		// the amount of (total/standard/reflected) flow entering each vertex
	vector<vector<GRBLinExpr>> fOut;				// the amount of flow leaving each vertex
	vector<vector<GRBLinExpr>> fSizeUsed(inst.C, vector<GRBLinExpr>(inst.W + 1)); // the number of arcs used of each item type 
	vector<GRBLinExpr> xAssigned(inst.J);			// the number of items of type j assigned to a bin / assigned to arcflow
	vector<GRBLinExpr> nbBinsPerColor(inst.C);		// the number of bins used for each color
	vector<GRBLinExpr> nbColorsPerBin(inst.B);		// the number of colors present in each bin
	vector<GRBLinExpr> fReflect(inst.C);			// the number of reflected arcs used per color
	vector<GRBLinExpr> binLoads(Bplus);				// the total load in each bin
	vector<vector<GRBLinExpr>> binLoadsPerColor(inst.C, vector<GRBLinExpr>(Bplus));	// the total load in each bin for each color
	GRBLinExpr nbBinsUsed = 0;						// the total number of bins used
	GRBLinExpr colorFragmentations = 0;				// the total number of color fragmentations
	int R, tail, head, c, w;

	if (reflect) {
		R = ceil(inst.W / double(2));
		fInS.resize(inst.C, vector<GRBLinExpr>(R + 1, 0));
		fInR.resize(inst.C, vector<GRBLinExpr>(R + 1, 0));
		fOut.resize(inst.C, vector<GRBLinExpr>(R + 1, 0));
	}
	else {
		fIn.resize(inst.C, vector<GRBLinExpr>(inst.W + 1, 0));
		fOut.resize(inst.C, vector<GRBLinExpr>(inst.W + 1, 0));
	}

	for (int c = 0; c < inst.C; c++) {					// loop over all colors
		for (int a = 0; a < G[c].A.size(); a++) {		// loop over all arcs
			tail = G[c].A[a][0], head = G[c].A[a][1], idx = G[c].A[a][2]; // decompose the arc

			if (reflect) {								// if reflect is used
				if (a == G[c].A.size() - 1) {			// final special arc
					fInR[c][R] += f[c][a];				// inflow of reflected arcs
					fOut[c][R] += f[c][a];				// outflow of all arcs
					fReflect[c] += f[c][a];				// the number of reflected arcs used
					nbBinsPerColor[c] += f[c][a];		// each pair of flows corresponds to one bin
				}
				else {									// for all other arcs
					fOut[c][tail] += f[c][a];     		// outflow
					if (idx > 0) {						// for item arcs
						fInS[c][head] += f[c][a];       // inflow using standard arcs (item + loss)
						w = inst.groupedItems[c][idx - 1][0];
						fSizeUsed[c][w] += f[c][a];		// number of items used of certain type
					}
					else if (idx == 0) {					// for loss arcs
						fInS[c][head] += f[c][a];        	// inflow using standard arcs (item + loss)
					}
					else {									// for reflected arcs
						fInR[c][head] += f[c][a];        	// inflow using reflected arcs
						w = inst.groupedItems[c][-idx - 1][0];
						fSizeUsed[c][w] += f[c][a];			// number of items used of certain type
						fReflect[c] += f[c][a];				// the number of reflected arcs used
						nbBinsPerColor[c] += f[c][a];		// each pair of flows corresponds to one bin
					}
				}
			}

			else {									// if reflect is not used
				fIn[c][head] += f[c][a];			// inflow
				fOut[c][tail] += f[c][a];			// outflow
				if (idx >= 0) {						// for item arcs
					w = inst.groupedItems[c][idx][0];
					fSizeUsed[c][w] += f[c][a];		// number of items used of certain type
				}
			}
		}

		if (not reflect) {
			nbBinsPerColor[c] = fOut[c][0];		// each flow corresponds to one monochromatic bin
		}
		nbBinsUsed += nbBinsPerColor[c];

		for (int b = 0; b < Bplus; b++) {		// loop over all traditional bins
			nbBinsPerColor[c] += y[c][b];		// colors in traditional bins add to color fragmentations as well
			nbColorsPerBin[b] += y[c][b];		// number of colors present in each bin
		}
	}

	for (int b = 0; b <= Bplus; b++) {			// loop over all traditional bins
		if (b != Bplus) { nbBinsUsed += z[b]; }	// the number of bins depends on the number of used traditional bins
		for (int j = 0; j < inst.J; j++) {		// loop over all item types
			xAssigned[j] += x[j][b];			// number of items of the type assigned to a bin / assigned to arcflow
			if (b != Bplus) {
				binLoads[b] += inst.countedItems[j][1] * x[j][b]; // the total load in each bin
				binLoadsPerColor[inst.countedItems[j][0]][b] += inst.countedItems[j][1] * x[j][b]; // the load per color
			}
		}
	}

	for (int c = 0; c < inst.C; c++) {				// loop over all colors
		colorFragmentations += nbBinsPerColor[c];	// add all color fragmentations per color
	}

	model.update();

	// set the objective: minimize the number of color fragmentations
	model.setObjective(colorFragmentations, GRB_MINIMIZE);

	// constraint 1: limited number of bins
	model.addConstr(nbBinsUsed <= inst.B);

	// constraints 2: assignment constraints
	for (int j = 0; j < inst.J; j++) {
		model.addConstr(xAssigned[j] == inst.countedItems[j][2], "assignment(" + to_string(j) + ")");
	}

	// constraints 3: flow conservation
	for (int c = 0; c < inst.C; c++) {			// loop over all colors
		if (not reflect) {						// if reflect is not used
			for (int v = 1; v < inst.W; v++) {  // loop over all vertices except 0 and W
				if (G[c].V[v]) {
					model.addConstr(fIn[c][v] == fOut[c][v]); // inflow = outflow
				}
			}
			// model.addConstr(fIn[c][inst.W] == fOut[c][0]); // bordery condition
		}
		else {									// if reflect is used
			for (int v = 1; v <= R; v++) {		// loop over all vertices
				if (G[c].V[v]) {				// only proceed if the vertex exists
					model.addConstr(fInS[c][v] == fInR[c][v] + fOut[c][v]); // standard inflow = reflected inflow + outflow
				}
			}
			model.addConstr(fOut[c][0] == 2 * fReflect[c]);	// each bin consists of 2 flows, one of which uses a reflected arc
		}
	}

	// constraints 4: item type quantities
	for (int j = 0; j < inst.J; j++) {								// loop over all item types
		c = inst.countedItems[j][0]; w = inst.countedItems[j][1];	// find color and size
		model.addConstr(fSizeUsed[c][w] == x[j][Bplus]);			// demand met
	}

	// constraints 5: bin loads
	for (int b = 0; b < Bplus; b++) {					// loop over all traditional bins
		model.addConstr(binLoads[b] <= inst.W * z[b]);	// bin capacity must be respected
	}

	// constraints 6: y-constraints
	for (int j = 0; j < inst.J; j++) {		// loop over all item types
		for (int b = 0; b < Bplus; b++) {	// loop over all traditional bins
			model.addConstr(x[j][b] <= min(inst.countedItems[j][2], inst.W / inst.countedItems[j][1]) * y[inst.countedItems[j][0]][b]);
		}
	}

	// constraints 7: z-constraints
	for (int c = 0; c < inst.C; c++) { // loop over all colors
		for (int b = 0; b < Bplus; b++) {	// loop over all traditional bins
			model.addConstr(y[c][b] <= z[b]);
		}
	}
	// (optional) multi-chromatic bins have at least two colors
	for (int b = 0; b < Bplus; b++) {
		model.addConstr(nbColorsPerBin[b] >= 2 * z[b]);
	}


	// (optional) constraints 8: bin load per color
	for (int c = 0; c < inst.C; c++) {		// loop over all colors
		for (int b = 0; b < Bplus; b++) {	// loop over all traditional bins
			model.addConstr(binLoadsPerColor[c][b] <= inst.W * y[c][b]);	// bin capacity must be respected for each color
		}
	}

	// (optional) constraints 9: lower bound on number of bins per color
	if (constrL2) {
		for (int c = 0; c < inst.C; c++) {						// loop over all colors
			model.addConstr(nbBinsPerColor[c] >= inst.LBs[c]);	// mimimum number of bins per color
		}
	}

	// (optional) constraints 10: symmetry breaking
	if (constrSym) {
		for (int b = 0; b < Bplus - 1; b++) {
			model.addConstr(z[b] >= z[b + 1]);
		}
	}

	model.update();

	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	// if a warm start is used, start with the provided feasible solution
	if (warmStart && inst.feasSol.size() > 0) {

		// find the colors per bin
		vector<vector<int>> colorsPerBin(inst.feasSol.size());
		for (int b = 0; b < inst.feasSol.size(); b++) {	// loop over all bins
			for (int c = 0; c < inst.C; c++) {			// loop over all colors
				if (inst.feasSol[b][c].size() > 0) {	// if the color is contained in the bin
					colorsPerBin[b].push_back(c);
				}
			}
		}

		// set the initial values for the multichromatic bins
		vector<vector<int>> xInit(inst.J, vector<int>(Bplus + 1, 0));
		vector<vector<int>> yInit(inst.C, vector<int>(Bplus, 0));
		vector<int> zInit(Bplus, 0);
		int q;

		// find a mapping from an items color and size to its index
		vector<vector<int>> cw2j(inst.C, vector<int>(inst.W+1, -1));
		for (int j = 0; j < inst.J; j++) {
			c = inst.countedItems[j][0]; w = inst.countedItems[j][1];
			cw2j[c][w] = j;
		}

		int multiBinIdx = -1;
		for (int b = 0; b < inst.feasSol.size(); b++) { // loop over all multichromatic bins
			if (colorsPerBin[b].size() > 1) {
				multiBinIdx++;							// open a new multichromatic bin
				zInit[multiBinIdx] = 1;
				for (const int& c : colorsPerBin[b]) { // loop over all the colors in the bin
					yInit[c][multiBinIdx] = 1;
					for (int j = 0; j < inst.feasSol[b][c].size(); j++) {
						w = inst.feasSol[b][c][j][0]; q = inst.feasSol[b][c][j][1];
						xInit[cw2j[c][w]][multiBinIdx] = q;
					}
				}
			}
		}

		// find out how many units of each type are packed through arcflow
		for (int j = 0; j < inst.J; j++) {
			xInit[j][Bplus] = inst.countedItems[j][2];
			for (int b = 0; b < Bplus; b++) {
				xInit[j][Bplus] -= xInit[j][b];
			}
		}

		// set the initial values
		for (int b = 0; b <= Bplus; b++) {
			if (b != Bplus) {
				z[b].set(GRB_DoubleAttr_Start, zInit[b]);
				for (int c = 0; c < inst.C; c++) {
					y[c][b].set(GRB_DoubleAttr_Start, yInit[c][b]);
				}
			}
			for (int j = 0; j < inst.J; j++) {
				x[j][b].set(GRB_DoubleAttr_Start, xInit[j][b]);
			}
		}

		// set the values for the monochromatic bins
		vector<vector<int>> fInit(inst.C);

		int vCur, vNew, a;
		vector<int> I1, I2;
		int loadI1, loadI2;
		int collidePoint;
		bool reflectUsed, reflectForI1;

		for (int c = 0; c < inst.C; c++) { // loop over all colors

			// initialize the array holding the initial values of the variables
			fInit[c].resize(G[c].A.size(), 0);

			// Find a mapping from the arcs (i.e., the two endpoints) to their indices
			vector<vector<int>> AS(R + 1, vector<int>(R + 1, -1));	// standard arcs
			vector<vector<int>> AR(R + 1, vector<int>(R + 1, -1));	// reflected arcs
			int zeroLossI;											// loss arc from 0 to 0
			vector<int> AL(R + 1);									// other loss arcs 
			for (int a = 0; a < G[c].A.size(); a++) {				// loop over all arcs
				tail = G[c].A[a][0], head = G[c].A[a][1], idx = G[c].A[a][2]; // decompose the arc
				if (idx > 0) { AS[tail][head] = a; }				// standard arc
				else if (idx < 0) { AR[tail][head] = a; }			// reflected arc
				else if (head == 0) { zeroLossI = a; }				// loss arc from 0 to 0
				else { AL[tail] = a; }								// other loss arcs
			}

			for (int b = 0; b < inst.feasSol.size(); b++) { // loop over all monochromatic bins of color c
				if (colorsPerBin[b].size() == 1 && colorsPerBin[b][0] == c) {
					// divide the items over two half bins
					I1.clear(); I2.clear(); loadI1 = 0; loadI2 = 0;
					for (int j = 0; j < inst.feasSol[b][c].size(); j++) {			// for every item type
						w = inst.feasSol[b][c][j][0]; q = inst.feasSol[b][c][j][1]; // decompose the item
						for (int k = 0; k < q; k++) {								// for all copies of the same length
							if (loadI1 >= R) {										// If I1 is full, add to I2
								I2.push_back(w); loadI2 += w;
							}
							else if (loadI2 >= R) {									// if I2 is full, add to I1
								I1.push_back(w); loadI1 += w;
							}
							else if (2 * loadI1 + w <= inst.W) {					// if at least half of the item fits in I1, add to I1
								I1.push_back(w); loadI1 += w;
							}
							else {													// otherwise, add to I2
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
							a = AS[vCur][vNew];
						}
						else {
							vNew = inst.W - loadI1;
							a = AR[vCur][vNew];
							reflectUsed = true;
						}
						fInit[c][a]++;	// add 1 to the flow on that arc
						vCur = vNew;	// move to the new point
					}

					if (not reflectForI1 && collidePoint == 0) {
						a = zeroLossI;
						fInit[c][a]++;
					}

					while (vCur != collidePoint) {	// take loss arcs until reaching the collide point
						a = AL[vCur];				// take a loss arc
						fInit[c][a]++;				// add 1 to the flow on that arc
						vCur = G[c].A[a][1];		// move to the head of the loss arc
					}
					if (reflectForI1 && not reflectUsed) {	// if reflect must still be used for I1
						a = AL[R];							// take the final reflected loss arc
						fInit[c][a]++;						// add 1 to the flow on that arc
					}

					// Halfbin I2
					vCur = 0;
					for (int i = 0; i < I2.size(); i++) {
						if (vCur + I2[i] <= R) {
							vNew = vCur + I2[i];
							a = AS[vCur][vNew];
						}
						else {
							vNew = inst.W - loadI2;
							a = AR[vCur][vNew];
							reflectUsed = true;
						}
						fInit[c][a]++;	// add 1 to the flow on that arc
						vCur = vNew;	// move to the new point
					}

					if (reflectForI1 && collidePoint == 0) {
						a = zeroLossI;
						fInit[c][a]++;
					}

					while (vCur != collidePoint) {	// take loss arcs until reaching the collide point
						a = AL[vCur];				// take a loss arc
						fInit[c][a]++;				// add 1 to the flow on that arc
						vCur = G[c].A[a][1];		// move to the head of the loss arc
					}
					if (not reflectForI1 && not reflectUsed) {  // if reflect must still be used for I2
						a = AL[R];								// take the final reflected loss arc
						fInit[c][a]++;							// add 1 to the flow on that arc
					}
				}
			}
		}

		// set the initial values
		for (int c = 0; c < inst.C; c++) {
			for (int a = 0; a < G[c].A.size(); a++) {
				f[c][a].set(GRB_DoubleAttr_Start, fInit[c][a]);
			}
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
		// it could be that the current instance is infeasible because the number of allowed multi-chromatic bins is too low
		if (Bplus < solUB.UB - inst.B) {
			sol.feas = -1;
		}
		else {
			sol.feas = 0;
		}
	}

	// if a solution has been found within the time limit
	else if (model.get(GRB_IntAttr_SolCount) >= 1) {
		sol.feas = 1;
		sol.UB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		// if the solution is optimal for the current value of Bplus
		if (sol.LB == sol.UB) {
			// optimality is only certain is the UB equals the MHLB bound or if Bplus is high enough
			if (sol.UB == inst.LBtot || Bplus >= solUB.UB - inst.B || Bplus == inst.B) {
				sol.opt = 1;
			}
			else {
				sol.opt = 0;
			}
		}

		// first decompose the mono-chromatic bins color by color
		vector<vector<vector<int>>> AByTail;
		int fval, nFlows;
		vector<int> a, sizesInBin;
		vector<vector<vector<int>>> currentBin;
		vector<int> halfBin, fullBin;							// The items packed in the current half bin / full bin 
		vector<vector<vector<int>>> standardHalfBinsPerEnd;		// The half bins not containing a reflected arc per end vertex 
		vector<vector<vector<int>>> reflectedHalfBinsPerEnd;	// The half bins containing a reflected arc per end vertex
		if (reflect) { standardHalfBinsPerEnd.resize(R + 1); reflectedHalfBinsPerEnd.resize(R + 1); }
		bool reflectUsed; // indicates if a reflected arc has been used on the path or not

		for (int c = 0; c < inst.C; c++) { // loop over all colors
			// reset the arc list
			AByTail.clear();
			if (reflect) AByTail.resize(R + 1);
			else AByTail.resize(inst.W + 1);

			// construct the new list
			for (int a = 0; a < G[c].A.size(); a++) {					// loop over all arcs
				fval = ceil(f[c][a].get(GRB_DoubleAttr_X) - EPSILON);	// find variable value
				for (int k = 0; k < fval; k++) {						// loop over all copies of the same arc
					AByTail[G[c].A[a][0]].push_back(G[c].A[a]);			// add arc to set of arcs with same tail
				}
			}

			// create a path one-by-one
			if (not reflect) {						// if reflect is not used
				nFlows = AByTail[0].size();
				for (int b = 0; b < nFlows; b++) {
					sizesInBin.clear();				// reset the bin
					currentBin.clear(); currentBin.resize(inst.C);
					tail = 0;						// start at 0
					while (tail != inst.W) {		// end after reaching W
						a = AByTail[tail].back();	// select an arc
						head = a[1];				// find the head of the arc
						idx = a[2];					// find the index of the arc
						AByTail[tail].pop_back();	// remove the used arc
						tail = head;				// move to the new point
						if (idx >= 0) {				// for item arcs: 
							sizesInBin.push_back(inst.groupedItems[c][idx][0]); // add item to the bin
						}
					}
					currentBin[c] = sortAndCount(sizesInBin);
					sol.binPacking.push_back(currentBin);
				}
			}

			else { // if reflect is used
				nFlows = AByTail[0].size();
				for (int b = 0; b < nFlows; b++) {  // loop over the to-be-filled bins
					halfBin.clear();				// reset the half bin
					tail = 0;						// start at 0
					reflectUsed = false;			// no reflect arc has been used yet
					while (AByTail[tail].size() > 0 && !reflectUsed) { // while still a leaving arc and no relected arc used
						a = AByTail[tail].back();	// select an arc
						AByTail[tail].pop_back();	// remove the used arc
						head = a[1];				// find the head of the arc
						idx = a[2];					// find the index of the arc
						if (idx > 0) {				// for item arcs
							halfBin.push_back(inst.groupedItems[c][idx - 1][0]);
						}
						else if (idx < 0) {			// for reflected arcs
							halfBin.push_back(inst.groupedItems[c][-idx - 1][0]);
							reflectUsed = true;
						}
						else if (tail == R) {		// for the loss/reflected arc from R to R
							reflectUsed = true;
						}
						tail = head;				// move to the new point
					}
					if (reflectUsed) reflectedHalfBinsPerEnd[head].push_back(halfBin);
					else standardHalfBinsPerEnd[head].push_back(halfBin);
				}

				// combine the half bins
				for (int end = 0; end <= R; end++) {
					while (reflectedHalfBinsPerEnd[end].size() > 0) {
						currentBin.clear();  currentBin.resize(inst.C);
						fullBin = reflectedHalfBinsPerEnd[end].back();
						halfBin = standardHalfBinsPerEnd[end].back();
						fullBin.insert(fullBin.end(), halfBin.begin(), halfBin.end());
						reflectedHalfBinsPerEnd[end].pop_back();
						standardHalfBinsPerEnd[end].pop_back();
						currentBin[c] = sortAndCount(fullBin);
						sol.binPacking.push_back(currentBin);
					}
				}
			}
		}

		// second decompose the multi-chromatic bins
		int xval, c;
		for (int b = 0; b < Bplus; b++) {									// for every multi-chromatic bin
			if (ceil(z[b].get(GRB_DoubleAttr_X) - EPSILON) > 0) {			// if the bin is used
				currentBin.clear(); currentBin.resize(inst.C);				// reset the bin
				for (int j = 0; j < inst.J; j++) {							// loop over all items
					xval = ceil(x[j][b].get(GRB_DoubleAttr_X) - EPSILON);	// find variable value
					if (xval > 0) {											// if the variable is used
						c = inst.countedItems[j][0];						// find the item color
						currentBin[c].push_back({ inst.countedItems[j][1], xval }); // add the item to the bin
					}
				}
				sol.binPacking.push_back(currentBin);						// add the bin to the set of bins
			}
		}
	}
	sol.timeT += getCPUTime() - start;	 // save the total time

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