#include "AFCSP.h"

void Graph::print() {
	// This function is used to print the arcflow graph

	cout << "Graph G = (V, A):\n";
	//cout << "number of vertices: " << sol.G.Nvert << "\n";
	//cout << "number of arcs: " << sol.G.Narcs << "\n";
	printVector(V, "- V");
	print2DVector(A, "- A", "arc");
}

Graph graphConstructionCSP(const vector<vector<int>>& wq, int W, bool terminalZero, bool useLoss) {
	// This function is used to construct the AF graph for a CSP instance
	//  with items wq (each row gives the size of an item type and its corresponding quantity)
	//  and bin capacity W

	Graph G; G.V.resize(W + 1, false); G.V[0] = true; // initialize graph
	int w, q;										  // current item size (w) and quantity (q)
	// Construct vertices and item arcs
	for (int j = 0; j < wq.size(); j++) {	          // loop over item types (ordered from large to small size)
		w = wq[j][0], q = wq[j][1];		              // find size (w) and quantity (q) of current item type
		for (int v = W - w; v >= 0; v--) {            // loop over vertices in backwards order
			if (G.V[v]) {							  // only process vertices that were already there
				for (int k = 0; k < q; k++) {		  // loop over copies of the same item type
					if (v + (k + 1) * w <= W) {
						G.A.push_back({ v + k * w,v + (k + 1) * w, j }); // add the new arc
						if (not(G.V[v + (k + 1) * w])) {  // continue as long as you reach vertices
							G.V[v + (k + 1) * w] = true;  // add the new vertex
						}
						else break;						  // stop if the new vertex was already there
					}
					else break;
				}
			}
		}
	}

	// Note that node W is not created here if it wasn't there yet: needed for AF2

	if (terminalZero) {
		for (int v = 1; v <= W; v++) {						// loop over vertices
			if (G.V[v]) { G.A.push_back({ v, 0, -1 }); }	// add an arc to 0 for each existing vertex
		}
	}
	else if (useLoss) {
		for (int v = 1; v < W; v++) {						// loop over vertices
			if (G.V[v]) { G.A.push_back({ v, W, -1 }); }	// add an arc to W for each existing vertex
		}
	}

	// Return the graph
	return G;
}

Solution solveAFCSP(const Instance& inst, const vector<vector<int>>& wq, int L2, int c, int early) {
	// This function solves a CSP instance using AF, given:
	// - items wq (each row contains the size w and corresponding quanity q for an item type) 
	// - bin capacity W
	// - # bins B
	// - it is assumed that all items have the same color c
	//		if color doesn't matter, it can be set to 0

	int timeLimit = 1800;    	   // time limit in seconds
	double start = getCPUTime();   // starting time
	Solution sol;				   // initialize Solution struct
	sol.method = "AFCSP";		   // method declaration

	// Create graph
	Graph G = graphConstructionCSP(wq, inst.W);
	// G.print(); // print the graph

	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();                       	// remove Gurobi message
	env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);         // create a new model

	// declaration of the variables for the model
	vector<GRBVar> f(G.A.size());
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
	// constraints 3 (optional): L2-based constraints
	if (L2 > 0) model.addConstr(fOut[0] >= L2);

	model.update();
	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	if (early > 0) {
		model.getEnv().set(GRB_DoubleParam_BestObjStop, early);
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
		if (sol.LB == sol.UB || sol.UB <= early) {
			sol.opt = 1;

			// store the items assigned to each bin (sorted by color)
			vector<int> binItemSizes;
			vector<vector<vector<int>>> currentBin;

			// first store the used arcs by tail
			int fval;
			vector<vector<vector<int>>> AByTail(inst.W + 1);
			for (int a = 0; a < G.A.size(); a++) {
				fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
				for (int k = 0; k < fval; k++) {
					AByTail[G.A[a][0]].push_back(G.A[a]); // add arc to set of arcs with same tail
				}
			}
			// then fill the bins one by one by finding a 0-W path every iteration
			int tail, head, idx;
			vector<int> a;
			for (int b = 0; b < sol.UB; b++) { // loop over the to-be-filled bins
				tail = 0;
				while (tail != inst.W) {
					a = AByTail[tail].back();	// select an arc
					head = a[1];				// find the head of the arc
					idx = a[2];					// find the index of the arc
					if (idx != -1) {			// add the item to the bin
						binItemSizes.push_back(wq[idx][0]);
					}
					AByTail[tail].pop_back();	// remove the used arc
					tail = head;				// move to the new point
				}
				currentBin.resize(inst.C);					// the bin contains a seperate vector per color
				currentBin[c] = sortAndCount(binItemSizes); // group and count the items in the bin
				sol.binPacking.push_back(currentBin);		// add the last bin to the set of all bins
				binItemSizes.clear();						// start a new empty bin
				currentBin.clear();							// start a new empty bin
			}
		}
	}
	sol.timeT = getCPUTime() - start; // save the total time
	return sol;                                 // return the solution
}

Solution solveAFCSPConcentration(const Instance& inst, const vector<vector<int>>& wq, int L2, int c, bool startWithReflect, int searchStrategy, bool lexicographic) {
	// This function solves an instance of the CSP with loss concentration
	
	bool printing = false;
	int timeLimit = 1800;    	   // time limit in seconds
	double start = getCPUTime();   // starting time
	Solution sol;				   // initialize Solution struct
	sol.method = "AFCSP";		   // method declaration

	// Create graph
	Graph G = graphConstructionCSP(wq, inst.W, true);
	// G.print(); // print the graph

	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();                       	// remove Gurobi message
	env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);         // create a new model

	// declaration of the variables for the model
	vector<GRBVar> f(G.A.size());
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
	vector<GRBLinExpr> typeUsed(wq.size(), 0); 	// the number of arcs used of each item type 
	vector<GRBLinExpr> lossAtLeast(inst.W, 0); // the number of loss arcs leaving at least a certain amount of loss space

	// calculate the linear expressions	
	for (int a = 0; a < G.A.size(); a++) {   				// loop over all arcs
		fIn[G.A[a][1]] += f[a];        						// inflow
		fOut[G.A[a][0]] += f[a];       						// outflow
		if (G.A[a][2] >= 0) {				// for item arcs
			typeUsed[G.A[a][2]] += f[a];	// number of items used of certain type
		}
		else {
			for (int l = 1; l <= inst.W - G.A[a][0]; l++) {
				lossAtLeast[l] += f[a];
			}
		}
	}
	model.update();

	// set the objective: minimize the number of flow going out of 0
	model.setObjective(fOut[0], GRB_MINIMIZE);
	// constraints 1: flow conservation
	for (int v = 0; v <= inst.W; v++) {       	// loop over all vertices
		if (G.V[v])
			model.addConstr(fIn[v] == fOut[v]); // inflow = outflow
	}
	// constraints 2: item type quantities
	for (int j = 0; j < wq.size(); j++) {			// loop over all item types
		model.addConstr(typeUsed[j] == wq[j][1]);	// demand met
	}
	// constraints 3 (optional): L2-based constraints
	if (L2 > 0) model.addConstr(fOut[0] >= L2);

	model.update();
	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));


	if (startWithReflect) {
		// start with the Reflect solution
		Solution solReflect = solveReflectCSP(inst, wq, L2, c);
		// solReflect.print();

		// find a mapping from the arcs to their indices
		vector<vector<int>> AN2I(inst.W + 1, vector<int>(inst.W + 1, -1));
		for (int a = 0; a < G.A.size(); a++) {   // loop over all arcs
			AN2I[G.A[a][0]][G.A[a][1]] = a;
		}

		// set the values of the initial solution
		vector<int> fInit(G.A.size(), 0);
		int cur_v, new_v, a;
		for (int b = 0; b < solReflect.binPacking.size(); b++) {// for every bin
			cur_v = 0;
			// item arcs
			for (int j = 0; j < solReflect.binPacking[b][c].size(); j++) {
				for (int k = 0; k < solReflect.binPacking[b][c][j][1]; k++) {
					new_v = cur_v + solReflect.binPacking[b][c][j][0];
					a = AN2I[cur_v][new_v];
					fInit[a]++;
					cur_v = new_v;
				}
			}
			// loss arc
			a = AN2I[cur_v][0];
			fInit[a]++;
		}
		for (int a = 0; a < G.A.size(); a++) {
			f[a].set(GRB_DoubleAttr_Start, fInit[a]);
		}
		model.getEnv().set(GRB_DoubleParam_BestObjStop, solReflect.UB);
		model.set(GRB_IntParam_Presolve, 0);
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

		if (startWithReflect) { 
			sol.LB = sol.UB; 
			model.set(GRB_IntParam_Presolve, -1);
		}

		// if the solution is optimal
		if (sol.LB == sol.UB) {
			sol.opt = 1;
			int minBins = sol.UB;
			GRBLinExpr emptyObjective = 0;
			model.setObjective(emptyObjective, GRB_MINIMIZE);
			model.addConstr(fOut[0] == minBins);
			GRBConstr tempConstr = model.addConstr(emptyObjective == 0); // add a mock constraint
			int fval, spaceGoal, mostSpace, totSpace;
			vector<vector<vector<int>>> AByTail(inst.W + 1);
			bool first = true;
			int LBspace, UBspace;

			// find the minimum weight
			int wMin = wq[wq.size() - 1][0];
			// find the sum of the weights
			int wSum = 0;
			for (int j = 0; j < wq.size(); j++) {
				wSum += wq[j][1] * wq[j][0];
			}

			if (not lexicographic) {
				// partially decompose the solution and find the maximum remaining space
				mostSpace = 0;
				for (int a = 0; a < G.A.size(); a++) {
					fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
					for (int k = 0; k < fval; k++) {
						AByTail[G.A[a][0]].push_back(G.A[a]); // add arc to set of arcs with same tail
					}
					if (fval > 0 && G.A[a][2] == -1) {
						mostSpace = max(mostSpace, inst.W - G.A[a][0]);
					}
				}
				totSpace = minBins * inst.W - wSum;

				if (printing) {
					// cout << "-----------------------------------------------------------\n";
					cout << "\n";
					cout << "Color " << c << ":\n";
					cout << "Minimum number of bins: " << sol.UB << ", total remaining space: " << totSpace << ", W-w_min = " << inst.W - wMin << "\n";
					cout << "Solution found with maximum remaining space: " << mostSpace << "\n";
				}

				// initialize the LB and UB for the search
				LBspace = mostSpace;
				UBspace = min(totSpace, inst.W - wq[wq.size() - 1][0]);

				// continue the search while the LB and UB do not match
				while (LBspace != UBspace) {

					// depending on the search strategy, set the next value to try
					if (searchStrategy == 0) { // binary search
						spaceGoal = max(LBspace + 1, LBspace + (UBspace - LBspace) / 2);
					}
					else if (searchStrategy == -1) { // descending search
						spaceGoal = UBspace;
					}
					else if (searchStrategy == 1) { // ascending search
						spaceGoal = LBspace + 1;
					}

					// remove the old constraint
					if (first) { first = false; }
					else { model.remove(tempConstr); }

					// add the new constraint
					tempConstr = model.addConstr(lossAtLeast[spaceGoal] >= 1);

					// resolve the model
					model.optimize();

					// if a solution is found, partially decompose it and update the LB
					if (model.get(GRB_IntAttr_SolCount) >= 1) {
						// partially decompose the solution and find the maximum remaining space
						mostSpace = 0;
						for (int a = 0; a < G.A.size(); a++) {
							fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON);
							for (int k = 0; k < fval; k++) {
								AByTail[G.A[a][0]].push_back(G.A[a]);
							}
							if (fval > 0 && G.A[a][2] == -1) {
								mostSpace = max(mostSpace, inst.W - G.A[a][0]);
							}
						}
						if (printing) { cout << "Solution found with maximum remaining space: " << mostSpace << "\n"; }
						LBspace = mostSpace;
					}
					// if no solution is found, update the UB
					else {
						if (printing) { cout << "No solution found with maximum remaining space: " << spaceGoal << "\n"; }
						for (int v = inst.W - spaceGoal + 1; v < inst.W; v++) {
							if (G.V[v]) {
								UBspace = inst.W - v;
								break;
							}
						}
					}
				}
			}

			else { // if solved lexicographically

				vector<int> spacePerBin;
				totSpace = minBins * inst.W - wSum;
				int spaceFixed = 0;

				// Lexicographically maximize the loss per bin
				for (int b = 0; b < minBins; b++) {	// loop over the bins

					if (b==0) {
						// find the space per bin
						for (int a = 0; a < G.A.size(); a++) {
							fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
							if (fval > 0 && G.A[a][2] == -1) {
								for (int k = 0; k < fval; k++) {
									spacePerBin.push_back(inst.W - G.A[a][0]);
								}
							}
						}
						sort(spacePerBin.begin(), spacePerBin.end(), greater<int>());
					}

					if (printing) {
						// cout << "-----------------------------------------------------------\n";
						cout << "\n";
						cout << "Color " << c << ", bin " << b << ":\n";
						cout << "Minimum number of bins: " << sol.UB << ", total remaining unfixed space: " << totSpace - spaceFixed << ", W-w_min = " << inst.W - wMin << "\n";
						printVector(spacePerBin, "spacePerBin");
						cout << "Solution found with remaining space in bin " << b << " of: " << spacePerBin[b] << "\n";
					}

					// initialize the LB and UB for the search
					LBspace = spacePerBin[b];
					UBspace = min(totSpace - spaceFixed, inst.W - wq[wq.size() - 1][0]);
					if (b > 0) { UBspace = min(UBspace, spacePerBin[b - 1]); }

					// continue the search while the LB and UB do not match
					while (LBspace != UBspace) {

						// depending on the search strategy, set the next value to try
						if (searchStrategy == 0) { // binary search
							spaceGoal = max(LBspace + 1, LBspace + (UBspace - LBspace) / 2);
						}
						else if (searchStrategy == -1) { // descending search
							spaceGoal = UBspace;
						}
						else if (searchStrategy == 1) { // ascending search
							spaceGoal = LBspace + 1;
						}

						// remove the old constraint
						model.remove(tempConstr);

						// add the new constraint
						tempConstr = model.addConstr(lossAtLeast[spaceGoal] >= b+1);

						// resolve the model
						model.optimize();

						// if a solution is found, partially decompose it and update the LB
						if (model.get(GRB_IntAttr_SolCount) >= 1) {
							// partially decompose the solution and find the maximum remaining space
							spacePerBin.clear();
							for (int a = 0; a < G.A.size(); a++) {
								fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
								if (fval > 0 && G.A[a][2] == -1) {
									for (int k = 0; k < fval; k++) {
										spacePerBin.push_back(inst.W - G.A[a][0]);
									}
								}
							}
							sort(spacePerBin.begin(), spacePerBin.end(), greater<int>());
							if (printing) { cout << "Solution found with remaining space in bin " << b << " of: " << spacePerBin[b] << "\n"; }
							LBspace = spacePerBin[b];
						}
						// if no solution is found, update the UB
						else {
							if (printing) { cout << "No solution found with remaining space in bin " << b << " of: " << spaceGoal << "\n"; }
							for (int v = inst.W - spaceGoal + 1; v < inst.W; v++) {
								if (G.V[v]) {
									UBspace = inst.W - v;
									break;
								}
							}
						}
					}
					if (LBspace <= 0) { break; }
					if (printing) { cout << "Maximum remaining space in bin " << b << ": " << LBspace << "\n"; }
					spaceFixed += LBspace;
					// fix the loss in the considered bin from now on
					model.addConstr(lossAtLeast[LBspace] >= b + 1);
				}
			}

			// Fully decompose the best solution that was found
			// store the items assigned to each bin (sorted by color), by finding a 0-0 cycle every iteration

			//// WHY WAS THIS THERE??? Possibly for the lexicographic stuff??
			//model.remove(tempConstr);
			//model.optimize();
			//for (int a = 0; a < G.A.size(); a++) {
			//	fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
			//	for (int k = 0; k < fval; k++) {
			//		AByTail[G.A[a][0]].push_back(G.A[a]); // add arc to set of arcs with same tail
			//	}
			//}
			vector<int> binItemSizes;
			vector<vector<vector<int>>> currentBin;
			int tail, head;
			vector<int> a;
			for (int b = 0; b < minBins; b++) { // loop over the to-be-filled bins
				tail = 0;
				while (true) {
					a = AByTail[tail].back();	// select an arc
					head = a[1];				// find the head of the arc
					idx = a[2];					// find the index of the arc
					AByTail[tail].pop_back();	// remove the used arc
					if (idx != -1) {			
						binItemSizes.push_back(wq[idx][0]); // add the item to the bin
						tail = head;					    // move to the new point
					}
					else {break;}
				}
				currentBin.resize(inst.C);					// the bin contains a seperate vector per color
				currentBin[c] = sortAndCount(binItemSizes); // group and count the items in the bin
				sol.binPacking.push_back(currentBin);		// add the last bin to the set of all bins
				binItemSizes.clear();						// start a new empty bin
				currentBin.clear();							// start a new empty bin
			}
		}
	}
	sol.timeT = getCPUTime() - start; // save the total time
	return sol;                                 // return the solution
}

Graph graphConstructionCSPReflect(const vector<vector<int>>& wq, int W) {
	// This function is used to construct the REFLECT AF graph for a CSP instance
	//  with items wq (each row gives the size of an item type and its corresponding quantity)
	//  and bin capacity W
	// 3rd index: 1,...,j = item arc, -1,...,-j = reflected arc, 0 = loss arc

	int R = ceil(W / double(2));                      // find halfway point
	Graph G; G.V.resize(R + 1, false); G.V[0] = true; // initialize graph
	int w, q;										  // current item size (w) and quantity (q)
	vector<bool> VRef(R + 1, false);                  // indicates which vertices can be reached with a reflected arc

	// construct vertices and item arcs
	for (int j = 0; j < wq.size(); j++) {	          // loop over item types (ordered from large to small size)
		w = wq[j][0], q = wq[j][1];		              // find size (w) and quantity (q) of current item type
		for (int v = R - 1; v >= 0; v--) {            // loop over vertices in backwards order
			if (G.V[v]) {							  // only process vertices that were already there
				for (int k = 0; k < q; k++) {		  // loop over copies of the same item type
					if (v + (k + 1) * w <= W / 2) {     // if the head still fits within the first half of the bin
						G.A.push_back({ v + k * w,v + (k + 1) * w, j + 1 }); // add the new standard arc
						if (not(G.V[v + (k + 1) * w])) {  // continue as long as you reach new vertices
							G.V[v + (k + 1) * w] = true;  // add the new vertex
						}
						else break;						  // stop if the new vertex was already there
					}
					else if (2 * (v + k * w) <= W - w) {    // if otherwise the reflected head is at least as large as the tail
						G.A.push_back({ v + k * w, W - (v + (k + 1) * w), -(j + 1) }); // add the new reflected arc
						// G.V[W - (v + (k + 1) * w)] = true;  // DON'T add the new vertex yet (since it is not a possible tail)
						VRef[W - (v + (k + 1) * w)] = true; // the vertex can be reached using a reflected arc
						break;                              // move to the next possible tail
					}
					else break;         // if the head doesn't fit at all, move to the next possible tail
				}
			}
		}
	}
	// activate the terminal vertex
	G.V[R] = true, VRef[R] = true;

	// add a 0 --> 0 loss arc in case of an item with w = W
	if (wq[0][0] == W) G.A.push_back({ 0,0,0 });

	// add the loss arcs
	// standard way: from every vertex to the next vertex
	int vCur = 0;
	for (int v = 1; v <= R; v++) {
		if (VRef[v]) G.V[v] = true;  // first activate vertices that are only reachable by reflected arcs
		if (G.V[v]) {
			G.A.push_back({ vCur, v, 0 });
			vCur = v;
		}
	}
	// alternative way to add loss arcs: from every vertex to the next vertex that is the head of at least 1 reflected arc
	// THIS METHOD DOES NOT YET SOLVE THE BUG THAT I HAD
	//vector<int> vTails = {0};
	//for (int v = 1; v <= R; v++) {
	//    if (VRef[v]) {
	//        for (const int& vTail : vTails) G.A.push_back({ vTail, v, 0 });
	//        vTails.clear();
	//    }
	//    if (G.V[v]) vTails.push_back(v);
	//}

	// add the final refect/loss arc (intentially given the last index in A)
	G.A.push_back({ R, R, 0 });

	// return the graph
	return G;
}

Solution solveReflectCSP(const Instance& inst, const vector<vector<int>>& wq, int L2, int c, int early) {
	// This function solves a CSP instance using REFLECT AF, given:
	// - items wq (each row contains the size w and corresponding quanity q for an item type) 
	// - bin capacity W
	// - # bins B
	// - it is assumed that all items have the same color c
	//		if color doesn't matter, it can be set to 0

	int timeLimit = 1800;    	   // time limit in seconds
	double start = getCPUTime();   // starting time
	Solution sol;				   // initialize Solution struct
	sol.method = "REFLECTCSP";	   // method declaration

	// Create graph
	Graph G = graphConstructionCSPReflect(wq, inst.W);
	// G.print(); // print the graph

	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();                       	// remove Gurobi message
	env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);         // create a new model

	// declaration of the variables for the model
	vector<GRBVar> f(G.A.size());
	// for a in A: fa = # arc a is used
	int tail, head, idx; // 3rd index : 1, ..., j = item arc, -1, ..., -j = reflected arc, 0 = loss arc
	for (int a = 0; a < G.A.size() - 1; a++) {
		idx = G.A[a][2];
		if (idx >= 1) {
			f[a] = model.addVar(0, wq[idx - 1][1], 0, GRB_INTEGER);
		}
		else if (idx <= -1) {
			f[a] = model.addVar(0, wq[-idx - 1][1], 0, GRB_INTEGER);
		}
		else {
			f[a] = model.addVar(0, 2 * inst.B, 0, GRB_INTEGER);
		}
	}
	f[G.A.size() - 1] = model.addVar(0, inst.B, 0, GRB_INTEGER);

	model.update();

	// declare linear expressions
	int R = ceil(inst.W / double(2));			// find halfway point
	vector<GRBLinExpr> fInS(R + 1, 0);    		// the amount of flow entering each vertex using standard arcs (item + loss)
	vector<GRBLinExpr> fInR(R + 1, 0);    		// the amount of flow entering each vertex using reflected arcs
	vector<GRBLinExpr> fOut(R + 1, 0);   		// the amount of flow leaving each vertex
	vector<GRBLinExpr> typeUsed(wq.size(), 0);	// the number of arcs used of each item type 
	GRBLinExpr fReflect(0);						// the number of reflected arcs used

	// calculate the linear expressions	
	for (int a = 0; a < G.A.size() - 1; a++) {   				 // loop over all arcs (treat final arc as special case)
		tail = G.A[a][0], head = G.A[a][1], idx = G.A[a][2]; // decompose the arc
		if (idx > 0) {					// for item arcs
			fInS[head] += f[a];        	// inflow using standard arcs (item + loss)
			typeUsed[idx - 1] += f[a];	// number of items used of certain type
		}
		else if (idx == 0) {			// for loss arcs
			fInS[head] += f[a];        	// inflow using standard arcs (item + loss)
		}
		else {							// for reflected arcs
			fInR[head] += f[a];        	// inflow using reflected arcs
			typeUsed[-idx - 1] += f[a];	// number of items used of certain type
			fReflect += f[a];			// number of reflected arcs
		}
		fOut[tail] += f[a];       		// outflow
	}
	// for the final loss arc:
	int aSpecial = G.A.size() - 1;
	fInR[R] += f[aSpecial], fOut[R] += f[aSpecial], fReflect += f[aSpecial];

	model.update();

	// set the objective: minimize the number of reflect arcs
	model.setObjective(fReflect, GRB_MINIMIZE);
	// constraints 1: flow conservation
	for (int v = 1; v <= R; v++) {       	// loop over all vertices
		if (G.V[v])
			model.addConstr(fInS[v] == fInR[v] + fOut[v]); // standard inflow = reflected inflow + outflow
	}
	// constraints 2: boundary conditions
	model.addConstr(fOut[0] == 2 * fReflect);	// each bin consists of 2 flows, one of which uses a relect arc
	// constraints 3: item type quantities
	for (int j = 0; j < wq.size(); j++) {			// loop over all item types
		model.addConstr(typeUsed[j] == wq[j][1]);	// demand met
	}
	// constraints 4 (optional): L2-based constraints
	if (L2 > 0) model.addConstr(fReflect >= L2);

	model.update();
	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	if (early > 0) {
		model.getEnv().set(GRB_DoubleParam_BestObjStop, early);
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
		if (sol.LB == sol.UB || sol.UB <= early) {
			sol.opt = 1;
			// store the items assigned to each bin (sorted by color)

			// first store the used arcs by tail
			int fval;
			vector<vector<vector<int>>> AByTail(R + 1);
			for (int a = 0; a < G.A.size(); a++) {
				fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
				for (int k = 0; k < fval; k++) {
					AByTail[G.A[a][0]].push_back(G.A[a]); // add arc to set of arcs with same tail
				}
			}
			// print3DVector(AByTail, "AByTail", "Tail");

			// then fill the half bins one by one by finding a path every iteration
			vector<int> halfBin; // The items packed in the current half bin 
			vector<vector<vector<int>>> standardHalfBinsPerEnd(R + 1); // The half bins not containing a reflected arc per end vertex 
			vector<vector<vector<int>>> reflectedHalfBinsPerEnd(R + 1); // The half bins containing a reflected arc per end vertex
			int tail, head, idx;
			vector<int> a;
			bool reflectUsed; // indicates if a reflected arc has been used on the path or not
			for (int b = 0; b < 2 * sol.UB; b++) { // loop over the to-be-filled bins
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

			vector<int> fullBin; // a combination of 2 half bins
			vector<vector<vector<int>>> currentBin;
			for (int end = 0; end <= R; end++) {
				while (reflectedHalfBinsPerEnd[end].size() > 0) {
					currentBin.resize(inst.C);
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

	sol.timeT = getCPUTime() - start;		    // save the total time
	return sol;                                 // return the solution
}