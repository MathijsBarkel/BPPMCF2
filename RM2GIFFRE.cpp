#include "RM2GIFFRE.h"

void GraphRM2GIFF::print() {
	cout << "Graph G = (V, A):\n";
	cout << "number of vertices: " << Nvert << "\n";
	cout << "number of arcs: " << Narcs << "\n";
	print2DVector(V_I2N, "- V_I2N", "idx");
	print3DVector(V_N2I, "- V_N2I:\n   vertex (-1,0,0): 0, rest", "layer/item index");
	print2DVector(A, "- A (1D vertices)", "arc");

	// also print the arcs with all vertices in their (j,v,gamma) representation
	cout << "- A (3D vertices):\n";
	int j1, j2, v1, v2, gamma1, gamma2, k;
	for (int a = 0; a < A.size(); a++) {
		j1 = V_I2N[A[a][0]][0];
		j2 = V_I2N[A[a][1]][0];
		v1 = V_I2N[A[a][0]][1];
		v2 = V_I2N[A[a][1]][1];
		gamma1 = V_I2N[A[a][0]][2];
		gamma2 = V_I2N[A[a][1]][2];
		k = A[a][2];
		cout << "   arc " << a << ": [" << "(" << j1 << ", " << v1 << ", " << gamma1 << "), (" << j2 << ", " << v2 << ", " << gamma2 << "), " << k << "]\n";
	}
	cout << "\n";
}

GraphRM2GIFF graphConstructionRM2GIFF(const Instance& inst, bool terminal) {
	// This function is used to construct the RM2GIFF graph for the method RM2GIFF

	// initialize a graph
	GraphRM2GIFF G;

	// initialize G.V_N2I: size JxWx2
	// - if vertex (j,v,gamma) does exist:    G.V_N2I[j][v][gamma] gives a 1D index for vertex (j,v,gamma). 
	// - If vertex v_c doesn't exist: G.V_N2I[j][v][gamma] = -1.
	// - The starting vertex gets index 0
	G.V_N2I.resize(inst.J, vector<vector<int>>(inst.W+1, vector<int>(2, -1)));

	// initialize G.V_I2N: size ?x3 (push_back each time to determine the # rows)
	// - G.ITN[idx] gives the 3-D (j,v, gamma) representation of the idx'th vertex
	G.V_I2N.push_back({ -1,0,0 });
	G.Nvert = 1; // keeps track of the number of nodes added

	// construct vertices and item arcs layer by layer (i.e., item type by item type)
	int c, w, q;						// current item color (c), size (w) and quantity (q)
	int idxTail, idxHead, idxItem = -1; // the 1D index of the current tail and current head of a considered arc and the index of the item type
	bool lastOfColor;					// a boolean indicating whether the current item is the last of its color
	int vNew, gammaNew;

	for (int j = 0; j < inst.J; j++) {
		c = inst.countedItems[j][0]; w = inst.countedItems[j][1]; q = inst.countedItems[j][2];	// find color, size and quantity
		lastOfColor = ((j == inst.J - 1) || (inst.countedItems[j + 1][0] != c));				// determine if this is the last item of its color 
		for (int v = 0; v < inst.W; v++) {														// loop over all load values
			for (int gamma = 0; gamma <= 1; gamma++) {											// for gamma = 0 and gamma = 1
				if (j == 0) {																	// find the index of the tail
					if (v == 0 && gamma == 0) { idxTail = 0; }
					else { idxTail = -1; }
				}
				else {idxTail = G.V_N2I[j - 1][v][gamma]; }
				if (idxTail != -1) {															// if the tail exists
					for (int k = 0; k <= q; k++) {												// loop over copies of the same item type
						if ((k == 0) && (j == inst.J - 1)) { continue; }						// skip the loss arcs in the final layer
						vNew = v + k * w;														// find the new load
						if (vNew <= inst.W) {													// if the new vertex still 'fits'
							gammaNew = ((k > 0) || (gamma == 1)) * (1 - lastOfColor);			// find the new value of gamma
							idxHead = G.V_N2I[j][vNew][gammaNew];								// find the 1D-idx of the current head
							if (idxHead == -1) { 												// if the head didn't exist yet: create the vertex and arc, and try adding another copy of the same item type
								idxHead = G.Nvert;												// introduce a new vertex index
								G.V_N2I[j][vNew][gammaNew] = G.Nvert;							// add the vertex
								G.V_I2N.push_back({ j,vNew, gammaNew });
								G.Nvert++;														// increment the number of vertices
							}
							G.A.push_back({ idxTail, idxHead, k });								// add the arc (even if the head existed already)
						}
						else break;																// go to the next vertex if the head would not fit anymore
					}
				}
			}
		}
	}

	if (terminal) {
		// Add terminal node: (J, W, 0)
		G.V_I2N.push_back({ inst.J,inst.W,0 });
		G.Nvert++;
		for (int i = 1; i < G.Nvert - 1; i++) {
			G.A.push_back({ i, G.Nvert - 1, 0 });
		}
	}


	// Return the graph
	G.Narcs = G.A.size();
	return G;
}

GraphRM2GIFF graphConstructionRM2GIFFUnitLoss(const Instance& inst) {
	// This function is used to construct the RM2GIFF graph for the method RM2GIFF, using unit loss arcs and a terminal node
	// The differences are:
	// - all nodes in the final layer are grouped into 1 terminal node
	// - we also use loss arcs to the terminal node
	// - we also use loss arcs from vertices with first coordinate W


	// initialize a graph
	GraphRM2GIFF G;

	// initialize G.V_N2I: size JxWx2
	// - if vertex (j,v,gamma) does exist:    G.V_N2I[j][v][gamma] gives a 1D index for vertex (j,v,gamma). 
	// - If vertex v_c doesn't exist: G.V_N2I[j][v][gamma] = -1.
	// - The starting vertex gets index 0
	G.V_N2I.resize(inst.J, vector<vector<int>>(inst.W + 1, vector<int>(2, -1)));

	// initialize G.V_I2N: size ?x3 (push_back each time to determine the # rows)
	// - G.ITN[idx] gives the 3-D (j,v, gamma) representation of the idx'th vertex
	G.V_I2N.push_back({ -1,0,0 });
	G.Nvert = 1; // keeps track of the number of nodes added

	// construct vertices and item arcs layer by layer (i.e., item type by item type)
	int c, w, q;						// current item color (c), size (w) and quantity (q)
	int idxTail, idxHead, idxItem = -1; // the 1D index of the current tail and current head of a considered arc and the index of the item type
	bool lastOfColor;					// a boolean indicating whether the current item is the last of its color
	int vNew, gammaNew;

	for (int j = 0; j < inst.J; j++) {
		c = inst.countedItems[j][0]; w = inst.countedItems[j][1]; q = inst.countedItems[j][2];		// find color, size and quantity
		lastOfColor = ((j == inst.J - 1) || (inst.countedItems[j + 1][0] != c));					// determine if this is the last item of its color 
		if (j == inst.J - 1) {																		// for the last layer
			idxHead = G.Nvert;																		// introduce the terminal node: (J-1,W,0)
			G.V_N2I[inst.J - 1][inst.W][0] = G.Nvert;												// add the vertex
			G.V_I2N.push_back({ inst.J - 1,inst.W,0 });
			G.Nvert++;
		}
		for (int v = 0; v <= inst.W; v++) {															// loop over all load values
			for (int gamma = 0; gamma <= 1; gamma++) {												// for gamma = 0 and gamma = 1
				if (j == 0) {																		// find the index of the tail
					if (v == 0 && gamma == 0) { idxTail = 0; }
					else { idxTail = -1; }
				}
				else { idxTail = G.V_N2I[j - 1][v][gamma]; }
				if (idxTail != -1) {																// if the tail exists
					for (int k = 0; k <= q; k++) {													// loop over copies of the same item type
						vNew = v + k * w;															// find the new load
						if (vNew <= inst.W) {														// if the new vertex still 'fits'
							if (j == inst.J - 1) { G.A.push_back({ idxTail, G.Nvert - 1, k }); }
							else {
								gammaNew = ((k > 0) || (gamma == 1)) * (1 - lastOfColor);			// find the new value of gamma
								idxHead = G.V_N2I[j][vNew][gammaNew];								// find the 1D-idx of the current head
								if (idxHead == -1) { 												// if the head didn't exist yet: create the vertex and arc, and try adding another copy of the same item type
									idxHead = G.Nvert;												// introduce a new vertex index
									G.V_N2I[j][vNew][gammaNew] = G.Nvert;							// add the vertex
									G.V_I2N.push_back({ j,vNew, gammaNew });
									G.Nvert++;														// increment the number of vertices
								}
								G.A.push_back({ idxTail, idxHead, k });								// add the arc (even if the head existed already)
							}
						}
						else break;																	// go to the next vertex if the head would not fit anymore
					}
				}
			}
		}
	}

	// Return the graph
	G.Narcs = G.A.size();
	return G;
}

Solution solveRM2GIFFRE(const Instance& inst, bool warmStart, double timeLimit, bool terminal, bool unitLoss) {
	// This function solves an instance using our reimplemenation of the RM2-GIFF method by Mehrani et al.

	double start = getCPUTime();     // starting time
	Solution sol;                    // initialize Solution struct

	// method declaration
	sol.method = "RM2GIFFRE";
	if (unitLoss) { terminal = false; }
	if (terminal) sol.method += "(terminal)";
	if (unitLoss) sol.method += "(unitLoss)";
	if (warmStart) sol.method += "(W)";

	// create graph
	GraphRM2GIFF G;
	if (unitLoss) {G = graphConstructionRM2GIFFUnitLoss(inst);}
	else {G = graphConstructionRM2GIFF(inst, terminal); }
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
	int j, k, kMax = 0;
	for (int a = 0; a < G.A.size(); a++) {
		j = G.V_I2N[G.A[a][1]][0];	// find corresponding item type
		k = G.A[a][2];				// find the quantity of the item to be packed
		kMax = max(kMax, k);		// update the maximum k
		if (k > 0) { // for item arcs
			f[a] = model.addVar(0, inst.countedItems[j][2]/k, 0, GRB_INTEGER);
		}
		else if (k == 0) { // for loss arcs
			f[a] = model.addVar(0, inst.B, 0, GRB_INTEGER);
		}
	}
	model.update();

	// declare linear expressions
	vector<GRBLinExpr> fIn(G.Nvert, 0);    											 	 // the amount of flow entering each vertex
	vector<GRBLinExpr> fOut(G.Nvert, 0);   												 // the amount of flow leaving each vertex
	vector<GRBLinExpr> typeUsed(inst.countedItems.size(), 0); 							 // the amount of arcs used of each item type 
	GRBLinExpr colorFragmentation = 0; 													 // the amount of color fragmentation = the number of arcs with k>0 from a node with gamma=0
	vector<GRBLinExpr> fragmentationsPerColor(inst.C, 0);								 // the number of color fragmentations per color = the number of arcs with k>0 from a node with gamma=0 to a node with c[j] = that color

	// calculate the linear expressions	
	int head, tail, c, gammaTail;
	for (int a = 0; a < G.A.size(); a++) {   					// loop over all arcs
		tail = G.A[a][0], head = G.A[a][1], k = G.A[a][2];		// find tail, head and multiplicity
		fIn[head] += f[a];        								// inflow
		fOut[tail] += f[a];		    							// outflow
		j = G.V_I2N[head][0];									// item type
		if (j != inst.J) { c = inst.countedItems[j][0]; }		// item color
		gammaTail = G.V_I2N[tail][2];							// gamma value of tail
		if (k > 0) {											// if the arc is an item arc
			typeUsed[j] += k * f[a];							// number of items used of certain type
			if (gammaTail == 0) {								// if the tail has gamma=0
				colorFragmentation += f[a];						// amount of color fragmentation
				fragmentationsPerColor[c] += f[a];				// amount of transitions for color
			}
		}
	}
	model.update();

	// set the objective: minimize the amount of color fragmentation
	model.setObjective(colorFragmentation, GRB_MINIMIZE);

	// constraints 1: flow conservation
	for (int v = 1; v < G.Nvert; v++) {       		// loop over all vertices
		if (terminal || unitLoss) {
			if (v != G.Nvert - 1) {
				model.addConstr(fIn[v] == fOut[v]);			// inflow = outflow
			}
		}
		else {
			model.addConstr(fIn[v] >= fOut[v]);			// inflow >= outflow
		}
	}

	// constraints 2: limited number of bins
	model.addConstr(fOut[0] <= inst.B);

	// constraints 3: item type quantities
	for (int j = 0; j < inst.countedItems.size(); j++) {		// loop over all item types
		model.addConstr(typeUsed[j] == inst.countedItems[j][2]);
	}

	// (optional) constraints 4: lower bound on transitioning into layers for each color
	for (int c = 0; c < inst.C; c++) {								// loop over all colors
		model.addConstr(fragmentationsPerColor[c] >= inst.LBs[c]);	// mimimum number of fragmentations per color
	}

	model.update();

	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	// start with an initial feasible solution if warmStart = true
	if (warmStart && inst.feasSol.size() > 0 && not terminal && not unitLoss) {

		// Find a mapping from an items color and size to its index
		vector<vector<int>> cw2j(inst.C, vector<int>(inst.W + 1, -1));
		int c, w;
		for (int j = 0; j < inst.J; j++) {
			c = inst.countedItems[j][0]; w = inst.countedItems[j][1];
			cw2j[c][w] = j;
		}

		// Find a mapping that gives the 1D-index of an arc based on (jHead (=jItem), vTail, gammaTail, k)
		vector<vector<vector<vector<int>>>> arcInfoToIdx(inst.J, vector<vector<vector<int>>>(inst.W, vector<vector<int>>(2, vector<int>(kMax+1,-1))));
		int vTail, jHead;
		for (int a = 0; a < G.A.size(); a++) {   												// loop over all arcs
			tail = G.A[a][0], head = G.A[a][1], k = G.A[a][2];									// find tail, head and multiplicity
			jHead = G.V_I2N[head][0];  vTail = G.V_I2N[tail][1]; gammaTail = G.V_I2N[tail][2];	// find info relevant to arc
			arcInfoToIdx[jHead][vTail][gammaTail][k] = a;										// add the arc to the map
		}

		// Initialize a vector with initial values
		vector<int> fInit(G.A.size(), 0);

		int q, jTail, cHead, vHead, gammaHead, a;
		for (int b = 0; b < inst.feasSol.size(); b++) {						// for every bin

			// Temporary
			// cout << b << "\n";

			// print3DVector(inst.feasSol[b], "Bin " + to_string(b), "color");
			jTail = -1; vTail = 0; gammaTail = 0; jHead = 0;				// start the bin
			for (int c = 0; c < inst.C; c++) {								// for every color
				if (inst.feasSol[b][c].size() > 0) {						// if the color is present in the bin
					for (int i = 0; i < inst.feasSol[b][c].size(); i++) {	// loop over the items of that color in the bin (from large to small)
						w = inst.feasSol[b][c][i][0];						// find item size
						j = cw2j[c][w];										// find item type
						q = inst.feasSol[b][c][i][1];						// find the number of copies of that item in the bin

						// Add loss arcs if jHead isn't yet j
						while (jHead != j) {
							//if (jTail == -1) { tail = 0; }
							//else { tail = G.V_N2I[jTail][vTail][gammaTail]; }
							//head = G.V_N2I[jHead][vHead][gammaHead];
							a = arcInfoToIdx[jHead][vTail][gammaTail][0];
							// cout << "Arc added: " << a << "\n";
							fInit[a]++;

							// go to the next vertex
							vHead = vTail;
							cHead = inst.countedItems[jHead][0];
							gammaHead = gammaTail;
							if (jHead == inst.J - 1 || cHead != inst.countedItems[jHead + 1][0]) {
								gammaHead = 0;
							}
							jTail = jHead;
							vTail = vHead;
							gammaTail = gammaHead;
							jHead++;
						}

						// Add item arc
						//if (jTail == -1) { tail = 0; }
						//else { tail = G.V_N2I[jTail][vTail][gammaTail]; }
						//head = G.V_N2I[jHead][vHead][gammaHead];
						a = arcInfoToIdx[jHead][vTail][gammaTail][q];
						// cout << "Arc added: " << a << "\n";
						fInit[a]++;

						// go to the next vertex
						vHead = vTail + q * w;
						cHead = inst.countedItems[jHead][0];
						gammaHead = 1;
						if (jHead == inst.J - 1 || cHead != inst.countedItems[jHead + 1][0]) {
							gammaHead = 0;
						}
						jTail = jHead;
						vTail = vHead;
						gammaTail = gammaHead;
						jHead++;
					}
				}
			}
		}

		// set the initial values
		for (int a = 0; a < G.A.size(); a++) {
			f[a].set(GRB_DoubleAttr_Start, fInit[a]);
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

		// store the items assigned to each bin (sorted by color)
		vector<vector<vector<int>>> currentBin;

		// first store the used arcs by tail
		int fval;
		vector<vector<vector<int>>> AByTail(G.Nvert);
		for (int a = 0; a < G.A.size(); a++) {
			fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
			for (int k = 0; k < fval; k++) {
				AByTail[G.A[a][0]].push_back(G.A[a]); // add arc to set of arcs with same tail
			}
		}
		// print3DVector(AByTail, "AByTail", "Tail");
		// then fill the bins one by one by finding a path leaving (0,0,0) every iteration
		int tail, head, k, w;
		vector<int> a;
		int nBinsUsed = ceil(fOut[0].getValue() - EPSILON); // find the number of bins used
		for (int b = 0; b < nBinsUsed; b++) { // loop over the to-be-filled bins
			vector<vector<int>> itemsInBinPerColor(inst.C);
			tail = 0;										// start at (0,0,0)
			while (AByTail[tail].size() > 0) {				// while you can still move on
				a = AByTail[tail].back();					// select an arc
				head = a[1];								// find the head of the arc
				k = a[2];									// find the multiplicity of the arc
				if (k > 0) {								// for item arcs: add item to the bin
					j = G.V_I2N[head][0];					// find item type
					c = inst.countedItems[j][0];			// find item color
					w = inst.countedItems[j][1];			// find item size
					for (int i = 0; i < k; i++) {			// for every copy of the item
						itemsInBinPerColor[c].push_back(w);	// add the item
					}
				}
				AByTail[tail].pop_back();					// remove the used arc
				tail = head;								// move to the new point
			}

			currentBin.resize(inst.C);						// the bin contains a seperate vector per color
			for (int c = 0; c < inst.C; c++) {				// for every color
				if (itemsInBinPerColor[c].size() > 0) {		// if the color is contained in the bin
					currentBin[c] = sortAndCount(itemsInBinPerColor[c]); // group and count the items of the current color in the bin
				}
			}
			sol.binPacking.push_back(currentBin);			// add the last bin to the set of all bins
			itemsInBinPerColor.clear();						// start a new empty bin
			currentBin.clear();								// start a new empty bin
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

	// Print the LP relaxation solution
	//double fval;
	//for (int a = 0; a < G.A.size(); a++) {
	//	fval = modelRelaxed.getVars()[a].get(GRB_DoubleAttr_X); // find variable value
	//	if (fval > 0) {
	//		cout << "arc " << a << ": " << fval << "\n";
	//	}
	//}

	return sol;		// return the solution
}

Solution solveRM2GIFFSetupMehrani(const Instance& inst, bool warmStart, double timeLimit) {
	// This function solves an instance using our reimplemenation of the RM2-GIFF method by Mehrani et al.,
	// Like in the implementation of Mehrani et al., we use unit loss arcs here and all nodes in the final layer are combined in a single terminal node
	// Furthermore, for the warm start we assume that colors were not re-ordered and that item weights are ordered from smallest to largest

	double start = getCPUTime();     // starting time
	Solution sol;                    // initialize Solution struct

	// method declaration
	sol.method = "RM2GIFFSetupMehrani";
	if (warmStart) sol.method += "(W)";

	// create graph
	GraphRM2GIFF G;
	G = graphConstructionRM2GIFFUnitLoss(inst); 
	// G.print(); // print the graph

	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();                       	// remove Gurobi message
	// env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);         // create a new model

	// declaration of the variables for the model
	vector<GRBVar> f(G.A.size());
	// for a in A: fa = # arc a is used
	int j, k, kMax = 0;
	for (int a = 0; a < G.A.size(); a++) {
		j = G.V_I2N[G.A[a][1]][0];	// find corresponding item type
		k = G.A[a][2];				// find the quantity of the item to be packed
		kMax = max(kMax, k);		// update the maximum k
		if (k > 0) { // for item arcs
			f[a] = model.addVar(0, inst.countedItems[j][2] / k, 0, GRB_INTEGER);
		}
		else if (k == 0) { // for loss arcs
			f[a] = model.addVar(0, inst.B, 0, GRB_INTEGER);
		}
	}
	model.update();

	// declare linear expressions
	vector<GRBLinExpr> fIn(G.Nvert, 0);    											 	 // the amount of flow entering each vertex
	vector<GRBLinExpr> fOut(G.Nvert, 0);   												 // the amount of flow leaving each vertex
	vector<GRBLinExpr> typeUsed(inst.countedItems.size(), 0); 							 // the amount of arcs used of each item type 
	GRBLinExpr colorFragmentation = 0; 													 // the amount of color fragmentation = the number of arcs with k>0 from a node with gamma=0
	vector<GRBLinExpr> fragmentationsPerColor(inst.C, 0);								 // the number of color fragmentations per color = the number of arcs with k>0 from a node with gamma=0 to a node with c[j] = that color

	// calculate the linear expressions	
	int head, tail, c, gammaTail;
	for (int a = 0; a < G.A.size(); a++) {   					// loop over all arcs
		tail = G.A[a][0], head = G.A[a][1], k = G.A[a][2];		// find tail, head and multiplicity
		fIn[head] += f[a];        								// inflow
		fOut[tail] += f[a];		    							// outflow
		j = G.V_I2N[head][0];									// item type
		if (j != inst.J) { c = inst.countedItems[j][0]; }		// item color
		gammaTail = G.V_I2N[tail][2];							// gamma value of tail
		if (k > 0) {											// if the arc is an item arc
			typeUsed[j] += k * f[a];							// number of items used of certain type
			if (gammaTail == 0) {								// if the tail has gamma=0
				colorFragmentation += f[a];						// amount of color fragmentation
				fragmentationsPerColor[c] += f[a];				// amount of transitions for color
			}
		}
	}
	model.update();

	// set the objective: minimize the amount of color fragmentation
	model.setObjective(colorFragmentation, GRB_MINIMIZE);

	// constraints 1: flow conservation
	for (int v = 1; v < G.Nvert-1; v++) {       	// loop over all vertices
		model.addConstr(fIn[v] == fOut[v]);			// inflow = outflow
	}

	// constraints 2: limited number of bins
	model.addConstr(fOut[0] <= inst.B);

	// constraints 3: item type quantities
	for (int j = 0; j < inst.countedItems.size(); j++) {		// loop over all item types
		model.addConstr(typeUsed[j] == inst.countedItems[j][2]);
	}

	// (optional) constraints 4: lower bound on transitioning into layers for each color
	for (int c = 0; c < inst.C; c++) {								// loop over all colors
		model.addConstr(fragmentationsPerColor[c] >= inst.LBs[c]);	// mimimum number of fragmentations per color
	}

	model.update();

	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	// start with an initial feasible solution if warmStart = true
	if (warmStart && inst.feasSol.size() > 0) {

		// Find a mapping from an items color and size to its index
		vector<vector<int>> cw2j(inst.C, vector<int>(inst.W + 1, -1));
		int c, w;
		for (int j = 0; j < inst.J; j++) {
			c = inst.countedItems[j][0]; w = inst.countedItems[j][1];
			cw2j[c][w] = j;
		}

		// Find a mapping that gives the 1D-index of an arc based on (jHead (=jItem), vTail, gammaTail, k)
		vector<vector<vector<vector<int>>>> arcInfoToIdx(inst.J, vector<vector<vector<int>>>(inst.W+1, vector<vector<int>>(2, vector<int>(kMax + 1, -1))));
		int vTail, jHead;
		for (int a = 0; a < G.A.size(); a++) {   												// loop over all arcs
			tail = G.A[a][0], head = G.A[a][1], k = G.A[a][2];									// find tail, head and multiplicity
			jHead = G.V_I2N[head][0];  vTail = G.V_I2N[tail][1]; gammaTail = G.V_I2N[tail][2];	// find info relevant to arc
			arcInfoToIdx[jHead][vTail][gammaTail][k] = a;										// add the arc to the map
		}

		// Initialize a vector with initial values
		vector<int> fInit(G.A.size(), 0);

		int q, jTail, cHead, vHead, gammaHead, a;
		for (int b = 0; b < inst.feasSol.size(); b++) {						// for every bin

			// Temporary
			// cout << b << "\n";

			// print3DVector(inst.feasSol[b], "Bin " + to_string(b), "color");
			jTail = -1; vTail = 0; gammaTail = 0; jHead = 0;				// start the bin
			for (int c = 0; c < inst.C; c++) {								// for every color
				if (inst.feasSol[b][c].size() > 0) {						// if the color is present in the bin
					for (int i = inst.feasSol[b][c].size() - 1; i << inst.feasSol[b][c].size() >= 0; i--) {	// loop over the items of that color in the bin (from small to large)
						w = inst.feasSol[b][c][i][0];						// find item size
						j = cw2j[c][w];										// find item type
						q = inst.feasSol[b][c][i][1];						// find the number of copies of that item in the bin

						// Add loss arcs if jHead isn't yet j
						while (jHead != j) {
							//if (jTail == -1) { tail = 0; }
							//else { tail = G.V_N2I[jTail][vTail][gammaTail]; }
							//head = G.V_N2I[jHead][vHead][gammaHead];
							a = arcInfoToIdx[jHead][vTail][gammaTail][0];
							// cout << "Arc added: " << a << "\n";
							fInit[a]++;

							// go to the next vertex
							vHead = vTail;
							cHead = inst.countedItems[jHead][0];
							gammaHead = gammaTail;
							if (jHead == inst.J - 1 || cHead != inst.countedItems[jHead + 1][0]) {
								gammaHead = 0;
							}
							jTail = jHead;
							vTail = vHead;
							gammaTail = gammaHead;
							jHead++;
						}

						// Add item arc
						a = arcInfoToIdx[jHead][vTail][gammaTail][q];
						// cout << "Arc added: " << a << "\n";
						fInit[a]++;

						// go to the next vertex
						vHead = vTail + q * w;
						cHead = inst.countedItems[jHead][0];
						gammaHead = 1;
						if (jHead == inst.J - 1 || cHead != inst.countedItems[jHead + 1][0]) {
							gammaHead = 0;
						}

						jTail = jHead;
						vTail = vHead;
						gammaTail = gammaHead;

						jHead++;
					}
				}
			}
			// Go to the terminal node
			while (jHead <= inst.J-1) {
				a = arcInfoToIdx[jHead][vTail][gammaTail][0];
				// cout << "Loss arc added: " << a << "\n";
				fInit[a]++;

				// go to the next vertex
				vHead = vTail;
				cHead = inst.countedItems[jHead][0];
				gammaHead = gammaTail;
				if (jHead == inst.J - 1 || cHead != inst.countedItems[jHead + 1][0]) {
					gammaHead = 0;
				}
				jTail = jHead;
				vTail = vHead;
				gammaTail = gammaHead;
				jHead++;
			}
		}

		// set the initial values
		for (int a = 0; a < G.A.size(); a++) {
			f[a].set(GRB_DoubleAttr_Start, fInit[a]);
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

		// store the items assigned to each bin (sorted by color)
		vector<vector<vector<int>>> currentBin;

		// first store the used arcs by tail
		int fval;
		vector<vector<vector<int>>> AByTail(G.Nvert);
		for (int a = 0; a < G.A.size(); a++) {
			fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
			for (int k = 0; k < fval; k++) {
				AByTail[G.A[a][0]].push_back(G.A[a]); // add arc to set of arcs with same tail
			}
		}
		// print3DVector(AByTail, "AByTail", "Tail");
		// then fill the bins one by one by finding a path leaving (0,0,0) every iteration
		int tail, head, k, w;
		vector<int> a;
		int nBinsUsed = ceil(fOut[0].getValue() - EPSILON); // find the number of bins used
		for (int b = 0; b < nBinsUsed; b++) { // loop over the to-be-filled bins
			vector<vector<int>> itemsInBinPerColor(inst.C);
			tail = 0;										// start at (0,0,0)
			while (AByTail[tail].size() > 0) {				// while you can still move on
				a = AByTail[tail].back();					// select an arc
				head = a[1];								// find the head of the arc
				k = a[2];									// find the multiplicity of the arc
				if (k > 0) {								// for item arcs: add item to the bin
					j = G.V_I2N[head][0];					// find item type
					c = inst.countedItems[j][0];			// find item color
					w = inst.countedItems[j][1];			// find item size
					for (int i = 0; i < k; i++) {			// for every copy of the item
						itemsInBinPerColor[c].push_back(w);	// add the item
					}
				}
				AByTail[tail].pop_back();					// remove the used arc
				tail = head;								// move to the new point
			}

			currentBin.resize(inst.C);						// the bin contains a seperate vector per color
			for (int c = 0; c < inst.C; c++) {				// for every color
				if (itemsInBinPerColor[c].size() > 0) {		// if the color is contained in the bin
					currentBin[c] = sortAndCount(itemsInBinPerColor[c]); // group and count the items of the current color in the bin
				}
			}
			sol.binPacking.push_back(currentBin);			// add the last bin to the set of all bins
			itemsInBinPerColor.clear();						// start a new empty bin
			currentBin.clear();								// start a new empty bin
		}
	}

	sol.timeT = getCPUTime() - start;	// save the total time

	return sol;		// return the solution
}
