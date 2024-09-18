#include "LF.h"

void Graph2D::print() {
	cout << "Graph G = (V, A):\n";
	//cout << "number of vertices: " << sol.G.Nvert << "\n";
	//cout << "number of arcs: " << sol.G.Narcs << "\n";
	// print2DVector(V_I2N, "- V_I2N", "idx");
	print2DVector(V_N2I, "- V_N2I:\n   vertex s: 0, rest", "layer/color");
	// print2DVector(A, "- A (1D vertices)", "arc");

	// also print the arcs with all vertices in their v_c representation
	cout << "- A (2D vertices):\n";
	int u, v, c, d, idx;
	for (int a = 0; a < A.size(); a++) {
		u = V_I2N[A[a][0]][1];
		c = V_I2N[A[a][0]][0];
		v = V_I2N[A[a][1]][1];
		d = V_I2N[A[a][1]][0];
		idx = A[a][2];
		if (c == -1) {
			cout << "   arc " << a << ": (" << "s" << ", " << v << "_" << d << ", " << idx << ")\n";
		}
		else if (d == -1) {
			cout << "   arc " << a << ": (" << u << "_" << c << ", " << "t" << ", " << idx << ")\n";
		}
		else {
			cout << "   arc " << a << ": (" << u << "_" << c << ", " << v << "_" << d << ", " << idx << ")\n";
		}
	}
	cout << "\n";
}

vector<int> determineColorOrder(const Instance& inst) {
	vector<vector<bool>> V(inst.C, vector<bool>(inst.W, false));
	vector<vector<int>> wq;
	int w, q;

	// create the standard AF graph for each color seperately
	for (int c = 0; c < inst.C; c++) {
		wq = inst.groupedItems[c];     // 2D array holding the sizes and quantities of all items of the current color										// current item size (w) and quantity (q)
		V[c][0] = true;				   // start with vertex 0 active
		for (int j = 0; j < wq.size(); j++) {	          // loop over item types (ordered from large to small size)
			w = wq[j][0], q = wq[j][1];		              // find size (w) and quantity (q) of current item type
			for (int v = inst.W - w; v >= 0; v--) {       // loop over vertices in backwards order
				if (V[c][v]) {							  // only process vertices that were already there
					for (int k = 0; k < q; k++) {		  // loop over copies of the same item type
						if (v + (k + 1) * w < inst.W) {	  // if the new vertex still 'fits'
							if (!V[c][v + (k + 1) * w]) { // if the head did not exist yet
								V[c][v + (k + 1) * w] = true; // add the new vertex
							}
							else break;					  // go to the next item type if the head did exist already  
						}
						else break;                       // also go to the next item type if the next copy does not fit anymore
					}
				}
			}
		}
	}

	// Greedy idea: Start by making the standard graph for each color seperately.
	// Assign a score to each layer. Pick the best one. That is layer 1.
	// Then shift all other layers to obtain their activated nodes given the new startpoints due to transition arcs from layer 1 
	// Give all layers a score again, and pick the one with the best score. This will be layer 2, etc.
	// If at one point in the algorithm all nodes have been activated --> terminate

	// Main variables
	vector<int> colorOrder;
	vector<bool> alreadyPlaced(inst.C, false);
	vector<int> shifts;
	vector<int> scores;
	int best;
	vector<vector<bool>> Vnew(inst.C, vector<bool>(inst.W, false));
	vector<bool> full(inst.W, true);
	int maxScore = inst.W * (inst.W - 1) / 2;

	for (int l = 0; l < inst.C-1; l++) {

		// 1. Give each layer a score: Sum_{v=1, ..., W-1 | v activated} (W-v)
		scores.clear(); scores.resize(inst.C, 0);
		for (int c = 0; c < inst.C; c++) {
			for (int v = 1; v < inst.W; v++) {
				if (V[c][v]) { scores[c] += inst.W - v; }
			}
		}

		// 2. Choose the best layer based on the weighted node score
		best = distance(scores.begin(), min_element(scores.begin(), scores.end()));
		if (scores[best] == maxScore) break;
		alreadyPlaced[best] = true;
		colorOrder.push_back(best);

		// 3. Find the shift amounts
		shifts.clear();
		for (int v = 1; v < inst.W; v++) {
			if (V[best][v]) {
				shifts.push_back(v);
			}
		}

		// 4. Find the new activated vertices per layer
		for (int c = 0; c < inst.C; c++) {
			if (alreadyPlaced[c]) { // set already processed layers to all 1's --> never the lowest score anymore
				Vnew[c] = full;
			}
			else {
				Vnew[c] = V[c];
				for (int v = 0; v < inst.W; v++) {
					if (V[c][v]) {
						for (auto& s : shifts) {
							if (v + s < inst.W) {
								Vnew[c][v + s] = true;
							}
							else break;
						}
					}
				}
			}
		}
		V = Vnew;
	}
	// add the final color (or final colors in case of terminating early)
	for (int c = 0; c < inst.C; c++) {
		if (!alreadyPlaced[c]) colorOrder.push_back(c);	
	}

	return colorOrder;
}

Graph2D graphConstructionLF(const Instance& inst, int transitionType) {
    // This function is used to construct the AF graph for the method LayerFlow
	// Depending on whether transitionType is set to 0, 1, or 2,
	//   different types of transition arcs are constructed

    // initialize a graph
    Graph2D G;

    // initialize G.V_N2I: size CxW
    // - if vertex v_c does exist:    G.V_N2I[c][v] gives a 1D index for vertex v_c. 
    // - If vertex v_c doesn't exist: G.V_N2I[c][v] = -1.
    // - The special vertex s is not encoded in G.V_N2I
    G.V_N2I.resize(inst.C, vector<int>(inst.W+1, -1));

    // initialize G.V_I2N: size ?x2 (push_back each time to determine the # rows)
    // - G.ITN[idx] gives the 2-D (c,v) representation of the idx'th vertex
    // vertex s (= (-1,0)) is node 0
    G.V_I2N.push_back({ -1,0 });
    G.Nvert = 1; // keeps track of the number of nodes added

    // initialize V: size W
    // - V[v] = true iff the nodes of size v have been activated in the current layer
    //                                          (and thus also all subsequent layers)
    vector<bool> V(inst.W+1, false); V[0] = true;

    // Construct vertices and item arcs layer by layer (i.e., color by color)
    vector<vector<int>> wq;			 // 2D array holding the sizes and quantities of all items of the current color
    int w, q;							// current item size (w) and quantity (q)
    int idxTail, idxHead, idxItem = -1; // the 1D index of the current tail and current head of a considered arc and the index of the item type

    for (int c = 0; c < inst.C; c++) {                // construct the graph per layer (i.e., per color)
        // set all vertices with sizes already activated to active in the new layer
        for (int v = 0; v < V.size(); v++) {
            if (V[v]) {
                G.V_N2I[c][v] = G.Nvert;
                G.V_I2N.push_back({ c,v });
                G.Nvert++;
            }
        }
        wq = inst.groupedItems[c];						// select the items of the current color

        for (int j = 0; j < wq.size(); j++) {	          // loop over item types (ordered from large to small size)
            idxItem++;                                    // update the index of the current item
            w = wq[j][0], q = wq[j][1];		              // find size (w) and quantity (q) of current item type
            for (int v = inst.W - w; v >= 0; v--) {       // loop over vertices in backwards order
                idxTail = G.V_N2I[c][v];                  // find the 1D idx of the current tail
                if (idxTail != -1) {				      // only process vertices that were already there
                    for (int k = 0; k < q; k++) {		  // loop over copies of the same item type
                        // Case 1: the new vertex still 'fits'
                        if (v + (k + 1) * w <= inst.W) {
                            idxHead = G.V_N2I[c][v + (k + 1) * w]; // find the 1D idx of the current head
                            // if the head didn't exist yet: create the vertex and arc, and try adding another copy of the same item type
                            if (idxHead == -1) {          
                                V[v + (k + 1) * w] = true;
                                idxHead = G.Nvert;
                                G.A.push_back({ idxTail, idxHead, idxItem });
                                G.V_N2I[c][v + (k + 1) * w] = G.Nvert;
                                G.V_I2N.push_back({ c,v + (k + 1) * w });
                                G.Nvert++;
                                idxTail = idxHead; // move to the next point
                            }
                            // if the head did exist yet: create just the arc and move on to the next item type
                            else {
                                G.A.push_back({ idxTail, idxHead, idxItem });
                                break; 
                            }
                        }

                        // Case 2: the new vertex does not 'fit' anymore --> go to the next item type
                        else break;
                    }
                }
            }
        }
    }

	// If we add a transition arc from every vertex to the vertices of the same value in all subsequent layers
	if (transitionType == 0) {
		// Add special arcs
		for (int c = 0; c < inst.C; c++) {
			G.A.push_back({ 0, G.V_N2I[c][0], -1 }); // starting arc (= special case of transition arc)
			for (int v = 1; v < inst.W; v++) {
				if (G.V_N2I[c][v] != -1) {
					// G.A.push_back({ G.V_N2I[c][v], 1, -2 }); // loss arc
					for (int d = c+1; d < inst.C; d++) {
						G.A.push_back({ G.V_N2I[c][v], G.V_N2I[d][v], -1 }); // transition arc
					}
				}
			}
		}
	}
	// If we add a transition arc from every vertex to and from the transition vertex of the same value 
	else if (transitionType == 1) {
		// Add the transition layer
		G.V_N2I.push_back(vector<int>(inst.W, -1));
		for (int v = 1; v < inst.W; v++) {
			if (V[v]){
				G.V_N2I[inst.C][v] = G.Nvert;
				G.V_I2N.push_back({inst.C, v});
				G.Nvert += 1;
			}
		}
		// Add special arcs
		for (int c = 0; c < inst.C; c++) {
			G.A.push_back({ 0, G.V_N2I[c][0], -1 }); // starting arc (= special case of transition arc)
			for (int v = 1; v < inst.W; v++) {
				if (G.V_N2I[c][v] != -1) {
					// G.A.push_back({ G.V_N2I[c][v], 1, -2 }); // loss arc
					if (c != inst.C - 1) {
						G.A.push_back({ G.V_N2I[c][v], G.V_N2I[inst.C][v], -3 }); // transition arc 1 (doesn't count towards objective)
					}
					if (c != 0) {
						G.A.push_back({G.V_N2I[inst.C][v], G.V_N2I[c][v], -1 }); // transition arc 2 (counts towards objective)
					}
				}
			}
		}
	}
	// If we add a transition arc from every regular vertex to and from the unique associated transition vertex, 
	// ... and from every transition vertex to the transition vertex of the next color and same value 
	else if (transitionType == 2) {
		// Add the transition vertices
		// The transition vertex corresponding to regular vertex (c,v) is the vertex (C + c, v)
		for (int c = 0; c < inst.C; c++) {
			G.V_N2I.push_back(vector<int>(inst.W, -1));
			for (int v = 1; v < inst.W; v++) {
				if (G.V_N2I[c][v] != -1){
					G.V_N2I[inst.C + c][v] = G.Nvert;
					G.V_I2N.push_back({inst.C + c, v });
					G.Nvert += 1;
				}
			}
		}
		// Add special arcs
		for (int c = 0; c < inst.C; c++) {
			G.A.push_back({ 0, G.V_N2I[c][0], -1 }); // starting arc (= special case of transition arc)
			for (int v = 1; v < inst.W; v++) {
				if (G.V_N2I[c][v] != -1) {
					// G.A.push_back({ G.V_N2I[c][v], 1, -2 }); // loss arc
					if (c != 0) {
						G.A.push_back({ G.V_N2I[inst.C + c][v], G.V_N2I[c][v], -1 }); // transition arc 1 (counts towards objective)
					}
					if (c != inst.C - 1) {
						G.A.push_back({ G.V_N2I[c][v], G.V_N2I[inst.C + c][v], -3 }); // transition arc 2 (doesn't count)
						G.A.push_back({ G.V_N2I[inst.C + c][v], G.V_N2I[inst.C + c + 1][v], -4 }); // transition arc 3 (doesn't count)
					}
				}
			}
		}
	}

    // Return the graph
    return G;
}

Solution solveLF(const Instance& inst, bool warmStart, double timeLimit, int transitionType, bool constrL2, bool constrSym) {
	// This function solves an instance using the LF method.
	// Depending on whether transitionType is set to 0, 1, or 2, the transition arcs in the graph are constructed differently.
	//	(Only version 0 appears in the paper.)

	double start = getCPUTime();     // starting time
	Solution sol;                    // initialize Solution struct

	// method declaration
	sol.method = "LF";
	if (transitionType == 1) sol.method += "(A)";
	else if (transitionType == 2) sol.method += "(B)";
	if (!constrL2) sol.method += "(no L2 constr)";
	if (!constrSym) sol.method += "(no sym constr)";
	if (warmStart) sol.method +="(W)";

	// create graph
	Graph2D G = graphConstructionLF(inst, transitionType);
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
	int idx;
	for (int a = 0; a < G.A.size(); a++) {
		idx = G.A[a][2];
		if (idx >= 0) { // for item arcs
			f[a] = model.addVar(0, inst.countedItems[idx][2], 0, GRB_INTEGER);
		}
		else if (idx == -1) { // for transition arcs
			f[a] = model.addVar(0, inst.B * inst.C, 0, GRB_INTEGER);
		}
		else { // for loss arcs (idx == -2), or for type 1, transition arcs from transition vertex (idx == -3))
			f[a] = model.addVar(0, inst.B, 0, GRB_INTEGER);
		}
	}
	model.update();

	// declare linear expressions
	vector<GRBLinExpr> fIn(G.Nvert, 0);    											 	 // the amount of flow entering each vertex
	vector<GRBLinExpr> fOut(G.Nvert, 0);   												 // the amount of flow leaving each vertex
	vector<GRBLinExpr> typeUsed(inst.countedItems.size(), 0); 							 // the amount of arcs used of each item type 
	GRBLinExpr colorFragmentation = 0; 													 // the amount of color fragmentation = the amount of transition arcs used
	vector<GRBLinExpr> transitionsPerColor(inst.C, 0);									 // the amount of transition arcs going into each layer
	vector<vector<GRBLinExpr>> higherTransitions(inst.C, vector<GRBLinExpr>(inst.W, 0)); // the number of higher transition vertices going into the same corresponding transition vertex

	// calculate the linear expressions	
	int head, tail, d;
	for (int a = 0; a < G.A.size(); a++) {   					// loop over all arcs
		tail = G.A[a][0], head = G.A[a][1], idx = G.A[a][2];	// find tail, head and index
		fIn[head] += f[a];        								 // inflow
		fOut[tail] += f[a];		    							 // outflow
		if (idx >= 0) typeUsed[idx] += f[a];					 // number of items used of certain type
		else if (idx == -1) {
			colorFragmentation += f[a];							// amount of color fragmentation
			d = G.V_I2N[head][0];								// color of layer of head
			transitionsPerColor[d] += f[a];						// amount of transitions into layer c
		}
		else if (transitionType == 1 && idx == -3) {			// in case of type 1 transitions, we look for transition arcs from real vertices to corresponding transition vertices
			d = G.V_I2N[tail][0];								// color of layer of tail
			for (int c = d + 1; c < inst.C; c++) {				// for all higher colors
				higherTransitions[c][G.V_I2N[tail][1]] += f[a]; // number of higher transitions to same transition vertex
			}
		}
	}
	model.update();

	// set the objective: minimize the amount of color fragmentation
	model.setObjective(colorFragmentation, GRB_MINIMIZE);

	// constraints 1: flow conservation
	for (int v = 1; v < G.Nvert; v++) {       	// loop over all vertices
		if (fOut[v].size() > 0) { model.addConstr(fIn[v] >= fOut[v]); }		// inflow >= outflow
	}
	
	// constraints 2: limited number of bins
	model.addConstr(fOut[0] <= inst.B);
	
	// constraints 3: item type quantities
	for (int j = 0; j < inst.countedItems.size(); j++) {					// loop over all item types	
		model.addConstr(typeUsed[j] == inst.countedItems[j][2]);			// demand met
	}
	// (optional) constraints 4: lower bound on transitioning into layers for each color
	if (constrL2) {
		for (int c = 0; c < inst.C; c++) {							// loop over all colors
			model.addConstr(transitionsPerColor[c] >= inst.LBs[c]);	// mimimum number of transition arcs
		}
	}

	// (optional) constraints 5 (when transitionType = 1 only): reduce symmetry by partially ordering colors
	if (transitionType == 1 && constrSym) {
		int c, v;
		for (int a = 0; a < G.A.size(); a++) {
			if (G.A[a][2] == -1 && G.A[a][0] != 0 && G.A[a][1] != 1 ) {
				c = G.V_I2N[G.A[a][1]][0];
				if (c != 0) {
					v = G.V_I2N[G.A[a][1]][1];
					model.addConstr(f[a] <= higherTransitions[c][v]);
				}
			}
		}
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
		vector<vector<int>> AN2I(G.Nvert, vector<int>(G.Nvert, -1));
		vector<int> AT2I(G.Nvert);

		for (int a = 0; a < G.A.size(); a++) {   					// loop over all arcs
			tail = G.A[a][0], head = G.A[a][1], idx = G.A[a][2];	// find tail, head and index
			AN2I[tail][head] = a;
		}

		int cur_v, new_v, w, q, v1, v2, a, cur_c;
		vector<int> fInit(G.A.size(), 0);
		for (int b = 0; b < inst.feasSol.size(); b++) {		// for every bin
			cur_c = -1;
			for (int c = 0; c < inst.C; c++) {				// for every color
				if (inst.feasSol[b][c].size() > 0) {
					if (cur_c == -1) {						// starting arc
						v1 = 0;
						v2 = G.V_N2I[c][0];
						cur_v = 0;
						a = AN2I[v1][v2];					// find 1D index of arc
						fInit[a]++;							// add 1 to the flow on that arc
					}
					else {									// transition arc
						v1 = G.V_N2I[cur_c][cur_v];
						v2 = G.V_N2I[c][cur_v];
						a = AN2I[v1][v2];					// find 1D index of arc
						fInit[a]++;							// add 1 to the flow on that arc
					}

					cur_c = c;
					for (int j = 0; j < inst.feasSol[b][c].size(); j++) {			// for every item type
						w = inst.feasSol[b][c][j][0], q = inst.feasSol[b][c][j][1];
						for (int k = 0; k < q; k++) {								// for every copy of the item type
							new_v = cur_v + w;										// find new value
							v1 = G.V_N2I[c][cur_v];									// find 1D indices of vertices
							v2 = G.V_N2I[c][new_v];
							cur_v = new_v;											// move to the new point
							a = AN2I[v1][v2];										// find 1D index of arc
							fInit[a]++;												// add 1 to the flow on that arc
						}
					}
				}
			}
			//if (cur_v != inst.W && cur_c != -1) {	// add a loss arc to t if t wasn't reached yet
			//	v1 = G.V_N2I[cur_c][cur_v];			// find 1D index of the current vertex
			//	a = AT2I[v1];						// find 1D index of the loss arc
			//	fInit[a]++;							// add 1 to the flow on that arc
			//}
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

		// then fill the bins one by one by finding an s-t path every iteration
		int tail, head, idx;
		vector<int> a;
		int NbinsUsed = ceil(fOut[0].getValue() - EPSILON); // find the number of bins used
		for (int b = 0; b < NbinsUsed; b++) { // loop over the to-be-filled bins
			vector<vector<int>> itemsInBinPerColor(inst.C);
			tail = 0;						// start at s
			while (AByTail[tail].size() > 0) {	// keep going as long as there is an outgoing arc
				a = AByTail[tail].back();	// select an arc
				head = a[1];				// find the head of the arc
				idx = a[2];					// find the index of the arc
				if (idx >= 0) {				// for item arcs: add item to the bin
					itemsInBinPerColor[inst.countedItems[idx][0]].push_back(inst.countedItems[idx][1]);
				}
				AByTail[tail].pop_back();	// remove the used arc
				tail = head;				// move to the new point
			}

			currentBin.resize(inst.C);					// the bin contains a seperate vector per color
			for (int c = 0; c < inst.C; c++) {			// for every color
				if (itemsInBinPerColor[c].size() > 0) { // if the color is contained in the bin
					currentBin[c] = sortAndCount(itemsInBinPerColor[c]); // group and count the items of the current color in the bin
				}
			}
			sol.binPacking.push_back(currentBin);		// add the last bin to the set of all bins
			itemsInBinPerColor.clear();			        // start a new empty bin
			currentBin.clear();							// start a new empty bin
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

	//double fval;
	//for (int a = 0; a < G.A.size(); a++) {
	//	fval = modelRelaxed.getVars()[a].get(GRB_DoubleAttr_X); // find variable value
	//	if (fval > 0) {
	//		cout << "arc " << a << ": " << fval << "\n";
	//	}
	//}

	return sol;										// return the solution
}

Solution solveModelAfterOrdering(Instance inst, bool warmStart, double timeLimit, string model, int transitionType, bool constrL2, bool constrSym) {
	// This function is used to call either RM2-GIFF or LayerFlow after first reordering the colors. 
	// Afterwards, the original color order is restored, and the results are presented with respect to this original order.


	double start = getCPUTime();    // starting time
	double timePO;

	vector<int> colorOrder = determineColorOrder(inst);
	// printVector(colorOrder, "color order");
		
	// Relabel the colors so that the color index corresponds to the layer
	vector<int> old2newColor(inst.C);
	for (int idx = 0; idx < inst.C; idx++) {
		old2newColor[colorOrder[idx]] = idx;
	}
	for (int i = 0; i < inst.I; i++) {
		inst.items[i][0] = old2newColor[inst.items[i][0]];
	}
	inst.preprocessing(true);   // reset the instance attributes
	inst.findL2();
	// inst.print();		// print information about the instance

	if (warmStart) {
		vector<vector<vector<vector<int>>>> feasSol2(inst.B, vector<vector<vector<int>>>(inst.C));
		int cNew;
		for (int c = 0; c < inst.C; c++) {
			cNew = old2newColor[c];
			for (int b = 0; b < inst.feasSol.size(); b++) {
				feasSol2[b][cNew] = inst.feasSol[b][c];
			}
		}
		inst.feasSol = feasSol2;
	}

	// save the additional preprocessing time due to the color reordering
	timePO = getCPUTime() - start;

	Solution sol;
	if (model == "LF") {
		// Apply the LF method
		sol = solveLF(inst, warmStart, timeLimit, transitionType, constrL2, constrSym);
	}
	else if (model == "RM2GIFF") {
		// Apply the RM2GIFF method
		sol = solveRM2GIFFRE(inst, warmStart, timeLimit);
	}

	// Revert the colors back
	vector<vector<vector<vector<int>>>> binPacking(sol.binPacking.size(), vector<vector<vector<int>>>(inst.C));
	for (int b = 0; b < sol.binPacking.size(); b++) {
		for (int c = 0; c < inst.C; c++) {
			if (sol.binPacking[b][c].size() > 0) {
				binPacking[b][colorOrder[c]] = sol.binPacking[b][c];
			}
		}
	}
	sol.binPacking = binPacking;
	if (warmStart) {
		vector<vector<vector<vector<int>>>> feasSol(inst.B, vector<vector<vector<int>>>(inst.C));
		int cOld;
		for (int c = 0; c < inst.C; c++) {
			cOld = colorOrder[c];
			for (int b = 0; b < inst.B; b++) {
				feasSol[b][cOld] = inst.feasSol[b][c];
			}
		}
		inst.feasSol = feasSol;
	}

	sol.timeP += timePO;			  // update the pre-processing time

	sol.timeT = getCPUTime() - start; // save the total time

	return sol;
}
