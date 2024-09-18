#include "IP2RE.h"

Solution solveIP2RE(const Instance& inst, bool warmStart, double timeLimit, bool constrL2, int constrSym) {
    // This function solves a BPPMCF instance using a reimplementation of the IP2 model by Mehrani et al. (2022)

    double start = getCPUTime();     // starting time
    Solution sol;                    // initialize Solution struct

    // method declaration
    sol.method = "IP2RE";
    if (warmStart) sol.method += "(W)";
    if (!constrL2) sol.method += " (-L2constr)";
    if (constrSym == 1) sol.method += " (+symCols)";
    else if (constrSym == 2) sol.method += " (+symItems)";
    else if (constrSym == 3) sol.method += " (+symLoads)";

    // create a model
    GRBEnv env = GRBEnv();              	// create an environment
    removeLine();                           // remove Gurobi message
    // env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
    env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
    GRBModel model = GRBModel(env);         // create a new model

    // declaration of the variables for the model
    // for j = 0,...,J-1, b = 0,...,B-1: xjb = # items of type j (i.e., color c and size w) packed in bin b
    vector<vector<GRBVar> > x;
    x.resize(inst.J, vector<GRBVar>(inst.B));
    // for c = 0,...,C-1, b = 0,...,B-1: ycb = 1 if color c contained in bin b, 0 otherwise
    vector<vector<GRBVar>> y;
    y.resize(inst.C, vector<GRBVar>(inst.B));

    // initialization of the variables for the model
    for (int b = 0; b < inst.B; b++) {
        for (int j = 0; j < inst.J; j++) {
            x[j][b] = model.addVar(0, inst.countedItems[j][2], 0, GRB_INTEGER);
        }
        for (int c = 0; c < inst.C; c++) {
            y[c][b] = model.addVar(0, 1, 0, GRB_BINARY);
        }
    }
    model.update();

    // add the objective function and the contraints
    // the objective: minimize color fragmentation
    GRBLinExpr obj = 0;                 	// the amount of color fragmentation
    for (int c = 0; c < inst.C; c++) {      // loop over all colors
        for (int b = 0; b < inst.B; b++) {  // loop over all bins
            obj += y[c][b];
        }
    }
    model.setObjective(obj, GRB_MINIMIZE);

    // constraints 1: all copies of each item type must be packed
    vector<GRBLinExpr> assigned(inst.J, 0);   // the number of bins that each item type is assigned to
    for (int j = 0; j < inst.J; j++) {        // loop over all item types
        for (int b = 0; b < inst.B; b++) {    // loop over all bins
            assigned[j] += x[j][b];
        }
        model.addConstr(assigned[j] == inst.countedItems[j][2]);
    }

    // constraints 2: the bin capacities must be respected
    vector<GRBLinExpr> totalSize(inst.B, 0);    // the sum of item sizes packed in each bin
    for (int b = 0; b < inst.B; b++) {          // loop over all bins
        for (int j = 0; j < inst.J; j++) {      // loop over all item types
            totalSize[b] += inst.countedItems[j][1] * x[j][b];
        }
        model.addConstr(totalSize[b] <= inst.W);
    }

    // constraints 3: the y-variables must be set to the correct values
    int maxNumberOfItems;
    for (int j = 0; j < inst.J; j++) {                                                              // loop over all item types
        for (int b = 0; b < inst.B; b++) {                                                          // loop over all bins
            // maxNumberOfItems = inst.countedItems[j][2];                                          // original version
            maxNumberOfItems = min(inst.countedItems[j][2], inst.W / inst.countedItems[j][1]);      // improved version
            model.addConstr(x[j][b] <= maxNumberOfItems * y[inst.countedItems[j][0]][b]);
        }
    }

    // constraints 4 (optional): fraction-based LB on the y-variables
    vector<vector<GRBLinExpr>> totalSizePerColor;  // the sum of item sizes of certain color packed in each bin
    totalSizePerColor.resize(inst.C, vector<GRBLinExpr>(inst.B, 0));
    for (int b = 0; b < inst.B; b++) {          // loop over all bins
        for (int j = 0; j < inst.J; j++) {      // loop over all item types
            totalSizePerColor[inst.countedItems[j][0]][b] += inst.countedItems[j][1] * x[j][b];
        }
        for (int c = 0; c < inst.C; c++) {      // loop over all colors
            model.addConstr(totalSizePerColor[c][b] <= inst.W * y[c][b]);
        }
    }
    // constraints 5 (optional): L2-based LB on the number of bins required for each color
    if (constrL2) {
        vector<GRBLinExpr> binsPerColor(inst.C, 0);   // the number of bins allowed to contain each color 
        for (int c = 0; c < inst.C; c++) {            // loop over all colors
            for (int b = 0; b < inst.B; b++) {        // loop over all bins
                binsPerColor[c] += y[c][b];
            }
            model.addConstr(binsPerColor[c] >= inst.LBs[c]);
        }
    }
    // constraints 6A (optional): non-increasing number of colors per bin
    if (constrSym == 1) {
        vector<GRBLinExpr> colorsPerBin(inst.B, 0); // counts the number of colors in each bin
        for (int b = 0; b < inst.B; b++) {          // loop over all bins
            for (int c = 0; c < inst.C; c++) {      // loop over all colors
                colorsPerBin[b] += y[c][b];
            }
            if (b != 0) {
                model.addConstr(colorsPerBin[b - 1] >= colorsPerBin[b]);
            }
        }
    }
    // constraints 6B (optional): non-increasing number of items per bin
    else if (constrSym == 2) {
        vector<GRBLinExpr> itemsPerBin(inst.B, 0);  // counts the number of items in each bin
        for (int b = 0; b < inst.B; b++) {          // loop over all bins
            for (int j = 0; j < inst.J; j++) {      // loop over all item types
                itemsPerBin[b] += x[j][b];
            }
            if (b != 0) {
                model.addConstr(itemsPerBin[b - 1] >= itemsPerBin[b]);
            }
        }
    }
    // constraints 6C (optional): non-increasing load per bin
    else if (constrSym == 3) {
        for (int b = 1; b < inst.B; b++) {          // loop over all bins
            model.addConstr(totalSize[b - 1] >= totalSize[b]);
        }
    }

    model.update();

    sol.timeP = getCPUTime() - start;     // save the preprocessing time

    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    // model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
    model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	// if a warm start is used, start with the provided feasible solution
	if (warmStart && inst.feasSol.size() > 0) {

        // find a mapping from an items color and size to its index
        vector<vector<int>> cw2j(inst.C, vector<int>(inst.W + 1, -1));
        int c, w;
        for (int j = 0; j < inst.J; j++) {
            c = inst.countedItems[j][0]; w = inst.countedItems[j][1];
            cw2j[c][w] = j;
        }
        vector<vector<int>> xInit(inst.J, vector<int>(inst.B, 0));
        vector<vector<bool>> yInit(inst.C, vector<bool>(inst.B, false));
        int q, j;
		for (int b = 0; b < inst.feasSol.size(); b++) {	// loop over all bins
			for (int c = 0; c < inst.C; c++) {			// loop over all colors
				if (inst.feasSol[b][c].size() > 0) {	// if the color is contained in the bin
                    yInit[c][b] = true;
                    for (int i = 0; i < inst.feasSol[b][c].size(); i++) {
                        w = inst.feasSol[b][c][i][0];
                        q = inst.feasSol[b][c][i][1];
                        j = cw2j[c][w];
                        xInit[j][b] += q;
                    }
				}
			}
		}

		// set the initial values
		for (int b = 0; b < inst.B; b++) {
            for (int c = 0; c < inst.C; c++) {
                y[c][b].set(GRB_DoubleAttr_Start, yInit[c][b]);
            }
			for (int j = 0; j < inst.J; j++) {
				x[j][b].set(GRB_DoubleAttr_Start, xInit[j][b]);
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
        int x_jb, c;
        vector<int> item;
        for (int b = 0; b < inst.B; b++) {
            currentBin.resize(inst.C);
            for (int j = 0; j < inst.J; j++) {
                x_jb = ceil(x[j][b].get(GRB_DoubleAttr_X) - EPSILON);
                if (x_jb > 0) {
                    c = inst.countedItems[j][0];
                    item = { inst.countedItems[j][1], x_jb };
                    currentBin[c].push_back(item);
                }
            }
            sol.binPacking.push_back(currentBin);
            currentBin.clear();
        }
    }

    sol.timeT = getCPUTime() - start; // save the total time

    // Seperately solve the LP relaxation
    if (not warmStart) {
        GRBModel modelRelaxed = model.relax();
        model.reset(1);
        modelRelaxed.optimize();
        if (modelRelaxed.get(GRB_IntAttr_Status) == GRB_OPTIMAL) { sol.LPrel = modelRelaxed.get(GRB_DoubleAttr_ObjVal); }
        sol.timeLP = getCPUTime() - start - sol.timeT;
    }

    return sol;                       // return the solution
}
