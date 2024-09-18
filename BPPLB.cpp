#include "BPPLB.h"

Solution solveBPPLB(const Instance& inst, vector<int> &LBs, bool concentration, bool reflect, bool lexicographic) {
	// This function applies the algorithm BPP-LB to the given instance

	int timeLimit = 1800;		           // time limit in seconds
	Solution sol;						   // initialize solution object			
	vector<Solution> subsolutions(inst.C); // we will solve a seperate CSP problem for each color
	vector<int> binColors;				   // store the color of the items in each bin
	sol.method = "BPP-LB";				   // declare the name of the method
	if (!concentration) sol.method += "(NoConcentration)";
	if (!reflect) sol.method += "(NoReflect)";
	if (lexicographic) sol.method += "(Lexicographic)"; // NEW
	sol.opt = 1;
	sol.feas = 1;

	LBs.resize(inst.C); // an array that will hold the L^*_c bounds

	// Phase 1: Solve a CSP subproblem for each color
	for (int c = 0; c < inst.C; c++) { 	// for each color
		vector<vector<int>> itemsSubset = inst.groupedItems[c]; // select only items of current color

		// solve a CSP using AF for only the items of the current color
		// also pass to the function the L2 LB for the current color
		if (concentration && reflect) subsolutions[c] = solveAFCSPConcentration(inst, itemsSubset, inst.LBs[c], c, true, -1, lexicographic);
		else if (concentration && not reflect) subsolutions[c] = solveAFCSPConcentration(inst, itemsSubset, inst.LBs[c], c, false, -1, lexicographic);
		else if (not concentration && reflect) subsolutions[c] = solveReflectCSP(inst, itemsSubset, inst.LBs[c], c);
		else if (not concentration && not reflect) subsolutions[c] = solveAFCSP(inst, itemsSubset, inst.LBs[c], c);		
		
		// subsolutions[c].print(true); // print the solutions per color

		// update the L^2_c bounds to the L^*_c bounds
		if (subsolutions[c].opt) LBs[c] = subsolutions[c].LB;

		// combine all results
		sol.feas = min(sol.feas, subsolutions[c].feas);
		sol.opt = min(sol.opt, subsolutions[c].opt);
		sol.LB += subsolutions[c].LB;
		sol.timeP += subsolutions[c].timeP;
		sol.timeT += subsolutions[c].timeT;
		sol.Nvar += subsolutions[c].Nvar;
		sol.Nconstr += subsolutions[c].Nconstr;
		sol.Ncoeff += subsolutions[c].Ncoeff;
		for (int b = 0; b < subsolutions[c].LB; b++) {
			binColors.push_back(c);
		}

		// If the time limit is already reached after a subproblem, break inmediately
		if (sol.timeT > timeLimit) {
			sol.feas = -1, sol.opt = 0, sol.Nvar = 0, sol.Nconstr = 0, sol.Ncoeff = 0;
			return sol;
		}
	}

	// Phase 2: (try to) combine the uni-colored bins in order to use at most B bins
	// if all subsolutions were feasible and optimal
	if (sol.feas == 1 && sol.opt) {
		// put all bins next to each other 
		sol.binPacking.resize(sol.LB);
		vector<int> binLoads(sol.LB, 0);					// a vector that is to contain the total size in each bin
		vector<vector<vector<vector<int>>>> subPackingc;	// store the packing for the subinstance with only color c
		int b = 0;
		for (int c = 0; c < inst.C; c++) {								// for every color
			subPackingc = subsolutions[c].binPacking;					// the packing for the subinstance with only color c
			for (int bc = 0; bc < subPackingc.size(); bc++) {			// for every bin in the current sub-packing
				sol.binPacking[b + bc] = subPackingc[bc];				// add the bin to the set of all bins
				for (int j = 0; j < subPackingc[bc][c].size(); j++) {	// for every item type used in the bin
					binLoads[b + bc] += subPackingc[bc][c][j][1] * subPackingc[bc][c][j][0]; // add the size to the total load
				}
			}
			b += subPackingc.size();
		}
		// printVector(binLoads, "binLoads");

		// if we respect the bin limit, the solution is optimal
		if (sol.binPacking.size() <= inst.B) {
			sol.feas = 1; sol.UB = sol.LB;
		}
		// otherwise, we try to combine bins. The bin loads become the items
		else {
			// first sort and count the items
			vector<vector<int>> binItems = sortAndCount(binLoads);
			// solve the CSP using AF with the bin loads as items
			Solution solCombined;
			if (reflect) solCombined = solveReflectCSP(inst, binItems, 0, 0, inst.B);
			else solCombined = solveAFCSP(inst, binItems, 0, 0, inst.B);

			// update the attributes of sol (note that the LB remains unchanged)
			sol.feas = min(sol.feas, solCombined.feas);
			sol.opt = min(sol.opt, solCombined.opt);
			sol.timeP += solCombined.timeP;
			sol.timeT += solCombined.timeT;
			sol.Nvar += solCombined.Nvar;
			sol.Nconstr += solCombined.Nconstr;
			sol.Ncoeff += solCombined.Ncoeff;
			vector<vector<vector<vector<int>>>> combinedPacking = solCombined.binPacking;

			// if we still exceed the bin limit we failed. However, the solution does give a valid LB
			if (combinedPacking.size() > inst.B) {
				sol.opt = 0, sol.feas = -1;
			}
			// if we now do respect the bin limit, the solution is optimal
			else {
				sol.feas = 1; sol.opt = 1; sol.UB = sol.LB;
				// we must unmerge the bins again
				// first we find out which original bins correspond to each load value
				vector<vector<int>> binsPerLoad(inst.W + 1); // each entry consistst the bin indices of bins containing that load
				for (int b = 0; b < binLoads.size(); b++) {  // loop over all bins
					binsPerLoad[binLoads[b]].push_back(b);   // add the bin to the list of bins corresponding to the load of the bin
				}
				// declare variables
				vector<vector<vector<vector<int>>>> newBinPacking;
				vector<vector<vector<int>>> currentBin(inst.C);
				vector<vector<int>> currentCombinedBin;
				int w, q, b;
				// now unmerge the bins
				for (int new_b = 0; new_b < combinedPacking.size(); new_b++) { // for each newly created bin
					currentBin.resize(inst.C);								   // the bin consists of a seperate vector for each color
					currentCombinedBin = combinedPacking[new_b][0];			   // extract the bin's contents: size and quantity of each 'bin' item.
					for (int j = 0; j < currentCombinedBin.size(); j++) {      // loop over all 'bin' items
						w = currentCombinedBin[j][0], q = currentCombinedBin[j][1]; // find the size and the quantity
						for (int k = 0; k < q; k++) {								// repeat for every copy of the 'bin' item with the current size
							b = binsPerLoad[w].back();								// find a bin with load w
							// note that bin b has color binColors[b]
							currentBin[binColors[b]] = sol.binPacking[b][binColors[b]]; //add the item to the new bin under the correct color
							binsPerLoad[w].pop_back();								// remove the 'bin' item from the list
						}
					}
					newBinPacking.push_back(currentBin); // add the bin to the set of all bins
					currentBin.clear();					 // start a new bin
				}
				sol.binPacking = newBinPacking;			 // save the packing
			}
		}
	}
	// return the solution
	return sol;
}