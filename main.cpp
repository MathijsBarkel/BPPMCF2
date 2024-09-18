#include "main.h"

int main() {
	// -------------------------------------------------------------------------
	// Part 1: Instance selection and pre-processing

	// Select an instance
	// string filename = "InstancesBPPMCF/Dataset 1/70-8/70-8-1.txt";
	// string filename = "InstancesBPPMCF/Dataset 2/120-2/120-2-1.txt";
	// string filename = "InstancesBPPMCF/Dataset 3/10-100-4/10-100-4-1.txt";
	// string filename = "InstancesBPPMCF/Dataset 4/50-400-3/50-400-3-1.txt";
	string filename = "InstancesBPPMCF/Triplets/t60_00.txt";
	Instance inst = readInstance(filename);
	// Instance inst = createInstance("Test_instance", 10, 4, 3, 6, { {0,4},{0,3},{0,1},{1,3},{1,2},{1,2},{1,1},{2,3},{2,2},{2,1} });
	// Instance inst = createInstance("Test_instance", 19, 9, 4, 24, { {0,16}, {0,15}, {0,13}, {0,13}, {1,15}, {1,14}, {1,2}, {2,8}, {2,7}, {2,3}, {2,3}, {2,1}, {3,8}, {3,8}, {3,6}, {3,6}, {3,3}, {3,3}, {3,3} });
	// Instance inst = createInstance("Weird_instance", 18, 6, 4, 24, { {1,14},{3,7},{2,3},{0,13},{3,8},{2,3},{0,15}, {2,6},{3,3},{0,13},{2,8},{3,3},{0,16},{2,6},{1,2},{1,15},{2,8},{3,1} });

	// Do pre-processing of the instance
	inst.preprocessing();

	// Find the minimum number of required bins
	string filenameBmin = "InstancesBPPMCF/minNumberOfBinsPerInstance.txt";
	inst.findMinNumberOfBins(filenameBmin);

	// Set the number of bins equal to the minimum number of required bins
	if (inst.Bmin > 0) { inst.B = inst.Bmin; } // change the number of bins

	// Find the L2-lower bounds
	inst.findL2();

	// Print information about the instance
	inst.print();

	// --------------------------------------------------------------
	// Part 2: Bounding procedures

	// Apply BPP-LB
	vector<int> LBs(inst.C);
	Solution solBPPLB = solveBPPLB(inst, LBs);
	solBPPLB.print(true);

	// Apply BPP-UB
	Solution solBPPUB = solveBPPUB(inst);
	solBPPUB.print(true);

	// Apply TS
	int maxNIterNoImprovement = 1000;	// Change the maximum number of iterations without improvement after which the algorithm terminates
	Solution solTS = solveTS(inst, solBPPUB, 50, 40, maxNIterNoImprovement, true, 1000000);
	solTS.print(true);

	// --------------------------------------------------------------
	// Part 3: Methods without warm start

	// Set a time limit
	double timeLimit = 1800;

	// Apply a reimplementation of IP2
	Solution solIP2RE = solveIP2RE(inst, false, timeLimit);
	solIP2RE.print(true);

	// Apply a reimplementation of RM2-GIFF
	Solution solRM2GIFFRE = solveModelAfterOrdering(inst, false, timeLimit, "RM2GIFF");
	solRM2GIFFRE.print(true);

	// Apply LayerFlow
	Solution solLF = solveModelAfterOrdering(inst, false, timeLimit);
	solLF.print(true);

	// Apply HierarchyFlow
	Solution solHF = solveHF(inst, false, timeLimit);
	solHF.print(true);

	// Apply MonoFlow-MultiBin
	Solution solMFMB = solveMFMB(inst, solTS, false, timeLimit - solTS.timeT);
	solMFMB.print(false);

	// --------------------------------------------------------------
	// Part 4: Methods with warm start
	if (not solBPPLB.opt && not solTS.opt) {

		// Update the L2 bounds to the L* bounds
		inst.LBs = LBs;
		inst.LBtot = 0; for (int c = 0; c < inst.C; c++) { inst.LBtot += LBs[c]; }

		// Save the feasible solution found by TS in the instance object
		inst.feasSol = solTS.binPacking;

		// Find the time that was already spent on the bounding procedures
		double preTime = solBPPLB.timeT + solTS.timeT;

		// Apply a reimplementation of IP2(W)
		Solution solIP2REW = solveIP2RE(inst, true, timeLimit - preTime);
		solIP2REW.timeT += preTime;
		solIP2REW.print(true);

		// Apply a reimplementation of RM2-GIFF(W)
		Solution solRM2GIFFREW = solveModelAfterOrdering(inst, true, timeLimit - preTime, "RM2GIFF");
		solRM2GIFFREW.timeT += preTime;
		solRM2GIFFREW.print(true);

		// Apply a reimplementation of RM2-GIFF(W), using the setup of Mehrani et al.
		Solution solRM2GIFFREMehraniW = solveRM2GIFFSetupMehrani(inst, true, timeLimit - preTime);
		solRM2GIFFREMehraniW.timeT += preTime;
		solRM2GIFFREMehraniW.print(true);

		// Apply LayerFlow(W)
		Solution solLFW = solveModelAfterOrdering(inst, true, timeLimit - preTime);
		solLFW.timeT += preTime;
		solLFW.print(true);

		// Apply HierarchyFlow(W)
		Solution solHFW = solveHF(inst, true, timeLimit - preTime);
		solHFW.timeT += preTime;
		solHFW.print(true);

		// Apply MonoFlow-MultiBin(W)
		Solution solMFMBW = solveMFMB(inst, solTS, true, timeLimit - preTime);
		solMFMBW.timeT += preTime;
		solMFMBW.print(true);
	}
}