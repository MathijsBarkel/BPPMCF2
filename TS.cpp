#include "TS.h"

Solution solveTS(const Instance& inst, Solution solBPPUB, int tabuTime, int focusIter, int maxNIterNoImprovement, bool diversification, int maxCFIncrease) {
    // This function applies the algorithm TS to the given instance

    // determine whether progress information is printed or not
    bool printing = false;

    // set the time limit in seconds
    int timeLimit = 1800;  		  
    double start = getCPUTime();  // starting time

    // set a random seed
    srand(53);                      // option 1: every time the same
    // srand((unsigned)time(NULL)); // option 2: every time different

    // first apply BPP-UB
    Solution sol = solBPPUB;
    timeLimit -= sol.timeT;
    sol.method = "TS";

    // initialize variables for the local search
    double TfindBestMove = 0, TimplementMove = 0, TupdateScores = 0;
    double TupdateScoresSuperSwap = 0, TupdateScoresSwap = 0, TupdateScoresSuperTransfer = 0, TupdateScoresTransfer = 0;
    bool lossConcentration = true;

    // if the BPP-UB solution is already optimal  --> return that solution
    // otherwise, try to improve it using local search
    if (sol.feas == 1 && not sol.opt) {
        
        // declare main storage variables
        vector<vector<vector<int>>> curPacking(inst.B);                     // curPacking[b][j]={c,w,idx} gives of the j'th item in the b'th bin its color c, size w, and index idx 
        vector<vector<vector<int>>> bestPacking;                            // same format as curPacking
        vector<int> binLoads(inst.B);                                       // binLoads[b] gives the total overall size of the items packed in bin b
        vector<vector<int>> binLoadsPerColor(inst.C, vector<int>(inst.B));  // binLoadsPerColor[c][b] gives the total overall size of the items of color c packed in bin b
        vector<vector<int>> itemInfo(inst.I);                               // itemInfo[idx]={c,w,b} gives the color c, size w and current bin b of the item with index idx 
        vector<int> emptyBins;                                              // emptyBins is a list of the indices of bins that are empty
        
        // transform the solution to a more workable form
        int idx = 0;
        for (int b = 0; b < sol.binPacking.size(); b++) {                                               // loop over all bins b
            for (int c = 0; c < inst.C; c++) {                                                          // loop over all colors c
                for (int j = 0; j < sol.binPacking[b][c].size(); j++) {                                 // loop over all item types of color c in bin b
                    binLoads[b] += sol.binPacking[b][c][j][1] * sol.binPacking[b][c][j][0];
                    binLoadsPerColor[c][b] += sol.binPacking[b][c][j][1] * sol.binPacking[b][c][j][0];
                    for (int k = 0; k < sol.binPacking[b][c][j][1]; k++) {                              // loop over all copies of the same item type
                        curPacking[b].push_back({ c, sol.binPacking[b][c][j][0], idx });
                        itemInfo[idx] = { c,sol.binPacking[b][c][j][0],b };
                        idx++;
                    }
                }
            }
        }
        for (int b = sol.binPacking.size(); b < inst.B; b++) {  // loop over all remaining (= empty) bins
            emptyBins.push_back(b);
        }

        bestPacking = curPacking;   // the initial best packing is simply the initial packing

        int iter = 0;

        if (printing) {
            cout << "----------------------------------------------------------\n";
            cout << "Start metaheuristic (#CF = " << solBPPUB.UB << ")" << "\n";
            // printSolution(curPacking, itemInfo, binLoads, binLoadsPerColor, iter);
        }

        // Start the local search 
        int iX, cX, wX, bX, lcXbX, lcXbY;
        int iY, cY, wY, bY, lcYbX, lcYbY;
        vector<int> iXBestTransfer, bYBestTransfer;
        vector<int> iXBestSwap, iYBestSwap;
        vector<int> cXBestSuperTransfer, bXBestSuperTransfer, bYBestSuperTransfer;
        vector<int> cXBestSuperSwap, bXBestSuperSwap, cYBestSuperSwap, bYBestSuperSwap;
        double scoreBestTransfer, scoreBestSwap, scoreBestSuperTransfer, scoreBestSuperSwap, time;
        int CF, CFBest, tabuIdx, nbItemsInBin, selectedMove;
        vector<vector<double>> transferScores(inst.I, vector<double>(inst.B, -double('inf')));
        vector<vector<double>> swapScores(inst.I, vector<double>(inst.I, -double('inf')));
        vector<vector<vector<double>>> superTransferScores(inst.C, vector<vector<double>>(inst.B, vector<double>(inst.B, -double('inf'))));
        vector<vector<vector<vector<double>>>> superSwapScores(inst.C, vector<vector<vector<double>>>(inst.B, vector<vector<double>>(inst.C, vector<double>(inst.B, -double('inf')))));
        vector<int> iToUpdate, iXToUpdate, iYToUpdate, bToUpdate;
        vector<bool> iInAffectedB;
        vector<vector<vector<int>>> tabuList(inst.C, vector<vector<int>>(inst.W + 1, vector<int>(inst.B, -1)));
        bool isTabu, tabuResetAllowed = true, CFIncreaseAllowed = false, backToBest = false;
        vector<vector<int>> tabuCandidates;
        int lastImprovement = 0;
        int lastFocusChange = -focusIter - 1;
        bool changeMode = false;
        vector<int> multiBins, multiBinColors, emptyTransfer; int selectedMultiBin, selectedColor;

        // Find the initial amount of color fragmentation
        CF = solBPPUB.UB; CFBest = CF;

        time = getCPUTime();

        // Repeat until the solution is probably optimal, until there are maxNIterNoImprovement iterations without iprovement, or until a time limit is exceeded
        while (CF > inst.LBtot && iter - lastImprovement <= maxNIterNoImprovement && getCPUTime()-start < timeLimit) {

            // change to loss concentration if 
            // - we were in color concentration mode and the best solution hasn't been improved in the past focusIter iterations
            // - or we reset the tabu list while in color concentration mode
            // change to color concentration if
            // - we were in loss concentration mode and the best solution hasn't been improved in the past focusIter iterations
            // - we improved the incumbent while in loss concentration mode
            // - or at the very start of the algorithm (for initialization purposes)
            if (backToBest || iter == 0 || diversification && (iter - lastFocusChange > focusIter || changeMode)) {

                lastFocusChange = iter;


                // There is a feature such that we go back to the best solution if the current solution is more than maxCFIncrease color fragmentations away from the best one
                // However, this feature did not significantly help, so the default value of maxCFIncrease is 1000000 (essentially turning off the feature)
                if (backToBest) {
                    // go back to the best solution found so far
                    curPacking = bestPacking;
                    CF = CFBest;
                    backToBest = false;

                    binLoads.clear(); binLoads.resize(inst.B);
                    binLoadsPerColor.clear(); binLoadsPerColor.resize(inst.C, vector<int>(inst.B));
                    for (int b = 0; b < curPacking.size(); b++) {
                        for (int j = 0; j < curPacking[b].size(); j++) {
                            binLoads[b] += curPacking[b][j][1];
                            binLoadsPerColor[curPacking[b][j][0]][b] += curPacking[b][j][1];
                            itemInfo[curPacking[b][j][2]][2] = b;
                        }
                    }
                    if (printing) { cout << "Solution too far away from the best one --> return to the best solution (#CF = " << CF << ")"; }
                }

                else {
                    lossConcentration = not(lossConcentration);
                    if (printing){
                        if (not changeMode && not iter == 0) { cout << "No improvement in last " << focusIter << " iterations\n"; }
                        cout << "--> change of focus to: ";
                        if (lossConcentration) { cout << "loss concentration"; }
                        else { cout << "color concentration"; }
                    }
                }
                if (printing) { cout << "\n\n"; }
                changeMode = false;
                
                // (Re-)initialize all superSwapScores
                for (int cX = 0; cX < inst.C; cX++) {
                    for (int bX = 0; bX < inst.B - 1; bX++) {
                        if (binLoadsPerColor[cX][bX] > 0) {
                            for (int cY = 0; cY < inst.C; cY++) {
                                for (int bY = bX + 1; bY < inst.B; bY++) {
                                    if (binLoadsPerColor[cY][bY] > 0) {
                                        superSwapScores[cX][bX][cY][bY] = findSuperSwapScore(cX, bX, cY, bY, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                                    }
                                }
                            }
                        }
                    }
                }
                // (Re-)initialize all swapScores
                for (int iX = 0; iX < inst.I - 1; iX++) {
                    for (int iY = iX + 1; iY < inst.I; iY++) {
                        swapScores[iX][iY] = findSwapScore(iX, iY, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                    }
                }
                // (Re-)initialize all superTransferScores
                for (int cX = 0; cX < inst.C; cX++) {
                    for (int bX = 0; bX < inst.B; bX++) {
                        if (binLoadsPerColor[cX][bX] > 0) {
                            for (int bY = 0; bY < inst.B; bY++) {
                                superTransferScores[cX][bX][bY] = findSuperTransferScore(cX, bX, bY, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                            }
                        }
                    }
                }
                // (Re-)initialize all transferScores
                for (int iX = 0; iX < inst.I; iX++) {
                    for (int bY = 0; bY < inst.B; bY++) {
                        transferScores[iX][bY] = findTransferScore(iX, bY, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                    }
                }
            }

            // Evaluate all moves
            emptyTransfer = {};
            scoreBestSuperSwap = -100;
            scoreBestSwap = -100;
            scoreBestSuperTransfer = -100;
            scoreBestTransfer = -100;

            // If there is an empty bin and at least one bin with multiple colors and the focus is CC --> make a random super transfer from a bin with multiple colors to that empty bin
            if (emptyBins.size() > 0 && CF > inst.B - emptyBins.size() && not lossConcentration ) {
                // find the bins with multiple colors
                multiBins = {};
                for (int b = 0; b < inst.B; b++) {
                    if (binLoads[b] > 0) {
                        for (int c = 0; c < inst.C-1; c++) {
                            if (binLoadsPerColor[c][b] > 0 && binLoadsPerColor[c][b] < binLoads[b]) {
                                multiBins.push_back(b);
                                break;
                            }
                        }
                    }
                }
                // randomly select one of the multibins to split up
                selectedMultiBin = multiBins[rand() % multiBins.size()]; 
                multiBinColors = {};
                for (int c = 0; c < inst.C; c++) {
                    if (binLoadsPerColor[c][selectedMultiBin] > 0) {
                        multiBinColors.push_back(c);
                    }
                }
                selectedColor = multiBinColors[rand() % multiBinColors.size()];
                emptyTransfer = { selectedColor, selectedMultiBin, emptyBins.back() };
                emptyBins.pop_back();
                scoreBestSuperTransfer = 0;
            }

            // Otherwise --> do a regular move
            else {
                // First, evaluate all super swap moves
                for (int cX = 0; cX < inst.C; cX++) {
                    for (int bX = 0; bX < inst.B - 1; bX++) {
                        if (binLoadsPerColor[cX][bX] > 0) {
                            for (int cY = 0; cY < inst.C; cY++) {
                                for (int bY = bX + 1; bY < inst.B; bY++) {
                                    if (binLoadsPerColor[cY][bY] > 0) {
                                        if (superSwapScores[cX][bX][cY][bY] >= scoreBestSuperSwap) {
                                            // ignore tabu if the super transfer swap decreases the CF by 1 or 2 and the solution is better than ever before
                                            if ((superSwapScores[cX][bX][cY][bY] > 5 && CF - 2 < CFBest) || (superSwapScores[cX][bX][cY][bY] > 2 && CF - 1 < CFBest)) {
                                                scoreBestSuperSwap = superSwapScores[cX][bX][cY][bY];
                                                cXBestSuperSwap = { cX }; bXBestSuperSwap = { bX }; cYBestSuperSwap = { cY }; bYBestSuperSwap = { bY };
                                            }
                                            // an item can't be moved back to a bin that it left in one of the last tabuTime iterations (unless the CF will be better than ever before)
                                            else {
                                                // determine if the move is tabu
                                                isTabu = false;
                                                for (int j = 0; j < curPacking[bX].size(); j++) {
                                                    iX = curPacking[bX][j][2];
                                                    if (itemInfo[iX][0] == cX && tabuList[itemInfo[iX][0]][itemInfo[iX][1]][bY] >= iter) {
                                                        isTabu = true;
                                                        break;
                                                    }
                                                }
                                                if (not isTabu) {
                                                    for (int j = 0; j < curPacking[bY].size(); j++) {
                                                        iY = curPacking[bY][j][2];
                                                        if (itemInfo[iY][0] == cY && tabuList[itemInfo[iY][0]][itemInfo[iY][1]][bX] >= iter) {
                                                            isTabu = true;
                                                            break;
                                                        }
                                                    }
                                                }

                                                if (not isTabu) {
                                                    if (superSwapScores[cX][bX][cY][bY] > scoreBestSuperSwap) {
                                                        scoreBestSuperSwap = superSwapScores[cX][bX][cY][bY];
                                                        cXBestSuperSwap = { cX }; bXBestSuperSwap = { bX }; cYBestSuperSwap = { cY }; bYBestSuperSwap = { bY };
                                                    }
                                                    else {
                                                        cXBestSuperSwap.push_back(cX); bXBestSuperSwap.push_back(bX); cYBestSuperSwap.push_back(cY); bYBestSuperSwap.push_back(bY);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Second, evaluate all swap moves
                if (scoreBestSuperSwap < 5) {
                    for (int iX = 0; iX < inst.I - 1; iX++) {
                        for (int iY = iX + 1; iY < inst.I; iY++) {
                            if (swapScores[iX][iY] >= scoreBestSwap) {

                                // ignore tabu if the swap decreases the CF by 1 and the solution is better than ever before
                                if (swapScores[iX][iY] > 2 && CF - 1 < CFBest) {
                                    scoreBestSwap = swapScores[iX][iY];
                                    iXBestSwap = { iX }; iYBestSwap = { iY };
                                }
                                // an item can't be moved back to a bin that it left in one of the last tabuTime iterations (unless the CF will be better than ever before)

                                else if (tabuList[itemInfo[iX][0]][itemInfo[iX][1]][itemInfo[iY][2]] < iter && tabuList[itemInfo[iY][0]][itemInfo[iY][1]][itemInfo[iX][2]] < iter) {
                                    if (swapScores[iX][iY] > scoreBestSwap) {
                                        scoreBestSwap = swapScores[iX][iY];
                                        iXBestSwap = { iX }; iYBestSwap = { iY };
                                    }
                                    else {
                                        iXBestSwap.push_back(iX); iYBestSwap.push_back(iY);
                                    }
                                }
                            }
                        }
                    }
                }

                // Third, evaluate all super transfer moves
                if (scoreBestSuperSwap < 4 && scoreBestSwap < 4) {
                    for (int cX = 0; cX < inst.C; cX++) {
                        for (int bX = 0; bX < inst.B; bX++) {
                            if (binLoadsPerColor[cX][bX] > 0) {
                                for (int bY = 0; bY < inst.B; bY++) {
                                    if (superTransferScores[cX][bX][bY] > scoreBestSuperTransfer) {
                                        // ignore tabu if the super transfer decreases the CF by 1 and the solution is better than ever before
                                        if (superTransferScores[cX][bX][bY] > 2 && CF - 1 < CFBest) {
                                            scoreBestSuperTransfer = superTransferScores[cX][bX][bY];
                                            cXBestSuperTransfer = { cX }; bXBestSuperTransfer = { bX };  bYBestSuperTransfer = { bY };
                                        }
                                        // an item can't be moved back to a bin that it left in one of the last tabuTime iterations (unless the CF will be better than ever before)
                                        else {
                                            // determine if the move is tabu
                                            isTabu = false;
                                            for (int j = 0; j < curPacking[bX].size(); j++) {
                                                iX = curPacking[bX][j][2];
                                                if (itemInfo[iX][0] == cX && tabuList[itemInfo[iX][0]][itemInfo[iX][1]][bY] >= iter) {
                                                    isTabu = true;
                                                    break;
                                                }
                                            }

                                            if (not isTabu) {
                                                if (superTransferScores[cX][bX][bY] > scoreBestSuperTransfer) {
                                                    scoreBestSuperTransfer = superTransferScores[cX][bX][bY];
                                                    cXBestSuperTransfer = { cX }; bXBestSuperTransfer = { bX }; bYBestSuperTransfer = { bY };
                                                }
                                                else {
                                                    cXBestSuperTransfer.push_back(cX); bXBestSuperTransfer.push_back(bX); bYBestSuperTransfer.push_back(bY);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Fourth, evaluate all transfer moves
                if (scoreBestSuperSwap < 1 && scoreBestSwap < 1 && scoreBestSuperTransfer < 1) {
                    for (int iX = 0; iX < inst.I; iX++) {
                        for (int bY = 0; bY < inst.B; bY++) {
                            if (transferScores[iX][bY] >= scoreBestTransfer) {
                                // an item can't be moved back to a bin that it left in one of the last tabuTime iterations
                                if (tabuList[itemInfo[iX][0]][itemInfo[iX][1]][bY] < iter) {
                                    if (transferScores[iX][bY] > scoreBestTransfer) {
                                        scoreBestTransfer = transferScores[iX][bY];
                                        iXBestTransfer = { iX }; bYBestTransfer = { bY };
                                    }
                                    else {
                                        iXBestTransfer.push_back(iX); bYBestTransfer.push_back(bY);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            TfindBestMove += getCPUTime() - time;
            time = getCPUTime();

            // If there is no feasible non-tabu move that does not add color fragmentation
            if (emptyTransfer.size() == 0 && scoreBestSuperSwap < -2 && scoreBestSwap < -2 && scoreBestSuperTransfer < -2 && scoreBestTransfer < -2) {
                
                if (scoreBestSuperSwap < -10 && scoreBestSwap < -10 && scoreBestSuperTransfer < -10 && scoreBestTransfer < -10) {
                    if (printing) { cout << "Termination because no more feasible non-tabu moves\n"; }
                    break; // fully terminate
                }
                
                // If allowed, reset the tabuList
                else if (tabuResetAllowed) {
                    for (int c = 0; c < inst.C; c++) {
                        for (int w = 0; w <= inst.W; w++) {
                            for (int b = 0; b < inst.B; b++) {
                                tabuList[c][w][b] = -1;
                            }
                        }
                    }
                    if (printing) { cout << "Reset of tabu list because no more feasible non-CF-increasing non-tabu moves\n"; }
                    tabuResetAllowed = false; // you can't reset again, unless in the mean time the #CF is improved
                    if (not lossConcentration) { changeMode = true; } // switch to focussing on loss concentration
                    else if (printing) { cout << "\n"; }
                    continue; // try again to find the best move, now that the tabus have been lifted and the mode has been changed
                }
                // Otherwise start allowing for CF-increasing moves
                else if (not CFIncreaseAllowed) {
                    CFIncreaseAllowed = true;
                    if (printing) { cout << "No more feasible non-CF-increasing non-tabu moves: allow for CF-increasing moves now\n\n"; }
                }
            }

            // Otherwise, implement the best move 
            tabuCandidates.clear();

            // Case I) The best move is a super swap move
            if (scoreBestSuperSwap >= scoreBestSwap && scoreBestSuperSwap >= scoreBestSuperTransfer && scoreBestSuperSwap >= scoreBestTransfer) {
                selectedMove = rand() % cXBestSuperSwap.size();
                cX = cXBestSuperSwap[selectedMove]; bX = bXBestSuperSwap[selectedMove]; cY = cYBestSuperSwap[selectedMove]; bY = bYBestSuperSwap[selectedMove];
                lcXbX = binLoadsPerColor[cX][bX]; lcXbY = binLoadsPerColor[cX][bY];
                lcYbX = binLoadsPerColor[cY][bX]; lcYbY = binLoadsPerColor[cY][bY];
                nbItemsInBin = curPacking[bX].size();
                iXToUpdate.clear();
                for (int j = nbItemsInBin - 1; j >= 0; j--) {
                    iX = curPacking[bX][j][2];
                    if (itemInfo[iX][0] == cX) {
                        itemInfo[iX][2] = bY;
                        iXToUpdate.push_back(iX);
                        curPacking[bX].erase(curPacking[bX].begin() + j);
                        tabuCandidates.push_back({ iX, bX });
                    }
                }
                nbItemsInBin = curPacking[bY].size();
                iYToUpdate.clear();
                for (int j = nbItemsInBin - 1; j >= 0; j--) {
                    iY = curPacking[bY][j][2];
                    if (itemInfo[iY][0] == cY) {
                        itemInfo[iY][2] = bX;
                        iYToUpdate.push_back(iY);
                        curPacking[bY].erase(curPacking[bY].begin() + j);
                        tabuCandidates.push_back({ iY, bY });
                    }
                }

                for (const int& iX : iXToUpdate) {
                    curPacking[bY].push_back({ itemInfo[iX][0], itemInfo[iX][1], iX });
                }
                for (const int& iY : iYToUpdate) {
                    curPacking[bX].push_back({ itemInfo[iY][0], itemInfo[iY][1], iY });
                }

                binLoads[bX] += lcYbY - lcXbX;
                binLoads[bY] += lcXbX - lcYbY;
                binLoadsPerColor[cX][bX] -= lcXbX; binLoadsPerColor[cX][bY] += lcXbX;
                binLoadsPerColor[cY][bX] += lcYbY; binLoadsPerColor[cY][bY] -= lcYbY;

                // randomly choose one of the moved items for the tabu
                tabuIdx = rand() % tabuCandidates.size();
                tabuList[itemInfo[tabuCandidates[tabuIdx][0]][0]][itemInfo[tabuCandidates[tabuIdx][0]][1]][tabuCandidates[tabuIdx][1]] = iter + tabuTime;

                if (scoreBestSuperSwap > 5) {
                    CF -= 2;
                    if (CF < CFBest) {
                        CFBest = CF;
                        bestPacking = curPacking;
                        tabuResetAllowed = true;
                        lastImprovement = iter; lastFocusChange = iter;
                        if (lossConcentration) { changeMode = true; }
                    }
                }
                else if (scoreBestSuperSwap > 2) {
                    CF -= 1;
                    if (CF < CFBest) {
                        CFBest = CF;
                        bestPacking = curPacking;
                        tabuResetAllowed = true;
                        lastImprovement = iter; lastFocusChange = iter;
                        if (lossConcentration) { changeMode = true; }
                    }
                }
                if (printing) {
                    cout << "iter " << iter + 1 << " - swap all items of color " << cX << " from bin " << bX << " [";
                    for (int j = 0; j < iXToUpdate.size(); j++) {
                        cout << iXToUpdate[j] << "(" << itemInfo[iXToUpdate[j]][1] << ")";
                        if (j != iXToUpdate.size() - 1) { cout << ", "; }
                    }
                    cout << "] with all items of color " << cY << " from bin " << bY << " [";
                    for (int j = 0; j < iYToUpdate.size(); j++) {
                        cout << iYToUpdate[j] << "(" << itemInfo[iYToUpdate[j]][1] << ")";
                        if (j != iYToUpdate.size() - 1) { cout << ", "; }
                    }
                    cout << "] (score: " << scoreBestSuperSwap << ", #CF = " << CF << ")" << "\n";
                }
            }

            // Case I) The best move is a swap move
            else if (scoreBestSwap >= scoreBestSuperSwap && scoreBestSwap >= scoreBestSuperTransfer && scoreBestSwap >= scoreBestTransfer) {
                selectedMove = rand() % iXBestSwap.size();
                iX = iXBestSwap[selectedMove], cX = itemInfo[iX][0]; wX = itemInfo[iX][1]; bX = itemInfo[iX][2];
                iY = iYBestSwap[selectedMove], cY = itemInfo[iY][0]; wY = itemInfo[iY][1]; bY = itemInfo[iY][2];
                itemInfo[iX][2] = bY;
                itemInfo[iY][2] = bX;
                binLoads[bX] -= wX; binLoads[bX] += wY;
                binLoads[bY] += wX; binLoads[bY] -= wY;
                binLoadsPerColor[cX][bX] -= wX; binLoadsPerColor[cX][bY] += wX;
                binLoadsPerColor[cY][bX] += wY; binLoadsPerColor[cY][bY] -= wY;
                curPacking[bY].push_back({ cX, wX, iX });
                for (int j = 0; j < curPacking[bX].size(); j++) {
                    if (curPacking[bX][j][2] == iX) {
                        curPacking[bX].erase(curPacking[bX].begin() + j);
                        break;
                    }
                }
                curPacking[bX].push_back({ cY, wY, iY });
                for (int j = 0; j < curPacking[bY].size(); j++) {
                    if (curPacking[bY][j][2] == iY) {
                        curPacking[bY].erase(curPacking[bY].begin() + j);
                        break;
                    }
                }
                if (scoreBestSwap > 2) {
                    CF -= 1;
                    if (CF < CFBest) {
                        CFBest = CF;
                        bestPacking = curPacking;
                        tabuResetAllowed = true;
                        lastImprovement = iter; lastFocusChange = iter;
                        if (lossConcentration) { changeMode = true; }
                    }
                }
                else if (scoreBestSwap < -5) {
                    CF += 2;
                    if (CF > CFBest + maxCFIncrease) {
                        backToBest = true;
                    }
                }
                else if (scoreBestSwap < -2) {
                    CF += 1;
                    if (CF > CFBest + maxCFIncrease) {
                        backToBest = true;
                    }
                }

                if (printing) { cout << "iter " << iter + 1 << " - swap item " << iX << " = (" << cX << ", " << wX << ")" << " from bin " << bX << " with item " << iY << " = (" << cY << ", " << wY << ")" << " from bin " << bY << " (score: " << scoreBestSwap << ", #CF = " << CF << ")" << "\n"; }

                // Update the tabu list
                tabuList[itemInfo[iX][0]][itemInfo[iX][1]][bX] = iter + tabuTime;
                tabuList[itemInfo[iY][0]][itemInfo[iY][1]][bY] = iter + tabuTime;
            }

            // Case III) The best move is a super transfer move
            else if (scoreBestSuperTransfer >= scoreBestSuperSwap && scoreBestSuperTransfer >= scoreBestSwap && scoreBestSuperTransfer >= scoreBestTransfer) {

                if (emptyTransfer.size() > 0) {
                    cX = emptyTransfer[0]; bX = emptyTransfer[1]; bY = emptyTransfer[2];
                }
                else {
                    selectedMove = rand() % cXBestSuperTransfer.size();
                    cX = cXBestSuperTransfer[selectedMove]; bX = bXBestSuperTransfer[selectedMove]; bY = bYBestSuperTransfer[selectedMove];
                }
                lcXbX = binLoadsPerColor[cX][bX];

                nbItemsInBin = curPacking[bX].size();
                iToUpdate.clear();
                for (int j = nbItemsInBin - 1; j >= 0; j--) {
                    iX = curPacking[bX][j][2];
                    if (itemInfo[iX][0] == cX) {
                        itemInfo[iX][2] = bY;
                        iToUpdate.push_back(iX);
                        curPacking[bY].push_back({ cX, curPacking[bX][j][1], iX });
                        curPacking[bX].erase(curPacking[bX].begin() + j);
                        tabuCandidates.push_back({ iX, bX });
                    }
                }

                binLoads[bX] -= lcXbX; binLoads[bY] += lcXbX;
                binLoadsPerColor[cX][bX] -= lcXbX; binLoadsPerColor[cX][bY] += lcXbX;

                if (binLoads[bX] == 0) {
                    emptyBins.push_back(bX);
                }

                // randomly choose one of the moved items for the tabu
                tabuIdx = rand() % tabuCandidates.size();
                tabuList[itemInfo[tabuCandidates[tabuIdx][0]][0]][itemInfo[tabuCandidates[tabuIdx][0]][1]][tabuCandidates[tabuIdx][1]] = iter + tabuTime;

                if (scoreBestSuperTransfer > 2) {
                    CF--;
                    if (CF < CFBest) {
                        CFBest = CF;
                        bestPacking = curPacking;
                        tabuResetAllowed = true;
                        lastImprovement = iter; lastFocusChange = iter;
                        if (lossConcentration) { changeMode = true; }
                    }
                }
                if (printing) {
                    cout << "iter " << iter + 1 << " - move all items of color " << cX << " from bin " << bX << " [";
                    for (int j = 0; j < iToUpdate.size(); j++) {
                        cout << iToUpdate[j] << "(" << itemInfo[iToUpdate[j]][1] << ")";
                        if (j != iToUpdate.size() - 1) { cout << ", "; }
                    }
                    cout << "] to ";
                    if (emptyTransfer.size() > 0) { cout << "empty(!) "; }
                    cout << "bin " << bY << " (score: " << scoreBestSuperTransfer << ", #CF = " << CF << ")" << "\n";
                }
            }

            // Case IV) The best move is a transfer move
            else {
                selectedMove = rand() % iXBestTransfer.size();
                iX = iXBestTransfer[selectedMove], cX = itemInfo[iX][0]; wX = itemInfo[iX][1]; bX = itemInfo[iX][2];
                bY = bYBestTransfer[selectedMove];
                itemInfo[iX][2] = bY;
                binLoads[bX] -= wX; binLoads[bY] += wX;
                binLoadsPerColor[cX][bX] -= wX; binLoadsPerColor[cX][bY] += wX;
                curPacking[bY].push_back({ cX, wX, iX });
                for (int j = 0; j < curPacking[bX].size(); j++) {
                    if (curPacking[bX][j][2] == iX) {
                        curPacking[bX].erase(curPacking[bX].begin() + j);
                        break;
                    }
                }

                if (scoreBestTransfer < -2) {
                    CF++;
                    if (CF > CFBest + maxCFIncrease) {
                        backToBest = true;
                    }
                }
                if (printing) { cout << "iter " << iter + 1 << " - move item " << iX << " = (" << cX << ", " << wX << ")" << " from bin " << bX << " to bin " << bY << " (score: " << scoreBestTransfer << ", #CF = " << CF << ")" << "\n"; }

                // Update the tabu list
                tabuList[itemInfo[iX][0]][itemInfo[iX][1]][bX] = iter + tabuTime;
            }

            TimplementMove += getCPUTime() - time;
            time = getCPUTime();

            // Print the new solution
            iter++;
            // printSolution(curPacking, itemInfo, binLoads, binLoadsPerColor, iter);

            // Update the transferScores, swapScores and superTransferScores
            if (not diversification || not changeMode) {
                bToUpdate = { bX, bY };
                for (int iX = 0; iX < inst.I; iX++) {
                    for (const int& bY : bToUpdate) {
                        transferScores[iX][bY] = findTransferScore(iX, bY, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                    }
                }

                TupdateScoresTransfer += getCPUTime() - time;

                iToUpdate.clear();
                iInAffectedB.clear(); iInAffectedB.resize(inst.I, false);
                for (const int& b : bToUpdate) {
                    for (int j = 0; j < curPacking[b].size(); j++) {
                        iToUpdate.push_back(curPacking[b][j][2]);
                        iInAffectedB[curPacking[b][j][2]] = true;
                    }
                }

                time = getCPUTime();

                for (const int& iX : iToUpdate) {
                    for (int bY = 0; bY < inst.B; bY++) {
                        transferScores[iX][bY] = findTransferScore(iX, bY, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                    }
                }

                TupdateScoresTransfer += getCPUTime() - time;
                time = getCPUTime();

                for (int iX = 0; iX < inst.I - 1; iX++) {
                    for (int iY = iX + 1; iY < inst.I; iY++) {
                        if (iInAffectedB[iX] || iInAffectedB[iY]) {
                            swapScores[iX][iY] = findSwapScore(iX, iY, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                        }
                    }
                }

                TupdateScoresSwap += getCPUTime() - time;
                time = getCPUTime();

                for (int bToRemove = 0; bToRemove < inst.B; bToRemove++) {
                    for (int bToAdd = 0; bToAdd < inst.B; bToAdd++) {
                        if (bToRemove != bToAdd && (bToRemove == bX || bToRemove == bY || bToAdd == bX || bToAdd == bY)) {
                            for (int cX = 0; cX < inst.C; cX++) {
                                superTransferScores[cX][bToRemove][bToAdd] = findSuperTransferScore(cX, bToRemove, bToAdd, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                            }
                        }
                    }
                }

                TupdateScoresSuperTransfer += getCPUTime() - time;
                time = getCPUTime();

                for (int b1 = 0; b1 < inst.B - 1; b1++) {
                    for (int b2 = b1 + 1; b2 < inst.B; b2++) {
                        if (b1 != b2 && (b1 == bX || b1 == bY || b2 == bX || b2 == bY)) {
                            for (int cX = 0; cX < inst.C; cX++) {
                                for (int cY = 0; cY < inst.C; cY++) {
                                    superSwapScores[cX][b1][cY][b2] = findSuperSwapScore(cX, b1, cY, b2, lossConcentration, inst.W, itemInfo, binLoads, binLoadsPerColor);
                                }
                            }
                        }
                    }
                }
                TupdateScoresSuperSwap += getCPUTime() - time;
                time = getCPUTime();
            }
        }
        // printSolution(curPacking, itemInfo, binLoads, binLoadsPerColor, iter);

        if (getCPUTime() >= timeLimit) {
            if (printing) { cout << "Termination because the time limit has been reached\n"; }
        }
        else if (iter - lastImprovement > maxNIterNoImprovement) {
            if (printing) { cout << "Termination because the maximum number of iterations without improvement (" << maxNIterNoImprovement << ") is reached\n"; }
        }

        // return the solution
        sol.feas = 1;
        sol.opt = 0;
        int nbColorFragmentations = 0;
        vector<vector<int>> curBin;
        vector<vector<vector<int>>> curCountedBin;
        vector<vector<vector<vector<int>>>> binPacking;
        for (int b = 0; b < inst.B; b++) {
            if (bestPacking[b].size() > 0) {
                curBin.clear(); curBin.resize(inst.C);
                for (int j = 0; j < bestPacking[b].size(); j++) {
                    curBin[bestPacking[b][j][0]].push_back(bestPacking[b][j][1]);
                }
                curCountedBin.clear(); curCountedBin.resize(inst.C);
                for (int c = 0; c < inst.C; c++) {
                    if (curBin[c].size() > 0) {
                        curCountedBin[c] = sortAndCount(curBin[c]);
                        nbColorFragmentations++;
                    }
                }
                binPacking.push_back(curCountedBin);
            }
        }
        //while (binPacking.size() != inst.B) {
        //    curCountedBin.clear(); curCountedBin.resize(inst.C);
        //    binPacking.push_back(curCountedBin);
        //}
        sol.feas = 1;
        sol.UB = nbColorFragmentations;
        if (sol.UB == sol.LB) {
            sol.opt = 1;
        }
        sol.binPacking = binPacking;
        // We save the number of iterations in the Nvar attribute
        sol.Nvar = iter, sol.Nconstr = 0, sol.Ncoeff = 0;
    }

    sol.timeT += getCPUTime() - start; // save the total time

    TupdateScores = TupdateScoresSuperSwap + TupdateScoresSwap + TupdateScoresSuperTransfer + TupdateScoresTransfer;
    //if (printing) {
    //    cout << "\n";
    //    cout << "Total time spent: " << sol.timeT << "s\n";
    //    cout << "- Time spent by BPP-UB: " << solBPPUB.timeT << "s\n";
    //    cout << "- Time spent finding the best move: " << TfindBestMove << "s\n";
    //    cout << "- Time spent implementing the best move: " << TimplementMove << "s\n";
    //    cout << "- Time spent updating the scores: " << TupdateScores << "s\n";
    //    cout << "--> Super swaps: " << TupdateScoresSuperSwap << "s\n";
    //    cout << "--> Single swaps: " << TupdateScoresSwap << "s\n";
    //    cout << "--> Super transfers: " << TupdateScoresSuperTransfer << "s\n";
    //    cout << "--> Single transfers: " << TupdateScoresTransfer << "s\n";
    //}

    return sol;
}

void printSolution(const vector<vector<vector<int>>>& curPacking, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor, int iter) {
    // This function is used to print the current solution from TS
    cout << "---------------------------------------------------------------\n";
    cout << "Solution after iter " << iter << ":\n";
    print3DVector(curPacking, "curPacking", "bin");
    print2DVector(itemInfo, "itemInfo [c,w,b]", "item");
    printVector(binLoads, "binLoads");
    print2DVector(binLoadsPerColor, "binLoadsPerColor", "color");
    cout << "---------------------------------------------------------------\n";
}

double findSuperSwapScore(int cX, int bX, int cY, int bY, bool lossConcentration, int W, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor) {
    /* Score ranges in mode CC / mode LC / overall:
    CF decreases by 2: (6,8]   / (5,7)   / (5, )
    CF decreases by 1: (3,4]   / (2,4)   / (2,5)
    CF unchanged:      0       / (-1,1)  / (-2,2)
    CF increases by 1: -       / -       / (-5,-2) 
    CF decreases by 2: -       / -       / ( ,-5)
    */
    
    int lX, lcXbX, lcYbX, lY, lcXbY, lcYbY;
    double score;
    lX = binLoads[bX]; lcXbX = binLoadsPerColor[cX][bX]; lcXbY = binLoadsPerColor[cX][bY];
    lY = binLoads[bY]; lcYbX = binLoadsPerColor[cY][bX]; lcYbY = binLoadsPerColor[cY][bY];

    // if the swap is feasible, give a score. Otherwise give score -inf
    if (bY != bX && lcXbX > 0 && lcYbY > 0 && lY + lcXbX - lcYbY <= W && lX - lcXbX + lcYbY <= W && (lX != lcXbX || lY != lcYbY)) {
        score = 0;
        if (cX != cY) {
            if (lcXbY > 0) {
                score += 3;
                if (not lossConcentration) {
                    score += (lcXbX + lcXbY - abs(lcXbX - lcXbY)) / double(W);
                }
            }
            if (lcYbX > 0) {
                score += 3;
                if (not lossConcentration) {
                    score += (lcYbX + lcYbY - abs(lcYbX - lcYbY)) / double(W);
                }
            }
        }
        if (lossConcentration) {
            score += (abs(lX - lY - 2 * lcXbX + 2 * lcYbY) - abs(lX - lY)) / double(W);
        }
    }
    else { score = -double('inf'); }

    // TEMPORARILY: Don't allow moves with score 0
    if (score == 0) { score = -1000; }

    return score;
}

double findSwapScore(int iX, int iY, bool lossConcentration, int W, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor) {
    /* Score ranges in mode CC / mode LC / overall:
    CF decreases by 2: -       / -       / -      <-- this would be a super swap
    CF decreases by 1: (2,5)   / (2,4)   / (2,5)
    CF unchanged:      (-2,2)  / (-1,1)  / (-2,2)
    CF increases by 1: (-5,-2) / (-4,-2) / (-5,-2)
    CF decreases by 2: [-8,6)  / (-7,-5) / ( ,-5)
    */

    // if the swap is feasible, give a score. Otherwise give score -inf
    int bX = itemInfo[iX][2]; int bY = itemInfo[iY][2];
    if (bX == bY) { return -double('inf'); }
    int wX = itemInfo[iX][1]; int lX = binLoads[bX];
    int wY = itemInfo[iY][1]; int lY = binLoads[bY];
    if (lY + wX - wY > W || lX - wX + wY > W) { return -double('inf'); }
    int cX = itemInfo[iX][0]; int cY = itemInfo[iY][0];
    if (wX == wY && cX == cY) { return -double('inf'); }
    int lcXbX = binLoadsPerColor[cX][bX]; int lcYbY = binLoadsPerColor[cY][bY];
    if (wX == lcXbX && wY == lcYbY) { return -double('inf'); }
    int lcXbY = binLoadsPerColor[cX][bY];  int lcYbX = binLoadsPerColor[cY][bX];

    int dCF = 0; double score = 0;

    // if the colors are different
    if (cX != cY) {
        if (lcXbX == wX) { dCF--; }
        if (lcXbY == 0) { dCF++; }
        if (lcYbY == wY) { dCF--; }
        if (lcYbX == 0) { dCF++; }
        score += -3 * dCF;
        if (lossConcentration) {
            score += (abs(lX - lY - 2 * wX + 2 * wY) - abs(lX - lY)) / double(W);
        }
        else {
            score += (abs(lcXbX - lcXbY - 2 * wX) - abs(lcXbX - lcXbY)) / double(W);
            score += (abs(lcYbX - lcYbY + 2 * wY) - abs(lcYbX - lcYbY)) / double(W);
        }
    }
    // if the colors are the same
    else {
        if (lossConcentration) {
            score += (abs(lX - lY - 2 * wX + 2 * wY) - abs(lX - lY)) / double(W);
        }
        else {
            score += (abs(lcXbX - lcXbY - 2 * wX + 2 * wY) - abs(lcXbX - lcXbY)) / double(W);
        }
    }
    // TEMPORARILY: Don't allow moves with score 0
    // if (score == 0) { score = -1000; }

    return score;
}

double findSuperTransferScore(int cX, int bX, int bY, bool lossConcentration, int W, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor) {
    /* Score ranges in mode CC / mode LC / overall:
    CF decreases by 2: -       / -       / (5, ) 
    CF decreases by 1: (3,4]   / (3,4]   / (2,5)
    CF unchanged:      0       / [-1,1]  / (-2,2)
    CF increases by 1: -       / -       / (-5,-2)
    CF decreases by 2: -       / -       / ( ,-5)
    */

    int lX, lY, lcXbX, lcXbY;
    double score = 0;
    lX = binLoads[bX]; lY = binLoads[bY]; lcXbX = binLoadsPerColor[cX][bX];

    // if the super transfer move is feasible, give a score. Otherwise, the score is -inf
    if (bY != bX && lcXbX > 0 && lY + lcXbX <= W && (lY != 0 || lX != lcXbX)) {
        lcXbY = binLoadsPerColor[cX][bY];
        score = 0;
        if (lcXbY > 0) { 
            score += 3; 
            if (not lossConcentration) { // if the focus is on color concentration
                score += (lcXbX + lcXbY - abs(lcXbX - lcXbY)) / double(W);
            }
        }
        if (lossConcentration) { // if the focus is on loss concentration
            score += (abs(lX - lY - 2 * lcXbX) - abs(lX - lY)) / double(W);
        }
    }
    else { score = -double('inf'); }

    // TEMPORARILY: Don't allow moves with score 0
    // if (score == 0) { score = -1000; }


    return score;
}

double findTransferScore(int iX, int bY, bool lossConcentration, int W, const vector<vector<int>>& itemInfo, const vector<int>& binLoads, const vector<vector<int>>& binLoadsPerColor) {
    /* Score ranges in mode CC / mode LC / overall:
    CF decreases by 2: -       / -       / (5, ) 
    CF decreases by 1: -       / -       / (2,5)   <-- this would be a super transfer
    CF unchanged:      (-1,1)  / [-1,1]  / (-2,2)
    CF increases by 1: (-4,-2) / (-4,-2) / (-5,-2)
    CF decreases by 2: -       / -       / ( ,-5)
    */
    
    int cX, wX, bX, lX, lcXbX, lY, lcXbY;
    double score = 0;
    wX = itemInfo[iX][1]; bX = itemInfo[iX][2]; lY = binLoads[bY];
    cX = itemInfo[iX][0]; lcXbX = binLoadsPerColor[cX][bX];

    // if the transfer move is feasible, give a score. Otherwise, the score is -inf
    if (bY != bX && wX != lcXbX && lY + wX <= W) {
        lcXbY = binLoadsPerColor[cX][bY];
        if (lcXbY == 0) { score -= 3; }
        if (lossConcentration) { // if we focus on loss concentration
            lX = binLoads[bX];
            score += (abs(lX - lY - 2 * wX) - abs(lX - lY)) / double(W);
        }
        else { // if we focus on color concentration
            score += (abs(lcXbX - lcXbY - 2 * wX) - abs(lcXbX - lcXbY)) / double(W);
        }
    }
    else { score = -double('inf'); }

    // TEMPORARILY: Don't allow moves with score 0
    // if (score == 0) { score = -1000; }

    return score;
}
