#include "helper_functions.h"

void Instance::print() {
    // This function prints information about the Instance struct
    cout << "\n------------------------------------------------------\n";
    cout << "Instance information (\"" << name << "\"):\n";
    cout << "- number of bins: B = " << B << '\n';
    cout << "- bin capacity: W = " << W << '\n';
    cout << "- number of colors: C = " << C << '\n';
    cout << "- number of items: I = " << I << '\n';
    cout << "- number of item types: J = " << J << '\n';
    cout << "- mimimum number of required bins: Bmin = " << Bmin << '\n';
    cout << "- L2-LB: " << LBtot;
    printVector(LBs, "");
    cout << "- sorted item sizes and demands per color: \n";
    for (int c = 0; c < C; c++) {
        cout << " --> " << "color " << c << ": ";
        if (groupedItems[c].size() == 0) {
            cout << "-\n";
            continue;
        }
        cout << "[";
        for (size_t i = 0; i < groupedItems[c].size(); i++) {
            cout << groupedItems[c][i][0] << "(" << groupedItems[c][i][1] << ")";
            if (i != groupedItems[c].size() - 1) cout << ", ";
        }
        cout << "]\n";

    }
    cout << "\n";
}

void Instance::preprocessing(bool reset, bool reverse) {
    // this function adds convenient attributes to the instance object
    // for every attribute items are sorted first acc. to color (lowest index first), then acc. to size (highest size first)
    // - items: list containing sorted uncounted tuples (c,w), possibly containing duplicates
    // - countedItems: list containing sorted counted tuples (c,w,q)
    // - groupedItems: nested list, each 'layer' corresponds to a color and contains a list of sorted counted tuples (w,q)
    // - colorlessItems: list containing sorted counted colorless tuples (w,q) (item colors are ignored)
    // also computes J: the number of unique item types
    // and scales the bin capacity and item sizes by 2 if the bin capacity is odd

    // ensure that the bin capacity is even (needed for AF2RR)
    //if (W % 2 != 0) {
    //    W *= 2;
    //    for (int i = 0; i < items.size(); i++) {
    //        items[i][1] *= 2;
    //    }
    //}

    // sort the list of all items
    if (reverse) { sort(items.begin(), items.end(), reverseItemOrdering); }
    else { sort(items.begin(), items.end(), itemOrdering); }

    // reset the attributes if the function is called for a second time
    if (reset) {
        countedItems.clear();
        groupedItems.clear();
    }

    // initialize
    int c, w, cPrev = 0, wPrev = items[0][1], q = 0;
    vector<int> colorlessItemsUnsorted;
    groupedItems.resize(C);
    J = 1;

    // set the attributes
    for (int i = 0; i < I; i++) { // loop over all uncounted sorted items
        c = items[i][0]; // color of current item
        w = items[i][1]; // size of current item
        colorlessItemsUnsorted.push_back(w);

        // if we switch colors or size
        if (w != wPrev || c != cPrev) {
            countedItems.push_back({ cPrev,wPrev,q });
            groupedItems[cPrev].push_back({ wPrev, q });
            cPrev = c, wPrev = w, q = 1, J++;
        }
        // if the item is the same
        else {q += 1;}
    }
    // finish processing the last item
    countedItems.push_back({ cPrev,wPrev,q });
    groupedItems[cPrev].push_back({ wPrev, q });

    // sort the list of colorless items
    colorlessItems = sortAndCount(colorlessItemsUnsorted);
 }

void Instance::findL2() {
    // This function finds the L2 lower bound by Martello and Toth for each color-subinstance
    // it's stored in the attribute LBs (and the sum of them all is stored in LBtot)

    // declare variables
    LBs.resize(C);
    LBtot = 0;
    vector<vector<int>> colorSubset;
    int S, F2, F3, K, LK, L1, L2, n1, n2, n3, S2, S3, changeF2, constantn1pn2;

    for (int c = 0; c < C; c++) { // we derive an LB for each color seperately

        // reinitialize the variables for the new color
        colorSubset = groupedItems[c];
        S = 0, n1 = 0, n2 = 0, n3 = 0, S2 = 0, S3 = 0, L2 = 0, F2 = -1, F3 = colorSubset.size();

        // first calculate L1 (which is L(K) for K=0, so subpart of L2)
        for (int j = 0; j < colorSubset.size(); j++) {
            S += colorSubset[j][1] * colorSubset[j][0];
        }
        L1 = ceil(S / double(W));

        // To compute L2, first determine max j such that w_j <= W/2.
        for (int j = 0; j < colorSubset.size(); j++) {
            if (colorSubset[j][0] <= W / 2) {
                F3 = j; // index of first item in N3
                break;
            }
        }

        // trivial case: only small items
        if (F3 == 0) {
            for (int j = 0; j < colorSubset.size(); j++) {
                S3 += colorSubset[j][1] * colorSubset[j][0];
            }
            LBs[c] = ceil(S3 / double(W));
        }

        // trivial case: only large items
        else if (F3 == colorSubset.size()) {
            for (int j = 0; j < colorSubset.size(); j++) {
                n1 += colorSubset[j][1];
            }
            LBs[c] = n1;
        }

        // non-trivial case:
        else {
            // first compute L(w_F3)
            K = colorSubset[F3][0];
            n3 = colorSubset[F3][1];
            S3 = colorSubset[F3][1] * colorSubset[F3][0];
            for (int j = F3 - 1; j >= 0; j--) {
                if (colorSubset[j][0] <= W - K) {
                    n2 += colorSubset[j][1];
                    S2 += colorSubset[j][1] * colorSubset[j][0];
                }
                else {
                    F2 = j + 1; // index of first item in N2
                    break;
                }
            }
            for (int j = F2 - 1; j >= 0; j--) {
                n1 += colorSubset[j][1];
                // S1 += colorSubset[j][1] * colorSubset[j][0];
            }

            constantn1pn2 = n1 + n2;
            L2 = constantn1pn2 + max(double(0), ceil(S3 / double(W) - n2 + S2 / double(W)));
            // now compute also L(K) for K = w_{F3+1}, ..., w_{end}
            for (int j = F3 + 1; j < colorSubset.size(); j++) {
                K = colorSubset[j][0];
                n3 += colorSubset[j][1];
                S3 += colorSubset[j][1] * colorSubset[j][0];
                changeF2 = 0;
                if (n2 < constantn1pn2) {
                    for (int i = F2 - 1; i >= 0; i--) {
                        if (colorSubset[i][0] <= W - K) {
                            changeF2++;
                            n2 += colorSubset[i][1];
                            S2 += colorSubset[i][1] * colorSubset[i][0];
                            // n1 -= colorSubset[i][1];
                            // S1 -= colorSubset[i][1] * colorSubset[i][0];
                        }
                        else {
                            F2 -= changeF2;
                            break;
                        }
                    }
                }
                LK = constantn1pn2 + max(double(0), ceil(S3 / double(W) - n2 + S2 / double(W)));
                L2 = max(L2, LK);
            }

            // save the final bound
            LBs[c] = max(L1, L2);
        }
        // update the sum of the L2 lower bounds
        LBtot += LBs[c];
    }
}

void Instance::findMinNumberOfBins(string filename) {
    // open the file
    ifstream file(filename);

    // read the file
    if (file.is_open()) { //if the file is open
        string line;
        bool first = true;
        while (true) {
            getline(file, line, ';');
            if (first) {
                line = line.substr(3);
                first = false;
            }
            if (line == name) {
                getline(file, line, '\n');
                Bmin = stoi(line);
                break;
            }
            else if (line == "") {
                Bmin = 0;
                break;
            }
            else {
                getline(file, line, '\n');
            }
        }

        // close the file
        file.close();
    }
}

vector<vector<int>> Instance::findColorSubset(const vector<bool> &subset) {
    // This function selects all items of the given colors according to size (high to low)
    vector<int> allItems;
    for (int i = 0; i < I; i++) {
        if (subset[items[i][0]]) { allItems.push_back(items[i][1]); }
    }
    return sortAndCount(allItems);
}

void Instance::draw() {
    // This function prints code that can be copied to LaTeX for drawing an instance
    cout << "\\begin{figure}[h]\n";
    cout << "\\centering\n";
    cout << "\\caption{Place a nice caption here}\n";
    double xscale, yscale;
    xscale = min(0.5, double(18)/(1.5*I + 0.5));
    yscale = min(0.5, double(18)/W);
    cout << "\\begin{tikzpicture}[xscale = " << xscale << ", yscale = "<< yscale << "]\n";

    // add horizontal grid lines
    for (int v = 0; v <= W; v++) {
        cout << "\\draw[gray, thin](-0.5, " <<v << ")--(" << 1.5 * I << ", " << v << "); \n";
    }
    int c, w, q, val, i = -1;
    for (int j = 0; j < countedItems.size(); j++) {     // for every item type
        c = countedItems[j][0];                         // find the item color
        w = countedItems[j][1];                         // find the item size
        q = countedItems[j][2];                         // find the item quantity
        for (int k = 0; k < q; k++) {                   // for every duplicate
            i++;                                        // add 1 to the number of items
            val = (c + 1) / double(C + 1) * 100;        // determine the color: lower c = lighter color
            // draw the item              
            cout << "\\draw[thick, fill = black!" << val << "!white](" << 1.5*i << ", " << 0 << ")--(" << 1.5*i + 1 << ", " << 0 << ")--(" << 1.5*i + 1 << ", " << w << ")--(" << 1.5*i << ", " << w << ")--(" << 1.5*i << ", " << 0 << "); \n";
        }
    }

    cout << "\\end{tikzpicture}\n";
    cout << "\\label{Place a nice label here}\n";
    cout << "\\end{figure}\n\n";
}

Instance readInstance(string filename) {
    // This function reads an instance from a .txt file and saves it as an Instance struct

    // Note: I also use getline to skip lines, maybe there is a more efficient way to do so. 
    // In particular, I can skip the 0-matrix in the middle of each file.

    // Instances can either be standard, or a triplet instance":
    /* Structure of each standard:
    1
    number of bins (B)
    bin capacity (W)
    1 blank line
    0-matrix of size BxW
    3 blank lines
    number of colors (C)
    number of items (I)
    1 blank line
    matrix of size Ix2: - first column: color (c), second column: size (w)
    */

    /* Structure of each triplet file:
    number of items (I)
    c1, w1
    c2, w2
    ...
    cI, wI
    */

    // initialize Instance object
    Instance inst;
    inst.name = filename;

    // open the file
    ifstream file(filename);

    // read the file
    if (file.is_open()) { //if the file is open
        string line;
        getline(file, line, '\n');
        inst.I = stoi(line);
        if (inst.I == 1) { // If the first element is a 1, then it's a standard instance
            getline(file, line, '\n');
            inst.B = stoi(line);
            getline(file, line, '\n');
            inst.W = stoi(line);
            for (int i = 0; i < inst.B + 4; i++) {
                getline(file, line, '\n');
            }
            getline(file, line, '\n');
            inst.C = stoi(line);
            getline(file, line, '\n');
            inst.I = stoi(line);
            getline(file, line, '\n');
            inst.items.resize(inst.I, vector<int>(2, 0));
            for (int i = 0; i < inst.I; i++) {
                getline(file, line, '\t');
                inst.items[i][0] = stoi(line);
                getline(file, line, '\n');
                inst.items[i][1] = stoi(line);
            }
        }
        else { // If the first element is not 1, then it's a triplet instance
            inst.W = 1000;
            inst.C = 3;
            inst.B = inst.I / 3;
            inst.items.resize(inst.I, vector<int>(2, 0));
            for (int i = 0; i < inst.I; i++) {
                getline(file, line, ',');
                inst.items[i][0] = stoi(line);
                getline(file, line, '\n');
                inst.items[i][1] = stoi(line);
            }
        }
        // close the file
        file.close();
    }
    else {
        // if the file cannot be opened: print error and return default Instance
        cout << "Unable to open file";
    }
    return inst;
}

Instance createInstance(string name, int I, int B, int C, int W, const vector<vector<int>>& items) {
    // This function creates an Instance object with the specified attributes
    Instance inst;
    inst.name = name, inst.I = I, inst.B = B, inst.C = C, inst.W = W, inst.items = items;
    return inst;
}

Instance generateInstance(string name, int I, int B, int C, int W, int r) {
    srand(r);
    Instance inst;
    inst.name = name, inst.I = I, inst.B = B, inst.C = C, inst.W = W;
    inst.items = {};
    int c, w;
    for (int i = 0; i < inst.I; i++) {
        if (i < inst.C) {
            c = i;
        }
        else {
            c = rand() % C;
        }
        w = rand() % (W-1) + 1;
        inst.items.push_back({ c,w });
    }

    return inst;
}

void Solution::print(bool full) {
    // print model name
    cout << "\n----------------------------------------\n";
    cout << "Result of algorithm "<< method <<":\n";

    // print solution status
    cout << "Solution status:\n- ";
    if (feas == 0) {cout << "Infeasible instance";}
    else if (opt) {cout << "Optimal solution\n- obj. value: " << LB;}
    else if (feas == -1 && (method.rfind("BPP-LB", 0) == 0)) { cout << "No easy combined solution\n- LB: " << LB; }
    else if (feas == -1) { cout << "No solution (time limit)"; }
    else if (method.rfind("BPP-UB",0) == 0 || method.rfind("TS",0) == 0) {cout << "Feasible solution\n- UB: " << UB;}
    else { cout << "Suboptimal solution (time limit)\n- LB: " << LB << "\n- UB: " << UB; }
    if ((feas == 1) && ((method.rfind("LF", 0) == 0) || (method.rfind("HF", 0) == 0) || (method.rfind("MFMB", 0) == 0) || (method.rfind("IP2RE", 0) == 0) || (method.rfind("RM2GIFFRE", 0) == 0) || (method.rfind("EM", 0) == 0))) {
        cout << "\n- LP rel. value: " << LPrel;
    }
    cout << "\n";

    // print running time
    cout << "Running time:\n";
    cout << "- total: " << timeT << "s\n";
    cout << "- preprocessing: " << timeP << "s\n";

    if (method.rfind("TS", 0) == 0 && Nconstr == 0) {
        // print number of iterations
        cout << "- number of iterations: " << Nvar << "\n";
    }
    else {
        // print model size
        cout << "Model size:\n";
        cout << "- number of variables: " << Nvar << "\n";
        cout << "- number of constraints: " << Nconstr << "\n";
        cout << "- number of non-zero coefficients: " << Ncoeff << "\n";
    }

    // print feasible solution if it was found and you want to print it
    if (full && feas == 1) {
        vector<int> colorsInBin;
        int col;
        if (opt) cout << "Optimal solution:\n";
        else cout << "Feasible solution:\n";
        for (int b = 0; b < binPacking.size(); b++) {
            cout << "- bin " << b << ": ";
            colorsInBin.clear();
            for (int c = 0; c < binPacking[b].size(); c++) {
                if (binPacking[b][c].size() > 0) {
                    colorsInBin.push_back(c);
                }
            }
            if (colorsInBin.size() == 0) cout << "- \n";
            else {
                for (int cIdx = 0; cIdx < colorsInBin.size(); cIdx++) {
                    col = colorsInBin[cIdx];
                    cout << "color " << col << ": [";
                    for (int j = 0; j < binPacking[b][col].size(); j++) {
                        cout << binPacking[b][col][j][0] << "(" << binPacking[b][col][j][1] << ")";
                        if (j != binPacking[b][col].size() - 1) cout << ", ";
                    }
                    cout << "]";
                    if (cIdx != colorsInBin.size() - 1) cout << ", ";
                }
                cout << "\n";
            }
        }
        cout << "\n";
    }
}

void Solution::draw(const Instance& inst) {
    // This function prints code that can be copied to LaTeX for drawing a solution

    cout << "\\begin{figure}[h]\n";
    cout << "\\centering\n";
    cout << "\\caption{Place a nice caption here}\n";
    double xscale, yscale;
    xscale = min(0.5, double(18)/inst.B);
    yscale = min(0.5, double(18)/inst.W);
    cout << "\\begin{tikzpicture}[xscale = " << xscale << ", yscale = " << yscale << "]\n";
    cout << "\\draw[gray, thin](0, 0) grid(" << inst.B << ", " << inst.W << ");\n";
    int v, w, val;
    for (int b = 0; b < binPacking.size(); b++) {                   // for each bin
        v = 0;                                                      // the current load
        for (int c = 0; c < inst.C; c++) {                          // for every color
            for (int j = 0; j < binPacking[b][c].size(); j++) {     // for every item type
                w = binPacking[b][c][j][0];                         // find the item size
                for (int k = 0; k < binPacking[b][c][j][1]; k++) {  // for every duplicate
                    val = (c + 1) / double(inst.C + 1) * 100;       // determine the color: lower c = lighter color
                    // draw the item              
                    cout << "\\draw[thick, fill = black!" << val << "!white](" << b << ", " << v << ")--(" << b + 1 << ", " << v << ")--(" << b + 1 << ", " << v + w << ")--(" << b << ", " << v + w << ")--(" << b << ", " << v << "); \n";
                    v += w;                                         // increase the bin's load
                }
            }
        }
    }
    cout << "\\end{tikzpicture}\n";
    cout << "\\label{Place a nice label here}\n";
    cout << "\\end{figure}\n\n";
}

string Solution::exportSol() {
    string text = ", " + method + ", " + to_string(max(feas,0)) + ", " + to_string(opt) + ", " + to_string(LPrel) + ", " + to_string(LB) + ", " + to_string(UB) + ", " + to_string(timeT) + ", " + to_string(timeP) + ", " + to_string(Nvar) + ", " + to_string(Nconstr) + ", " + to_string(Ncoeff);
    return text;
}

bool itemOrdering(const vector<int>& v1, const vector<int>& v2) {
    // An item is ordered earlier if its color index is lower
    // or if its size is larger in case that the colors are the same
    return 100000 * v1[0] - v1[1] < 100000 * v2[0] - v2[1];
}

vector<vector<int>> sortAndCount(vector<int> &vec) {
    // This function sorts the entries in a vector from largest to smallest and counts occurences of each value
    vector<vector<int>> result;
    sort(vec.begin(), vec.end(), greater<int>());
    int q = 0; int w = vec[0];
    for (int i = 0; i < vec.size(); i++) {  // loop over all items
        if (vec[i] == w) {                  // if the item size is the same
            q++;
        }
        else {                              // if the item size is different
            result.push_back({ w, q });
            w = vec[i];
            q = 1;
        }
    }
    // finish processing the last item type
    result.push_back({ w, q });
    return result;
}

double getCPUTime() {
    return (double)clock() / CLOCKS_PER_SEC;
}

void printVector(const vector<int>& vec, string name) {
    // This function prints a 1-dimensional vector
    cout << name << ": ";

    if (vec.size() == 0) {
        cout << "[]\n";
    }
    else {
        cout << "[";
        for (size_t i = 0; i < vec.size(); i++) {
            cout << vec[i];
            if (i != vec.size() - 1) {
                cout << " ";
            }
        }
        cout << "]\n";
    }
}

void printVector(const vector<bool>& vec, string name) {
    // This function prints a 1-dimensional vector
    cout << name << ": ";

    if (vec.size() == 0) {
        cout << "[]\n";
    }
    else {
        cout << "[";
        for (size_t i = 0; i < vec.size(); i++)
        {
            cout << vec[i];
            if (i != vec.size() - 1) cout << " ";
        }
        cout << "]\n";
    }
}

void printVector(const vector<double>& vec, string name) {
    // This function prints a 1-dimensional vector
    cout << name << ": ";

    if (vec.size() == 0) {
        cout << "[]\n";
    }
    else {
        cout << "[";
        for (size_t i = 0; i < vec.size(); i++)
        {
            cout << vec[i];
            if (i != vec.size() - 1) cout << " ";
        }
        cout << "]\n";
    }
}

void print2DVector(const vector<vector<int>>& vec, string name, string rowname) {
    // This function prints a 2-dimensional vector
    cout << name << ": \n";
    for (size_t i = 0; i < vec.size(); i++)
    {
        printVector(vec[i], "   " + rowname + " " + to_string(i));
    }
}

void print2DVector(const vector<vector<bool>>& vec, string name, string rowname) {
    // This function prints a 2-dimensional vector
    cout << name << ": \n";
    for (size_t i = 0; i < vec.size(); i++)
    {
        printVector(vec[i], "   " + rowname + " " + to_string(i));
    }
}

void print2DVector(const vector<vector<double>>& vec, string name, string rowname) {
    // This function prints a 2-dimensional vector
    cout << name << ": \n";
    for (size_t i = 0; i < vec.size(); i++)
    {
        printVector(vec[i], "   " + rowname + " " + to_string(i));
    }
}

void print3DVector(const vector<vector<vector<int>>>& vec, string name, string rowname) {
    // This function prints a 3-dimensional vector
    cout << name << ": \n";
    for (size_t i = 0; i < vec.size(); i++) {
        cout << "   " << rowname << " " << i << ": ";
        if (vec[i].size() == 0) {
            cout << "-\n";
            continue;
        }
        cout << "[";
        for (size_t j = 0; j < vec[i].size(); j++) {
            cout << "(";
            for (size_t k = 0; k < vec[i][j].size(); k++) {
                cout << vec[i][j][k];
                if (k != vec[i][j].size() - 1) cout << ",";
            }
            cout << ")";
            if (j != vec[i].size() - 1) cout << ", ";
        }
        cout << "]\n";
    }
}

void print3DVector(const vector<vector<vector<bool>>>& vec, string name, string rowname) {
    // This function prints a 3-dimensional vector
    cout << name << ": \n";
    for (size_t i = 0; i < vec.size(); i++) {
        cout << "   " << rowname << " " << i << ": ";
        if (vec[i].size() == 0) {
            cout << "-\n";
            continue;
        }
        cout << "[";
        for (size_t j = 0; j < vec[i].size(); j++) {
            cout << "(";
            for (size_t k = 0; k < vec[i][j].size(); k++) {
                cout << vec[i][j][k];
                if (k != vec[i][j].size() - 1) cout << ",";
            }
            cout << ")";
            if (j != vec[i].size() - 1) cout << ", ";
        }
        cout << "]\n";
    }
}

void print3DVector(const vector<vector<vector<double>>>& vec, string name, string rowname) {
    // This function prints a 3-dimensional vector
    cout << name << ": \n";
    for (size_t i = 0; i < vec.size(); i++) {
        cout << "   " << rowname << " " << i << ": ";
        if (vec[i].size() == 0) {
            cout << "-\n";
            continue;
        }
        cout << "[";
        for (size_t j = 0; j < vec[i].size(); j++) {
            cout << "(";
            for (size_t k = 0; k < vec[i][j].size(); k++) {
                cout << vec[i][j][k];
                if (k != vec[i][j].size() - 1) cout << ",";
            }
            cout << ")";
            if (j != vec[i].size() - 1) cout << ", ";
        }
        cout << "]\n";
    }
}

void removeLine() {
    #ifdef _WIN32
        cout << "\x1b[1F\x1b[2K"; // remove Gurobi message
        cout << "\x1b[1F\x1b[2K"; // remove Gurobi message
    #endif
}