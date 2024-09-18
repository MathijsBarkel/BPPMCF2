#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <sstream> 
#include "gurobi_c++.h"

const double EPSILON = 0.00001; // small constant
const int M = 1000000;			// big constant

struct Instance {
	// instance attributes
	string name;							  // name of the instance
	int I; 									  // number of items
	int B; 									  // number of bins
	int C;									  // number of colors
	int W;									  // bin capacity
	int J;		   							  // number of different item types					
	vector<vector<int>> items;				  // color and weight of each item
	vector<vector<int>> colorlessItems;		  // sizes and quantities of all items (ignoring color)
	vector<vector<int>> countedItems;		  // items counted acc. to color and size
	vector<vector<vector<int>>> groupedItems; // items grouped acc. to color
	vector<int> LBs;						  // the L2 lower bound for each color seperately
	int LBtot;								  // the sum of the L2 bounds
	int Bmin;							      // the minimum number of required bins
	vector<vector<vector<vector<int>>>> feasSol; // a feasible solution

	// instance methods
	void preprocessing(bool reset = false, bool reverse = false);
	void findL2();
	void findMinNumberOfBins(string filename);
	vector<vector<int>> findColorSubset(const vector<bool>& subset);
	void print();
	void draw();

	//void sortItems();
	//void countItems();
	//void groupItems();
	//void discardColorSort();
};

struct Solution {
	string method;
	bool opt = false;
	int feas = -1, LB = 0, UB = 1000000, Nvar = 0, Nconstr = 0, Ncoeff = 0;
	double LPrel = 0;
	double timeP = 0, timeLP = 0, timeT = 0;
	vector<vector<vector<vector<int>>>> binPacking;
	void print(bool full=false);
	void draw(const Instance& inst);
	string exportSol();
};

Instance readInstance(string filename);
Instance createInstance(string name, int I, int B, int C, int W, const vector<vector<int>>& items);
Instance generateInstance(string name, int I, int B, int C, int W, int r=0);
bool itemOrdering(const vector<int>& v1, const vector<int>& v2);
vector<vector<int>> sortAndCount(vector<int> &vec);
double getCPUTime();
void printVector(const vector<int>& vec, string name);
void printVector(const vector<bool>& vec, string name);
void printVector(const vector<double>& vec, string name);
void print2DVector(const vector<vector<int>>& vec, string name, string rowname);
void print2DVector(const vector<vector<bool>>& vec, string name, string rowname);
void print2DVector(const vector<vector<double>>& vec, string name, string rowname);
void print3DVector(const vector<vector<vector<int>>>& vec, string name, string rowname);
void print3DVector(const vector<vector<vector<bool>>>& vec, string name, string rowname);
void print3DVector(const vector<vector<vector<bool>>>& vec, string name, string rowname);
void removeLine();
#endif