/*
	A fast clustering algorithm for PPs guided by a tree
*/


#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <numeric>
#include <regex>
#include <algorithm>
#include <iomanip>
#include "Tree.h"

class CCluster{
public:
	static CCluster *Instance();				// Singleton access
	// Initiation
	bool AddNames(std::vector <string> Names) { if(!_names.empty()) { return false; }_names = Names; return true; };
												// First initiation function that reads the sequence names
	bool AddTree(std::string Tree) { if(_tree.NoSeq() > 0) { return false; } _tree = CTree(Tree,_names); MakePairs(); return true; }
												// Second initiation function that gets the splits from the tree
	// Access for the pairs needed for calculation
	vector <vector <int> > PairsToCalculate();							// Gets the list of pairs to calculate
	int NoPairs() { return _all_pairs.size(); }
	// The clustering methods
	std::vector <std::vector <int> > OutputClusters(vector <double> PPs, string seq, double threshold, int clusterMethod = 0);		// The result of the clustering

	// Whether to output a warning or not
	bool Warning() { return _warningNoInfo; }

private:
	bool _ready = false;
	// Some heuristic stuff
	// 1. Similarity check -- accept columns if there's a certain proportion of the same character; usually a bad idea
	bool _doSimilarityCheck = false;			// Whether to do a similarity check rather than use just PPs
	double _similarityCutOff = 0.8;				// The cutoff used
	// 2. Accept with no info -- accept a split if there's no PP to support or refute it; important in sparse alignments
	bool _acceptNoInfo = false;
	bool _warningNoInfo = false;
	// 3. Number threshold
	double _tightThreshold = 0.95;				// The tight threshold considered guarantee of relationship
	int _numberPastThreshold = -1;				// Simple number that must pass threshold (extension of single linkage clustering); -1 used mean
	// 4. The number of pairwise comparisons to make for each split
	int _approxNumber = 15;
	// Variables
	static CCluster * _cluster;					// Singleton
	std::vector <string> _names;				// Sequence names (defines order in tree and MSA
	std::vector <SSplit> _splits;				// The tree splits (not including trivial splits)
	std::vector <std::vector<std::vector<int> > > _pairs;
												// The set of pairs used to assess each split
	std::vector <vector <int> > _all_pairs;		// The full set of non-redundant pairs that need to be calculated

	CTree _tree;

	// Functions
	int NoSeq() { return _names.size(); }
	// Values linked to the clustering
	void MakePairs();																			// Calculates the pairs needed for clustering
	vector <vector <int> > GetPairs(int splitNum);												// Get the pairs for a specific split number
	bool TestSplit(int split2Test, vector <vector <int> > &curSplit, string seq, double threshold, vector <double> &PPs, int testMethod = 0);
																								// Function that uses PPs to test a specific split
	bool TestSubsplit(int split2test, vector <int> &testSplit);									// Test whether the subsplit affected by split2get
	vector <vector <int> > AddSplit(int split2Add, vector <vector <int> > &curSplit);			// Adds the _split(split2Add) to the current set of splits

};
