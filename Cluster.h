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
	bool AddTree(std::string Tree) {
		if(Instance()->_tree.NoSeq() > 0) { return false; }
		Instance()->_tree = CTree(Tree,_names);
		MakePairs();
		return true;
	}
												// Second initiation function that gets the splits from the tree
	void SetOptions(bool acceptNoInfo, int approxNumber, bool partialFilter) {		// Function for setting the approximation stuff
		_acceptNoInfo = acceptNoInfo;
		_approxNumber = approxNumber;
		_doDivvying = !partialFilter;
	}
	// Access for the pairs needed for calculation
	vector <vector <int> > PairsToCalculate();							// Gets the list of pairs to calculate
	int NoPairs() { return _all_pairs.size(); }
	// The clustering methods
	std::vector <std::vector <int> > OutputClusters(vector <double> PPs, string seq, double threshold);		// The result of the clustering

	// Whether to output a warning or not
	bool Warning() { return _warningNoInfo; }

private:
	bool _ready = false;
	// Type of filtering
	bool _doDivvying = true;						// TRUE = divvying; FALSE = partialFilter;
	// Some heuristic stuff
	// 1. Accept with no info -- accept a split if there's no PP to support or refute it; important in sparse alignments
	bool _acceptNoInfo = false;
	bool _warningNoInfo = false;
	// 2. The number of pairwise comparisons to make for each split
	int _approxNumber = 10;
	// Variables
	static CCluster * _cluster;					// Singleton
	vector <tuple< SSplit, vector <vector <int> > > > _splitPairs;	// The splits and the sets of pairs that define them
	std::vector <string> _names;				// Sequence names (defines order in tree and MSA)
	std::vector <vector <int> > _all_pairs;		// The full set of non-redundant pairs that need to be calculated

	CTree _tree;

	// Functions
	int NoSeq() { return _names.size(); }
	// Values linked to the clustering
	void MakePairs();																			// Calculates the pairs needed for clustering
	vector <vector <int> > GetPairs(int splitNum);												// Get the pairs for a specific split number
	bool TestSplit(int split2Test, vector <vector <int> > &curSplit, string seq, double threshold, vector <double> &PPs);
																								// Function that uses PPs to test a specific split
	double ScoreSplit(tuple <SSplit, vector <vector <int> > > split2Test, vector <vector <int> > &curSplit, string seq, vector <double> &PPs);
	bool TestSubsplit(int split2test, vector <int> &testSplit);									// Test whether the subsplit affected by split2get
	bool TestSubsplit(SSplit split2test,vector <int> &testSplit);
	vector <vector <int> > AddSplit(int split2Add, vector <vector <int> > &curSplit);			// Adds the _split(split2Add) to the current set of splits

	// Alternative clustering methods
	void SmartDivisive(vector <vector <int> > &RetClusters, vector <double> &PPs, string seq, double threshold);
};
