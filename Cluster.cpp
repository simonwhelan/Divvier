/*
	A fast clustering algorithm for PPs guided by a tree
*/


#include "Cluster.h"
#include "Tree.h"

using namespace::std;
#include <unordered_set>

// Some initialisation
CCluster *CCluster::_cluster = NULL;
CCluster * CCluster::Instance() {
	if(!_cluster) {
		_cluster = new CCluster();
	}
	return _cluster;
}

// MakePairs()
// ---
// Computes the distance matrix from the tree and uses it to select a set of pairs for testing during divvying
// Gets the splits from the tree
void CCluster::MakePairs() {
	assert(_tree.NoSeq() == NoSeq());
	assert(_splitPairs.empty());
	// Now midpoint root the tree
	if(_doUPGMA) { _tree.MidpointRoot(); }
	// Get the splits
	vector <SSplit> splits = _tree.BuildSplits();
	sort(splits.begin(),splits.end(),[](auto const &a, auto const &b) {
		if(a.Left.size() == 1) { assert(a.Left[0] == 0); return true; }
		if(a.Right.size() == 1 && b.Right.size() == 1) { return a.Right[0] < b.Right[0]; }
		return my_max(a.Left.size(),a.Right.size()) > my_max(b.Left.size(),b.Right.size());
	});
	// Distance matrix
	vector <double> distances = _tree.GetTreePW();
	vector <int> big, small;
	vector <tuple <int , int , double> > distSet;			// The set of pairwise comparisons (int i, int j) and their distance (double)
	vector <vector <int> > pairs2add;
	for(SSplit &split : splits) {
		// Get the full list of pairwise comparisons in the splits
		if(split.Left.size() > split.Right.size()) {
			big = split.Left; small = split.Right;
		} else {
			big = split.Right; small = split.Left;
		}
		// Create the list of pairwise comparisons sorted by distance. x always small; y alway big
		for(int &x : small) {
			for(int &y : big) {
				distSet.push_back(make_tuple(x,y,distances[(y*NoSeq()) + x]));
			}
		}
		sort(begin(distSet), end(distSet), [](auto const &t1, auto const &t2) {
			return get<2>(t1) < get<2>(t2);
		});
		// Obtain samples. The rules are:
		// 1. Unique sequences from the small data set are the priority, each with a different partner
		vector <int> small_count(NoSeq(),0),big_count(NoSeq(),0);
		int small_max = ceil( ( (double) _approxNumber / (double) small.size() ) +1 );
		int big_max = ceil( ( (double) _approxNumber / (double) big.size() ) + 1);
		for(int i = 0; i < distSet.size(); i++) {
			if(small_count[get<0>(distSet[i])] >= small_max) { continue; }
			if(big_count[get<1>(distSet[i])] >= big_max) { continue; }
			pairs2add.push_back(vector<int>{get<0>(distSet[i]),get<1>(distSet[i])});
			small_count[get<0>(distSet[i])]++;
			big_count[get<1>(distSet[i])]++;
			if(pairs2add.size() >= _approxNumber) { break; }	// Finish when we have the full list
		}
		_splitPairs.push_back( tuple<SSplit ,vector <vector <int> > >(split , pairs2add) );
		// Clean up
		pairs2add.clear();
		distSet.clear();
		big.clear();
		small.clear();
	}

	// Now extend the pair sets to maximise coverage
	_all_pairs = PairsToCalculate();
	vector <vector <int> > workingPairs;
	for(auto &split : _splitPairs) {
		get<1>(split).clear();
		small = get<0>(split).Left; big = get<0>(split).Right;
		for(auto &pair : _all_pairs) {
			if((find(small.begin(), small.end(),pair[0]) != small.end() && find(big.begin(), big.end(), pair[1]) != big.end()) ||
			 (find(small.begin(), small.end(),pair[1]) != small.end() && find(big.begin(), big.end(), pair[0]) != big.end())) {
				get<1>(split).push_back(pair);
			}
		}
	}
	_ready = true;
}

// Get the set of pairs for testing split numbered splitNum
vector <vector <int>> CCluster::GetPairs(int splitNum) {
	return get<1>(_splitPairs[splitNum]);
}

vector <vector <int> > CCluster::PairsToCalculate() {

	if(!_all_pairs.empty()) { return _all_pairs; }
	vector <vector <int> > retVec;
	// Get the list. Note only upper 1/2 of diagonal
	for(auto x : _splitPairs) {
		for(auto y : get<1>(x)) {
			if(y[0] > y[1]) {
				retVec.push_back( vector<int>{y[1],y[0]} );
			} else {
				retVec.push_back( vector<int>{y[0],y[1]} );
			}
		}
	}
	if(retVec.size() == 0) { cout << "\nCCluster::PairsToCalculate() :: Need to have pairs to calculate... Have you initialised properly? Is the Tree wrong?\n"; exit(-1); }
	// Sort it
	sort(retVec.begin(), retVec.end(), [](auto x, auto y) {
		return (x[1] < y[1] && x[0] < y[0]);
	});
	sort(retVec.begin(), retVec.end(), [](vector <int> x, vector <int> y) {
		if(x[0] < y[0]) { return true; }
		if(x[0] > y[0]) { return false; }
		if(x[1] < y[1]) { return true; }
		return false;
	});
	// Remove redundancy
	for(int i = 0 ; i < retVec.size() - 1; i++) {
		if(retVec[i+1][0] == retVec[i][0] && retVec[i+1][1] == retVec[i][1]) {
			retVec.erase(retVec.begin() + i + 1);
			i--;
		}
	}
	return retVec;
}

#define DEBUG_TESTSPLIT 0
//vector < vector <int> > CCluster::OutputClusters(vector <double> PPs, double threshold, int clusterMethod) {
vector < vector <int> > CCluster::OutputClusters(vector <double> PPs, string seq, double threshold, int clusterMethod) {
	assert(_ready);
	vector <vector <int> > retSplits;
	// Do appropriate clustering
	if(_doUPGMA) {
		PseudoUPGMA(retSplits,PPs,seq,threshold);
	} else {
		SmartDivisive(retSplits,PPs,seq,threshold);
	}
	// Sort them so they're in a nice order; the structure of hte tree splits mean they might not be
	sort(retSplits.begin(), retSplits.end(),[](const vector<int>& a, const vector<int>& b) {
	  return a[0] < b[0];
	});
	return retSplits;
}

void CCluster::SmartDivisive(vector <vector <int> > &retSplits, vector <double> &PPs, string seq, double threshold) {
	vector <int> starter(NoSeq(),0);
	for(int i = 0; i < NoSeq(); i++) { starter[i] = i; }
	retSplits.push_back(starter);
	// Find what clusters to make. Can change this function
	for(int i = 0 ; i < _splitPairs.size(); i++) {
#if DEBUG_TESTSPLIT == 1
		cout << "\n---\nCurrent splits\n";
		for(auto &out : retSplits) {
			cout << " | ";
			for(auto &v : out) { cout << "[" << v << "]" << seq[v] << " " << flush; }
		}
#endif
		// Debug code to do filtering approach
//		if(_splits[i].Left.size() != 1 && _splits[i].Right.size() != 1) { continue; }
		if(!TestSplit(i,retSplits, seq, threshold,PPs)) {
			retSplits = AddSplit(i,retSplits);
		}
	}
}

void CCluster::PseudoUPGMA(vector <vector <int> > &retSplits, vector <double> &PPs, string seq, double threshold) {
	assert(_ready);
	vector <int> starter(NoSeq(),0);
	for(int i = 0; i < NoSeq(); i++) { starter[i] = i; }
	retSplits.push_back(starter);
	// Find what splits to perform. True means the split will be applied
	vector <bool> DoSplits(_splitPairs.size(),true);
	for(int i = 0 ; i < _splitPairs.size(); i++) {
#if DEBUG_TESTSPLIT == 1
		cout << "\n---\nTrying split";
		cout << "\n["<<i<<"]: " << get<0>(_splitPairs[i]).Left << " | " << get<0>(_splitPairs[i]).Right;
#endif
		// Only consider groupings of 2 or more for UPGMA
		if(get<0>(_splitPairs[i]).Left.size() == 1 || get<0>(_splitPairs[i]).Right.size() == 1) {
			if(get<0>(_splitPairs[i]).Left.size() == 1) { assert(get<0>(_splitPairs[i]).Left[0] == i); }
			else { assert(get<0>(_splitPairs[i]).Right[0] == i); }
			continue;
		}

		if(TestSplit(i,retSplits, seq, threshold,PPs)) {
			DoSplits[i] = false;
			// Do the sequences below on the descendant part of split
			vector <int> seq2Join;
			if(get<0>(_splitPairs[i]).rootLeft) { seq2Join = get<0>(_splitPairs[i]).Right; }
			else { seq2Join = get<0>(_splitPairs[i]).Left; }
			for(int x : seq2Join) { DoSplits[x] = false; }
		}
	}
	// Add the existing splits
	for(int i = 0 ; i < _splitPairs.size(); i++) {
//		if(get<0>(_splitPairs[i]).Left.size() == 1 || get<0>(_splitPairs[i]).Right.size() == 1) { continue; }
		if(DoSplits[i]) {
			retSplits = AddSplit(i,retSplits);
		}
	}
}

// Calculates a test statistic based on PPs and compares it to threshold. If greater it passes and returns true
bool CCluster::TestSplit(int split2Test, vector <vector <int> > &curSplit, string seq, double threshold, vector <double> &PPs, int testMethod) {
	double testStat = 0;
	int similarity_count = 0;	// Proportion of pairs sharing same character

	// Statistics relating to the PPs used to assess the PP
	vector <double> splitPPs;	// The PPs for this split
	int numPast = 0;			// The number past the threshold
	// New stuff testing only the values in the present split when possible
	bool testActive = true;
	vector <int> activeSplit;
	vector <double> activePPs;
	for	(auto &v : curSplit) {
		if(TestSubsplit(split2Test,v)) {activeSplit = v; break;}
	}
	// Check all gaps
	bool leftOkay = false, rightOkay = false;
	for(auto  &v : get<0>(_splitPairs[split2Test]).Left) {
		if(seq[v] == '-') { continue; }
		leftOkay = true;
	}
	for(auto  &v : get<0>(_splitPairs[split2Test]).Right) {
			if(seq[v] == '-') { continue; }
			rightOkay = true;
		}
	if(!(leftOkay && rightOkay)) { return true;}

#if DEBUG_TESTSPLIT == 1
		cout << "\nseq: " << seq;
		cout << "\nTesting split["<< split2Test << "] " << get<0>(_splitPairs[split2Test]).Left << " | " << get<0>(_splitPairs[split2Test]).Right;
		cout << "\nActive split["<<activeSplit.size()<<"]: "<< activeSplit;
		cout << "\nPP";
		for(vector <int>  &v : get<1>(_splitPairs[split2Test])) {
			if(seq[v[0]] == '-' || seq[v[1]] == '-') { continue; }
			cout << "  [" << v[0] << seq[v[0]]<< "," << v[1] << seq[v[1]] << "]" << PPs[(v[0] * NoSeq()) + v[1]];
		}
#endif

	// Collect the PPs
	for(vector <int>  &v : get<1>(_splitPairs[split2Test])) {
		assert(v.size() == 2);
		if(seq[v[0]] == '-' || seq[v[1]] == '-') { continue; }
		// Normal statistic
		splitPPs.push_back(PPs[(v[0] * NoSeq()) + v[1]]);
		// Statistic where only the active split is considered
		if(testActive)  {
			if(find(activeSplit.begin(),activeSplit.end(),v[0]) != activeSplit.end() && find(activeSplit.begin(),activeSplit.end(),v[1]) != activeSplit.end())
				{ activePPs.push_back(PPs[(v[0] * NoSeq()) + v[1]]); }
		}
		// Similarity statistic
		if(seq[v[0]] == seq[v[1]] ) { similarity_count++; }
	}
	if(splitPPs.size() == 0) { // If there's no PP pairs then no evidence either way and go with MSA
		_warningNoInfo = true;
#if DEBUG_TESTSPLIT == 1
		cout << "\nCount = 0 return";
#endif
		return _acceptNoInfo;
	}
	// Similarity check
	if(_doSimilarityCheck) {
		if((double) similarity_count / (double) splitPPs.size() > _similarityCutOff) { return true; } 			// High similarity columns always returned true
	}


	switch(testMethod) {
	case 0: 				// Compute the mean and compare to threshold
		numPast = 0;
		for(auto  &v : splitPPs) {
			if(v > _tightThreshold) { numPast++; }
			testStat += v;
		}
		if(splitPPs.size() == 0) { // If there's no PP pairs then no evidence either way and go with MSA
			_warningNoInfo = true;
#if DEBUG_TESTSPLIT == 1
			cout << "\nCount = 0 return";
#endif
			return _acceptNoInfo;
		}
#if DEBUG_TESTSPLIT == 1
		cout << "\nNaiveStat= " << testStat / (double) get<1>(_splitPairs[split2Test]).size() << " ; CorrectStat: " << testStat / (double)splitPPs.size();
		if(!activePPs.empty()) { cout << " ; activeStat: " << (double) Sum(&activePPs) / (double) activePPs.size(); }
#endif

		if(numPast >= _numberPastThreshold && _numberPastThreshold > 0) {
#if DEBUG_TESTSPLIT == 1
			cout << "\nThreshold (" << numPast << " / " << _numberPastThreshold << ") return";
#endif
			return true;
		}
		// Decision between the active or the normal stat
		if(activePPs.size() >= activeSplit.size() - 1 || activePPs.size() > 2) {
			if( ( (double) Sum(&activePPs) / (double) activePPs.size() ) + DBL_EPSILON >= threshold) { return true; }
		} else {
			if((testStat / (double) splitPPs.size()) + DBL_EPSILON >= threshold) { return true; }
		}
		break;

	default:
		cout << "\nUnknown testMethod passed to CCluster::TestSplit(...)\n"; exit(-1);
	};
#if DEBUG_TESTSPLIT == 1
			cout << "\nsplit rejected return";
#endif
	return false;
}

vector <vector <int> > CCluster::AddSplit(int split2Add, vector <vector <int> > &curSplit) {
	vector <vector <int> > retSplit;
	// Find the split in curSplit that has elements from both Left and Right _split(split2Add) present
	for(vector <int> &split : curSplit) {
		if(!TestSubsplit(split2Add,split)) {
			retSplit.push_back(split);
			continue;
		}
		vector <int> newSplit;
		for(int &left : get<0>(_splitPairs[split2Add]).Left) {
			if(find(split.begin(),split.end(),left) != split.end()) { newSplit.push_back(left); }
//			if(split_set.find(left) != split_set.end()) { newSplit.push_back(left); }
		}
		assert(newSplit.size() > 0);
		retSplit.push_back(newSplit);		// Note: no sorting required because they're pre-sorted
		newSplit.clear();
		for(int &right : get<0>(_splitPairs[split2Add]).Right) {
			if(find(split.begin(),split.end(),right) != split.end()) { newSplit.push_back(right); }
//			if(split_set.find(right) != split_set.end()) { newSplit.push_back(right); }
		}
		assert(newSplit.size() > 0);
		retSplit.push_back(newSplit);		// Note: no sorting required because they're pre-sorted
	}
	return retSplit;
}
/*
vector <vector <int> > CCluster::JoinSplit(int split2Join, vector <vector <int> > &curSplit) {
	vector <vector <int> > retSplit;
	for(vector <int> &split : curSplit) {

	}
}
*/
// Function that returns the current set of sequences affected by split2get
bool CCluster::TestSubsplit(int split2test, vector <int> &testSplit) {
	unordered_set<int> split_set ( testSplit.begin(),testSplit.end() );
	bool inLeft = false, inRight = false;
	for(int &left : get<0>(_splitPairs[split2test]).Left) {
		if(split_set.find(left) != split_set.end()) { inLeft = true; break; }
	}
	for(int &right : get<0>(_splitPairs[split2test]).Right) {
		if(split_set.find(right) != split_set.end()) { inRight = true; break; }
	}
	if(!inLeft || !inRight) { return false; }
	return true;
}


