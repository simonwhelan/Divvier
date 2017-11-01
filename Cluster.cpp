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

#define DEBUG_TESTSPLIT 0
#define DEBUG_SCOREDETAIL 0

// MakePairs()
// ---
// Computes the distance matrix from the tree and uses it to select a set of pairs for testing during divvying
// Gets the splits from the tree
void CCluster::MakePairs() {
	assert(_tree.NoSeq() == NoSeq());
	assert(_splitPairs.empty());
	// Get the splits; no sorting required because of greedy algorithm
	vector <SSplit> splits = _tree.BuildSplits();
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
		// 		Unique sequences from the small data set are the priority, each with a different partner
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
	sort(retVec.begin(), retVec.end(), [](vector <int> x, vector <int> y) { // Sort both columns
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

//vector < vector <int> > CCluster::OutputClusters(vector <double> PPs, double threshold, int clusterMethod) {
vector < vector <int> > CCluster::OutputClusters(vector <double> PPs, string seq, double threshold) {
	assert(_ready);
	vector <vector <int> > retSplits;
	// Do clustering
	SmartDivisive(retSplits,PPs,seq,threshold);
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
	// Greedily remove splits, worst first
	vector <bool> splitDone(_splitPairs.size(),false);
#if DEBUG_TESTSPLIT == 1
	cout << "\n--- Greedy splitting ----";
#endif
	while(true) {
		tuple <int,double> minScore(-1,1.0);
		// Get score for remaining splits
		for(int i = 0; i < _splitPairs.size(); i++) {
			if(splitDone[i]) { continue; }
			double score = ScoreSplit(_splitPairs[i],retSplits,seq,PPs);
			if(score + DBL_EPSILON < threshold && score < get<1>(minScore)) {
				get<0>(minScore) = i; get<1>(minScore) = score;
			}
		}
		if(get<0>(minScore) < 0) { break; } // No minimum score found
		splitDone[get<0>(minScore)] = true;
		retSplits = AddSplit(get<0>(minScore),retSplits);
#if DEBUG_TESTSPLIT == 1
		cout << "\nAdded split[" << get<0>(minScore) << "] == " << get<1>(minScore) << ": " << get<0>(_splitPairs[get<0>(minScore)]).Left << " | " << get<0>(_splitPairs[get<0>(minScore)]).Right;
		cout << "\nNew splits ";
		for(auto &out : retSplits) {
			cout << " | ";
			for(auto &v : out) { cout << "[" << v << "]" << seq[v] << " " << flush; }
		}
#endif
	}
}

// Calculates a test statistic based on PPs and compares it to threshold. If greater it passes and returns true
bool CCluster::TestSplit(int split2Test, vector <vector <int> > &curSplit, string seq, double threshold, vector <double> &PPs) {
	double score = ScoreSplit(_splitPairs[split2Test],curSplit,seq,PPs);
	if(score + DBL_EPSILON > threshold) { return true; }
	return false;
}

double CCluster::ScoreSplit(tuple <SSplit, vector <vector <int> > > split2Test, vector <vector <int> > &curSplit, string seq, vector <double> &PPs) {
	// Statistics relating to the PPs used to assess the PP (backup statistic)
	vector <double> splitPPs;	// The PPs for this split
	vector <int> activeSplit;	// The split actually being made
	vector <double> activePPs;	// The active PPs for this split
	for	(auto &v : curSplit) {
		if(TestSubsplit(get<0>(split2Test),v)) { activeSplit = v; break;}
	}
	// Check all gaps
	bool leftOkay = false, rightOkay = false;
	for(auto  &v : get<0>(split2Test).Left) {
		if(seq[v] == '-') { continue; }
		leftOkay = true;
	}
	for(auto  &v : get<0>(split2Test).Right) {
			if(seq[v] == '-') { continue; }
			rightOkay = true;
		}
	if(!(leftOkay && rightOkay)) { return 1.0;}

#if DEBUG_SCOREDETAIL == 1
		cout << "\nseq: " << seq;
		cout << "\nTesting split " << get<0>(split2Test).Left << " | " << get<0>(split2Test).Right;
		cout << "\nActive split["<<activeSplit.size()<<"]: "<< activeSplit;
		cout << "\nThere are " << get<1>(split2Test).size() << " pairs";
		cout << "\nPP";
		for(vector <int>  &v : get<1>(split2Test)) {
			if(seq[v[0]] == '-' || seq[v[1]] == '-') { continue; }
			cout << "  [" << v[0] << seq[v[0]]<< "," << v[1] << seq[v[1]] << "]" << PPs[(v[0] * NoSeq()) + v[1]];
		}
#endif

	// Collect the PPs
	for(vector <int>  &v : get<1>(split2Test)) {
		assert(v.size() == 2);
		if(seq[v[0]] == '-' || seq[v[1]] == '-') { continue; }
		// Normal statistic
		splitPPs.push_back(PPs[(v[0] * NoSeq()) + v[1]]);
		// Statistic where only the active split is considered
		if(find(activeSplit.begin(),activeSplit.end(),v[0]) != activeSplit.end() && find(activeSplit.begin(),activeSplit.end(),v[1]) != activeSplit.end())
			{ activePPs.push_back(PPs[(v[0] * NoSeq()) + v[1]]); }
	}
	// If there's no PP pairs then no evidence either way and go with _acceptNoInfo
	if(splitPPs.size() == 0) {
		_warningNoInfo = true;
		if(_acceptNoInfo) { return 1.0; } else { return 0.0; }
	}
#if DEBUG_SCOREDETAIL == 1
	cout << "\nCompleteStat: " << Sum(&splitPPs) / (double)splitPPs.size();
	if(!activePPs.empty()) { cout << " ; activeStat: " << (double) Sum(&activePPs) / (double) activePPs.size(); }
#endif
	// Decision between the active or the normal stat
	if(activePPs.size() >= activeSplit.size() - 1 || activePPs.size() > 2) {
		return ( (double) Sum(&activePPs) / (double) activePPs.size() );
	}
	return (Sum(&splitPPs) / (double) splitPPs.size());
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
		}
		assert(newSplit.size() > 0);
		retSplit.push_back(newSplit);		// Note: no sorting required because they're pre-sorted
		newSplit.clear();
		for(int &right : get<0>(_splitPairs[split2Add]).Right) {
			if(find(split.begin(),split.end(),right) != split.end()) { newSplit.push_back(right); }
		}
		assert(newSplit.size() > 0);
		retSplit.push_back(newSplit);		// Note: no sorting required because they're pre-sorted
	}
	return retSplit;
}

// Function that returns the current set of sequences affected by split2get
bool CCluster::TestSubsplit(int split2test, vector <int> &testSplit) {
	return TestSubsplit(get<0>(_splitPairs[split2test]),testSplit);
}
bool CCluster::TestSubsplit(SSplit split2test,vector <int> &testSplit) {
	unordered_set<int> split_set ( testSplit.begin(),testSplit.end() );
	bool inLeft = false, inRight = false;
	for(int &left : split2test.Left) {
		if(split_set.find(left) != split_set.end()) { inLeft = true; break; }
	}
	for(int &right : split2test.Right) {
		if(split_set.find(right) != split_set.end()) { inRight = true; break; }
	}
	if(!inLeft || !inRight) { return false; }
	return true;
}


