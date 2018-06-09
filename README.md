## Divvier: a program for removing MSA uncertainty by Simon Whelan

Divvier in the beta stages of development. It is not released and the methodology is is in the process of being written up for publication. If you wish to use it then please contact me and I can ensure that you use it correctly. (It's currently very easy to interpret its output incorrectly.) 

For now the recommended usage is along the lines of:
```
./divvier -partial -mincol 5 myMSAfile.fas
```
This will give you a partially filtered MSA (individual characters removed; _not_ columns) with at least 5 characters in each column. 

The following options are available in divvier

Clustering options:
	-divvy       : do standard divvying (DEFAULT)
	-partial     : do partial filtering by testing removal of individual characters
	-thresh X    : set the threshold for divvying to X (DEFAULT = 0.801)

Approximation options: 
	-approx X    : minimum number of characters tested in a split during divvying (DEFAULT X = 10)
	-checksplits : go through sequence and ensure there's a pair for every split. Can be slow
	-HMMapprox   : Do the pairHMM bounding approximation (DEFAULT)
	-HMMexact    : Do the full pairHMM and ignore bounding

Output options: 
	-mincol X    : Minimum number of characters in a column to output when divvying/filtering (DEFAULT X = 2)
	-divvygap    : Output a gap instead of the static * character so divvied MSAs can be used in phylogeny

