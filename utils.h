/////////////////////////////////////////////////////////////////
// utils.h
//
// Default constants for use in AMAP.  The emission
// probabilities were computed using the program used to build
// the BLOSUM62 matrix from the BLOCKS 5.0 dataset.  Transition
// parameters were obtained via unsupervised EM training on the
// BALIBASE 2.0 benchmark alignment database.
/////////////////////////////////////////////////////////////////



int readSeq(char *inFile);
void error(char *fmt, ... );

EXTERN int Nseq;
EXTERN char **align;
EXTERN char **sequence;
EXTERN char **names;
EXTERN int *lens;
EXTERN int alen;
EXTERN double TOT_DIST;
EXTERN double **dists;
EXTERN int uguide; // User provided guide tree
EXTERN char guidetree[200];
#define MIN(X,Y) ((X) < (Y) ? : (X) : (Y))
