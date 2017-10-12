

#define ADD 1	// New
#define LONG
#define ENDTRANS
#define NCSTATE
#define WEIGHT
#define WEIGHTING 1 // New
#define verbose 1 // new


EXTERN double initDistribDefault[];
EXTERN int NUMSTATE;

#ifndef LONG 
EXTERN double gapOpenDefault;
EXTERN double gapExtendDefault;
#else
EXTERN double gapOpenDefault;
EXTERN double gapOpenDefaultL;
EXTERN double gapExtendDefault;
EXTERN double gapExtendDefaultL;
#endif



EXTERN double selfM;
EXTERN double MtoXY;
EXTERN double selfXY;
EXTERN double XYtoM;
#ifdef LONG
EXTERN double MtoXYL;
EXTERN double selfXYL;
EXTERN double XYLtoM;
#endif
#ifdef NCSTATE
EXTERN double MtoC;
EXTERN double NtoM;
EXTERN double selfN;
EXTERN double selfC;
#endif

EXTERN double emitSingleDefault[20];
EXTERN double emitPairsDefault[20][20];



void initHMM(int len);
void calc_posterior(int len);
