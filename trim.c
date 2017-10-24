#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#define EXTERN
#include "trim.h"
#undef EXTERN

#define EXTERN extern
#include "utils.h"
#include "hmm.h"
#include "matrices.h"
#undef EXTERN


void ReadAndPrepZorro(char * file){

  // Hard coded options
  JTT = 0;
  PMB = 0;
  PAM = 1;
  MATRICES = 0;

  readSeq(file);

  initHMM(alen);
  calc_prep(alen);
  // calc_posterior(alen);

  
}
