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


int main(int argc,char **argv){ 
  if(argc < 2){
    printf("Input file not specified\n");
    exit(EXIT_FAILURE);
  }

  clock_t start, end;
  double cpu_time_used;
  start = clock();


  // Hard coded options
  JTT = 0;
  PMB = 0;
  PAM = 1;
  MATRICES = 0;

  readSeq(argv[argc-1]);

  printf("%d  %d\n",Nseq,alen);
  for(int i = 0 ; i < Nseq ; i++) {
	  printf("%s\n",names[i]);
  }

  initHMM(alen);
  calc_posterior(alen);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  fprintf(stderr,"CPU time used: %.2f seconds\n",cpu_time_used);
  return 1;
  
}
