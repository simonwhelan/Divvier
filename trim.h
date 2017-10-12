#define addLogProb(x,y) ( (y) > (x) ) ? ((y) += log(1 + exp((x)-(y)))) : ((y)=(x)+log(1 + exp((y)-(x))))

#define MAX_LINE_LEN 2000
#define MAX_FILENAME_LEN 80
#define N_BASES 4
#define N_PEPT 20
#define MAX_LABEL_LEN 20

