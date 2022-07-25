#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#define FILE_NAME "input.txt"
#define MAX_STR_LEN 5000
#define NUM_OF_WEIGHTS 4

char *getStrExactName(FILE *fp);
char *getDynStr(char *str);
char *readWeightsAndSeq1FromFile(FILE *fp, double *weights);
char **readSeq2FromFile(FILE *fp, int *num_of_seq2);
int *getAllMutIndex(int seqLen, int numOfMut);
int getNumOfMut(int seq2Len);
int *getSeqBestMsAndOffsetOnGpu(char *seq2, char *seq1, int seq1Len, int seq2Len, int *allIndexsToMakeMut, double numOfMut, int numOfAlignment, double *weights, int *bestOffset);
void writhResultOnFile(int numProc, int **bestMsForAllSeq,int *bestOffsetForAllSeq);