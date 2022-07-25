#include "func.h"

char *getDynStr(char *str)
{
    char *theStr;
    theStr = (char *)malloc((strlen(str) + 1) * sizeof(char));
    if (!theStr)
        return NULL;

    strcpy(theStr, str);
    return theStr;
}

int *getAllMutIndex(int seqLen, int numOfMut)
{
    int index = 0;
    int *indexs = (int *)malloc(sizeof(int) * (numOfMut * 2));
    for (int i = 0; i < seqLen - 1; i++)
    {
        for (int j = i + 1; j < seqLen; j++)
        {
            //  printf("(%d,%d)\n",i,j);
            indexs[index] = i;
            index++;
            indexs[index] = j;
            index++;
        }
    }
    return indexs;
}
char *readWeightsAndSeq1FromFile(FILE *fp, double *weights)
{
    // get the weights
    for (int i = 0; i < NUM_OF_WEIGHTS; i++)
        fscanf(fp, "%lf", &weights[i]);

    // get seq1
    char temp[MAX_STR_LEN];
    fscanf(fp, "%s", temp);
    char *seq1 = getDynStr(temp);

    // printf("%s\n", seq1);
    return seq1;
}

int getNumOfMut(int seq2Len)
{
    int numOfMut = 0;
    for (int i = 0; i < seq2Len - 1; i++)
    {
        for (int j = i + 1; j < seq2Len; j++)
        {
            numOfMut++;
        }
    }
    return numOfMut;
}

char **readSeq2FromFile(FILE *fp, int *num_of_seq2)
{

    // get num of seq2
    fscanf(fp, "%d", num_of_seq2);
    char **seq2 = (char **)malloc(sizeof(char *) * *num_of_seq2);
    // Get all of seq2
    for (int i = 0; i < *num_of_seq2; i++)
    {
        char temp[MAX_STR_LEN];
        fscanf(fp, "%s", temp);
        seq2[i] = getDynStr(temp);
    }
    return seq2;
}



void writhResultOnFile(int numProc, int **bestMsForAllSeq, int *bestOffsetForAllSeq)
{
    FILE *fptr;
    fptr = fopen("output", "w");
    if (fptr == NULL)
    {
        printf("Error!");
        exit(1);
    }
    char seqNum[12] = {'S', 'e', 'q', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ':', '\n'};
    for (int i = 0; i < numProc; i++)
    {
        fprintf(fptr, "Seq number: %d MS=(%d,%d) best offset=%d\n", i + 1, bestMsForAllSeq[i][0], bestMsForAllSeq[i][1], bestOffsetForAllSeq[i]);
    }
}