#include "func.h"

int main(int argc, char *argv[])
{
    double t1, t2;

    int my_rank;       /* rank of process */
    int numProc;       /* number of processes */
    MPI_Status status; /* return status for receive */

    /* start up MPI */
    MPI_Init(&argc, &argv);

    /* find out process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    t1 = MPI_Wtime();

    double weights[4];
    FILE *fp;
    if ((fp = fopen(FILE_NAME, "r")) == 0)
    {
        printf("cannot open file %s for reading\n", FILE_NAME);
        exit(0);
    }
    char *seq1 = readWeightsAndSeq1FromFile(fp, weights);
    int **bestMsForAllSeq = (int **)malloc(sizeof(int *) * numProc);
    for (int i = 0; i < numProc; i++)
    {
        bestMsForAllSeq[i] = (int *)malloc(sizeof(int) * 2);
    }
    int *bestOffsetForAllSeq = (int *)malloc(sizeof(int) * numProc);

    if (my_rank == 0)
    {
        int num_of_seq2, offset;
        int seq2len, bestOffset, seq1Len;

        char **seq2 = readSeq2FromFile(fp, &num_of_seq2);
        for (int i = 1; i < numProc; i++)
        {
            seq2len = strlen(seq2[i]);

            MPI_Send(&seq2len, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(seq2[i], seq2len, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
        seq2len = strlen(seq2[0]);
        seq1Len = strlen(seq1);

        double numOfMut = getNumOfMut(seq2len);
        int *indexToMakeMut = getAllMutIndex(seq2len, numOfMut);
        int *bestMs = getSeqBestMsAndOffsetOnGpu(seq2[0], seq1, seq1Len, seq2len, indexToMakeMut, numOfMut, seq1Len - seq2len + 2, weights, &bestOffset);

        bestMsForAllSeq[0] = bestMs;
        bestOffsetForAllSeq[0] = bestOffset;
    }
    else
    {

        int seq2Len;
        int seq1Len = strlen(seq1);
        MPI_Recv(&seq2Len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        char *mySeq2 = (char *)malloc(sizeof(char) * seq2Len);
        MPI_Recv(mySeq2, seq2Len, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        double numOfMut = getNumOfMut(seq2Len);

        int *indexToMakeMut = getAllMutIndex(seq2Len, numOfMut);
        int bestOffset, findMs = 0;
        ;
        int *bestMs = getSeqBestMsAndOffsetOnGpu(mySeq2, seq1, seq1Len, seq2Len, indexToMakeMut, numOfMut, seq1Len - seq2Len + 2, weights, &bestOffset);
        double bestScore = 0;

        MPI_Send(bestMs, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&bestOffset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    if (my_rank == 0)
    {
        for (int i = 1; i < numProc; i++)
        {
            MPI_Recv(bestMsForAllSeq[i], 2, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&bestOffsetForAllSeq[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
        }

        writhResultOnFile(numProc, bestMsForAllSeq, bestOffsetForAllSeq);

        t2 = MPI_Wtime();
        printf("runing time is %f\n", t2 - t1);
    }

    MPI_Finalize();
    return 0;
}
