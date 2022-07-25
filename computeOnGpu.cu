#include <stdio.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <stdlib.h>
#include "omp.h"
#define NUM_OF_WEIGHTS 4
#define SIZE_OF_CONS 9
#define SIZE_OF_S_CONS 11
__constant__ char conservativeGrope[SIZE_OF_CONS][5] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
__constant__ char semiConservativeGroup[SIZE_OF_S_CONS][7] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};

__device__ int onGroup(char c1, char c2, int grpIND)
{
    int j, count;
    char currLetter;

    if (grpIND == 0)
    {

        for (int i = 0; i < SIZE_OF_CONS; i++)
        {
            j = 0;
            count = 0;
            currLetter = conservativeGrope[i][j];

            while (currLetter != '\0')
            {
                if (c1 == currLetter || c2 == currLetter)
                    count++;
                if (count == 2)
                    return 1;
                j++;
                currLetter = conservativeGrope[i][j];
            }
        }
    }
    else
    {
        for (int i = 0; i < SIZE_OF_S_CONS; i++)
        {
            j = 0;
            count = 0;
            currLetter = semiConservativeGroup[i][j];

            while (currLetter != '\0')
            {
                if (c1 == currLetter || c2 == currLetter)
                    count++;
                if (count == 2)
                    return 1;
                j++;
                currLetter = semiConservativeGroup[i][j];
            }
        }
    }

    return 0;
}
__global__ void getArrBestScorePermut(const char *seq2, int seq2Len, const char *Seq1, int *IndexsToMakeMut, int numOfmut, int numOfAlignment, double *weights, double *BestScorePerMut, int *bestAlignment)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numOfmut)
    {
        BestScorePerMut[i] = 0;

        for (int j = 0; j < numOfAlignment; j++)
        {
            int corectseq1Aligment = 0;
            double currScore = 0;
            for (int k = 0; k < seq2Len; k++)
            {
                if (k == IndexsToMakeMut[i + i] || k == IndexsToMakeMut[i + i + 1])

                    corectseq1Aligment++;

                else
                {
                    if (seq2[k] == Seq1[k + j - corectseq1Aligment])
                        currScore += weights[0];

                    else if (onGroup(seq2[k], Seq1[k + j - corectseq1Aligment], 0) == 1)
                        currScore -= weights[1];

                    else if (onGroup(seq2[k], Seq1[k + j - corectseq1Aligment], 1) == 1)
                        currScore -= weights[2];

                    else
                        currScore -= weights[3];
                }
            }
            if (currScore > BestScorePerMut[i])
            {

                BestScorePerMut[i] = currScore;
                bestAlignment[i] = j;
            }
        }
    }
}

int *getSeqBestMsAndOffsetOnGpu(char *seq2, char *seq1, int seq1Len, int seq2Len, int *allIndexsToMakeMut, double numOfMut, int numOfAlignment, double *weights, int *bestOffset)
{

    cudaError_t err = cudaSuccess;

    char *d_seq1, *d_seq2;
    int *d_allIndexsToMakeMut, *d_bestAlignmentPerMut;
    double *d_weights, *d_BestScorePerMut;

    // Allocate the device input seq2
    err = cudaMalloc((void **)&d_seq2, seq2Len * sizeof(char));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device seq2 (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Allocate the device input mutant indexs
    err = cudaMalloc((void **)&d_allIndexsToMakeMut, (numOfMut * 2) * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device mutant indexs (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Allocate the device input seq1
    err = cudaMalloc((void **)&d_seq1, seq1Len * sizeof(char));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device seq1 (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Allocate the device input best alignment array
    err = cudaMalloc((void **)&d_bestAlignmentPerMut, numOfMut * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device best alignment array (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Allocate the device input weights
    err = cudaMalloc((void **)&d_weights, NUM_OF_WEIGHTS * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device weights array (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Allocate the device input best score per mutant
    err = cudaMalloc((void **)&d_BestScorePerMut, numOfMut * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device best score per mutant array (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Copy the host input seq2  in host memory to the device input vectors in device memory
    err = cudaMemcpy(d_seq2, seq2, seq2Len, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy seq2 from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Copy the host input indexs array  in host memory to the device input vectors in device memory
    err = cudaMemcpy(d_allIndexsToMakeMut, allIndexsToMakeMut, (numOfMut * 2) * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy indexs array from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Copy the host input seq1 in host memory to the device input vectors in device memory
    err = cudaMemcpy(d_seq1, seq1, seq1Len, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy seq1 from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Copy the host input weights in host memory to the device input vectors in device memory
    err = cudaMemcpy(d_weights, weights, NUM_OF_WEIGHTS * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy weights from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    int threadsPerBlock = 1024;
    int blocksPerGrid = (numOfMut + threadsPerBlock - 1) / threadsPerBlock;
    double *h_BestScorePerMut = (double *)malloc(sizeof(double) * numOfMut);
    int *h_bestOffset = (int *)malloc(sizeof(int) * numOfMut);
    getArrBestScorePermut<<<blocksPerGrid, threadsPerBlock>>>(d_seq2, seq2Len, d_seq1, d_allIndexsToMakeMut, numOfMut, numOfAlignment, d_weights, d_BestScorePerMut, d_bestAlignmentPerMut);

    // Copy the device input best score per mutant array in device memory to the host input vectors in host memory
    err = cudaMemcpy(h_BestScorePerMut, d_BestScorePerMut, numOfMut * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy BestScorePerMut array from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Copy the device input best offset per mutant array in device memory to the host input vectors in host memory
    err = cudaMemcpy(h_bestOffset, d_bestAlignmentPerMut, numOfMut * sizeof(int), cudaMemcpyDeviceToHost);

    double bs = 0;
    int indexOfBestMsAndOffset;

#pragma omp parallel for
    for (int i = 0; i < numOfMut; i++)
    {
        if (h_BestScorePerMut[i] > bs)
        {
            bs = h_BestScorePerMut[i];
            indexOfBestMsAndOffset = i;
        }
    }
    *bestOffset = h_bestOffset[indexOfBestMsAndOffset];
    int *MS = (int *)malloc(sizeof(int) * 2);
    MS[0] = allIndexsToMakeMut[indexOfBestMsAndOffset + indexOfBestMsAndOffset] + 1;
    MS[1] = allIndexsToMakeMut[indexOfBestMsAndOffset + indexOfBestMsAndOffset + 1] + 1;
    
    // Free device global memory
    err = cudaFree(d_allIndexsToMakeMut);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device allIndexsToMakeMut  (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaFree(d_bestAlignmentPerMut);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device bestAlignmentPerMut  (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaFree(d_BestScorePerMut);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device BestScorePerMut  (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaFree(d_seq1);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device seq1  (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaFree(d_seq2);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device seq2  (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaFree(d_weights);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device weights  (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Free host memory
    free(h_bestOffset);
    free(h_BestScorePerMut);

    return MS;
}