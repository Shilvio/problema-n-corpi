#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int numberBody, seed, maxTime = 3;
char fileInput[] = "../../Generate/particle.txt";
__constant__ double G = 6.67384E-11;
// double const G = 1;
int statGPU()
{
    int numberGPU;
    cudaGetDeviceCount(&numberGPU);
    if (numberGPU < 1)
    {
        printf("non sono state rilevate GPU adeguate per esegiure il programma");
        exit(1);
    } 

    cudaDeviceProp pr;
    cudaGetDeviceProperties(&pr, 0); // thread per blocco 877

    printf("gpu stats : \n name: %s, \n global memory size: %zu,\n shared memory per block: %zu,\n registers per block: %d,\n warp size: %d,\n memory pitch: %zu,\n max thread per block: %d,\n max thread dimention: %d,\n max grid size: %d,\n total constant memory: %zd,\n clock rate: %d,\n multiprocessor count: %d",
           pr.name, pr.totalGlobalMem, pr.sharedMemPerBlock, pr.regsPerBlock, pr.warpSize, pr.memPitch, pr.maxThreadsPerBlock,
           *pr.maxThreadsDim, *pr.maxGridSize, pr.totalConstMem, pr.clockRate, pr.multiProcessorCount);
    // massima dim memoria per blocco/grandezza struct particella
    // printf("\n%d\n",f);
}

int main()
{
    // apro il file dove si trovano tutte le particelle
    statGPU();
    exit(1);
}