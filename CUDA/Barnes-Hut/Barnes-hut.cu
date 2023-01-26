#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

// costanti e variabili host
int maxCells, numberBody, seed, maxTime = 1;
char fileInput[] = "../../Generate/particle.txt";
double *x, *y, *m, *velX, *velY, *forceX, *forceY;
double maxSize = 6.162025e+070;

// costanti e variabili gpu
__constant__ double G = 6.67384E-11; // costante gravitazione universale
__constant__ double THETA = 0.5;     // thetha per il calcolo delle forze su particell
__device__ int pPointer;

///////////////////////////////////////////GPU ERRORCHECK///////////////////////////////////////////////////////////////
#define gpuErrchk(ans)                        \
    {                                         \
        gpuAssert((ans), __FILE__, __LINE__); \
    }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort)
            exit(code);
    }
}
__device__ int h = 0;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// funzioni gpu

// funzione kernel per trovare la particella da inserire
__device__ int findCell(int x, int y)
{
    printf("ppointer:%d\n", pPointer);
}
// funzione kernel per creare l'albero
__global__ void createTree(double *xP, double *yP, double *up, double *down, double *left, double *right, int *child)
{

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int cell = findCell(xP[id], yP[id]);
}
// funzione kernel per inizializzare la variabile globale puntatore
__global__ void setPointer(int num)
{
    pPointer = num;
}

// funzioni host

void getInput(FILE *file)
{
    x = (double *)malloc(sizeof(double) * numberBody);
    y = (double *)malloc(sizeof(double) * numberBody);
    m = (double *)malloc(sizeof(double) * numberBody);

    velX = (double *)malloc(sizeof(double) * numberBody);
    velY = (double *)malloc(sizeof(double) * numberBody);
    forceX = (double *)malloc(sizeof(double) * numberBody);
    forceY = (double *)malloc(sizeof(double) * numberBody);
    // prendo i dati per tutti i corpi
    printf("\n");
    for (int i = 0; i < numberBody; i++)
    {
        // prendo i dati dal file

        fscanf(file, "%lf%lf%lf%lf%lf", &x[i], &y[i], &m[i], &velX[i], &velY[i]);

        // imposto le forze iniziali a zero
        forceX[i] = 0;
        forceY[i] = 0;
        printf("particle %d xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", i, x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }
    printf("\n");
    // chiudo il file
    fclose(file);
}

FILE *initial()
{

    // mi apro il file in lettura
    FILE *file = fopen(fileInput, "r");
    // prendo il seed
    fscanf(file, "%d", &seed);
    printf("\n");
    printf("seed: %d\n", seed);
    // prendo il numero di corpi
    fscanf(file, "%d", &numberBody);
    printf("numero particelle: %d\n", numberBody);
    // calcolo max cell offset
    maxCells = ((numberBody * 2 + 12000) * 4);
    return file;
}

void compute(int time)
{
    printf("entro in compute\n");
    double *xP, *yP, *up, *down, *left, *right;
    int *child;

    printf("inizio l'allocazione su device: \n");
    printf("\n");
    cudaMalloc((void **)&xP, sizeof(double) * numberBody);
    printf("malloc 0 funzionante\n");
    gpuErrchk(cudaMalloc((void **)&yP, sizeof(double) * numberBody));
    printf("malloc 1 funzionante\n");
    gpuErrchk(cudaMalloc((void **)&up, sizeof(double)));
    printf("malloc 2 funzionante\n");
    gpuErrchk(cudaMalloc((void **)&down, sizeof(double)));
    printf("malloc 3 funzionante\n");
    gpuErrchk(cudaMalloc((void **)&left, sizeof(double)));
    printf("malloc 4 funzionante\n");
    gpuErrchk(cudaMalloc((void **)&right, sizeof(double)));
    printf("malloc 5 funzionante\n");
    gpuErrchk(cudaMalloc((void **)&child, sizeof(int) * maxCells * 4));
    printf("malloc 6 funzionante\n");
    printf("\n");
    // copio array delle posizioni x e y delle particelle
    cudaMemcpy(xP, x, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    printf("array x particelle copiato \n");
    cudaMemcpy(yP, y, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    printf("array y particelle copiato \n");
    printf("\n");
    // alloco le 4 posizioni iniziali della griglia
    cudaMemset(up, maxSize, sizeof(double));
    cudaMemset(down, -maxSize, sizeof(double));
    cudaMemset(left, -maxSize, sizeof(double));
    cudaMemset(right, maxSize, sizeof(double));
    printf("griglia copiata \n");
    printf("\n");
    // setto array dei figli a -1 (null)
    cudaMemset(&child[maxCells - 1], -1, sizeof(int));
    printf("array childs inizializzato \n");
    printf("\n");
    // invoco la funzione per settarre la variabile puntatore globale nel device
    setPointer<<<1, 1>>>(maxCells - 1);
    printf("puntatore settato \n");
    printf("\n");
    gpuErrchk(cudaDeviceSynchronize());
    printf("sincronizzo kernel");
    printf("\n");

    // eseguo funzioni cuda
    for (int i = 0; i < time; i++)
    {

        // funzione che genera l'albero
        createTree<<<4, 1>>>(xP, yP, up, down, left, right, child);
        cudaDeviceSynchronize();
        printf("albero generato, e kernel sincronizzati \n");

        // calculateCenterMass<<<?>>>(?);
        // cudaDeviceSynchronize();
        // calculateMove<<<?>>>(?);
        // cudaDeviceSynchronize();
    }

    // libero memoria
    cudaFree(up);
    cudaFree(down);
    cudaFree(left);
    cudaFree(right);
    cudaFree(child);
    cudaFree(xP);
    cudaFree(yP);
    printf("memoria liberata sul device \n");
}

// stampa le particelle
void printer()
{
    printf("\n");
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle %d xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", i, x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }

    printf("\n");
}

int main()
{
    // avvio getInput
    getInput(initial());
    // avvio compute
    printf("avvio compute\n");
    compute(maxTime);
    // stampo i risultati del calcolo
    printer();
}