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
__device__ int findCell(int body, double x, double y, double up, double down, double left, double right, int *child,int cell,int numBody)
{   
    int father=cell;
    int childPath = 0;
    //    _________________________
    //   |          3 |          1 |
    //   |    (NW)    |    (NE)    |
    //   |            |            |
    //   |     -+     |     ++     |
    //   |____________|____________|
    //   |          2 |          0 |
    //   |     --     |     +-     |
    //   |            |            |
    //   |    (SW)    |    (SE)    |
    //   |____________|____________|
    //
    if(x<=0.5*(left+right)){
        //+2
        childPath +=2;
        right = 0.5*(left+right);
    } else {
        //+0
        left = 0.5*(left+right);
    }
    if(y>0.5*(up+down)){
        //+1
        childPath +=1;
        down = 0.5*(up+down);
    }else{
        //+0
        up = 0.5*(up+down);
    }  

    cell=child[cell-childPath];

    // ciclo fino a che non trovo una foglia
    while(cell >= numBody){
        
        father = cell;
        childPath=0;

        if(x<=0.5*(left+right)){
            //+2
            childPath +=2;
            right = 0.5*(left+right);
        } else {
            //+0
            left = 0.5*(left+right);
        }
        if(y>0.5*(up+down)){
            //+1
            childPath +=1;
            down = 0.5*(up+down);
        }else{
            //+0
            up = 0.5*(up+down);
        }  

        //Possibbile creazione di centro di massa
        cell = child[cell - childPath];
    }
    
    if (cell != -2){
        
        if(atomicCAS(&child[father],cell,-2)==cell){
            if(cell == -1){
                child[father] = cell;
            }else{
                
            }
        }
    }



}
// funzione kernel per creare l'albero
__global__ void createTree(double *xP, double *yP, double up, double down, double left, double right, int *child,int dimMax,int numBody)
{

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    // uccido il thread che non deve inserire particelle
    if(id>numBody){
        return;
    }
    int cell = findCell(id,xP[id], yP[id], up, down, left, right, child, dimMax, numBody);
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
    double *xP, *yP, up=maxSize, down=-maxSize, left=-maxSize, right=maxSize;
    int *child;

    printf("inizio l'allocazione su device: \n");
    printf("\n");
    cudaMalloc((void **)&xP, sizeof(double) * numberBody);
    printf("malloc 0 funzionante\n");
    gpuErrchk(cudaMalloc((void **)&yP, sizeof(double) * numberBody));
    printf("malloc 1 funzionante\n");
    gpuErrchk(cudaMalloc((void **)&child, sizeof(int) * maxCells * 4));
    printf("malloc 2 funzionante\n");
    printf("\n");
    // copio array delle posizioni x e y delle particelle
    cudaMemcpy(xP, x, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    printf("array x particelle copiato \n");
    cudaMemcpy(yP, y, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    printf("array y particelle copiato \n");
    printf("\n");
    // setto array dei figli a -1 (null)
    cudaMemset(&child[ maxCells - 1], -1, sizeof(int));
    cudaMemset(&child[ maxCells - 2], -1, sizeof(int));
    cudaMemset(&child[ maxCells - 3], -1, sizeof(int));
    cudaMemset(&child[ maxCells - 4], -1, sizeof(int));
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
        createTree<<<4, 1>>>(xP, yP, up, down, left, right, child, maxCells-1, numberBody);
        cudaDeviceSynchronize();
        printf("albero generato, e kernel sincronizzati \n");

        // calculateCenterMass<<<?>>>(?);
        // cudaDeviceSynchronize();
        // calculateMove<<<?>>>(?);
        // cudaDeviceSynchronize();
    }

    // libero memoria
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
// stampa i risultati su solution.txt
void printerFile(){
    FILE* solution=fopen("solution.txt","w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution,"%e,%e,%e,%e,%e,%e,%e\n", x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }
    fclose(solution);
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