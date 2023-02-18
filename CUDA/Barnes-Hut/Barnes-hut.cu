#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

// costanti e variabili host
int maxCells, numberBody, seed, maxTime = 5;       // numero massimo di celle, numero particelle, seed di generazione, tempo massimo di esecuzione
char fileInput[] = "../../Generate/particle.txt";  // file di input contenente le particelle iniziali
double *x, *y, *m, *velX, *velY, *forceX, *forceY; // futuri puntatori agli array dei field delle particelle (posizioni: x,y,  masse, velocità: x,y,  forze: x,y)

// costanti e variabili gpu
__device__ const double G = 6.67384E-11; // costante gravitazione universale
__device__ const double THETA = 0.75;    // theta per il calcolo delle forze su particell
__device__ const int stackSize = 16;     // size dello stack per la gestione della ricorsione e per la pila delle profondità
__device__ const int blockSize = 256;    // dimensione dei bocchi, usata per gestire le memorie shared
__device__ int pPointer;                 // puntatore alla prima cella libera dell' array delle celle
__device__ const int deltaTime = 1;      // delta time

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

// calcolo la bounding box delle particelle, applicando tecniche di riduzione gpu
__global__ void boundingBox(double *xP, double *yP, int numBody, double *up, double *down, double *left, double *right, int *lock)
{

    // id del body gestito da i vari blocchi
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int stride = blockDim.x * gridDim.x;

    // controllo se ho una particella da controllare
    if (id >= numBody)
    {
        return;
    }

    // creiamo due tipi variabili, una locale e un array (shared) per ogniuno dei 4 valori
    float xMin = xP[id];
    float xMax = xP[id];
    float yMin = yP[id];
    float yMax = yP[id];

    __shared__ float leftCache[blockSize];
    __shared__ float rightCache[blockSize];
    __shared__ float upCache[blockSize];
    __shared__ float downCache[blockSize];

    int offset = stride;

    // finche mi trovo tra le particelle, cerco i valori minimi e massimi
    while (id + offset < numBody)
    {
        xMin = fminf(xMin, xP[id + offset]);
        xMax = fmaxf(xMax, xP[id + offset]);
        yMin = fminf(yMin, yP[id + offset]);
        yMax = fmaxf(yMax, yP[id + offset]);
        offset += stride;
    }

    // salvo i valori nella memoria shared relativa al thread
    leftCache[threadIdx.x] = xMin;
    rightCache[threadIdx.x] = xMax;
    upCache[threadIdx.x] = yMax;
    downCache[threadIdx.x] = yMin;

    __syncthreads();

    // applico la riduzione dimezzando ogni volta i thread, ottimizzando l'utilizzo dei warp
    int i = blockDim.x / 2;
    while (i != 0)
    {
        if (threadIdx.x < i)
        {
            // confronto i valori vari con quelli in shared e li sovrascrivo con i rispettivi minimi e massimi
            leftCache[threadIdx.x] = fminf(leftCache[threadIdx.x], leftCache[threadIdx.x + i]);
            rightCache[threadIdx.x] = fmaxf(rightCache[threadIdx.x], rightCache[threadIdx.x + i]);
            upCache[threadIdx.x] = fmaxf(upCache[threadIdx.x], upCache[threadIdx.x + i]);
            downCache[threadIdx.x] = fminf(downCache[threadIdx.x], downCache[threadIdx.x + i]);
        }
        __syncthreads();
        i /= 2;
    }

    // il thread 0, esegue la funzione di master e confronta il risultato con gli altri bloccchi
    if (threadIdx.x == 0)
    {
        // utiliziamo un mutex per accedere alla memoria globale evitando concorrenza
        while (atomicCAS(lock, 0, 1) != 0);

        *left = fminf(*left, leftCache[0]);
        *right = fmaxf(*right, rightCache[0]);
        *up = fmaxf(*up, upCache[0]);
        *down = fminf(*down, downCache[0]);

        atomicExch(lock, 0);
    }
}

// espandiamo la bounding box, per allinearci con i calcoli svolti nell'implementazione single thread
__global__ void boundingBoxExpander(double *up, double *down, double *left, double *right)
{
    *right = *right + 1;
    *left = *left - 1;
    *up = *up + 1;
    *down = *down - 1;
                                                                        //printf("\n\nBounding box: up: %e,down: %e,left: %e,right: %e\n\n",*up,*down,*left,*right);
}

__global__ void calculateMovement(double *xP, double *yP, double *mP,double *forceX, double *forceY, double *velX, double *velY,int numBody){

    int body = threadIdx.x + blockDim.x * blockIdx.x;

    if(body>=numBody){
        return;
    };

    xP[body] += deltaTime * velX[body];
    yP[body] += deltaTime * velY[body];
    velX[body] += deltaTime / mP[body] * forceX[body];
    velY[body] += deltaTime / mP[body] * forceY[body];
}

// funzione che calcola il movimento delle particelle e le rispettive forze
__global__ void calculateForce(int *child, double *xP, double *yP, double *mP, int point, int numBody, double *forceX, double *forceY, double *left, double *right)
{
    

    int body = threadIdx.x + blockDim.x * blockIdx.x;

    if(body>=numBody){
        return;
    };
    
    double size = (*right - *left);
    double forceXb = forceX[body];
    double forceYb = forceY[body];
    double xPb = xP[body];
    double yPb = yP[body];
    double mPb = mP[body];

    // utlizziamo due stack, uno per il calcolo del padre e uno per le profondità
    __shared__ int stack[stackSize * blockSize];
    __shared__ int depths[stackSize * blockSize];
    //il puntatore parte a -1, utilizzato come valore nullo
    int stackPoint = -1;

    // ciclo per quelli che sono sicuro siano i suoi figli e aggiorno i rispettivi stack
    for (int i = 0; i < 4; i++)
    {
        int cell = child[point - i];
        if (cell != -1)
        {   
            stackPoint++;
            stack[blockDim.x * stackPoint + threadIdx.x] = cell;
            depths[blockDim.x * stackPoint + threadIdx.x] = 1;
        }
    }

    // ciclo finchè ho filgi da analizzare, carico il figlio e aggiorno lo stack pointer
    while (stackPoint >= 0)
    {
        int cell = stack[blockDim.x * stackPoint + threadIdx.x];
        int depth = depths[blockDim.x * stackPoint + threadIdx.x];
        stackPoint--;
        double dist = sqrtf(pow(xPb - xP[cell], 2) + pow(yPb - yP[cell], 2));

        // controllo di non star confrontando la particella con se stessa
        if (dist == 0)
        {
            continue;
        }
        // se sto guardando una particella calcolo le forze
        if (cell < numBody)
        {
                                                                                    // printf("size: %e\n",(size/pow(2,depth)));
            double xDiff = xPb - xP[cell];                             // calcolo la distanza tra la particella 1 e la 2
            double yDiff = yPb - yP[cell];                             // (il centro di massa del nodo = particella)
            double cubeDist = dist * dist * dist;                           // elevo al cubo la distanza e applico la formula di newton
            forceXb -= ((G * mPb * mP[cell]) / cubeDist) * xDiff; // per il calcolo della forza sui 2 assi
            forceYb -= ((G * mPb * mP[cell]) / cubeDist) * yDiff;
                                                                                    //printf("dist:%e mass:%e xPcell: %e cell: %d\n",mPb,mP[cell],yP[cell],cell);
                                                                                    // printf("body %d, cell %d\n",body,cell);
        }
        else
        {
            // se va oltre il THETA calcolo approssimo, usiamo solo la x per il calcolo del tetha,  \
                (si potrebbe usare un abs max, per vedere chi è il massimo tra x e y) 
            if (((size / pow(2, depth)) / dist < THETA))
            {   
                double xDiff = xPb - xP[cell];                                 // calcolo la distanza tra la particella 1 e la 2
                double yDiff = yPb - yP[cell];                                 // (il centro di massa del nodo = particella)
                double cubeDist = dist * dist * dist;                          // elevo al cubo la distanza e applico la formula di newton
                forceXb -= ((G * mPb * mP[cell]) / cubeDist) * xDiff;          // per il calcolo della forza sui 2 assi
                forceYb -= ((G * mPb * mP[cell]) / cubeDist) * yDiff;
                                                                                    
                                                                                    //printf("ciao\n");
                                                                                    // printf("body %d, size %d",body,((G * mP[body] * mP[cell]) / cubeDist) * xDiff);
            }
            else
            {
                // aggiungo i figli allo stack 
                for (int i = 0; i < 4; i++)
                {
                    int newCell = child[cell - i];
                    if (newCell != -1)
                    {
                        stackPoint++;
                        stack[blockDim.x * stackPoint + threadIdx.x] = newCell;
                        depths[blockDim.x * stackPoint + threadIdx.x] = depth + 1;
                    }
                }
            }
        }
    }

                                                                                // printf("body %d, X %e, Y %e\n",body,forceX[body],forceY[body]);
    // aggiorno i valori delle particelle relative al delta-time
    forceX[body] = forceXb;
    forceY[body] = forceYb;
}

// funzione di calcolo dei centri di massa
__global__ void calculateCenterMass(int *child, double *xP, double *yP, double *mP, int point)
{
    
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int stride = blockDim.x * gridDim.x;
    point -= 4;
 
    for (int i = point - (4 * id); i > pPointer; i -= (4 * stride))
    {
        xP[i] /= mP[i];
        yP[i] /= mP[i];
    }
}

//resetta gli array per il prossimo timestep
__global__ void resetArray(double *xP, double *yP, double *massP, int point)
{

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int stride = blockDim.x * gridDim.x;
    point -= 4;

    for (int i = point - (4 * id); i > pPointer; i -= (4 * stride))
    {

                                                            // printf("\n(%d) ",i);
        xP[i] = 0;
        yP[i] = 0;
        massP[i] = 0;
                                                            // printf("(%d)mass: %e x:%e y:%e\n",i,mP[i],xP[i],yP[i]);
    }
}

// funzione per la creazione dell'albero
__global__ void createTree(double *x, double *y, double *mass, double *upP, double *downP, double *leftP, double *rightP, int *child, int cell, int numBody)
{
    int body = threadIdx.x + blockDim.x * blockIdx.x;
    // uccido il thread che non deve inserire particelle
    if (body > numBody)
    {
        return;
    }
    int father = cell;
    bool newBody = true;
    bool finish = false;
    int childPath;

    double up = *upP;
    double down = *downP;
    double left = *leftP;
    double right = *rightP;

    while (!finish)
    {

        // se inserisco una nuova particella
        if (newBody)
        {
            newBody = false;
            childPath = 0;

            // assegno i path ai figli
            if (x[body] <= ((right - left) / 2) + left)
            {
                //+2                                    
                childPath += 2;                         
                right = ((right - left) / 2) + left;    //    _________________________
            }                                           //   |          3 |          1 |
            else                                        //   |    (NW)    |    (NE)    |
            {                                           //   |            |            |
                //+0                                    //   |     -+     |     ++     |
                left = ((right - left) / 2) + left;     //   |____________|____________|
            }                                           //   |          2 |          0 |
            if (y[body] > ((up - down) / 2) + down)     //   |     --     |     +-     |
            {                                           //   |            |            |
                //+1                                    //   |    (SW)    |    (SE)    |
                childPath += 1;                         //   |____________|____________|
                down = ((up - down) / 2) + down;        
            }
            else
            {
                //+0
                up = ((up - down) / 2) + down;
            }
        }
        cell = child[father - childPath];
        // ciclo fino a che non trovo una foglia e assegno i path
        while (cell >= numBody)
        {

            father = cell;
            childPath = 0;
            if (x[body] <= ((right - left) / 2) + left)
            {
                //+2
                childPath += 2;
                right = ((right - left) / 2) + left;
            }
            else
            {
                //+0
                left = ((right - left) / 2) + left;
            }
            if (y[body] > ((up - down) / 2) + down)
            {
                //+1
                childPath += 1;
                down = ((up - down) / 2) + down;
            }
            else
            {
                //+0
                up = ((up - down) / 2) + down;
            }
            atomicAdd(&x[father], mass[body] * x[body]);
            atomicAdd(&y[father], mass[body] * y[body]);
            atomicAdd(&mass[father], mass[body]);

            cell = child[father - childPath];
        }
        // controllo se la cella è libera
        if (cell != -2)
        {
            int lock = father - childPath;
            // blocco la cella per lavoraci, utilizzando una funzione atomica
            if (atomicCAS(&child[lock], cell, -2) == cell)
            {
                if (cell == -1)
                {
                    child[body] = father;
                    child[lock] = body;
                    finish = true;
                }
                else
                {
                    while (cell >= 0 && cell < numBody)
                    {

                        // scalo al prossimo indice con cella libera
                        int newCell = atomicAdd(&pPointer, -4);
                        if (newCell - 3 < numBody)
                        {
                            printf("\nNon ho spazio disponibile\n");
                            return;
                        }
                        // assegno ai figli il valore -1, ovvero puntatore a null
                        child[newCell] = -1;
                        child[newCell - 1] = -1;
                        child[newCell - 2] = -1;
                        child[newCell - 3] = -1;

                        // inserisco la vecchia particella
                        childPath = 0;

                        if (x[cell] <= ((right - left) / 2) + left)
                        {
                            //+2
                            childPath += 2;
                        }
                        if (y[cell] > ((up - down) / 2) + down)
                        {
                            //+1
                            childPath += 1;
                        }

                        x[newCell] += mass[cell] * x[cell];
                        y[newCell] += mass[cell] * y[cell];
                        mass[newCell] += mass[cell];

                        child[cell] = newCell;
                        child[newCell - childPath] = cell;

                        // vedo dove inserire una nuova particella
                        childPath = 0;
                        father = newCell;
                        if (x[body] <= ((right - left) / 2) + left)
                        {
                            //+2
                            childPath += 2;
                            right = ((right - left) / 2) + left;
                        }
                        else
                        {
                            //+0
                            left = ((right - left) / 2) + left;
                        }
                        if (y[body] > ((up - down) / 2) + down)
                        {
                            //+1
                            childPath += 1;
                            down = ((up - down) / 2) + down;
                        }
                        else
                        {
                            //+0
                            up = ((up - down) / 2) + down;
                        }

                        x[father] += mass[body] * x[body];
                        y[father] += mass[body] * y[body];
                        mass[father] += mass[body];

                        cell = child[newCell - childPath];

                        child[newCell - childPath] = -2;

                        __threadfence();
                        child[lock] = newCell;

                        lock = newCell - childPath;
                    }

                    child[body] = father;
                    child[lock] = body;
                    finish = true;
                }
            }
        }
        cell = child[father - childPath];
    }
}

// funzione kernel per inizializzare la variabile globale puntatore
__global__ void setPointer(int num)
{
    pPointer = num - 5;
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
        printf("particle %d xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n",
               i, x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
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
    maxCells = ((numberBody * 2 + 50) * 4);
    // maxCells = ((numberBody * 2 + 12000) * 4);
    return file;
}

                                                                                                                                        __global__ void set0(int *child)
                                                                                                                                        {
                                                                                                                                            child[pPointer] = 0;
                                                                                                                                            child[pPointer - 1] = 0;
                                                                                                                                        }

// funzione grafica per stampare l'albero creato da crateTree()
void printerTree(int *array, int state, int max, int point)
{
    if (state == 0)
    {
        int counter = 0;
        printf("(%d) ", point);
        for (int i = point; i >= 0; i--)
        {
            printf("%d , ", array[i]);
            counter++;
            if (counter % 4 == 0)
            {
                if (array[i - 4] == 0 && array[i - 5] == 0)
                {
                    break;
                }
                printf("\n(%d) ", i - 1);
            }
        }
        printf("\n");
        return;
        printf("\n\nPosizione dei body: ");
        int counter2 = max;
        for (int i = max - 1; i >= 0; i--)
        {
            counter2--;
            printf("(%d) %d , ", counter2, array[i]);
        }
        printf("\n%d count %d", point, counter);
        printf("\n\n");
        printf("1\n");
        printerTree(array, state + 1, max, point - 1);
        printf("0\n");
        printerTree(array, state + 1, max, point);
        printf("2\n");
        printerTree(array, state + 1, max, point - 2);
        printf("3\n");
        printerTree(array, state + 1, max, point - 3);
        return;
    }

    for (int i = 0; i < state; i++)
    {
        printf("\t");
    }
    if (array[point] < -1)
    {
        printf("error");
        return;
    }
    if (array[point] < max)
    {
        if (array[point] == -1)
        {
            printf("void\n");
        }
        else
        {
            printf("%d ", point);
            printf("point: %d\n", array[point]);
        }
        return;
    }

    printf("1\n");
    printerTree(array, state + 1, max, array[point] - 1);
    for (int i = 0; i < state; i++)
    {
        printf("\t");
    }
    printf("0\n");
    printerTree(array, state + 1, max, array[point]);
    for (int i = 0; i < state; i++)
    {
        printf("\t");
    }
    printf("2\n");
    printerTree(array, state + 1, max, array[point] - 2);
    for (int i = 0; i < state; i++)
    {
        printf("\t");
    }
    printf("3\n");
    printerTree(array, state + 1, max, array[point] - 3);
}

// riporto i valori da kernel a host
void returnCuda(double *xP, double *yP, double *velXP, double *velYP, double *forceXP, double *forceYP)
{

    gpuErrchk(cudaMemcpy(x, xP, sizeof(double) * numberBody, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(y, yP, sizeof(double) * numberBody, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(velX, velXP, sizeof(double) * numberBody, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(velY, velYP, sizeof(double) * numberBody, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(forceX, forceXP, sizeof(double) * numberBody, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(forceY, forceYP, sizeof(double) * numberBody, cudaMemcpyDeviceToHost));
}

// funzione di esecuzione dei vari kernell
void compute(int time)
{
    double *xP, *yP, *massP;

    double *up, *down, *left, *right;
    int *child, *lock;

    double *forceXP, *forceYP, *velXP, *velYP;

    // alloco la memoria dei vari parametrio sul device
    gpuErrchk(cudaMalloc((void **)&xP, sizeof(double) * maxCells));
    gpuErrchk(cudaMalloc((void **)&yP, sizeof(double) * maxCells));
    gpuErrchk(cudaMalloc((void **)&child, sizeof(int) * maxCells));
    gpuErrchk(cudaMalloc((void **)&massP, sizeof(double) * maxCells));
    gpuErrchk(cudaMalloc((void **)&forceXP, sizeof(double) * numberBody));
    gpuErrchk(cudaMalloc((void **)&forceYP, sizeof(double) * numberBody));
    gpuErrchk(cudaMalloc((void **)&velXP, sizeof(double) * numberBody));
    gpuErrchk(cudaMalloc((void **)&velYP, sizeof(double) * numberBody));
    gpuErrchk(cudaMalloc((void **)&up, sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&down, sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&left, sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&right, sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&lock, sizeof(int)));

    // copio array delle posizioni x, y e masse delle particelle
    gpuErrchk(cudaMemcpy(xP, x, sizeof(double) * numberBody, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(yP, y, sizeof(double) * numberBody, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(massP, m, sizeof(double) * numberBody, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(velXP, velX, sizeof(double) * numberBody, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(velYP, velY, sizeof(double) * numberBody, cudaMemcpyHostToDevice));
                                                                                                                int *childH = (int *)malloc(sizeof(int) * maxCells);
    cudaDeviceSynchronize();
    // eseguo funzioni cuda
    for (int i = 0; i < time; i++)
    {
        // invoco la funzione per settarre la variabile puntatore globale nel device
        setPointer<<<1, 1>>>(maxCells);
        boundingBox<<<1, 4>>>(xP, yP, numberBody, up, down, left, right, lock);     //sizeBlock

        cudaDeviceSynchronize();

        boundingBoxExpander<<<1, 1>>>(up, down, left, right);

        // setto array dei figli a -1 (null)
        gpuErrchk(cudaMemset(&child[maxCells - 1], -1, sizeof(int)));
        gpuErrchk(cudaMemset(&child[maxCells - 2], -1, sizeof(int)));
        gpuErrchk(cudaMemset(&child[maxCells - 3], -1, sizeof(int)));
        gpuErrchk(cudaMemset(&child[maxCells - 4], -1, sizeof(int)));

        cudaDeviceSynchronize();
        // genero l'albero
        createTree<<<1, 4>>>(xP, yP, massP, up, down, left, right, child, maxCells - 1, numberBody); //precisa
        cudaDeviceSynchronize();
        // sincronizzo i kernel a fine esecuzione
                                                                                                                    //set0<<<1, 1>>>(child);
                                                                                                                    cudaMemcpy(childH, child, sizeof(int) * maxCells, cudaMemcpyDeviceToHost);
                                                                                                                    // ritorno l'albero a l'host per la stampa e lo stampo
                                                                                                                    //printerTree(childH, 0, numberBody, maxCells - 1);

        // calcolo centri di massa

        calculateCenterMass<<<2, 2>>>(child, xP, yP, massP, maxCells - 1); 
        cudaDeviceSynchronize();

        // calcolo spostamento particelle
        calculateForce<<<2, 2>>>(child, xP, yP, massP, maxCells - 1, numberBody, forceXP, forceYP, left, right);    //precisa sizeBlock
        cudaDeviceSynchronize();

        calculateMovement<<<2,2>>>(xP, yP, massP,forceXP, forceYP, velXP, velYP, numberBody); //precisa
        cudaDeviceSynchronize();


        resetArray<<<2, 2>>>(xP, yP, massP, maxCells - 1);

        gpuErrchk(cudaMemset(lock, 0, sizeof(int)));

        gpuErrchk(cudaMemset(up, 0, sizeof(double)));
        gpuErrchk(cudaMemset(down, 0, sizeof(double)));
        gpuErrchk(cudaMemset(left, 0, sizeof(double)));
        gpuErrchk(cudaMemset(right, 0, sizeof(double)));

        cudaDeviceSynchronize();
    }

    //ritorno i dati
    returnCuda(xP, yP, velXP, velYP, forceXP, forceYP);

    // libero memoria
    free(childH);
    cudaFree(child);
    cudaFree(xP);
    cudaFree(yP);
    cudaFree(massP);
                                                                // printf("memoria liberata sul device \n");
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
void printerFile()
{
    FILE *solution = fopen("solution.txt", "w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution, "%e,%e,%e,%e,%e,%e,%e\n", x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }
    fclose(solution);
}

int main()
{
    // avvio getInput
    getInput(initial());
    // avvio compute
    compute(maxTime);
    // stampo i risultati del calcolo
    printer();
}