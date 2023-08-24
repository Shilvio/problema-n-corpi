#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

// costanti e variabili host
int maxCells, numberBody, seed, maxTime = 5;
char fileInput[] = "./particle.txt";
double *x, *y, *m, *velX, *velY, *forceX, *forceY;
int error_h = 0;

cudaDeviceProp pr;

// costanti e variabili gpu
__device__ const double G = 6.67384E-11; // costante gravitazione universale
__device__ const double THETA = 1.2;     // theta per il calcolo delle forze su particell
__device__ const int blockSize = 256;    // dimensione dei bocchi, usata per gestire le memorie shared
__device__ int pPointer;                 // puntatore alla prima cella libera dell' array delle celle
__device__ const int deltaTime = 1;      // delta time
__device__ int error = 0;

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

// FUNZIONI GPU

// calcolo spotamento particelle
__global__ void calculateMovement(double *xP, double *yP, double *mP, double *forceX, double *forceY, double *velX, double *velY, int numBody)
{

    int body = threadIdx.x + blockDim.x * blockIdx.x;

    if (body >= numBody)
    {
        return;
    };

    xP[body] += deltaTime * velX[body];
    yP[body] += deltaTime * velY[body];
    velX[body] += deltaTime / mP[body] * forceX[body];
    velY[body] += deltaTime / mP[body] * forceY[body];
}

// calcolo delle forze
__global__ void calculateForce(int *child, double *xP, double *yP, double *mP, int point, int numBody, double *forceX, double *forceY, double *left, double *right)
{
    int body = threadIdx.x + blockDim.x * blockIdx.x;

    if (body >= numBody)
    {
        return;
    };

    double size = (*right - *left); // dimensione cella
    double forceXb = forceX[body];  // forze xy della particella
    double forceYb = forceY[body];
    double xPb = xP[body]; // posizioni e massa xy della partiella
    double yPb = yP[body];
    double mPb = mP[body];

    int depth = 0; // profondita albero

    int cell = point; // cella dell' array di analisi / origine
    int pre = cell;   // cella da cui proveniamo

    double dist = sqrtf(pow(xPb - xP[cell], 2) + pow(yPb - yP[cell], 2)); // distanza tra la particella e il centro di massa

    if (dist == 0)
    {
        return;
    }

    if (((size / pow(2, depth)) / dist < THETA))
    {
        double xDiff = xPb - xP[cell];                        // calcolo la distanza tra la particella 1 e la 2
        double yDiff = yPb - yP[cell];                        // (il centro di massa del nodo = particella)
        double cubeDist = dist * dist * dist;                 // elevo al cubo la distanza e applico la formula di newton
        forceXb -= ((G * mPb * mP[cell]) / cubeDist) * xDiff; // per il calcolo della forza sui 2 assi
        forceYb -= ((G * mPb * mP[cell]) / cubeDist) * yDiff;

        forceX[body] = forceXb;
        forceY[body] = forceYb;
        return;
    }

    // ciclo finchè ho filgi da analizzare, carico il figlio e aggiorno lo stack pointer
    while (true)
    {
        // la cella non ha figli esco
        if (cell == -1)
        {
            break;
        }
        // la cella ha figli e si analizza il primo
        if (pre == cell)
        {
            cell = child[pre];

            for (int i = 1; i < 4; i++)
            {
                if (cell == -1)
                {
                    cell = child[pre - i];
                }
                else
                {
                    break;
                }
            }

            // aumento la profondità dell' albero
            depth++;
            pre = cell;

            // se la cella non ha figli o particelle, esco
            if (cell == -1)
            {
                break;
            }
            // trovo il prossimo figlio
        }
        else
        {

            int i;
            // trovo il figlio che ho analizzato
            for (i = 0; i < 4; i++)
            {
                if (child[cell - i] == pre)
                {
                    pre = cell;
                    break;
                }
            }

            // non ha altri figli quindi torno al padre
            if (i == 3)
            {
                depth--;
                cell = child[pre - 4];
                continue;
            }
            // se non sono ancora uscito
            i++;                   // incremento il puntatore per vedere il valore del prossimo valore dell'array
            cell = child[pre - i]; // inserisco in cell il valore
            i++;                   // incremento ancora per controllare

            for (; i < 4; i++)
            {
                // se il quadrante in cui mi trovo è vuoto, pendo il valore successivo dell'array ed esco
                if (cell == -1)
                {
                    cell = child[pre - i];
                }
                else
                {
                    break;
                }
            }

            // se la cella comunque non ha figli torno al padre
            if (cell == -1)
            {
                depth--;
                cell = child[pre - 4]; // -4 e il puntatore al padre nell'array
                continue;
            }
            // scendo e aumento la profondità
            pre = cell;
            depth++;
        }

        double dist = sqrtf(pow(xPb - xP[cell], 2) + pow(yPb - yP[cell], 2));

        // controllo di non star confrontando la particella con se stessa
        if (dist == 0)
        {
            cell = child[pre];
            depth--;
            continue;
        }

        // se sto guardando una particella calcolo le forze
        if (cell < numBody)
        {
            double xDiff = xPb - xP[cell];                        // calcolo la distanza tra la particella 1 e la 2
            double yDiff = yPb - yP[cell];                        // (il centro di massa del nodo = particella)
            double cubeDist = dist * dist * dist;                 // elevo al cubo la distanza e applico la formula di newton
            forceXb -= ((G * mPb * mP[cell]) / cubeDist) * xDiff; // per il calcolo della forza sui 2 assi
            forceYb -= ((G * mPb * mP[cell]) / cubeDist) * yDiff;

            cell = child[pre];
            depth--;
            continue;
        }
        // se va oltre il THETA calcolo approssimo, usiamo solo la x per il calcolo del tetha, (si potrebbe usare un abs max, per vedere chi è il massimo tra x e y)

        if (((size / pow(2, depth)) / dist < THETA))
        {
            double xDiff = xPb - xP[cell];                        // calcolo la distanza tra la particella 1 e la 2
            double yDiff = yPb - yP[cell];                        // (il centro di massa del nodo = particella)
            double cubeDist = dist * dist * dist;                 // elevo al cubo la distanza e applico la formula di newton
            forceXb -= ((G * mPb * mP[cell]) / cubeDist) * xDiff; // per il calcolo della forza sui 2 assi
            forceYb -= ((G * mPb * mP[cell]) / cubeDist) * yDiff;

            cell = child[pre - 4];
            depth--;
        }
    }

    // aggiorno i valori delle particelle relative al delta-time
    forceX[body] = forceXb;
    forceY[body] = forceYb;
}

// setto x y iniziali della bounding box
__global__ void initialPosition(double *up, double *down, double *left, double *right, double *x, double *y)
{

    *up = y[0];
    *down = y[0];
    *left = x[0];
    *right = x[0];
}

// bounding box delle particelle, applicando tecniche di riduzione gpu
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
            if (id + i < numBody)
            {
                // confronto i valori vari con quelli in shared e li sovrascrivo con i rispettivi minimi e massimi
                leftCache[threadIdx.x] = fminf(leftCache[threadIdx.x], leftCache[threadIdx.x + i]);
                rightCache[threadIdx.x] = fmaxf(rightCache[threadIdx.x], rightCache[threadIdx.x + i]);
                upCache[threadIdx.x] = fmaxf(upCache[threadIdx.x], upCache[threadIdx.x + i]);
                downCache[threadIdx.x] = fminf(downCache[threadIdx.x], downCache[threadIdx.x + i]);
            }
        }
        __syncthreads();
        i /= 2;
    }

    // il thread 0, esegue la funzione di master e confronta il risultato con gli altri bloccchi
    if (threadIdx.x == 0)
    {
        // utiliziamo un mutex per accedere alla memoria globale evitando concorrenza
        while (atomicCAS(lock, 0, 1) != 0)
            ;
        __threadfence();
        *left = fminf(*left, leftCache[0]);
        *right = fmaxf(*right, rightCache[0]);
        *up = fmaxf(*up, upCache[0]);
        *down = fminf(*down, downCache[0]);
        atomicExch(lock, 0);
    }
}

// espandiamo la bounding box,
// per allinearci con i calcoli svolti nell'implementazione single thread
__global__ void boundingBoxExpander(double *up, double *down, double *left, double *right)
{

    *right = *right + 1;
    *left = *left - 1;
    *up = *up + 1;
    *down = *down - 1;
}

// calcolo dei centri di massa
__global__ void calculateCenterMass(int *child, double *xP, double *yP, double *mP, double numBody)
{

    int body = threadIdx.x + blockDim.x * blockIdx.x;
    if (body >= numBody)
    {
        return;
    }
    int cell = child[body];
    while (true)
    {
        double mass = 0;
        double mcX = 0;
        double mcY = 0;
        // itero per i 4 quadranti componenti del quadrante in questione
        for (int i = 0; i < 4; i++)
        {
            int childCell = child[cell - i];
            //  controllo che non sia vuoto
            if (childCell == -1)
            {
                continue;
            }
            double massCell = mP[childCell]; // centro di massa della cella figlio
            // se non è stata calcolata la massa allora un altro thread ci sta lavorando su
            if (massCell == 0)
            {
                return;
            }
            // aggiorno i valori
            mass += massCell;
            mcX += massCell * xP[childCell];
            mcY += massCell * yP[childCell];
        }
        // calcolo del centro di massa
        xP[cell] = mcX / mass;
        yP[cell] = mcY / mass;
        mP[cell] = mass;

        // controllo il padre
        cell = child[cell - 4];
        if (cell == -1)
        {
            return;
        }
    }
}

// creazione dell'albero
__global__ void createTree(double *x, double *y, double *mass, double *upA, double *downA, double *leftA, double *rightA, int *child, int cell, int numBody)
{
    double up = *upA, down = *downA, left = *leftA, right = *rightA;
    int body = threadIdx.x + blockDim.x * blockIdx.x;

    // uccido il thread che non deve inserire particelle
    if (body >= numBody)
    {
        return;
    }

    int father = cell;
    bool newBody = true;
    bool finish = false;
    int childPath;

    while (!finish)
    {

        // se inserisco una nuova particella
        if (newBody)
        {
            newBody = false;
            childPath = 0;

            // assegno i path ai figli
            if (x[body] <= 0.5 * (left + right))
            {
                //+2
                childPath += 2;
                right = 0.5 * (left + right);
            }
            else
            {
                //+0
                left = 0.5 * (left + right);
            }
            if (y[body] > 0.5 * (up + down))
            {
                //+1
                childPath += 1;
                down = 0.5 * (up + down);
            }
            else
            {
                //+0
                up = 0.5 * (up + down);
            }
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
        }

        // quadrante di inizio
        cell = child[father - childPath];

        // ciclo fino a che non trovo una foglia e assegno i path
        while (cell >= numBody)
        {
            father = cell;
            childPath = 0;
            if (x[body] <= 0.5 * (left + right))
            {
                //+2
                childPath += 2;
                right = 0.5 * (left + right);
            }
            else
            {
                //+0
                left = 0.5 * (left + right);
            }
            if (y[body] > 0.5 * (up + down))
            {
                //+1
                childPath += 1;
                down = 0.5 * (up + down);
            }
            else
            {
                //+0
                up = 0.5 * (up + down);
            }

            // quadrante scelto
            cell = child[father - childPath];
        }

        // controllo se nessun thread sta operando sulla cella
        if (cell != -2)
        {
            int lock = father - childPath;
            // blocco la cella per lavoraci, utilizzando una funzione atomica
            if (atomicCAS(&child[lock], cell, -2) == cell)
            {
                // se il quadrante è vuoto
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
                        int newCell = atomicAdd(&pPointer, -5);
                        // controllo se ho spazio disponibile nell'array per operare
                        if (newCell - 5 < numBody)
                        {
                            error = 1;
                            child[lock] = -1;
                            return;
                        }
                        // assegno ai figli il valore -1, ovvero puntatore a null e salvo il puntatore al padre
                        child[newCell] = -1;
                        child[newCell - 1] = -1;
                        child[newCell - 2] = -1;
                        child[newCell - 3] = -1;
                        child[newCell - 4] = father;

                        // inserisco la vecchia particella
                        childPath = 0;

                        if (x[cell] <= 0.5 * (left + right))
                        {
                            //+2
                            childPath += 2;
                        }
                        if (y[cell] > 0.5 * (up + down))
                        {
                            //+1
                            childPath += 1;
                        }

                        child[cell] = newCell;             // aggiorno il quadrante della vecchia particella
                        child[newCell - childPath] = cell; // inserisco la vecchia particella nel suo nuovo quadrante

                        // vedo dove inserire una nuova particella
                        childPath = 0;
                        father = newCell;
                        if (x[body] <= 0.5 * (left + right))
                        {
                            //+2
                            childPath += 2;
                            right = 0.5 * (left + right);
                        }
                        else
                        {
                            //+0
                            left = 0.5 * (left + right);
                        }

                        if (y[body] > 0.5 * (up + down))
                        {
                            //+1
                            childPath += 1;

                            down = 0.5 * (up + down);
                        }
                        else
                        {
                            //+0
                            up = 0.5 * (up + down);
                        }

                        cell = child[newCell - childPath]; // inserisco il contenuto del quadrante
                        child[newCell - childPath] = -2;   // blocco il quadrante
                        __threadfence();
                        child[lock] = newCell; // inserendo il nuovo quadrante libero il padre

                        lock = newCell - childPath; // aggiorno il valore del lock

                        // se in cell c'è la vecchia particella ricomincio
                    }

                    // aggiorno i valori della nuova particella
                    child[body] = father;
                    child[lock] = body;
                    finish = true;
                }
            }
            //__syncthreads();
        }
        cell = child[father - childPath];
    }
}

// funzione kernel per inizializzare la variabile globale puntatore al primo punto libero dell'array child
__global__ void setPointer(int num)
{
    pPointer = num - 6;
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
    for (int i = 0; i < numberBody; i++)
    {
        // prendo i dati dal file
        fscanf(file, "%lf%lf%lf%lf%lf", &x[i], &y[i], &m[i], &velX[i], &velY[i]);

        // imposto le forze iniziali a zero
        forceX[i] = 0;
        forceY[i] = 0;
    }
    // chiudo il file
    fclose(file);
}

FILE *initial()
{
    // mi apro il file in lettura
    FILE *file = fopen(fileInput, "r");
    // prendo il seed
    fscanf(file, "%d", &seed);
    // prendo il numero di corpi
    fscanf(file, "%d", &numberBody);
    // calcolo max cell offset
    maxCells = ((numberBody * 50 + 12000) * 5);
    return file;
}

__global__ void checkError(int *er)
{
    *er = error;
}

__global__ void resetArray(double *xP, double *yP, double *massP, int point)
{
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int stride = blockDim.x * gridDim.x;
    point -= 4;

    for (int i = point - (4 * id); i > pPointer; i -= (4 * stride))
    {
        xP[i] = 0;
        yP[i] = 0;
        massP[i] = 0;
    }
}

// stampa l'albero a terminale creato da crateTree()
void printerTree(int *array, int state, int max, int point, double up, double down, double left, double right, double *x, double *y)
{
    if (state == 0)
    {
        int counter = 0;
        /*printf("(%d) ",point);
        for(int i=point;i>=0;i--){
            printf("%d , ",array[i]);
            counter++;
            if(counter%5==0){
                if(array[i-5]==0){
                    break;
                }
                printf("\n(%d) ",i-1);
            }
        }
        printf("\n\nPosizione dei body: ");
        int counter2=max;
        for(int i=max-1;i>=0;i--){
            counter2--;
            printf("(%d) %d , ",counter2,array[i]);
        }
        printf("\n");*/

        freopen("output.txt", "a+", stdout);
        // return;
        printf("\n%d count %d", point, counter);
        printf("\n\n");
        printf("0\n");
        printerTree(array, state + 1, max, point, (up + down) * 0.5, down, (left + right) * 0.5, right, x, y);
        printf("1\n");
        printerTree(array, state + 1, max, point - 1, up, (up + down) * 0.5, (left + right) * 0.5, right, x, y);
        printf("2\n");
        printerTree(array, state + 1, max, point - 2, (up + down) * 0.5, down, left, (left + right) * 0.5, x, y);
        printf("3\n");
        printerTree(array, state + 1, max, point - 3, up, (up + down) * 0.5, left, (left + right) * 0.5, x, y);
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

            printf("body: %d x:%e y:%e\n", array[point], x[array[point]], y[array[point]]);
        }
        return;
    }
    printf("0 up: %e down: %e left %e right %e\n\n", up, down, left, right);
    printerTree(array, state + 1, max, array[point], (up + down) * 0.5, down, (left + right) * 0.5, right, x, y);
    for (int i = 0; i < state; i++)
    {
        printf("\t");
    }
    printf("1 up: %e down: %e left %e right %e\n\n", up, down, left, right);
    printerTree(array, state + 1, max, array[point] - 1, up, (up + down) * 0.5, (left + right) * 0.5, right, x, y);
    for (int i = 0; i < state; i++)
    {
        printf("\t");
    }
    printf("2 up: %e down: %e left %e right %e\n\n", up, down, left, right);
    printerTree(array, state + 1, max, array[point] - 2, (up + down) * 0.5, down, left, (left + right) * 0.5, x, y);
    for (int i = 0; i < state; i++)
    {
        printf("\t");
    }
    printf("3 up: %e down: %e left %e right %e\n\n", up, down, left, right);
    printerTree(array, state + 1, max, array[point] - 3, up, (up + down) * 0.5, left, (left + right) * 0.5, x, y);
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

// stampa delle particelle
void printer()
{
    printf("\n");
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle %d xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", i, x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }

    printf("\n");
}

// funzione di esecuzione dei vari kernell
void compute(int time)
{
    int *er;
    cudaMalloc((void **)&er, sizeof(int));
    double *xP, *yP, *massP;
    double *up, *down, *left, *right;
    int *child, *lock;
    double *forceXP, *forceYP, *velXP, *velYP;

    // variabili di ottimizzazione GPU
    int boundingNumBlocks = (numberBody / blockSize) + 1;
    int preciseNumThread = pr.maxThreadsPerBlock;
    int preciseNumBlocks;
    int preciseNumBlockSize = (numberBody / blockSize) + 1;

    if (numberBody < preciseNumThread)
    {
        preciseNumBlocks = 1;
        preciseNumThread = numberBody;
    }
    else
    {
        if (pr.maxThreadsPerMultiProcessor % preciseNumThread != 0)
        {
            preciseNumThread = pr.maxThreadsPerMultiProcessor / 2;
        }
        preciseNumBlocks = (numberBody / preciseNumThread) + 1;
    }

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
    cudaMemcpy(xP, x, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    cudaMemcpy(yP, y, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    cudaMemcpy(massP, m, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    gpuErrchk(cudaMemcpy(velXP, velX, sizeof(double) * numberBody, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(velYP, velY, sizeof(double) * numberBody, cudaMemcpyHostToDevice));

    gpuErrchk(cudaDeviceSynchronize());
    // eseguo funzioni cuda
    for (int i = 0; i < time; i++)
    {
        // invoco la funzione per settarre la variabile puntatore globale nel device
        setPointer<<<1, 1>>>(maxCells);

        initialPosition<<<1, 1>>>(up, down, left, right, xP, yP);
        cudaDeviceSynchronize();

        boundingBox<<<boundingNumBlocks, blockSize>>>(xP, yP, numberBody, up, down, left, right, lock);
        cudaDeviceSynchronize();

        boundingBoxExpander<<<1, 1>>>(up, down, left, right);

        // setto array dei figli a -1 (null)
        cudaMemset(&child[maxCells - 1], -1, sizeof(int));
        cudaMemset(&child[maxCells - 2], -1, sizeof(int));
        cudaMemset(&child[maxCells - 3], -1, sizeof(int));
        cudaMemset(&child[maxCells - 4], -1, sizeof(int));
        cudaMemset(&child[maxCells - 5], -1, sizeof(int));

        // genero l'albero
        createTree<<<preciseNumBlocks, preciseNumThread>>>(xP, yP, massP, up, down, left, right, child, maxCells - 1, numberBody);
        cudaDeviceSynchronize();
        // sincronizzo i kernel a fine esecuzione

        checkError<<<1, 1>>>(er);
        cudaDeviceSynchronize();
        gpuErrchk(cudaMemcpy(&error_h, er, sizeof(int), cudaMemcpyDeviceToHost));
        if (error_h != 0)
        {
            printf("\nNon ho spazio disponibile per creare l'albero\n");
            break;
        }

        // calcolo centri di massa
        calculateCenterMass<<<preciseNumBlocks, preciseNumThread>>>(child, xP, yP, massP, numberBody);
        cudaDeviceSynchronize();

        calculateForce<<<preciseNumBlockSize, blockSize>>>(child, xP, yP, massP, maxCells - 1, numberBody, forceXP, forceYP, left, right);
        cudaDeviceSynchronize();

        // calcolo spostamento particelle
        calculateMovement<<<preciseNumBlocks, preciseNumThread>>>(xP, yP, massP, forceXP, forceYP, velXP, velYP, numberBody);
        cudaDeviceSynchronize();

        resetArray<<<preciseNumBlocks, preciseNumThread>>>(xP, yP, massP, maxCells - 1);
        gpuErrchk(cudaMemset(lock, 0, sizeof(int)));
        cudaDeviceSynchronize();
    }
    returnCuda(xP, yP, velXP, velYP, forceXP, forceYP);
    // libero memoria
    cudaFree(child);
    cudaFree(xP);
    cudaFree(yP);
    cudaFree(massP);
}

// stampa i risultati su file
void printerFile()
{
    FILE *solution = fopen("CUDA_Barnes-hut.txt", "w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution, "%e,%e,%e,%e,%e,%e,%e\n", x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }
    fclose(solution);
}

int statGPU()
{
    int numberGPU;
    cudaGetDeviceCount(&numberGPU);
    if (numberGPU < 1)
    {
        printf("Non sono state rilevate GPU Cuda per esegiure il programma");
        exit(1);
    }
    cudaGetDeviceProperties(&pr, 0);
}

int main()
{

    // verifico le stat della gpu e la sua presenza
    statGPU();
    // avvio getInput
    getInput(initial());
    // avvio compute
    compute(maxTime);
    // stampo i risultati del calcolo
    if (error_h != 0)
    {
        printf("\nNon completabile\n");
        return 0;
    }
    // printer();
    printerFile();
    return 0;
}