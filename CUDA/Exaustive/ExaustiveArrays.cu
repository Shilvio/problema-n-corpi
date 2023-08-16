#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// numero di corpi, seed di generazione, numero di iterazioni globali
int numberBody, seed, maxTime = 3;
char fileInput[] = "./particle.txt";
// costante gravitazione universale
__constant__ double G = 6.67384E-11;

// puntatori ad array di parametri delle particelle
double *x, *y, *m, *velX, *velY, *forceX, *forceY;

// check presenza gpu e assegnazione thread massimi per blocco
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
    cudaGetDeviceProperties(&pr, 0);
    int f = pr.sharedMemPerBlock / (sizeof(double) * 3);
    /* f = massima dim memoria per blocco/grandezza struct particella
    f sono i thread per blocco massimi calcolati */
    if (f > pr.maxThreadsPerBlock)
    {
        f = pr.maxThreadsPerBlock;
    }

    if (pr.maxThreadsPerMultiProcessor % f)
    {
        int h = pr.maxThreadsPerMultiProcessor;
        while (h > f)
        {
            h = h / 2;
        }

        f = h;
    }
    return f;
}

// calcolo della posizione per una particella di interesse per un dato intervallo di tempo e velocità
__device__ void calculatePosition(double xMe, double yMe, double mMe, double velXMe, double velYMe, double forceXMe, double forceYMe, int time, int ID, double *xEnd, double *yEnd, double *velXEnd, double *velYEnd, double *forceXEnd, double *forceYEnd)
{
    xMe += time * velXMe;
    yMe += time * velYMe;
    velXMe += time / mMe * forceXMe;
    velYMe += time / mMe * forceYMe;

    xEnd[ID] = xMe;
    yEnd[ID] = yMe;
    velXEnd[ID] = velXMe;
    velYEnd[ID] = velYMe;
    forceXEnd[ID] = forceXMe;
    forceYEnd[ID] = forceYMe;
}

// calcolo tutte le forze relative ad una particella di interesse in un intervallo
__global__ void calculateTotalForce(double *xStart, double *yStart, double *mStart, double *velXStart, double *velYStart, double *forceXStart, double *forceYStart, double *xEnd, double *yEnd, double *velXEnd, double *velYEnd, double *forceXEnd, double *forceYEnd, int tot)
{

    int sizeMAx = blockDim.x;
    extern __shared__ double temp[];
    int tId = blockIdx.x * blockDim.x + threadIdx.x;
    int ID = threadIdx.x;

    // prendo le variabili della particella dall' array iniziale
    if (tId >= tot)
    {
        return;
    }
    double xMe = xStart[tId];
    double yMe = yStart[tId];
    double mMe = mStart[tId];
    double velXMe = velXStart[tId];
    double velYMe = velYStart[tId];
    double forceXMe = forceXStart[tId];
    double forceYMe = forceYStart[tId];

    // seleziono i blocchi di particelle e itero per esse e calcolo le forze per un unità di tempo
    for (int i = 0; i < (tot / sizeMAx) + 1; i++)
    {
        // carico in temp la mia particella del blocco
        if (i * sizeMAx + ID < tot)
        {
            temp[ID] = xStart[i * sizeMAx + ID];
            temp[sizeMAx + ID] = yStart[i * sizeMAx + ID];
            temp[sizeMAx + sizeMAx + ID] = mStart[i * sizeMAx + ID];
        }
        __syncthreads();
        // calcolo le forze utilizzando le sole particelle presenti in temp
        for (int j = 0; j < sizeMAx; j++)
        {
            // salto il calcolo con la stessa particella
            if (tId == i * sizeMAx + j)
            {
                continue;
            }
            if (i * sizeMAx + j >= tot)
            {
                break;
            }

            // calcolo delle forze
            double xDiff = xMe - temp[j];
            double yDiff = yMe - temp[sizeMAx + j];
            double dist = sqrt(xDiff * xDiff + yDiff * yDiff);
            double cubeDist = dist * dist * dist;
            forceXMe -= ((G * mMe * temp[sizeMAx + sizeMAx + j]) / cubeDist) * xDiff;
            forceYMe -= ((G * mMe * temp[sizeMAx + sizeMAx + j]) / cubeDist) * yDiff;
        }
        __syncthreads();
    }
    // calcolo la posizione per un unità di tempo
    calculatePosition(xMe, yMe, mMe, velXMe, velYMe, forceXMe, forceYMe, 1, tId, xEnd, yEnd, velXEnd, velYEnd, forceXEnd, forceYEnd);
}

// funzione usata per inizializzare il programma e avviare i calcoli per tutte le unità di tempo
void compute(int time)
{
    int thread = statGPU();
    int block = (numberBody / thread) + 1;

    double *xStart, *yStart, *mStart, *velXStart, *velYStart, *forceXStart, *forceYStart;
    double *xEnd, *yEnd, *velXEnd, *velYEnd, *forceXEnd, *forceYEnd;

    // allocazione array di inizio, da cui prendo le particelle
    cudaMalloc((void **)&xStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&yStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&mStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&velXStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&velYStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&forceXStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&forceYStart, sizeof(double) * numberBody);

    // allocazione array dei risultati
    cudaMalloc((void **)&xEnd, sizeof(double) * numberBody);
    cudaMalloc((void **)&yEnd, sizeof(double) * numberBody);
    cudaMalloc((void **)&velXEnd, sizeof(double) * numberBody);
    cudaMalloc((void **)&velYEnd, sizeof(double) * numberBody);
    cudaMalloc((void **)&forceXEnd, sizeof(double) * numberBody);
    cudaMalloc((void **)&forceYEnd, sizeof(double) * numberBody);

    // copio array delle posizioni x e y delle particelle
    cudaMemcpy(xEnd, x, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    cudaMemcpy(yEnd, y, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    cudaMemcpy(mStart, m, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    cudaMemcpy(velXEnd, velX, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    cudaMemcpy(velYEnd, velY, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    cudaMemcpy(forceXEnd, forceX, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
    cudaMemcpy(forceYEnd, forceY, sizeof(double) * numberBody, cudaMemcpyHostToDevice);

    // avvio i calcoli e scambio gli array di inizio con quello di fine
    for (int i = 0; i < time; i++)
    {
        double *temp = xStart;
        xStart = xEnd;
        xEnd = temp;

        temp = yStart;
        yStart = yEnd;
        yEnd = temp;

        temp = velXStart;
        velXStart = velXEnd;
        velXEnd = temp;

        temp = velYStart;
        velYStart = velYEnd;
        velYEnd = temp;

        temp = forceXStart;
        forceXStart = forceXEnd;
        forceXEnd = temp;

        temp = forceYStart;
        forceYStart = forceYEnd;
        forceYEnd = temp;

        calculateTotalForce<<<block, thread, sizeof(double) * thread * 3>>>(xStart, yStart, mStart, velXStart, velYStart, forceXStart, forceYStart, xEnd, yEnd, velXEnd, velYEnd, forceXEnd, forceYEnd, numberBody);
        cudaDeviceSynchronize();
    }

    // copia finale dei risultati
    cudaMemcpy(x, xEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(y, yEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(velX, velXEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(velY, velYEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(forceX, forceXEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(forceY, forceYEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);

    // libero la memoria a fine programma
    cudaFree(xStart);
    cudaFree(yStart);
    cudaFree(mStart);
    cudaFree(velXStart);
    cudaFree(velYStart);
    cudaFree(forceXStart);
    cudaFree(forceYStart);

    cudaFree(xEnd);
    cudaFree(yEnd);
    cudaFree(velXEnd);
    cudaFree(velYEnd);
    cudaFree(forceXEnd);
    cudaFree(forceYEnd);
}

// popolo l'array con le particelle nel file
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

// aprire il file e prendere i primi valori (seed e numero di corpi)
FILE *initial()
{
    // mi apro il file in lettura
    FILE *file = fopen(fileInput, "r");
    // prendo il seed
    fscanf(file, "%d", &seed);
    // prendo il numero di corpi
    fscanf(file, "%d", &numberBody);
    return file;
}

// stamap dei risultati su file
void printerFile()
{
    FILE *solution = fopen("CUDA_Exaustive.txt", "w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution, "%e,%e,%e,%e,%e,%e,%e\n", x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }
    fclose(solution);
}

// printer di debug
void printer()
{
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }
}

int main()
{
    // apro il file dove si trovano tutte le particelle
    // alloco la memoria per l'array che contierrà tutte le particelle (p1)
    // popolo l'array
    getInput(initial());

    // calcolo il movimento delle particelle nel tempo richiesto
    compute(maxTime);

    printerFile();
    // printer(); // per debug

    // libero la memoria
    free(x);
    free(y);
    free(m);
    free(velX);
    free(velY);
    free(forceX);
    free(forceY);
    exit(1);
}