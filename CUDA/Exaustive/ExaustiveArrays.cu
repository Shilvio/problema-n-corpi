#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int numberBody, seed, maxTime = 3;
char fileInput[] = "../../Generate/particle.txt";
__constant__ double G = 6.67384E-11;
double *x, *y, *m, *velX, *velY, *forceX, *forceY;
/*
__global__ void test(double* xStart,double* yStart,double* mStart,double* velXStart,double* velYStart,double* forceXStart,double* forceYStart,double* xEnd,double* yEnd,double* mEnd,double* velXEnd,double* velYEnd,double* forceXEnd,double* forceYEnd,int tot){
    printf("ciao test\n");
}*/

int statGPU() {
    int numberGPU;
    cudaGetDeviceCount(&numberGPU);
    if(numberGPU<1){
        printf("non sono state rilevate GPU adeguate per esegiure il programma");
        exit(1);
    }
    cudaDeviceProp pr;
    cudaGetDeviceProperties(&pr,0);//thread per blocco 877
    int f = pr.sharedMemPerBlock/(sizeof(double)*3); //massima dim memoria per blocco/grandezza struct particella 
    //f sono i thread per blocco massimi
    //printf("\n%d\n",f);
    if(f > pr.maxThreadsPerBlock){
        f=pr.maxThreadsPerBlock;
    }
    if(pr.maxThreadsPerMultiProcessor%f){
        int h=pr.maxThreadsPerMultiProcessor;
        while (h>f)
        {
            h=h/2;
        }
        
        f=h;
    }
    printf("\n%d\n",f);
    return f;
}
/*
__device__ void printerKer(particle *p1, int numberBody)
{
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
}
*/
// calcolo del cambio della posizione per una particella di interesse per un dato intervallo di tempo e velocità
__device__ void calculatePosition(double xMe,double yMe,double mMe,double velXMe,double velYMe,double forceXMe,double forceYMe, int time,int ID, double* xEnd,double* yEnd,double* velXEnd,double* velYEnd,double* forceXEnd,double* forceYEnd)
{
    xMe += time * velXMe;
    yMe += time * velYMe;
    velXMe += time / mMe * forceXMe;
    velYMe += time / mMe * forceYMe;

    xEnd[ID]=xMe;
    yEnd[ID]=yMe;
    velXEnd[ID]=velXMe;
    velYEnd[ID]=velYMe;
    forceXEnd[ID]=forceXMe;
    forceYEnd[ID]=forceYMe;
}

// calcolo tutte le forze relative ad una particella di interesse
__global__ void calculateTotalForce(double* xStart,double* yStart,double* mStart,double* velXStart,double* velYStart,double* forceXStart,double* forceYStart,double* xEnd,double* yEnd,double* velXEnd,double* velYEnd,double* forceXEnd,double* forceYEnd,int tot)
{
    
    int sizeMAx=blockDim.x;
    extern __shared__ double temp[];
    int tId=blockIdx.x*blockDim.x+threadIdx.x;
    int ID=threadIdx.x;

    if(tId>=tot){
        return;
    }
    double xMe=xStart[tId];
    double yMe=yStart[tId];
    double mMe=mStart[tId];
    double velXMe=velXStart[tId];
    double velYMe=velYStart[tId];
    double forceXMe=forceXStart[tId];
    double forceYMe=forceYStart[tId];
    
    for(int i=0;i<(tot/sizeMAx)+1;i++){//da scalare (+1)
        if(i*sizeMAx+ID<tot){
            temp[ID]=xStart[i*sizeMAx+ID];
            temp[sizeMAx+ID]=yStart[i*sizeMAx+ID];
            temp[sizeMAx+sizeMAx+ID]=mStart[i*sizeMAx+ID];
        }
        __syncthreads();
        for (int j = 0; j < sizeMAx; j++)//da scalare
        {
            if (tId == i*sizeMAx+j)
            {
                continue;
            }
            //aggiunto 
            if(i*sizeMAx+j>=tot){
                break;
            }
            // calcolo delle forze
            double xDiff = xMe - temp[j];
            double yDiff = yMe - temp[sizeMAx+j];
            double dist = sqrt(xDiff * xDiff + yDiff * yDiff);
            double cubeDist = dist * dist * dist;
            forceXMe -= ((G * mMe * temp[sizeMAx+sizeMAx+j]) / cubeDist) * xDiff;
            forceYMe -= ((G * mMe * temp[sizeMAx+sizeMAx+j]) / cubeDist) * yDiff;
            //printf("%e\n",mStart[i*sizeMAx+ID]);
                                                                                            //printf("%d con %d\n",tId,i*sizeMAx+j);
                                                                                            //printf("\n%e",((G * me.mass * temp[i].mass) / cubeDist) * yDiff);
                                                                                            //printf("\n px=%e py=%e fX=%e fY=%e \n",me.x,me.y,me.forceX,me.forceY);
                                                                                            //printf("px=%e py=%e \n",temp[i].x,temp[i].y);
        }
        __syncthreads();
    }
    
    calculatePosition(xMe,yMe,mMe,velXMe,velYMe,forceXMe,forceYMe, 1,tId,xEnd,yEnd,velXEnd,velYEnd,forceXEnd,forceYEnd);
}
// calcolo il movimento delle particelle nel tempo richiesto
void compute(int time)
{
    int thread=statGPU();
    int block=(numberBody/thread)+1;
    
    double *xStart, *yStart, *mStart, *velXStart, *velYStart, *forceXStart, *forceYStart;
    double *xEnd, *yEnd, *velXEnd, *velYEnd, *forceXEnd, *forceYEnd;

    cudaMalloc((void **)&xStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&yStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&mStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&velXStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&velYStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&forceXStart, sizeof(double) * numberBody);
    cudaMalloc((void **)&forceYStart, sizeof(double) * numberBody);

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
    // alloco le 4 posizioni iniziali della griglia
    // setto array dei figli a -1 (null)
    for(int i=0;i<time;i++){
        
        double* temp=xStart;
        xStart=xEnd;
        xEnd=temp;

        temp=yStart;
        yStart=yEnd;
        yEnd=temp;

        temp=velXStart;
        velXStart=velXEnd;
        velXEnd=temp;

        temp=velYStart;
        velYStart=velYEnd;
        velYEnd=temp;

        temp=forceXStart;
        forceXStart=forceXEnd;
        forceXEnd=temp;

        temp=forceYStart;
        forceYStart=forceYEnd;
        forceYEnd=temp;

        calculateTotalForce<<<block,thread,sizeof(double)*thread*3>>>(xStart,yStart,mStart,velXStart,velYStart,forceXStart,forceYStart,xEnd,yEnd,velXEnd,velYEnd,forceXEnd,forceYEnd,numberBody);
        cudaDeviceSynchronize();
        
                                                                                            //printf("\ncambio\n");
    }

    cudaMemcpy(x ,xEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(y, yEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(velX, velXEnd,  sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(velY, velYEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(forceX, forceXEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);
    cudaMemcpy(forceY, forceYEnd, sizeof(double) * numberBody, cudaMemcpyDeviceToHost);

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

// aprire il file e prendere i primi valori (seed e numero di corpi)
FILE *initial()
{
    // mi apro il file in lettura
    FILE *file = fopen(fileInput, "r");
    // prendo il seed
    fscanf(file, "%d", &seed);
    printf("%d\n", seed);
    // prendo il numero di corpi
    fscanf(file, "%d", &numberBody);
    printf("%d\n", numberBody);
    return file;
}

void printerFile(){
    FILE* solution=fopen("solutionArray.txt","w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution,"%e,%e,%e,%e,%e,%e,%e\n", x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }
    fclose(solution);
}

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
    printf("\n");
    
    printerFile();
    printer();

    free(x);
    free(y);
    free(m);
    free(velX);
    free(velY);
    free(forceX);
    free(forceY);
    exit(1);
}