#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int numberBody, seed, maxTime = 3;
char fileInput[] = "../../Generate/particle.txt";
__constant__ double G = 6.67384E-11;
//double const G = 1;

typedef struct particle
{
    double x;
    double y;
    double mass;
    double forceX;
    double forceY;
    double velX;
    double velY;
} particle;

int statGPU() {
    int numberGPU;
    cudaGetDeviceCount(&numberGPU);
    if(numberGPU<1){
        printf("non sono state rilevate GPU adeguate per esegiure il programma");
        exit(1);
    }

    cudaDeviceProp pr;
    cudaGetDeviceProperties(&pr,0);//thread per blocco 877
    int f = pr.sharedMemPerBlock/sizeof(particle); //massima dim memoria per blocco/grandezza struct particella 
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

__device__ void printerKer(particle *p1, int numberBody)
{
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
}

// calcolo del cambio della posizione per una particella di interesse per un dato intervallo di tempo e velocità
__device__ void calculatePosition(particle p, int time,int ID, particle* p1End)
{
    p.x += time * p.velX;
    p.y += time * p.velY;
    p.velX += time / p.mass * p.forceX;
    p.velY += time / p.mass * p.forceY;

    p1End[ID]=p;
}

// calcolo tutte le forze relative ad una particella di interesse
__global__ void calculateTotalForce(particle *p1Start,particle* p1End,int tot)
{
    int sizeMAx=blockDim.x;
    extern __shared__ particle temp[];
    int tId=blockIdx.x*blockDim.x+threadIdx.x;

    if(tId>=tot){
        return;
    }

    int ID=threadIdx.x;
    particle me=p1Start[tId];

    //printf("%d, %d\n",blockIdx.x,ID);

    for(int i=0;i<(tot/sizeMAx)+1;i++){//da scalare (+1)
        if(i*sizeMAx+ID<tot){
            temp[ID]=p1Start[i*sizeMAx+ID];
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
            double xDiff = me.x - temp[j].x;
            double yDiff = me.y - temp[j].y;
            double dist = sqrt(xDiff * xDiff + yDiff * yDiff);
            double cubeDist = dist * dist * dist;
            me.forceX -= ((G * me.mass * temp[j].mass) / cubeDist) * xDiff;
            me.forceY -= ((G * me.mass * temp[j].mass) / cubeDist) * yDiff;
                                                                                            //printf("%d con %d\n",tId,i*sizeMAx+j);
                                                                                            //printf("\n%e",((G * me.mass * temp[i].mass) / cubeDist) * yDiff);
                                                                                            //printf("\n px=%e py=%e fX=%e fY=%e \n",me.x,me.y,me.forceX,me.forceY);
                                                                                            //printf("px=%e py=%e \n",temp[i].x,temp[i].y);
        }
        __syncthreads();

    }

    calculatePosition(me, 1,tId,p1End);
}

void printerFile(particle *p1)
{
    FILE* solution=fopen("solution.txt","w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution,"%e,%e,%e,%e,%e,%e,%e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
    fclose(solution);
}

void printer(particle *p1)
{
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
}

// calcolo il movimento delle particelle nel tempo richiesto
void compute(int time, particle *p1)
{
    int thread=statGPU();
    int block=(numberBody/thread)+1;
    particle *p1Dstart,*p1Dend;
    cudaMalloc((void**)&p1Dstart,sizeof(particle) * numberBody);
    cudaMalloc((void**)&p1Dend,sizeof(particle) * numberBody);
    cudaMemcpy(p1Dend,p1,sizeof(particle) * numberBody,cudaMemcpyHostToDevice);

    for(int i=0;i<time;i++){
        
        particle* temp=p1Dstart;
        p1Dstart=p1Dend;
        p1Dend=temp;

        calculateTotalForce<<<block,thread,sizeof(particle)*thread>>>(p1Dstart,p1Dend,numberBody);
        cudaDeviceSynchronize();
        
                                                                                            //printf("\ncambio\n");
    }

    cudaMemcpy(p1,p1Dend,sizeof(particle) * numberBody,cudaMemcpyDeviceToHost);
    cudaFree(p1Dstart);
    cudaFree(p1Dend);
    
    

}

// popolo l'array con le particelle nel file
void getInput(FILE *file, particle *p1)
{
    // prendo i dati per tutti i corpi
    for (int i = 0; i < numberBody; i++)
    {
        // prendo i dati dal file
        fscanf(file, "%lf%lf%lf%lf%lf", &p1[i].x, &p1[i].y, &p1[i].mass, &p1[i].velX, &p1[i].velY);
        // imposto le forze iniziali a zero
        p1[i].forceX = 0;
        p1[i].forceY = 0;
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
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
    printf("%d\n", seed);
    // prendo il numero di corpi
    fscanf(file, "%d", &numberBody);
    printf("%d\n", numberBody);
    return file;
}

int main()
{
    // apro il file dove si trovano tutte le particelle
    FILE *file = initial();
    // alloco la memoria per l'array che contierrà tutte le particelle (p1)
    particle *p1 = (particle*) malloc(sizeof(particle) * numberBody);
    // popolo l'array
    getInput(file, p1);

    // calcolo il movimento delle particelle nel tempo richiesto
    compute(maxTime, p1);
    printf("\n");
    printer(p1);

    printerFile(p1);
    fclose(file);
    free(p1);
    exit(1);
}