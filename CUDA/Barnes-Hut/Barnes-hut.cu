#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

int numberBody, seed, maxTime = 3;
char fileInput[] = "../../Generate/particle.txt";

__constant__ double G = 6.67384E-11; // costante gravitazione universale
__constant__ double THETA = 0.5; // thetha per il calcolo delle forze su particell



double maxSize = 6.162025e+070;
// double maxSize = 100;
// int count = 0;
double *x,*y,*m;
 //&p1[i].x, &p1[i].y, &p1[i].mass, &p1[i].velX, &p1[i].velY

//struct particella
typedef struct particle
{
    double x;      // posizione x
    double y;      // posizione y
    double mass;   // massa
    double forceX; // forza applicata alla particella sull' asse x
    double forceY; // forza applicata alla particella sull' asse y
    double velX;   // velocità sull' asse x
    double velY;   // velocità sull' asse y
    
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
    /*
    int thread=statGPU();
    int block=(numberBody/thread)+1;
    */
    //int sizeTree=numberBody*2+12000;

    double *xP,*yP,*up,*down,*left,*right;
    int *child;

    // allocazione della memoria a device
    // gpuErrchk(); da aggiungere
    cudaMalloc((void**)&xP,sizeof(particle) * numberBody);
    cudaMalloc((void**)&yP,sizeof(particle) * numberBody);
    cudaMalloc((void**)&up,sizeof(double));
    cudaMalloc((void**)&down,sizeof(double));
    cudaMalloc((void**)&left,sizeof(double));
    cudaMalloc((void**)&right,sizeof(double));
    cudaMalloc((void**)&child,sizeof(int)*(numberBody*2+12000)*4);    
    
    cudaMemcpy(xP,p1,sizeof(double) * numberBody,cudaMemcpyHostToDevice);
    cudaMemcpy(yP,p1,sizeof(double) * numberBody,cudaMemcpyHostToDevice);
    cudaMemset(up,maxSize,sizeof(double));
    cudaMemset(down,-maxSize,sizeof(double));
    cudaMemset(left,-maxSize,sizeof(double));
    cudaMemset(right,maxSize,sizeof(double));
    cudaMemset(&child[((numberBody*2+12000)*4)-1],-1,sizeof(int));

    
    for(int i=0;i<time;i++){

        //createTree<<<1,4>>>(tree,p1Dstart,sizeTree);
        cudaDeviceSynchronize();
        //calculateCenterMass<<<?>>>(?);
        cudaDeviceSynchronize();
        //calculateMove<<<?>>>(?);
        cudaDeviceSynchronize();
                                                                                            //printf("\ncambio\n");
    }
    

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
        x[i]= p1[i].x;
        y[i]= p1[i].y;
        m[i]= p1[i].mass;
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
    // alloco memoria per variabili host
    particle *p1 = (particle*) malloc(sizeof(particle) * numberBody);
    x= (double*) malloc(sizeof(int)*numberBody);
    y= (double*) malloc(sizeof(int)*numberBody);
    m= (double*) malloc(sizeof(int)*numberBody);
    //inizializzo array di indirizzi child
    // popolo gli array
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