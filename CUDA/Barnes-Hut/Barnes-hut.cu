#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) getchar();
   }
}

char fileInput[] = "../../Generate/particle.txt";
double *x,*y,*m,*velX,*velY,*forceX,*forceY;
double maxSize = 6.162025e+070;
int numberBody, seed, maxTime = 3;

__constant__ double G = 6.67384E-11; // costante gravitazione universale
__constant__ double THETA = 0.5; // thetha per il calcolo delle forze su particell
__device__ int ppointer;


// double maxSize = 100;
// int count = 0;

 //&p1[i].x, &p1[i].y, &p1[i].mass, &p1[i].velX, &p1[i].velY



__device__ int findCell(int x,int y){
    printf("ppointer:%d\n",ppointer);
}

__global__ void createTree(double* xP,double* yP,double* up,double* down,double* left,double* right,int* child){
    
    int id=threadIdx.x+blockDim.x*blockIdx.x;
    int cell=findCell(xP[id],yP[id]);
}

__global__ void setppointer(int num){
    ppointer=num;
}

int statGPU() {
    int numberGPU;
    cudaGetDeviceCount(&numberGPU);
    if(numberGPU<1){
        printf("non sono state rilevate GPU adeguate per esegiure il programma");
        exit(1);
    }

    cudaDeviceProp pr;
    cudaGetDeviceProperties(&pr,0);//thread per blocco 877
    int f = pr.sharedMemPerBlock/sizeof(double); //massima dim memoria per blocco/grandezza struct particella 
    //printf("\n%d\n",f);

    if(pr.maxThreadsPerMultiProcessor%f){

        int h=pr.maxThreadsPerMultiProcessor;

        while (h>f)
        {
            h=h/2;
        }
        
        f=h;
    }
    //printf("\n%d\n",f);
    return f;
}

/*void printerFile(particle *p1)
{
    FILE* solution=fopen("solution.txt","w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution,"%e,%e,%e,%e,%e,%e,%e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
    fclose(solution);
}*/

void printer()
{
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle xPos= %e, yPos= %e, mass= %e\n", x[i], y[i], m[i]);//, forceX, forceY, velX, velY); , forceX= %e, forceY= %e, velX= %e, velY= %e
    }
}

// calcolo il movimento delle particelle nel tempo richiesto
void compute(int time)
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
    //cudaGetLastError
    printf("ciao\n");
    
    gpuErrchk(cudaMalloc((void**)&xP,sizeof(double) * numberBody));
    cudaMalloc((void**)&yP,sizeof(double) * numberBody);
    cudaMalloc((void**)&up,sizeof(double));
    cudaMalloc((void**)&down,sizeof(double));
    cudaMalloc((void**)&left,sizeof(double));
    cudaMalloc((void**)&right,sizeof(double));
    cudaMalloc((void**)&child,sizeof(int)*(numberBody*2+12000)*4);    
    
    cudaMemcpy(xP,x,sizeof(double) * numberBody,cudaMemcpyHostToDevice);
    cudaMemcpy(yP,y,sizeof(double) * numberBody,cudaMemcpyHostToDevice);
    cudaMemset(up,maxSize,sizeof(double));
    cudaMemset(down,-maxSize,sizeof(double));
    cudaMemset(left,-maxSize,sizeof(double));
    cudaMemset(right,maxSize,sizeof(double));
    cudaMemset(&child[((numberBody*2+12000)*4)-1],-1,sizeof(int));

    setppointer<<<1,1>>>(((numberBody*2+12000)*4)-1);
    cudaDeviceSynchronize();
    
    for(int i=0;i<time;i++){
        
        createTree<<<4,1>>>(xP,yP,up,down,left,right,child);
        cudaDeviceSynchronize();
        //calculateCenterMass<<<?>>>(?);
        //cudaDeviceSynchronize();
        //calculateMove<<<?>>>(?);
        //cudaDeviceSynchronize();
                                                                                            //printf("\ncambio\n");
    }
    cudaFree(xP);
    cudaFree(yP);
    cudaFree(up);
    cudaFree(down);
    cudaFree(left);
    cudaFree(right);
    cudaFree(child);

}

// popolo l'array con le particelle nel file
void getInput(FILE *file)
{
    x= (double*) malloc(sizeof(double)*numberBody);
    y= (double*) malloc(sizeof(double)*numberBody);
    m= (double*) malloc(sizeof(double)*numberBody);
    
    velX= (double*) malloc(sizeof(double)*numberBody);
    velY= (double*) malloc(sizeof(double)*numberBody);
    forceX= (double*) malloc(sizeof(double)*numberBody);
    forceY= (double*) malloc(sizeof(double)*numberBody);
    // prendo i dati per tutti i corpi
    for (int i = 0; i < numberBody; i++)
    {   
        // prendo i dati dal file
        fscanf(file, "%lf%lf%lf%lf%lf", &x[i], &y[i], &m[i], &velX[i], &velY[i]);
        // imposto le forze iniziali a zero
        forceX[i]=0;
        forceY[i]=0;
        //printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
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
    //inizializzo array di indirizzi child
    // popolo gli array
    getInput(file);

    // calcolo il movimento delle particelle nel tempo richiesto
    compute(maxTime);
    printf("\n");
    //printer(p1);

    //printerFile(p1);
    fclose(file);
    exit(1);
}