#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

int numberBody, seed, maxTime = 3;
char fileInput[] = "../../Generate/particle.txt";
double const G = 6.67384E-11; // costante gravitazione universale
// double const G=1;
double const THETA = 0.5; // thetha per il calcolo delle forze su particell
double maxSize = 6.162025e+070;
// double maxSize = 100;
// int count = 0;

// struct particella
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

typedef struct massCenter
{
    double x;    // poszione x del centro di massa
    double y;    // posizione y del centro di massa
    double mass; // massa totale al centro di massa
} massCenter;

typedef struct quadTree
{
    char id[20];         // index del nodo utile solo al debug
    double x;            // x del centro dell' albero
    double y;            // y del centro del' albero
    double s;            // dimensione
    particle *p;         // puntatore a una particella
    massCenter *mc;      // centro di massa per quadrante
    bool div;            // check della divisione dell' albero
    struct quadTree *nw; // ramo nord ovest dell' albero guardando la suddivisione del quandrante
    struct quadTree *ne; // ramo nord est dell' albero guardando la suddivisione del quandrante
    struct quadTree *sw; // ramo sud ovest dell' albero guardando la suddivisione del quandrante
    struct quadTree *se; // ramo sud est dell' albero guardando la suddivisione del quandrante
                         //    _________________________
                         //   |          4 |          1 |
                         //   |    (NW)    |    (NE)    |
                         //   |            |            |
                         //   |     -+     |     ++     |
                         //   |____________|____________|
                         //   |          3 |          2 |
                         //   |     --     |     +-     |
                         //   |            |            |
                         //   |    (SW)    |    (SE)    |
                         //   |____________|____________|

} quadTree;

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

        //calculateTotalForce<<<block,thread,sizeof(particle)*thread>>>(p1Dstart,p1Dend,numberBody);
        cudaDeviceSynchronize();
        
                                                                                            //printf("\ncambio\n");
    }

    cudaMemcpy(p1,p1Dend,sizeof(particle) * numberBody,cudaMemcpyDeviceToHost);
    

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