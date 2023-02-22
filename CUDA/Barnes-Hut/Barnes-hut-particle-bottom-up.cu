#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

// costanti e variabili host
int maxCells, numberBody, seed, maxTime = 1;
char fileInput[] = "../../Generate/particle.txt";
double *x, *y, *m, *velX, *velY, *forceX, *forceY;
double maxSize=20;
//double maxSize = 6.162025e+070;

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

//funzione di calcolo dei centri di massa
__global__ void calculateCenterMass(int* child,double* xP,double* yP,double* mP,double numBody){
    
    int body=threadIdx.x + blockDim.x * blockIdx.x;
    if(body>=numBody){
        return;
    }
    int cell=child[body];
    while (true){
        double mass=0;
        double mcX=0;
        double mcY=0;

        for(int i=0;i<4;i++){
            int childCell=child[cell-i];
            if(childCell==-1){
                continue;
            }
            double massCell=mP[childCell];
            if(massCell==0){
                return;
            }
            mass+=massCell;
            mcX+=massCell*xP[childCell];
            mcY+=massCell*yP[childCell];
        }
        xP[cell]=mcX/mass;
        yP[cell]=mcY/mass;
        mP[cell]=mass;
                                                                                printf("cell %d, mass %e\n",cell,mass);
        cell=child[cell-4];
        if(cell==-1){
            return;
        }
    }
}

// funzione per la creazione dell'albero
__global__ void createTree(double* x, double* y,double* mass, double up, double down, double left, double right, int *child,int cell,int numBody)
{  
    int body = threadIdx.x + blockDim.x * blockIdx.x;
    // uccido il thread che non deve inserire particelle
    if(body>numBody){
        return;
    }
    int father=cell;
    bool newBody=true;
    bool finish=false;
    int childPath;
    while(!finish){

        //se inserisco una nuova particella
        if(newBody){
            newBody = false;
            childPath = 0;

            //assegno i path ai figli
            if(x[body]<=0.5*(left+right)){
                //+2                            //    _________________________
                childPath +=2;                  //   |          3 |          1 |
                right = 0.5*(left+right);       //   |    (NW)    |    (NE)    |
            } else {                            //   |            |            |
                //+0                            //   |     -+     |     ++     |
                left = 0.5*(left+right);        //   |____________|____________|
            }                                   //   |          2 |          0 |
            if(y[body]>0.5*(up+down)){          //   |     --     |     +-     |
                //+1                            //   |            |            |
                childPath +=1;                  //   |    (SW)    |    (SE)    |
                down = 0.5*(up+down);           //   |____________|____________|
            }else{
                //+0
                up = 0.5*(up+down);
            }  
        }
        cell=child[father-childPath];
        // ciclo fino a che non trovo una foglia e assegno i path
        while(cell >= numBody){
            
            father = cell;
            childPath=0;
            if(x[body]<=0.5*(left+right)){
                //+2
                childPath +=2;
                right = 0.5*(left+right);
            } else {
                //+0
                left = 0.5*(left+right);
            }
            if(y[body]>0.5*(up+down)){
                //+1
                childPath +=1;
                down = 0.5*(up+down);
            }else{
                //+0
                up = 0.5*(up+down);
            }

            cell = child[father - childPath];
        }
                                                                                //printf("cell: %d\n",cell);
        //controllo se la cella è libera
        if (cell != -2){
            int lock=father-childPath;
                                                                                //printf("cell2: %d\n",cell);
            //blocco la cella per lavoraci, utilizzando una funzione atomica
            if(atomicCAS(&child[lock],cell,-2)==cell){
                if(cell == -1){
                                                                                //printf("lock:%d id:%d d %f, u %f, r %f, l %f\n",lock,body,down,up,right,left);
                    //child[body]=lock;
                    child[body]=father;
                    child[lock] = body;
                    finish=true;     
                }else{
                    while(cell>=0 && cell<numBody){

                        //scalo al prossimo indice con cella libera
                        int newCell = atomicAdd(&pPointer,-5);
                        if(newCell-4<numBody){
                            printf("\nNon ho spazio disponibile\n");
                            return;
                        }
                        //assegno ai figli il valore -1, ovvero puntatore a null
                        child[newCell]=-1;
                        child[newCell-1]=-1;
                        child[newCell-2]=-1;
                        child[newCell-3]=-1;
                        child[newCell-4]=father;
                        
                        //inserisco la vecchia particella
                        childPath=0;
                                                                                        //double down2=down,up2=up,left2=left,right2=right;

                        if(x[cell]<=0.5*(left+right)){
                            //+2
                            childPath +=2;
                            //right = 0.5*(left+right);
                        }else{
                            //left = 0.5*(left+right);
                        }
                        if(y[cell]>0.5*(up+down)){
                            //+1
                            childPath +=1;
                            //down = 0.5*(up+down);
                        }else{
                            //up = 0.5*(up+down);
                        }
                                                                                        //printf("move lock:%d id:%d d %f, u %f, r %f, l %f\n",newCell-childPath,cell,down2,up2,right2,left2);
                        //child[cell]=newCell-childPath;

                        child[cell]=newCell;
                        child[newCell-childPath]=cell;

                        //vedo dove inserire una nuova particella
                        childPath=0;
                        father = newCell;
                        if(x[body]<=0.5*(left+right)){
                            //+2
                            childPath +=2;
                            right = 0.5*(left+right);
                        } else {
                            //+0
                            left = 0.5*(left+right);
                        }
                        if(y[body]>0.5*(up+down)){
                            //+1
                            childPath +=1;
                            down = 0.5*(up+down);
                        }else{
                            //+0
                            up = 0.5*(up+down);
                        }
                        
                        cell=child[newCell-childPath];
                        
                        child[newCell-childPath]=-2;
                        

                        __threadfence();
                        child[lock]=newCell;

                        lock= newCell-childPath;

                    }
                                                                                        //printf("lock:%d id:%d d %f, u %f, r %f, l %f\n",lock,body,down,up,right,left);
                    //child[body]=lock;

                    child[body]=father;
                    child[lock]=body;
                    finish=true;
                }
            }
            //__syncthreads();
        }
        cell = child[father - childPath];

    }
                                                                                //printf("%d",cell);
}

// funzione kernel per inizializzare la variabile globale puntatore
__global__ void setPointer(int num)
{
    pPointer = num-6;
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
        printf("particle %d xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n",\
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
    //maxCells = ((numberBody * 2 + 12000) * 4);
    return file;
}

                                                                                                    __global__ void set0(int* child){
                                                                                                        child[pPointer-4]=0;
                                                                                                    }

//funzione grafica per stampare l'albero creato da crateTree()
void printerTree(int* array, int state, int max,int point){
    if(state==0){
        int counter=0;
        printf("(%d) ",point);
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
        printf("\n%d count %d",point,counter);
        printf("\n\n");
        printf("0\n");
        printerTree(array,state+1,max,point);
        printf("1\n");
        printerTree(array,state+1,max,point-1);
        printf("2\n");
        printerTree(array,state+1,max,point-2);
        printf("3\n");
        printerTree(array,state+1,max,point-3);
        return;
    }

    for(int i=0;i<state;i++){
        printf("\t");
    }
    if(array[point]<-1){
        printf("error");
        return;        
    }
    //printf("%d numero",array[point]);
    if(array[point]<max){
        if(array[point]==-1){
            printf("void\n");
        }else{
            printf("%d ",point);
            printf("point: %d\n",array[point]);
        }
        return;
    }
    printf("0\n");
    printerTree(array,state+1,max,array[point]);
    for(int i=0;i<state;i++){
        printf("\t");
    }
    printf("1\n");
    printerTree(array,state+1,max,array[point]-1);
    for(int i=0;i<state;i++){
        printf("\t");
    }
    printf("2\n");
    printerTree(array,state+1,max,array[point]-2);
    for(int i=0;i<state;i++){
        printf("\t");
    }
    printf("3\n");
    printerTree(array,state+1,max,array[point]-3);

}

//funzione di esecuzione dei vari kernell
void compute(int time)
{
    double *xP, *yP, *massP;
    
    double up=maxSize, down=-maxSize, left=-maxSize, right=maxSize,*massR;
    int *child;
    

    //alloco la memoria dei vari parametrio sul device
    gpuErrchk(cudaMalloc((void **)&xP, sizeof(double) * maxCells * 4));
    gpuErrchk(cudaMalloc((void **)&yP, sizeof(double) * maxCells * 4));   
    gpuErrchk(cudaMalloc((void **)&child, sizeof(int) * maxCells * 4));   
    gpuErrchk(cudaMalloc((void **)&massP, sizeof(double) * maxCells * 4));
    // copio array delle posizioni x, y e masse delle particelle
    cudaMemcpy(xP, x, sizeof(double) * numberBody, cudaMemcpyHostToDevice);   
    cudaMemcpy(yP, y, sizeof(double) * numberBody, cudaMemcpyHostToDevice);  
    cudaMemcpy(massP, m, sizeof(double) * numberBody, cudaMemcpyHostToDevice);
                                                                                                int* childH=(int*) malloc( sizeof(int) * maxCells * 4);
    
    gpuErrchk(cudaDeviceSynchronize());
    // eseguo funzioni cuda
    for (int i = 0; i < time; i++)
    {
        // invoco la funzione per settarre la variabile puntatore globale nel device
        setPointer<<<1,1>>>(maxCells);
        // setto array dei figli a -1 (null)
        cudaMemset(&child[ maxCells - 1], -1, sizeof(int));
        cudaMemset(&child[ maxCells - 2], -1, sizeof(int));
        cudaMemset(&child[ maxCells - 3], -1, sizeof(int));
        cudaMemset(&child[ maxCells - 4], -1, sizeof(int));
        cudaMemset(&child[ maxCells - 5], -1, sizeof(int));

        // genero l'albero
        createTree<<<1, 4>>>(xP, yP, massP, up, down, left, right, child, maxCells-1, numberBody);
        cudaDeviceSynchronize();
        // sincronizzo i kernel a fine esecuzione
                                                                                                set0<<<1,1>>>(child);
                                                                                                cudaMemcpy(childH,child,sizeof(int) * maxCells * 4,cudaMemcpyDeviceToHost);
                                                                                                // ritorno l'albero a l'host per la stampa e lo stampo
                                                                                                printerTree(childH,0,numberBody,maxCells-1);
         
        
        // calcolo centri di massa
        calculateCenterMass<<<4,1>>>(child,xP,yP,massP,numberBody);
        cudaDeviceSynchronize();
        
        // calcolo spostamento particelle
        // calculateMove<<<?>>>(?);

        //gpuErrchk(cudaMalloc((void **)&massR, sizeof(double) * maxCells * 4));
        //cudaMemcpy(massR, massP, sizeof(double) * numberBody, cudaMemcpyDeviceToDevice);
        // cudaDeviceSynchronize();

        cudaFree(massP);
        massP=massR;
        
    }
    // libero memoria
                                                                                                free(childH);
    cudaFree(child);
    cudaFree(xP);
    cudaFree(yP);
    cudaFree(massP);
    //printf("memoria liberata sul device \n");
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
    compute(maxTime);
    // stampo i risultati del calcolo
    //printer();
}