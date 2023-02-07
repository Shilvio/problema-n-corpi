#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

// costanti e variabili host
int maxCells, numberBody, seed, maxTime = 1;
char fileInput[] = "../../Generate/particle.txt";
double *x, *y, *m, *velX, *velY, *forceX, *forceY;
double maxSize=100;
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
__global__ void boundingBox(double *xP, double* yP,int numBody,double *up, double *down, double *left, double *right,int* lock){

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int stride = blockDim.x*gridDim.x;
    float xMin = xP[id];
	float xMax = xP[id];
	float yMin = yP[id];
	float yMax = yP[id];

    __shared__ float leftCache[256];
    __shared__ float rightCache[256];
    __shared__ float upCache[256];
    __shared__ float downCache[256];

    int offset = stride;
	while(id + offset < numBody){
		xMin = fminf(xMin, xP[id + offset]);
		xMax = fmaxf(xMax, xP[id + offset]);
		yMin = fminf(yMin, yP[id + offset]);
		yMax = fmaxf(yMax, yP[id + offset]);
		offset += stride;
	}

	leftCache[threadIdx.x] = xMin;
	rightCache[threadIdx.x] = xMax;
	upCache[threadIdx.x] = yMax;
	downCache[threadIdx.x] = yMin;

	__syncthreads();

	// assumes blockDim.x is a power of 2!
	int i = blockDim.x/2;
	while(i != 0){
		if(threadIdx.x < i){
			leftCache[threadIdx.x] = fminf(leftCache[threadIdx.x], leftCache[threadIdx.x + i]);
			rightCache[threadIdx.x] = fmaxf(rightCache[threadIdx.x], rightCache[threadIdx.x + i]);
			upCache[threadIdx.x] = fmaxf(upCache[threadIdx.x], upCache[threadIdx.x + i]);
			downCache[threadIdx.x] = fminf(downCache[threadIdx.x], downCache[threadIdx.x + i]);
		}
		__syncthreads();
		i /= 2;
	}

	if(threadIdx.x == 0){
		while (atomicCAS(lock, 0 ,1) != 0); // lock
		*left = fminf(*left, leftCache[0])-1;
		*right = fmaxf(*right, rightCache[0])+1;
		*up = fmaxf(*up, upCache[0])+1;
		*down = fminf(*down, downCache[0])-1;
		atomicExch(lock, 0); // unlock
                                                                            //printf("Bounding box: up: %e,down: %e,left: %e,right: %e",*up,*down,*left,*right);
	}
}

__global__ void calculateMovement(int* child,double* xP,double* yP,double* mP,int point,int numBody,double* forceX, double* forceY, double* velX, double* velY,double size){
    
    int body = threadIdx.x + blockDim.x * blockIdx.x;
    int id=threadIdx.x;
    int depth=0;
    //double dist = sqrt(pow(xP[body] - t->mc->x, 2) + pow(yP[body] - t->mc->y, 2));
    extern __shared__ int stack[];
    int stakPoint=0;
    for(int i=0;i<4;i++){
        int cell=child[point-i];
        if(cell!=-1){
            stack[blockDim.x*stakPoint+threadIdx.x]=cell;
            stakPoint++;
        }
    }
    depth++;
    while (stakPoint!=0){
        int cell = stack[blockDim.x*stakPoint+threadIdx.x];
        double dist = rsqrt(pow(xP[body] - xP[cell], 2) + pow(yP[body] - yP[cell], 2));
        stakPoint--;
        if(cell<numBody){
            double xDiff = xP[body] - xP[cell];                                // calcolo la distanza tra la particella 1 e la 2
            double yDiff = yP[body] - yP[cell];                                // (il centro di massa del nodo = particella)
            double cubeDist = dist * dist * dist;                              // elevo al cubo la distanza e applico la formula di newton
            forceX[body] -= ((G * mP[body] * mP[cell]) / cubeDist) * xDiff;    // per il calcolo della forza sui 2 assi
            forceY[body] -= ((G * mP[body] * mP[cell]) / cubeDist) * yDiff;
        }else{
            if((size/pow(2,depth))){

            }
        }
    }
    


}

//funzione di calcolo dei centri di massa
__global__ void calculateCenterMass(int* child,double* xP,double* yP,double* mP,int point){

    int body=threadIdx.x + blockDim.x * blockIdx.x;
    int stride= blockDim.x*gridDim.x;
    point-=4;

    //printf("ppointer: %d",point-(4*body));

    for(int i=point-(4*body);i>pPointer;i-=(4*stride)){

        //printf("\n(%d) ",i);
        xP[i]/=mP[i];
        yP[i]/=mP[i];
        printf("(%d)mass: %e x:%e y:%e\n",i,mP[i],xP[i],yP[i]);

    }
}

// funzione per la creazione dell'albero
__global__ void createTree(double* x, double* y,double* mass, double *upP, double *downP, double *leftP, double *rightP, int *child,int cell,int numBody)
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
    
    double up=*upP;
    double down=*downP;
    double left=*leftP;
    double right=*rightP;
                                                                        //printf("Bounding box2: up: %e,down: %e,left: %e,right: %e\n",up,down,left,right);
    while(!finish){

        //se inserisco una nuova particella
        if(newBody){
            newBody = false;
            childPath = 0;

            //assegno i path ai figli
            if(x[body]<=((right-left)/2)+left){
                //+2                            //    _________________________
                childPath +=2;                  //   |          3 |          1 |
                right = ((right-left)/2)+left;  //   |    (NW)    |    (NE)    |
            } else {                            //   |            |            |
                //+0                            //   |     -+     |     ++     |
                left = ((right-left)/2)+left;   //   |____________|____________|
            }                                   //   |          2 |          0 |
            if(y[body]>((up-down)/2)+down){     //   |     --     |     +-     |
                //+1                            //   |            |            |
                childPath +=1;                  //   |    (SW)    |    (SE)    |
                down = ((up-down)/2)+down;      //   |____________|____________|
            }else{
                //+0
                up = ((up-down)/2)+down;
            }  
        }
        cell=child[father-childPath];
        // ciclo fino a che non trovo una foglia e assegno i path
        while(cell >= numBody){
            
            father = cell;
            childPath=0;
            if(x[body]<=((right-left)/2)+left){
                //+2
                childPath +=2;
                right = ((right-left)/2)+left;
            } else {
                //+0
                left = ((right-left)/2)+left;
            }
            if(y[body]>((up-down)/2)+down){
                //+1
                childPath +=1;
                down = ((up-down)/2)+down;
            }else{
                //+0
                up = ((up-down)/2)+down;
            }
                                                                                //printf("lavoro su %d\n",father);
            atomicAdd(&x[father],mass[body]*x[body]);
            atomicAdd(&y[father],mass[body]*y[body]);
            atomicAdd(&mass[father],mass[body]);

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
                        int newCell = atomicAdd(&pPointer,-4);
                        if(newCell-3<numBody){
                            printf("\nNon ho spazio disponibile\n");
                            return;
                        }
                        //assegno ai figli il valore -1, ovvero puntatore a null
                        child[newCell]=-1;
                        child[newCell-1]=-1;
                        child[newCell-2]=-1;
                        child[newCell-3]=-1;
                        
                        //inserisco la vecchia particella
                        childPath=0;
                                                                                        //double down2=down,up2=up,left2=left,right2=right;

                        printf("x cell %e\n",x[cell]);

                        if(x[cell]<=((right-left)/2)+left){
                            //+2
                            childPath +=2;
                            //right2 = ((right-left)/2)+left;
                        }else{
                            //left2 = ((right-left)/2)+left;
                        }
                        if(y[cell]>((up-down)/2)+down){
                            //+1
                            childPath +=1;
                            //down2 = ((up-down)/2)+down;
                        }else{
                            //up2 = ((up-down)/2)+down;
                        }
                                                                                        //printf("move lock:%d id:%d d %f, u %f, r %f, l %f\n",newCell-childPath,cell,down2,up2,right2,left2);
                        //child[cell]=newCell-childPath;

                                                                                        //printf("lavoro su %d\n",newCell);
                        x[newCell]+=mass[cell]*x[cell];
                        y[newCell]+=mass[cell]*y[cell];
                        mass[newCell]+=mass[cell];

                        child[cell]=newCell;
                        child[newCell-childPath]=cell;

                        //vedo dove inserire una nuova particella
                        childPath=0;
                        father = newCell;
                        if(x[body]<=((right-left)/2)+left){
                            //+2
                            childPath +=2;
                            right = ((right-left)/2)+left;
                        } else {
                            //+0
                            left = ((right-left)/2)+left;
                        }
                        if(y[body]>((up-down)/2)+down){
                            //+1
                            childPath +=1;
                            down = ((up-down)/2)+down;
                        }else{
                            //+0
                            up = ((up-down)/2)+down;
                        }
                                                                                        //printf("lavoro su %d\n",father);
                        x[father]+=mass[body]*x[body];
                        y[father]+=mass[body]*y[body];
                        mass[father]+=mass[body];

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
    pPointer = num-5;
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
                                                                                                        child[pPointer]=0;
                                                                                                        child[pPointer-1]=0;
                                                                                                    }

//funzione grafica per stampare l'albero creato da crateTree()
void printerTree(int* array, int state, int max,int point){
    if(state==0){
        int counter=0;
        printf("(%d) ",point);
        for(int i=point;i>=0;i--){
            printf("%d , ",array[i]);
            counter++;
            if(counter%4==0){
                if(array[i-4]==0 && array[i-5]==0){
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
        printf("1\n");
        printerTree(array,state+1,max,point-1);
        printf("0\n");
        printerTree(array,state+1,max,point);
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

    printf("1\n");
    printerTree(array,state+1,max,array[point]-1);
    for(int i=0;i<state;i++){
        printf("\t");
    }
    printf("0\n");
    printerTree(array,state+1,max,array[point]);
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
    double *xR, *yR, *massR;
    
    double *up, *down, *left, *right;
    int *child,*lock;

    double *foceXP,*foceYP,*velXP,*velYP;
    

    //alloco la memoria dei vari parametrio sul device
    gpuErrchk(cudaMalloc((void **)&xP, sizeof(double) * maxCells * 4));
    gpuErrchk(cudaMalloc((void **)&yP, sizeof(double) * maxCells * 4));   
    gpuErrchk(cudaMalloc((void **)&child, sizeof(int) * maxCells * 4));   
    gpuErrchk(cudaMalloc((void **)&massP, sizeof(double) * maxCells * 4));
    gpuErrchk(cudaMalloc((void **)&foceXP, sizeof(double) * numberBody));
    gpuErrchk(cudaMalloc((void **)&foceYP, sizeof(double) * numberBody));
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
    cudaMemcpy(velXP, velX, sizeof(double) * numberBody, cudaMemcpyHostToDevice);   
    cudaMemcpy(velYP, velY, sizeof(double) * numberBody, cudaMemcpyHostToDevice);  
                                                                                                int* childH=(int*) malloc( sizeof(int) * maxCells * 4);
    
    gpuErrchk(cudaDeviceSynchronize());
    // eseguo funzioni cuda
    for (int i = 0; i < time; i++)
    {
        // invoco la funzione per settarre la variabile puntatore globale nel device
        setPointer<<<1,1>>>(maxCells);
        boundingBox<<<1, 4>>>(xP,yP,numberBody,up, down,left, right,lock);
        // setto array dei figli a -1 (null)
        cudaMemset(&child[ maxCells - 1], -1, sizeof(int));
        cudaMemset(&child[ maxCells - 2], -1, sizeof(int));
        cudaMemset(&child[ maxCells - 3], -1, sizeof(int));
        cudaMemset(&child[ maxCells - 4], -1, sizeof(int));
        cudaMemset(&child[ maxCells - 5], -1, sizeof(int));

        cudaDeviceSynchronize();
        // genero l'albero
        createTree<<<1, 4>>>(xP, yP, massP, up, down, left, right, child, maxCells-1, numberBody);
        cudaDeviceSynchronize();
        // sincronizzo i kernel a fine esecuzione
                                                                                                set0<<<1,1>>>(child);
                                                                                                cudaMemcpy(childH,child,sizeof(int) * maxCells * 4,cudaMemcpyDeviceToHost);
                                                                                                // ritorno l'albero a l'host per la stampa e lo stampo
                                                                                                printerTree(childH,0,numberBody,maxCells-1);
         
        
        // calcolo centri di massa
        
        calculateCenterMass<<<2,2>>>(child,xP,yP,massP,maxCells-1);
        cudaDeviceSynchronize();
        
        // calcolo spostamento particelle
        calculateMovement<<<2,2,sizeof(int)*64*2>>>(child,xP,yP,massP,maxCells-1, numberBody,foceXP,foceYP,velXP,velYP,maxSize);

        //gpuErrchk(cudaMalloc((void **)&xR, sizeof(double) * maxCells * 4));

        //gpuErrchk(cudaMalloc((void **)&yR, sizeof(double) * maxCells * 4));

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