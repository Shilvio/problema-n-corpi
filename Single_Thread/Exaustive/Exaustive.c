#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int numberBody,seed,maxTime=2;
char fileInput[]="../../Generate/particle.txt";
double const G=6.67384E-11;


typedef struct particle{
    float x;
    float y;
    float mass;
    double forceX;
    double forceY;
    float velX;
    float velY;
}particle;

void calculatePosition(particle* p,int time){
    p->x+=time*p->velX;
    p->y+=time*p->velY;
    p->velX+=time/p->mass*p->forceX;
    p->velY+=time/p->mass*p->forceY;
}

void calculateTotalForce(particle* p1,int j){
    for(int i=0; i<numberBody;i++){
        if (i==j){
            continue;
        }
        float xDiff=p1[j].x-p1[i].x;
        float yDiff=p1[j].y-p1[i].y;
        float dist=sqrt(xDiff*xDiff+yDiff*yDiff);
        float cubeDist=dist*dist*dist;
        p1[j].forceX-=G*(double)p1[j].mass*(double)p1[i].mass/(double)cubeDist*xDiff;
        p1[j].forceY-=G*(double)p1[j].mass*(double)p1[i].mass/(double)cubeDist*yDiff;
    }
}

void stampa(particle* p1){
    for(int i=0;i<numberBody;i++){
        printf("particle %f, %f, %f, %f, %f, %f, %f\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].forceX,p1[i].forceY,p1[i].velX,p1[i].velY);
    }
}

// da finire
void compute(int time,particle* p1){
    for(int i=0; i<time; i++){
        for(int j=0;j<numberBody;j++){
            calculateTotalForce(p1,j);
        }
        for(int j=0;j<numberBody;j++){
            calculatePosition(&p1[j],1);
        }
    }
}

//riempo l'array con le particelle nel file
void getInput(FILE* file,particle* p1){
    //for numberBody
    for(int i=0;i<numberBody;i++){
        fscanf(file,"%f%f%f%f%f",&p1[i].x,&p1[i].y,&p1[i].mass,&p1[i].velX,&p1[i].velY);
        p1[i].forceX=0;
        p1[i].forceY=0;
        printf("particle %.2f, %f, %f, %f, %f, %f, %f\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].forceX,p1[i].forceY,p1[i].velX,p1[i].velY);
    }
    fclose(file);
}

FILE* initial(){
    FILE* file=fopen(fileInput,"r");
    //prendo e stampo il seed
    fscanf(file,"%d",&seed);
    printf("%d\n",seed);
    //prendo e stampo il numero di corpi
    fscanf(file,"%d",&numberBody);
    printf("%d\n",numberBody);
    return file;
}

int main(){
    FILE* file=initial();
    //array che contiene tutte le particelle (p1)
    particle* p1=malloc(sizeof(particle)*numberBody);
    getInput(file,p1);
    /*printf("\n");
                                                                    stampa(p1);
    compute(maxTime,p1);
    printf("\n");
                                                                    stampa(p1);*/
}