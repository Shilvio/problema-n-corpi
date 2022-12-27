#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int numberBody,seed,maxTime;
char fileInput[]="../../Generate/particle.txt";
float const G=(float)6.67384E-11;;

typedef struct particle{
    float x;
    float y;
    float mass;
    float forceX;
    float forceY;
}particle;

/**/
void calculateTotalForce(particle p,particle* p1,int j){
    for(int i=0; i<numberBody;i++){
        if (i==j){
            continue;
        }
        float xDiff=p1[j].x-p1[i].x;
        float yDiff=p1[j].y-p1[i].y;
        float dist=sqrt(xDiff*xDiff+yDiff*yDiff);
        float cubeDist=dist*dist*dist;
        p1[j].forceX-=G*p1[j].mass*p1[i].mass/cubeDist*xDiff;
        p1[j].forceY-=G*p1[j].mass*p1[i].mass/cubeDist*yDiff;
    }
}

// da finire
void compute(int time,particle* p1){
    for(int i=0; i<time; i++){
        for(int j=0;j<numberBody;j++){
            calculateTotalForce(p1[j],p1,j);
        }
        for(int j=0;j<numberBody;j++){
            calculatePositionVel(p1[j]);
        }
    }
}

void getInput(FILE* file,particle* p1){
    //mi creo l'arry che dove poi metterò tutte le particelle
    //for numberBody
    for(int i=0;i<numberBody;i++){
        fscanf(file,"%f%f%f%f%f",&p1[i].x,&p1[i].y,&p1[i].mass,&p1[i].forceX,&p1[i].forceY);
        printf("particle %f, %f, %f, %f, %f\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].forceX,p1[i].forceY);
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
    particle* p1=malloc(sizeof(particle)*numberBody);
    getInput(file,p1);
    printf("\n");
    for(int i=0;i<numberBody;i++){
        printf("particle %f, %f, %f, %f, %f\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].forceX,p1[i].forceY);
    }
    //compute(maxTime,p1);
}