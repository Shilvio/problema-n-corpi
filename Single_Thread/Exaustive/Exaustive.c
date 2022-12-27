#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int numberBody,seed;

typedef struct particle{
    float x;
    float y;
    float mass;
    float vel_x;
    float vel_y;
}particle;


void compute(int T){
    for(int i=0; i<T; i++){
        printf("ciao\n");
    }
}

void getInput(FILE* file,particle* p1){
    //mi creo l'arry che dove poi metterò tutte le particelle
    //for numberBody
    for(int i=0;i<numberBody;i++){
        fscanf(file,"%f%f%f%f%f",&p1[i].x,&p1[i].y,&p1[i].mass,&p1[i].vel_x,&p1[i].vel_y);
        printf("particle %f, %f, %f, %f, %f\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].vel_x,p1[i].vel_y);
    }
    fclose(file);
}

FILE* initial(){
    FILE* file=fopen("../../Generate/particle.txt","r");
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
        printf("particle %f, %f, %f, %f, %f\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].vel_x,p1[i].vel_y);
    }
    //compute(T);
}