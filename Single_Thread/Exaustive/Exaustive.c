#include <stdio.h>
#include <math.h>

int numberBody;

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

particle* getInput(){
    FILE* file=fopen("../../Generate/particle.txt","r");
    int seed;
    //prendo e stampo il seed
    fscanf(file,"%d",&seed);
    printf("%d\n",seed);
    ////prendo e stampo il numero di corpi
    fscanf(file,"%d",&numberBody);
    printf("%d\n",numberBody);
    char temp[10];
    particle particles[numberBody];
    //for numberBody
    for(int i=0;i<numberBody;i++){
        fscanf(file,"%f%f%f%f%f",&particles[i].x,&particles[i].y,&particles[i].mass,&particles[i].vel_x,&particles[i].vel_y);
        printf("particle %f, %f, %f, %f, %f\n",particles[i].x,particles[i].y,particles[i].mass,particles[i].vel_x,particles[i].vel_y);
    }
    fclose(file);
    return particles;
}

int main(){
    int T=3;
    particle* p1=getInput();
    printf("\n");
    for(int i=0;i<numberBody;i++){
        printf("particle %f, %f, %f, %f, %f\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].vel_x,p1[i].vel_y);
    }
    //compute(T);
}