#include <stdio.h>
#include <math.h>

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

void getInput(){
    FILE* file=fopen("../../Generate/particle.txt","r");
    int seed;
    fscanf(file,"%d",&seed);
    printf("%d\n",seed);
    int numberBody;
    fscanf(file,"%d",&numberBody);
    printf("%d\n",numberBody);
    char temp[10];
    particle particles[numberBody];
    //for numberBody
    particle p;
    fscanf(file,"%s%f%s%f%s%f%s%f%s%f",&temp,&p.x,&temp,&p.y,&temp,&p.mass,&temp,&p.vel_x,&temp,&p.vel_y);
    particles[0]=p;
    fscanf(file,"%s%f%s%f%s%f%s%f%s%f",&temp,&p.x,&temp,&p.y,&temp,&p.mass,&temp,&p.vel_x,&temp,&p.vel_y);
    particles[1]=p;
    printf("particle %f, %f, %f, %f, %f",particles[0].x,p.y,p.mass,p.vel_x,p.vel_y);
}

int main(){
    int T=3;
    getInput();
    //compute(T);
}