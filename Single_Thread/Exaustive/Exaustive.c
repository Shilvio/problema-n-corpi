#include <stdio.h>
#include <math.h>

int numberN_corpi =5;



struct particle{
    float y;
    float x;
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
}

int main(){
    int T=3;
    getInput();
    //compute(T);
}