#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int const NUMBER_BODY = 4;
float const dim=100;
float const mass=100;
float const vel=200;
unsigned int const seed=0;

int main(){
    FILE* file=fopen("particle.txt","w");
    srand(seed);
    fprintf(file,"%d %d\n",seed,NUMBER_BODY);
    for(int i=0;i<NUMBER_BODY;i++)
        fprintf(file,"%f %f %f %f %f\n",rand()/(RAND_MAX/dim),rand()/(RAND_MAX/dim),\
        rand()/(RAND_MAX/mass),rand()/(RAND_MAX/vel)-vel/2,rand()/(RAND_MAX/vel)-vel/2);
    return 0;
}