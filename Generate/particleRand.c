#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int const NUMBER_BODY = 6;
double const dim=100;
double const mass=5.000e+010;
double const vel=10;
unsigned int const seed=0;

int main(){
    FILE* file=fopen("particle.txt","w");
    srand(seed);
    fprintf(file,"%d %d\n",seed,NUMBER_BODY);
    for(int i=0;i<NUMBER_BODY;i++)
        fprintf(file,"%f %f %f %f %f\n",rand()/(RAND_MAX/dim)-dim/2,rand()/(RAND_MAX/dim)-dim/2,\
        rand()/(RAND_MAX/mass),rand()/(RAND_MAX/vel)-vel/2,rand()/(RAND_MAX/vel)-vel/2);
    fclose(file);
    return 0;
}