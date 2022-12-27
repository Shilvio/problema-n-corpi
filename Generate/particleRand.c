#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int const NUMBER_BODY = 4;
float const dim=100;
float const mass=100;
float const vel=100;
unsigned int const seed=0;

int main(){
    FILE* file=fopen("particle.txt","w");
    srand(seed);
    fprintf(file,"%d\n",seed);
    for(int i=0;i<NUMBER_BODY;i++)
        fprintf(file,"x=%f, y=%f, m=%f, vx=%f, vy=%f\n",(float)rand()/((float)RAND_MAX/dim),(float)rand()/((float)RAND_MAX/dim),\
        (float)rand()/((float)RAND_MAX/mass),(float)rand()/((float)RAND_MAX/vel),(float)rand()/((float)RAND_MAX/vel));
    return 0;
}