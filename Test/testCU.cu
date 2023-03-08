#include <stdio.h>
#include <stdlib.h>
#include <math.h>

__device__ int h=50;

__global__ void fun(){
    double a = -1.1111111111;
    double b = -1.2222222222;
    double max = fmaxf(a,b);
    printf("\n il massimo è :%e ", max);
}

int main(){
    fun<<<1,1>>>();
}