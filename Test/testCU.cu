#include <stdio.h>
#include <stdlib.h>
#include <math.h>

__device__ int h=50;

__global__ void fun(){

    int point=atomicAdd(&h,-1);

    printf("io occupo %d\n",point);
}

int main(){
    fun<<<3,2>>>();
}