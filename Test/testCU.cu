#include <stdio.h>
#include <stdlib.h>
#include <math.h>

__global__ void fun(){

    __shared__ int h;

    h=0;

    int point=atomicAdd(&h,1);

    printf("io occupo %d\n",point);
}

int main(){
    fun<<<1,65>>>();
}