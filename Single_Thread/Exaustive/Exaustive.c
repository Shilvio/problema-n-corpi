#include <stdlib.h>
#include <math.h>



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

int main(){
    int T=3;
    compute(T);
}