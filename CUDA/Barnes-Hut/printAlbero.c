#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

void printerTree(int* array, int state, int max,int point){
    if(state==0){
        printf("0\n");
        printerTree(array,state+1,max,point);
        printf("1\n");
        printerTree(array,state+1,max,point-1);
        printf("2\n");
        printerTree(array,state+1,max,point-2);
        printf("3\n");
        printerTree(array,state+1,max,point-3);
        return;
    }

    for(int i=0;i<state;i++){
        printf("\t");
    }
    if(array[point]<-1){
        printf("error");
        return;        
    }
    //printf("%d numero",array[point]);
    if(array[point]<max){
        if(array[point]==-1){
            printf("void\n");
        }else{
            printf("point: %d\n",array[point]);
        }
        return;
    }
    printf("0\n");
    printerTree(array,state+1,max,array[point]);
    for(int i=0;i<state;i++){
        printf("\t");
    }
    printf("1\n");
    printerTree(array,state+1,max,array[point]-1);
    for(int i=0;i<state;i++){
        printf("\t");
    }
    printf("2\n");
    printerTree(array,state+1,max,array[point]-2);
    for(int i=0;i<state;i++){
        printf("\t");
    }
    printf("3\n");
    printerTree(array,state+1,max,array[point]-3);

}


int main()
{
    int array[]={0,1,2,3,   -1,-1,    -1,3,2,-1,     -1,9,0,-1,     1,-1,13,-1 };
    printerTree(array,0,4,17);
    return 0;
}