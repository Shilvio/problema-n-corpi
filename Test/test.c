#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int numberBody,seed,maxTime=2;
char fileInput[]="../../Generate/particle.txt";
double const G=6.67384E-11;

int main(){
    char c[10];
    FILE* f=fopen("test.txt","r");
    fscanf(f,"%s",&c);
    char *stringEnd;
    double num = strtod(c, &stringEnd); 
    printf("%e",num);
}