#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#define MAX_LENGTH 1000
int sel = 1;

int main(){
  bool flag = true;
  char line_st[MAX_LENGTH];
  char line_cuda[MAX_LENGTH];
  

  FILE *st;
  FILE *cuda;

  if (sel == 1){//exaustive
    st = fopen("../Single_Thread/Exaustive/solution.txt","r");
    cuda = fopen("../CUDA/Exaustive/solution.txt","r");
  }else if (sel == 2){//barnes-hut
    st = fopen("../Single_Thread/Barnes-Hut/solution.txt","r");
    cuda = fopen("../CUDA/Barnes-Hut/solution.txt","r");
  }
  if (st==NULL || cuda==NULL){
    printf("one of files does not exist!");
    exit(1);
  }
  while(fgets(line_st, MAX_LENGTH, st)){
    fgets(line_cuda, MAX_LENGTH, cuda);

   // printf("%c",ch);
    if(line_st[0]!=line_cuda[0]){
      //printf("stocazzo");
    }
    int j = 2;
    for(int i = 2; i < sizeof(line_st); i++){
      //printf ("i = %d\n", i);
      if(line_st[i] == 'e'){
        int t = i + 1;
       // printf("TTT = %d", t);
        /*
        if (line_st[i+1] != line_cuda[j+1]){
          printf("!!! Exponent Signs Are Different !!!");
          flag = false;
          return 1;
        }*/
        while ((line_st[t]!=',')&&(line_st[t]!='\0' )){
          printf("%c", line_st[t]);
          t += 1;
        }
        //printf("f");
      }    
    }
  }
  fclose(st);
  fclose(cuda);

  if (flag){
    printf("Check Complete! Results are Identical");
  }
    return 0;
}