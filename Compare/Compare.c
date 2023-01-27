/*
funziona con i file scritti dalla seguente funzione:

void printerFile(particle *p1)
{
    FILE* solution=fopen("solution.txt","w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution,"%e,%e,%e,%e,%e,%e,%e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
    fclose(solution);
}
void printerFile()
{
    FILE* solution=fopen("solution.txt","w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution,"%e,%e,%e,%e,%e,%e,%e\n", x[i], y[i], m[i], forceX[i], forceY[i], velX[i], velY[i]);
    }
    fclose(solution);
}
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
//massima lunghezza stringa
#define MAX_LENGTH 1000 
//massimo numero di iterazioni sulle righe
#define MAX_ITERATIONS 100 
//flag di verifica
bool flag = true;
//array per contenere le stringhe
char reference_line[MAX_LENGTH]="";
char compared_line[MAX_LENGTH]="";
//caratteri di confronto se i risultati sono diversi
char ref_char = ' ';
char comp_char = ' ';
//puntatori per i file
FILE *reference_file;
FILE *compared_file;
//path delle varie cartelle
char cuda[]="../CUDA";
char mpi[]="../MPI";
char st[]="../Single_Thread";
char barnes_hut[]="/Barnes-Hut/solution.txt";
char exaustive[]="/Exaustive/solution.txt";
char exaustiveArrays[]="/Exaustive/solutionArray.txt";

//restituisce 1 se i risultati combaciano, 0 altrimenti
int main(){
  //caricamento dei file
  reference_file = fopen(strcat(cuda,exaustive),"r"); //selezione path per il file di riferimento (inseriteli a mano e non rompete i coglioni)
  compared_file = fopen(strcat(cuda,exaustiveArrays),"r");//selezione path per il file da verificare
  //verifa esistenza file
  if (reference_file==NULL || compared_file==NULL){ //se uno dei file non esiste interrompi esecuzione
    printf("-> One of Files does Not Exist ! <-");
    return 0;
  }
  //lettura riga per riga
  int i,j,line,section=0;//uso 2 contatori separati per le righe poiche sono di dimensioni diverse, piu un contatore di righe e sezioni
  while(fgets(reference_line, MAX_LENGTH, reference_file)){ //finche esistono le righe nel primo file
    fgets(compared_line, MAX_LENGTH, compared_file); //leggo le righe nel secondo file
    if(flag && line<=MAX_ITERATIONS){//finche i risultati sono uguali continuo e le iterazioni minori del massimo
      i=j=0;//azzero i contatori per ogni riga
      //lettura carattere per carattere
      line++;
      while(reference_line[i]!='\0'&&compared_line[j]!='\0'){//scorro la riga carattere per carattere finche non incontro il carattere di fine riga
        if(reference_line[i]!=compared_line[j]){//se i caratteri sono diversi
          ref_char=reference_line[i];//memorizzo i caratteri per poi stamparli
          comp_char=compared_line[j];
          flag = false;//allora i risultati sono diversi
          break;
        }
        //controllo degli esponenti
        if(reference_line[i] == 'e'){//se siamo arrivati all'esponente nella prima riga
          i += 1; // aggiorno i contatori per leggere il segno dell'esponente
          j += 1;
          section++;
          //esponente riga di riferimento
          char ref_expo[6]=""; //creo ed inizializzo l'array per contere gli esponenti
          int e = 0; // contatore per aggiornare gli array
          while (reference_line[i]!=','&&reference_line[i]!='\0' ){ //finche non arriviamo alla virgola o al carattere di fine riga
            ref_expo[e] = reference_line[i];//salvo l'esponente
            e += 1;//aggiorno i contatori
            i += 1;
          }
          //esponente riga da confrontare
          char comp_expo[6]="";
          e = 0; 
          while (compared_line[j]!=','&&compared_line[j]!='\0' ){ 
            comp_expo[e] = compared_line[j];
            e += 1;
            j += 1;
          }
          //confronto numerico tra gli esponenti
          if (atoi(ref_expo)!=atoi(comp_expo)){//se gli esponenti sono diversi
            printf("\n!!! Exponents Are Different !!!\n\n");
            printf("riga %d, section %d, ref_i %d\nexp_ref = %d\nexp_com = %d\n%s%s\n",line,section-((line-1)*7),i,atoi(ref_expo),atoi(comp_expo),reference_line,compared_line);
            return 0;//interrompi
          }
        }else{   
          //aggiorno i contatori per ogni carattere
          i++;
          j++; 
        }
      }
    }else{//se i risultati sono diversi, interrompo
      break;
    }
  }
  //chiusura dei file
  fclose(reference_file);
  fclose(compared_file);
  //verifica finale tramite flag
  if (flag){//se true i risultati combaciano
    printf("\nCheck Complete!\nResults are Identical\n");
    return 1;//ritorna 1 in caso positivo
  }//altrimenti no
  printf("\n!!! Results Are Different !!!\n\n");
  //indico riga, sezione, caratteri e le righe stesse dove si verifica la differenza
  printf("line %d, section %d\nrefe_line[%d]='%c'\ncomp_line[%d]='%c'\n%s%s\n",line,1+section-((line-1)*7),i,ref_char,j,comp_char,reference_line,compared_line);
  return 0;
}