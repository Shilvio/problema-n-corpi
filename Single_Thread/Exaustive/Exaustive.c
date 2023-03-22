#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//variabili per il calcolo del tempo di esecuzione


int numberBody, seed, maxTime = 5;
char fileInput[] = "../../Generate/particle.txt";
double const G = 6.67384E-11;
//double const G = 1;
double const deltaTime=1;

typedef struct particle
{
    double x;
    double y;
    double mass;
    double forceX;
    double forceY;
    double velX;
    double velY;
} particle;

// calcolo del cambio della posizione per una particella di interesse per un dato intervallo di tempo e velocità
void calculatePosition(particle *p, int time)
{
    p->x += time * p->velX;
    p->y += time * p->velY;
    p->velX += time / p->mass * p->forceX;
    p->velY += time / p->mass * p->forceY;
}

// calcolo tutte le forze relative ad una particella di interesse
void calculateTotalForce(particle *p1, int j)
{
    // itero per tutte le particelle tranne quella di interesse
    for (int i = 0; i < numberBody; i++)
    {
        if (i == j)
        {
            continue;
        }
        // calcolo delle forze
        double xDiff = p1[j].x - p1[i].x;
        double yDiff = p1[j].y - p1[i].y;
        double dist = sqrt(xDiff * xDiff + yDiff * yDiff);
        double cubeDist = dist * dist * dist;
        p1[j].forceX -= ((G * p1[j].mass * p1[i].mass) / cubeDist) * xDiff;
        p1[j].forceY -= ((G * p1[j].mass * p1[i].mass) / cubeDist) * yDiff;
                                                                                        //printf("\n%e",((G * p1[j].mass * p1[i].mass) / cubeDist) * yDiff);
                                                                                        //printf("\n px=%e py=%e fX=%e fY=%e \n",p1[j].x,p1[j].y,p1[j].forceX,p1[j].forceY);
                                                                                        //printf("px=%e py=%e \n",p1[i].x,p1[i].y);
    }
}

void printerFile(particle *p1)
{
    FILE* solution=fopen("solution.txt","w");
    for (int i = 0; i < numberBody; i++)
    {
        fprintf(solution,"%e,%e,%e,%e,%e,%e,%e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
    fclose(solution);
}

void printer(particle *p1)
{
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
}

// calcolo il movimento delle particelle nel tempo richiesto
void compute(int time, particle *p1)
{
    // itero per il tempo tichiesto
    for (int i = 0; i < time; i++)
    {
        // itero per ogni corpo
        for (int j = 0; j < numberBody; j++)
        {
            // calcolo tutte le forze relative ad una particella di interesse
            calculateTotalForce(p1, j);
        }
        // itero per ogni corpo
        for (int j = 0; j < numberBody; j++)
        {
            // calcolo del cambio della posizione per una particella di interesse
            calculatePosition(&p1[j], deltaTime);
        }
    }
}

// popolo l'array con le particelle nel file
void getInput(FILE *file, particle *p1)
{
    // prendo i dati per tutti i corpi
    for (int i = 0; i < numberBody; i++)
    {
        // prendo i dati dal file
        fscanf(file, "%lf%lf%lf%lf%lf", &p1[i].x, &p1[i].y, &p1[i].mass, &p1[i].velX, &p1[i].velY);
        // imposto le forze iniziali a zero
        p1[i].forceX = 0;
        p1[i].forceY = 0;
        //printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
    // chiudo il file
    fclose(file);
}

// aprire il file e prendere i primi valori (seed e numero di corpi)
FILE *initial()
{
    // mi apro il file in lettura
    FILE *file = fopen(fileInput, "r");
    // prendo il seed
    fscanf(file, "%d", &seed);
    printf("%d\n", seed);
    // prendo il numero di corpi
    fscanf(file, "%d", &numberBody);
    printf("%d\n", numberBody);
    return file;
}

int main()
{   
    clock_t start, end;
    double cpu_time_used;   
    // apro il file dove si trovano tutte le particelle
    FILE *file = initial();
    // alloco la memoria per l'array che contierrà tutte le particelle (p1)
    particle *p1 = malloc(sizeof(particle) * numberBody);
    // popolo l'array
    getInput(file, p1);
    // calcolo il movimento delle particelle nel tempo richiesto
    start = clock();  
    
                                           //avvio il clock per calcolo del tempo di esecuzione
    compute(maxTime, p1);
    end = clock();                                             //fermo il clock per il calcolo del tempo di esecuzione
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC; //converto in secondi
    //cpu_time_used = (end - start);                    //converto in clock time
    printf("\nla funzione ha richiesto: %e secondi\n", cpu_time_used); 

    printf("\n");
    printer(p1);

    printerFile(p1);

    fclose(file);
    free(p1);
    exit(1);
}