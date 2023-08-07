#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// numero di corpi, seed di generazione, numero di iterazioni globali
int numberBody, seed, maxTime = 5;
char fileInput[] = "../../Generate/particle.txt";
// costante di gravitazione universale
double const G = 6.67384E-11;
// numero di unità di tempo per ogni iterazione
double const deltaTime = 1;

// numero processi
int nProcess;
// id processo
int id;

// numero di corpi minimo per processo
int bodyForProcess;
// corpi rimanenti dopo la divisione intera dei corpi/processi
int remainBodies;
// punto in cui partono i body del processo
int myPoint;
// numero di corpi da leggere
int bodyToRead;

// struct particella in MPI
MPI_Datatype MPIParticle;

// struct particelle in c
typedef struct particle
{
  double x;      // valore di x della particella
  double y;      // valore di y della particella
  double mass;   // massa della particella
  double forceX; // forza sull'asse delle x
  double forceY; // forza sull'asse delle y
  double velX;   // velocità sull' asse delle x
  double velY;   // velocità sull' asse delle y
} particle;

// funzione calcolo delle forze
void calculateTotalForce(particle *p1, int j)
{
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
  }
}

// calcola posizioni e velocità delleparticelle
void calculatePosition(particle *p, int time)
{
  p->x += time * p->velX;
  p->y += time * p->velY;
  p->velX += time / p->mass * p->forceX;
  p->velY += time / p->mass * p->forceY;
}

// funzione che stampa i risultati su file
void printerFile(particle *p1)
{
  FILE *solution = fopen("solution.txt", "w");
  for (int i = 0; i < numberBody; i++)
  {
    fprintf(solution, "%e,%e,%e,%e,%e,%e,%e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
  }
  fclose(solution);
}

// funzione MPI per il brodcast dei risultati post iterazione
void multiBrodcast(particle *p1)
{
  // definisce le posizioni dei corpi rispettiva al numero di processi
  int count = 0;
  for (int i = 0; i < nProcess; i++)
  {
    if (i < remainBodies)
    {
      MPI_Bcast(&p1[count], bodyForProcess + 1, MPIParticle, i, MPI_COMM_WORLD);
      count += bodyForProcess + 1;
    }
    else
    {
      MPI_Bcast(&p1[count], bodyForProcess, MPIParticle, i, MPI_COMM_WORLD);
      count += bodyForProcess;
    }
  }
}

// funzione di carimaneto dati dal file
void getInput(FILE *file, particle *p1)
{
  // per tutti i corpi
  for (int i = 0; i < numberBody; i++)
  {
    // caricamento dei dati dal file
    fscanf(file, "%lf%lf%lf%lf%lf", &p1[i].x, &p1[i].y, &p1[i].mass, &p1[i].velX, &p1[i].velY);
    // imposto le forze iniziali a zero
    p1[i].forceX = 0;
    p1[i].forceY = 0;
  }
  // chiusura il file
  fclose(file);
}

//
FILE *initial()
{
  // mi apro il file in lettura
  FILE *file = fopen(fileInput, "r");
  // prendo il seed
  fscanf(file, "%d", &seed);
  // prendo il numero di corpi
  fscanf(file, "%d", &numberBody);

  // setup for MPI
  bodyForProcess = numberBody / nProcess;
  remainBodies = numberBody % nProcess;
  bodyToRead = bodyForProcess;

  // trova il punto dell' array da passare per essere usato dal thread e ripartiziona le particelle in più
  myPoint = id * bodyForProcess;
  if (id < remainBodies)
  {
    myPoint += id;
    bodyToRead++;
  }
  else
  {
    myPoint += remainBodies;
  }

  return file;
}

// setup per MPI, creazione del comunicatore
void mpiStart()
{
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  MPI_Type_contiguous(sizeof(particle), MPI_BYTE, &MPIParticle);
  MPI_Type_commit(&MPIParticle);
}

// funzione che stampa i dati dell' array di particelle
void printer(particle *p1)
{
  for (int i = 0; i < numberBody; i++)
  {
    printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
  }
}

int main()
{

  // inizializza il programma MPI
  mpiStart();

  // prende i dati dal file delle particelle
  FILE *file = initial();

  // crea array particelle
  particle *p1 = malloc(sizeof(particle) * numberBody);

  if (id == 0)
  {
    getInput(file, p1);
  }

  // broadcast dell' array di particelle a tutti i thread usati
  MPI_Bcast(p1, numberBody, MPIParticle, 0, MPI_COMM_WORLD);

  // per tutti i time step
  for (int i = 0; i < maxTime; i++)
  {
    // calcolo delle forze
    for (int i = 0; i < bodyToRead; i++)
      calculateTotalForce(p1, myPoint + i);
    // calcolo delle posizioni
    for (int i = 0; i < bodyToRead; i++)
      calculatePosition(&p1[myPoint + i], deltaTime);

    multiBrodcast(p1);
  }
  // stampa dei risultati da parte del processo 0
  if (id == 0)
  {
    // printer(p1); // print per debug
    printerFile(p1);
  }
  fclose(file);
  free(p1);
  // ferma il programma mpi
  MPI_Finalize();
}
