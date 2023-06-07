#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int numberBody, seed, maxTime = 1;
char fileInput[] = "../../Generate/particle.txt";
double const G = 6.67384E-11;
// double const G = 1;
double const deltaTime = 1;

int nProcess;
int id;

int bodyForProcess;
int remainBodies;
int myPoint;
int bodyToRead;

MPI_Datatype MPIParticle;
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

void calculateTotalForce(particle *p1)
{
  for (int i = myPoint; i < myPoint + bodyToRead; i++)
  {
    // printf("id = %d body: %d end %d\n", id, i, myPoint + bodyToRead);
    p1[i].mass = id;
  }
}

void multiBrodcast(particle *p1)
{
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
    // printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
  }
  // chiudo il file
  fclose(file);
}

FILE *initial()
{
  // mi apro il file in lettura
  FILE *file = fopen(fileInput, "r");
  // prendo il seed
  fscanf(file, "%d", &seed);
  if (id == 0)
    printf("%d\n", seed);
  // prendo il numero di corpi
  fscanf(file, "%d", &numberBody);
  if (id == 0)
    printf("%d\n", numberBody);

  // setup for MPI
  bodyForProcess = numberBody / nProcess;
  remainBodies = numberBody % nProcess;
  bodyToRead = bodyForProcess;

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

void mpiStart()
{
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  MPI_Type_contiguous(sizeof(particle), MPI_BYTE, &MPIParticle);
  MPI_Type_commit(&MPIParticle);
}

void printer(particle *p1)
{
  for (int i = 0; i < numberBody; i++)
  {
    printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
  }
}

int main()
{

  mpiStart();

  FILE *file = initial();

  particle *p1 = malloc(sizeof(particle) * numberBody);
  if (id == 0)
  {
    getInput(file, p1);
  }
  MPI_Bcast(p1, numberBody, MPIParticle, 0, MPI_COMM_WORLD);

  for (int i = 0; i < maxTime; i++)
  {
    calculateTotalForce(p1);
    multiBrodcast(p1);
  }
  if (id == 0)
    printer(p1);

  MPI_Finalize();
}
