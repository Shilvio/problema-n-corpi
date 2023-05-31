#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int numberBody, seed, maxTime = 5;
char fileInput[] = "../../Generate/particle.txt";
double const G = 6.67384E-11;
// double const G = 1;
double const deltaTime = 1;

int nProcess;
int id;

int bodyForProcess;
int remainBodies;

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

void getInput(FILE *file, particle *p1)
{

  // ci spostiamo nel punto giusto del file--------------------------------------------------------------------------------------

  for (int i = 0; i < id * bodyForProcess; i++)
  {
    // getline(&buffer,&size,stdin);
    fscanf(file, "%lf%lf%lf%lf%lf", &p1[0].x, &p1[0].y, &p1[0].mass, &p1[0].velX, &p1[0].velY);
    // pritnf()
    //  printf("id:%d ciao\n", id);
  }

  int pre;
  int bodyToRead = bodyForProcess;

  if (id < remainBodies)
  {
    pre = id;
    bodyToRead++;
  }
  else
  {
    pre = remainBodies;
  }

  for (int i = 0; i < pre; i++)
  {
    fscanf(file, "%lf%lf%lf%lf%lf", &p1[0].x, &p1[0].y, &p1[0].mass, &p1[0].velX, &p1[0].velY);
  }

  // prendo i dati per tutti i corpi
  for (int i = 0; i < bodyToRead; i++)
  {
    // prendo i dati dal file
    fscanf(file, "%lf%lf%lf%lf%lf", &p1[i].x, &p1[i].y, &p1[i].mass, &p1[i].velX, &p1[i].velY);
    // imposto le forze iniziali a zero
    p1[i].forceX = 0;
    p1[i].forceY = 0;
    printf("id:%d particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", id, p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
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

int main()
{

  mpiStart();
  FILE *file = initial();

  bodyForProcess = numberBody / nProcess;
  remainBodies = numberBody % nProcess;

  particle *p1 = malloc(sizeof(particle) * (bodyForProcess + 1));
  getInput(file, p1);

  MPI_Finalize();
}