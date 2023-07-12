#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int numberBody, seed, maxTime = 5;
char fileInput[] = "../../Generate/particle.txt";
double const G = 6.67384E-11;
// double const G = 1;
double const deltaTime = 1;

int nProcess;
int id;

double size[4];

int bodyForProcess;
int remainBodies;
int myPoint;
int bodyToRead;

#define min(i, j) (((i) < (j)) ? (i) : (j))
#define max(i, j) (((i) > (j)) ? (i) : (j))

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

// calcola forze
void calculateTotalForce(particle *p1, int j)
{
  for (int i = 0; i < numberBody; i++)
  {
    /*if (i == j)
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
    // printf("\n%e",((G * p1[j].mass * p1[i].mass) / cubeDist) * yDiff);
    // printf("\n px=%e py=%e fX=%e fY=%e \n",p1[j].x,p1[j].y,p1[j].forceX,p1[j].forceY);
    // printf("px=%e py=%e \n",p1[i].x,p1[i].y);*/
  }
}

// calcolo la bounding box delle particelle, applicando tecniche di riduzione gpu
void boundingBox(particle *p1)
{

  double up = p1[0].y, down = p1[0].y, left = p1[0].x, right = p1[0].x;
  // controllo se ho più thread che particelle
  if (id < numberBody)
  { // ciclo il mo settore di particelle nell' array globale
    for (int i = 0; i < bodyToRead; i++)
    {
      up = max(up, p1[i + myPoint].y);
      down = min(down, p1[i + myPoint].y);
      left = min(left, p1[i + myPoint].x);
      right = max(right, p1[i + myPoint].x);
    }

    up += 1;
    down -= 1;
    left -= 1;
    right += 1;
  }

  // possibile barrier
  MPI_Allreduce(&up, &up, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&down, &down, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&right, &right, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&left, &left, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD);

  // printf("\n\nBounding box: up: %e,down: %e,left: %e,right: %e\n\n", up, down, left, right);
  size[0] = up;
  size[1] = down;
  size[2] = left;
  size[3] = right;
}

// calcola posizioni e velocità delleparticelle
void calculatePosition(particle *p, int time)
{
  p->x += time * p->velX;
  p->y += time * p->velY;
  p->velX += time / p->mass * p->forceX;
  p->velY += time / p->mass * p->forceY;
}

void printerFile(particle *p1)
{
  FILE *solution = fopen("solution.txt", "w");
  for (int i = 0; i < numberBody; i++)
  {
    fprintf(solution, "%e,%e,%e,%e,%e,%e,%e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
  }
  fclose(solution);
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

  boundingBox(p1);
  printf("\n\nBounding box: up: %e,down: %e,left: %e,right: %e\n\n", size[0], size[1], size[2], size[3]);
  // createTree();
  // calculateCenterMass();
  if (id == 0)
  {
    // printer(p1);
    // printerFile(p1);
  }
  fclose(file);
  free(p1);
  // ferma il programma mpi server
  MPI_Finalize();
}
