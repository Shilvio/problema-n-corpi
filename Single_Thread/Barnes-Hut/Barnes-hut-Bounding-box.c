#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

int numberBody, seed, maxTime = 1;
char fileInput[] = "../../Generate/particle.txt";
double const G = 6.67384E-11; // costante gravitazione universale
double const THETA = 0.75;     // thetha per il calcolo delle forze su particell
double maxSize;
double up,down,left,right;

#define min(i, j) (((i) < (j)) ? (i) : (j))
#define max(i, j) (((i) > (j)) ? (i) : (j))

// int count = 0;

// struct particella
typedef struct particle
{
    double x;      // posizione x
    double y;      // posizione y
    double mass;   // massa
    double forceX; // forza applicata alla particella sull' asse x
    double forceY; // forza applicata alla particella sull' asse y
    double velX;   // velocità sull' asse x
    double velY;   // velocità sull' asse y
} particle;

typedef struct massCenter
{
    double x;    // poszione x del centro di massa
    double y;    // posizione y del centro di massa
    double mass; // massa totale al centro di massa
} massCenter;

typedef struct quadTree
{
    char id[20];         // index del nodo utile solo al debug
    double up;           // dimensione superiore del quadrante
    double down;         // dimensione inferiore del quadrante
    double left;         // dimensione sinistra del quadrante
    double right;        // dimensione destra del quadrante
    particle *p;         // puntatore a una particella
    massCenter *mc;      // centro di massa per quadrante
    struct quadTree *nw; // ramo nord ovest dell' albero guardando la suddivisione del quandrante
    bool div;            // check della divisione dell' albero
    struct quadTree *ne; // ramo nord est dell' albero guardando la suddivisione del quandrante
    struct quadTree *sw; // ramo sud ovest dell' albero guardando la suddivisione del quandrante
    struct quadTree *se; // ramo sud est dell' albero guardando la suddivisione del quandrante
                         //    _________________________
                         //   |          4 |          1 |
                         //   |    (NW)    |    (NE)    |
                         //   |            |            |
                         //   |     -+     |     ++     |
                         //   |____________|____________|
                         //   |          3 |          2 |
                         //   |     --     |     +-     |
                         //   |            |            |
                         //   |    (SW)    |    (SE)    |
                         //   |____________|____________|

} quadTree;

void boundingBox(particle *p1)
{
    up = p1[0].y;
    down = p1[0].y;
    left = p1[0].x;
    right = p1[0].x;
    for (int i = 0; i < numberBody; i++)
    {
        right = max(p1[i].x, right);
        left = min(p1[i].x, left);
        up = max(p1[i].y, up);
        down = min(p1[i].y, down);
    }

    right += 1;
    left -= 1;
    up += 1;
    down -= 1;
    printf("boundingbox: up: %e, down: %e,left: %e,right: %e \n\n",up,down,left,right);
}

// funzione che analizza la presenza della particella in un dato quadrante
bool contains(quadTree *t, particle *p)
{

    return (p->x <=  t->right &&
            p->x > t->left &&
            p->y <=  t->up &&
            p->y > t->down );
}

// funzione che inizializza l' albero o nuovi nodi di esso
struct quadTree *newNode(double up,double down,double left,double right, char *idF, char *son)
{
    quadTree *t = (quadTree *)malloc(sizeof(quadTree));
    if (t == NULL)
    {
        // controllo di memori
        printf("out of memory");
        exit(1);
    }
    t->id[0] = '\0';
    strcat(t->id, idF);
    strcat(t->id, ":");
    strcat(t->id, son);

    t->up = up;
    t->down = down;
    t->left = left;
    t->right = right;
    t->div = false;
    t->p = NULL;

    t->ne = NULL;
    t->se = NULL;
    t->nw = NULL;
    t->sw = NULL;
    t->mc = NULL;

    return t;
}

// funzione che divide l' albero nelle sue 4 ramificazioni
void divide(quadTree *t)
{
    if (t->div)
    {
        return;
    }
    // printf("quand 1 %f %f %f\n",t->s/2, t->x + t->s/4, t->y + t->s/4,t->id,"1");
    // printf("quand 2 %f %f %f\n",t->s/2, t->x + t->s/4, t->y - t->s/4,t->id,"2");
    // printf("quand 3 %f %f %f\n",t->s/2, t->x - t->s/4, t->y - t->s/4,t->id,"3");
    // printf("quand 4 %f %f %f\n",t->s/2, t->x - t->s/4, t->y + t->s/4,t->id,"4");
    double horizontalSize=(t->right-t->left)/2;
    double verticalSize=(t->up-t->down)/2;

    double horizontalCenter=horizontalSize+t->left;
    double verticalCenter=verticalSize+t->down;
    
    horizontalSize/=2;
    verticalSize/=2;
    
    t->ne = newNode(t->up, verticalCenter, horizontalCenter, t->right, t->id, "1");
    t->se = newNode(verticalCenter,t->down,horizontalCenter, t->right, t->id, "2");
    t->sw = newNode(verticalCenter,t->down,t->left,horizontalCenter, t->id, "3");
    t->nw = newNode(t->up, verticalCenter,t->left,horizontalCenter, t->id, "4");
    t->div = true; // assegno true al flag di divisione
    return;
}

// funzione che inserisce le particelle in un determinato quadrante (max 1 per quadrante)
void insert(particle *p, quadTree *t)
{
    // printf("%d\n",count-1);
    if (!contains(t, p))
    {
        // printf("uscito\n");
        return;
    }
    if (t->p == NULL)
    {
        t->p = p;
        // printf("no particle in quadrant");
        return;
    }

    if (!t->div)
    {

        divide(t);
        insert(t->p, t->nw);
        insert(t->p, t->ne);
        insert(t->p, t->sw);
        insert(t->p, t->se);
    }

    insert(p, t->nw);
    insert(p, t->ne);
    insert(p, t->sw);
    insert(p, t->se);
    return;
}

// funzione che costruisce l' albero
void buildquadTree(particle *p1, quadTree *t)
{
    for (int i = 0; i < numberBody; i++)
    {
        if (contains(t, &p1[i]))
        {
            printf("out rangeparticle\n");
            exit(1);
        }
        insert(&p1[i], t);
    }
}

// funzione che stampa l' albero per eventuale debug e visualizzazione
void printer(quadTree *t, int i)
{
    if (t->p == NULL)
    {
        for (int j = 0; j < i; j++)
        {
            printf("     ");
        }
        printf("vuoto\n");
        return;
    }

    for (int j = 0; j < i; j++)
    {
        printf("     ");
    }

    if (t->mc == NULL)
    {
        printf("id=%s\n", t->id);
    }
    else
    {
        printf("id=%s cm(x= %e, y= %e, mass= %e)\n", t->id, t->mc->x, t->mc->y, t->mc->mass);
    }

    i += 1;
    if (t->div)
    {
        for (int j = 0; j < i; j++)
        {
            printf("     ");
        }
        printf("ne\n");
        printer(t->ne, i + 1);
        for (int j = 0; j < i; j++)
        {
            printf("     ");
        }
        printf("se\n");
        printer(t->se, i + 1);
        for (int j = 0; j < i; j++)
        {
            printf("     ");
        }
        printf("sw\n");
        printer(t->sw, i + 1);
        for (int j = 0; j < i; j++)
        {
            printf("     ");
        }
        printf("nw\n");
        printer(t->nw, i + 1);
    }
    else
    {
        for (int j = 0; j < i; j++)
        {
            printf("     ");
        }
        printf(" pX %e pY %e mass %e\n", t->p->x, t->p->y, t->p->mass);
    }
    return;
}

// funzione che stampa i dati dell' array di particelle
void printerAlt(particle *p1)
{
    for (int i = 0; i < numberBody; i++)
    {
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
}

// funzione che legge il file di inpiut
void getInput(FILE *file, particle *p1)
{ // popolo l'array con le particelle nel file
    for (int i = 0; i < numberBody; i++)
    {                                                                                               // prendo i dati per tutti i corpi
        fscanf(file, "%lf%lf%lf%lf%lf", &p1[i].x, &p1[i].y, &p1[i].mass, &p1[i].velX, &p1[i].velY); // prendo i dati dal file
        p1[i].forceX = 0;                                                                           // imposto le forze iniziali a zero
        p1[i].forceY = 0;
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n", p1[i].x, p1[i].y, p1[i].mass, p1[i].forceX, p1[i].forceY, p1[i].velX, p1[i].velY);
    }
    fclose(file); // chiudo il file
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
// calcolo il centro di massa e la massa totale per le particelle in ogni sottoquadrato
massCenter *centerMass(quadTree *c)
{
    massCenter *mc = malloc(sizeof(massCenter));
    if (!c->div)
    { // se non è stato diviso
        if (c->p == NULL)
        { // se la particella non c'e
            mc->x = 0;
            mc->y = 0;
            mc->mass = 0;
        }
        else
        {                    // se c'è una particella
            mc->x = c->p->x; // allora la massa e il centro di massa sono la posizione e la massa della particella stessa
            mc->y = c->p->y;
            mc->mass = c->p->mass;
        }

        c->mc = mc;
        // printf("id=%s\n",c->id);
        // printf("id=%s cm(x= %e, y= %e, mass= %e)\n",c->id,c->mc->x,c->mc->y,c->mc->mass);
        // printf(" pX %e pY %e mass %e\n",c->p->x,c->p->y,c->p->mass);
        return mc;
    }
    // altrimenti per tutti i 4 figli
    massCenter *ne = centerMass(c->ne); // calcolo centro di massa alla radice del quadrante nord est (1)
    massCenter *se = centerMass(c->se); // calcolo centro di massa alla radice del quadrante sud est (2)
    massCenter *sw = centerMass(c->sw); // calcolo centro di massa alla radice del quadrante sud ovest (3)
    massCenter *nw = centerMass(c->nw); // calcolo centro di massa alla radice del quadrante  ord ovest (4)
    // la massa di un nodo è la somma delle masse dei
    // figli = mass(1) + mass(2) + mass(3) + mass(4)
    mc->mass = ne->mass + se->mass + sw->mass + nw->mass;
    // il centro di massa di un nodo è la somma pesata dei centri di massa dei
    // poizione = (mass(1)*cm(1) + mass(2)*cm(2) + mass(3)*cm(3) + mass(4)*cm(4)) / mass
    mc->x = (ne->mass * ne->x + nw->mass * nw->x + se->mass * se->x + sw->mass * sw->x) / mc->mass;
    mc->y = (ne->mass * ne->y + nw->mass * nw->y + se->mass * se->y + sw->mass * sw->y) / mc->mass;

    c->mc = mc;
    // printf("id=%s\n",c->id);
    // printf("mod id=%s cm(x= %e, y= %e, mass= %e)\n",c->id,c->mc->x,c->mc->y,c->mc->mass);

    return mc;
}

// funzione per applicare le forze alle particelle
void threeForce(quadTree *t, particle *p)
{
    // printf("id=%s",t->id);
    // printf(" cmX=%f cmY=%f\n",t->mc->x,t->mc->y);
    double dist = sqrt(pow(p->x - t->mc->x, 2) + pow(p->y - t->mc->y, 2));
    if (dist == 0)
    {
        return;
    }
    if (!t->div) // se non e diviso
    {
        if (t->p != NULL) // se c'è una sola particella
        {
            double xDiff = p->x - t->mc->x;                                // calcolo la distanza tra la particella 1 e la 2
            double yDiff = p->y - t->mc->y;                                // (il centro di massa del nodo = particella)
            double cubeDist = dist * dist * dist;                          // elevo al cubo la distanza e applico la formula di newton
            p->forceX -= ((G * p->mass * t->mc->mass) / cubeDist) * xDiff; // per il calcolo della forza sui 2 assi
            p->forceY -= ((G * p->mass * t->mc->mass) / cubeDist) * yDiff;
            // printf("%e\n",((G * p->mass * t->mc->mass) / cubeDist) * yDiff);
            // printf("x=%e y=%e px=%e py=%e\n",t->mc->x,t->mc->y,p->x,p->velY);
        }
        return;
    }

    if ((t->right-t->left) / dist < THETA) // uso il theta come filtro per decidere la scala di approssimazione
    {
        double xDiff = p->x - t->mc->x; // calcolo le forze come espresso sopra
        double yDiff = p->y - t->mc->y;
        double cubeDist = dist * dist * dist;
        p->forceX -= ((G * p->mass * t->mc->mass) / cubeDist) * xDiff;
        p->forceY -= ((G * p->mass * t->mc->mass) / cubeDist) * yDiff;
        return;
    }
    threeForce(t->ne, p); // applico la stessa funzione attraverso figli del nodo
    threeForce(t->se, p);
    threeForce(t->sw, p);
    threeForce(t->nw, p);
    return;
}
// funzione che calcola la nuova posizione e la nuova velocità della particella
void calculatePosition(particle *p, int time)
{
    p->x += time * p->velX;
    p->y += time * p->velY;
    p->velX += time / p->mass * p->forceX;
    p->velY += time / p->mass * p->forceY;
}
// funzione per distruggere l' albero a ogni iterazione dell' algoritmo
// scandaglio l' albero e libero mano mano la memoria
void destroyTree(quadTree *c)
{
    if (!c->div)
    {
        free(c);
        return;
    }
    destroyTree(c->ne);
    destroyTree(c->se);
    destroyTree(c->sw);
    destroyTree(c->nw);
    free(c);
    return;
}

// funzione per l' inizializzazione dell' albero
void initialQuad(quadTree *c)
{
    c->up=up;
    c->down=down;
    c->left=left;
    c->right=right;
    c->id[0] = '1';
    c->id[1] = '\0';
    c->div = false;
    c->p = NULL;
    c->ne = NULL;
    c->se = NULL;
    c->nw = NULL;
    c->sw = NULL;
    c->mc = NULL;
}

// funzione di calcolo e coordinamento del programma
void compute(particle *p1, int time)
{

    for (int j = 0; j < time; j++)
    {
        boundingBox(p1);
        quadTree *c = (quadTree *)malloc(sizeof(quadTree));

        initialQuad(c);
        for (int i = 0; i < numberBody; i++)
        {
            // printf("%d number=%d x=%e y=%e \n",count,i,p1[i].x,p1[i].y);
            // count++;
            if (!contains(c, &p1[i]))
            { // se la particella esce dal quadrato principale
                printf("%d out range particle\n", i);
                exit(1);
            }
            insert(&p1[i], c);
        }
        printer(c,0);
        centerMass(c);
        // printf("\ncalcolato il centro di massa\n\n");
        // printer(c,0);
        for (int i = 0; i < numberBody; i++)
        {
            // printf("\n");
            threeForce(c, &p1[i]);

            calculatePosition(&p1[i], 1);
        }
        // printf("\ncalcolato lo spostamento\n\n");
        // printerAlt(p1);
        destroyTree(c);
        // printer(c,0);
    }
}

int main()
{
    FILE *file = initial();                               // apro il file
    particle *p1 = malloc(sizeof(particle) * numberBody); // alloco la memoria per l'array che conterrà tutte le particelle (p1)
    getInput(file, p1);                                   // lancio lo script per recuperare le particelle
    printf("\n");

    compute(p1, maxTime); // applico la funzione compute per gestire il programma
    //printerAlt(p1);       // stampo i risultati
    return 0;
}