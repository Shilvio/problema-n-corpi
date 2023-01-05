#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

int numberBody,seed,maxTime=3;
char fileInput[]="../../Generate/particle.txt";
double const G=6.67384E-11; // costante gravitazione universale
//double const G=1;
double const THETA= 0.5; //thetha per il calcolo delle forze su particell
double maxSize =6.162025e+070;
//double maxSize = 100;
int count=0;

typedef struct particle{
    double x;      // oszione x della particella
    double y;      // posizione y della particella
    double mass;   // massa della particella
    double forceX; // forza applicata alla particella sull' asse x
    double forceY; // forza applicata alla particella sull' asse y
    double velX;   // velocità della particella a un dato momento sull' asse x
    double velY;   // velocità della particella a un dato momento sull' asse y
}particle;

typedef struct massCenter{
    double x;      // poszione x del centro di massa
    double y;      // posizione y del centro di massa
    double mass;   // massa totale al centro di massa
}massCenter;

typedef struct quadTree{
    char id[20];
    double x;            // x del centro dell' albero
    double y;            // y del centro del' albero
    double s;            // dimensione
    particle* p;         // puntatore a una particella
    massCenter* mc;      // centro di massa per quadrante
    struct quadTree* nw; // ramo nord ovest dell' albero guardando la suddivisione del quandrante
    bool div;            // check della divisione dell' albero
    struct quadTree* ne; // ramo nord est dell' albero guardando la suddivisione del quandrante
    struct quadTree* sw; // ramo sud ovest dell' albero guardando la suddivisione del quandrante
    struct quadTree* se; // ramo sud est dell' albero guardando la suddivisione del quandrante
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
    
}quadTree;


// funzione che analizza la presenza della particella in un dato quadrante
bool contains(quadTree* t, particle* p){
                                                                        //printf("x %f, y %f",t->x,t->y);
        return (p->x <= (t->x+t->s/2) &&
        p->x > (t->x-t->s/2) &&
        p->y <= (t->y+t->s/2) &&
        p->y > (t->y-t->s/2));       
}

// funzione che inizializza l' albero 
struct quadTree* newNode(double s,double x, double y,char* idF,char* son){  
    quadTree* t = (quadTree*) malloc(sizeof(quadTree)); 
    if(t==NULL){
        printf("out of memory");
        exit(1);
    }
    t->id[0]='\0';
    strcat(t->id,idF);
    strcat(t->id,":");
    strcat(t->id,son);
    
    t->s=s;
    t->x=x;
    t->y=y;
    t->div=false;
    t->p=NULL;

    t->ne=NULL;
    t->se=NULL;
    t->nw=NULL;
    t->sw=NULL;
    t->mc=NULL;
    

    return t;
}

// funzione che divide l' albero nelle sue 4 ramificazioni
void divide(quadTree* t){
    if(t->div){
        return;
    }
                                                                        //printf("quand 1 %f %f %f\n",t->s/2, t->x + t->s/4, t->y + t->s/4,t->id,"1");
                                                                        //printf("quand 2 %f %f %f\n",t->s/2, t->x + t->s/4, t->y - t->s/4,t->id,"2");
                                                                        //printf("quand 3 %f %f %f\n",t->s/2, t->x - t->s/4, t->y - t->s/4,t->id,"3");
                                                                        //printf("quand 4 %f %f %f\n",t->s/2, t->x - t->s/4, t->y + t->s/4,t->id,"4");
    t->ne = newNode(t->s/2, t->x + t->s/4, t->y + t->s/4,t->id,"1");    
    t->se = newNode(t->s/2, t->x + t->s/4, t->y - t->s/4,t->id,"2");   
    t->sw = newNode(t->s/2, t->x - t->s/4, t->y - t->s/4,t->id,"3");  
    t->nw = newNode(t->s/2, t->x - t->s/4, t->y + t->s/4,t->id,"4");    
    t->div = true;  //assegno true al flag di divisione
    return;
}

// funzione che inserisce le particelle in un determinato quadrante (max 1 per quadrante)
void insert(particle* p,quadTree* t){
                                                                            //printf("%d\n",count-1);
    if(!contains(t,p)){  
                                                                            //printf("uscito\n");
        return;
    }
                                                                            //printf("ciao");   
    if(t->p==NULL){
        t->p=p;
                                                                            //printf("no particle in quadrant");
        return;
    }
    
    
    if(!t->div){     
                                                                            //printf("loop\n");
        divide(t);
        insert(t->p,t->nw);
        insert(t->p,t->ne);
        insert(t->p,t->sw);
        insert(t->p,t->se);
    }   
    
    insert(p,t->nw);
    insert(p,t->ne);
    insert(p,t->sw);
    insert(p,t->se);
    return;
}

// funzione che costruisce l' albero
void buildquadTree(particle* p1,quadTree* t){ 
    for(int i=0;i<numberBody;i++){
        if(contains(t,&p1[i])){
            printf("out rangeparticle\n");
            exit(1);
        }
        insert(&p1[i],t);
    }
}

// funzione che stampa l' albero per eventuale debug e visualizzazione
void printer(quadTree* t,int i){
    if(t->p==NULL){
        for(int j=0;j<i;j++){
        printf("     ");
        }
        printf("vuoto\n");
        return;
    }
    
    for(int j=0;j<i;j++){
        printf("     ");
    }
    
    if(t->mc==NULL){
        printf("id=%s\n",t->id);
    }
    else{
        printf("id=%s cm(x= %e, y= %e, mass= %e)\n",t->id,t->mc->x,t->mc->y,t->mc->mass);
    }

    i+=1;
    if(t->div){
        for(int j=0;j<i;j++){
            printf("     ");
        }
        printf("ne\n");
        printer(t->ne,i+1);
        for(int j=0;j<i;j++){
            printf("     ");
        }
        printf("se\n");
        printer(t->se,i+1);
        for(int j=0;j<i;j++){
            printf("     ");
        }
        printf("sw\n");
        printer(t->sw,i+1);
        for(int j=0;j<i;j++){
            printf("     ");
        }
        printf("nw\n");
        printer(t->nw,i+1);

    }else{
        for(int j=0;j<i;j++){
            printf("     ");
        }
        printf(" pX %e pY %e mass %e\n",t->p->x,t->p->y,t->p->mass);
    }
    return;
}

void printerAlt(particle* p1){
    for(int i=0;i<numberBody;i++){
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].forceX,p1[i].forceY,p1[i].velX,p1[i].velY);
    }
}

void getInput(FILE* file,particle* p1){     //popolo l'array con le particelle nel file
    for(int i=0;i<numberBody;i++){          //prendo i dati per tutti i corpi
        fscanf(file,"%lf%lf%lf%lf%lf",&p1[i].x,&p1[i].y,&p1[i].mass,&p1[i].velX,&p1[i].velY);//prendo i dati dal file
        p1[i].forceX=0;                     //imposto le forze iniziali a zero
        p1[i].forceY=0;
        printf("particle xPos= %e, yPos= %e, mass= %e, forceX= %e, forceY= %e, velX= %e, velY= %e\n",p1[i].x,p1[i].y,p1[i].mass,p1[i].forceX,p1[i].forceY,p1[i].velX,p1[i].velY);
    }
    fclose(file);//chiudo il file
}

//aprire il file e prendere i primi valori (seed e numero di corpi)
FILE* initial(){
    //mi apro il file in lettura
    FILE* file=fopen(fileInput,"r");
    //prendo il seed
    fscanf(file,"%d",&seed);
                                                                        printf("%d\n",seed);
    //prendo il numero di corpi
    fscanf(file,"%d",&numberBody);
                                                                        printf("%d\n",numberBody);
    return file;
}
// calcolo il centro di massa e la massa totale per le particelle in ogni sottoquadrato
massCenter* centerMass(quadTree* c){
    massCenter* mc=malloc(sizeof(massCenter)); 
    if(!c->div){ //se non è stato diviso
        if(c->p==NULL){ //se la particella non c'e
            mc->x=0; 
            mc->y=0;
            mc->mass=0;
            
        }else{  //se c'è una particella
            mc->x=c->p->x; //allora la massa e il centro di massa sono la posizione e la massa della particella stessa
            mc->y=c->p->y;
            mc->mass=c->p->mass;
        }
                                                                            
        c->mc=mc;
                                                                            //printf("id=%s\n",c->id);
                                                                            //printf("id=%s cm(x= %e, y= %e, mass= %e)\n",c->id,c->mc->x,c->mc->y,c->mc->mass);
                                                                            //printf(" pX %e pY %e mass %e\n",c->p->x,c->p->y,c->p->mass);
        return mc;
    }   
    //altrimenti per tutti i 4 figli 
    massCenter* ne=centerMass(c->ne); //calcolo centro di massa alla radice del quadrante nord est (1)
    massCenter* se=centerMass(c->se); //calcolo centro di massa alla radice del quadrante sud est (2)
    massCenter* sw=centerMass(c->sw); //calcolo centro di massa alla radice del quadrante sud ovest (3)
    massCenter* nw=centerMass(c->nw); //calcolo centro di massa alla radice del quadrante  ord ovest (4)
    //la massa di un nodo e la somma delle masse dei figli = mass(1) + mass(2) + mass(3) + mass(4)
    mc->mass= ne->mass + se->mass + sw->mass + nw->mass;
    // il centro di massa di un nodo e la somma pesata dei centri di massa dei figli = (mass(1)*cm(1) + mass(2)*cm(2) + mass(3)*cm(3) + mass(4)*cm(4)) / mass
    mc->x= (ne->mass*ne->x + nw->mass*nw->x + se->mass*se->x + sw->mass*sw->x) / mc->mass;
    mc->y= (ne->mass*ne->y + nw->mass*nw->y + se->mass*se->y + sw->mass*sw->y) / mc->mass;

    c->mc=mc;
                                                                            //printf("id=%s\n",c->id);
                                                                            //printf("mod id=%s cm(x= %e, y= %e, mass= %e)\n",c->id,c->mc->x,c->mc->y,c->mc->mass);
    
    return mc;
}

void threeForce(quadTree* t,particle* p){
    //printf("id=%s",t->id);
    //printf(" cmX=%f cmY=%f\n",t->mc->x,t->mc->y);
    double dist = sqrt(pow(p->x - t->mc->x,2) + pow(p->y - t->mc->y,2));
    if(dist==0){
        return;
    }
    if(!t->div){ //se non e diviso
        if(t->p!=NULL){ //se c'è una particella
            double xDiff = p->x - t->mc->x ;
            double yDiff = p->y - t->mc->y ;
            double cubeDist = dist * dist * dist;
            p->forceX -= ((G * p->mass * t->mc->mass) / cubeDist) * xDiff;
            p->forceY -= ((G * p->mass * t->mc->mass) / cubeDist) * yDiff;
                                                                                        //printf("%e\n",((G * p->mass * t->mc->mass) / cubeDist) * yDiff);
                                                                                        //printf("x=%e y=%e px=%e py=%e\n",t->mc->x,t->mc->y,p->x,p->velY);
        }
                                                                                        //printf("ciao2");
        return;
    }
    
    if (t->s/dist < THETA){
                                                                                        //printf("ciao");
        double xDiff = p->x - t->mc->x ;
        double yDiff = p->y - t->mc->y ;
        double cubeDist = dist * dist * dist;
        p->forceX -= ((G * p->mass * t->mc->mass) / cubeDist) * xDiff;
        p->forceY -= ((G * p->mass * t->mc->mass) / cubeDist) * yDiff;
        return;
    } 
    //printf("ciao3");
    threeForce(t->ne,p);
    threeForce(t->se,p);
    threeForce(t->sw,p);
    threeForce(t->nw,p);
    return;
}

void calculatePosition(particle* p,int time){
    p->x+=time*p->velX;
    p->y+=time*p->velY;
    p->velX+=time/p->mass*p->forceX;
    p->velY+=time/p->mass*p->forceY;
}
// funzione per distruggere l' albero a ogni iterazione dell' algoritmo
void destroyTree(quadTree* c){
    if(!c->div){
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

void initialQuad(quadTree* c){
    c->x=0;
    c->y=0;
    c->s=maxSize;
    c->id[0]='1';
    c->id[1]='\0';
    c->div=false;
    c->p=NULL;
    c->ne=NULL;
    c->se=NULL;
    c->nw=NULL;
    c->sw=NULL;
    c->mc=NULL;
}

void compute(particle* p1,int time){
    
    for(int j=0;j<time;j++){
        quadTree* c=(quadTree*)malloc(sizeof(quadTree));
        
        initialQuad(c);
        for(int i=0;i<numberBody;i++){
                                                                //printf("%d number=%d x=%e y=%e \n",count,i,p1[i].x,p1[i].y);
            count++;
            if(!contains(c,&p1[i])){                             //se la particella esce dal quadrato principale
                printf("%d out range particle\n",i);
                exit(1);
            }
            insert(&p1[i],c);
        }
                                                                //printer(c,0);
        centerMass(c);
                                                                //printf("\ncalcolato il centro di massa\n\n");
                                                                //printer(c,0);
                                                                //printer(c,0);
        for(int i=0;i<numberBody;i++){
                                                                //printf("\n");            
            threeForce(c,&p1[i]);

            calculatePosition(&p1[i],1);
        }
                                                                //printf("\ncalcolato lo spostamento\n\n");
                                                                //printerAlt(p1);
        destroyTree(c);
        //printer(c,0);
    }
}

int main(){
    FILE* file=initial();       //alloco la memoria per l'array che conterrà tutte le particelle (p1)
    particle* p1=malloc(sizeof(particle)*numberBody);   //popolo l'array
    getInput(file,p1);
    
    
    //printer(c,0);
                                                        printf("\n");
    
    compute(p1,maxTime);
                                                        printerAlt(p1);
    return 0;
}