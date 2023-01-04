#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

int numberBody,seed,maxTime=2;
char fileInput[]="../../Generate/particle.txt";
double const G=6.67384E-11; // costante gravitazione universale

typedef struct particle{
    double x;      // oszione x della particella
    double y;      // posizione y della particella
    double mass;   // massa della particella
    double forceX; // forza applicata alla particella sull' asse x
    double forceY; // forza applicata alla particella sull' asse y
    double velX;   // velocità della particella a un dato momento sull' asse x
    double velY;   // velocità della particella a un dato momento sull' asse y
}particle;

typedef struct quadTree{
    char id[20];
    double x;            // x del centro dell' albero
    double y;            // y del centro del' albero
    double s;            // dimensione
    particle* p;         // puntatore a una particella
    bool div;            // check della divisione dell' albero
    struct quadTree* nw; // ramo nord ovest dell' albero guardando la suddivisione del quandrante
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
    if(!contains(t,p)){  
        return;
    }
                                                                            //printf("ciao");   
    if(t->p==NULL){
        t->p=p;
                                                                            //printf("no particle in quadrant");
        return;
    }
    
    if(!t->div){     
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
    printf("id=%s\n",t->id);
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
        printf(" pX %f pY %f\n",t->p->x,t->p->y);
    }
    return;
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

void initialQuad(quadTree* c){
    c->x=0;
    c->y=0;
    c->s=1000;
    c->id[0]='1';
    c->id[1]='\0';
    c->div=false;
    c->p=NULL;
    c->ne=NULL;
    c->se=NULL;
    c->nw=NULL;
    c->sw=NULL;
}

int main(){
    FILE* file=initial();       //alloco la memoria per l'array che conterrà tutte le particelle (p1)
    particle* p1=malloc(sizeof(particle)*numberBody);   //popolo l'array
    getInput(file,p1);
    
    quadTree* c=(quadTree*)malloc(sizeof(quadTree));
    initialQuad(c);
    
    for(int i=0;i<numberBody;i++){
        particle* p=&p1[i];
        insert(p,c);
    }
    printer(c,0);
    return 0;
}