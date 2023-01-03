#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

int numberBody,seed,maxTime=2;
char fileInput[]="../../Generate/particle.txt";
double const G=6.67384E-11;

typedef struct particle{
    double x;
    double y;
    double mass;
    double forceX;
    double forceY;
    double velX;
    double velY;
}particle;

typedef struct quadTree{
    double x;
    double y;
    double s; //dimensione
    particle p;
    bool div;
    struct quadTree* nw;
    quadTree* ne;
    quadTree* sw;
    quadTree* se;
}quadTree;

bool contains(quadTree* t, particle p){
        return (p.x < (t->x+t->s/2) &&
        p.x > (t->x-t->s/2) &&
        p.y < (t->y+t->s/2) &&
        p.y > (t->y-t->s/2));
        
}

quadTree* newNode(double s,double x, double y){
    quadTree* t = (quadTree*) malloc(sizeof(quadTree));

    if(t==NULL){
        printf("out of memory");
        exit(1);
    }
    t->s=s;
    t->x=x;
    t->y=y;
    t->div=false;
    t->p=NULL;

    t->ne=NULL;
    t->se=NULL;
    t->nw=NULL;
    t->sw=NULL,
    
    return t;
}

void divide(quadTree* t){
    if(t->div){
        return;
    }
    
    t->ne = newNode(t->s/2, t->x + t->s/4, t->y + t->s/4); 
    t->se = newNode(t->s/2, t->x + t->s/4, t->y - t->s/4); 
    t->sw = newNode(t->s/2, t->x - t->s/4, t->y - t->s/4); 
    t->nw = newNode(t->s/2, t->x - t->s/4, t->y + t->s/4); 

    t->div = true;
}

void insert(particle p,quadTree* t){
    if(!contains(t,p)){
        return;
    }

    if(t->p==NULL){
        t->p=p;
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


void buildquadTree(particle* p1,quadTree* t){
    fot(int i=0;i<numberBody;i++){
        insert(p1[i],t);
    }
}

int main(){
    printf("puzzi");
    
}