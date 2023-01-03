#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

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
    char id[20];
    double x;
    double y;
    double s; //dimensione
    particle* p;
    bool div;
    struct quadTree* nw;
    struct quadTree* ne;
    struct quadTree* sw;
    struct quadTree* se;
}quadTree;

bool contains(quadTree* t, particle* p){
                                                                        //printf("x %f, y %f",p->x,p->y);
        return (p->x < (t->x+t->s/2) &&
        p->x > (t->x-t->s/2) &&
        p->y < (t->y+t->s/2) &&
        p->y > (t->y-t->s/2));
        
}

struct quadTree* newNode(double s,double x, double y,char* idF,char* son){
    
    quadTree* t = (quadTree*) malloc(sizeof(quadTree));
    
    if(t==NULL){
        printf("out of memory");
        exit(1);
    }
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
    
    t->div = true;
}

void insert(particle* p,quadTree* t){
    if(!contains(t,p)){
        return;
    }
    
    if(t->p==NULL){
        t->p=p;
                                                                            printf("no particle in quadrant");
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
    for(int i=0;i<numberBody;i++){
        insert(&p1[i],t);
    }
}

void printer(quadTree* t,int i){
    if(t->p==NULL){
        return;
    }
    if(i>5){
        return;
    }
    for(int j=0;j<i;j++){
        printf("\t");
    }
    printf("id=%s ",t->id);
    if(t->div){
        for(int j=0;j<i;j++){
            printf("\t");
        }
        printf("ne\n");
        printer(t->ne,i+1);
        for(int j=0;j<i;j++){
            printf("\t");
        }
        printf("se\n");
        printer(t->se,i+1);
        for(int j=0;j<i;j++){
            printf("\t");
        }
        printf("sw\n");
        printer(t->sw,i+1);
        for(int j=0;j<i;j++){
            printf("\t");
        }
        printf("nw\n");
        printer(t->nw,i+1);

    }
    return;
}

int main(){
    quadTree* c=(quadTree*)malloc(sizeof(quadTree));
    particle* p=(particle*)malloc(sizeof(particle));
    
    p->x=4;
    p->y=4;
    
    c->x=0;
    c->y=0;
    c->s=20;
    c->id[0]="1";
    c->div=false;
    c->p=NULL;
    c->ne=NULL;
    c->se=NULL;
    c->nw=NULL;
    c->sw=NULL;
    
    insert(p,c);
                                                                        //printf("ciao");
    //printer(c,0);
    p->x=2;
    p->y=2;
    insert(p,c);
                                                                        //printf("ciao");
    //printer(c,0);
    return 0;
}