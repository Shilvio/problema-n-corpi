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
                                                                        //printf("x %f, y %f",t->x,t->y);
        return (p->x <= (t->x+t->s/2) &&
        p->x > (t->x-t->s/2) &&
        p->y <= (t->y+t->s/2) &&
        p->y > (t->y-t->s/2));
        
}

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

    return;
}

int insert(particle* p,quadTree* t){
    if(!contains(t,p)){
        return 0;
    }

                                                                            //printf("ciao");
    
    if(t->p==NULL){
        t->p=p;
                                                                            //printf("no particle in quadrant");
        return 1;
    }else{

        if(!t->div){
            
            divide(t);
            if(!insert(t->p,t->nw)){
                if(!insert(t->p,t->ne)){
                    if(!insert(t->p,t->sw)){
                        insert(t->p,t->se);
                    }
                }
            }
        }

        if(insert(p,t->nw))
            return 1;
        if(insert(p,t->ne))
            return 1;
        if(insert(p,t->sw))
            return 1;
        if(insert(p,t->se))
            return 1;

    }
    return 0;
}


void buildquadTree(particle* p1,quadTree* t){
    for(int i=0;i<numberBody;i++){
        insert(&p1[i],t);
    }
}

void printer(quadTree* t,int i){
    if(t->p==NULL){
        for(int j=0;j<i;j++){
        printf("     ");
        }
        printf("vuoto\n");
        return;
    }
    /*if(i>5){
        return;
    }*/
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

int main(){
    quadTree* c=(quadTree*)malloc(sizeof(quadTree));
    particle* p=(particle*)malloc(sizeof(particle));
    
    p->x=4.3;
    p->y=4.3;
    
    c->x=0;
    c->y=0;
    c->s=20;
    c->id[0]='1';
    c->id[1]='\0';
    c->div=false;
    c->p=NULL;
    c->ne=NULL;
    c->se=NULL;
    c->nw=NULL;
    c->sw=NULL;
    particle* d=(particle*)malloc(sizeof(particle));
    insert(p,c);
                                                                        //printf("ciao");
    //printer(c,0);
    d->x=-0.2;
    d->y=-0.2;
    insert(d,c);
                                                                        //printf("ciao");
    printer(c,0);
    return 0;
}