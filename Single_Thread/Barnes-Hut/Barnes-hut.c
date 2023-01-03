#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

int numberBody,seed,maxTime=2;
char fileInput[]="../../Generate/particle.txt";
double const G=6.67384E-11;

typedef struct point{
    double x;
    double y;
}point;


typedef struct particle{
    point pos;
    double mass;
    double forceX;
    double forceY;
    double velX;
    double velY;
}particle;

typedef struct quadTree{
    int data; //key
    double s; //are di interesse
    point center; //centro del quadrato
    bool div; //verifica se è stato diviso

    particle* p;//particella

    quadTree* nw; //nodi figli
    quadTree* ne;
    quadTree* se;
    quadTree* sw;

}quadTree;

quadTree* newNode (int data, double s, double x, double y) 
{ 
    //alloca memoria 
    quadTree* quad = malloc(sizeof(quadTree)); 

    // se non c'è spazio chiudi
    if(quad == NULL){
        printf ("memoria non allocata\n");
        exit(1);
    }
    // assegno i dati al quadTree, s è la dimensione del quantrato, center il centro
    // p è nullo dato che non indica ancora una particella, l' albero resta non diviso, quindi divsta a false.
    quad->data = data;
    quad->s = s; 
    quad->center.x = x; quad->center.y = y;
    quad->p = NULL;
    quad->div = false;

    // inizializzo i nodi figli a null
    quad->ne = NULL; 
    quad->se = NULL;
    quad->sw = NULL;
    quad->nw = NULL;

    return (struct quadTree*)quad; 
} 

void printTree(quadTree* nd){
    if (nd == NULL)
        return;
    // mostra i dati del nodo
    printf(" %*c(%d) \n %*c[%f, %f] \n %*c[%f]\n\n",50, ' ', nd->data, 42,' ', nd->center.x, nd->center.y, 47, ' ', nd->s);
    
    // mostra i dati dei sottonodi formattati
    if(nd->ne != NULL)
        printf("%*c|NE:%d|  ",34,' ',nd->ne->data);
    if(nd->se != NULL)
        printf("|SE:%d|  ",nd->se->data);
    if(nd->sw != NULL)
        printf("|SW:%d|  ",nd->sw->data);
    if(nd->nw != NULL)
        printf("|NW:%d|",nd->nw->data);
    printf("\n\n");
    
    // va in ricorsione fino alle foglie
    display_tree(nd->ne);
    display_tree(nd->se);
    display_tree(nd->sw);
    display_tree(nd->nw);
}

//distrugge l'albero
void deconstruct_tree(quadTree* root){
    if(root != NULL){
        deconstruct_tree(root->ne);
        deconstruct_tree(root->se);
        deconstruct_tree(root->sw);  
        deconstruct_tree(root->nw);
        if(root->data<0){
            free(root->p);
        }
        free(root);
    }
}

//suddivide l'albero
void subdivide(quadTree* nd, int* track){ 
    
    if(nd == NULL){
        return;
    }

    // chiama newNode function per ogni filgio che non era Null ripetto al padre, e gli assegna un blocco di memoria di tipo(quadTree)
    // un id intero (track)  ai nodi che sono rami, ovvero quelli che non sono folgie, anche alle celle vuote viende assegnato un id definito da (track) 
    // univoco a qualsiasi albero generato
    //    _________________________
    //   |            |            |
    //   |    (NW)    |    (NE)    |
    //   |            |            |
    //   |     -+     |     ++     |
    //   |____________|____________|
    //   |            |            |
    //   |     --     |     +-     |
    //   |            |            |
    //   |    (SW)    |    (SE)    |
    //   |____________|____________|

    nd->ne = newNode(*track-1, nd->s/2, nd->center.x + nd->s/4, nd->center.y + nd->s/4); 
    nd->se = newNode(*track-2, nd->s/2, nd->center.x + nd->s/4, nd->center.y - nd->s/4); 
    nd->sw = newNode(*track-3, nd->s/2, nd->center.x - nd->s/4, nd->center.y - nd->s/4); 
    nd->nw = newNode(*track-4, nd->s/2, nd->center.x - nd->s/4, nd->center.y + nd->s/4); 
    
      
    nd->div = true; 
    *track = *track-4;
    // printf("Track: %i\n", *track);

}

//verifico che la particella sia contenuta nel nodo
bool contains(quadTree* nd, struct point p){
        return (
        p.x < (nd->center.x+nd->s/2) &&
        p.x > (nd->center.x-nd->s/2) &&
        p.y < (nd->center.y+nd->s/2) &&
        p.y > (nd->center.y-nd->s/2));
        
}
//costruzione albero in bottom up
int insert(quadTree* nd, particle* p, int *index, int* track){

    // se il quadrante non può contenerla
    if (!contains(nd,p->pos)){ 
        return 0; //particella non contenuta
    } 

    if(nd->p==NULL){ //assegno un puntatore dell' albero alla particella nel quadrante iteressato(teniamo la capacità del quadrante a 1)
        nd->p = p;
        nd->data = *index; //assegno il numero del corpo preso dall' array dei corpi, aiuta a ritornare i dati in caso il corpo sia in una foglia
        // printf("Pointer to %i\n", nd->data);
        return 0; // una volta trovato il corpo ritorniamo 0 e usciamo dalla ricorsiona
    } 
    else{
    
        if(nd->div!=true){ // controlliamo se il quadrante non isa stato diviso
            subdivide(nd, track); // se non è diviso chiamaiamo la funzione divide
            // controlla in quale figlio inserire la particella già presente 
            int temp = nd->data;
            nd->data = *track+4;
            *track = *track-1;
            insert(nd->ne, nd->p, &temp, track);  
            insert(nd->se, nd->p, &temp, track);  
            insert(nd->sw, nd->p, &temp, track); 
            insert(nd->nw, nd->p, &temp, track);
        }
    }

    return insert(nd->ne, p, index, track)|| // Since insert is an int function, the return statement here returns 0 if 
    insert(nd->se, p, index, track)||  // a node is found to point the body and thus save the body at. So for example if we
    insert(nd->sw, p, index, track)|| // have 2 bodies and the first one is always saved at root ( being empthy and not divided)
    insert(nd->nw, p, index, track); // then the second Body 1 if it is at NW subcell after division when it goes to the return state-
        // -ment it will go through the OR terms as the contain function will not let the first 3 OR terms ( being insert(nd->NE,...) || ... 
        // to insert(nd->SW,...) and pick the insert(nd->NW,...) to assign. The same process is repeated every time the function calls itself
        // and checks for where to put the body.
}




int main(){

}