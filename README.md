<img src="https://user-images.githubusercontent.com/63518174/210822356-9cfed022-b237-43fe-af06-ab5266f65b31.png" alt="drawing" width="200"/>

Demetryo Grad

Matteo Bernardo

Silvio Gori

# Problema N-Corpi


Il problema N-corpi o "N-body problem" dall' inglese, è uno dei problemi più incisivi della fisica matematica, basato sulla seconda legge del moto e la legge di gravitazione universale stipulate da Newton.

Viene posto il quesito del calcolo del moto di **n** particelle, soggette alla legge di gravitazione universale date: le loro masse, velocità e posizione nello spazio e il loro mutare nel tempo.

## Metodo Naive

Utilizzando la formula di gravitazione universale unita alla seconda legge del moto  e definendo 2 particelle **q** e **k** tali che:

- **m<b><sub>q </sub></b>** **m<b><sub>k </sub></b>** siano le loro rispettive masse;
- **s<b><sub>q </sub></b>** **s<b><sub>k </sub></b>**  la loro posizione nello spazio;
- **t** un determinato istante di tempo;
- ed essendo **G** la costante di gravitazione universale;

la forza **f<b><sub>qk </sub></b>** che esercita **k** sulla particella **q** è data da:

$$f_{qk} = {G * m_q * m_k \over |s_q( t)-s_k( t)|^3} *{[s_q( t)s_k( t)]}$$
Per calcolare le forze che tutte le particelle esercitsno su **q** , definito come: **F<b><sub>q </sub></b>**, possiamo usare la sommatoria:
	 $$F_q = \sum^{n-1}_{k=0 \space\land \space k!=q}{f_qk}$$
applicando tale sommatoria a tutte le particelle si può creare un **algoritmo** che calcola ogni volta le forze relative alle particelle, viene definito come approccio **naive**, ergo calcolando in bruteforce tutte le forze delle particelle e successivamente la posizione e velocità relativa. Ciò aumenta notevolmente i tempi di simulazione specialmente per grandi numeri di particelle, avendo un costo computazionale elevato.


## Algoritmo di Barnes-Hut

Un altro approccio per ridurre drasticamente il numero di calcoli che deve svolgere il programma è quello di utilizzare l'**algoritmo di Barnes-Hut**. 

In queste simulazioni verrà usato un albero a **4 nodi** per radice, detto **quadTree** e le simulazioni verranno eseguite su piano bidimensionale.

l' algoritmo cponsiste nel dividere la regione di spazio dove sono contenute le particelle in 4 parti, continuare ricorsivamente la divisione fino a raggiungere la capacità di una particella per quandrante.

Si procede calcolando il centro di massa di ogni quadrante dell' albero partendo dalle foglie fino a alla radice, per ogni centro di massa **cm**.
Ogni **cm** è coposto da un valore che indica la mass e 2 valori che indicano la posizione del suo epicentro (nel caso di una simulazione bidimensionale), per ogni quadrante verrà calcolato come:

- in caso ci sia una particella sola nel quadrante ergo nel caso si arrivi ad una foglia : il **cm** prende come parametri la massa e i valori di posizione della particella stessa.
   
- in caso si stesse esplorando una radice e quindi ci siano più cm per quandrnte, perchè avente figli:

  -  per determinare la massa, vengono presi i **cm** deggli **n** figli e sommati tra loro.
    
  -  per ottenerne posizioni di esse:
  si fa la somma del prodotto tra le masse (**mass<b><sub>f0-(n-1)</sub></b>**) e posizioni di tutti i figli(**X<b><sub>f0-(n-1) </sub></b>** ) e si divide per la massa totale del ** cm ** (**mass<b><sub>cm </sub></b>**).
  Per esempio, in una griglia bidimensionale per trovare la x del epicentro del cm ( **X<b><sub>cm </sub></b>** )si può usare la formula:
	
$$X_{cm} = ({\sum^{i=0}_{n-1}mass_{fi} * X_{fi}})\space/\space{mass_{cm}}$$

Una volta calcolati i centri di massa si passa al calcolo delle forze che agiscono sulle particelle, si esamina l' albero dalla radice alle foglie, guardando prima i figli ( post-order traversal) e per ridurre i calcoli si passa tutto a un filtro, una costante **θ**, di solito dal valore compreso tra 0 e 1.5, che determina fino a quale grado avviene l' approssimazione:
- se **θ** = 0, la velocità dell' algoritmo sarà uguale a quello dell' algoritmo naive.
-  Per valori **θ** maggiori di 1 si inizia a vedere l' impatto dell' approssimazione sulla velocità di calcolo anche su poche particelle, aumentando naturalmente l' errore d' approssimazione sui vari calcoli.
  
Il calcolo delle forze delle partielle procede con le stesse formule del metodo naive, ma riducendo i calcoli usando la massa dei centri di massa al posto delle varie masse delle altre particelle, i centri di massa da utilizzare vengono scelti in base al filtro del **θ**, riducendo l' accesso all' albero e quindi la quantità di calcoli per paritcella, per ogni particella.


# Implementazione seriale

## Metodo Naive



## Algoritmo di Barnes-Hut




## Fonti:
