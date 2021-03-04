/**************************************************************************************
*
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
*
* Progetto dell'algoritmo di Regressione
* in linguaggio assembly x86-32 + SSE
*
* Fabrizio Angiulli, aprile 2019
*
**************************************************************************************/

/*
*
* Software necessario per l'esecuzione:
*
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
*
* entrambi sono disponibili come pacchetti software
* installabili mediante il packaging tool del sistema
* operativo; per esempio, su Ubuntu, mediante i comandi:
*
*    sudo apt-get install nasm
*    sudo apt-get install gcc
*
* potrebbe essere necessario installare le seguenti librerie:
*
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
*
* Per generare il file eseguibile:
*
* nasm -f elf64 regression64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o regression64.o regression46c.c -o regression64c -lm && ./regression64c $pars
*
* oppure
*
* ./runregression64
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>
#include <omp.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*
#define eps         1e-8d

void somm_SGD(VECTOR stimaTheta, MATRIX Dbatch, VECTOR yBatch, int j, int t, VECTOR ret);
void calcolaTheta(MATRIX Dbatch, MATRIX G, int v, VECTOR yBatch, VECTOR stimaTheta, type eta, int t, int adagrad,VECTOR sommatoria);
void sommatoriaSGD(VECTOR stimaTheta, MATRIX Dbatch, VECTOR yBatch, int indexSommatoria, int t, int mt, VECTOR sommatoria);
void sommatoria_interna(VECTOR stimaTheta, MATRIX Dbatch, type G[], VECTOR yBatch, int j, int t, VECTOR sommatoria);
void sommatoriaAdagrad(VECTOR stimaTheta, MATRIX Dbatch, type G[], VECTOR yBatch, int indexSommatoria, int t, int mt, VECTOR sommatoria);

typedef struct {
    MATRIX x; //data set
    VECTOR y; //label set
    MATRIX xast; //data set convertito
    int n; //numero di punti del data set
    int d; //numero di dimensioni del data set
    int k; //dimensione del batch
    int degree; //grado del polinomio
    type eta; //learning rate
    //STRUTTURE OUTPUT
    VECTOR theta; //vettore dei parametri
    int t; //numero di parametri, dimensione del vettore theta
    int iter; //numero di iterazioni
    int adagrad; //accelerazione adagrad
    int silent; //silenzioso
    int display; //stampa risultati
} params;

void calcola_adagrad(params * input);
void calcola_sgd(params * input);
void calcolaTheta_adagrad(MATRIX Dbatch, MATRIX G, int v, VECTOR yBatch, VECTOR stimaTheta, type eta, int t, int mt, VECTOR sommatoria);
void calcolaTheta_sgd(MATRIX Dbatch, int v, VECTOR yBatch, VECTOR stimaTheta, type eta, int t, int mt, VECTOR sommatoria);
void numeroComb(int deg, int d, int *numComb, int *indice);


/*
*
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate
* 	mediante un array (double*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere
* 	memorizzate mediante array di array (double**).
*
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
*
*/

void* get_block(int size, int elements) {
    return _mm_malloc(elements*size,32);
}

void free_block(void* p) {
    _mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
    return (MATRIX) get_block(sizeof(type),rows*cols);
}

int* alloc_matrix_int(int rows, int cols) {
    return (int*) get_block(sizeof(int),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
    free_block(mat);
}

void dealloc_matrix_int(int* mat) {
    free_block(mat);
}

/*
*
* 	load_data
* 	=========
*
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
*
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
*
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice.
*****************************************************************************
*
*/
MATRIX load_data(char* filename, int *n, int *k) {
    FILE* fp;
    int rows, cols, status, i;

    fp = fopen(filename, "rb");

    if (fp == NULL){
        printf("'%s': bad data file name!\n", filename);
        exit(0);
    }

    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);

    MATRIX data = alloc_matrix(rows,cols);
    status = fread(data, sizeof(type), rows*cols, fp);
    fclose(fp);

    *n = rows;
    *k = cols;

    return data;
}

/*
* 	save_data
* 	=========
*
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
*
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
*/
void save_data(char* filename, void* X, int n, int k) {
    FILE* fp;
    int i;
    fp = fopen(filename, "wb");
    if(X != NULL){
        fwrite(&k, 4, 1, fp);
        fwrite(&n, 4, 1, fp);
        for (i = 0; i < n; i++) {
            fwrite(X, sizeof(type), k, fp);
            //printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
            X += sizeof(type)*k;
        }
    }
    fclose(fp);
}

// PROCEDURE ASSEMBLY
extern void prova(params* input);
/**
 * Inizializza il vettore x ricevuto come parametro al valore c.
 * @param x , vettore da inizializzare
 * @param n , lunghezza del vettore
 * @param c , valore con cui inizializzare il vettore
 */
extern void inizializza(VECTOR x, int n, type c);
extern void inserisciComb(int deg, int nCH, int degh, MATRIX x, int* mComb, MATRIX xastSomm); //assegna valori a x*
extern void init_padd(VECTOR ar, MATRIX x, int mt, int n); //inizializza la prima componente di ogni riga (sempre pari a 1) e le componenti di padding (poste a 0).
/**
 * La funzione calcola un vettore risultate a seguito del calcolo della seguente formula:
 * for(i from 1 to len)
 * vector_dest[i] = vector_dest[i] + fact_mul * vector_src[i];
 * @param fact_mul, fattore moltiplicativo
 * @param vector_src, vettore che funge da sorgente
 * @param vector_dest vettore che funge da destinazione
 * @param len, lunghezza del vettore
 */
extern void calcolaVettore(type fact_mul, VECTOR vector_src, VECTOR vector_dest, int len);
/**
 * La funzione calcola il prodotto scalare tra due vettori
 * @param vector_1, primo vettore
 * @param vector_2, secondo vettore
 * @param len ,lunghezza dei vettori
 * @param ret, indirizzo di memoria in cui caricare il risultato del prodotto scalare
 */
extern void prodottoScalare(VECTOR vector_1, VECTOR vector_2, int len, type* ret);
extern void sommatoriaInterna(type delta, VECTOR Dbatch,VECTOR G, int t , VECTOR sommatoria);

/*
* 	creaNum
* 	=========
*
*	crea il numero composto da h cifre uguali e pari a cifra
*	es: h=3 e cifra=4, creaNum restituisce 444
*/
int creaNum(int cifra, int h) {
    return pow(10, h)/9*cifra;
}
/**
 * Il metodo effettua il cambiamento di base di un numero i in base bp, e lo converte in base ba
 * @param i, numero di partenza
 * @param bp, base di partenza del numero i
 * @param ba, base in cui convertire il numero i
 */
int cambioBase(int i, int bp, int ba) {
	int ret=0, c=0;
	int divi=i;
	while (divi!=0) {
		ret += (divi % ba) * pow(bp, c);
		divi = divi/ba;
		c++;
	}
	return ret;
}

/*
 * combinazioni2
 * =============
 * genera la matrice delle combinazioni di grado massimo (deg)
 */
void combinazioni2(int* ris, int h, int d) {
	int i, j, divi;
	int cur=0, flag=0, start=0;

	for(i=0;!flag; i++, start++){
	    divi=start;
	    cur = divi % d; //valore della variabile
		divi/=d; //avanzamento alla variabile più a sinistra
	    if(cur == d-1)// se l'ultimo numero  uguale a d-1 al prossimo start ci sarà un numero non ordinato
	        flag=1;
		ris[i*h +h-1]=cur;
		for(j=h-2; j>=0; j--){
		    cur = divi % d;
			ris[i*h + j] = cur;
			divi/=d;
			if(flag && cur != d-1){ //il primo numero != d-1 � il numero a cui devono arrivare tutti gli 0 a destra di esso (e.g.1144->1200->+22)
        			start += cambioBase(creaNum(cur+1, h-j-1), d, 10); //viene creato il numero in base a quanti zeri (e.g. j=2) e cambiato di base per essere sommato a start (-1 perch� sopra start++)
        			flag=0;
			}
		}
	}
}

void convert_data(params* input){
    int deg=input->degree;
    int d=input->d;
    int n=input->n;
    MATRIX x=input->x;
    int t, h;
    int numComb[deg+1];
    int indice[deg+1];
    numeroComb(deg, d, numComb, indice);
    t= indice[deg];
    int mt = (t%4!=0)? t+ 4-t%4: t;
    MATRIX ret=alloc_matrix(n, mt);
    //inizializzazione matrice ret
    int delta = mt - t;
    ret[0]=1.0;
    if(delta != 0){
        type ar[4] = {0.0, 0.0, 0.0, 0.0};
        ar[delta]++;
        init_padd(ar, ret+t, mt, n);
        for(int i=t; i<mt; i++)
            ret[(n-1)*mt+i]=0.0;
    }
    else
        for (int i = 1; i < n; ++i)
            ret[i*mt]=1.0;

    int* mComb = alloc_matrix_int(numComb[deg], deg);
    combinazioni2(mComb, deg, d);
    int oss;

    for( h=1; h<=deg; h++) {
        #pragma omp parallel num_threads(3)
        {
        #pragma omp for
            for (oss = 0; oss < n; oss++)
                inserisciComb(deg, numComb[h], deg-h, x+oss*d, mComb, ret+indice[h-1]+oss*mt);
        }
        }
    input->t=t;
    input->xast=ret;
}

/*
 * numeroComb
 * ==========
 * calcola il numero di combinazioni per i dati d e deg e restituisce il numero di colonne di xast
 */
void numeroComb(int deg, int d, int *numComb, int *indice) {
    int nc=d;
    numComb[0]=1;
    numComb[1]=d;
    indice[0]=1;
    indice[1]=d+1;
    for(int h=2; h <=deg; h++){
        nc*=d+h-1;
        nc/=h;
        numComb[h]= nc;
        indice[h] = indice[h-1]+nc;
    }
}

void sgd(params * input){
    if(input->adagrad)
        calcola_adagrad(input);
    else
        calcola_sgd(input);
}


void calcola_sgd(params * input){
    MATRIX D=input->xast;
    int t = input->t;
    type eta = input->eta; // learning rate
    int k = input->k; // dimensione batch
    int iter = input->iter; //numero massimo di iterazioni
    VECTOR y = input->y; // label set
    int n = input->n; //numero righe della matrice D;
    int mt = (t%4!=0)? t+ 4-t%4: t;

    //estendo stima di theta e sommatoria
    VECTOR stimaTheta = alloc_matrix(mt, 1);
    inizializza(stimaTheta, mt, 0.0);
    int it = 0;
    int v;

    VECTOR sommatoria  = alloc_matrix(mt,1);
    while (it < iter){
        for(int i = 0 ; i < n; i+= k){
            v = (n-i>=k)? k:n-i;
            calcolaTheta_sgd((D + i * mt), v, (y + i), stimaTheta, eta, t, mt, sommatoria);
        }
        it++;
    }
    input->theta=stimaTheta;
}

void calcola_adagrad(params * input){
    MATRIX D=input->xast;
    int t = input->t;
    type eta = input->eta; // learning rate
    int k = input->k; // dimensione batch
    int iter = input->iter; //numero massimo di iterazioni
    type* y = input->y; // label set
    int n = input->n; //numero righe della matrice D;
    int mt = (t%4!=0)? t+ 4-t%4: t;

    MATRIX G = alloc_matrix(k, mt);
    inizializza(G, k * mt, eps);

    VECTOR stimaTheta = alloc_matrix(mt,1);
    inizializza(stimaTheta, mt, 0.0);
    int it = 0;
    int v;
    VECTOR sommatoria  = alloc_matrix(mt, 1);
    while (it < iter){
        for(int i = 0 ; i < n; i+= k){
            v = (n-i>=k)? k:n-i;
            calcolaTheta_adagrad((D + i * mt), G, v, (y + i), stimaTheta, eta, t, mt, sommatoria);
        }
        it++;
    }
   input->theta=stimaTheta;
}

/**
 * La funzione calcola il nuovo theta a partire dall'osservazione precedente
 * @param Dbatch, sottoinsieme della matrice D in cui sono prese solo le righe che rientrano nell'batch
 * quindi: puntatore iniziale + dimensione batch * indiceBloccoOsservazioni * osservazioniPerRiga
 * @param v, dimensione del batch corrente
 * @param yBatch, sottoinsieme della vettore y considerando il batch
 * @param stimaTheta, theta all'iterazione precedente
 * @param eta valore del learning rate
 * @param t, dimensione di theta
 */
void calcolaTheta_adagrad(MATRIX Dbatch, MATRIX G, int v, VECTOR yBatch, VECTOR stimaTheta, type eta, int t, int mt, VECTOR sommatoria) {
    sommatoriaAdagrad(stimaTheta, Dbatch, G, yBatch, v, t, mt, sommatoria);
    type mul =(eta * (-1))/ (type) v;
    calcolaVettore(mul, sommatoria, stimaTheta, t);
}


void calcolaTheta_sgd(MATRIX Dbatch, int v, VECTOR yBatch, VECTOR stimaTheta, type eta, int t, int mt, VECTOR sommatoria) {
    sommatoriaSGD(stimaTheta, Dbatch, yBatch, v, t, mt, sommatoria);
    type mul = (eta * (-1)) / (type) v ;
    calcolaVettore(mul, sommatoria, stimaTheta, t);
}

/**
 * Calcola il valore della sommatoria SGD
 * @param stimaTheta theta precedente
 * @param Dbatch, sottoinsieme della matrice D in cui sono prese solo le righe che rientrano nell'batch
 * @param yBatch, sottoinsieme del vettore y considerando il batch
 * @param indexSommatoria , indice corrente della sommatoria
 * @param t, dimensione di theta ed anche di una riga della matrice Dbatch
 */
void sommatoriaSGD(VECTOR stimaTheta, MATRIX Dbatch, VECTOR yBatch, int indexSommatoria, int t, int mt, VECTOR sommatoria) {
    inizializza(sommatoria, mt, 0.0);
    for(int j= 0; j< indexSommatoria; j++ ){
        type ps;
        prodottoScalare(stimaTheta, (Dbatch + j*mt), t, &ps);
        type delta = ps - yBatch[j];
        calcolaVettore(delta,(Dbatch + j * mt), sommatoria, t);
    }
}

/**
 * La funzione calcola il valore della sommatoria all'interno della funzione adagrad
 * @param stimaTheta , stima di theta precedente
 * @param Dbatch ,
 * @param G
 * @param yBatch
 * @param indexSommatoria
 * @param t
 * @param sommatoria
 */
void sommatoriaAdagrad(VECTOR stimaTheta, MATRIX Dbatch, type G[], VECTOR yBatch, int indexSommatoria, int t, int mt, VECTOR sommatoria) {
    inizializza(sommatoria, mt, 0.0);
    for(int j= 0; j< indexSommatoria; j++ ){
        type gjv;
        type ps;
        prodottoScalare(stimaTheta, (Dbatch + j*mt), t, &ps);
        type delta = ps - yBatch[j];
        sommatoriaInterna(delta, (Dbatch + j * mt) ,(G + j * mt), t, sommatoria);
    }
}

int main(int argc, char** argv) {

    char fname[256];
    char* dsname;
    char* filename;
    int i, j, k;
    clock_t t;
    float time;
    int yd = 1;

    //
    // Imposta i valori di default dei parametri
    //

    params* input = malloc(sizeof(params));

    input->x = NULL;
    input->y = NULL;
    input->xast = NULL;
    input->n = 0;
    input->d = 0;
    input->k = -1;
    input->degree = -1;
    input->eta = -1;
    input->iter = -1;
    input->adagrad = 0;
    input->theta = NULL;
    input->t = 0;
    input->adagrad = 0;
    input->silent = 0;
    input->display = 0;

    //
    // Visualizza la sintassi del passaggio dei parametri da riga comandi
    //

    if(argc <= 1){
        printf("%s D -batch <k> -degree <deg> -eta <eta> -iter <it> [-adagrad]\n", argv[0]);
        printf("\nParameters:\n");
        printf("\tD: il nome del file, estensione .data per i dati x, estensione .labels per le etichette y\n");
        printf("\t-batch <k>: il numero di campini nel batch\n");
        printf("\t-degree <deg>: il grado del polinomio\n");
        printf("\t-eta <eta>: il learning rate\n");
        printf("\t-iter <it>: il numero di iterazioni\n");
        printf("\t-adagrad: l'acceleratore AdaGrad\n");
        exit(0);
    }

    //
    // Legge i valori dei parametri da riga comandi
    //

    int par = 1;
    while (par < argc) {
        if (par == 1) {
            filename = argv[par];
            par++;
        } else if (strcmp(argv[par],"-s") == 0) {
            input->silent = 1;
            par++;
        } else if (strcmp(argv[par],"-d") == 0) {
            input->display = 1;
            par++;
        } else if (strcmp(argv[par],"-batch") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing batch dimension value!\n");
                exit(1);
            }
            input->k = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-degree") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing degree value!\n");
                exit(1);
            }
            input->degree = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-eta") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing eta value!\n");
                exit(1);
            }
            input->eta = atof(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-iter") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing iter value!\n");
                exit(1);
            }
            input->iter = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-adagrad") == 0) {
            input->adagrad = 1;
            par++;
        } else{
            printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
            par++;
        }
    }

    //
    // Legge i dati e verifica la correttezza dei parametri
    //

    if(filename == NULL || strlen(filename) == 0){
        printf("Missing input file name!\n");
        exit(1);
    }

    dsname = basename(strdup(filename));
    sprintf(fname, "%s.data", filename);
    input->x = load_data(fname, &input->n, &input->d);


    sprintf(fname, "%s.labels", filename);
    input->y = load_data(fname, &input->n, &yd);

    if(input->k < 0){
        printf("Invalid value of batch dimension parameter!\n");
        exit(1);
    }

    if(input->degree < 0){
        printf("Invalid value of degree parameter!\n");
        exit(1);
    }

    if(input->eta < 0){
        printf("Invalid value of eta parameter!\n");
        exit(1);
    }

    if(input->iter < 0){
        printf("Invalid value of iter parameter!\n");
        exit(1);
    }

    //
    // Visualizza il valore dei parametri
    //

    if(!input->silent){
        printf("Input data name: '%s.data'\n", filename);
        printf("Input label name: '%s.labels'\n", filename);
        printf("Data set size [n]: %d\n", input->n);
        printf("Number of dimensions [d]: %d\n", input->d);
        printf("Batch dimension: %d\n", input->k);
        printf("Degree: %d\n", input->degree);
        printf("Eta: %f\n", input->eta);
        if(input->adagrad)
            printf("Adagrad enabled\n");
        else
            printf("Adagrad disabled\n");
    }
    //
    // Conversione Dati
    //

    t = clock();
    convert_data(input);
    t = clock() - t;
    time = ((float)t)/CLOCKS_PER_SEC;

    sprintf(fname, "%s.xast", dsname);

    if(!input->silent)
        printf("Conversion time = %f secs\n", time);
    else
        printf("%.3f\n", time);

    //
    // Regressione
    //

    t = clock();
    sgd(input);
    t = clock() - t;
    time = ((float)t)/CLOCKS_PER_SEC;

    if(!input->silent) {
        printf("Regression time = %f secs\n", time);
    }
    else {
        printf("%.3f\n", time);
    }

    //
    // Salva il risultato di theta
    //

    if(!input->adagrad) {
        sprintf(fname, "%s.theta.sgdomp", dsname);
    }
    else
        sprintf(fname, "%s.theta.adagradomp", dsname);
    save_data(fname, input->theta, input->t, 1);

    if(input->display){
        printf("theta: [");
        for(i=0; i<input->t-1; i++)
            printf("%f,", input->theta[i]);
        printf("%f]\n", input->theta[i]);
    }
    if(!input->silent)
        printf("\nDone.\n");
    return 0;
}