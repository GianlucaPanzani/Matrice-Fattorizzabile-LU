#include <stdio.h>
#include <stdlib.h>
#include <termios.h>
#include <string.h>

static struct termios term, oterm;
int REAL_DIM, LU = 0, PLU = 0;

int getch() {
    int c = 0;
    tcgetattr(0, &oterm);
    memcpy(&term, &oterm, sizeof(term));
    term.c_lflag &= ~(ICANON | ECHO);
    term.c_cc[VMIN] = 1;
    term.c_cc[VTIME] = 0;
    tcsetattr(0, TCSANOW, &term);
    c = getchar();
    tcsetattr(0, TCSANOW, &oterm);
    return c;
}

void stampaMatrice(char name, double** m, int N) {
    if(m == 0) {
        printf("error: stampa: Empty Matrix");
        return;
    }
    for(int i = 0; i < N; ++i) {
        if(i == N/2)
            printf("%c = |", name);
        else
            printf("    |");
        for(int j = 0; j < N; ++j)
            if(m[i][j] >= 0)
                printf("   %0.2f", m[i][j]);
            else
                printf("  %0.2f", m[i][j]);
        printf("   |\n");
    }
}

void stampaVettore(double* v, int N) {
    if(v == 0) {
        printf("error: stampa: Empty Matrix");
        return;
    }
    for(int i = 0; i < N; ++i) {
        if(i == N/2)
            printf("v = |");
        else
            printf("    |");
        if(v[i] >= 0)
            printf("  %0.2f  ", v[i]);
        else
            printf(" %0.2f  ", v[i]);
        printf("|\n");
    }
}

// sfrutta la TECNICA PIVOTING (o la TECNICA DEL MASSIMO PIVOT)
void shuffle(double** m, int N, int k) {
    LU--;
    PLU++;
    if(m == 0) {
        printf("error: stampa: Empty Matrix");
        return;
    }
    if(k == N-1) {
        REAL_DIM--;
        return;
    }

    int max = k+1;
    for(int i = k+1; i < N; ++i) {
        if(m[i][k] != 0 && m[max][k] < m[i][k])
            max = i;
    }

    if(m[max][k] == 0) {
        REAL_DIM--;
        return;
    }

    double d;
    for(int j = k; j < N; ++j) {
        d = m[max][j];
        m[max][j] = m[k][j];
        m[k][j] = d;
    }
}

int main(int argc, char* argv[]) { // eseguire con "./exe 1" per vedere i vari output durante l'esecuzione
    int x, N, flag = 0;

    if(argc > 1)
        flag = atoi(argv[1]);
    
    printf("\nCALCOLO MATRICE TRIANGOLARE SUPERIORE A PARTIRE DA MATRICE NXN\n\nInserire la dimensione della matrice N: ");
    scanf("%d", &N);
    if(N < 2) {
        printf("N > 1");
        return 0;
    }
    REAL_DIM = N;

    // A(k) = E(k)*A(k-1)
    // E(k) = I - v*e(k)
    double** a = (double**) malloc(N*sizeof(double*));
    double** L = (double**) malloc(N*sizeof(double*));
    printf("Esempio: inserimento di una matrice 2x2:\t  x y z w --> corrisponde alla matrice:\n x  y\n z  w\n\n");
    printf("Adesso inserisci gli elementi della matrice:\n");
    for(int i = 0; i < N; ++i) {
        a[i] = (double*) malloc(N*sizeof(double));
        L[i] = (double*) malloc(N*sizeof(double));
        for(int j = 0; j < N; ++j) {
            scanf("%d", &x);
            a[i][j] = (double) x;
            if(i == j)
                L[i][j] = 1;
            else
                L[i][j] = 0;
        }
    }
        
    int k = 0; // k --> passo di Gauss
    double* v = (double*) malloc(N*sizeof(double));
    while(k != N) {
        v[k] = 0;
        if(a[k][k] == 0)
            shuffle(a, N, k);
        for(int i = k+1; i < N; ++i) {
            v[i] = L[i][k] = a[i][k]/a[k][k];
            a[i][k] = 0;
        }
        for(int i = k+1; i < N; ++i) {
            for(int j = k+1; j < N; ++j) {
                a[i][j] = a[i][j]-(v[i]*a[k][j]);
            }
        }
        if(flag) {
            printf("\nMatrice A, Matrice L e Vettore v (al passo k = %d):\n", k);
            stampaMatrice('A', a, N);
            stampaMatrice('L', L, N);
            stampaVettore(v, N);
            printf("Premi qualsiasi tasto per continuare...\n");
            getch();
        }
        k++;
    }

    printf("\n Matrice TRIANGOLARE SUPERIORE ottenuta con GAUSS:\n");
    stampaMatrice('A', a, N);
    printf("\nProprieta' Matrice:\n");
    printf("\tdim(A) = %d\n", REAL_DIM);
    if(LU < 0 && REAL_DIM < N) {
        printf("\tSINGOLARE (NON INVERTIBILE)\n\tNON FATTORIZZABILE LU/PLU\n\t#SOLUZIONI = +inf\n");
        return 0;
    } else if(LU < 0 && REAL_DIM == N) {
        printf("\tNON SINGOLARE (INVERTIBILE)\n\tFATTORIZZABILE PLU\n\t#SOLUZIONI FINITE\n");
        return 0;
    } else {
        printf("\tNON SINGOLARE (INVERTIBILE)\n\tFATTORIZZABILE LU\n\t#SOLUZIONI FINITE\n\n");
    }
    stampaMatrice('L', L, N);
    
    
    // Possibile aggiunta:
    // calcolare l'inversa di L e trovare U cosi: U = L^(-1)*A

    // inversa di L
    double** Linv = (double**) malloc(N*sizeof(double*));
    for(int i = 0; i < N; ++i) {
        Linv[i] = (double*) malloc(N*sizeof(double));
        for(int j = 0; j < N; ++j)
            if(i == j)
                Linv[i][j] = 1;
            else
                Linv[i][j] = 0;
    }
    // to be continued....


    free(v);
    for(int i = 0; i < N; ++i) {
        free(a[i]);
        free(L[i]);
    }
    free(a);
    free(L);
    return 0;
}