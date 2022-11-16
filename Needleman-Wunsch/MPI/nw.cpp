#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>

#include <mpi.h>
#define TAM_MAX 20000
#define ROOT_ID 0
using namespace std;

void inicializacion(char primeraSecuencia[], char segundaSecuencia[]);
void matrizDeScore(char primeraSecuencia[], char segundaSecuencia[], int id);
void printMatriz(char primeraSecuencia[], char segundaSecuencia[]);
int MAYOR(int a, int b);

int matriz[TAM_MAX][TAM_MAX];
int match = 1;
int missmatch = 1;
int gap = 1;
string A,B;
int n, m;

void print_matrix(int matrix[][TAM_MAX], int n, int m)
{
    for(int loop = 0; loop < n; loop++){
        for(int loop2 = 0; loop2 < m; loop2++){
            printf(" %d ", matrix[loop][loop2]);
        }
        printf("\n");
    }
}

inline pair<string, string> get_optimal_alignment()
{
    string retA, retB;
    stack<char> SA, SB;
    int ii = n, jj = m;
    while (ii != 0 || jj != 0)
    {
        if (ii == 0)
        {
            SA.push('-');
            SB.push(B[jj-1]);
            jj--;
        }
        else if (jj == 0)
        {
            SA.push(A[ii-1]);
            SB.push('-');
            ii--;
        }
        else
        {
            int S = (A[ii-1] == B[jj-1]) ? match : -missmatch;
            if (matriz[ii][jj] == matriz[ii-1][jj-1] + S)
            {
                SA.push(A[ii-1]);
                SB.push(B[jj-1]);
                ii--; jj--;
            }
            else if (matriz[ii-1][jj] > matriz[ii][jj-1])
            {
                SA.push(A[ii-1]);
                SB.push('-');
                ii--;
            }
            else
            {
                SA.push('-');
                SB.push(B[jj-1]);
                jj--;
            }
        }
    }
    while (!SA.empty())
    {
        retA += SA.top();
        retB += SB.top();
        SA.pop();
        SB.pop();
    }
    return make_pair(retA, retB);
}



int main(){	
	double initialTime = MPI_Wtime();

	int id, num_processos, i;

	MPI_Init(NULL, NULL);
	MPI_Status status;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processos);

	printf("Proceso %d entre %d\n", id, num_processos);

	FILE *arq;
	char primeraSecuencia[TAM_MAX]; 
    char segundaSecuencia[TAM_MAX];
    
	arq = fopen("./input1.txt", "rt");
	if (arq == NULL) {
		printf("Problemas abriendo archivo\n");
	 	return 0;
	}
	
	while (!feof(arq)){	  
  		fgets(primeraSecuencia, TAM_MAX, arq);
  		primeraSecuencia[strcspn(primeraSecuencia, "\n")] = 0;	  
	}
	
	fclose(arq);
	
	arq = fopen("./input2.txt", "rt");
	if (arq == NULL) {
		printf("Problemas abriendo archivo\n");
	 	return 0;
	}
	
	while (!feof(arq)){	  
  		fgets(segundaSecuencia, TAM_MAX, arq);
	  	segundaSecuencia[strcspn(segundaSecuencia, "\n")] = 0;
	}
	
	fclose(arq);

	if (id == ROOT_ID) {
		printf("Primera Sequencia: %s \n", primeraSecuencia);
    	printf("Segunda Sequencia: %s \n", segundaSecuencia);
        A = primeraSecuencia;
        B = segundaSecuencia;
        n = A.length(); // length of gene1
        m = B.length(); // length of gene2
	
	}	
	
    inicializacion(primeraSecuencia, segundaSecuencia);

	MPI_Barrier(MPI_COMM_WORLD);

    matrizDeScore(primeraSecuencia, segundaSecuencia, id);

	if (id == ROOT_ID) {
		//printMatriz(primeraSecuencia, segundaSecuencia);	
        //print_matrix(matriz, n , m );
        pair<string, string> alignment = get_optimal_alignment();

        printf("\n %s \n %s \n", alignment.first.c_str(), alignment.second.c_str());

	}

	MPI_Barrier(MPI_COMM_WORLD);
	double endTime = MPI_Wtime();
	
	if (id == ROOT_ID) {
		printf("Time elapsed: %fs\n", endTime - initialTime);
	}
	
	MPI_Finalize();

    
	return 0;
}

void inicializacion(char primeraSecuencia[], char segundaSecuencia[]) {
    int sizeprimeraSecuencia = strlen(primeraSecuencia);
    int sizesegundaSecuencia = strlen(segundaSecuencia);
    
	matriz[0][0] = 0;
    
	for (int i = 0; i <= sizeprimeraSecuencia ; i++) {    	
        matriz[i][0] = -i*(gap); 
	}	
		
	for (int j = 0; j <= sizesegundaSecuencia ; j++) {    	
		matriz[0][j] = -j*(gap);
	}
}

void matrizDeScore(char primeraSecuencia[], char segundaSecuencia [], int id) {
    int sizeprimeraSecuencia = strlen(primeraSecuencia);
    int sizesegundaSecuencia = strlen(segundaSecuencia);
	MPI_Status status;

	int num_processos;

	MPI_Comm_size(MPI_COMM_WORLD, &num_processos);
    
    for (int i = id + 1; i < sizeprimeraSecuencia + 1; i++) {
        for (int j = 1; j < sizesegundaSecuencia + 1; j++) {			
			if (id != ROOT_ID) {
				int linhaRecebida[TAM_MAX];
				MPI_Recv(&linhaRecebida, sizeprimeraSecuencia, MPI_INT, id - 1, i - 1, MPI_COMM_WORLD, &status);				
				memcpy(&matriz[i - 1], &linhaRecebida, sizeof(int)*TAM_MAX);
			}
        
            int valorDiagonal = 0;
			
            if (primeraSecuencia[i-1] == segundaSecuencia[j-1]) {
			    valorDiagonal = matriz[i - 1][j - 1] + match;
            } else {
			    valorDiagonal = matriz[i - 1][j - 1] + (-missmatch);
            }
    
            int valorIzq = matriz[i][j - 1] - gap;
            int valorCima =  matriz[i - 1][j] - gap;
            int maximoScore = MAYOR(valorDiagonal, MAYOR( valorIzq, valorCima) );
            
            matriz[i][j] = maximoScore;

			if (id != num_processos - 1) {
				MPI_Send(&matriz[i], sizeprimeraSecuencia, MPI_INT, id + 1, i, MPI_COMM_WORLD);
			}
        }
		if (id == ROOT_ID) {
			for (int k = 1; k < num_processos; k++) {

				if (i + k < sizeprimeraSecuencia) {
					int linhaRecebida[TAM_MAX];
					MPI_Recv(&linhaRecebida, sizeprimeraSecuencia + 1, MPI_INT, id + k, i + k, MPI_COMM_WORLD, &status);					
					memcpy(&matriz[i + k], &linhaRecebida, sizeof(int)*TAM_MAX);					
				}
			}
		} else {
			MPI_Send(&matriz[i], sizeprimeraSecuencia + 1, MPI_INT, ROOT_ID, i, MPI_COMM_WORLD);
		}
    }
}

int MAYOR (int a, int b) {
	if (a > b) {
		return a;
	}
	
	return b;
}

void printMatriz(char primeraSecuencia[], char segundaSecuencia []){
    int sizeprimeraSecuencia = strlen(primeraSecuencia);
    int sizesegundaSecuencia = strlen(segundaSecuencia);
    
    printf("\t\t");
    for (int i = 0; i < sizesegundaSecuencia; i++) {		    
		printf("%c\t", segundaSecuencia[i]);
	}
    
    for (int i = 0; i < sizeprimeraSecuencia + 1; i++){
        printf("\n");
        
        if (i > 0) {
        	printf("%c\t", primeraSecuencia[i - 1]);	
		} else {
			printf("\t");
		}
		
        for (int j = 0; j < sizesegundaSecuencia + 1; j++){
        	if (matriz[i][j] < 0) {
        		printf("%d\t",matriz[i][j]);	
			} else {
				printf(" %d\t",matriz[i][j]);
			}
            
        }

    }
	printf("\n");
}

