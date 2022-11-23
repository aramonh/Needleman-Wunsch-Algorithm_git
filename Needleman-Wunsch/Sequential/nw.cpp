/*
  Needleman-Wunsch algorithm  
  by Adriano Ramón Hernández
  Universidad Nacional de Colombia
*/

#include <stdio.h>
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


#define MAX_N 10001

using namespace std;

int n, m;
int match_score, mismatch_score, gap_score;
string A, B;
int dp[MAX_N][MAX_N];

/*
 Needleman-Wunsch algorithm para determinar la alineación óptima entre dos secuencias
 suponiendo una puntuación dada por matches , gaps and mismatches.
 Complexity: O(n * m) time, O(n * m) memory
*/

void print_matrix(int matrix[][MAX_N], int n, int m)
{
    for(int loop = 0; loop < n; loop++){
        for(int loop2 = 0; loop2 < m; loop2++){
            printf(" %d ", matrix[loop][loop2]);
        }
        printf("\n");
    }
}

/* 
Rellena la matriz segun las 2 secuncial iniciales, 
Acorde a los score en matchs, mismatches y gaps
*/
inline int needleman_wunsch()
{
    
    for (int i=0;i<=n;i++) dp[i][0] = dp[0][i] = -i * gap_score;

    for (int j=-m ; j<=n ;j++)
    {
        for (int i=max(1,j) ; i<=min(n, j+m ) ; i++)
        {
          
            int j1 = i-j;
            int S = (A[i-1] == B[j1-1]) ? match_score : -mismatch_score;
            dp[i][j1] = max(dp[i-1][j1-1] + S, max(dp[i-1][j1] - gap_score, dp[i][j1-1] - gap_score));
        }
    }

    return dp[n][m];
}

/*
Recorre la matriz generada,
buscando el alineamiento mas optimo
y generando las nuevas secuencias.
*/
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
            int S = (A[ii-1] == B[jj-1]) ? match_score : -mismatch_score;
            if (dp[ii][jj] == dp[ii-1][jj-1] + S)
            {
                SA.push(A[ii-1]);
                SB.push(B[jj-1]);
                ii--; jj--;
            }
            else if (dp[ii-1][jj] > dp[ii][jj-1])
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

inline void nw(){

    match_score = 1, mismatch_score = 1, gap_score = 1; // Constants de score quemadas

    A = "CGATGCTAGCGTATCGTAGTCTATCGTAC";
    B = "ACGATGCTAGCGTTTCGTATCATCGTA";

    n = A.length(); // length of gene1
  	m = B.length(); // length of gene2
	
    needleman_wunsch();
    pair<string, string> alignment = get_optimal_alignment();

    printf("\n %s \n %s \n", alignment.first.c_str(), alignment.second.c_str());


}

int main()
{

    clock_t start = clock();
    nw();
    printf("Time elapsed: %.6fs \n", (double)(clock() - start) / CLOCKS_PER_SEC);

    return 0;
}