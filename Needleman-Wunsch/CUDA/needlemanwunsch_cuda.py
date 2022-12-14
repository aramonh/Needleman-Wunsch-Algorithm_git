# -*- coding: utf-8 -*-
"""needlemanWunsch_CUDA.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1IMxSu_LYkuPzHkDrR1DQO7XopHSuxfoF
"""

!apt-get --purge remove cuda nvidia* libnvidia-*
!dpkg -l | grep cuda- | awk '{print $2}' | xargs -n1 dpkg --purge
!apt-get remove cuda-*
!apt autoremove
!apt-get update

!wget  --no-clobber https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-repo-ubuntu1804_10.0.130-1_amd64.deb
#install CUDA kit dpkg
!dpkg -i cuda-repo-ubuntu1804_10.0.130-1_amd64.deb
!sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
!apt-get update
!apt-get install cuda-10-0

# Commented out IPython magic to ensure Python compatibility.
!nvcc --version
!pip install git+https://github.com/andreinechaev/nvcc4jupyter.git
# %load_ext nvcc_plugin

# Commented out IPython magic to ensure Python compatibility.
# #@title NeedlemanWunsch
# 
# %%cu
# /************* NeedlemanWunsch ******************************************************/
# #include <stdio.h>
# // For the CUDA runtime routines (prefixed with "cuda_")
# #include <cuda.h>
# #include <cuda_runtime.h>
# #include <stdio.h>
# #include <math.h>
# #include <string.h>
# #include <time.h>
# #include <iostream>
# #include <vector>
# #include <list>
# #include <string>
# #include <algorithm>
# #include <queue>
# #include <stack>
# #include <set>
# #include <map>
# #include <complex>
# #include <cstring>
# #define MAX_N 1001
# #include <string.h>
# using namespace std;
# // Threads per CTA dimension
# int THREADS = 1024;
# 
# // Blocks per grid dimension (assumes THREADS divides N evenly)
# int BLOCKS = 1;
# int n, m;
# int match_score, mismatch_score, gap_score;
# string A, B;
# int dp[MAX_N][MAX_N];
# 
# void print_matrix(int matrix[][MAX_N], int n, int m)
# {
#     for(int loop = 0; loop < n; loop++){
#         for(int loop2 = 0; loop2 < m; loop2++){
#             printf(" %d ", matrix[loop][loop2]);
#         }
#         printf("\n");
#     }
# }
# 
# 
# inline int needleman_wunsch()
# {
#     for (int j=-m ; j<=n ;j++)
#     {
#         for (int i=max(1,j) ; i<=min(n, j+m ) ; i++)
#         {
#           
#             int j1 = i-j;
#             int S = (A[i-1] == B[j1-1]) ? match_score : -mismatch_score;
#             dp[i][j1] = max(dp[i-1][j1-1] + S, max(dp[i-1][j1] - gap_score, dp[i][j1-1] - gap_score));
#         }
#     }
# 
#     return dp[n][m];
# }
# 
# inline pair<string, string> get_optimal_alignment()
# {
#     string retA, retB;
#     stack<char> SA, SB;
#     int ii = n, jj = m;
#     while (ii != 0 || jj != 0)
#     {
#         if (ii == 0)
#         {
#             SA.push('-');
#             SB.push(B[jj-1]);
#             jj--;
#         }
#         else if (jj == 0)
#         {
#             SA.push(A[ii-1]);
#             SB.push('-');
#             ii--;
#         }
#         else
#         {
#             int S = (A[ii-1] == B[jj-1]) ? match_score : -mismatch_score;
#             if (dp[ii][jj] == dp[ii-1][jj-1] + S)
#             {
#                 SA.push(A[ii-1]);
#                 SB.push(B[jj-1]);
#                 ii--; jj--;
#             }
#             else if (dp[ii-1][jj] > dp[ii][jj-1])
#             {
#                 SA.push(A[ii-1]);
#                 SB.push('-');
#                 ii--;
#             }
#             else
#             {
#                 SA.push('-');
#                 SB.push(B[jj-1]);
#                 jj--;
#             }
#         }
#     }
#     while (!SA.empty())
#     {
#         retA += SA.top();
#         retB += SB.top();
#         SA.pop();
#         SB.pop();
#     }
#     return make_pair(retA, retB);
# }
#  __global__ void initializate(const char *A,const char *B, int n , int m , int *C , int gap_score){
#     int i = blockDim.x * blockIdx.x + threadIdx.x;
#     if(i<=(n*m)){
#         if(i<=(m-1) ){
#             C[i] = -i * gap_score;
#         }else{
#             if((i)%m==0){
#                C[i] = (-i/m) * gap_score; 
#             }
#         }
#     }    
# }
#  __global__ void fill(const char *A,const char *B, int n , int m , int *C , int match_score , int mismatch_score){
#     // Compute each thread's global row and column index
#     int k = blockDim.x * blockIdx.x + threadIdx.x;
#     int j = 0;
#     int i = 0;
#     if(k<n*m){
#         //dp[i][j1] = max(dp[i-1][j1-1] + S, max(dp[i-1][j1] - gap_score, dp[i][j1-1] - gap_score)); 
#         if((k)%m==0){
#             i = (k/m);
#             j = 0;
#         }else{
#             j++;
#         }
#         int S = (A[i-1] == B[i-1]) ? match_score : -mismatch_score;
#     }
# 
# }
# 
# inline void nw(){
#     
#     cudaError_t err = cudaSuccess;
# 
#     match_score = 1, mismatch_score = 1, gap_score = 1; // Constants de score quemadas
# 
#     A = "AGGGCT";
#     B = "AGGCA";
#     //A = "CGATGCTAGCGTATCGTAGTCTATCGTAC";
#     //B = "ACGATGCTAGCGTTTCGTATCATCGTA";
# 
#     n = A.length(); // length of gene1
# 	  m = B.length(); // length of gene2
# 
#     char *h_A, *h_B;
#     int *h_C;
# 
#     size_t size_A = sizeof(char)*n;
#     size_t size_B = sizeof(char)*m;
#     size_t size_C = sizeof(int)*m*n;
# 
#     h_A = (char*) malloc(size_A);
#     h_B = (char*) malloc(size_B);
#     h_C = (int*) malloc(size_C);
# 
#     h_A = "AGGGCT";
#     h_B = "AGGCA";
# 
#     // Allocate the device input vector A
#     char *d_A = NULL;
#     err = cudaMalloc((void **)&d_A, size_A);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to allocate device vector A (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     // Allocate the device input vector B
#     char *d_B = NULL;
#     err = cudaMalloc((void **)&d_B, size_B);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to allocate device vector B (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     // Allocate the device input vector C
#     int *d_C = NULL;
#     err = cudaMalloc((void **)&d_C, size_C);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to allocate device vector C (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     err = cudaMemcpy(d_A, h_A, size_A, cudaMemcpyHostToDevice);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to copy vector A from host to device (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     err = cudaMemcpy(d_B, h_B, size_B, cudaMemcpyHostToDevice);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to copy vector B from host to device (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
# 
#     initializate<<<BLOCKS, THREADS>>>(d_A,d_B,n,m,d_C, gap_score);
#     // Wait for GPU to finish before accessing on host
#     cudaDeviceSynchronize();
#     err = cudaGetLastError();
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     fill<<<BLOCKS, THREADS>>>(d_A,d_B,n,m,d_C , match_score, mismatch_score);
#     // Wait for GPU to finish before accessing on host
#     cudaDeviceSynchronize();
#     err = cudaGetLastError();
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     err = cudaMemcpy(h_C, d_C, size_C, cudaMemcpyDeviceToHost);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to copy vector C from device to host (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     // Free device global memory
#     err = cudaFree(d_A);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to free device vector A (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
#     err = cudaFree(d_B);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to free device vector B (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     err = cudaFree(d_C);
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to free device vector C (error code %s)!\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
# 
#     // Reset the device and exit
#     err = cudaDeviceReset();
# 
#     if (err != cudaSuccess)
#     {
#         fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
#         exit(EXIT_FAILURE);
#     }
# 
#     int j,i = 0;
#     for(int k = 0; k<(n*m) ; k++){
#         if((k)%m==0){
#             i = (k/m);
#             j = 0;
#         }else{
#             j++;
#         }
#         dp[i][j] = h_C[k];
#     }
# 
# 
# 
# 
#     needleman_wunsch();
#     //print_matrix(dp,n,m);
#     pair<string, string> alignment = get_optimal_alignment();
# 
#     printf("\n %s \n %s \n", alignment.first.c_str(), alignment.second.c_str());
# 
# }
# 
# 
# 
# int main()
# {
# 
#     clock_t start = clock();
# 
#     nw();
# 
#     printf("Time elapsed: %.6fs \n", (double)(clock() - start) / CLOCKS_PER_SEC);
# 
#     return 0;
# }