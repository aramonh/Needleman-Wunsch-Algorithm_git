
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
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
#include <cstring>
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
 
#define MAX_SOURCE_SIZE (0x100000)
 #define MAX_N 1001
using namespace std;

string A,B;
int n, m;
int match_score, mismatch_score, gap_score;
int dp[MAX_N][MAX_N];



void print_matrix(int matrix[][MAX_N], int n, int m)
{
    for(int loop = 0; loop < n; loop++){
        for(int loop2 = 0; loop2 < m; loop2++){
            printf(" %d ", matrix[loop][loop2]);
        }
        printf("\n");
    }
}


inline int needleman_wunsch()
{
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


void nw(void) {

    match_score = 1, mismatch_score = 1, gap_score = 1; // Constants de score quemadas

    A = "CGATGCTAGCGTATCGTAGTCTATCGTAC";
    B = "ACGATGCTAGCGTTTCGTATCATCGTA";

    n = A.length(); // length of gene1
	  m = B.length(); // length of gene2

    // Create the two input vectors
    char *A_k = (char*)malloc(sizeof(char)*n);
    char *B_k = (char*)malloc(sizeof(char)*m);

    strcpy(A_k, A.c_str());
    strcpy(B_k, B.c_str());
 
    // Load the kernel source code into the array source_str
    FILE *fp;
    char *source_str;
    size_t source_size;
 
    fp = fopen("nw.cl", "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );
 
    // Get platform and device information
    cl_platform_id platform_id = NULL;
    cl_device_id device_id = NULL;   
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_DEFAULT, 1, 
            &device_id, &ret_num_devices);


   // Create an OpenCL context
    cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);
 
    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
 
    // Create memory buffers on the device for each vector 
    cl_mem a_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            n * sizeof(char), NULL, &ret);
    cl_mem b_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY,
            m * sizeof(char), NULL, &ret);
    cl_mem c_mem_obj = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
            n*m * sizeof(int), NULL, &ret);
 
    // Copy the lists A and B to their respective memory buffers
    ret = clEnqueueWriteBuffer(command_queue, a_mem_obj, CL_TRUE, 0,
            n * sizeof(char), A_k, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, b_mem_obj, CL_TRUE, 0, 
            m * sizeof(char), B_k, 0, NULL, NULL);

    // Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, 
            (const char **)&source_str, (const size_t *)&source_size, &ret);
 
    // Build the program
    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);

     // Create the OpenCL kernel
    cl_kernel kernel = clCreateKernel(program, "nw", &ret);
 
    // Set the arguments of the kernel
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&a_mem_obj);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&b_mem_obj);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&c_mem_obj);
    ret = clSetKernelArg(kernel, 3, sizeof(n), &n);
    ret = clSetKernelArg(kernel, 4, sizeof(m), &m);
    ret = clSetKernelArg(kernel, 5, sizeof(gap_score), &gap_score);

     // Execute the OpenCL kernel on the list
    size_t global_item_size = n*m; // Process the entire lists
    size_t local_item_size = 29; // Divide work items into groups of 64
    ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, 
            &global_item_size, &local_item_size, 0, NULL, NULL);
 
    // Read the memory buffer C on the device to the local variable C
    int *C = (int*)malloc(sizeof(int)*n*m);
    ret = clEnqueueReadBuffer(command_queue, c_mem_obj, CL_TRUE, 0, 
            n*m * sizeof(int), C, 0, NULL, NULL);

    int j,i = 0;
    for(int k = 0; k<(n*m) ; k++){
        if((k)%m==0){
            i = (k/m);
            j = 0;
        }else{
            j++;
        }
        dp[i][j] = C[k];
    }

    //print_matrix(dp,n,m);

    // Clean up
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(a_mem_obj);
    ret = clReleaseMemObject(b_mem_obj);
    ret = clReleaseMemObject(c_mem_obj);
    ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);
    free(A_k);
    free(B_k);
    free(C);

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