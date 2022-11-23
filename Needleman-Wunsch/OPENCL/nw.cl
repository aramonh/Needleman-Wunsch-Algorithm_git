__kernel void nw(__global char *A, __global char *B, __global int *C, int n , int m , int gap_score ) {
 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    if(i<=(n*m)){
        if(i<=(m-1) ){
            C[i] = -i * gap_score;
        }else{
            if((i)%m==0){
               C[i] = (-i/m) * gap_score; 
            }
        }
    } 
}