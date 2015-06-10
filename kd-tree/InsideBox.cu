//
//  InsideBox.cu
//  
//
//  Created by Ann-Christine Vossberg on 6/3/15.
//
//

#include <stdio.h>


// Here you can set the device ID that was assigned to you
#define MYDEVICE 0

// Simple utility function to check for CUDA runtime errors
void checkCUDAError(const char *msg);



__global__ void insideBox( int *trees, int numBlocks )
{
    int warpSize = 32;
    int warpIdx = threadIdx.x / warpSize;
    int i = warpIdx;
    int j = threadIdx.x % 32; //=0 bis 32
    
    if( ((trees[i][j].x >= start.x && trees[i][j].x <= end.x) || (end.x == 0 && start.x == 0))  && ((trees[i][j].y >= start.y && trees[i][j].y <= end.y) || (end.y == 0 && start.y == 0)) && ((trees[i][j].z >= start.z && trees[i][j].z <= end.z) || (end.z == 0 && start.z == 0))){
    }
    
    
    
    int i = numBlocks*blockIdx.x + threadIdx.x;
    d_a[i] = blockIdx.x + threadIdx.x;
    
}

void cudaMain(){
    cudaSetDevice(MYDEVICE);



}
