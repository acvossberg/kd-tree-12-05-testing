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



__global__ void testKernel( int *trees, int numBlocks )
{
    int i = numBlocks*blockIdx.x + threadIdx.x;
    d_a[i] = blockIdx.x + threadIdx.x;
    
}

void cudaMain(){
    cudaSetDevice(MYDEVICE);



}
