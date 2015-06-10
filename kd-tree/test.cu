// includes, system
#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include "test.hpp"

// Here you can set the device ID that was assigned to you
#define MYDEVICE 0

// Simple utility function to check for CUDA runtime errors
void checkCUDAError(const char *msg);

// Part 3 of 5: implement the kernel
__global__ void myFirstKernel( int *d_a, int numBlocks )
{
    int i = numBlocks*blockIdx.x + threadIdx.x;
    d_a[i] = blockIdx.x + threadIdx.x;
    
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
void test( )
{
    cudaSetDevice(MYDEVICE);
    // pointer for host memory
    int *h_a;
    
    // pointer for device memory
    int *d_a;
    
    // define grid and block size
    int numBlocks = 8;
    int numThreadsPerBlock = 8;
    
    // Part 1 of 5: allocate host and device memory
    size_t memSize = numBlocks * numThreadsPerBlock * sizeof(int);
    h_a = (int *) malloc(memSize);
    cudaMalloc((void **)&d_a, memSize);
    
    // Part 2 of 5: configure and launch kernel
    dim3 dimGrid( 8 );
    dim3 dimBlock( 8 );
    myFirstKernel<<< 8 , 8 >>>(d_a, numBlocks);
    
    // block until the device has completed
    cudaThreadSynchronize();
    
    // check if kernel execution generated an error
    checkCUDAError("kernel execution");
    
    // Part 4 of 5: device to host copy
    cudaMemcpy(h_a, d_a, memSize, cudaMemcpyDeviceToHost );
    
    // Check for any CUDA errors
    checkCUDAError("cudaMemcpy");
    
    // Part 5 of 5: verify the data returned to the host is correct
    for (int i = 0; i < numBlocks  ; i++)
    {
        for (int j = 0; j < numThreadsPerBlock  ; j++)
        {
            assert(h_a[i * numThreadsPerBlock + j] == i + j);
        }
    }
    
    // free device memory
    cudaFree(d_a);
    
    // free host memory
    free(h_a);
    
    // If the program makes it this far, then the results are correct and
    // there are no run-time errors.  Good work!
    printf("Correct!\n");
    
}

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }                         
}
