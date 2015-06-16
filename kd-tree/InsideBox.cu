//
//  InsideBox.cu
//  
//
//  Created by Ann-Christine Vossberg on 6/3/15.
//
//

#include <stdio.h>
#include <iostream>
//#include "cuPrintf.cu" braucht es?
#include "InsideBox.hpp"


// Here you can set the device ID that was assigned to you
#define MYDEVICE 0

// Simple utility function to check for CUDA runtime errors
void checkCUDAError(const char *msg);


//think about: threadIDx.y!! index s.d. threadIDx.y die treeArray_y bearbeitet -
//for this would have to change treeArray_x & treeArray_y etc --> tree..
//no nested if's no recursive.. -
__global__
void insideBox( int *treeArray_x, int *treeArray_y, int *treeArray_z, int *treeArray_ID, int *box)
{
    //TODO: change to get warpSize
    int warpSize = 32;
    int warpIdx = threadIdx.x / warpSize;
    int i = warpIdx;
    int j = threadIdx.x % 32; //=0 bis 32
    int treeSize = warpSize-1;
    int index = i*treeSize+j;
    
    
    
    //evaluate with another if(which is the next, then you exit from the if.. and you follow the next
    if( ((treeArray_x[index] >= box[0] && treeArray_x[index] <= box[1]) || (box[0] == 0 && box[1] == 0))  && ((treeArray_y[index] >= box[2] && treeArray_y[index] <= box[3]) || (box[2] == 0 && box[3] == 0)) && ((treeArray_z[index] >= box[4] && treeArray_z[index] <= box[5]) || (box[4] == 0 && box[5] == 0))){
        //inside box
        
    }
        else{
            //printf("not inside box");
        treeArray_ID[index] = -1;
    }
}

/*__device__
void traverseTree(int *treeArray_x, int *treeArray_y, int *treeArray_z, int *treeArray_ID, int *box, int pos, int startOfTree, int endOfTree){
    
    //check if inside box:
    if( ((treeArray_x[pos] >= box[0] && treeArray_x[pos] <= box[1]) || (box[0] == 0 && box[1] == 0))  && ((treeArray_y[pos] >= box[2] && treeArray_y[pos] <= box[3]) || (box[2] == 0 && box[3] == 0)) && ((treeArray_z[pos] >= box[4] && treeArray_z[pos] <= box[5]) || (box[4] == 0 && box[5] == 0))){
        //inside box
        printf("I am thread nr. %d", threadIdx.x);
    }
    else{
        printf("not inside box");
        treeArray_ID[pos] = -1;
    }
    
    //Abbruchkriterium:
    //TODO: < oder <= ???
    if(pos <= endOfTree){
        
       //left child:
        pos = pos + pos % startOfTree; //caution: nicht *= 2; weil
       traverseTree(treeArray_x, treeArray_y, treeArray_z, treeArray_ID, box, pos, startOfTree, endOfTree);
        
       //right child:
       pos+=1;
       traverseTree(treeArray_x, treeArray_y, treeArray_z, treeArray_ID, box, pos, startOfTree, endOfTree);
    }
}


//Each thread starts at Node 0 of it's "own" tree. Traverses tree and changes treeArray_ID to -1, if not inside box.
__global__
void Insidebox(int *treeArray_x, int *treeArray_y, int *treeArray_z, int *treeArray_ID, int *box, int tree_size){
    
    //for each thread has it's own tree starting here (+1, because of 2*i, 2*i+1 sonst = 0)
    int startOfTree = threadIdx.x * tree_size + 1; //1, 32, 64, ...
    int endOfTree = startOfTree + tree_size;
    
    traverseTree(treeArray_x, treeArray_y, treeArray_z, treeArray_ID, *box, startOfTree);

}
*/


//TODO: template
void cudaMain(int number_of_trees, int tree_size, int treeArray_x[], int treeArray_y[], int treeArray_z[], int treeArray_ID[], int box[]){
    
    cudaSetDevice(MYDEVICE);
    std::cout << "number of trees: " << number_of_trees << std::endl;
    std::cout << "tree size: " << tree_size << std::endl;
    
    //TODO: int ----> num_t
    int size_of_forest = number_of_trees*tree_size*sizeof(int);
    int *d_treeArray_x;
    int *d_treeArray_y;
    int *d_treeArray_z;
    int *d_treeArray_ID;
    int *d_box;
    
    
    //allocate memory
    cudaMalloc(&d_treeArray_x, size_of_forest);
    cudaMalloc(&d_treeArray_y, size_of_forest);
    cudaMalloc(&d_treeArray_z, size_of_forest);
    cudaMalloc(&d_treeArray_ID, size_of_forest);
    //TODO: generic
    cudaMalloc(&d_box, 6*sizeof(int));
    
    //send trees to gpu
    cudaMemcpy(d_treeArray_x, treeArray_x, size_of_forest, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_y, treeArray_y, size_of_forest, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_z, treeArray_z, size_of_forest, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_ID, treeArray_ID, size_of_forest, cudaMemcpyHostToDevice);
    //TODO: generic
    cudaMemcpy(d_box, box, 6*sizeof(int), cudaMemcpyHostToDevice);
    
    
    //search forest for points inside box_dimensions - returns all treeArray_ID's which are inside box - rest are filled with -1
    insideBox<<<1,1024>>>(d_treeArray_x, d_treeArray_y, d_treeArray_z, d_treeArray_ID, d_box);
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Error: %s\n", cudaGetErrorString(err));
    
    cudaMemcpy(treeArray_ID, d_treeArray_ID, size_of_forest, cudaMemcpyDeviceToHost);
    
    
    std::cout << "Size of forest: " << size_of_forest << std::endl;
    //print out ID's which are in box:
    for(int i = 0; i< number_of_trees*tree_size; i++){
        std::cout << "ID: " << treeArray_ID[i]<< std::endl;
    }
    
    
    //free space
    cudaFree(d_treeArray_x);
    cudaFree(d_treeArray_y);
    cudaFree(d_treeArray_z);
    cudaFree(d_treeArray_ID);
    cudaFree(d_box);

}
