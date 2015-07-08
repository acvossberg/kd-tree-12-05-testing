//
//  InsideBox.cu
//  
//
//  Created by Ann-Christine Vossberg on 6/3/15.
//
//

#include <stdio.h>
#include <iostream>
#include "InsideBox.hpp"


// Here you can set the device ID that was assigned to you
#define MYDEVICE 0

// Simple utility function to check for CUDA runtime errors
void checkCUDAError(const char *msg);


//think about: threadIDx.y!! index s.d. threadIDx.y die treeArray_y bearbeitet -
//for this would have to change treeArray_x & treeArray_y etc --> tree..
//no nested if's no recursive.. -
template <typename T>
__global__
void insideBox_test( T *treeArray_x, T *treeArray_y, T *treeArray_z, int *treeArray_ID, T *box)
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

template <typename T>
__device__
void traverseTree( T *treeArray_x, T *treeArray_y, T *treeArray_z, int *treeArray_ID, T *box, int pos, int startOfTree, int endOfTree){
    
    //printf("\n threadIdx: %d startOfTree %d, endOfTree %d", threadIdx.x, startOfTree, endOfTree);
    
    if(startOfTree + pos -1 <= endOfTree){
    
        ///CHECK HERE!!!! STARTOFTREE+POS == INDEX  AND POS = i , i->2*i etc.
        if( ((treeArray_x[startOfTree+pos] >= box[0] && treeArray_x[startOfTree+pos] <= box[1]) || (box[0] == 0 && box[1] == 0))  && ((treeArray_y[startOfTree+pos] >= box[2] && treeArray_y[startOfTree+pos] <= box[3]) || (box[2] == 0 && box[3] == 0)) && ((treeArray_z[startOfTree+pos] >= box[4] && treeArray_z[startOfTree+pos] <= box[5]) || (box[4] == 0 && box[5] == 0))){
            //inside box
        }
        else{
            //printf("\n not inside box at position %d with thread nr: %d ", startOfTree+pos, threadIdx.x);
            treeArray_ID[startOfTree+pos] = -1;
        }
    
        //Abbruchkriterium:
        //TODO: < oder <= ???
        //if(startOfTree+pos < endOfTree){
        
            //left child:
            pos *= 2;
            traverseTree(treeArray_x, treeArray_y, treeArray_z, treeArray_ID, box, pos, startOfTree, endOfTree);
        
            //right child:
            pos += 1;
            traverseTree(treeArray_x, treeArray_y, treeArray_z, treeArray_ID, box, pos, startOfTree, endOfTree);
        //}
    }
}


//Each thread starts at Node 0 of it's "own" tree. Traverses tree and changes treeArray_ID to -1, if not inside box.
template <typename T>
__global__
void insideBox(T *treeArray_x, T *treeArray_y, T *treeArray_z, int *treeArray_ID, T *box, int tree_size){
    
    //for each thread has it's own tree starting here
    int startOfTree = threadIdx.x * tree_size ;
    int endOfTree = startOfTree + tree_size - 1;
    traverseTree(treeArray_x, treeArray_y, treeArray_z, treeArray_ID, box, 1, startOfTree, endOfTree);
}



//more generic approach!!

template <typename T>
__device__
void traverseTree( T *treeArray_values, int *treeArray_ID, T *box, int pos, int startOfTree, int endOfTree){
    
    //printf("\n threadIdx: %d startOfTree %d, endOfTree %d", threadIdx.x, startOfTree, endOfTree);
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    if(startOfTree + pos -1 <= endOfTree){
        ///CHECK HERE!!!! STARTOFTREE+POS == INDEX  AND POS = i , i->2*i etc.
       if( ((treeArray_values[row][startOfTree+pos] >= box[row] && treeArray_values[row][startOfTree+pos] <= box[row]) || (box[row] == 0 && box[row] == 0))  && ((treeArray_values[row][startOfTree+pos] >= box[row] && treeArray_values[row][startOfTree+pos] <= box[row]) || (box[row] == 0 && box[row] == 0)) && ((treeArray_values[row][startOfTree+pos] >= box[row] && treeArray_values[row][startOfTree+pos] <= box[row]) || (box[row] == 0 && box[row] == 0))){
            //inside box
        }
        else{
            //printf("\n not inside box at position %d with thread nr: %d ", startOfTree+pos, threadIdx.x);
            treeArray_ID[startOfTree+pos] = -1;
        }
        
        //Abbruchkriterium:
        //TODO: < oder <= ???
        //TODO: if hier weg
        if(startOfTree+pos < endOfTree){
            
            //left child:
            pos *= 2;
            traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree);
            
            //right child:
            pos += 1;
            traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree);
        }
    }
}

template <typename T>
__global__
void insideBox(T *treeArray_values, int *treeArray_ID, T *box, int tree_size, int number_of_dimensions){
    
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    //for each thread has it's own tree starting here
    int startOfTree = threadIdx.x * tree_size ;
    int endOfTree = startOfTree + tree_size - 1;
    printf("\n threadIdx: %d startOfTree %d, endOfTree %d, row %d, col %d", threadIdx.x, startOfTree, endOfTree, row, col);
    //traverseTree(treeArray_values, treeArray_ID, box, 1, startOfTree, endOfTree);
}


template <typename T>
//void Cuda_class<T>::cudaMain(int number_of_trees, int tree_size, T treeArray_x[], T treeArray_y[], T treeArray_z[], int treeArray_ID[], T box[]){
void Cuda_class<T>::cudaMain(int number_of_trees, int tree_size, T *treeArray_values, int *treeArray_ID, T box[],  int number_of_dimensions){
    
    cudaSetDevice(MYDEVICE);
    std::cout << "number of trees: " << number_of_trees << std::endl;
    std::cout << "tree size: " << tree_size << std::endl;
    std::cout << "number of dimensions: " << number_of_dimensions << std::endl;
    
    //TODO: int ----> num_t
    int size_of_forest = number_of_trees*tree_size*sizeof(int);
    T *d_treeArray_values;
    //T *d_treeArray_x;
    //T *d_treeArray_y;
    //T *d_treeArray_z;
    int *d_treeArray_ID;
    T *d_box;
    
    
    //allocate memory
    /*cudaMalloc(&d_treeArray_x, size_of_forest);
    cudaMalloc(&d_treeArray_y, size_of_forest);
    cudaMalloc(&d_treeArray_z, size_of_forest);*/
    cudaMalloc(&d_treeArray_values, size_of_forest);
    cudaMalloc(&d_treeArray_ID, size_of_forest);
    //TODO: generic
    cudaMalloc(&d_box, number_of_dimensions*2*sizeof(T));
    
    //send trees to gpu
    /*cudaMemcpy(d_treeArray_x, treeArray_x, size_of_forest, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_y, treeArray_y, size_of_forest, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_z, treeArray_z, size_of_forest, cudaMemcpyHostToDevice);*/
    cudaMemcpy(d_treeArray_values, treeArray_values, size_of_forest*number_of_dimensions, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_ID, treeArray_ID, size_of_forest, cudaMemcpyHostToDevice);
    //TODO: generic
    cudaMemcpy(d_box, box, number_of_dimensions*2*sizeof(T), cudaMemcpyHostToDevice);
    
    
    //search forest for points inside box_dimensions - returns all treeArray_ID's which are inside box - rest are filled with -1
    //TODO: do not change treeArray_ID's - make separate array.
    dim3 dimBlock(number_of_trees, number_of_dimensions);
    insideBox<T><<<1,dimBlock>>>(d_treeArray_values, d_treeArray_ID, d_box, tree_size, number_of_dimensions);
    //YourKernel<<<dimGrid, dimBlock>>>(d_A,d_B); //Kernel invocation
    
    /*
    //test wether insideBox works
    int *d_treeArray_ID_copy;
    int test_ID[number_of_trees*tree_size];
    //int test_treeArray_ID = std::copy(treeArray_ID);
    cudaMalloc(&d_treeArray_ID_copy, size_of_forest);
    cudaMemcpy(d_treeArray_ID_copy, treeArray_ID, size_of_forest, cudaMemcpyHostToDevice);
    //insideBox_test<<<1,1024>>>(d_treeArray_x, d_treeArray_y, d_treeArray_z, d_treeArray_ID_copy, d_box);
    cudaMemcpy(test_ID, d_treeArray_ID_copy, size_of_forest, cudaMemcpyDeviceToHost);
    //finish test
     */
    
    cudaMemcpy(treeArray_ID, d_treeArray_ID, size_of_forest, cudaMemcpyDeviceToHost);
    
    /*
    bool correctID=true;
    for(int i = 0; i<number_of_trees*tree_size; i++){
        correctID = correctID && (treeArray_ID[i] == test_ID[i]);
    }
    printf("\n All ID's found in box are %d", correctID );
     */
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Error: %s\n", cudaGetErrorString(err));
    
    std::cout << "\n Size of forest: " << size_of_forest << std::endl;
//
//    //print out ID's which are in box:
//    for(int i = 0; i< number_of_trees*tree_size; i++){
//        std::cout << "ID: " << treeArray_ID[i]<< std::endl;
//    }
    
    
    //free space
    /*cudaFree(d_treeArray_x);
    cudaFree(d_treeArray_y);
    cudaFree(d_treeArray_z);*/
    cudaFree(d_treeArray_values);
    cudaFree(d_treeArray_ID);
    cudaFree(d_box);

}
template class Cuda_class<int>;
