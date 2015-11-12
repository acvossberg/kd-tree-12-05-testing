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
void traverseTree( T *treeArray_values, int *treeArray_ID, T *box, int pos, int startOfTree, int endOfTree, int number_of_dimensions){
    
    //int row = blockIdx.y * blockDim.y + threadIdx.y;
    //printf("\n first value: startOfTree+pos*blokDim.y + row %d, row %d", (startOfTree+pos)*blockDim.y+row, row);
    //printf("\n box[0] %d, box[1] %d, box[2] %d, box[3] %d, box[4] %d, box[5] %d, thread %d", box[0], box[1], box[2], box[3], box[4], box[5], threadIdx.x);
    
    //TODO: erstes/ erst hits in tree wird nicht processed
    if(startOfTree + pos - 1 <= endOfTree){
        
        for( int i=0; i<number_of_dimensions; i++){
            if( treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] >= box[2*i] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] <= box[2*i+1] ){
                //inside box
                //printf("\n value: %d , box-min: %d, box-max: %d, ID: %d, thread %d", treeArray_values[startOfTree+number_of_dimensions*(pos-1)+i], box[2*i], box[2*i+1], treeArray_ID[startOfTree+pos-1], threadIdx.x);
            }
            else{
                //printf("\n value: %d , box-min: %d, box-max: %d, ID: %d, thread %d NOOOOOT", treeArray_values[startOfTree+number_of_dimensions*(pos-1)+i], box[2*i], box[2*i+1], treeArray_ID[startOfTree+pos-1], threadIdx.x);
                treeArray_ID[startOfTree+pos-1] = -1;
            }
            
        }
        
        //Abbruchkriterium:
        //TODO: < oder <= ???
        //TODO: if hier weg
        if(startOfTree+pos < endOfTree){
            
            //left child:
            pos *= 2;
            traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree, number_of_dimensions);
            
            //right child:
            pos += 1;
            traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree, number_of_dimensions);
        }
    }
}


template <typename T>
__global__
void insideBox(T *treeArray_values, int *treeArray_ID, T *box, int tree_size, int number_of_dimensions){
    
    //for each thread has it's own tree starting here
    //TODO: STARTOFTREE falsch.. ist die gesamte position, ohne berücksichtigung der number_of_dimensions. Die müssen berücksichtigt werden!!!
    int startOfTree = threadIdx.x * tree_size;
    int endOfTree = startOfTree + (tree_size - 1);
    for(int i = startOfTree; i<endOfTree; i++ ){
        //printf("\n (%d, %d, %d) \t ID: %d,\t startOfTree: %d, \t position of point: %d, thread: %d", treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*i+0], treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*i+1], treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*i+2], treeArray_ID[startOfTree+i], startOfTree*number_of_dimensions, startOfTree*number_of_dimensions+number_of_dimensions*i, threadIdx.x);
    }
    //printf("\n threadIdx: %d startOfTree %d, endOfTree %d, row %d, col %d, blockDim %d", threadIdx.x, startOfTree, endOfTree, row, col, blockDim.y);
    traverseTree(treeArray_values, treeArray_ID, box, 1, startOfTree, endOfTree, number_of_dimensions);
}


template <typename T>
//void Cuda_class<T>::cudaMain(int number_of_trees, int tree_size, T treeArray_x[], T treeArray_y[], T treeArray_z[], int treeArray_ID[], T box[]){
void Cuda_class<T>::cudaMain(int number_of_trees, int tree_size, T *treeArray_values, int *treeArray_ID, T box[],  int number_of_dimensions){
    
    cudaSetDevice(MYDEVICE);
    std::cout << "number of trees: " << number_of_trees << std::endl;
    std::cout << "tree size: " << tree_size << std::endl;
    std::cout << "number of dimensions: " << number_of_dimensions << std::endl;
    std::cout << "box: " << box[0] << " " << box[1] << " " << box[2] << " " << box[3] << " " << box[4] << " " << box[5] << std::endl;
    //TODO: int ----> num_t
    int size_of_forest = number_of_trees*tree_size*sizeof(int);
    T *d_treeArray_values;
    int *d_treeArray_ID;
    T *d_box;
    
    
    //allocate memory
    //TODO: do outside of cudaMain
    cudaMalloc(&d_treeArray_values, size_of_forest*number_of_dimensions);
    cudaMalloc(&d_treeArray_ID, size_of_forest);
    cudaMalloc(&d_box, number_of_dimensions*2*sizeof(T));
    
    //send trees to gpu
    cudaMemcpy(d_treeArray_values, treeArray_values, size_of_forest*number_of_dimensions, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_ID, treeArray_ID, size_of_forest, cudaMemcpyHostToDevice);
    cudaMemcpy(d_box, box, number_of_dimensions*2*sizeof(T), cudaMemcpyHostToDevice);
    
    
    //search forest for points inside box_dimensions - returns all treeArray_ID's which are inside box - rest are filled with -1
    //TODO: do not change treeArray_ID's - make separate array.
    
    //insideBox<T><<<Anzahl benutzte Blöcke, Anzahl Threads>>> = <<<Anzahl benutzte Blöcke, Anzahl Baeume >>>
    //weil ein Thread == ein Baum
    insideBox<T><<<1,number_of_trees>>>(d_treeArray_values, d_treeArray_ID, d_box, tree_size, number_of_dimensions);
    //YourKernel<<<dimGrid, dimBlock>>>(d_A,d_B); //Kernel invocation
    
    
    cudaMemcpy(treeArray_ID, d_treeArray_ID, size_of_forest, cudaMemcpyDeviceToHost);
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Error: %s\n", cudaGetErrorString(err));
    
    std::cout << "\n Size of forest: " << size_of_forest << std::endl;

    //print out ID's which are in box:
//    for(int i = 0; i< number_of_trees*tree_size; i++){
//        std::cout << "ID: " << treeArray_ID[i]<< std::endl;
//    }
    
    
    //free space
    cudaFree(d_treeArray_values);
    cudaFree(d_treeArray_ID);
    cudaFree(d_box);

}

template class Cuda_class<int>;
