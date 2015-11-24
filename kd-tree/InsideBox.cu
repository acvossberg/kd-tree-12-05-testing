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

//more generic approach!!

template <typename T>
__device__
void traverseTree( T *treeArray_values, int *treeArray_ID, T *box, int pos, int startOfTree, int endOfTree, int number_of_dimensions){
    
    //printf("\n first value: startOfTree+pos*blokDim.y + row %d, row %d", (startOfTree+pos)*blockDim.y+row, row);
    //printf("\n box[0] %d, box[1] %d, box[2] %d, box[3] %d, box[4] %d, box[5] %d, thread %d", box[0], box[1], box[2], box[3], box[4], box[5], threadIdx.x);
    bool flag_not_found_invalid = true;
    int lastLevel = ceil(log2(double(endOfTree+1))-1);
    
    
    if(startOfTree + pos - 1 <= endOfTree){
        //calculate which tree level we are on to know which dimension was sorted
        int level = ceil(log2(double(pos+1))-1);
        int level_of_dimension = level%number_of_dimensions;
        //a mod b = a - floor(a / b) * b
        
        //printf("\n level %d, lastlevel %d,levelOfDimension %d, ID %d, pos %d, startOfTree %d, thrad %d", level, lastLevel,level_of_dimension, treeArray_ID[startOfTree+pos-1], pos, startOfTree, threadIdx.x);
        
        
        
        //check wether invalid encountered, continue search:
        if((treeArray_ID[startOfTree+pos-1] != -1)){
            
            //if node has sorted dimension in box, continue both branches:
            if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] >= box[2*level_of_dimension] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] <= box[2*level_of_dimension+1]){
                
                //check wether the node is inside the box:
                for(int i=0; i<number_of_dimensions; i++){
                    if(i == level_of_dimension) continue;
                    if( treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] >= box[2*i] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] <= box[2*i+1] ){
                    //entirely inside box for all dimensions
                    }
                    else{
                        //not totally inside box
                        //printf("\n thread %d is changing ID %d of tree starting at %d, exact position: %d", threadIdx.x, treeArray_ID[startOfTree+pos-1], startOfTree, startOfTree+pos-1);
                        treeArray_ID[startOfTree+pos-1] = -2;
                    }
                }
                
                //continue both branches:
                if( ! (level == lastLevel-1 && flag_not_found_invalid==false)){
                    //left child:
                    pos *= 2;
                    traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
                if( ! (level == lastLevel-1 && flag_not_found_invalid==false)){
                    //right child:
                    pos += 1;
                    traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            }
            //if sorted dimension is larger than box follow branch of smaller child = left child
            else if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level] > box[2*level+1]){
                treeArray_ID[startOfTree+pos-1] = -2;
                if( ! (level == lastLevel-1 && flag_not_found_invalid==false)){
                    //left child:
                    pos *= 2;
                    traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            }
            //if sorted dimension is smaller than box, follow branch of larger child = right child
            else{
                treeArray_ID[startOfTree+pos-1] = -2;
                if( ! (level == lastLevel-1 && flag_not_found_invalid==false)){
                    //right child:
                    pos += 1;
                    traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            }
        }
        //if invalid is encountered, all future lastLevel checks don't have to be done
        //flag is set
        else{
            flag_not_found_invalid = false;
        }
            
            
            
            
            
            /*
            if( treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] >= box[2*i] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] <= box[2*i+1] ){
                
                //inside box
                //printf("\n value: %d , box-min: %d, box-max: %d, ID: %d, thread %d", treeArray_values[startOfTree+number_of_dimensions*(pos-1)+i], box[2*i], box[2*i+1], treeArray_ID[startOfTree+pos-1], threadIdx.x);
            }
            else{
                //printf("\n value: %d , box-min: %d, box-max: %d, ID: %d, thread %d NOOOOOT", treeArray_values[startOfTree+number_of_dimensions*(pos-1)+i], box[2*i], box[2*i+1], treeArray_ID[startOfTree+pos-1], threadIdx.x);
                treeArray_ID[startOfTree+pos-1] = -1;
                
                //left child:
                pos *= 2;
                traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree, number_of_dimensions);
                
                //right child:
                pos += 1;
                traverseTree(treeArray_values, treeArray_ID, box, pos, startOfTree, endOfTree, number_of_dimensions);c
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
        }*/
    }
}
template <typename T>
__device__
void traverseTreeRecursiveIF( T *treeArray_values, int *treeArray_ID, int *treeArray_results, T *box, int pos, int startOfTree, int endOfTree, int number_of_dimensions){
    
    if(startOfTree + pos - 1 <= endOfTree){
        
        //calculate which tree level we are on to know which dimension was sorted
        int level = ceil(log2(double(pos+1))-1);
        int level_of_dimension = level%number_of_dimensions;
        int lastLevel = ceil(log2(double(endOfTree+1))-1);
        //a mod b = a - floor(a / b) * b
        
        //cout << "global position: " << startOfTree+pos -1 << " und ID: " << treeArray_ID[startOfTree+pos-1] << " level: " << level <<  " levelofDimension: " << level_of_dimension << endl;
        
        
        //check wether invalid encountered, continue search:
        if(treeArray_ID[startOfTree+pos-1] != 0){
            
            //if node has sorted dimension in box, continue both branches:
            //cout << box[2*level_of_dimension] <<  " <= " << treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] << " <= " << box[2*level_of_dimension+1] << endl;
            
            if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] >= box[2*level_of_dimension] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] <= box[2*level_of_dimension+1]){
                bool inside = true;
                //cout << "number of nodes traversed " << number_of_nodes_traversed << endl;
                
                //check wether the node is inside the box:
                for(int i=0; i<number_of_dimensions; i++){
                    if(i == level_of_dimension) continue;
                    //cout << treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] << " >= " << box[2*i] <<  " && " << treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] << " <= " << box[2*i+1] << endl;
                    
                    if( treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] >= box[2*i] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] <= box[2*i+1] ){
                        //entirely inside box for all dimensions
                        
                    }
                    else{
                        //not totally inside box
                        //printf("\n thread %d is changing ID %d of tree starting at %d, exact position: %d", threadIdx.x, treeArray_ID[startOfTree+pos-1], startOfTree, startOfTree+pos-1);
                        //cout << "not inside \t ID: " << treeArray_ID[startOfTree+pos-1] <<" bei pos: " << startOfTree+pos-1 << " weg! BEIDE BRANCHES" << endl;
                        inside = false;
                    }
                }
                if(inside){
                    treeArray_results[startOfTree+pos-1] = treeArray_ID[startOfTree+pos-1];
                    //cout << "yes inside \t ID: " << treeArray_ID[startOfTree+pos-1] << " BEIDE BRANCHES" << endl;
                }
                
                //continue both branches:
                if(level != lastLevel){
                    //cout << "ID: " << treeArray_ID[startOfTree+pos-1] <<" bei pos: " << startOfTree+pos-1 << " BEIDE BRANCHES" << endl;
                    //left child:
                    pos *= 2;
                    traverseTreeRecursiveIF(treeArray_values, treeArray_ID,treeArray_results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                    
                    //right child:
                    pos += 1;
                    traverseTreeRecursiveIF(treeArray_values, treeArray_ID,treeArray_results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
                
            }
            //if sorted dimension is larger than box follow branch of smaller child = left child
            else if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] > box[2*level_of_dimension+1]){
                //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " LINKES KIND" << endl;
                //cout << "number of nodes traversed " << number_of_nodes_traversed << endl;
                //left child:
                if(level != lastLevel){
                    pos *= 2;
                    traverseTreeRecursiveIF(treeArray_values, treeArray_ID,treeArray_results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            }
            //if sorted dimension is smaller than box, follow branch of larger child = right child
            else{
                //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " RECHTES KIND" << endl;
                //cout << "number of nodes traversed " << number_of_nodes_traversed << endl;
                if(level != lastLevel){
                    //right child:
                    pos *= 2;
                    pos += 1;
                    traverseTreeRecursiveIF(treeArray_values, treeArray_ID, treeArray_results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            }
        }
    }
}

template <typename T>
__device__
void traverseTreeIterative( T *treeArray_values, int *treeArray_ID, int *treeArray_results, T *box, int *queue, int pos, int startOfTree, int endOfTree, int number_of_dimensions){
    
    int lastLevel = ceil(log2(double(endOfTree+1))-1);
    
    
    queue[startOfTree] = pos;
    int queueFront = 0;
    int queueRear = 1;
    int queueSize = 1;
    int numberOfMightHits = 0;
    
    
    while(queueSize != 0){
        queueSize--;
        pos = queue[queueFront++];
        
        int level = ceil(log2(double(pos+1))-1);
        int level_of_dimension = level%number_of_dimensions;
        
        //cout << "position " << pos-1 << endl;
        
        //if sorted dimension inside box continue with both branches
        if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] >= box[2*level_of_dimension] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] <= box[2*level_of_dimension+1]){
            //put left and right child in queue:
            //cout << "ID: " << treeArray_ID[startOfTree+pos-1] <<" bei pos: " << startOfTree+pos-1 << " BEIDE BRANCHES" << endl;
            
            //nested ifs - bad!
            if(level != lastLevel){
                queue[queueRear++] = pos*2;
                queue[queueRear++] = pos*2+1;
                queueSize+=2;
            }
            
            //and check if node is totally inside box - then right it to results
            //possibility(?) can this be checked after tree traversed? Can an extra thread check this? 2 threads per tree?
            //write pos to array? What about size of array on GPU? Size of tree must be allocated? Is there so much space?
            //per warp needed memory: 3*TreeSize - NO !!! can be put into results :)
            //toCheckIfInside.push_back(pos);
            //these results have to be checked again!!!
            treeArray_results[startOfTree+numberOfMightHits] = pos; //treeArray_ID[startOfTree+pos-1];
            numberOfMightHits++;
            
        }
        //else if sorted dimensions > inside box continue with left child
        else if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] > box[2*level_of_dimension+1]
                && level != lastLevel){
            //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " LINKES KIND" << endl;
            queue[queueRear++] = pos*2;
            queueSize++;
        }
        //else sorted dimension < inside box continue with right child
        else if(level != lastLevel){
            //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " RECHTES KIND" << endl;
            queue[queueRear++] = pos*2+1;
            queueSize++;
        }
    }
    
    //check nodes, that might be inside box:
    for(int j = 0; j<=numberOfMightHits;j++){
        pos = treeArray_results[startOfTree+j];
        treeArray_results[startOfTree+j] = 0;
        
        bool inside = true;
        for(int i=0; i<number_of_dimensions; i++){
            if( treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] >= box[2*i] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] <= box[2*i+1] ){
                //entirely inside box for all dimensions
            }
            else{
                inside = false;
            }
        }
        if(inside){
            treeArray_results[startOfTree+pos-1] = treeArray_ID[startOfTree+pos-1];
        }
    }
}

template <typename T>
__global__
void insideBox(T *treeArray_values, int *treeArray_ID, int *treeArray_results, T *box, int *queue, int tree_size, int number_of_dimensions){
    
    //for each thread has it's own tree starting here
    //TODO: STARTOFTREE falsch.. ist die gesamte position, ohne berücksichtigung der number_of_dimensions. Die müssen berücksichtigt werden!!!
    int startOfTree = threadIdx.x * tree_size;
    int endOfTree = startOfTree + (tree_size - 1);
    /*for(int i = startOfTree; i<endOfTree; i++ ){
        printf("\n (%d, %d, %d) \t ID: %d,\t startOfTree: %d, \t position of point: %d, thread: %d", treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*i+0], treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*i+1], treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*i+2], treeArray_ID[startOfTree+i], startOfTree*number_of_dimensions, startOfTree*number_of_dimensions+number_of_dimensions*i, threadIdx.x);
    }
    printf("\n threadIdx: %d startOfTree %d, endOfTree %d, row %d, col %d, blockDim %d", threadIdx.x, startOfTree, endOfTree, row, col, blockDim.y);
     */
    //traverseTree(treeArray_values, treeArray_ID, box, 1, startOfTree, endOfTree, number_of_dimensions);
    //traverseTreeRecursiveIF(treeArray_values, treeArray_ID, treeArray_results, box, 1, startOfTree, endOfTree, number_of_dimensions);
    traverseTreeIterative(treeArray_values, treeArray_ID, treeArray_results, box, queue, 1, startOfTree, endOfTree, number_of_dimensions);
}

template <typename T>
void Cuda_class<T>::cudaInsideBox(int number_of_trees, int tree_size, int number_of_dimensions, T *treeArray_values, int *treeArray_ID, int *treeArray_results, T box[], int* queue){
    
    //insideBox<T><<<Anzahl benutzte Blöcke, Anzahl Threads>>> = <<<Anzahl benutzte Blöcke, Anzahl Baeume >>>
    //weil ein Thread == ein Baum
    int max_threads_per_block = 1024;
    int number_of_blocks = (number_of_trees + max_threads_per_block -1) / max_threads_per_block ;
    if(number_of_trees > 1024){number_of_trees = 1024;}
    
    insideBox<T><<<number_of_blocks,number_of_trees>>>(d_treeArray_values, d_treeArray_ID, d_treeArray_results, d_box, d_queue, tree_size, number_of_dimensions);
    //YourKernel<<<dimGrid, dimBlock>>>(d_A,d_B); //Kernel invocation
}

template <typename T>
void Cuda_class<T>::cudaCopyToDevice(int number_of_trees_, int tree_size_, T *treeArray_values, int *treeArray_ID, int* treeArray_results, T box[],int *queue, int number_of_dimensions_){
    number_of_trees = number_of_trees_;
    number_of_dimensions = number_of_dimensions_;
    tree_size = tree_size_;

    cudaSetDevice(MYDEVICE);
    std::cout << "number of trees: " << number_of_trees << std::endl;
    std::cout << "tree size: " << tree_size << std::endl;
    std::cout << "number of dimensions: " << number_of_dimensions << std::endl;
    std::cout << "box: " << box[0] << " " << box[1] << " " << box[2] << " " << box[3] << " " << box[4] << " " << box[5] << std::endl;
    //TODO: int ----> num_t
    size_of_forest =  number_of_trees*tree_size*sizeof(int);
    
    //allocate memory
    //TODO: do outside of cudaMain
    cudaMalloc(&d_treeArray_values, size_of_forest*number_of_dimensions);
    cudaMalloc(&d_treeArray_ID, size_of_forest);
    cudaMalloc(&d_treeArray_results, size_of_forest);
    cudaMalloc(&d_box, number_of_dimensions*2*sizeof(T));
    cudaMalloc(&d_queue, size_of_forest);
    
    //send trees to gpu
    cudaMemcpy(d_treeArray_values, treeArray_values, size_of_forest*number_of_dimensions, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_ID, treeArray_ID, size_of_forest, cudaMemcpyHostToDevice);
    cudaMemcpy(d_treeArray_results, treeArray_results, size_of_forest, cudaMemcpyHostToDevice);
    cudaMemcpy(d_box, box, number_of_dimensions*2*sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(d_queue, queue, size_of_forest, cudaMemcpyHostToDevice);
}

template <typename T>
void Cuda_class<T>::cudaCopyToHost(int* treeArray_results){
    
    cudaMemcpy(treeArray_results, d_treeArray_results, size_of_forest, cudaMemcpyDeviceToHost);
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Error: %s\n", cudaGetErrorString(err));
    
    std::cout << "\n Size of forest: " << size_of_forest << std::endl;
    
//    //print out ID's which are in box:
//        for(int i = 0; i< number_of_trees*tree_size; i++){
//            std::cout << "ID: " << treeArray_ID[i]<< std::endl;
//        }
    
    
    //free space
    cudaFree(d_treeArray_values);
    cudaFree(d_treeArray_ID);
    cudaFree(d_treeArray_results);
    cudaFree(d_box);

}

template class Cuda_class<int>;
