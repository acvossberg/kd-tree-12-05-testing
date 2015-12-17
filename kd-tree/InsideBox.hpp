//
//  hi.h
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 6/10/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#ifndef __kd_tree_12_05__hi__
#define __kd_tree_12_05__hi__

#include <stdio.h>
#include <iostream>
#define MYDEVICE 0


#endif /* defined(__kd_tree_12_05__hi__) */

template <typename T>
__global__
void insideBox(T *treeArray_x, T *treeArray_y, T *treeArray_z, int *treeArray_ID, T *box, int tree_size);

template <typename T>
__global__
void insideBox(T *treeArray_values, int *treeArray_ID, int *treeArray_results, T *box, int tree_size, int number_of_dimensions);
//void insideBox(T *treeArray_values, int *treeArray_ID, int *treeArray_results, T *box, int *queue, int tree_size, int number_of_dimensions);
template <class T>
class Cuda_class{
public:
    
    //void cudaMain(int number_of_trees, int tree_size, T *treeArray_values, int *treeArray_ID, T box[],  int number_of_dimensions);
    
    void cudaInsideBox(int number_of_trees, int tree_size, int number_of_dimensions, T *treeArray_values, int *treeArray_ID,int *treeArray_results, T box[], int numberOfHits);
    //void cudaInsideBox(int number_of_trees, int tree_size, int number_of_dimensions, T *treeArray_values, int *treeArray_ID,int *treeArray_results, T box[], int *queue, int numberOfHits);
    void cudaCopyToDevice(int number_of_trees, int tree_size, T *treeArray_values, int *treeArray_ID, int *treeArray_results, T box[], int number_of_dimensions, int numberOfHits);
    //void cudaCopyToDevice(int number_of_trees, int tree_size, T *treeArray_values, int *treeArray_ID, int *treeArray_results, T box[], int *queue, int number_of_dimensions, int numberOfHits);
    void cudaCopyToHost(int* treeArray_results, int numberOfHits);
    
private:
    int tree_size;
    int number_of_trees;
    int number_of_dimensions;
    int size_of_forest;
    T *d_treeArray_values;
    int *d_treeArray_ID;
    int *d_treeArray_results;
    T *d_box;
    int *d_queue;
};