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
void insideBox(T *treeArray_values, int *treeArray_ID, T *box, int tree_size, int number_of_dimensions);

template <class T>
class Cuda_class{
    //(int tree_size_, int number_of_trees_, int number_of_dimensions_, T *d_treeArray_values_, int *d_treeArray_ID_, T *d_box_):tree_size(tree_size_),  number_of_trees(number_of_trees_), number_of_dimensions(number_of_dimensions_), d_treeArray_ID(d_treeArray_ID_), d_box(d_box_){
public:
    
    //void cudaMain(int number_of_trees, int tree_size, T *treeArray_values, int *treeArray_ID, T box[],  int number_of_dimensions);
    void cudaInsideBox(int number_of_trees, int tree_size, int number_of_dimensions, T *treeArray_values, int *treeArray_ID, T box[]);
    void cudaCopyToDevice(int number_of_trees, int tree_size, T *treeArray_values, int *treeArray_ID, T box[],  int number_of_dimensions);
    void cudaCopyToHost(int* treeArray_ID);
private:
    int tree_size;
    int number_of_trees;
    int number_of_dimensions;
    int size_of_forest;// = number_of_trees*tree_size*sizeof(int);
    T *d_treeArray_values;
    int *d_treeArray_ID;
    T *d_box;
};

