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
void insideBox(T **treeArray_values, int *treeArray_ID, T *box, int tree_size, int number_of_dimensions);

template <class T>
class Cuda_class{
public:
    
    //void cudaMain(int number_of_trees, int tree_size, T treeArray_x[], T treeArray_y[], T treeArray_z[], int treeArray_ID[], T box[]);
    void cudaMain(int number_of_trees, int tree_size, T **treeArray, int *treeArray_ID, T box[], int number_of_dimensions);

};

