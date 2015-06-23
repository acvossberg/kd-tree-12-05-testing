//
//  cuda_class.cu
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 6/18/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#include <stdio.h>
#include <iostream>

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER  __device__ __host__
#else
#define CUDA_CALLABLE_MEMBER
#endif



template <typename T>
extern void cudaMain(int number_of_trees, int tree_size, T treeArray_x[], T treeArray_y[], T treeArray_z[], T treeArray_ID[], T box[]);
    



