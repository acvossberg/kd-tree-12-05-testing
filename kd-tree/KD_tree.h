//
//  KD_tree.h
//  kd-tree-08-05
//
//  Created by Ann-Christine Vossberg on 5/8/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#ifndef __kd_tree_08_05__KD_tree__
#define __kd_tree_08_05__KD_tree__

#include <stdio.h>
using namespace std;

template <typename num_t>
struct Point{
    num_t  x,y;//,z; //dim = 1, 2, 3
};
template <typename num_t>
struct PointCloud
{
    std::vector<Point<num_t>>  pts;
};

template <class T>
class KD_tree{
public:
    KD_tree(vector<Point<T>> &data, vector<int> &dimensions);
    void printTree();
    void KD_tree_recursive(int left, int right, int k);
    
    
private:
    void selectMedian(int dim, int median, int left, int right);
    vector<Point<T>> data;
    vector<int> dim;
    
};


#endif /* defined(__kd_tree_08_05__KD_tree__) */
