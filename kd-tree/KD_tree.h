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
//#include "simple_kd_tree.h"
using namespace std;
template <typename num_t>
struct Point{
    num_t  x,y,z,ID; //dim = 1, 2, 3
};

template < class num>
class SimpleKDtree;

template <class T>
class KD_tree{
public:
    KD_tree(vector<Point<T>> &data, vector<int> dimensions);
    void printTree();
    void KD_tree_recursive(int left, int right, int k, int pos);
    bool testTree( SimpleKDtree<T> *simpleTree);
    vector<Point<T>> get_tree_as_vector();
    
private:
    void printData();
    void selectMedian(int d, int median, int left, int right, int pos);
    void original_order_median(int median, int d, int left, int right);
    T get_value(int d, Point<T> val );
    vector<Point<T>> data;
    const vector<int> dim;
    vector<Point<T>> result;

    
};


#endif /* defined(__kd_tree_08_05__KD_tree__) */
