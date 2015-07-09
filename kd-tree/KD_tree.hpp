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
#include "sorter.hpp"
#include <vector>

using namespace std;

template < class num>
class SimpleKDtree;

template <class T>
class KD_tree{
public:
    KD_tree(vector<Point<T>> &data, vector<vector<T>> &vectorData_, vector<int> &dataID_,int begin_, int end_, vector<int>& dimensions, T **transformable_trees_,int *treesArray_ID, int& offset_);
    void printTree();
    void KD_tree_recursive(int left, int right, int k, int pos);
    bool testTree( SimpleKDtree<T> *simpleTree);
    vector<Point<T>> get_tree_as_vector();
    T get_value(int d, Point<T> val );

    
private:
    void printData();
    void selectMedian(int d, int median, int left, int right, int pos);
    void original_order_median(int median, int d, int left, int right);
    void original_order_median_new(int median, int d, int left, int right);
    void original_order_median2D(int median_position, int d, int left, int right);
    
    T qselect(int k, int li, int hi, int dimension_offset);
    T qselect2D(int k, int li, int hi, int dimension_offset);
    
    
    vector<Point<T>>& data;
    vector<vector<T>>& vectorData;
    vector<int> dataID;
    int begin;
    int end;
    const vector<int>& dim;
    vector<Point<T>> result;
    T** transformable_trees;
    int* treesArray_ID;
    int& offset;

    
};


#endif /* defined(__kd_tree_08_05__KD_tree__) */
