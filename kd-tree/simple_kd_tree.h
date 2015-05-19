//
//  simple_kd_tree.h
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 5/12/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#ifndef __kd_tree_12_05__simple_kd_tree__
#define __kd_tree_12_05__simple_kd_tree__

#include <stdio.h>¨
#include <vector>
#include "KD_tree.h"

using namespace std;

template <class T>
struct KDnode{
    Point<T> values;
    KDnode<T> *left, *right;
    
   KDnode(Point<T> val) {
        this->values = val;
    }
    KDnode(){
        
    }
    
    KDnode(Point<T> val, KDnode<T> left, KDnode<T> right) {
        this->values = val;
        this->left = left;
        this->right = right;
    }
};

template < class T>
class SimpleKDtree{
    

public:
    SimpleKDtree(vector<int> &dimensions);
    void make_SimpleKDtree(vector<Point<T>> cloud, int left, int right, int dim);
    KDnode<T> *newSimpleKDtreeNode(Point<T> data);
    
private:
    vector<int> dim;
    void addHelper(KDnode<T> *root, Point<T> val);
    void add(Point<T> val);
    
    KDnode<T> *root;
    
    
};

#endif /* defined(__kd_tree_12_05__simple_kd_tree__) */


