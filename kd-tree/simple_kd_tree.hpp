//
//  simple_kd_tree.h
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 5/12/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#ifndef __kd_tree_12_05__simple_kd_tree__
#define __kd_tree_12_05__simple_kd_tree__

#include <stdio.h>
#include <vector>
#include "sorter.hpp"
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
    SimpleKDtree(const vector<int> &dimensions);
    void make_SimpleKDtree(vector<Point<T>> cloud, int left, int right, int d);
    KDnode<T> *newSimpleKDtreeNode(Point<T> data);
    bool sameTree(vector<Point<T>> Kdtree, int i);
    vector<Point<T>> get_as_vector();
    
    ~SimpleKDtree();
    
private:
    void addHelper(KDnode<T> *root, Point<T> val);
    void add(Point<T> val);
    bool sameTreeHelper(KDnode<T> *n, int i, vector<Point<T>> Kdtree);
    void make_SimpleKDtreeHelper(vector<Point<T>> cloud, KDnode<T> *root, int left, int right, int k);
    bool compare_tree_node(KDnode<T> root, Point<T> Kdtree);
    //T get_value(int d, KDnode<T> *val );
    T get_value(int d, Point<T> val );
    
    
    int ID_val = 4;
    vector<int> dim;
    KDnode<T> *root;
    int height;
    int number_nodes;
    Point<T> zero_point;
    vector<Point<T>> res_vector;
};

#endif /* defined(__kd_tree_12_05__simple_kd_tree__) */


