//
//  simple_kd_tree.cpp
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 5/12/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#include "simple_kd_tree.h"
#include "KD_tree.h"
#include <vector>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>


using namespace std;

template < class T>
struct KDnode<T> *SimpleKDtree<T>::newSimpleKDtreeNode(Point<T> data)
{
    //FIXME: BAD_ACCESS
    KDnode<T> *node;
    node->values = data ;
    node->left = NULL;
    node->right = NULL;
    
    return node;
}

template < class T>
SimpleKDtree<T>::SimpleKDtree(vector<int> &dimensions){
    dim = dimensions;
}
template < class T>
void SimpleKDtree<T>::addHelper(KDnode<T> *root, Point<T> val) {
    if (root->value > val) {
        if (!root->left) {
            root->left = new KDnode<T>(val);
        } else {
            addHelper(root->left, val);
        }
    } else {
        if (!root->right) {
            root->right = new KDnode<T>(val);
        } else {
            addHelper(root->right, val);
        }
    }
}
template <class T>
void SimpleKDtree<T>::add(Point<T> val) {
    if (root) {
        this->addHelper(root, val);
    } else {
        root = new KDnode<T>(val);
    }
}
template <class T>
void SimpleKDtree<T>::make_SimpleKDtree(vector<Point<T>> cloud, int left, int right, int k){
    int mid = cloud.size()/2;
    
    if(cloud.size() > 1){
        
        //check ob korrekt gerundet wird
        k = k%dim.size();
        
        //printTree();
        sort(cloud.begin(), cloud.end(), sorter(dim[k]));
        
        
        //this->values = cloud[mid];
        if(root){
            root->values = cloud[mid];
        }
        else{
            root = new KDnode<T>(cloud[mid]);
        }
        
        vector<Point<T>> cloud_left(cloud.begin(), cloud.begin() + mid-1);
        vector<Point<T>> cloud_right(cloud.begin()+mid, cloud.end());
        k++;
        root->left = new KDnode<T>;
        root->right = new KDnode<T>;
        //left:
        make_SimpleKDtree(cloud_left, left, mid-1, k);
        //right:
        make_SimpleKDtree(cloud_right, mid+1, right, k);
    }
    if(cloud.size() == 1) root->values = cloud[0];
}


