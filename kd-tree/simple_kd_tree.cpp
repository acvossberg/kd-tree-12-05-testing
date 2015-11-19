//
//  simple_kd_tree.cpp
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 5/12/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#include "simple_kd_tree.hpp"
#include "KD_tree.hpp"
#include <vector>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>


using namespace std;
template < class T>
SimpleKDtree<T>::~SimpleKDtree() { // destructor implementation
    res_vector.clear();
    if(root){
        delete[] root;}
    dim.clear();
}

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
SimpleKDtree<T>::SimpleKDtree(vector<int> const &dimensions){
    dim = dimensions;
    zero_point.x = zero_point.y = zero_point.z = 0;

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
void SimpleKDtree<T>::make_SimpleKDtreeHelper(vector<Point<T>> cloud, KDnode<T> *root, int left, int right, int k){
    int mid = cloud.size()/2;
    height--;
    if(cloud.size() > 1){

        k = k%dim.size();
        //s.d. stable_sort die richtige Reihenfolge hat, jedes mal nach indices sortieren:
        sort(cloud.begin(), cloud.end(), sortern<T>(ID_val));
        stable_sort(cloud.begin(), cloud.end(), sortern<T>(dim[k]));
        
        if(root){
            root->values = cloud[mid];
        }
        else{
            root = new KDnode<T>(cloud[mid]);
        }
        
        vector<Point<T>> cloud_left(cloud.begin(), cloud.begin() + mid);
        vector<Point<T>> cloud_right(cloud.begin()+mid+1, cloud.end());
        root->left = new KDnode<T>(zero_point);
        root->right = new KDnode<T>(zero_point);
        
        k++;
        //left:
        make_SimpleKDtreeHelper(cloud_left,root->left, left, mid-1, k);
        height++;
        //right:
        make_SimpleKDtreeHelper(cloud_right,root->right, mid+1, right, k);
        height++;
    }
    if(cloud.size() == 1){
        root->values = cloud[0];
        if(height == 1){
            Point<T> zero_point;
            zero_point.x = zero_point.y = zero_point.z = 0;
            //cout << "put extra children here: " << cloud[0].x << " " << cloud[0].y << endl;
            root->left = new KDnode<T>(zero_point);
            root->right = new KDnode<T>(zero_point);
        }
    }

}

template <class T>
void SimpleKDtree<T>::make_SimpleKDtree(vector<Point<T>> cloud, int left, int right, int k){
    height = floor(log2(cloud.size()))+1;
    number_nodes = pow(2,floor(log2(cloud.size())) +1) - 1;
    Point<T> placeholder;
    placeholder.x = 0; //TODO: some better marker/placeholder needed!
    placeholder.y = 0;
    placeholder.z = 0;
    res_vector.resize(number_nodes,placeholder);
    
    if(root){
        root = new KDnode<T>();
        this->make_SimpleKDtreeHelper(cloud, root, 0, cloud.size()-1, 0);
    }
    else{
        root = new KDnode<T>();
        this->make_SimpleKDtreeHelper(cloud, root, 0, cloud.size()-1, 0);
    }
}
template <class T>
bool SimpleKDtree<T>::compare_tree_node(KDnode<T> root, Point<T> Kdtree){
    
    //std::cout << "simpletree: (" << get_value(1, root.values) << ", " << get_value(2, root.values) << ", " << get_value(3, root.values) << ") ID: "<< get_value(4, root.values) << "\t realtree: (" << get_value(1, Kdtree) << ", " << get_value(2, Kdtree) << ", " << get_value(3, Kdtree) << ") ID:" << get_value(4, Kdtree) << std::endl;
    
    for(int i = 0; i<dim.size()-1; ++i){
        int d = dim[i];
        if(get_value(d, root.values) != get_value(d, Kdtree)){return false;}
    }
    return true;
}

template <class T>
bool SimpleKDtree<T>::sameTreeHelper(KDnode<T> *root, int i, vector<Point<T>> Kdtree){
    
    if(i >= Kdtree.size()) return true;
    if (!root){
        return false;
     }
    
    if(compare_tree_node(root->values, Kdtree[i-1])){
        return sameTreeHelper(root->left, 2*i, Kdtree) && sameTreeHelper(root->right, 2*i+1, Kdtree);
    }
    
    else if(root->values.x == -1 && root->values.y == -1  && Kdtree[i-1].x == 0 && Kdtree[i-1].y == 0 ){
        //reached end-node
        return true;
    }
    else{
        cout << "hups! not the same " << endl;
        return false;
    }
    
    return true;

}


template <class T>
bool SimpleKDtree<T>::sameTree(vector<Point<T>> Kdtree, int i){
    return sameTreeHelper(this->root, i, Kdtree);
}

template < class T>
T SimpleKDtree<T>::get_value(int d, Point<T> val ){
    switch (d){
        case 1:
            return val.x;
        case 2:
            return val.y;
        case 3:
            return val.z;
        case 4:
            return val.ID;
        default:
            cout << "error no correct comparison-value" <<endl;
            return val.x;
            
    }
}
