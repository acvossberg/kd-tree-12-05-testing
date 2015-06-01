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
template < class T>
void SimpleKDtree<T>::make_SimpleKDtreeHelper(vector<Point<T>> cloud, KDnode<T> *root, int left, int right, int k){
    int mid = cloud.size()/2;
    number_nodes--;
    height--;
    if(cloud.size() > 1){
        //check ob korrekt gerundet wird
        k = k%dim.size();
        //s.d. stable_sort die richtige Reihenfolge hat, jedes mal nach indices sortieren:
        sort(cloud.begin(), cloud.end(), sorter(ID_val));
        stable_sort(cloud.begin(), cloud.end(), sorter(dim[k]));
        
        /*cout<< "CLOUD STABLE SORT: " <<endl;
        for(int i = 0; i< cloud.size(); i++){
            cout <<"("<< cloud[i].x << ", " << cloud[i].y << ", " << cloud[i].z << ")" << endl;
        }
        cout << "\n \n" << endl;*/
        
        //this->values = cloud[mid];
        if(root){
            root->values = cloud[mid];
        }
        else{
            root = new KDnode<T>(cloud[mid]);
        }
        
        vector<Point<T>> cloud_left(cloud.begin(), cloud.begin() + mid);
        vector<Point<T>> cloud_right(cloud.begin()+mid+1, cloud.end());
        k++;
        root->left = new KDnode<T>(zero_point);
        root->right = new KDnode<T>(zero_point);
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
            cout << "put extra children here: " << cloud[0].x << " " << cloud[0].y << endl;
            root->left = new KDnode<T>(zero_point);
            root->right = new KDnode<T>(zero_point);
        }
        if(number_nodes <= 3){
            
            cout << "put extra children here: " << cloud[0].x << " " << cloud[0].y << endl;
            root->left = new KDnode<T>(zero_point);
            root->right = new KDnode<T>(zero_point);
        }
    }
    //if(height > 0 && cloud.size() == 0)

}

template <class T>
void SimpleKDtree<T>::make_SimpleKDtree(vector<Point<T>> cloud, int left, int right, int k){
    height = log2(cloud.size())+1;
    number_nodes = pow(2,floor(log2(cloud.size())) +1) - 1;
    if(root){
        this->make_SimpleKDtreeHelper(cloud, root, 0, cloud.size()-1, 0);
    }
    else{
        root = new KDnode<T>();
        this->make_SimpleKDtreeHelper(cloud, root, 0, cloud.size()-1, 0);

    }
    
}


template <class T>
bool SimpleKDtree<T>::sameTreeHelper(KDnode<T> *root, int i, vector<Point<T>> Kdtree){
    
    if(i >= Kdtree.size()) return true;
    //cout << root->values.x << " " << root->values.y << endl;
    if (!root){
     cout<< "could not access the child-node at " << i-1 << " tree: "  << Kdtree[i-1].x << ", " << Kdtree[i-1].y << ", " << Kdtree[i-1].z <<  ")" << endl;
     return false;
     }
    //TODO: change to also compare .z
    if(root->values.x == Kdtree[i-1].x && root->values.y == Kdtree[i-1].y){
        //cout << "pair " << i-1 << " is the same! (" << Kdtree[i-1].x << ", " << Kdtree[i-1].y << ")" << endl;
        return sameTreeHelper(root->left, 2*i, Kdtree) && sameTreeHelper(root->right, 2*i+1, Kdtree);
    }
    else if(root->values.x == 0 && root->values.y == 0 && Kdtree[i-1].x == -1 && Kdtree[i-1].y == -1){
        //end-node
        return true;
    }
    else{
        cout << "at " << i-1 << " NOT same! (" << Kdtree[i-1].x << ", " << Kdtree[i-1].y << " " << Kdtree[i-1].z << ")" << endl;
        cout << "vs " << root->values.x << " " << root->values.y << " " << root->values.z << endl;
        cout << "hups! not the same " << endl;
        return false;
    }
    
    return true;

}
template <class T>
bool SimpleKDtree<T>::sameTree(vector<Point<T>> Kdtree, int i){
    return sameTreeHelper(this->root, i, Kdtree);
}
