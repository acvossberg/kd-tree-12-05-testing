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
    KDnode<T> *node;
    node->values = data ;
    node->left = NULL;
    node->right = NULL;
    
    return node;
}
template<class T>
bool SimpleKDtree<T>::same_value(Point<T> &a, Point<T> &b, int d){
    switch (d){
        case 1:
            return a.x == b.x;
        case 2:
            return a.y == b.y;
        case 3:
            return a.z == b.z;
        default:
            cout << "error at == - no correct dimension for comparison" << endl;
            return a.x == b.x;
            
    }
}
template < class T>
SimpleKDtree<T>::SimpleKDtree(vector<int> &dimensions){
    dim = dimensions;
    
}
template < class T>
void SimpleKDtree<T>::make_SimpleKDtreeHelper(vector<Point<T>> cloud, KDnode<T> *root, int left, int right, int k){
    int mid = cloud.size()/2;
    height--;
    number_nodes--;
    if(cloud.size() > 1){
        //check ob korrekt gerundet wird
        k = k%dim.size();
        
        //printTree();
        //was macht sort genau mit gleichen values?
        stable_sort(cloud.begin(), cloud.end(), sorter(dim[k]));
        
        int first = mid-1;
        while(same_value(cloud[mid], cloud[first], dim[k])){
            cout << "same value" << endl;
            first--;
        }
        int last = mid+1;
        while(same_value(cloud[mid], cloud[last], dim[k])){
            cout << "same value" << endl;
            last++;
        }
        //int ID_mid = first + (last - first)/2;
        //cout << "ID " << ID_mid << " mid " << mid << endl;
        //this->values = cloud[mid];
        
        
        if(root){
            root->values = cloud[mid];
            //stimmt sowieso nicht... ID_mid.. denn cloud wird jede iteration gekürzt!!!!
            cout << "(" << cloud[mid].x << ", " << cloud[mid].y << ", " << cloud[mid].z << ")" << endl;
        }
        else{
            root = new KDnode<T>(cloud[ID_mid]);
        }
        
        vector<Point<T>> cloud_left(cloud.begin(), cloud.begin() + ID_mid);
        vector<Point<T>> cloud_right(cloud.begin()+ID_mid+1, cloud.end());
        k++;
        root->left = new KDnode<T>;
        root->right = new KDnode<T>;
        //left:
        make_SimpleKDtreeHelper(cloud_left,root->left, left, ID_mid-1, k);
        //right:
        height++;
        make_SimpleKDtreeHelper(cloud_right,root->right, ID_mid+1, right, k);
        height++;
    }
    if(cloud.size() == 1){
        root->values = cloud[0];
        //if(number_nodes <= 3){
        if(height == 0){
            Point<T> zero_point;
            zero_point.x = zero_point.y = zero_point.z = 0;
            cout << "put extra children here: " << cloud[0].x << " " << cloud[0].y << endl;
            root->left = new KDnode<T>(zero_point);
            root->right = new KDnode<T>(zero_point);
        }
    }

}

template <class T>
void SimpleKDtree<T>::make_SimpleKDtree(vector<Point<T>> cloud, int left, int right, int k){
    height = floor(log2(cloud.size()));
    number_nodes = pow(2, height + 1) - 1;
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
        cout << "pair " << i-1 << " is the same! (" << Kdtree[i-1].x << ", " << Kdtree[i-1].y << ")" << endl;
        return sameTreeHelper(root->left, 2*i, Kdtree) && sameTreeHelper(root->right, 2*i+1, Kdtree);
    }
    else if(root->values.x == 0 && root->values.y == 0 && Kdtree[i-1].x == -1 && Kdtree[i-1].y == -1){
        //end-node
        return true;
    }
    else{
        cout << "hups! not the same at " << i << " root " << root->values.x << " " << root->values.y << " Kdtree " << Kdtree[i-1].x << " " << Kdtree[i-1].y << endl;
        return false;
    }
    
    return true;

}
template <class T>
bool SimpleKDtree<T>::sameTree(vector<Point<T>> Kdtree, int i){
    return sameTreeHelper(this->root, i, Kdtree);
}
