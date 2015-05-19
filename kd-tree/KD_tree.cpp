//
//  KD_tree.cpp
//  kd-tree-08-05
//
//  Created by Ann-Christine Vossberg on 5/8/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "KD_tree.h"

using namespace std;
class sorter {
private:
    short dim_;
public:
    sorter(short dim) : dim_(dim) {}
    bool operator()(Point<int> const &a, Point<int> const &b) const {
        switch (dim_){
            case 1:
                return a.x < b.x;
            case 2:
                return a.y < b.y;
            /*case 3:
                return a.z < b.z;*/
            default:
                cout << "error - no correct dimension for comparison" << endl;
                return a.x < b.x;
        }
    }
};
template < class T>
KD_tree<T>::KD_tree(vector<Point<T>> &cloud, vector<int> &dimensions){
    data = cloud;
    dim = dimensions;
}

template <class T>
void KD_tree<T>::selectMedian(int dim, int median, int left, int right)// dim = 1, 2 oder 3
{
    cout << "left " << left << " right " << right << " med " << median << endl;

    //nth_element sorts data left - right.
    cout << median << "-th element with dim = " << dim << endl;
    
    nth_element(data.begin()+left, data.begin() + median, data.begin()+right, sorter(dim));
    cout << "median " <<  data[median].x << " " << data[median].y << endl;

    //sorts element s.t. all smaller than median on the left and larger on right

}
template<class T>
void KD_tree<T>::printTree(){
    for(int i = 0; i<data.size(); i++){
        cout << "( " << data[i].x << ", " << data[i].y << ") \t";
        
    }
    cout << endl;
}
template <class T>
void KD_tree<T>::KD_tree_recursive(int left, int right, int k){
    int med = left + (right+1 - left)/2; //left + (right - left)/2;

    //return statement:
    if(left <right){
    
        //check ob korrekt gerundet wird
        k = k%dim.size();
    
        printTree();

        selectMedian(dim[k], med, left, right+1);
        k++;
    

        //left:
        KD_tree_recursive(left, med-1, k);
        //right:
        KD_tree_recursive(med+1, right, k);
    
    }
}


