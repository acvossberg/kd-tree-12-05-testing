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
#include "KD_tree.hpp"

using namespace std;

template < class T>
KD_tree<T>::KD_tree(vector<Point<T>> &cloud, vector<int> dimensions, vector<vector<T>> &transformable_trees_) : dim(dimensions){
    
    transformable_trees = transformable_trees_;
    //data is reference??? or copy??
    data = cloud;
    int height = floor(log2(data.size()));
    int max_number_nodes = pow(2,height+1) - 1;
    Point<T> placeholder;
    placeholder.x = -1; //TODO: some better marker/placeholder needed!
    placeholder.y = -1;
    placeholder.z = -1;
    result.resize(max_number_nodes, placeholder);
}
template <class T>
vector<Point<T>> KD_tree<T>::get_tree_as_vector(){
    return result;
}
template < class T>
T KD_tree<T>::get_value(int d, Point<T> val ){
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
template<class T>
void KD_tree<T>::printTree(){
    cout << "printing tree with number of nodes: " << result.size() << endl;
    for(int i = 0; i<result.size();i++){
        cout << "(" << result[i].x << ", " << result[i].y << ", " << result[i].z << ") with ID " << result[i].ID <<  endl;
    
    }
    cout << "done" << endl;
}
template <class T>
void KD_tree<T>::selectMedian(int d, int median, int left, int right, int pos)// dim = 1, 2 oder 3
{
    //nth_element sorts data left - right.
    //sorts element s.t. all smaller than median on the left and larger on right
    nth_element(data.begin()+left, data.begin() + median, data.begin()+right, sorter<T>(d));
    
    //cout << "after sorted with nth_element: von " << left << " bis " <<  right << " with dim= " << d << endl;
    //printData();
    
    //bring duplicate values of median in original order
    original_order_median(median, d, left, right);
    
    result[pos-1] = data[median];
}

template <class T>
void KD_tree<T>::original_order_median(int median_position, int d, int left, int right){
    T median_value = get_value(d, data[median_position]);
    
    //make median values contiguous
    int med_left = median_position;
    int med_right = median_position;
    int leftIt = left;
    int rightIt = right-1;
    
    while(leftIt < med_left){
        //find next element left of median, that is not median
        while(get_value(d, data[med_left]) == median_value && med_left >= left ){ med_left--; }
        //find first element from left, that is median
        while(get_value(d, data[leftIt]) != median_value){ leftIt++;}
        
        if(leftIt < med_left){
            //cout << "watch out - left-med switched " << leftIt << " " << med_left << endl;
    
            swap(data[leftIt], data[med_left]);}
    }
    
    if(med_right == rightIt) med_right++;
    while(med_right < rightIt){
        //find next element left of median, that is not median
        while(get_value(d, data[med_right]) == median_value && med_right < right){ med_right++; }
        //find first element from left, that is median
        while(get_value(d, data[rightIt]) != median_value){ rightIt--;}
        
        if(med_right < rightIt){
            //cout << "watch out - right-med switched " << rightIt << " " << med_right << endl;
            swap(data[rightIt], data[med_right]);
        }
    }
    
    sort(data.begin()+leftIt, data.begin()+med_right, sorter<T>(4));
}

template < class T>
void KD_tree<T>::printData(){
    for(int i = 0; i< data.size(); i++){
        cout << data[i].x << ", " << data[i].y << ", " << data[i].z << " ID: " << data[i].ID << endl;
    }
    cout << "\n" << endl;
}
template <class T>
void KD_tree<T>::KD_tree_recursive(int left, int right, int k, int pos){
    int med = left + (right+1 - left)/2; //left + (right - left)/2;
    //if even:
    if(left <= right){
    
        //check ob korrekt gerundet wird
        k = k%dim.size();
        int d = dim[k];
        selectMedian( d, med, left, right+1, pos);
        k++;
    
        //left:
        pos *= 2;
        KD_tree_recursive(left, med-1, k, pos);
        pos+=1;
        //right:
        KD_tree_recursive(med+1, right, k, pos);
    
    }
    
}

