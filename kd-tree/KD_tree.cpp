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
KD_tree<T>::KD_tree(vector<Point<T>> &cloudn, vector<Hit<T>> &cloud, vector<int>& dimensions, T *transformable_trees_, int *treesArray_ID_, int& offset_) : dim(dimensions), offset(offset_), transformable_trees(transformable_trees_), treesArray_ID(treesArray_ID_), data(cloud), datan(cloudn){
    
    //offset = offset_;
    //transformable_trees = transformable_trees_;
    //data is reference??? or copy??
    //data = cloud;
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
    nth_element(datan.begin()+left, datan.begin() + median, datan.begin()+right, sortern<T>(d+1));
    
    //cout << "after sorted with nth_element: von " << left << " bis " <<  right << " with dim= " << d << endl;
    
    
    //bring duplicate values of median in original order
    original_order_medianN(median, d+1, left, right);
    original_order_median(median, d, left, right);
    //printData();
    
    result[pos-1] = datan[median];
    
    
    //insert hit into KDtree:
    for(int i = 0; i<dim.size();i++){
        transformable_trees[offset*dim.size()+(pos-1)*dim.size()+i] = data[median].datapoints[i];
    }
    treesArray_ID[offset+pos-1] = data[median].ID;
    
}
template <class T>
void KD_tree<T>::original_order_medianN(int median_position, int d, int left, int right){
    T median_value = get_value(d, datan[median_position]);
    //make median values contiguous
    int med_left = median_position;
    int med_right = median_position;
    int leftIt = left;
    int rightIt = right-1;
    
    //cout << "med_left: " << med_left << " leftIt " << leftIt << " med_right: " << med_right << " rightIt " << rightIt << endl;
    //cout << "size of data: " << datan.size() << endl;
    while(leftIt < med_left){
        //find next element left of median, that is not median
        while(get_value(d, datan[med_left]) == median_value && med_left >= left ){ med_left--; }
        //find first element from left, that is median
        while(get_value(d, datan[leftIt]) != median_value){ leftIt++;}
        
        if(leftIt < med_left){
            //cout << "watch out - left-med switched " << leftIt << " " << med_left << endl;
            swap(datan[leftIt], datan[med_left]);
        }
    }
    if(med_right == rightIt) med_right++;
    while(med_right < rightIt){
        //find next element left of median, that is not median
        while(get_value(d, datan[med_right]) == median_value && med_right < right){ med_right++; //cout << "med_right: " << med_right << " and value " << get_value(d, datan[med_right]) << endl;
        }
        //find first element from left, that is median
        while(get_value(d, datan[rightIt]) != median_value){ rightIt--;}
        
        if(med_right < rightIt){
            //cout << "watch out - right-med switched " << rightIt << " " << med_right << endl;
            swap(datan[rightIt], datan[med_right]);
        }
    }
    
    sort(datan.begin()+leftIt, datan.begin()+med_right, sortern<T>(4));
}

template <class T>
void KD_tree<T>::original_order_median(int median_position, int d, int left, int right){
    //cout << "d " << d << " left " << left << " right " << right << endl;
    T median_value = data[median_position].datapoints[d];
    //make median values contiguous
    int med_left = median_position;
    int med_right = median_position;
    int leftIt = left;
    int rightIt = right-1;
    //cout << "med_left: " << med_left << " leftIt " << leftIt << " med_right: " << med_right << " rightIt " << rightIt << endl;

    //cout << "size of data: " << data.size() << endl;
    while(leftIt < med_left){
        //find next element left of median, that is not median
        while(data[med_left].datapoints[d] == median_value && med_left >= left ){
            med_left--;
            if(med_left < 0){break;}
            //cout << "med_left: " << med_left << endl;
        }
        //find first element from left, that is median
        while(data[leftIt].datapoints[d] != median_value){ leftIt++; //cout << "leftIt: " << leftIt << endl;
        }
        
        if(leftIt < med_left){
            //cout << "watch out - left-med switched " << leftIt << " " << med_left << endl;
            swap(data[leftIt], data[med_left]);
        }
    }
    if(med_right == rightIt) med_right++;
    while(med_right < rightIt){
        //find next element left of median, that is not median
        while(data[med_right].datapoints[d] == median_value && med_right < right){
            med_right++;
            if(med_right > data.size()-1){break;}
             //cout << "med_right: " << med_right << " and value " << data[med_right].datapoints[d] << endl;
        }
        //find first element from left, that is median
        while(data[rightIt].datapoints[d] != median_value){ rightIt--;//cout << "rightIt: " << rightIt << endl;
        }
        
        if(med_right < rightIt){
            //cout << "watch out - right-med switched " << rightIt << " " << med_right << endl;
            swap(data[rightIt], data[med_right]);
        }
    }
    
    sort(data.begin()+leftIt, data.begin()+med_right, sorterID<T>());
}

template < class T>
void KD_tree<T>::printData(){
    for(int i = 0; i< data.size(); i++){
        cout <<"data " << data[i].datapoints[0] << ", " << data[i].datapoints[1] << ", " << data[i].datapoints[2] << " ID: " << data[i].ID << endl;
        cout << "datan " << datan[i].x << ", " << datan[i].y << ", " << datan[i].z << " ID: " << datan[i].ID << endl;    }
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
        selectMedian( d-1, med, left, right+1, pos);
        k++;
    
        //left:
        pos *= 2;
        KD_tree_recursive(left, med-1, k, pos);
        pos+=1;
        //right:
        KD_tree_recursive(med+1, right, k, pos);
    
    }
    
}

