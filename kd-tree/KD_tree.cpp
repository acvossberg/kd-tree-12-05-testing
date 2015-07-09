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
KD_tree<T>::KD_tree(vector<Point<T>> &cloud,vector<vector<T>> &data_,vector<int> &dataID_,int begin_, int end_, vector<int>& dimensions, T **transformable_trees_, int *treesArray_ID_, int& offset_) : dim(dimensions), offset(offset_), transformable_trees(transformable_trees_), treesArray_ID(treesArray_ID_), data(cloud), vectorData(data_), dataID(dataID_){
    
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
    begin = begin_;
    end = end_;
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
template < class T>
int KD_tree<T>::partition(vector<T> &input, int NodeP, int NodeR, int dimension_offset)
{
    int p = NodeP*dim.size()+dimension_offset;
    int r = NodeR*dim.size()+dimension_offset;
    int pivot = input[r];
    std::cout << "IN PARTITION: input.size " << input.size() <<" NodeP:" << NodeP << " NodeR:" << NodeR << " p:" << p << " r:" << r  << " dim-offset " << dimension_offset<< " with pivot: " << pivot <<  std::endl;
    
    
    
    while ( NodeP < NodeR )
    {
        //std::cout << " NodeP:" << NodeP << " NodeR:" << NodeR << " p:" << p << " r:" << r  << std::endl;
        while ( input[p] < pivot ){
            NodeP++;
            p+=dim.size();
            std::cout << "p is: " << p << " input[p] " << input[p] << std::endl;
        }
        while ( input[r] > pivot ){
            NodeR--;
            r-=dim.size();
            std::cout << "r is: " << r << " input[r] " << input[r] << std::endl;

        }
        /*if( input[p] == input[r] ){
            r+=dim.size();
            NodeR++;
        }*/
        if ( NodeP < NodeR ) {
            //swap dim.size number of elements:
            swap_ranges(input.begin() + p, input.begin() + p + dim.size(), input.begin() + r );
            //std::cout << "IN PARTITION" <<std::endl;
        }
    }

    return NodeR;
}

template < class T>
T KD_tree<T>::quick_select(vector<T> &input, int NodeP, int NodeR, int NodeK, int dimension_offset)
{
    int p = NodeP*dim.size()+dimension_offset;
    int r = NodeR*dim.size()+dimension_offset;
    
    std::cout << "input.size " << input.size() <<" NodeP:" << NodeP << " NodeR:" << NodeR << " p:" << p << " r:" << r  << " dim-offset " <<dimension_offset<< std::endl;
    
    
    if ( p == r ) return input[p];
    int NodeJ = partition(input, NodeP, NodeR, dimension_offset);
    int length = NodeJ - NodeP + 1;
    if ( length == NodeK ) return input[NodeJ*dim.size()+dimension_offset];
    else if ( NodeK < length ) return quick_select(input, NodeP, NodeJ - 1, NodeK, dimension_offset);
    else return quick_select(input, NodeJ + 1, NodeR, NodeK - length, dimension_offset);
}

template<class T>
void KD_tree<T>::printTree()
{
    cout << "printing tree with number of nodes: " << result.size() << endl;
    for(int i = 0; i<result.size();i++){
        cout << "(" << result[i].x << ", " << result[i].y << ", " << result[i].z << ") with ID " << result[i].ID <<  endl;
    
    }
    cout << "done" << endl;
}



template <class T>
T KD_tree<T>::qselect(int k, int li, int hi, int dimension_offset) {

    int ksmall = k*dim.size()+dimension_offset;
    
    if (hi - li <= 1){ return vectorData[ksmall];}
    int j = li;
    
    //std::swap(pArray[j], pArray[k]);
    std::swap_ranges(vectorData.begin() + j*dim.size(), vectorData.begin() + j*dim.size() + dim.size(), vectorData.begin() + k*dim.size());
    std::swap(dataID[j], dataID[k]);
    
    //cout << "not yet dumped" << endl;
    for (int i = j = li + 1; i < hi; i++){
        //if (pArray[i] < pArray[li])
        if (vectorData[i*dim.size()+dimension_offset] < vectorData[ li*dim.size()+dimension_offset]){
            //std::swap(pArray[j++], pArray[i]);
            std::swap_ranges(vectorData.begin() + j*dim.size(), vectorData.begin() + j*dim.size() + dim.size(), vectorData.begin() + i*dim.size());
            std::swap(dataID[j], dataID[i]);
            
            j++;
        }
    }
    //cout << "k " << k << " j " << j << endl;
    //std::swap(pArray[--j], pArray[li]);
    --j;
    std::swap_ranges(vectorData.begin() + j*dim.size(), vectorData.begin() + j*dim.size() + dim.size(), vectorData.begin() + li*dim.size());
    std::swap(dataID[j], dataID[li]);
    
    if (k < j) return qselect(k, li, j, dimension_offset);
    if (k > j) return qselect( k - j, j + 1, hi, dimension_offset);
    return vectorData[j];
}

template <class T>
void KD_tree<T>::original_order_median_new(int median_position, int d, int left, int right){
    //T median_value = get_value(d, data[median_position]);
    T median_value = vectorData[median_position*dim.size()+d];
    cout << median_value << " vectorData vs data " << get_value(d+1, data[median_position]) << endl;
    //make median values contiguous
    int med_left = median_position;
    int med_right = median_position;
    int leftIt = left;
    int rightIt = right-1;
    
    while(leftIt < med_left){
        //find next element left of median, that is not median
        //while(get_value(d, data[med_left]) == median_value && med_left >= left ){ med_left--; }
        cout << "data: " << get_value(d+1, data[med_left]) << " vs vectorData " << vectorData[med_left*dim.size()+d] << endl;
        while(vectorData[med_left*dim.size()+d] == median_value && med_left >= left){ med_left--; }
        //find first element from left, that is median
        //while(get_value(d, data[leftIt]) != median_value){ leftIt++;}
        while(vectorData[leftIt*dim.size()+d] != median_value){ leftIt++;}
        if(leftIt < med_left){
            //cout << "watch out - left-med switched " << leftIt << " " << med_left << endl;
            //swap(data[leftIt], data[med_left]);}
            
            swap(dataID[leftIt], dataID[med_left]);
            std::swap_ranges(vectorData.begin()+leftIt*dim.size(), vectorData.begin()+leftIt*dim.size()+dim.size(), vectorData.begin()+med_left*dim.size());
        }
    }

    if(med_right == rightIt) med_right++;
    while(med_right < rightIt){
        //find next element left of median, that is not median
        //while(get_value(d, data[med_right]) == median_value && med_right < right){ med_right++; }
        while(vectorData[med_right*dim.size()+d] == median_value && med_right < right){ med_right++; }
        //find first element from left, that is median
        //while(get_value(d, data[rightIt]) != median_value){ rightIt--;}
        while(vectorData[rightIt*dim.size()+d] != median_value){ rightIt--;}
            
        if(med_right < rightIt){
            //cout << "watch out - right-med switched " << rightIt << " " << med_right << endl;
            //std::swap(data[rightIt], data[med_right]);
                
            std::swap(dataID[rightIt], dataID[med_right]);
            std::swap_ranges(vectorData.begin()+rightIt*dim.size(), vectorData.begin()+rightIt*dim.size()+dim.size(), vectorData.begin()+med_right*dim.size());
        }
    }
        
    //sort(data.begin()+leftIt, data.begin()+med_right, sorter<T>(4));
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
    //original_order_median_new(median, d-1, left, right);
    
    result[pos-1] = data[median];
    
    //transformable_trees[0][offset+pos-1] = data[median].ID;
    treesArray_ID[offset+pos-1] = data[median].ID;
    transformable_trees[0][offset+pos-1] = data[median].x;
    transformable_trees[1][offset+pos-1] = data[median].y;
    transformable_trees[2][offset+pos-1] = data[median].z;
    
    //checking if all data[].ID == dataID TODO: take away!!!
    /*for(int i = 0; i<data.size(); i++){
        if(data[i].ID != dataID[i]){cout << " not the same" << endl;}
        else{cout << "THE SAME!!!!" << endl;}
    }*/
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
    
            swap(data[leftIt], data[med_left]);
        }
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

