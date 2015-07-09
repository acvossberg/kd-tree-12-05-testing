//
//  sorter.h
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 6/1/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#ifndef __kd_tree_12_05__sorter__
#define __kd_tree_12_05__sorter__

#include <stdio.h>
#include <vector>
#include <iostream>

template <typename num_t>
struct Point{
    num_t  x,y,z,ID; //dim = 1, 2, 3
};

template <typename num_t>
struct Hit{
    std::vector<num_t> datapoints;
    int ID;
};

template <class num_T>
class sortern {
private:
    short dim_;
public:
    sortern(short d) : dim_(d) {}
    
    bool operator()(Point<num_T> const &a, Point<num_T> const &b) const {
        switch (dim_){
            case 1:
                return a.x < b.x;
            case 2:
                return a.y < b.y;
            case 3:
                return a.z < b.z;
            case 4:
                return a.ID < b.ID;
            default:
                std::cout << "error - no correct dimension for comparison" << std::endl;
                return a.x < b.x;
        }
    }
};
    


template <class num_T>
class sorter {
private:
    short dim_;
public:
    sorter(short d) : dim_(d) {}
    
    bool operator()(Hit<num_T> const &a, Hit<num_T> const &b) const {
        return a.datapoints[dim_] < b.datapoints[dim_];
        /*switch (dim_){
         case 1:
         return a.x < b.x;
         case 2:
         return a.y < b.y;
         case 3:
         return a.z < b.z;
         case 4:
         return a.ID < b.ID;
         default:
         std::cout << "error - no correct dimension for comparison" << std::endl;
         return a.x < b.x;
         }*/
    }
};
template <class num_T>
class sorterID{
public:
    sorterID(){}
    bool operator()(Hit<num_T> const &a, Hit<num_T> const &b) const {
        return a.ID < b.ID;
    }
};


#endif /* defined(__kd_tree_12_05__sorter__) */
