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
#include <iostream>

template <typename num_t>
struct Point{
    num_t  x,y,z,ID; //dim = 1, 2, 3
};


template <class num_T>
class sorter {
private:
    short dim_;
public:
    sorter(short d) : dim_(d) {}
    
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
#endif /* defined(__kd_tree_12_05__sorter__) */
