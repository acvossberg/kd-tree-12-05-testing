//
//  main.cpp
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 5/12/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#include <stdio.h>
#include <vector>
#include "KD_tree.h"
#include "KD_tree.cpp"
#include "simple_kd_tree.h"
#include "simple_kd_tree.cpp"
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;


template <typename num_t>
void generateRandomPointCloud(vector<Point<num_t>> &point, const size_t N, const int max_range = 10)
{
    cout << "Generating "<< N << " point cloud...";
    point.resize(N);
    for (size_t i=0;i<N;i++)
    {
        point[i].x = max_range * (rand() % 1000) / num_t(1000);
        point[i].y = max_range * (rand() % 1000) / num_t(1000);
        //point[i].z = max_range * (rand() % 1000) / T(1000);
    }
    point[0].x = 2;
    point[0].y = 3;
    point[1].x = 5;
    point[1].y = 4;
    point[2].x = 9;
    point[2].y = 6;
    point[3].x = 4;
    point[3].y = 7;
    point[4].x = 8;
    point[4].y = 1;
    point[5].x = 7;
    point[5].y = 2;
    std::cout << "done\n";
}


int main()
{
    //: <int> to generic --> DONE
    //TODO: check, which tree should be built.. cut top or bottom?
    //TODO: make test compare to other tree
    //TODO: free kd_tree.. ?
    //: schön alles in class machen --> DONE

    //type to use:
    typedef int num_t;
    
    vector<Point<num_t>> cloud;
    
    // Generate points:
    generateRandomPointCloud(cloud, 6);
    
    vector<int> dim = {1,2};
    KD_tree<num_t> tree(cloud, dim);
    tree.KD_tree_recursive(0, cloud.size()-1, 0);
    
    /*SimpleKDtree<num_t> simple_kd_tree(dim);
    KDnode<num_t> *root = simple_kd_tree.newSimpleKDtreeNode(cloud[0]);
    simple_kd_tree.make_SimpleKDtree(cloud, root, 0, cloud.size()-1, 0);*/
    
    
    
    
    SimpleKDtree<num_t> *bst = new SimpleKDtree<num_t>(dim);
    bst->make_SimpleKDtree(cloud, 0, cloud.size()-1, 0);

    
    //TODO: hier testen ob beide gleich sind
    
    
    return 0;
}