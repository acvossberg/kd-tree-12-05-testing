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
    cout << "Generating "<< N << " point cloud...\n";
    point.resize(N);
    for (size_t i=0;i<N;i++)
    {
        point[i].x = max_range * (rand() % 1000) / num_t(1000);
        point[i].y = max_range * (rand() % 1000) / num_t(1000);
        point[i].z = max_range * (rand() % 1000) / num_t(1000);
        point[i].ID = i;
        cout << point[i].x << ", " << point[i].y << ", " << point[i].z << " ID: " << point[i].ID << endl;

    }
    std::cout << "done\n \n";
}

template <typename num_t>
bool test(KD_tree<num_t> KDtree, SimpleKDtree<num_t>* Tsimple){
    cout << "testing ..." << endl;
    
    vector<Point<num_t>> KDVect;
    KDVect = KDtree.get_tree_as_vector();
    
    if(Tsimple->sameTree(KDVect, 1)){
        cout << "yay! trees are the same" << endl;
        return true;
    }
    else{
        cout << "o oh - debug some more!" << endl;
        
        return false;
    }
}

int main()
{
    //: <int> to generic --> DONE
    //TODO: check, which tree should be built.. cut top or bottom?
    //: make test compare to other tree -> DONE
    //TODO: free kd_tree.. ?
    //: schön alles in class machen --> DONE
    //TODO: test for 3D
    //TODO: test for larger - random trees.
    //TODO: test for double/float
    
    
    //type to use:
    typedef int num_t;
    
    vector<Point<num_t>> cloud;
    
    // Generate points:
    generateRandomPointCloud(cloud, 101);
    
    vector<int> dimensions = {1,2};
    KD_tree<num_t> tree(cloud, dimensions);
    tree.KD_tree_recursive(0, cloud.size()-1, 0, 1);
    tree.printTree();
    
    SimpleKDtree<num_t> *bst = new SimpleKDtree<num_t>(dimensions);
    bst->make_SimpleKDtree( cloud, 0, cloud.size()-1, 0);
    cout << " \n \n" << endl;
    
    //: hier testen ob beide gleich sind -> DONE
    cout << test(tree, bst) << endl;
    cloud.clear();

    return 0;
}