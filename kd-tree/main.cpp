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
#include <future>

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
template <typename num_t>
void print_Pointvector(vector<Point<num_t>> a){
    for(int i = 0; i<a.size(); i++){
        cout << "( "<< a[i].x << ", " << a[i].y << ", " << a[i].z << ")" << endl;
    }
    cout << "\n" << endl;
}
template <typename num_t>
void make_tree(vector<Point<num_t>> cloud, vector<int> dimensions, vector<vector<Point<num_t>>> &trees, int Id){
    KD_tree<num_t> tree(cloud, dimensions);
    tree.KD_tree_recursive(0, cloud.size()-1, 0, 1);
    trees[Id] = tree.get_tree_as_vector();
    //print_Pointvector(trees[Id]);
}

template <typename num_t>
vector<vector<Point<num_t>>> make_forest(vector<Point<num_t>> cloud,vector<int> dimensions, int datapoints_per_tree, int nthreads){
    //vector with trees inside.
    vector<vector<Point<num_t>>> trees(nthreads);//( 2, std::vector<float> ( 2, 0 ) )
    vector<std::future<void>> futures;
    std::future<void> test_future;
    
    
    for(int id = 0; id < nthreads; ++id){
        //TODO:auch aufsplitten - das kann jeder thread selbst tun
        //TODO: maybe way to use part of vector without copying
        vector<Point<num_t>> threadcloud (cloud.begin()+id*datapoints_per_tree, cloud.begin()+(id+1)*datapoints_per_tree);
        //make_tree( threadcloud, dimensions, trees, id);
        //futures.push_back(async(make_tree, threadcloud, dimensions, std::ref(trees), id));
        futures.push_back(std::async(launch::async, make_tree<num_t>, threadcloud, dimensions, std::ref(trees), id));
    }
    
    //testing - while coding: muss wieder weg
    //vector<Point<num_t>> threadcloud (cloud.begin()+0*datapoints_per_tree, cloud.begin()+(0+1)*datapoints_per_tree);
    //test_future = async(make_tree<num_t>, threadcloud, dimensions, std::ref(trees), 0);

    for(auto &e : futures) {
        e.get();
    }
    
    return trees;
}

int main()
{
    //: <int> to generic --> DONE
    //: check, which tree should be built.. cut top or bottom? ->TOP ->DONE
    //: make test compare to other tree -> DONE
    //: schön alles in class machen --> DONE
    //: test for 3D -> DONE
    //: test for larger - random trees. -> DONE
    //: test for double/float -> DONE
    //TODO: free kd_tree.. ?
    
    //type to use:
    typedef int num_t;
    
    vector<Point<num_t>> cloud;
    
    // Generate points:
    generateRandomPointCloud(cloud, 100);
    
    //must be defined {1, 2, 3} = {x, y, z}
    vector<int> dimensions = {1,2,3};
    
    //get_size_of_tree from cuda_device --> #datapoints per thread.. = datapoints per tree
    int datapoints_per_tree = 20;
    
    //make_forest() in parallel using futures, promises and tasks async!
    int threads = cloud.size()/datapoints_per_tree;
    vector<vector<Point<num_t>>> trees = make_forest<num_t>(cloud, dimensions, datapoints_per_tree, threads);
    
    
    //print
    for(int i = 0; i< trees.size(); i++){
        print_Pointvector(trees[i]);
    }
    
    
    
    /*
    KD_tree<num_t> tree(cloud, dimensions);
    tree.KD_tree_recursive(0, cloud.size()-1, 0, 1);
    tree.printTree();
    
    SimpleKDtree<num_t> *bst = new SimpleKDtree<num_t>(dimensions);
    bst->make_SimpleKDtree( cloud, 0, cloud.size()-1, 0);
    cout << " \n \n" << endl;
    
    //: hier testen ob beide gleich sind -> DONE
    cout << test(tree, bst) << endl;
    */
     
    cloud.clear();

    return 0;
}