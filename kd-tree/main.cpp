//
//  main.cpp
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 5/12/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#include <stdio.h>
#include <vector>
#include "KD_tree.hpp"
#include "KD_tree.cpp"
#include "simple_kd_tree.hpp"
#include "simple_kd_tree.cpp"
#include "cuda_runtime.h"
#include "cuda.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <future>
#include "InsideBox.hpp"

#define MYDEVICE 0

using namespace std;


template <typename num_t>
void generateRandomPointCloud( vector<Point<num_t>> &point, const size_t N, const int max_range = 10)
{
    cout << "Generating "<< N << " point cloud...\n";
    point.resize(N);
    for (size_t i=0;i<N;i++)
    {
        point[i].x = max_range * (rand() % 1000) / num_t(1000);
        point[i].y = max_range * (rand() % 1000) / num_t(1000);
        point[i].z = max_range * (rand() % 1000) / num_t(1000);
        point[i].ID = i;
        //cout << point[i].x << ", " << point[i].y << ", " << point[i].z << " ID: " << point[i].ID << endl;
        
    }
    std::cout << "done\n \n";
}

template <typename num_t>
bool testing_trees(vector<Point<num_t>> treevector, SimpleKDtree<num_t>* Tsimple){
    
    if(Tsimple->sameTree(treevector, 1)){
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
void make_tree(vector<Point<num_t>> cloud, vector<int> dimensions, vector<vector<num_t>>& transformable_trees, int offset){
    KD_tree<num_t> tree(cloud, dimensions, transformable_trees, offset);
    tree.KD_tree_recursive(0, cloud.size()-1, 0, 1);
    //trees[Id] = tree.get_tree_as_vector();
}

//TODO: change all the copying around.. maybe use std::move
//TODO: check speedup by changing number of threads
template <typename num_t>
void make_forest(vector<Point<num_t>> &cloud,vector<int> &dimensions, int datapoints_per_tree, int nthreads, vector<vector<num_t>> &transformable_trees){
    //vector<vector<Point<num_t>>> trees(nthreads);
    vector<std::future<void>> futures;
    
    int nodes_per_tree = datapoints_per_tree;
    for(int id = 0; id < nthreads; ++id){
        //TODO: auch aufsplitten - das kann jeder thread selbst tun
        //TODO: maybe way to use part of vector without copying
        
        if(id == nthreads-1){
            datapoints_per_tree = cloud.size() -  datapoints_per_tree*id;
        }
        //cout << "offset in main: " << id*nodes_per_tree << " id: " << id << " nodes_per_tree: " << nodes_per_tree << endl;
        vector<Point<num_t>> threadcloud(cloud.begin()+id*datapoints_per_tree, cloud.begin()+(id+1)*datapoints_per_tree);
        futures.push_back(std::async(launch::async, make_tree<num_t>, threadcloud, dimensions, std::ref(transformable_trees), id*nodes_per_tree));
        
    }
    
    for(auto &e : futures) {
        e.get();
    }
    
    //return trees;
}

template <typename num_t>
vector<int> inBox(Point<num_t> start, Point<num_t> end, vector<vector<Point<num_t>>> &trees ){
    vector<int> result;
    
    for(int i=0; i<trees.size(); i++){
        for(int j = 0; j< trees[i].size();j++){
            
            if( ((trees[i][j].x >= start.x && trees[i][j].x <= end.x) || (end.x == 0 && start.x == 0))  && ((trees[i][j].y >= start.y && trees[i][j].y <= end.y) || (end.y == 0 && start.y == 0)) && ((trees[i][j].z >= start.z && trees[i][j].z <= end.z) || (end.z == 0 && start.z == 0))){
                result.push_back(trees[i][j].ID);
            }
        
        }
    
    }
    return result;
}

template<typename num_t>
vector<vector<Point<num_t>>> convertTree(vector<vector<num_t>> &tree, int number_of_trees, int datapoints_per_tree, int numberOfHits){
    vector<vector<Point<num_t>>> trees;
    trees.resize(number_of_trees, vector<Point<num_t>>(datapoints_per_tree));
    cout << "Number of trees: " << trees.size()<< endl;
    
    int tree_nr = -1;
    for(int i = 0; i< (number_of_trees)*datapoints_per_tree; i++){
        if(i%datapoints_per_tree == 0){
            tree_nr++;
            
        }
        //cout << "tree_nr: "<< tree_nr << " i " << i%datapoints_per_tree << endl;
        trees[tree_nr][i%datapoints_per_tree].ID = tree[0][i];
        trees[tree_nr][i%datapoints_per_tree].x = tree[1][i];
        trees[tree_nr][i%datapoints_per_tree].y = tree[2][i];
        trees[tree_nr][i%datapoints_per_tree].z = tree[3][i];
    }
    /*for(int i = 0; i<datapoints_per_tree; i++){
        cout << "trees[number_of_trees-1][i].x " << trees[number_of_trees-1][i].x << endl;
    }*/
    
    trees[number_of_trees-1].resize(pow(2,floor(log2(numberOfHits - datapoints_per_tree*(number_of_trees-1))) +1) - 1);
    return trees;
}

//check wether tree_array is correct by generating SimpleKDtree and comparing
template <typename num_t>
vector<vector<Point<num_t>>> test_correct_trees(vector<vector<num_t>> trees_array_transformable, int datapoints_per_tree, int threads, vector<int> dimensions, int numberOfHits, vector<Point<num_t>> &cloud){

    
    //make vector<vector<Point<num_t>>> trees from vector<vector<num_t>> s.t. comparable to SimpleKDtree:
    vector<vector<Point<num_t>>> trees = convertTree(trees_array_transformable, threads, datapoints_per_tree, numberOfHits);


    bool correctTree=true;
    int points_in_tree = datapoints_per_tree;
    for(int i = 0; i < threads; i++){
        SimpleKDtree<num_t> *bst = new SimpleKDtree<num_t>(dimensions);
        if(i == threads-1){
            points_in_tree = numberOfHits - datapoints_per_tree*i;
        }
        vector<Point<num_t>> threadcloud (cloud.begin()+i*points_in_tree, cloud.begin()+(i+1)*points_in_tree);
        bst->make_SimpleKDtree(threadcloud, 0, threadcloud.size()-1, 0);
        correctTree = correctTree && testing_trees(trees[i], bst);
        delete bst;
    }
    if(correctTree){cout << "\nAll tree's are correct from tree_arrays!!!!" << endl; }
    return trees;
}

void printDevProp(cudaDeviceProp devProp)
{
    //http://www.cs.fsu.edu/~xyuan/cda5125/examples/lect24/devicequery.cu
    printf("Major revision number:         %d\n",  devProp.major);
    printf("Minor revision number:         %d\n",  devProp.minor);
    printf("Name:                          %s\n",  devProp.name);
    printf("Total global memory:           %u\n",  devProp.totalGlobalMem);
    printf("Total shared memory per block: %u\n",  devProp.sharedMemPerBlock);
    printf("Total registers per block:     %d\n",  devProp.regsPerBlock);
    printf("Warp size:                     %d\n",  devProp.warpSize);
    printf("Maximum memory pitch:          %u\n",  devProp.memPitch);
    printf("Maximum threads per block:     %d\n",  devProp.maxThreadsPerBlock);
    for (int i = 0; i < 3; ++i)
        printf("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
    for (int i = 0; i < 3; ++i)
        printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
    printf("Clock rate:                    %d\n",  devProp.clockRate);
    printf("Total constant memory:         %u\n",  devProp.totalConstMem);
    printf("Texture alignment:             %u\n",  devProp.textureAlignment);
    printf("Concurrent copy and execution: %s\n",  (devProp.deviceOverlap ? "Yes" : "No"));
    printf("Number of multiprocessors:     %d\n",  devProp.multiProcessorCount);
    printf("Kernel execution timeout:      %s\n",  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
    
    return;
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
    //: free kd_tree.. ?
    
    //type to use:
    typedef int num_t;
    
    vector<Point<num_t>> cloud;
    
    // Generate points:
    int numberOfHits = 1000;
    generateRandomPointCloud(cloud, numberOfHits);
    
    //must be defined {1, 2, 3} = {x, y, z}
    vector<int> dimensions = {1,2,3};
    int number_of_dimensions = 3;
    
    //get_size_of_tree from cuda_device --> #datapoints per thread.. = datapoints per tree
    int device;
    cudaGetDevice(&device);
    std::cout<< "devices are: " << device << endl;
    
    cudaDeviceProp devProp;
    cudaGetDeviceProperties(&devProp, device);
    printDevProp(devProp);
    int max_threads = devProp.warpSize;
    //TODO: here calculate #nodes need
    
    cout << "number of warps " << max_threads << endl;
    
    int warp_size = max_threads;
    
    //???: will it be one tree per warp or per block?
    //formula: 2^(floor(log_2(64)+1))-1 is always max when one less than warp_size
    int datapoints_per_tree = warp_size-1;
    cout << "datapoints_per_tree: " << datapoints_per_tree << endl;
    
    
    //round up: q = (x + y - 1) / y;
    int threads = (numberOfHits+datapoints_per_tree-1)/datapoints_per_tree;
    vector<vector<num_t>> trees_array_transformable;
    trees_array_transformable.resize(number_of_dimensions+1, vector<num_t >(threads*datapoints_per_tree, -1));
    
    std::cout << "  trees_array_transformable.size() " << trees_array_transformable.size() << "  trees_array_transformable[0].size() " << trees_array_transformable[0].size() << endl;

    
    make_forest<num_t>(cloud, dimensions, datapoints_per_tree, threads, trees_array_transformable);
    //cout << "Number of trees: " << trees.size()<< endl;
    
    
    
    int* treeArray_ID_new = &trees_array_transformable[0][0];
    int* treeArray_x_new = &trees_array_transformable[1][0];
    int* treeArray_y_new = &trees_array_transformable[2][0];
    int* treeArray_z_new = &trees_array_transformable[3][0];
    
    
    /*for(int i = (threads-1)*datapoints_per_tree; i < threads*datapoints_per_tree; i++){
        std::cout << trees_array_transformable[1][i] << endl;
    }*/
    
    
    
    
    
    
    /*
    
    bool correctTree=true;
    int points_in_tree = datapoints_per_tree;
    cout << "threadcloud is made for " << 0 << " till "<< threads-1 << endl;
    for(int i = 0; i < threads; i++){
        SimpleKDtree<num_t> *bst = new SimpleKDtree<num_t>(dimensions);
        if(i == threads-1){
            points_in_tree = numberOfHits-  datapoints_per_tree*i;
        }
        vector<Point<num_t>> threadcloud (cloud.begin()+i*points_in_tree, cloud.begin()+(i+1)*points_in_tree);
        bst->make_SimpleKDtree(threadcloud, 0, threadcloud.size()-1, 0);
        correctTree = correctTree && testing_trees(trees[i], bst);
        delete bst;
    }
    if(correctTree){cout << "\nAll tree's are correct" << endl; }
    */
    
    
    vector<vector<Point<num_t>>> new_tree = test_correct_trees(trees_array_transformable, datapoints_per_tree, threads, dimensions, numberOfHits, cloud);
    
    
    //std::cout << "tree.size() " << trees.size() << " vs new: " << new_tree.size() << endl;
    //std::cout << "tree[0].size() " << trees[0].size() << " vs new: " << new_tree[0].size() << endl;
    
    
    
       // print_Pointvector(trees[trees.size()-1]);
        //cout << "vs new tree: \n" << endl;
        //print_Pointvector(new_tree[trees.size()-1]);
        


    
    //make trees into array (instead vector<vector< >> and copy this array over
    //: should be done while making trees and not converted afterwards ---> DONE
    /*
    int* treeArray_x = new int[trees.size()*trees[0].size()];
    int* treeArray_y = new int[trees.size()*trees[0].size()];
    int* treeArray_z = new int[trees.size()*trees[0].size()];
    int* treeArray_ID = new int[trees.size()*trees[0].size()];
    
    
    for(int i=0; i< trees.size() ; i++){
        for(int j = 0; j < trees[i].size(); j++){
            treeArray_x[i*trees[i].size()+j] = trees[i][j].x;
            treeArray_y[i*trees[i].size()+j] = trees[i][j].y;
            treeArray_z[i*trees[i].size()+j] = trees[i][j].z;
            treeArray_ID[i*trees[i].size()+j] = trees[i][j].ID;
        }
    }
    */
    
     
    
    //int size_of_forest = sizeof(int)*trees.size()*trees[0].size();
    
    //check array: - wieder weg!
    /*for(int i = 0; i < 992; i++){
        cout << treeArray_x[i] << endl;
    }*/
    
    //make box, in which should be searched for hits
    //set all other dimensions to zero, if not used:
    int box[6] = {2, 8, 0, 0, 0, 0};
    
    /*Cuda_class<num_t> p;
    p.cudaMain(threads, datapoints_per_tree, treeArray_x_new, treeArray_y_new, treeArray_z_new, treeArray_ID_new, box);
    //cudaMain<int>(trees.size(), trees[0].size(), treeArray_x, treeArray_y, treeArray_z, treeArray_ID, box);
    */
    cloud.clear();
    
    return 0;
}



