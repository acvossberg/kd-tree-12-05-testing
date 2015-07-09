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
void generateRandomPointCloud(vector<Point<num_t>> &pointn, vector<Hit<num_t>> &point, const size_t N, const int max_range = 10)
{
    cout << "Generating "<< N << " point cloud...\n";
    point.resize(N);
    pointn.resize(N);
    for (size_t i=0;i<N;i++)
    {
        pointn[i].x = max_range * (rand() % 1000) / num_t(1000);
        pointn[i].y = max_range * (rand() % 1000) / num_t(1000);
        pointn[i].z = max_range * (rand() % 1000) / num_t(1000);
        pointn[i].ID = i;
        //vector<num_t> p = {max_range * (rand() % 1000) / num_t(1000), max_range * (rand() % 1000) / num_t(1000), max_range * (rand() % 1000) / num_t(1000)};
        
        vector<num_t> p = {pointn[i].x, pointn[i].y,pointn[i].z};
        point[i].datapoints.resize(3);
        point[i].datapoints = p;
        point[i].ID = i;
        
        
        cout << point[i].datapoints[0] << ", " << point[i].datapoints[1] << ", " << point[i].datapoints[2] << " ID: " << point[i].ID << endl;
        cout << pointn[i].x << ", " << pointn[i].y << ", " << pointn[i].z << " ID: " << pointn[i].ID << endl;
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
void make_tree(vector<Point<num_t>> cloudn, vector<Hit<num_t>> cloud, vector<int> dimensions, num_t *transformable_trees, int *treesArray_ID, int offset){
    KD_tree<num_t> tree(cloudn, cloud, dimensions, transformable_trees, treesArray_ID, offset);
    tree.KD_tree_recursive(0, cloud.size()-1, 0, 1);
}

//TODO: change all the copying around.. maybe use std::move
//TODO: check speedup by changing number of threads
template <typename num_t>
void make_forest(vector<Point<num_t>> &cloudn, vector<Hit<num_t>> &cloud, vector<int> &dimensions, int datapoints_per_tree, int nthreads, num_t *transformable_trees, int *treesArray_ID){
    vector<std::future<void>> futures;
    
    int nodes_per_tree = datapoints_per_tree;
    for(int id = 0; id < nthreads; ++id){
        //TODO: auch aufsplitten - das kann jeder thread selbst tun
        //TODO: maybe way to use part of vector without copying
        
        if(id == nthreads-1){
            datapoints_per_tree = cloud.size() -  datapoints_per_tree*id;
        }
        //TODO: NO COPY!!
        vector<Hit<num_t>> threadcloud(cloud.begin()+id*datapoints_per_tree, cloud.begin()+(id+1)*datapoints_per_tree);
        vector<Point<num_t>> threadcloudn(cloudn.begin()+id*datapoints_per_tree, cloudn.begin()+(id+1)*datapoints_per_tree);
        
        
        //cout << "cloud.size() " << threadcloud.size() << " cloudn.size() " << threadcloudn.size() << endl;
        futures.push_back(std::async(launch::async, make_tree<num_t>, threadcloudn, threadcloud, dimensions, transformable_trees, treesArray_ID, id*nodes_per_tree));
        //make_tree<num_t>(threadcloudn, threadcloud, dimensions, transformable_trees, treesArray_ID, id*nodes_per_tree);
    }
    
    for(auto &e : futures) {
        e.get();
    }
    
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
vector<vector<Point<num_t>>> convertTree(num_t *tree, int *treesArray_ID, int number_of_trees, int datapoints_per_tree, int numberOfHits){
    vector<vector<Point<num_t>>> trees;
    trees.resize(number_of_trees, vector<Point<num_t>>(datapoints_per_tree));
    
    int tree_nr = -1;
    for(int i = 0; i< (number_of_trees)*datapoints_per_tree; i++){
        if(i%datapoints_per_tree == 0){
            tree_nr++;
            
        }
        trees[tree_nr][i%datapoints_per_tree].ID = treesArray_ID[i];
        trees[tree_nr][i%datapoints_per_tree].x = tree[0+3*i];
        trees[tree_nr][i%datapoints_per_tree].y = tree[1+3*i];
        trees[tree_nr][i%datapoints_per_tree].z = tree[2+3*i];
    }

    trees[number_of_trees-1].resize(pow(2,floor(log2(numberOfHits - datapoints_per_tree*(number_of_trees-1))) +1) - 1);
    return trees;
}

//check wether tree_array is correct by generating SimpleKDtree and comparing
template <typename num_t>
void test_correct_trees(num_t *trees_array_transformable, int *treesArray_ID,  int datapoints_per_tree, int threads, vector<int> dimensions, int numberOfHits, vector<Point<num_t>> &cloud){

    
    vector<vector<Point<num_t>>> trees = convertTree(trees_array_transformable, treesArray_ID, threads, datapoints_per_tree, numberOfHits);


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
    //TODO: remove Point completely
    //TODO: change 2D-array to 1D-array
    
    
    //type to use:
    typedef int num_t;
    
    vector<Hit<num_t>> cloud;
    vector<Point<num_t>> cloudn;
    // Generate points:
    int numberOfHits = 1000;
    generateRandomPointCloud(cloudn, cloud, numberOfHits);
    
    //must be defined {1, 2, 3} = {x, y, z}
    vector<int> dimensions = {1,2,3};
    int number_of_dimensions = dimensions.size();
    
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
    trees_array_transformable.resize(number_of_dimensions+1, vector<num_t>(threads*datapoints_per_tree, -1));
    
   
    //all not-used elements in treesArray are by default set to zero by compiler
    int *treesArray_ID = new int[threads*datapoints_per_tree];
    num_t *treesArray;
    treesArray = new num_t [threads*datapoints_per_tree*number_of_dimensions];//[number_of_dimensions+1][threads*datapoints_per_tree];
    
    //make 2D-array:
    /*num_t **treesArray;//[number_of_dimensions+1][threads*datapoints_per_tree];
    treesArray = new int *[number_of_dimensions];
    for(int i = 0; i <number_of_dimensions; i++){
        treesArray[i] = new int[threads*datapoints_per_tree];
    }
    */
    
    make_forest<num_t>(cloudn, cloud, dimensions, datapoints_per_tree, threads, treesArray, treesArray_ID);
    
    //test if trees made with make_forest are correct:
    test_correct_trees(treesArray,treesArray_ID, datapoints_per_tree, threads, dimensions, numberOfHits, cloudn);
    
   
    
    //make box, in which should be searched for hits
    //set all other dimensions to zero, if not used:
    int box[6] = {2, 8, 0, 0, 0, 0};
    
    Cuda_class<num_t> p;
    //p.cudaMain(threads, datapoints_per_tree, treeArray_x_new, treeArray_y_new, treeArray_z_new, treeArray_ID, box);
    p.cudaMain(threads, datapoints_per_tree, treesArray, treesArray_ID, box, number_of_dimensions);
    cloud.clear();
    return 0;
}



