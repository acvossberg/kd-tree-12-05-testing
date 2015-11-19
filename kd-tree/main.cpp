//
//  main.cpp
//  kd-tree-12-05
//
//  Created by Ann-Christine Vossberg on 5/12/15.
//  Copyright (c) 2015 Ann-Christine Vossberg. All rights reserved.
//

#include <stdio.h>
#include <vector>
#include <iterator>
#include "KD_tree.hpp"
#include "KD_tree.cpp"
#include "simple_kd_tree.hpp"
#include "simple_kd_tree.cpp"
//#include "cuda_runtime.h"
//#include "cuda.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <future>
//#include "InsideBox.hpp"
#include <fstream>
#include <iomanip>
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
        pointn[i].ID = i+1;
        //vector<num_t> p = {max_range * (rand() % 1000) / num_t(1000), max_range * (rand() % 1000) / num_t(1000), max_range * (rand() % 1000) / num_t(1000)};
        
        vector<num_t> p = {pointn[i].x, pointn[i].y,pointn[i].z};
        point[i].datapoints.resize(3);
        point[i].datapoints = p;
        point[i].ID = i+1;
        
        
        //cout << point[i].datapoints[0] << ", " << point[i].datapoints[1] << ", " << point[i].datapoints[2] << " ID: " << point[i].ID << endl;
        //cout << pointn[i].x << ", " << pointn[i].y << ", " << pointn[i].z << " ID: " << pointn[i].ID << endl;
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

//TODO: change all the copying around.. maybe use std::move --> DONE through references and pointers
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
            cout << "datapoints_per_tree for last tree: " << datapoints_per_tree << endl;
        }
        
        //TODO: NO COPY!!
        vector<Hit<num_t>> threadcloud(cloud.begin()+id*datapoints_per_tree, cloud.begin()+(id+1)*datapoints_per_tree);
        vector<Point<num_t>> threadcloudn(cloudn.begin()+id*datapoints_per_tree, cloudn.begin()+(id+1)*datapoints_per_tree);
        
        
        //cout << "id " << id << " von " << id*datapoints_per_tree << " bis " << (id+1)*datapoints_per_tree << endl;
        futures.push_back(std::async(launch::async, make_tree<num_t>, threadcloudn, threadcloud, dimensions, transformable_trees, treesArray_ID, id*nodes_per_tree));
        //make_tree<num_t>(threadcloudn, threadcloud, dimensions, transformable_trees, treesArray_ID, id*nodes_per_tree);
    }
    
    for(auto &e : futures) {
        e.get();
    }
    
}

template <typename num_t>
vector<int> inBox(int threads, int datapoints_per_tree, int box[], vector<vector<Point<num_t>>> &trees ){
    vector<int> result;
    
    for(int i=0; i<threads; i++){
        for(int j=0; j< datapoints_per_tree; j++){
            //cout << trees[i][j].x << " " << trees[i][j].y << " " << trees[i][j].z << " ID:" << trees[i][j].ID <<  endl;
            if( box[0] <= trees[i][j].x && box[1]>=trees[i][j].x && box[2] <=trees[i][j].y && box[3] >= trees[i][j].y && box[4] <= trees[i][j].z && box[5] >= trees[i][j].z){
                result.push_back(trees[i][j].ID);
            }
            else{ result.push_back(-1);}
        }
    }
    /*
    for(int i=0; i<trees.size(); i++){
        for(int j = 0; j< trees[i].size();j++){
            
            if( (trees[i][j].x >= start.x && trees[i][j].x <= end.x)  && (trees[i][j].y >= start.y && trees[i][j].y <= end.y) && (trees[i][j].z >= start.z && trees[i][j].z <= end.z) ){
                    result.push_back(trees[i][j].ID);
            }
            else{
                result.push_back(-1);
            }
        
        }
    
    }*/
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
vector<vector<Point<num_t>>> test_correct_trees(num_t *trees_array_transformable, int *treesArray_ID,  int datapoints_per_tree, int threads, vector<int> dimensions, int numberOfHits, vector<Point<num_t>> &cloud){

    //convertTree only for 3D tree. Convers 1D tree to trees.x, trees.y, trees.z, trees.ID
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
    return trees;
}

/*void printDevProp(cudaDeviceProp devProp)
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
}*/
int number_of_nodes_traversed = 0;
void traverseTreeCPU( vector<int> &treeArray_values, vector<int> &treeArray_ID, vector<int> &results, vector<int> &box, int pos, int startOfTree, int endOfTree, int number_of_dimensions){
    
    number_of_nodes_traversed++;

    if(startOfTree + pos - 1 <= endOfTree){
        
        //calculate which tree level we are on to know which dimension was sorted
        int level = ceil(log2(double(pos+1))-1);
        int level_of_dimension = level%number_of_dimensions;
        int lastLevel = ceil(log2(double(endOfTree+1))-1);
        //a mod b = a - floor(a / b) * b
        
        //cout << "global position: " << startOfTree+pos -1 << " und ID: " << treeArray_ID[startOfTree+pos-1] << " level: " << level <<  " levelofDimension: " << level_of_dimension << endl;
        
        
        //check wether invalid encountered, continue search:
        if(treeArray_ID[startOfTree+pos-1] != 0){
            
            //if node has sorted dimension in box, continue both branches:
            //cout << box[2*level_of_dimension] <<  " <= " << treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] << " <= " << box[2*level_of_dimension+1] << endl;
            
            if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] >= box[2*level_of_dimension] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] <= box[2*level_of_dimension+1]){
                bool inside = true;
                //cout << "number of nodes traversed " << number_of_nodes_traversed << endl;
                
                //check wether the node is inside the box:
                for(int i=0; i<number_of_dimensions; i++){
                    if(i == level_of_dimension) continue;
                    //cout << treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] << " >= " << box[2*i] <<  " && " << treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] << " <= " << box[2*i+1] << endl;
                    
                    if( treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] >= box[2*i] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] <= box[2*i+1] ){
                        //entirely inside box for all dimensions
                        
                    }
                    else{
                        //not totally inside box
                        //printf("\n thread %d is changing ID %d of tree starting at %d, exact position: %d", threadIdx.x, treeArray_ID[startOfTree+pos-1], startOfTree, startOfTree+pos-1);
                        //cout << "not inside \t ID: " << treeArray_ID[startOfTree+pos-1] <<" bei pos: " << startOfTree+pos-1 << " weg! BEIDE BRANCHES" << endl;
                        inside = false;
                    }
                }
                if(inside){
                    results[startOfTree+pos-1] = treeArray_ID[startOfTree+pos-1];
                    //cout << "yes inside \t ID: " << treeArray_ID[startOfTree+pos-1] << " BEIDE BRANCHES" << endl;
                }
                
                //continue both branches:
                if(level != lastLevel){
                    //cout << "ID: " << treeArray_ID[startOfTree+pos-1] <<" bei pos: " << startOfTree+pos-1 << " BEIDE BRANCHES" << endl;
                    //left child:
                    pos *= 2;
                    traverseTreeCPU(treeArray_values, treeArray_ID,results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                
                    //right child:
                    pos += 1;
                    traverseTreeCPU(treeArray_values, treeArray_ID,results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            
            }
            //if sorted dimension is larger than box follow branch of smaller child = left child
            else if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] > box[2*level_of_dimension+1]){
                //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " LINKES KIND" << endl;
                //cout << "number of nodes traversed " << number_of_nodes_traversed << endl;
                //left child:
                if(level != lastLevel){
                    pos *= 2;
                    traverseTreeCPU(treeArray_values, treeArray_ID,results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            }
            //if sorted dimension is smaller than box, follow branch of larger child = right child
            else{
                //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " RECHTES KIND" << endl;
                //cout << "number of nodes traversed " << number_of_nodes_traversed << endl;
                if(level != lastLevel){
                    //right child:
                    pos *= 2;
                    pos += 1;
                    traverseTreeCPU(treeArray_values, treeArray_ID, results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            }
        }
    }
}

void insideBoxCPU(vector<int> &treeArray_values, vector<int> &treeArray_ID, vector<int> &results, vector<int> &box, int tree_size, int number_of_dimensions){
    
    int threadIdx = 0;
    int startOfTree = threadIdx* tree_size;
    int endOfTree = startOfTree + (tree_size - 1);
    
    traverseTreeCPU(treeArray_values, treeArray_ID,results, box, 1, startOfTree, endOfTree, number_of_dimensions);
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
    vector<int> VectNumberOfHits = {100,100,100,1000,2000,3000, 4000, 5000, 10000,15000, 20000,30000};
    
    ofstream myThreadFile ("ThreadingTimes.txt");
    ofstream myCudaFile("CudaTimes.txt");
    ofstream myTreeSizeFile("TreeSizes.txt");
    ofstream myNodesTraversed("NodesTraversed.txt");
    ofstream myNodesInsideBox("InsideBox.txt");
    std::chrono::high_resolution_clock::time_point startMakingForestWithThreads;
    std::chrono::high_resolution_clock::time_point endMakingForestWithThreads;
    std::chrono::high_resolution_clock::time_point startInsideBox; 
    std::chrono::high_resolution_clock::time_point endInsideBox;
    std::chrono::high_resolution_clock::time_point startCopyToDevice;
    std::chrono::high_resolution_clock::time_point endCopyToDevice;
    //for(int i = 0; i < 12; i++){
    int numberOfHits = 100;
    vector<int> number_nodes;
    while( numberOfHits <= 30000){
        numberOfHits+=100;
        // Generate points:
        generateRandomPointCloud(cloudn, cloud, numberOfHits);
        number_of_nodes_traversed = 0;
        
        //must be defined {1, 2, 3} = {x, y, z}
        vector<int> dimensions = {1,2,3};
        int number_of_dimensions = dimensions.size();
        
        //get_size_of_tree from cuda_device --> #datapoints per thread.. = datapoints per tree
        /*int device;
        cudaGetDevice(&device);
        std::cout<< "devices are: " << device << endl;
        
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, device);
        printDevProp(devProp);
        int warp_size = devProp.warpSize;
        int max_threads_per_block = devProp.maxThreadsPerBlock;*/
        //TODO: here calculate #nodes need
        
        
        //cup-version
        int warp_size = 1;
        
        
        //???: will it be one tree per warp or per block?
        //formula: 2^(floor(log_2(64)+1))-1 is always max when one less than warp_size
        //int datapoints_per_tree = warp_size-1;
        int datapoints_per_tree = (numberOfHits+warp_size-1)/warp_size;
        int treeSize = pow(2, floor(log2(datapoints_per_tree)+1))-1;
        datapoints_per_tree = treeSize;
        //round up: q = (x + y - 1) / y;
        int threads = (numberOfHits+datapoints_per_tree-1)/datapoints_per_tree;
        //threads = warp_size-1;
        
        //for cpu:
        threads = 1;
        cout << "number of warps " << warp_size << endl;
        cout << "datapoints_per_tree: " << datapoints_per_tree << endl;
        cout << "treeSize " << treeSize << endl;
        cout << "number of threads e.g. number of trees " << threads << endl;
        
        //all not-used elements in treesArray are by default set to zero by compiler
        int *treesArray_ID = new int[threads*datapoints_per_tree]();
        //initializing results with 0
        int *treeResults_ID = new int[threads*datapoints_per_tree]();
        num_t *treesArray;
        treesArray = new num_t [threads*datapoints_per_tree*number_of_dimensions]();//[number_of_dimensions+1][threads*datapoints_per_tree];
        startMakingForestWithThreads = std::chrono::high_resolution_clock::now();
        make_forest<num_t>(cloudn, cloud, dimensions, datapoints_per_tree, threads, treesArray, treesArray_ID);
        endMakingForestWithThreads = std::chrono::high_resolution_clock::now();
        
        
        //test if trees made with make_forest are correct:
        vector<vector<Point<num_t>>> trees = test_correct_trees(treesArray, treesArray_ID, datapoints_per_tree, threads, dimensions, numberOfHits, cloudn);
        
        
        /*
        //print treesArray
        int c=0;
        for(int i=0; i<datapoints_per_tree*threads*number_of_dimensions-1; i+=number_of_dimensions){
            cout << treesArray[i+0] << " " << treesArray[i+1] << " " << treesArray[i+2] << " ID:" << treesArray_ID[c] << endl;
            c++;
            
        }
        for(int i=0; i<threads; i++){
            cout << " NEW TREE" << endl;
            
            for( int j=0; j< datapoints_per_tree; j++){
                cout << treesArray[i*datapoints_per_tree*number_of_dimensions+j*number_of_dimensions+0] << " " << treesArray[i*datapoints_per_tree*number_of_dimensions+j*number_of_dimensions+1] << " " << treesArray[i*datapoints_per_tree*number_of_dimensions+j*number_of_dimensions+2] << " ID:" << treesArray_ID[i*datapoints_per_tree+j] << endl;
                cout << trees[i][j].x << " " << trees[i][j].y << " " << trees[i][j].z << " ID:" << trees[i][j].ID <<  endl;
            }
        }
      */
    
        //make box, in which should be searched for hits
        //set all other dimensions to zero, if not used:
        int box[6] = {8, 8, 1, 1, 2, 2};
        
        //testing on CPU
        std::vector<int> VtreesArray(treesArray, treesArray + threads*datapoints_per_tree*number_of_dimensions);
        std::vector<int> VtreesArray_ID(treesArray_ID, treesArray_ID + threads*datapoints_per_tree);
        std::vector<int> VtreeResults_ID(treeResults_ID, treeResults_ID + threads*datapoints_per_tree);
        std::vector<int> Vbox(box, box + 6);
        startInsideBox = std::chrono::high_resolution_clock::now();
        insideBoxCPU(VtreesArray, VtreesArray_ID, VtreeResults_ID, Vbox, treeSize, number_of_dimensions);
        endInsideBox = std::chrono::high_resolution_clock::now();
        
        number_nodes.push_back(number_of_nodes_traversed);
        for(int i=0; i<number_nodes.size();i++){
            cout  << number_nodes[i] << ",";
        }
        
        //get unique ID's
        auto last = std::unique(VtreeResults_ID.begin(), VtreeResults_ID.end());
        VtreeResults_ID.erase(last, VtreeResults_ID.end());
        
        /*
        Cuda_class<num_t> tree;
    
        startCopyToDevice = std::chrono::high_resolution_clock::now();
        tree.cudaCopyToDevice(threads, datapoints_per_tree, treesArray, treesArray_ID, box, number_of_dimensions);
        endCopyToDevice = std::chrono::high_resolution_clock::now();
    
        startInsideBox = std::chrono::high_resolution_clock::now();
        tree.cudaInsideBox(threads, datapoints_per_tree, number_of_dimensions, treesArray, treesArray_ID, box);
        cudaDeviceSynchronize();
        endInsideBox = std::chrono::high_resolution_clock::now();
    
        tree.cudaCopyToHost(treesArray_ID);
         */
         
        
        //TESTING CORRECTNESS ------------------------------------------------------------------------------------------------------
        //test wether the resulting ID's are correct, and the only ones inside box:
        /*vector<int> dummyResult = inBox(threads, datapoints_per_tree, box,trees);
        
        for( int i = 0; i<threads*datapoints_per_tree; i++){
            if(VtreeResults_ID[i] != 0){
                std::cout << "ID real: " << VtreeResults_ID[i]<< " ID dummy: " << dummyResult[i]<<std::endl;
            }
            if(VtreeResults_ID[i] != dummyResult[i]){
                //std::cout << "ID real: " << VtreeResults_ID[i]<< " ID dummy: " << dummyResult[i]<<std::endl;
                if(VtreeResults_ID[i] == 0 && dummyResult[i] == -1) continue;
                    std::cout << "NOT same!!! ID real: " << VtreeResults_ID[i]<< " ID dummy: " << dummyResult[i]<< " cloudsize " << numberOfHits << std::endl;
            }
        }
        
        for(int i = 0; i< threads*datapoints_per_tree; i++){
            if(VtreeResults_ID[i] != 0){
            std::cout << "ID: " << VtreeResults_ID[i]<< std::endl;
            }
        }*/
        
        cloud.clear();
        
        myThreadFile << std::chrono::duration_cast<std::chrono::microseconds>(endMakingForestWithThreads-startMakingForestWithThreads).count() << ",";
        myCudaFile  << std::chrono::duration_cast<std::chrono::microseconds>(endInsideBox-startInsideBox).count() << ",";
        myTreeSizeFile << to_string(treeSize) << ",";
        myNodesTraversed << to_string(number_of_nodes_traversed) << ",";
        myNodesInsideBox << to_string(VtreeResults_ID.size()) << ",";
    }
    myThreadFile.close();
    myCudaFile.close();
    myTreeSizeFile.close();
    myNodesTraversed.close();
    
    return 0;
}



