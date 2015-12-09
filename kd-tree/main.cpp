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
#include "cuda_runtime.h"
#include "cuda.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <future>
#include "InsideBox.hpp"
#include <fstream>
#include <numeric>
#include <iomanip>
#include <random>
#include <math.h>
#define MYDEVICE 0

using namespace std;


template <typename num_t>
void generateRandomPointCloud(int datapoints_per_tree, vector<Point<num_t>> &pointn, vector<Hit<num_t>> &point, const size_t N, const int max_range = 10)
{
    //TODO: change, s.t. each tree has hit !
    //jetzt hat ein tree-forest genau ein hit inside box. Realistisch? nein! jeder tree sollte ein hit inside box haben!
    //in jeder threadcloud ein hit inside box
    //d.h. alle anzahl_hits_per_tree + 30 muss box-hit sein
    cout << "Generating "<< N << " point cloud...\n";
    point.resize(N);
    pointn.resize(N);
    
    const int range_from_x  = 0;
    const int range_to_x    = 10;
    const int range_from_y  = -50;
    const int range_to_y    = 50;
    const int range_from_z  = 0;
    const int range_to_z    = 2*M_PI;
    
    std::random_device rand_dev;
    
    for (size_t i=0;i<N;i++)
    {
        
        std::mt19937 generator(rand_dev());
        std::uniform_real_distribution<double>  distr_x(range_from_x, range_to_x);
        std::uniform_real_distribution<double>  distr_y(range_from_y, range_to_y);
        std::uniform_real_distribution<double>  distr_z(range_from_z, range_to_z);
        
        pointn[i].x = distr_x(generator);
        pointn[i].y = distr_y(generator);
        pointn[i].z = distr_z(generator);
        pointn[i].ID = i+1;
        
        
        
        
        
        //not needed anymore, since not one hit per tree should be inside box, but more
        /*
         int a = (int)pointn[i].x;
         int b = (int)pointn[i].y;
         int c = (int)pointn[i].z;
         
         while(a == specialPoint_x && b == specialPoint_y && c == specialPoint_z){
         //cout << " 15 14 5 ID:" << i+1 << endl;
         pointn[i].x = max_range * (rand() % 1000) / num_t(500);
         pointn[i].y = max_range * (rand() % 1000) / num_t(500);
         pointn[i].z = max_range * (rand() % 1000) / num_t(500);
         
         a = pointn[i].x;
         b = pointn[i].y;
         c = pointn[i].z;
         }
         
         if(i%datapoints_per_tree == 30){
         pointn[i].x = specialPoint_x;
         pointn[i].y = specialPoint_y;
         pointn[i].z = specialPoint_z;
         }*/
        
        
        
        
        //vector<num_t> p = {max_range * (rand() % 1000) / num_t(1000), max_range * (rand() % 1000) / num_t(1000), max_range * (rand() % 1000) / num_t(1000)};
        
        vector<num_t> p = {pointn[i].x, pointn[i].y,pointn[i].z};
        point[i].datapoints.resize(3);
        point[i].datapoints = p;
        point[i].ID = i+1;
        
        
        //cout << point[i].datapoints[0] << ", " << point[i].datapoints[1] << ", " << point[i].datapoints[2] << "\t ID: " << point[i].ID << endl;
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
vector<int> inBox(int threads, int datapoints_per_tree, num_t box[], vector<vector<Point<num_t>>> &trees ){
    vector<int> result;
    for(int i=0; i<threads; i++){
        for(int j=0; j< datapoints_per_tree; j++){
            if( box[0] <= trees[i][j].x && box[1]>=trees[i][j].x && box[2] <=trees[i][j].y && box[3] >= trees[i][j].y && box[4] <= trees[i][j].z && box[5] >= trees[i][j].z){
                result.push_back(trees[i][j].ID);
            }
            else{ result.push_back(-1);}
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

int number_of_nodes_traversed = 0;
void CPUtraverseTreeRecursiveIF( vector<int> &treeArray_values, vector<int> &treeArray_ID, vector<int> &results, vector<int> &box, int pos, int startOfTree, int endOfTree, int number_of_dimensions){
    
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
                    CPUtraverseTreeRecursiveIF(treeArray_values, treeArray_ID,results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                    
                    //right child:
                    pos += 1;
                    CPUtraverseTreeRecursiveIF(treeArray_values, treeArray_ID,results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
                
            }
            //if sorted dimension is larger than box follow branch of smaller child = left child
            else if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] > box[2*level_of_dimension+1]){
                //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " LINKES KIND" << endl;
                //cout << "number of nodes traversed " << number_of_nodes_traversed << endl;
                //left child:
                if(level != lastLevel){
                    pos *= 2;
                    CPUtraverseTreeRecursiveIF(treeArray_values, treeArray_ID,results, box, pos, startOfTree, endOfTree, number_of_dimensions);
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
                    CPUtraverseTreeRecursiveIF(treeArray_values, treeArray_ID, results, box, pos, startOfTree, endOfTree, number_of_dimensions);
                }
            }
        }
    }
}
void CPUtraverseTreeIterative(vector<int> &treeArray_values, vector<int> &treeArray_ID, vector<int> &results, vector<int> &box, vector<int> &queue,int pos, int startOfTree, int endOfTree, int number_of_dimensions){
    
    int lastLevel = ceil(log2(double(endOfTree+1))-1);
    
    
    //for GPU queue will be array (for debugging reasons vector here)
    
    queue[0]= pos;
    int queueFront = 0;
    int queueRear = 1;
    int queueSize = 1;
    int numberOfMightHits = 0;
    number_of_nodes_traversed = 0;
    
    
    while(queueSize != 0){
        number_of_nodes_traversed++;
        queueSize--;
        pos = queue[queueFront++];
        
        int level = ceil(log2(double(pos+1))-1);
        int level_of_dimension = level%number_of_dimensions;
        
        //cout << "position " << pos-1 << endl;
        
        //if sorted dimension inside box continue with both branches
        if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] >= box[2*level_of_dimension] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] <= box[2*level_of_dimension+1]){
            //put left and right child in queue:
            //cout << "ID: " << treeArray_ID[startOfTree+pos-1] <<" bei pos: " << startOfTree+pos-1 << " BEIDE BRANCHES" << endl;
            
            //nested ifs - bad!
            if(level != lastLevel){
                queue[queueRear++] = pos*2;
                queue[queueRear++] = pos*2+1;
                queueSize+=2;
            }
            
            //and check if node is totally inside box - then right it to results
            //possibility(?) can this be checked after tree traversed? Can an extra thread check this? 2 threads per tree?
            //write pos to array? What about size of array on GPU? Size of tree must be allocated? Is there so much space?
            //per warp needed memory: 3*TreeSize - NO !!! can be put into results :)
            //toCheckIfInside.push_back(pos);
            //these results have to be checked again!!!
            results[startOfTree+numberOfMightHits] = pos; //treeArray_ID[startOfTree+pos-1];
            numberOfMightHits++;
            
        }
        //else if sorted dimensions > inside box continue with left child
        else if(treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+level_of_dimension] > box[2*level_of_dimension+1]
                && level != lastLevel){
            //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " LINKES KIND" << endl;
            queue[queueRear++] = pos*2;
            queueSize++;
        }
        //else sorted dimension < inside box continue with right child
        else if(level != lastLevel){
            //cout << "ID: " << treeArray_ID[startOfTree+pos-1] << " bei pos: " << startOfTree+pos-1 << " RECHTES KIND" << endl;
            queue[queueRear++] = pos*2+1;
            queueSize++;
        }
    }
    
    //check nodes, that might be inside box:
    for(int j = 0; j<=numberOfMightHits;j++){
        pos = results[startOfTree+j];
        results[startOfTree+j] = 0;
        
        bool inside = true;
        for(int i=0; i<number_of_dimensions; i++){
            if( treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] >= box[2*i] && treeArray_values[startOfTree*number_of_dimensions+number_of_dimensions*(pos-1)+i] <= box[2*i+1] ){
                //entirely inside box for all dimensions
            }
            else{
                inside = false;
            }
        }
        if(inside){
            results[startOfTree+pos-1] = treeArray_ID[startOfTree+pos-1];
        }
    }
}


void insideBoxCPU(vector<int> &treeArray_values, vector<int> &treeArray_ID, vector<int> &results, vector<int> &box,vector<int> &queue, int tree_size, int number_of_dimensions){
    
    int threadIdx = 0;
    int startOfTree = threadIdx* tree_size;
    int endOfTree = startOfTree + (tree_size - 1);
    
    CPUtraverseTreeRecursiveIF(treeArray_values, treeArray_ID,results, box, 1, startOfTree, endOfTree, number_of_dimensions);
    //CPUtraverseTreeIterative(treeArray_values, treeArray_ID,results, box, queue, 1, startOfTree, endOfTree, number_of_dimensions);
    //CPUtraverseTreeBFSQueue(treeArray_values, treeArray_ID,results, box, 1, startOfTree, endOfTree, number_of_dimensions);
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
    typedef double num_t;
    
    vector<Hit<num_t>> InnerCloud;
    vector<Hit<num_t>> OuterCloud;
    vector<Point<num_t>> InnerCloud_simple;
    vector<Point<num_t>> OuterCloud_simple;
    
    fstream myOutputFile("./OutputIterative.txt", std::fstream::in | std::fstream::out | std::fstream::app);
    ofstream myThreadFile ("ThreadingTimesIterative.txt");
    ofstream myCudaFile("CudaTimesIterative.txt");
    ofstream myCudaStdev("CudaStdDevIterative.txt");
    ofstream myTreeSizeFile("TreeSizesIterative.txt");
    ofstream myNodesTraversed("NodesTraversedIterative.txt");
    ofstream myNodesStdev("NodesStdDevIterative.txt");
    ofstream myNodesInsideBox("InsideBoxIterative.txt");
    /*
     fstream myOutputFile("./OutputRecursive.txt", std::fstream::in | std::fstream::out | std::fstream::app);
     ofstream myThreadFile ("ThreadingTimes.txt");
     ofstream myCudaFile("CudaTimes.txt");
     ofstream myCudaStdev("CudaStdDev.txt");
     ofstream myTreeSizeFile("TreeSizes.txt");
     ofstream myNodesTraversed("NodesTraversed.txt");
     ofstream myNodesStdev("NodesStdDev.txt");
     ofstream myNodesInsideBox("InsideBox.txt");
     */
    std::chrono::high_resolution_clock::time_point startMakingForestWithThreads;
    std::chrono::high_resolution_clock::time_point endMakingForestWithThreads;
    std::chrono::high_resolution_clock::time_point startInsideBox;
    std::chrono::high_resolution_clock::time_point endInsideBox;
    std::chrono::high_resolution_clock::time_point startCopyToDevice;
    std::chrono::high_resolution_clock::time_point endCopyToDevice;
    //for(int i = 0; i < 12; i++){
    int numberOfHits = 0;
    vector<int> number_nodes;
    
    //while( numberOfHits <= 50000){
    numberOfHits+=200;
    int treeSize;
    std::vector<int> VtreeArray_results;
    
    //for deviation and median
    std::vector<double> CudaTimes;
    std::vector<int> NodesTraversed;
    for(int p = 0; p<=10;p++){
        number_of_nodes_traversed = 0;
        
        //must be defined {1, 2, 3} = {x, y, z}
        //set which dimensions
        vector<int> dimensions = {1,2,3};
        int number_of_dimensions = dimensions.size();
        
        //get_size_of_tree from cuda_device --> #datapoints per thread.. = datapoints per tree
        //get GPU device properties
        int device;
        cudaGetDevice(&device);
        std::cout<< "devices are: " << device << endl;
        
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, device);
        printDevProp(devProp);
        int warp_size = devProp.warpSize;
        int max_threads_per_block = devProp.maxThreadsPerBlock;
        //TODO: here calculate #nodes need
        
        
        //cup-version
        //int warp_size = 1;
        
        
        //???: will it be one tree per warp or per block?
        //formula: 2^(floor(log_2(64)+1))-1 is always max when one less than warp_size
        //int datapoints_per_tree = warp_size-1;
        int datapoints_per_tree = (numberOfHits+warp_size-1)/warp_size;
        int treeSize = pow(2, floor(log2(datapoints_per_tree)+1))-1;
        datapoints_per_tree = treeSize;
        //round up: q = (x + y - 1) / y;
        int threads = (numberOfHits+datapoints_per_tree-1)/datapoints_per_tree;
        //threads = warp_size-1;
        
        // Generate points:
        /*const int max_range = 10;
         int specialPoint_x = max_range * (rand() % 1000) / num_t(500);
         int specialPoint_y = max_range * (rand() % 1000) / num_t(500);
         int specialPoint_z = max_range * (rand() % 1000) / num_t(500);
         */
        generateRandomPointCloud(datapoints_per_tree, OuterCloud_simple, OuterCloud, numberOfHits);
        generateRandomPointCloud(datapoints_per_tree, InnerCloud_simple, InnerCloud, numberOfHits);
        
        
        //for cpu:
        //threads = 1;
        cout << "number of warps " << warp_size << endl;
        cout << "datapoints_per_tree: " << datapoints_per_tree << endl;
        cout << "treeSize " << treeSize << endl;
        cout << "number of threads e.g. number of trees " << threads << endl;
        
        //all not-used elements in treesArray are by default set to zero by compiler
        int *treesArray_ID = new int[threads*datapoints_per_tree]();
        //initializing results with 0
        int *treeArray_results = new int[threads*datapoints_per_tree]();
        int *queue = new int[threads*datapoints_per_tree]();
        num_t *treesArray;
        treesArray = new num_t [threads*datapoints_per_tree*number_of_dimensions]();//[number_of_dimensions+1][threads*datapoints_per_tree];
        startMakingForestWithThreads = std::chrono::high_resolution_clock::now();
        make_forest<num_t>(InnerCloud_simple, InnerCloud, dimensions, datapoints_per_tree, threads, treesArray, treesArray_ID);
        endMakingForestWithThreads = std::chrono::high_resolution_clock::now();
        
        
        //test if trees made with make_forest are correct:
        vector<vector<Point<num_t>>> trees = test_correct_trees(treesArray, treesArray_ID, datapoints_per_tree, threads, dimensions, numberOfHits, InnerCloud_simple);
        
        
        //make box, in which should be searched for hits
        //set all other dimensions to zero, if not used:
        //make each box, s.t. it has one hit inside tree NOT inside forest!
        //int box[warp_size*2*number_of_dimensions];
        num_t box[numberOfHits*6];
        
        for(int i = 0; i<numberOfHits; i+=6 ){
            
            box[i+0] = (OuterCloud[i].datapoints[0]-1 < 0)      ? 10-(OuterCloud[i].datapoints[0]-1)        : OuterCloud[i].datapoints[0]-1;
            box[i+1] = (OuterCloud[i].datapoints[0]+1 > 10)     ? OuterCloud[i].datapoints[0]+1-10          : OuterCloud[i].datapoints[0]+1;
            box[i+2] = (OuterCloud[i].datapoints[1]+10 > 50)    ? OuterCloud[i].datapoints[1]+10-100        : OuterCloud[i].datapoints[1]+10;
            box[i+3] = (OuterCloud[i].datapoints[1]-10 < -50)   ? 50-(OuterCloud[i].datapoints[i]-10)       : OuterCloud[i].datapoints[i]-10;
            box[i+4] = (OuterCloud[i].datapoints[2]+.1 > 2*M_PI)? (OuterCloud[i].datapoints[2]+.1)-2*M_PI   : OuterCloud[i].datapoints[2]+.1;
            box[i+5] = (OuterCloud[i].datapoints[2]-.1 < 0)     ? 2*M_PI-(OuterCloud[i].datapoints[2]-.1)   : OuterCloud[i].datapoints[2]-.1;
            
        }
        
        //testing on CPU
        /*std::vector<int> VtreesArray(treesArray, treesArray + threads*datapoints_per_tree*number_of_dimensions);
         std::vector<int> VtreesArray_ID(treesArray_ID, treesArray_ID + threads*datapoints_per_tree);
         std::vector<int> VtreeArray_results(treeArray_results, treeArray_results + threads*datapoints_per_tree);
         std::vector<int> Vbox(box, box + 6);*/
        
        
        vector<int> dummyResult = inBox(threads, datapoints_per_tree, box, trees);
        
        
        /*
         std::vector<int> Vqueue(treeSize);
         startInsideBox = std::chrono::high_resolution_clock::now();
         insideBoxCPU(VtreesArray, VtreesArray_ID, VtreeArray_results, Vbox, Vqueue, treeSize, number_of_dimensions);
         endInsideBox = std::chrono::high_resolution_clock::now();
         
         number_nodes.push_back(number_of_nodes_traversed);
         */
        
        /*for(int i=0; i<number_nodes.size();i++){
         cout  << number_nodes[i] << ",";
         }*/
        
        
        
        
        Cuda_class<num_t> tree;
        
        startCopyToDevice = std::chrono::high_resolution_clock::now();
        tree.cudaCopyToDevice(threads, datapoints_per_tree, treesArray, treesArray_ID, treeArray_results, box, queue, number_of_dimensions, numberOfHits);
        endCopyToDevice = std::chrono::high_resolution_clock::now();
        
        startInsideBox = std::chrono::high_resolution_clock::now();
        tree.cudaInsideBox(threads, datapoints_per_tree, number_of_dimensions, treesArray, treesArray_ID, treeArray_results, box, queue, numberOfHits);
        cudaDeviceSynchronize();
        endInsideBox = std::chrono::high_resolution_clock::now();
        
        tree.cudaCopyToHost(treeArray_results, numberOfHits);
        
        
        //TESTING CORRECTNESS ------------------------------------------------------------------------------------------------------
        //test wether the resulting ID's are correct, and the only ones inside box:
        //vector<int> dummyResult = inBox(threads, datapoints_per_tree, box,trees);
        /*bool resultsCorrect = true;
         for( int i = 0; i<threads*datapoints_per_tree; i++){
         if(treeArray_results[i] != 0){
         std::cout << "ID real: " << treeArray_results[i]<< " ID dummy: " << dummyResult[i]<<std::endl;
         cout << "Values: " << treesArray[i*number_of_dimensions+0] << " " << treesArray[i*number_of_dimensions+1] << " " << treesArray[i*number_of_dimensions+2] << endl;
         }
         if(treeArray_results[i] != dummyResult[i]){
         //std::cout << "ID real: " << VtreeArray_results[i]<< " ID dummy: " << dummyResult[i]<<std::endl;
         if(treeArray_results[i] == 0 && dummyResult[i] == -1) continue;
         std::cout << "NOT same!!! ID real: " << treeArray_results[i]<< " ID dummy: " << dummyResult[i]<< " cloudsize " << numberOfHits << std::endl;
         resultsCorrect = false;
         }
         }
         if(resultsCorrect) cout << "All inside Box are correct!" << endl;*/
        /*sort( treeArray_results.begin(), treeArray_results.end() );
         VtreeArray_results.erase( unique( treeArray_results.begin(), treeArray_results.end() ), treeArray_results.end());
         
         for(int i = 0; i< VtreeArray_results.size(); i++){
         std::cout << "FINAL ID: " << VtreeArray_results[i]<< std::endl;
         }*/
        
        OuterCloud.clear();
        InnerCloud.clear();
        OuterCloud_simple.clear();
        InnerCloud_simple.clear();
        
        CudaTimes.push_back(std::chrono::duration_cast<std::chrono::microseconds>(endInsideBox-startInsideBox).count());
        NodesTraversed.push_back(number_of_nodes_traversed);
        
        delete[] treeArray_results;
        delete[] queue;
        delete[] treesArray_ID;
        delete[] treesArray;
        
    }
    //getting mean and standard deviation
    double CudaSum = std::accumulate(CudaTimes.begin(), CudaTimes.end(), 0.0);
    double CudaMean = CudaSum / CudaTimes.size();
    double CudaSq_sum = std::inner_product(CudaTimes.begin(), CudaTimes.end(), CudaTimes.begin(), 0.0);
    double CudaStdev = std::sqrt(CudaSq_sum / CudaTimes.size() - CudaMean * CudaMean);
    
    double NodesSum = std::accumulate(NodesTraversed.begin(), NodesTraversed.end(), 0.0);
    double NodesMean = NodesSum / NodesTraversed.size();
    double NodesSq_sum = std::inner_product(NodesTraversed.begin(), NodesTraversed.end(), NodesTraversed.begin(), 0.0);
    double NodesStdev = std::sqrt(NodesSq_sum / NodesTraversed.size() - NodesMean * NodesMean);
    
    
    if(myOutputFile.is_open())
    {
        myOutputFile << CudaMean << ";" << CudaStdev << ";" << treeSize << ";" << numberOfHits << ";\n";
    }
    else
    {
        cout << "Error writing to ";
    }
    
    myThreadFile << std::chrono::duration_cast<std::chrono::microseconds>(endMakingForestWithThreads-startMakingForestWithThreads).count() << ",";
    myCudaFile  << CudaMean << ",";
    myCudaStdev << CudaStdev << ",";
    myTreeSizeFile << to_string(treeSize) << ",";
    myNodesTraversed << NodesMean << ",";
    myNodesStdev << NodesStdev << ",";
    
    //cout << "VtreeArray_results.size() " << VtreeArray_results.size() << endl;
    //myNodesInsideBox << to_string(treeArray_results.size()) << ",";
    //}
    myThreadFile.close();
    myCudaFile.close();
    myTreeSizeFile.close();
    myNodesTraversed.close();
    
    return 0;
}