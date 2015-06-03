all:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib
	g++ -std=gnu++11 -L/usr/local/cuda/lib64 -lcudart -I/usr/local/cuda/include main.cpp KD_tree.hpp simple_kd_tree.hpp sorter.hpp -lstdc++ -pthread -o test.out -L/usr/local/cuda/lib -lcudart -lcuda