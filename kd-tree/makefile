GENCODE_SM20	:= -arch=sm_20 
GENCODE_FLAGS	:= $(GENCODE_SM20)

LDFLAGS		:= -L/usr/local/cuda/lib64 -L/usr/local/cuda/lib -lcuda -lcudart
CCFLAGS		:= -pthread -std=gnu++11 -lstdc++

NVCCFLAGS	:= -rdc=true -c

INCLUDES	:= -I/usr/local/cuda/include 


NVCC := nvcc
GCC := g++


#g++ -std=gnu++11 -L/usr/local/cuda/lib64 -lcudart -I/usr/local/cuda/include InsideBox.o InsideBox_link.o main.cpp #KD_tree.hpp simple_kd_tree.hpp  sorter.hpp -lstdc++ -pthread -o test.out -L/usr/local/cuda/lib -lcudart -lcuda

all: build

build: KD_tree

InsideBox.o: InsideBox.cu
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) $(INCLUDES) -o $@ $<
	$(NVCC) -dlink  $(GENCODE_FLAGS)  -o InsideBox_link.o $@

KD_tree.o: main.cpp KD_tree.hpp simple_kd_tree.hpp sorter.hpp
	$(GCC) $(CCFLAGS) $(INCLUDES) -o $@ -c $<

KD_tree: InsideBox.o KD_tree.o InsideBox_link.o
	$(GCC) $(CCFLAGS) -o $@ $+ $(LDFLAGS) $(EXTRA_LDFLAGS)

run: build
	./KD_tree

#all:
#	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib
#	g++ -std=gnu++11 -L/usr/local/cuda/lib64 -lcudart -I/usr/local/cuda/include main.cpp KD_tree.hpp simple_kd_tree.hpp sorter.hpp -lstdc++ -pthread -o test.out -L/usr/local/cuda/lib -lcudart -lcuda


clean:
	rm -rf *.o