(1)compile
mpic++  -I ./ -I include/ calculategk.cpp platform.cpp simulator.cpp -fopenmp -g -o plt-mpi

(2)running
mpirun  -n x -host node1:1,node2:1,...noden:1 ./plt-mpi network sample

(3)Multithreaded setup
#define NUM_THREADS n (n is the thread number) in file platform.cpp

(4)test sample
refer to gendataset_guass.cpp

initialization, Self-defined functionality or Replacing the simulator please refer to LFQAP https://github.com/xzphys/LFQAP


The other github projects integrated in this project:
1. For full amplitude simulator, our platform integrated Quantum-Computing-Library, but do some refactor, trimmed some unused functionality. The github website of Quantum-Computing-Library is https://github.com/AbeerVaishnav13/Quantum-Computing-Library.
2. For tensor network simulator, our platform integrated qtorch, the source code in folder. The github website of qtorch is https://github.com/aspuru-guzik-group/qtorch.
