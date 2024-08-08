Developers:
Xin Zhang1, Xiaoyu Li2, Lianfu Wei3, Qinsheng Zhu4, Geng Chen4, Wenjie Shun5, Lianhui Yu6, Yuexian Hou1,*
1 (College of Intelligence and Computing, Tianjin University, Tianjin 300350, China)
2 (School of Information and Software Engineering, University of Electronic Science and Technology of China, Chengdu 610054)
3 (College of Information Science and Technology, Southwest Jiaotong University, Chengdu 610031)
4 (School of Computer Science and Engineering, University of Electronic Science and Technology of China, Chengdu 610000)
5 (School of Electronic Science and Engineering, University of Electronic Science and Technology of China, Chengdu 610000)
6 (School of Physics, University of Electronic Science and Technology of China, Chengdu 610000)
*yxhou@tju.edu.cn

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
