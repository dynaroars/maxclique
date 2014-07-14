# README #
This project implements an Ant-based algorithm to find a maximum clique from a given graph. It can be used with a compiler (e.g. GCC or Intel Compiler suite) that supports Open MP for parallelism.

Additional information can be found from these papers: 
* Bui, T., T. Nguyen, and J. Rizzo, “Parallel Shared Memory Strategies for Ant-Based Optimization Algorithms”, Genetic and Evolutionary Computation Conference (GECCO), 2009


## Setup ##
The source code of this program is released under the BSD license and can be downloaded using the command "hg clone https://nguyenthanhvuh@bitbucket.org/nguyenthanhvuh/maxclique/"

This code is written in C++ and can be compiled using GCC, e.g., gcc main.cc .
The program has been tested using the following setup:
* Debian Linux
* GNU C++ Compiler


Other related information
* **Converter**: convert from Dimacs binary to ascii and vice versa
* **Benchmark graphs**: http://www.cs.hbg.psu.edu/benchmarks/clique.html
* **Papers**: see these papers in my [publication](http://cs.unm.edu/~tnguyen/) section

