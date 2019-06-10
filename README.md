# Introduction #
This project implements an Ant-based algorithm to find a maximum clique from a given graph. It can be used with a compiler (e.g. GCC or Intel Compiler suite) that supports OpenMP for parallelism.

## Requires ##
* C++ Compiler

## Compile and Run ###
```
$ cd src
$ make

# run the program ../examples/brock200_1.clq.b with seed 1
$ ./maxclique  ../examples/brock200_1.clq.b 1
nthreads: 1 cliq_siz: 19 vertices: 200 edges: 14834 clique: (58,7,102,57,31,143,121,61,126,103,187,21,56,37,45,9,124,106,42) seed: 1
```

You might also be interested in

* **Converter**: convert from Dimacs binary to ascii and vice versa
* **Benchmark graphs**: http://www.cs.hbg.psu.edu/benchmarks/clique.html
* **Papers**: see these papers in my [publication](http://cse.unl.edu/~tnguyen/) section

   * Bui, T., T. Nguyen, and J. Rizzo, “Parallel Shared Memory Strategies for Ant-Based Optimization Algorithms”, Genetic and Evolutionary Computation Conference (GECCO), 2009
