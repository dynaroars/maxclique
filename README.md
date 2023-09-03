## Introduction

**AMC** implements a parallel Ant-based algorithm to find a maximum clique from a given graph. It can be used with a compiler (e.g. GCC) that supports OpenMP for parallelism.

## Requires
* C++ Compiler

## Compile and Run
```
$ cd src
$ make

# run the program ../examples/brock200_1.clq.b with seed 1
$ ./amc  ../examples/brock200_1.clq.b 1
nthreads: 1 cliq_siz: 19 vertices: 200 edges: 14834 clique: (58,7,102,57,31,143,121,61,126,103,187,21,56,37,45,9,124,106,42) seed: 1
```

## Miscs
You might also be interested in

* **Benchmark graphs**: [https://github.com/dynaroars/npbench/](https://github.com/dynaroars/npbench/)
* **[Converter](https://github.com/dynaroars/npbench/tree/master/instances/converter)**: convert from DIMACS binary to ASCII and vice versa

## Publications
   * Bui, T., T. Nguyen, and J. Rizzo, “Parallel Shared Memory Strategies for Ant-Based Optimization Algorithms”, Genetic and Evolutionary Computation Conference (GECCO), 2009
