# Neuromechanical model of forward and backward locomotion in C. elegans 

This repository contains all the classes necessary to evolve a neuromechanical model of C. elegans to produce forward and backward locomotion. 

We are using this to better understand how the same neural circuit can produce multiple behaviors (i.e., forward and backward) and moulate between them through command neurons.

Work in collaboration with Dr. Erick Olivares and Prof. Randall Beer.

## Instructions for use

1. Compile using the Makefile: 
```
$ make
```
2. Perform an evolutionary run (takes aprox. 2 minutes): 
```
$ ./main
```
3. Visualize the evolutionary progress: 
```
$ python viz.py
```
4. Visualize the movement of the body using Mathematica: 
```
$ Mathematica viz.nb
```
Optional. You can also visualize the dynamic neural and muscle activity by plotting "muscles.dat" and "neural.dat".
