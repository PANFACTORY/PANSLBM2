# PANSLBM2

A library of Topology Optimization with LBM.

# Usage

## Setup

Clone source codes into your environment.  
Although your environment is allowed Windows, Linux or MacOS, use Windows for example bellow.

```
git clone https://github.com/PANFACTORY/PANSLBM2.git && cd PANSLBM2
```

After that, working directory is ```PANSLBM2/```.  

## 2D heatsink

For example, run a program solving 2D steady-state heatsink problem.  
Use MinGW with AVX2 and OpenMP.  

```
g++ -mavx -fopenmp -o build/heatsink.exe production/heatsink.cpp
```
```
build/heatsink.exe
```
And results are output in ```result/```.  
Of course, you can also use other compilers like Intel C++ Compiler.

## 3D heatsink

For example, run a program solving 3D steady-state heatsink problem.  
Use MinGW with AVX2, OpenMP and MicrosoftMPI.   

You have to place two files in your working directory to use MicrosoftMPI.  
- libmsmpi.a
- msmpi.def

And add or uncomment ```#define _USE_MPI_DEFINES``` at 1st line in ```production/heatsink3D.cpp```.  
After that, build and execute with these commands.
```
g++ -mavx -fopenmp -lmsmpi -L . -o build/heatsink3D.exe production/heatsink3D.cpp
```
```
mpiexec -n 4 build/heatsink3D.exe 2 2 1
```
And results are output in ```result/```.  
Of course, you can also use other compilers like Intel C++ Compiler and other MPI libraries like Intel MPI. 

# Dependency

- AVX2
- OpenMP
- MPI

# License

[MIT](./LICENSE)