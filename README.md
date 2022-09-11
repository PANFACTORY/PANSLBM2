# PANSLBM2

# Usage

For example, run a program solving 2D steady-state heatsink.  
Use MinGW with AVX2 and openmp.

```
git clone https://github.com/PANFACTORY/PANSLBM2.git
```
```
cd PANSLBM2
```
```
g++ -mavx -fopenmp -o <your output directory>/<your output file> production/heatsink.cpp
```
```
<your output directory>/<your output file>
```