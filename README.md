This repository is for benchmarking ParSy vs. other libraries i.e., MKL Pardiso, Pardiso, and PaSTiX. 

Since ParSy's code is fixed for any matrix, the code is generated once and is used in this repository for testing. 

## Installation
### Library requirements
Suitesparse, MKL Pardiso, METIS, and optionally SCOTCH and PaSTiX libraries need to be installed and their corresponding variables 
need to be set in the CMakeLists.txt file in the root directory.


### Building the project

The first step is to set the environmental variables corresponding to each library. The following shows how the variables are set in bash.
```bash
export MKLROOT <path to MKL>
export SUITEROOT <path to Suitesparse>
export METISROOT <path to METIS> 
```
The environment variable can be also set manually in the cmake file.
The second step, the framework should be built using the following commands:
You should use cmake to build the :
```bash
cd where/you/cloned/parsy_bench
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make 
```

### Evaluation
For evaluating ParSy Code, you first need to download the matrices with matrix market format from Suitesparse matrix repository or any other sources. For symmetric matrices, only lower half of the matrix has to be stored i.e., similar to Suitesparse matrix repository. A program is created to generate lower half of the matrix if it is stored in full. 

After downloading the matrices, the following commands evaluate ParSy for all matrices.
```bash
cd where/you/cloned/parsy_bench/scripts
./eval ../build/examples/choleskyTest  path/to/matrices/folder
```
And for evaluating libraries the following commands can be used:
```bash
cd where/you/cloned/parsy_bench/scripts
./lib_eval ../build/libExample/MKLCholesky  path/to/matrices/folder
```