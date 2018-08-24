//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_UTIL_H
#define CHOLOPENMP_UTIL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include "mkl.h"
#include "def.h"

/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
bool readMatrix_old(std::string fName, int &n, int &NNZ, int* &col,
                    int* &row, double* &val){
 /*This function reads the input matrix from "fName" file and
  * allocate memory for matrix A, L and U.
  * - The input file is a coordinate version and e
  * ach row of the file shows (col, row, nnz)
  * - The matrices are zero-indexed
  */
 std::ifstream inFile;
 inFile.open(fName);
 inFile >> n;
 inFile >> n;
 inFile>>NNZ;
 int factorSize= (n * n) / 2;//Worst case assumption
 if(n <= 0 || NNZ <= 0)
  return false;
 col = new int[n + 1]();
 // colL = new int[n + 1]; colU = new int[n + 1];
 row = new int[NNZ];
 // rowL = new int[factorSize]; rowU = new int[factorSize];
 val = new double[NNZ];
 // valL = new double[factorSize]; valU = new double[factorSize];
 if(!val || !col || !row)
  return false;
 //Initializing the result vector
 int y, x, colCnt=0, nnzCnt=0;
 double value;

 col[0]=0;
 for (int i = 0; nnzCnt<NNZ; ) {//Reading from file row by row
  inFile>>x;x--;
  inFile>>y;y--;//zero indexing
  inFile>>value;
  if(y > n)
   return false;
  if(y==i){
   val[nnzCnt]=value;
   row[nnzCnt]=x;
   colCnt++; nnzCnt++;
  }
  else{//New col
   col[i+1]=col[i]+colCnt;
   i++;//next iteration
   colCnt=1;
   val[nnzCnt]=value;
   row[nnzCnt]=x;
   nnzCnt++;
  }

 }
 col[n]= col[n - 1] + colCnt;//last col

 return true;
}

/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
bool readMatrix(std::string fName, size_t &n, size_t &NNZ, int* &col,//FIXME change col type to size_t
                int* &row, double* &val){
 /*This function reads the input matrix from "fName" file and
  * allocate memory for matrix A, L and U.
  * - The input file is a coordinate version and e
  * ach row of the file shows (col, row, nnz)
  * - The matrices are zero-indexed
  */

 std::ifstream inFile;
 inFile.open(fName);
 std::string line,banner, mtx, crd, arith, sym;
 /*  File format:
  *    %%MatrixMarket matrix coordinate real general/symmetric/...
  *    % ...
  *    % (optional comments)
  *    % ...
  *    #rows    #non-zero
  *    Triplet in the rest of lines: row    col    value
  */
 std::getline(inFile,line);
 for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
 std::istringstream iss(line);
 if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
  std::cout<<"Invalid header (first line does not contain 5 tokens)\n";
  return false;
 }

 if(banner.compare("%%matrixmarket")) {
  std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")\n";
  return false;
 }
 if(mtx.compare("matrix")) {
  std::cout<<"Not a matrix; this driver cannot handle that.\"\n";
  return false;
 }
 if(crd.compare("coordinate")) {
  std::cout<<"Not in coordinate format; this driver cannot handle that.\"\n";
  return false;
 }
 if(arith.compare("real")) {
  if(!arith.compare("complex")) {
   std::cout<<"Complex matrix; use zreadMM instead!\n";
   return false;
  }
  else if(!arith.compare("pattern")) {
   std::cout<<"Pattern matrix; values are needed!\n";
   return false;
  }
  else {
   std::cout<<"Unknown arithmetic\n";
   return false;
  }
 }
 while (!line.compare(0,1,"%"))
 {
  std::getline(inFile, line);
 }
 std::istringstream issDim(line);
 if (!(issDim >> n >> n >> NNZ)){
  std::cout<<"The matrix dimension is missing\n";
  return false;
 }
 if(n <= 0 || NNZ <= 0)
  return false;
 col = new int[n + 1]();
 // colL = new int[n + 1]; colU = new int[n + 1];
 row = new int[NNZ];
 // rowL = new int[factorSize]; rowU = new int[factorSize];
 val = new double[NNZ];
 // valL = new double[factorSize]; valU = new double[factorSize];
 if(!val || !col || !row)
  return false;
 //Initializing the result vector
 int y, x, colCnt=0, nnzCnt=0;
 double value;

 col[0]=0;
 for (int i = 0; nnzCnt<NNZ; ) {//Reading from file row by row
  inFile>>x;x--;
  inFile>>y;y--;//zero indexing
  inFile>>value;
  if(y > n)
   return false;
  if(y==i){
   val[nnzCnt]=value;
   row[nnzCnt]=x;
   colCnt++; nnzCnt++;
  }
  else{//New col
   col[i+1]=col[i]+colCnt;
   i++;//next iteration
   colCnt=1;
   val[nnzCnt]=value;
   row[nnzCnt]=x;
   nnzCnt++;
  }

 }
 col[n]= col[n - 1] + colCnt;//last col

 return true;
}




/*
 * reading a ordering from PasTiX stored format.
 */
bool readOrdering(std::string fName, size_t n, size_t* perm){


 std::ifstream inFile;
 inFile.open(fName);
 std::string line,banner, mtx, crd, arith, sym;
 size_t dontKnow;
 /*  File format:
  *    # 0
  *    #??    #rows
  *    values values ...
  */
 std::getline(inFile,line);
 std::istringstream iss(line);

 while (!line.compare(0,1,"%")){
  std::getline(inFile, line);
 }


 if (!(iss >> n )){
  std::cout<<"The matrix dimension is missing\n";
  return false;
 }

 if(n <= 0 )
  return false;
  //Initializing the result vector
 int y, x, colCnt=0, nnzCnt=0;
 double value;

 for (int i = 0; colCnt<n;i++,colCnt++) {//Reading from file row by row
  inFile>>x;
  perm[i]=x;
 }
 assert(colCnt==n);
}

bool enableColdCache(int n, std::ifstream &f){
 /*
  * n specifies the size of data for double computation. It depends
  * on the cache size
  */
 //TODO check file during read
 assert(!f.fail());
 double curVal;
 double **waste=new double*[n];
 for (int i = 0; i < n; ++i) {
  waste[i] = new double[n];
 }
 for (int i = 0; i < n; ++i) {
  for (int j = 0; j < n; ++j) {
   f>>curVal;
   waste[i][j]=curVal;
  }
 }
 for (int i = 0; i < n; ++i) {
  for (int j = 0; j < n; ++j) {
   for (int k = 0; k < n; ++k) {
    waste[i][j] += waste[i][k]*waste[k][j];
   }
  }
 }
 for (int i = 0; i < n; ++i) {
  delete waste[i];
 }
 delete waste;
 return true;
}

/*
 *
 */

void rhsInit(int n, MKL_INT *Ap, MKL_INT *Ai, double *Ax, double *b){
 /*generating a rhs that produces a result of all 1 vector*/
 for (int j = 0; j < n; ++j) {
  b[j]=0;
 }
 for (int c = 0; c < n ; ++c) {
  for (MKL_INT cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
   b[Ai[cc]]+=Ax[cc];
  }
 }
}

/*
 * RHS initilization for blocked
 */

void rhsInitBlocked(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, double *Ax, double *b){
 /*generating a rhs that produces a result of all 1 vector*/
 for (int j = 0; j < n; ++j) {
  b[j]=0;
 }
 for (int c = 0; c < n ; ++c) {
  for (int cc = Ap[c], j=0; cc < Ap[c + 1]; ++cc, ++j) {
   size_t curRow = Ai[AiP[c]+j];
   b[curRow]+=Ax[cc];
  }
 }
}

/*
 * Testing lower triangular solve to make sure all x elements are ONE.
 */

int testTriangular(size_t n, const double *x) {//Testing
 int test=0;
 for (int i = 0; i < n; ++i) {
  if(1-x[i]<0.001)
   test++;
  //else
  // cout<<i<<";";
 }
 if(n-test>0){
  return false;
 }
 return true;
}

/*
 * Converting BCSC to CSC
 */
int bcsc2csc(
  //in
  size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
  int *sup2col, double *Ax,
  //out
  MKL_INT *Cp, MKL_INT *Ci, double *Cx
){
 size_t actualNNZ=0;
 Cp[0] = 0;
 for (int i = 0; i < nBlocks; ++i) {
  int curCol = sup2col[i];
  int nxtCol = sup2col[i+1];
  int supWdt = nxtCol-curCol;
  assert(supWdt>0);

  for (int j = curCol; j < nxtCol; ++j) {
   for (MKL_INT k = Ap[j] + (j-curCol), kk=AiP[curCol] + (j-curCol);
        k<Ap[j+1]; ++k, ++kk) {
    Cx[actualNNZ] = Ax[k];
    assert(Ai[kk] < n);
    Ci[actualNNZ]= Ai[kk];
    actualNNZ++;
   }
   Cp[j+1] = actualNNZ;
  }
 }
 return 1;
}


void swapWSet(std::vector<std::vector<int>>& List, int i, int j){
 std::vector<int> tmp;
 if(List[i].size()==0 || List[j].size()==0){
  std::cout<<"sssss\n";
  return;
 }
 tmp.insert(tmp.begin(),List[i].begin(),List[i].end());
 assert(tmp.size() == List[i].size());
 List[i].erase(List[i].begin(),List[i].end());
 assert(List[i].size()==0);
 assert(tmp.size()>0);
 List[i].insert(List[i].begin(),List[j].begin(),List[j].end());

 List[j].erase(List[j].begin(),List[j].end());
 List[j].insert(List[j].begin(),tmp.begin(),tmp.end());
}


/*
 * Check whether the matrix is stored in full
 * return the lower half if it is full
 * Makes sense only for symmetric matrices.
 */
CSC * computeLowerTriangular(CSC *A){
 CSC *lowerHalf = new CSC;
 size_t lNNZcnt=0;
 for (int i = 0; i < A->ncol; ++i) {
  for (int j = A->p[i]; j < A->p[i+1]; ++j) {
   if(A->i[j]>=i)
    lNNZcnt++;
  }
 }
 //lower is already stored
 if(lNNZcnt==A->nzmax || lNNZcnt==0)
  return NULL;
 allocateAC(lowerHalf,A->ncol,lNNZcnt,0,true);
 size_t curNNZ=0;
 lowerHalf->p[0] = 0;
 //copying into new space
 for (int i = 0; i < A->ncol; ++i) {
  for (int j = A->p[i]; j < A->p[i+1]; ++j) {
   if(A->i[j]>=i){
    lowerHalf->i[curNNZ] = i;
    curNNZ++;
   }
  }
  lowerHalf->p[i+1] = curNNZ;
 }
 return lowerHalf;
}


#endif //CHOLOPENMP_UTIL_H
