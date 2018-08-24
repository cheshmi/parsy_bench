//
// Created by kazem on 4/15/18.
//

#ifndef PROJECT_SEREIAL_CHOLESKY_01_H
#define PROJECT_SEREIAL_CHOLESKY_01_H
/*#define FLEXCOMPLEX
typedef long int  INT;
#define MKL_INT INT*/
#include <stdlib.h>
#include <cmath>

bool cholesky_left_01(int n, int* c, int* r, double* values,
                      MKL_INT* lC, MKL_INT* lR, double* &lValues, int *prunePtr, int *pruneSet) {
 /*
  * Performs a Cholesky decomposition on a given matrix (c, r, values), i.e.
  * stored in compressed column format, and produces L, which are
  * stored in column compressed format.
  * (n, c, r, values) : IN : input matrix
  * (lC, lR) : IN : The column and rwo sparsity patterns of L
  * (lValues) : OUT : the nonzero values of the L factor
  * (pruneSet, prunePtr) : IN : the row sparsity pattern of the factor L
  */
 double tmp;
 double *f = new double[n]();
 int *finger = new int[n];
 int spCol = 0;
 for (int colNo = 0; colNo < n; ++colNo) {
  //Uncompressing a col into a 1D array
  for (int nzNo = c[colNo]; nzNo < c[colNo + 1]; ++nzNo) {
   f[r[nzNo]] = values[nzNo];//Copying nonzero of the col
  }
  for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]; ++i) {
   spCol = pruneSet[i];
   bool sw=false;
   for (int l = lC[spCol]; l < lC[spCol + 1]; ++l) {
    if (lR[l] == colNo && !sw) {
     tmp = lValues[l];
     sw=true;
    }
    if(sw){
     f[lR[l]] -= lValues[l] * tmp;
    }
   }
  }
  if (f[colNo] <= 0)
   return false; //The matrix is not SPD
  double tmpSqrt = sqrt(f[colNo]);
  f[colNo] = 0;
  lValues[lC[colNo]] = tmpSqrt;
  for (int j = lC[colNo] + 1; j < lC[colNo + 1]; ++j) {
   lValues[j] = f[lR[j]] / tmpSqrt;
   f[lR[j]] = 0;
  }
 }
 delete[]finger;
 delete[]f;
 return true;
}

bool cholesky_left_04(int n, int* c, int* r, double* values,
                      MKL_INT* lC, MKL_INT* lR, double* &lValues, int *prunePtr, int *pruneSet) {
 /*
  * Performs a Cholesky decomposition on a given matrix (c, r, values), i.e.
  * stored in compressed column format, and produces L, which are
  * stored in column compressed format.
  * (n, c, r, values) : IN : input matrix
  * (lC, lR) : IN : The column and rwo sparsity patterns of L
  * (lValues) : OUT : the nonzero values of the L factor
  * (pruneSet, prunePtr) : IN : the row sparsity pattern of the factor L
  */

 double *f = new double[n]();
 int *finger = new int[n];
 //Determining the column pointer
 for (int k = 0; k < n; ++k) {
  finger[k] = lC[k];
 }
 int spCol = 0;
 for (int colNo = 0; colNo < n; ++colNo) {
  //Uncompressing a col into a 1D array
  for (int nzNo = c[colNo]; nzNo < c[colNo + 1]; ++nzNo) {
   f[r[nzNo]] = values[nzNo];//Copying nonzero of the col
  }
  for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]-1; ++i) {
   spCol = pruneSet[i];
   for (int l = lC[spCol]; l < lC[spCol + 1]; ++l) {
    if (lR[l] >= colNo) {
     f[lR[l]] -= lValues[l] * lValues[finger[spCol]];
    }
   }
   finger[spCol]++;
  }
  finger[colNo]++;//Skip diagonal row
  if (f[colNo] <= 0)
   return false; //The matrix is not SPD
  double tmpSqrt = sqrt(f[colNo]);
  f[colNo] = 0;
  lValues[lC[colNo]] = tmpSqrt;
  for (int j = lC[colNo] + 1; j < lC[colNo + 1]; ++j) {
   lValues[j] = f[lR[j]] / tmpSqrt;
   f[lR[j]] = 0;
  }
 }
 delete[]finger;
 delete[]f;
 return true;
}
#endif //PROJECT_SEREIAL_CHOLESKY_01_H
