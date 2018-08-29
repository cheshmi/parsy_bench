//
// Created by kazem on 3/2/18.
//
#include <stdio.h>
#define MKL_ILP64
#define FLEXCOMPLEX
typedef long int  INT;
#define MKL_INT INT
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <mkl.h>
#include "mkl_spblas.h"
#include "mkl_service.h"

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "omp.h"
#include "metis.h"
#include "../common/Util.h"
#include "../common/def.h"
#include "../cholesky/parallel_PB_Cholesky_05.h"
#include "../cholesky/Transpose.h"
#include "../../../programs/metis-5.1.0/include/metis.h"
#include "../cholesky/LSparsity.h"
#include "MKL_Utils.h"

//#define DEBUG 1




int main(int argc, char *argv[]) {

 std::string f1 = argv[1];
 int *colA, *rowA;
 double *valL;
 double *valA;
 int maxSupWid, maxCol;
 size_t n, nnzA;

 //std::vector<profilingInfo> piArray;

 if (!readMatrix(f1, n, nnzA, colA, rowA, valA))
  return -1;

/* if(wasteFile.fail())
  return -1;*/
 int numThread = atoi(argv[2]);
 int chunk = atoi(argv[3]);
 int costParam = atoi(argv[4]);//Inner parts
 int levelParam = atoi(argv[5]);// level distance
 int blasThreads = atoi(argv[6]);
 int finalSeqNode = atoi(argv[7]);

 /*
  * Calling Cholesky to generate blocked triangular matrix
  */

 omp_set_num_threads(numThread);

 MKL_Set_Num_Threads(numThread);
 // MKL_Set_Num_Threads(1);
 MKL_Domain_Set_Num_Threads(blasThreads, MKL_DOMAIN_BLAS);

 int *prunePtr, *pruneSet;
 int *levelPtr = NULL, *levelSet = NULL, *parPtr = NULL,
   *partition = NULL;
 int nLevels = 0, nPar = 0;

 double *timingChol = new double[4 + numThread]();//for time measurement
 double orderingTime = 0;
 int nrelax[3] = {4, 16, 48};//TODO
 double zrelax[3] = {0.8, 0.1, 0.05};
 int status = 0;
 CSC *Amat = new CSC;
 Amat->nzmax = nnzA;
 Amat->ncol = Amat->nrow = n;
 Amat->stype = -1;
 Amat->xtype = CHOLMOD_REAL;
 Amat->packed = TRUE;
 Amat->p = colA;
 Amat->i = rowA;
 Amat->x = valA;
 Amat->nz = NULL;
 Amat->sorted = TRUE;
 BCSC *L = analyze_p2(1, Amat, NULL, NULL, nrelax, zrelax,
                      n, prunePtr, pruneSet,
                      nLevels, levelPtr, levelSet,
                      nPar, parPtr, partition,
                      costParam, levelParam, finalSeqNode,
                      status, maxSupWid, maxCol, orderingTime);

 valL = new double[L->xsize]();
 delete[]L->pi;
 delete[]L->i;
 delete[]L->ColCount;

 CSC *A1 = ptranspose(Amat, 2, L->Perm, NULL, 0, status);
 CSC *A2 = ptranspose(A1, 2, NULL, NULL, 0, status);

 for (int i = 0; i < L->xsize; ++i) {
  valL[i] = 0.0;
 }


 cholesky_left_par_05(n, A2->p, A2->i, A2->x, L->p, L->s, L->i_ptr, valL,
                      L->super, L->nsuper, timingChol,
#ifndef PRUNE
                      L->sParent, A1->p, A1->i, L->col2Sup,
#else
   prunePtr,pruneSet,
#endif
                      nLevels, levelPtr, levelSet,
                      nPar, parPtr, partition,
                      chunk, numThread, maxSupWid + 1, maxCol + 1);

 allocateAC(Amat, 0, 0, 0, FALSE);
 allocateAC(A1, 0, 0, 0, FALSE);
 allocateAC(A2, 0, 0, 0, FALSE);





 //*******************************************************************************
 //    Declaration of local variables:
 //*******************************************************************************
 //double *x = new double[n];

 //// ------- Variables and making CSC matrix
 MKL_INT stat;
 MKL_INT num_calls = 1000;
 double *x = new double[n];
 long int nnzL = L->xsize;
 //Timing
 std::chrono::time_point<std::chrono::system_clock> start, end;
 double durationSimple = 0, duration2 = 0, duration1 = 0;
 double *y = new double[n]();
 double *y1 = new double[n]();
 double alpha = 1.0, beta = 0.0;
 double error = 0.0;
 //MKL_INT *ia = new MKL_INT[n + 1];
 MKL_INT *ia = (MKL_INT *) mkl_calloc(n+1, sizeof(MKL_INT), 64);

 MKL_INT *ja = (MKL_INT *) mkl_calloc(L->xsize, sizeof(MKL_INT), 64);
 double *a = (double *) mkl_calloc(L->xsize, sizeof(double), 64);;

 struct matrix_descr descrA;
 int iterNo = 5;
 std::chrono::duration<double> elapsed_seconds;
 //Converting BCSC to CSC
 bcsc2csc(n, L->nsuper, L->p, L->s, L->i_ptr, L->super, valL, ia, ja, a);

 delete []L->super;
 delete []L->sParent;
 delete []L->s;
 delete []L->col2Sup;
 delete []L->p;
 delete []L->i_ptr;
 delete []valL;
 if(levelPtr!=NULL)
  delete []levelPtr;
 if(levelPtr!=NULL)
  delete []levelSet;
 if(parPtr!=NULL )
  delete []parPtr;
 if(partition!=NULL)
  delete []partition;

 printf("%s; %d; ", f1.c_str(), numThread);
 ////--------- CSC
 printf(";*;");
 // Descriptor of main sparse matrix properties
 // Structure with sparse matrix stored in CSR format
 sparse_matrix_t cscA;
 //bcsc2csc()


#if 0 //
 rhsInitBlocked(n,L->nsuper,L->p,L->s,L->i_ptr,valL,y);
 for (int j = 0; j < n; ++j) {
  assert(x[j] == y[j]);
 }
#endif

 // Create handle with matrix stored in CSC format

 stat = mkl_sparse_d_create_csc(&cscA, SPARSE_INDEX_BASE_ZERO,
                         n,  // number of rows
                         n,  // number of cols
                         ia,
                         ia + 1,
                         ja,
                         a);
 if (stat != SPARSE_STATUS_SUCCESS) {
  printf("executor failed on CSC;");
 }
 // Create matrix descriptor
 descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
 descrA.mode = SPARSE_FILL_MODE_LOWER;
 descrA.diag = SPARSE_DIAG_NON_UNIT;

 //// Inspector is not supported for CSC
/* stat = mkl_sparse_set_sv_hint(cscA, SPARSE_OPERATION_NON_TRANSPOSE,
                               descrA, num_calls);
 if (stat != SPARSE_STATUS_SUCCESS) {
  printf("analysis failed with %d\n", stat);
  return -1;
 }*/


 //// Executor**********************************************************************
// omp_set_num_threads(numThread);
 for (int i = 0; i < iterNo; ++i) {
  for (int jj = 0; jj < n; ++jj) {
   y[jj] = 0;
   x[jj]=1;
  }
  rhsInit(n, ia, ja, a, x);
  start = std::chrono::system_clock::now();
  stat = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE,
                           alpha, cscA, descrA, x, y);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration1 = elapsed_seconds.count();

  if (stat != SPARSE_STATUS_SUCCESS) {
   printf("executor failed on CSC;");
   return -1;
  }
  if (!testTriangular(n, y)) {
   printf(" %f CSC TRNS failed;", error);
  } else {
   printf(" %f ;", duration1);
  }

 }


#if 1
 //// For timing y1 = A * y
 for (int jj = 0; jj < n; ++jj) {
  y1[jj] = 0;
 }
 start = std::chrono::system_clock::now();
 mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                 alpha,
                 cscA,
                 descrA,
                 y,
                 beta,
                 y1);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 duration1 = elapsed_seconds.count();
 //printf(" %f ;",duration1);

 //MY SPMV
 // For validation perform y1 = A * y
#if 0
 for (int jj = 0; jj < n; ++jj) {
  y1[jj] = 0;
 }
 start = std::chrono::system_clock::now();
 spmv_csc(n, ia, ja, a, y, y1);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 duration1 = elapsed_seconds.count();
 //printf(" %f ;",duration1);
 // Release matrix handle and deallocate matrix
#endif

 ////----------CSR
 printf(";*;");
 sparse_matrix_t csrA;
 MKL_INT job[7] = {1, 0, 0, 0, 0, 1, 1};
 MKL_INT info = 0;
 const MKL_INT size_csr = n;
 char transa;
 char matdescra[6];
 transa = 'n';
 matdescra[0] = 't';
 matdescra[1] = 'l';
 matdescra[2] = 'n';
 matdescra[3] = 'c';
 alpha = 1.0;
 double *csrVal = (double *) mkl_calloc(nnzL, sizeof(double), 64);
 MKL_INT *csrColInd = (MKL_INT *) mkl_calloc(nnzL, sizeof(MKL_INT), 64);
 MKL_INT *csrRowPtr = (MKL_INT *) mkl_calloc(n + 1, sizeof(MKL_INT), 64);
 for (int j = 0; j < n; ++j) {
  x[j] = 0.0;
  y[j] = 0.0;
 }

 mkl_dcsrcsc(job, &size_csr, csrVal, csrColInd, csrRowPtr, a, ja, ia, &info);


 // Create handle with matrix stored in CSR format
 stat = mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO,
                                size_csr,  // number of rows
                                size_csr,  // number of cols
                                csrRowPtr,
                                csrRowPtr + 1,
                                csrColInd,
                                csrVal);
 if (stat != SPARSE_STATUS_SUCCESS) {
  printf("analysis failed with %d;", stat);
//  return -1;
 }
 // Create matrix descriptor
 descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
 descrA.mode = SPARSE_FILL_MODE_LOWER;
 descrA.diag = SPARSE_DIAG_NON_UNIT;

 //omp_set_num_threads(numThread);



 //// Inspector -
 start = std::chrono::system_clock::now();
 stat = mkl_sparse_set_sv_hint(csrA, SPARSE_OPERATION_NON_TRANSPOSE,
                               descrA, num_calls);
 if (stat != SPARSE_STATUS_SUCCESS) {
  printf("CSR analysis failed with %d ;", stat);
 // return -1;
 }
 stat = mkl_sparse_optimize(csrA);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 duration1 = elapsed_seconds.count();

 if (stat != SPARSE_STATUS_SUCCESS) {
  printf("CSR Applying analysis failed with %d;", stat);
 // return -1;
 } else {
  printf(" %f ;", duration1);
 }

 for (int i = 0; i < iterNo; ++i) {

  rhsInit(n, ia, ja, a, x);

#if 0//For testing CSR
  for (int k = 0; k < n; ++k) {
   double tmp=0;
   for (int i = csrRowPtr[k]; i < csrRowPtr[k+1]; ++i) {
    tmp+=csrVal[i];
   }
   assert(tmp-x[k]< 0.001);
  }
#endif
  for (int jj = 0; jj < n; ++jj) {
   y[jj] = 0;
  }
  start = std::chrono::system_clock::now();
  // Compute y = alpha * A^{-1} * x
  stat = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE,
                           alpha,
                           csrA,
                           descrA,
                           x,
                           y);
  /* mkl_dcsrsv(&transa, &size_csr, &alpha, matdescra,
              csrVal,
              csrColInd,
              csrRowPtr,
              csrRowPtr+1,
              x,
              y);*/
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration1 = elapsed_seconds.count();

  if (stat != SPARSE_STATUS_SUCCESS) {
   printf("CSR trsv failed with %d ;", stat);
  // return -1;
  } else
   printf(" %f ;", duration1);
 }


 // For validation perform y1 = A * y
 start = std::chrono::system_clock::now();
 mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                 alpha,
                 csrA,
                 descrA,
                 y,
                 beta,
                 y1);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 duration1 = elapsed_seconds.count();

 //printf(" %f ;", duration1);


 /*if (!testTriangular(n,y)) {
  printf(" %f CSR TRNS failed;", error);
 }*/

 for (int i = 0; i < n; i++) {
  error += x[i] - y1[i];
 };
 if (error > 1e-8) {
  printf("; %f   VALIDATION FAILED for CSR;", error);
 }
 mkl_free(ia);
 mkl_free(ja);
 mkl_free(a);
 mkl_sparse_destroy(cscA);
#endif
#if 1
 printf(";*;");
 //// BCSR
 const int blockSize = 2;
 int ActualNNZ = ia[n];
 char sparsity = 'L';
 char twoArgs = 'N';
 CRSArrays csrL;
 BCRSArrays bcrL;
 sparse_matrix_t bsrA;
 int m = n;

 csrL.a = csrVal;
 csrL.ia = csrRowPtr;
 csrL.ja = csrColInd;
 csrL.m = n;
 csrL.nnz = ActualNNZ;

 double *x_v = new double[blockSize * n]();
 double *y_v = new double[blockSize * n]();
 double *y1_v = new double[blockSize * n]();
 CRS_to_BCRS(csrL, &bcrL, blockSize);

 mkl_sparse_destroy(csrA);
 mkl_free(csrVal);
 mkl_free(csrRowPtr);
 mkl_free(csrColInd);
/*
  int maxBlockNo = (m + blockSize -1) / blockSize;
 double *absr = new double[maxBlockNo*blockSize*blockSize]();
 int *jab = new int[maxBlockNo+1]();
 int *iab = new int[maxBlockNo+1];
 int ldabsr = blockSize*blockSize;
 MKL_INT job_bsr[6] = {0,//If job(1)=0, the matrix in the CSR format is converted to the BSR format;
                   0,//If job(2)=0, zero-based indexing for the matrix in CSR format is used;
                   0,//If job(3)=0, zero-based indexing for the matrix in the BSR format is used;
                   0,
                   0,
                   1 //If job(6)>0, all output arrays absr, jab, and iab are filled in for the output storage.
 };
 mkl_dcsrbsr(job_bsr,&m,&blockSize,&ldabsr,
             csrVal,csrColInd,csrRowPtr,
             absr,jab,iab,&info);*/
/* stat = mkl_sparse_d_create_bsr ( &bsrA, SPARSE_INDEX_BASE_ZERO,
                                  SPARSE_LAYOUT_ROW_MAJOR,
                                  n,n,
                                  blockSize,
                                  csrRowPtr,
                                  csrRowPtr+1,
                                  csrColInd,
                                  csrVal);*/
 start = std::chrono::system_clock::now();
 stat = mkl_sparse_d_create_bsr(&bsrA, SPARSE_INDEX_BASE_ZERO,
                                SPARSE_LAYOUT_ROW_MAJOR,
                                bcrL.m,  // number of rows
                                bcrL.m,  // number of cols
                                blockSize, // block size
                                bcrL.ia,
                                bcrL.ia + 1,
                                bcrL.ja,
                                bcrL.a);
 if (stat != SPARSE_STATUS_SUCCESS)
  printf("Handle not created for BSR! %d ;", stat);

 // Create matrix descriptor
 descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
 descrA.mode = SPARSE_FILL_MODE_LOWER;
 descrA.diag = SPARSE_DIAG_NON_UNIT;
 stat = mkl_sparse_set_sv_hint(bsrA, SPARSE_OPERATION_NON_TRANSPOSE,
                               descrA, num_calls);

 if (stat != SPARSE_STATUS_SUCCESS) {
  printf("CSR analysis failed with %d ;", stat);
  //return -1;
 }
 stat = mkl_sparse_optimize(bsrA);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 duration1 = elapsed_seconds.count();
 if (stat != SPARSE_STATUS_SUCCESS) {
  printf("BSR analysis failed %d ;", stat);
 } else
  printf(" %f ;", duration1);

 for (int i = 0; i < iterNo; ++i) {
// rhsInit(n,ia,ja,a,x);

  for (int jj = 0; jj < n; ++jj) {
   y[jj] = 0;
   x[jj] = 0;
   x_v[jj] = 1.0;
  }
  start = std::chrono::system_clock::now();
  // Compute y = alpha * A^{-1} * x
  stat = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE,
                           alpha,
                           bsrA,
                           descrA,
                           x_v,
                           y_v);
  /*mkl_dbsrtrsv(&sparsity,&twoArgs,&twoArgs, &bcrL.nbBlockRows, &bcrL.lb,
               bcrL.a, bcrL.ia, bcrL.ja,
               x, y);*/
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration1 = elapsed_seconds.count();
  if (stat != SPARSE_STATUS_SUCCESS) {
   printf("BSR trsv failed with %d ,", stat);
  } else
   printf(" %f ;", duration1);
 }

/* if (!testTriangular(n,y)) {
  printf(" %f BSR TRNS failed;", error);
 }*/
 // For validation perform y1 = A * y
 start = std::chrono::system_clock::now();
/* mkl_cspblas_dbsrgemv ( &twoArgs,
                               &bcrL.nbBlockRows , &bcrL.lb ,
                               bcrL.a , bcrL.ia , bcrL.ja , x ,y);*/
 mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE,
                   alpha,
                   bsrA,
                   descrA,
                   y_v,
                   beta,
                   y1_v );
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 duration1=elapsed_seconds.count();
 for (int i = 0; i < n; i++ )
 {
  error += x[i]-y1_v[i];
 };
 if (error > 1e-8)
 {
  printf("; %f   VALIDATION FAILED for BSR;", error);
 }


 //printf(" %f ;\n", duration1);
 mkl_sparse_destroy(bsrA);
  delete []y_v;
 delete []y1_v;
#endif



/* mkl_sparse_destroy(csrA);
 mkl_free(csrVal);
 mkl_free(csrRowPtr);
 mkl_free(csrColInd);*/



/* delete []csrVal;
 delete []csrRowPtr;
 delete []csrColInd;*/

 delete[]x;
 delete[]y;
 delete[]y1;

}
