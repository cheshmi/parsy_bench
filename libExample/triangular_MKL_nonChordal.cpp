//
// Created by kazem on 7/13/18.
//

#include <stdio.h>
#define MKL_ILP64
#define FLEXCOMPLEX
typedef int  INT;
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
#include "amd.h"

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


 /*
  * Calling Cholesky to generate blocked triangular matrix
  */

 omp_set_num_threads(numThread);

 MKL_Set_Num_Threads(numThread);
 // MKL_Set_Num_Threads(1);
 MKL_Domain_Set_Num_Threads(numThread, MKL_DOMAIN_BLAS);

 int *prunePtr, *pruneSet;
 int *levelPtr = NULL, *levelSet = NULL, *parPtr = NULL,
   *partition = NULL;
 int nLevels = 0, nPar = 0;
 int ncol=n;
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
 int *Perm = new int[n], *inPerm;


#ifdef GIVEN
 //pastix_data_t **pastix_data;
 L->ordering = CHOLMOD_METIS;
 for (int l = 0; l < A->nrow; ++l) {
  Perm[l] = inPerm[l];
 }

#elif SCOTCH
 L->ordering = CHOLMOD_METIS;
 CSC *ATrans;
 unsigned long nnzFull = A->nzmax*2;//Symmetric case
 ATrans = ptranspose(A, 0, NULL, NULL, 0, status) ;

 SCOTCH_Num baseVal=0;
 SCOTCH_Num          vertnbr;                    /* Number of vertices */
 SCOTCH_Graph        grafdat;                    /* Source graph       */
 SCOTCH_Ordering     ordedat;                    /* Graph ordering     */
 SCOTCH_Strat        stradat;                    /* Ordering strategy  */
 SCOTCH_Num          straval;
 //Making the graph for passing it to metis, it should have
 //both upper and lower parts
 //allocateAC(AFull,ncol,nnzFull,0,TRUE);
 SCOTCH_Num *permtab = new SCOTCH_Num[ncol]();                    /* Permutation array  */
 SCOTCH_Num *AFullp = new SCOTCH_Num[ncol+1]();
 SCOTCH_Num *AFulli = new SCOTCH_Num[nnzFull]();
 SCOTCH_Num ncolIDXT = ncol;

 AFullp[0]=0;
 for (int i = 0; i < ncol; ++i) {
  int nnzOfCurCol = ATrans->p[i+1]-ATrans->p[i]-1;
  nnzOfCurCol += A->p[i+1]-A->p[i]-1;
  AFullp[i+1] =(long int) AFullp[i]+nnzOfCurCol;
  //copying Upper part, ignoring diagonal since it is in L already
  int base=AFullp[i];
  for (int j = ATrans->p[i],k=0; j < ATrans->p[i+1]-1; ++j,++k) {
   AFulli[base+k] =(long int) ATrans->i[j];
  }
  //copying L part
  base+=ATrans->p[i+1]-ATrans->p[i]-1;
  for (int j = A->p[i]+1,k=0; j < A->p[i+1]; ++j, ++k) {
   AFulli[base+k] =(long int) A->i[j];
  }
 }
 if(SCOTCH_stratInit (&stradat))
  return NULL;
 char straStr[550];
 sprintf(straStr, SCOTCH_STRAT_DIRECT);
 if(SCOTCH_stratGraphOrder(&stradat,straStr))
  return NULL;
 nnzFull=AFullp[ncol];
 if(SCOTCH_graphBuild (&grafdat,baseVal,ncol,AFullp,NULL,NULL,NULL,
                   nnzFull,AFulli,NULL))
  return NULL;
 if(SCOTCH_graphCheck(&grafdat))
  return NULL;
 /*if(SCOTCH_graphOrderList(&grafdat, ncol, NULL, &stradat,
                          permtab, NULL, NULL, NULL, NULL)) { *//* Create ordering *//*
  return NULL;
 }*/
 /*if(SCOTCH_graphOrderInit    (&grafdat, &ordedat, permtab, NULL, NULL, NULL, NULL)) {  Create ordering
  return NULL;
 }
 if(SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat)) {  Perform ordering
  return NULL;
 }*/

 if(SCOTCH_graphOrderInit    (&grafdat, &ordedat, permtab, NULL, NULL, NULL, NULL)) {  // Create ordering
  return NULL;
 }
 if(SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat)) { // Perform ordering
  return NULL;
 }

 if (SCOTCH_graphOrderCheck (&grafdat, &ordedat) != 0)
  return NULL;
 for (int i = 0; i < ncol; ++i) {
  Lperm[i]=permtab[i];
  //std::cout<<Lperm[i];
 }
 SCOTCH_graphOrderExit (&grafdat, &ordedat);
 SCOTCH_stratExit      (&stradat);
 SCOTCH_graphExit      (&grafdat);

 allocateAC(ATrans,ATrans->nrow,ATrans->nzmax,ATrans->stype,false);
 delete [] AFullp;
 delete [] AFulli;
 delete [] permtab;

#elif METIS
 CSC *ATrans;
 unsigned long nnzFull = Amat->nzmax * 2;//Symmetric case
 ATrans = ptranspose(Amat, 0, NULL, NULL, 0, status);
#if 0
 for (int i = 0; i < ncol; ++i) {
  for (int j = A->p[i]; j < A->p[i+1]; ++j) {
   std::cout<<A->i[j]<<";";
  }
  std::cout<<"\n";
 }
 std::cout<<"---\n";
 for (int i = 0; i < ncol; ++i) {
  for (int j = ATrans->p[i]; j < ATrans->p[i+1]; ++j) {
   std::cout<<ATrans->i[j]<<";";
  }
  std::cout<<"\n";
 }
 std::cout<<"==\n";
#endif

 //Making the graph for passing it to metis, it should have
 //both upper and lower parts
 //allocateAC(AFull,ncol,nnzFull,0,TRUE);
 idx_t options1[METIS_NOPTIONS];
 METIS_SetDefaultOptions(options1);

 idx_t *AFullp = new idx_t[ncol + 1]();
 idx_t *AFulli = new idx_t[nnzFull]();
 idx_t ncolIDXT = ncol;
 idx_t *weigt = new idx_t[ncol];
 idx_t *LpermIDX = new idx_t[ncol];
 idx_t *ILpermIDX = new idx_t[ncol];
 for (int i = 0; i < ncol; ++i) {
  LpermIDX[i] = 0;
  ILpermIDX[i] = 0;
  weigt[i] = 1;
 }
 AFullp[0] = 0;
 for (int i = 0; i < ncol; ++i) {
  int nnzOfCurCol = ATrans->p[i + 1] - ATrans->p[i] - 1;
  nnzOfCurCol += Amat->p[i + 1] - Amat->p[i] - 1;
  AFullp[i + 1] = (long int) AFullp[i] + nnzOfCurCol;
  //copying Upper part, ignoring diagonal since it is in L already
  int base = AFullp[i];
  for (int j = ATrans->p[i], k = 0; j < ATrans->p[i + 1] - 1; ++j, ++k) {
   AFulli[base + k] = (long int) ATrans->i[j];
  }
  //copying L part
  base += ATrans->p[i + 1] - ATrans->p[i] - 1;
  for (int j = Amat->p[i] + 1, k = 0; j < Amat->p[i + 1]; ++j, ++k) {
   AFulli[base + k] = (long int) Amat->i[j];
  }
 }
#if 0
 for (int i = 0; i < ncol; ++i) {
  for (int j = AFullp[i]; j < AFullp[i+1]; ++j) {
   std::cout<<AFulli[j]<<";";
  }
  std::cout<<"\n";
 }
 std::cout<<"\n";
#endif

 int retMet = METIS_NodeND(&ncolIDXT, AFullp, AFulli, NULL, options1,
                           LpermIDX, ILpermIDX);
 assert(retMet == METIS_OK);
 if (retMet != METIS_OK) {
  std::cout << " " << retMet << "\n";
  exit(10);
 }
 for (int i = 0; i < ncol; ++i) {
  Perm[i] = LpermIDX[i];
  //std::cout<<Lperm[i];
 }
 allocateAC(ATrans, ATrans->nrow, ATrans->nzmax, ATrans->stype, false);
 METIS_Free(AFullp);
 METIS_Free(AFulli);
 METIS_Free(weigt);
 METIS_Free(LpermIDX);
 METIS_Free(ILpermIDX);
#else
 double info[20]={0};
 double Control[2];
 Control [0] = 10; //TODO check later //AMD_Dense
 Control [1] = TRUE; //AMD_AGGRESSIVE
 L->ordering = CHOLMOD_AMD;
 amd_order(ncol,colA,rowA,Perm,NULL,info);
#endif


#ifdef VERIFY
 auto checkOrder = new bool[ncol]();
 for (int i = 0; i < ncol; ++i) checkOrder[i] = false;
 for (int i = 0; i < ncol; ++i) {
  checkOrder[Perm[i]] = true;
 }
 for (int i = 0; i < ncol; ++i) {
  assert(checkOrder[i] == true);
 }
 delete checkOrder;
#endif
 CSC *A1 = ptranspose(Amat, 2, Perm, NULL, 0, status);
 CSC *A2 = ptranspose(A1, 2, NULL, NULL, 0, status);





 allocateAC(Amat, 0, 0, 0, FALSE);
 allocateAC(A1, 0, 0, 0, FALSE);






 //*******************************************************************************
 //    Declaration of local variables:
 //*******************************************************************************
 //double *x = new double[n];

 //// ------- Variables and making CSC matrix
 MKL_INT stat;
 MKL_INT num_calls = 1000;
 double *x = new double[n];
 int nzA=A2->p[ncol];
 //Timing
 std::chrono::time_point<std::chrono::system_clock> start, end;
 double durationSimple = 0, duration2 = 0, duration1 = 0;
 double *y = new double[n]();
 double *y1 = new double[n]();
 double alpha = 1.0, beta = 0.0;
 double error = 0.0;
 //MKL_INT *ia = new MKL_INT[n + 1];
 /*MKL_INT *ia = (MKL_INT *) mkl_calloc(n+1, sizeof(MKL_INT), 64);

 MKL_INT *ja = (MKL_INT *) mkl_calloc(nzA, sizeof(MKL_INT), 64);
 double *a = (double *) mkl_calloc(nzA, sizeof(double), 64);*/
 MKL_INT *ia = (MKL_INT *) A2->p;

 MKL_INT *ja = (MKL_INT *) A2->i;
 double *a = A2->x;

 struct matrix_descr descrA;
 int iterNo = 5;
 std::chrono::duration<double> elapsed_seconds;
 //Converting BCSC to CSC
 //bcsc2csc(n, L->nsuper, L->p, L->s, L->i_ptr, L->super, valL, ia, ja, a);


 printf("%s; %d; ", f1.c_str(), numThread);
 ////--------- CSC
 printf(";CSC;");
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
 printf(";CSR;");
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
 double *csrVal = (double *) mkl_calloc(nzA, sizeof(double), 64);
 MKL_INT *csrColInd = (MKL_INT *) mkl_calloc(nzA, sizeof(MKL_INT), 64);
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
#if 0

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
/* mkl_free(ia);
 mkl_free(ja);
 mkl_free(a);*/
 mkl_sparse_destroy(cscA);
#endif

#endif
#if 1
 printf(";BCSR;");
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
allocateAC(A2, 0, 0, 0, FALSE);
 delete[]x;
 delete[]y;
 delete[]y1;

}