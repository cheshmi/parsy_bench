/*******************************************************************************
* Copyright 2004-2017 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
*   Content : Intel(R) MKL PARDISO C example
*
********************************************************************************
*/
/* -------------------------------------------------------------------- */
/* Example program to show the use of the "PARDISO" routine */
/* on symmetric linear systems */
/* -------------------------------------------------------------------- */
/* This program can be downloaded from the following site: */
/* www.pardiso-project.org */
/* */
/* (C) Olaf Schenk, Department of Computer Science, */
/* University of Basel, Switzerland. */
/* Email: olaf.schenk@unibas.ch */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <mkl.h>
#include "mkl_service.h"

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "omp.h"
#include "metis.h"
#include "Util.h"
#include "def.h"
#include "Transpose.h"

//#define DEBUG 1




int main(int argc, char *argv[]) {
 int stat=0;
 if(argc < 2){
  printf("input args are missing! \n");
  return -1;
 }
 std::string inpath;
 inpath = argv[1];
 int tmp =1;
 int iterNo=1;
 if (argc >= 3)
  tmp = atoi(argv[2]);
 if(argc >=4)
  iterNo = atoi(argv[3]);
 if(argc>=5)
  std::string orderFileName = argv[4];
 int num_thread = tmp; //Default is serial

 /* Matrix data. */
 double  *a ;
 int    *ja ;
 int    *ia ;
 int n, nnz;
 size_t n_t, nnz_t;

 int job[7] = {1, 1, 0, 0, 0,1, 1};
 int info =0;
 if (!readMatrix(inpath,n_t,nnz_t,ia,ja,a)){
  printf("Input file is not loaded! \n");
  return -1;
 }
 printf("loaded!\n");
 n = n_t; nnz = nnz_t;
/*
    double  *a = (double*)mkl_calloc(nnz,sizeof(double),64);
    MKL_INT    *ja = (int*)mkl_calloc(nnz,sizeof(int),64);
    MKL_INT    *ia = (int*)mkl_calloc(n+1,sizeof(int),64);;
#ifdef DEBUG
    printf("Allocation done! \n");
#endif
    mkl_dcsrcsc(job, &n, a, ja, ia, val, row, col, &info);*/

#if 0
 printf("%d , %d, %d \n",n,ia[n-1], ia[n-2]);
    for (int k = 0; k < n; ++k) {
        for (int j = ia[k]; j < ia[k + 1]; ++j) {
            printf("A[%d][%d]=%f \n",k+1,ja[j],a[j]);
        }
    }

    printf("Conversion done! \n");
#endif

#if 0
 printf("%d , %d, %d, %d \n",n,ia[n], ia[n-1],nnz);
    for (int k = 0; k < n+1; ++k) {
        printf("%d => %d \n",k, ia[k]);
    }
    for (int k = 0; k < nnz; ++k) {
        printf("%d => %d=%f \n",k, ja[k],a[k]);
    }

    printf("Conversion done! \n");
#endif

 MKL_INT mtype = 2;       /* Real symmetric positive  matrix */
 /* RHS and solution vectors. */
 //double b[8], x[8];
 double *b = (double*)mkl_calloc(n,sizeof(double),64);
 double *x = (double*)mkl_calloc(n,sizeof(double),64);
 MKL_INT nrhs = 1;     /* Number of right hand sides. */
 /* Internal solver memory pointer pt, */
 /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
 /* or void *pt[64] should be OK on both architectures */
 long int pt[64];
 /* Pardiso control parameters. */
 MKL_INT iparm[64];
 MKL_INT maxfct, mnum, phase, error, msglvl;
 /* Auxiliary variables. */
 MKL_INT i;
 double ddum;          /* Double dummy */
 MKL_INT idum;         /* Integer dummy. */

 /*---------------------------------------------------
  * Ordering
  -------------------------------------------------*/
 int *perm = new int[n];
#if GIVEN
   readOrdering(orderFileName,n,perm);
#else
 int status=0;
 CSC *ATrans, *A=new CSC;
 A->nzmax=nnz;
 A->nrow=A->ncol=n;
 A->x = a;
 A->p=ia;
 A->i=ja;
 A->stype=-1;
 A->packed=1;
 A->nz=NULL;
 long nnzFull = A->nzmax*2;//Symmetric case
 ATrans = ptranspose(A, 0, NULL, NULL, 0, status) ;

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
// assert(METIS_SetDefaultOptions(options1) == METIS_OK);
 //options[METIS_OPTION_NO2HOP] = 1;
/* options1[METIS_OPTION_PTYPE]=METIS_PTYPE_KWAY;
 options1[METIS_OPTION_OBJTYPE]=METIS_OBJTYPE_VOL;
 options1[METIS_OPTION_CTYPE]   = METIS_CTYPE_SHEM;
 options1[METIS_OPTION_IPTYPE]  = METIS_IPTYPE_NODE;
 options1[METIS_OPTION_RTYPE]   = METIS_RTYPE_GREEDY;
 options1[METIS_OPTION_NCUTS] = 1;
 options1[METIS_OPTION_NSEPS] = 1;
 options1[METIS_OPTION_NUMBERING] = 0;
 options1[METIS_OPTION_NITER] = 10;
 options1[METIS_OPTION_SEED] = -1;
 options1[METIS_OPTION_MINCONN] = 0;
 options1[METIS_OPTION_NO2HOP] = 1;
 options1[METIS_OPTION_CONTIG]  = 0;
 options1[METIS_OPTION_COMPRESS] = 0;
 options1[METIS_OPTION_CCORDER] = 1;
 options1[METIS_OPTION_PFACTOR]=100;
 options1[METIS_OPTION_UFACTOR] = 40;
 options1[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
 */
/*
 for (int i = 0; i < METIS_NOPTIONS; ++i) {
  options1[i]=1;
 }
*/
 idx_t *AFullp = new idx_t[A->ncol+1]();
 idx_t *AFulli = new idx_t[nnzFull]();
 idx_t ncolIDXT = A->ncol;
 idx_t *weigt = new idx_t[A->ncol];
 idx_t *LpermIDX = new idx_t[A->ncol];
 idx_t *ILpermIDX = new idx_t[A->ncol];
 for (int i = 0; i < A->ncol; ++i) {
  LpermIDX[i]=0;ILpermIDX[i]=0;weigt[i]=1;
 }
 AFullp[0]=0;
 for (int i = 0; i < A->ncol; ++i) {
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
#if 0
 for (int i = 0; i < ncol; ++i) {
  for (int j = AFullp[i]; j < AFullp[i+1]; ++j) {
   std::cout<<AFulli[j]<<";";
  }
  std::cout<<"\n";
 }
 std::cout<<"\n";
#endif

 int retMet= METIS_NodeND(&ncolIDXT,AFullp,AFulli,NULL,options1,
                          LpermIDX,ILpermIDX);
 //assert(retMet==METIS_OK);
 if(retMet!=METIS_OK){
  std::cout<<" "<<retMet<<"\n";
  exit(10);
 }
 for (int i = 0; i < A->ncol; ++i) {
  perm[i]=LpermIDX[i];
  //std::cout<<Lperm[i];
 }
 allocateAC(ATrans,ATrans->nrow,ATrans->nzmax,ATrans->stype,false);
 METIS_Free(AFullp);
 METIS_Free(AFulli);
 METIS_Free(weigt);
 METIS_Free(LpermIDX) ;
 METIS_Free(ILpermIDX);
#endif
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
 for ( i = 0; i < 64; i++ )
 {
  iparm[i] = 0;
 }
 iparm[0] = 1;         /* No solver default */
 iparm[1] = 2;         /* Fill-in reordering from METIS */
 iparm[3] = 0;         /* No iterative-direct algorithm */
 //iparm[4] = 0;         /* No user fill-in reducing permutation */
 iparm[4] = 1;         /* No user fill-in reducing permutation */
 iparm[5] = 0;         /* Write solution into x */
 iparm[6] = 0;         /* Not in use */
 iparm[7] = 2;         /* Max numbers of iterative refinement steps */
 iparm[8] = 0;         /* Not in use */
 iparm[9] = 8;        /* Perturb the pivot elements with 1E-8 */
 iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */
 iparm[11] = 0;        /* A^TX=B */
 iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
 iparm[13] = 0;        /* Output: Number of perturbed pivots */
 iparm[14] = 0;        /* Not in use */
 iparm[15] = 0;        /* Not in use */
 iparm[16] = 0;        /* Not in use */
 iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
 iparm[18] = 1;       /* Output: Mflops for LU factorization */
 iparm[19] = 0;        /* Output: Numbers of CG Iterations */


 //
 iparm[26] = 1;
 //iparm[23] = 1; //TODO: Later enable to se if the parallelism is better
 iparm[34] = 1;
 //Because iparm[4]==0 so:
 iparm[30] = 0;
 iparm[35] = 0;
 maxfct = 1;           /* Maximum number of numerical factorizations. */
 mnum = 1;         /* Which factorization to use. */
 msglvl = 1;           /* Print statistical information in file */
 error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
 for ( i = 0; i < 64; i++ )
 {
  pt[i] = 0;
 }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
//    pardisoinit (pt,  &mtype, iparm);
/*    iparm[35] = 0;*/
 //iparm[1] = 0;

 omp_set_num_threads(num_thread);
 //MKL_Domain_Set_Num_Threads(4,MKL_DOMAIN_BLAS);
 printf("yaaaaaaaay\n");
 phase = 11;
 PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
          &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
 if ( error != 0 )
 {
  printf ("\nERROR during symbolic factorization: %d", error);
  exit (1);
 }
 printf ("\nReordering completed ... ");
 printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
 printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
 phase = 22;
 for (int l = 0; l < iterNo; ++l) {
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
           &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if ( error != 0 )
  {
   printf ("\nERROR during numerical factorization: %d", error);
   exit (2);
  }
 }
 printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
 phase = 33;
 iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
 /* Set right hand side to one. */
 for ( i = 0; i < n; i++ )
 {
  b[i] = 1;
 }
 PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
          &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);
 if ( error != 0 )
 {
  printf ("\nERROR during solution: %d", error);
  exit (3);
 }
 printf ("\nSolve completed ... \n");
/*    printf ("\nThe solution of the system is: ");
    for ( i = 0; i < n; i++ )
    {
        printf ("\n x [%d] = % f", i, x[i]);
    }
    printf ("\n");*/
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
 phase = -1;           /* Release internal memory. */
 PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
          &n, &ddum, ia, ja, perm, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error);

/*    mkl_free(row);
    mkl_free(col);
    mkl_free(val);*/
 //mkl_free(ja);
 //mkl_free(ia);
 //mkl_free(a);
 delete []ia;
 delete []ja;
 delete []a;
 delete []perm;
 return 0;
}
