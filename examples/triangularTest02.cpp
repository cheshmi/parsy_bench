//
// Created by kazem on 2/27/18.
//

#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <cholmod.h>
#include <cholmod_function.h>
#include <omp.h>
#include <Triangular_BCSC.h>
#include <Inspection_Level.h>
#include <Triangular_CSC.h>

#include "Ordering.h"
#include "Inspection_Prune.h"
#include "Inspection_Block.h"
#include "Util.h"
#include "PB_Cholesky.h"
#include "LSparsity.h"
#include "mkl.h"
#include "parallel_PB_Cholesky_05.h"
#include "Cholesky_Perf_Modeling.h"
#include "Cholesky_Perf_Comm_modeling.h"

#define CPUTIME (SuiteSparse_time ( ))
#define BCSC_TRNG
#undef DEBUG
//#define FLOPCNT


using namespace std;



int main(int argc, char *argv[]) {

 std::string f1 = argv[1];
 int *colA, *rowA;
 double *valL;
 double   *valA;
 int  maxSupWid, maxCol;
 size_t n, nnzA;

 std::vector<profilingInfo> piArray;

 if (!readMatrix(f1,n,nnzA,colA,rowA,valA))
  return -1;
 std::string waste="/home/kazem/UFDB/rajat21.mtx";
 ifstream wasteFile;
 wasteFile.open(waste);
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

 // MKL_Set_Num_Threads(1);
 MKL_Domain_Set_Num_Threads(blasThreads,MKL_DOMAIN_BLAS);

 int *prunePtr, *pruneSet;
 int *levelPtr = NULL, *levelSet = NULL, *parPtr = NULL,
   *partition =NULL;
 int nLevels=0, nPar=0;

 double *timingChol = new double[4+numThread]();//for time measurement
 double orderingTime=0;
 int nrelax[3] = {4,16,48};//TODO
 double zrelax[3] = {0.8,0.1,0.05};
 int status=0;
 CSC *Amat = new CSC;
 Amat->nzmax = nnzA; Amat->ncol=Amat->nrow=n;
 Amat->stype=-1;Amat->xtype=CHOLMOD_REAL;Amat->packed=TRUE;
 Amat->p = colA; Amat->i = rowA; Amat->x=valA; Amat->nz = NULL;
 Amat->sorted = TRUE;
 BCSC *L = analyze_p2(1,Amat,NULL,NULL,nrelax,zrelax,
                      n,prunePtr,pruneSet,
                      nLevels, levelPtr,levelSet,
                      nPar, parPtr, partition,
                      costParam,levelParam,finalSeqNode,
                      status, maxSupWid, maxCol, orderingTime);

 valL = new double[L->xsize]();
 delete []L->pi;
 delete []L->i;
 delete []L->ColCount;

 CSC *A1 = ptranspose(Amat,2,L->Perm,NULL,0,status);
 CSC *A2 = ptranspose(A1,2,NULL,NULL,0,status);

 for (int i = 0; i < L->xsize; ++i) {
  valL[i]=0.0;
 }


 cholesky_left_par_05(n,A2->p,A2->i,A2->x,L->p,L->s,L->i_ptr,valL,
                      L->super,L->nsuper, timingChol,
#ifndef PRUNE
                      L->sParent,A1->p, A1->i, L->col2Sup,
#else
   prunePtr,pruneSet,
#endif
                      nLevels,levelPtr,levelSet,
                      nPar, parPtr, partition,
                      chunk, numThread, maxSupWid+1,maxCol+1);

 allocateAC(Amat,0,0,0,FALSE);
 allocateAC(A1,0,0,0,FALSE);
 allocateAC(A2,0,0,0,FALSE);
/*
 * ********************* Triangular Solve
 */

 int *col, *row;
 double  *y, *val, *x;
 int nnz;
 std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;
 double duration4 = 0 ,duration3 = 0, duration2=0, duration1=0;

 x=new double[n]();

#ifdef FLOPCNT
 //***************Serial
 int *ia = new int[n + 1];
 int *ja = new int[L->xsize];
 double *a = new double[L->xsize];
 bcsc2csc(n, L->nsuper, L->p, L->s, L->i_ptr, L->super, valL, ia, ja, a);
 unsigned long counts=0;
 rhsInit(n,ia,ja,a,x);
 counts = flopCoutLSolve(n,ia,ja,a,x);
 std::cout<<L->xsize<<";"<<counts<<";";
 delete []ia;
 delete []ja;
 delete []a;
#endif


#ifdef CSC
 //***************Serial
    rhsInit(n,col,row,val,x);
    start = std::chrono::system_clock::now();
    lsolve(n,col,row,val,x);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    duration1=elapsed_seconds.count();
    std::cout<<duration1<<",";


    //****************Parallel CSC
    int *levelPtr, *levelSet, levels=0;
    levels= buildLevelSet_CSC(n, nnz, col, row, levelPtr, levelSet);
    rhsInit(n,col,row,val,x);
    start = std::chrono::system_clock::now();
    lsolvePar(n,col,row,val,x,levels,levelPtr,levelSet, chunk);
    //lsolvePar2(n,col,row,val,x);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    duration2=elapsed_seconds.count();
    std::cout<<duration2<<","<< levels <<",";
#endif

#ifdef BCSC_TRNG
 //*************** BCSC
 int *levelbPtr, *levelbSet, blevels=0;
 int *col2sup = L->col2Sup;
 int supNo=0, newNNZ=0;
 int *sup2col = L->super;
 size_t *newCol = L->p;
 int *newRow = L->s;
 double *newVal = valL;
 size_t *rowP = L->i_ptr;
 size_t nBlocks = L->nsuper;
 double max_tmp=0.0;
 int nChild=NULL;
 int ATreeHeight = getTreeHeightBruteForce(nBlocks,L->sParent);

 std::cout<<f1<<","<<levelParam<<","<<finalSeqNode<<","<<n<<","<< ATreeHeight <<","
          <<nBlocks<<","<<nnz <<",";
 int iterno=5;
 //*************** Serial Blocked
 for (int j = 0; j < iterno; ++j) {
  rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
  start = std::chrono::system_clock::now();
  blockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration2 = elapsed_seconds.count();
  std::cout << duration2 << ",";
 }
 std::cout<<"*:";

#if 0
 for (int j = 0; j < blevels; ++j) {
        for (int i = levelbPtr[j]; i < levelbPtr[j+1]; ++i) {
            std::cout<<levelbSet[i]<<",";
        }
        std::cout<<"\n";
    }
#endif
 //*************** Parallel Blocked
 //omp_set_num_threads(1);
 /*blevels = buildLevelSet_BCSC(n,nnz,col,rowP,newRow,nBlocks,

                             sup2col,col2sup,levelbPtr,levelbSet);*/
 levelbPtr = new int[nBlocks]();
 levelbSet = new int[nBlocks]();
 blevels = getLevelSet(nBlocks,L->sParent,levelbPtr,levelbSet);
 for (int j = 0; j < iterno; ++j) {
  rhsInitBlocked(n,nBlocks,newCol,newRow,rowP,newVal,x);
  start = std::chrono::system_clock::now();
  //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
  leveledBlockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,nBlocks,
                       x,blevels,levelbPtr,levelbSet,chunk);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  duration2=elapsed_seconds.count();
  if(testTriangular(n, x))
   std::cout<<duration2<<",";
  else
   std::cout<<"H1 failed,";
 }

 std::cout<<"*:";
 //*************** Parallel Blocked H2
 //omp_set_num_threads(1);
 for (int j = 0; j < iterno; ++j) {
  rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
  start = std::chrono::system_clock::now();
  //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
  H2LeveledBlockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks,
                         x, nLevels, levelPtr, levelSet,
                         nPar, parPtr, partition, chunk);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration2 = elapsed_seconds.count();
  if (testTriangular(n, x))
   std::cout << duration2 << ",";
  else
   std::cout << "H2 failed,";
 }
 std::cout<<"*:";

 //*************** Parallel Blocked H2 Peeled
 //omp_set_num_threads(1);
 for (int j = 0; j < iterno; ++j) {
  rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
  start = std::chrono::system_clock::now();
  //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
  H2LeveledBlockedLsolve_Peeled(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks,
                         x, nLevels, levelPtr, levelSet,
                         nPar, parPtr, partition, chunk,numThread);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration2 = elapsed_seconds.count();
  if (testTriangular(n, x))
   std::cout << duration2 << ",";
  else
   std::cout << "H2 failed,";
 }
 std::cout<<"*:";
#endif



//TODO HACK
 delete []L->super;
 delete []L->sParent;
 delete []L->s;
 delete []L->col2Sup;
 delete []L->p;
 delete []L->i_ptr;



#if DEBUG > 0
 for (int i = n-10; i < n; ++i) {
            std::cout<<i<<":\n";
            for (int m = colL[i],cnt=0; m < colL[i+1]; ++m, ++cnt) {
                if(!std::isfinite(valL[m])) {
                    std::cout << "Error in colA "<< i;
                    return -1;
                }
                if(rowL[li_ptr[i]+cnt] >= i )
                    std::cout<<valL[m]<<",";
            }
            std::cout<<"\n";
        }
        std::cout<<"\n";
#endif
// delete []col2sup;
#ifdef PRUNE
 delete []prunePtr; delete []pruneSet;
#endif
 if(levelPtr!=NULL)
  delete []levelPtr;
 if(levelPtr!=NULL)
  delete []levelSet;
 if(parPtr!=NULL )
  delete []parPtr;
 if(partition!=NULL)
  delete []partition;
 //delete []contribs;
 //delete []map;
 delete []valL;
 //delete []colL;
 //delete []li_ptr;
 delete []timingChol;

 return 0;
}

