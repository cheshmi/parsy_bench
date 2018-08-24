//
// Created by kazem on 1/29/18.
//

#ifndef CHOLOPENMP_CHOLESKY_PERF_MODELING_H
#define CHOLOPENMP_CHOLESKY_PERF_MODELING_H

#include <chrono>
#include <omp.h>
#include "MyBLAS.h"
#define TIMING
//#undef TIMING
#define TIMING1
#undef TIMING1
#define BLASTIMING
#undef BLASTIMING
//#define PRUNE
#undef PRUNE
#ifdef MKL
#include "mkl.h"
#include "Reach.h"
#include "performanceModel.h"

#endif
#ifdef OPENBLAS
#include "blas/cblas.h"
#endif

struct oneDBlock{
 double time;
 int supWdth,supRow,supSize;
 int contribNo, contribNNZs, allNNZs;
 long factFlops, contribFlop, allFlops;
 long flopsTruncated, NNZTruncated;
 oneDBlock():time(0.0),supRow(0),supWdth(0),supSize(0),
             contribNo(0),contribNNZs(0),factFlops(0),
             contribFlop(0),allNNZs(0),allFlops(0),
             NNZTruncated(0),flopsTruncated(0){}
 void sumAll(){allNNZs=supSize+contribNNZs;
  allFlops=factFlops+contribFlop;
 flopsTruncated=allFlops/100;
  NNZTruncated=allNNZs/100;
 }

 void addNoNNZFact(oneDBlock a){
  time+=a.time; supWdth+=a.supWdth;supRow+=a.supRow;allNNZs+=a.allNNZs;
  contribNo+=a.contribNo; contribNNZs+=a.contribNNZs; factFlops+=a.factFlops;
  contribFlop+=a.contribFlop; allFlops+=a.allFlops;
 }

 void divideByScalar(int num){
  time/=num; supWdth/=num;supRow/=num;allNNZs/=num;
  contribNo/=num; contribNNZs/=num; factFlops/=num;
  contribFlop/=num; allFlops/=num;
 }
};

bool cmpBlock1D(oneDBlock *a, oneDBlock *b){
 if(a->allNNZs > b->allNNZs )
  return true;
 return false;
}



bool cholesky_Perf_Model(int n, int* c, int* r, double* values,
                          size_t *lC, int* lR, size_t* Li_ptr, double* lValues,
                          int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                          int *aTree, int *cT, int *rT, int *col2Sup,
#else
  int *prunePtr, int *pruneSet,
#endif

                          int nLevels, int *levelPtr, int *levelSet,
                          int nPar, int *parPtr, int *partition,
                          int chunk, int threads, int super_max,
                          int col_max, double *nodCost=NULL) {
 /*
  * For timing using BLAS
  */
 MKL_Domain_Set_Num_Threads(1,MKL_DOMAIN_BLAS);
 int top = 0;
 std::vector<oneDBlock*> profiling;
 oneDBlock *curIter=NULL;
 // int *xi = new int[2*supNo]();
 //int col_max = n;
 int *map, *xi;
 double *contribs;
 double *blasTimePerThread = timing+4;
 int info;
 int thth=0;
 double one[2], zero[2];
 one[0] = 1.0;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
 one[1] = 0.;
 zero[0] = 0.;     /* BETA for *syrk, *herk, and *gemm */
 zero[1] = 0.;
 std::chrono::time_point<std::chrono::system_clock> start, end,startin,endin;
 std::chrono::duration<double> elapsed_seconds;
 double duration4 = 0 ,duration3 = 0, duration2=0, duration1=0;


 map = new int[n]();
 contribs = new double[super_max * col_max]();
 xi = new int[2*supNo]();
 for (int s = 1; s <= supNo; ++s) {
  curIter = new oneDBlock;
#ifdef TIMING
 start = std::chrono::system_clock::now();
#endif
   int curCol = s != 0 ? blockSet[s - 1] : 0;
   int nxtCol = blockSet[s];
   int supWdt = nxtCol - curCol;
   int nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
   for (int i = Li_ptr[curCol], cnt = 0; i < Li_ptr[nxtCol]; ++i) {
    map[lR[i]] = cnt++;//mapping L rows position to actual row idx
   }
   //copy the columns from A to L
   for (int i = curCol; i < nxtCol; ++i) {//Copy A to L
    int pad = i - curCol;
    for (int j = c[i]; j < c[i + 1]; ++j) {
     // if(r[j]>=i)//does not need to save upper part.
     lValues[lC[i] + map[r[j]]] = values[j];
     //   else
     //      printf("dddd\n");
    }
   }

   double *src, *cur = &lValues[lC[curCol]];//pointing to first element of the current supernode
#ifndef PRUNE
   top = ereach_sn(supNo,cT,rT,curCol,nxtCol,col2Sup, aTree,xi,xi+supNo);
  curIter->contribNo = supNo-top;
   for(int i = top; i < supNo; ++i){
    int lSN = xi[i];

#else
    for (int i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
      int lSN = pruneSet[i];
#endif
    int nSupRs = 0;
    int cSN = blockSet[lSN];//first col of current SN
    int cNSN = blockSet[lSN + 1];//first col of Next SN
    int Li_ptr_cNSN = Li_ptr[cNSN];
    int Li_ptr_cSN = Li_ptr[cSN];
    int nSNRCur = Li_ptr_cNSN - Li_ptr_cSN;
    int supWdts = cNSN - cSN;//The width of current src SN
    int lb = 0, ub = 0;
    bool sw = true;
    for (int j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
     //finding the overlap between curCol and curCol+supWdt in the src col
     if (lR[j] >= curCol && sw) {
      //src*transpose(row lR[j])
      lb = j - Li_ptr_cSN;
      sw = false;
     }
     if (lR[j] < curCol + supWdt && !sw) {
      ub = j - Li_ptr_cSN;
     }
    }
    nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
    curIter->contribNNZs+=(nSupRs*supWdts);
    int ndrow1 = ub - lb + 1;
    int ndrow3 = nSupRs - ndrow1;
    curIter->contribFlop+=(ndrow1*supWdts*ndrow3);
    src = &lValues[lC[cSN] + lb];//first element of src supernode starting from row lb
    double *srcL = &lValues[lC[cSN] + ub + 1];
#ifdef BLASTIMING
    startBlas = std::chrono::system_clock::now();
#endif
#ifdef MKL
    dsyrk("L", "N", &ndrow1, &supWdts, one, src, &nSNRCur, zero,
          contribs, &nSupRs);
#endif
#ifdef OPENBLAS
    dsyrk_("L","N",&ndrow1,&supWdts,one,src,&nSNRCur,zero,
                    contribs,&nSupRs);
#endif
#ifdef MYBLAS
    //TODO
#endif
    if (ndrow3 > 0) {
#ifdef MKL
     dgemm("N", "C", &ndrow3, &ndrow1, &supWdts, one, srcL, &nSNRCur,
           src, &nSNRCur, zero, &contribs[ndrow1], &nSupRs);
#endif
#ifdef OPENBLAS
     dgemm_("N","C",&ndrow3,&ndrow1,&supWdts,one,srcL,&nSNRCur,
                       src,&nSNRCur,zero,contribs+ndrow1,&nSupRs );
#endif
#ifdef MYBLAS
     //TODO
#endif
#ifdef BLASTIMING
     endBlas = std::chrono::system_clock::now();
     elapsed_seconds = (endBlas-startBlas);
     blasTimePerThread[threadID]+=elapsed_seconds.count();
#endif

    }
    //copying contrib to L
    for (int i = 0; i < ndrow1; ++i) {//Copy contribs to L
     int col = map[lR[Li_ptr_cSN + i + lb]];//col in the SN
     for (int j = i; j < nSupRs; ++j) {
      int cRow = lR[Li_ptr_cSN + j + lb];//corresponding row in SN
      //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
      cur[col * nSupR + map[cRow]] -= contribs[i * nSupRs + j];
     }
    }
   }
#ifdef BLASTIMING
   startBlas = std::chrono::system_clock::now();
#endif
#ifdef MKL
   dpotrf("L", &supWdt, cur, &nSupR, &info);
#endif
#ifdef OPENBLAS
   dpotrf_("L",&supWdt,cur,&nSupR,&info);
#endif
#ifdef MYBLAS
   Cholesky_col(nSupR,supWdt,cur);
#endif

   int rowNo = nSupR - supWdt;
#ifdef MKL
   dtrsm("R", "L", "C", "N", &rowNo, &supWdt, one,
         cur, &nSupR, &cur[supWdt], &nSupR);
#endif
#ifdef OPENBLAS
   dtrsm_("R", "L", "C", "N", &rowNo, &supWdt,one,
               cur,&nSupR,&cur[supWdt],&nSupR);
#endif
#ifdef MYBLAS
   for (int i = supWdt; i < nSupR; ++i) {
            lSolve_dense_col(nSupR,supWdt,cur,&cur[i]);
        }//TODO
#endif
#ifdef BLASTIMING
   endBlas = std::chrono::system_clock::now();
   elapsed_seconds = (endBlas-startBlas);
   blasTimePerThread[threadID]+=elapsed_seconds.count();
#endif

#ifdef TIMING
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 duration2=elapsed_seconds.count();
  curIter->time=duration2;
  curIter->supRow=nSupR;
  curIter->supWdth=supWdt;
  curIter->supSize=supWdt*nSupR;
  curIter->factFlops += OPS_TRSM(supWdt,rowNo);
  curIter->factFlops += OPS_PPF(supWdt);
  curIter->sumAll();
  profiling.push_back(curIter);
 /*std::cout<< s-1<<";"<<supWdt<<";"<<nSupR<<";"<<supWdt*nSupR<<";"<<numOfConts<<";"
          <<nnzOfContrib <<";"<<duration2<<"\n";*/

#endif
 }
 //Sorting profiling info
 std::vector<oneDBlock*> reducedProfiling;
 std::sort(profiling.begin(),profiling.end(),cmpBlock1D);
 oneDBlock *tmp;
 for (int l = 0; l < supNo; ++l) {
  //merge all profiling blocks with the same truncated values
  int count=0;
  tmp = new oneDBlock;
  tmp->NNZTruncated = profiling[l]->NNZTruncated;
  while(tmp->NNZTruncated == profiling[l]->NNZTruncated){
   tmp->addNoNNZFact(*profiling[l]);
   count++;l++;
   if(l==supNo)
    break;
  }
  tmp->divideByScalar(count); // averaging
  reducedProfiling.push_back(tmp);

 }
 for (int l = 0; l < reducedProfiling.size(); ++l) {
  std::cout<<reducedProfiling[l]->NNZTruncated<<";"<<reducedProfiling[l]->time<<
           ";"<<reducedProfiling[l]->allFlops<<";"<<reducedProfiling[l]->time<<";\n";
 }

 delete []xi;
 delete []contribs;
 delete []map;
 for (int k = 0; k < supNo; ++k) {
  delete profiling[k];
 }
 return true;
}


#endif //CHOLOPENMP_CHOLESKY_PERF_MODELING_H
