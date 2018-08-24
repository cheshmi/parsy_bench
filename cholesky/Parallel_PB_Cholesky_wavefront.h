//
// Created by kazem on 3/20/18.
//

#ifndef PROJECT_PARALLEL_PB_CHOLESKY_WAVEFRONT_H
#define PROJECT_PARALLEL_PB_CHOLESKY_WAVEFRONT_H

#undef PRUNE
#define MKL
bool cholesky_left_par_waveFront(int n, int* c, int* r, double* values,
                                 size_t *lC, int* lR, size_t* Li_ptr, double* lValues,
                          int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                                 int *aTree, int *cT, int *rT, int *col2Sup,
#else
  int *prunePtr, int *pruneSet,
#endif
                          int nLevels, int *levelPtr, int *levelSet,
                          int chunk, int threads,int super_max
                          ,int col_max) {
 /*
  * For timing using BLAS
  */
 int top=0;
 int *xi;
 int *map ;
 double *contribs ;
 int info;
 double one [2], zero [2];
 one [0] =  1.0 ;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
 one [1] =  0. ;
 zero [0] = 0. ;     /* BETA for *syrk, *herk, and *gemm */
 zero [1] = 0. ;

 for (int lev = 0; lev < nLevels; ++lev) {
#pragma omp parallel private(xi,map, contribs)
  {
   map = new int[n]();
   contribs = new double[super_max*col_max]();
   xi = new int[2*n+1]();
   #pragma omp for  schedule(dynamic,chunk)
   for (int lIter = levelPtr[lev]; lIter < levelPtr[lev+1]; ++lIter) {
    int s = levelSet[lIter]+1;
    //printf("Thread %d has completed iteration %d.\n", omp_get_thread_num( ), s);
    int curCol = s!=0 ? blockSet[s - 1] : 0;
    int nxtCol = blockSet[s];
    int supWdt = nxtCol-curCol;
    int nSupR = Li_ptr[nxtCol]-Li_ptr[curCol];//row size of supernode
    for (int i = Li_ptr[curCol],cnt=0; i < Li_ptr[nxtCol]; ++i) {
     map[lR[i]] = cnt++;//mapping L rows position to actual row idx
    }
    //copy the columns from A to L
    for (int i = curCol; i < nxtCol; ++i) {//Copy A to L
     int pad=i-curCol;
     for (int j = c[i]; j < c[i+1] ; ++j) {
      // if(r[j]>=i)//does not need to save upper part.
      lValues[lC[i]+map[r[j]]] = values[j];
      //   else
      //      printf("dddd\n");
     }
    }

    double *src, *cur=&lValues[lC[curCol]];//pointing to first element of the current supernode
#ifndef PRUNE
    top = ereach_sn(supNo,cT,rT,curCol,nxtCol,col2Sup, aTree,xi,xi+supNo);
    assert(top>=0);
    for(int i = top; i < supNo; ++i){
     int lSN = xi[i];

#else
     for (int i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
      int lSN = pruneSet[i];
#endif

     int nSupRs = 0;
     int cSN = blockSet[lSN];//first col of current SN
     int cNSN = blockSet[lSN+1];//first col of Next SN
     size_t Li_ptr_cNSN = Li_ptr[cNSN];
     size_t Li_ptr_cSN = Li_ptr[cSN];
     int nSNRCur=Li_ptr_cNSN-Li_ptr_cSN;
     int  supWdts=cNSN-cSN;//The width of current src SN
     int lb=0,  ub=0;
     bool sw=true;
     for (int j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
      //finding the overlap between curCol and curCol+supWdt in the src col
      if (lR[j] >= curCol && sw) {
       //src*transpose(row lR[j])
       lb=j-Li_ptr_cSN;
       sw=false;
      }
      if(lR[j] < curCol+supWdt && !sw){
       ub=j-Li_ptr_cSN;
      }
     }
     nSupRs=Li_ptr_cNSN-Li_ptr_cSN-lb;
     int ndrow1=ub-lb+1;
     int ndrow3 = nSupRs-ndrow1;
     src=&lValues[lC[cSN]+lb];//first element of src supernode starting from row lb
     double *srcL = &lValues[lC[cSN]+ub+1];
#ifdef MKL
     dsyrk("L","N",&ndrow1,&supWdts,one,src,&nSNRCur,zero,
           contribs,&nSupRs);
#endif
#ifdef OPENBLAS
     dsyrk_("L","N",&ndrow1,&supWdts,one,src,&nSNRCur,zero,
                            contribs,&nSupRs);
#endif
#ifdef MYBLAS
     //TODO
#endif
     if(ndrow3>0){
#ifdef MKL
      dgemm("N","C",&ndrow3,&ndrow1,&supWdts,one,srcL,&nSNRCur,
            src,&nSNRCur,zero,&contribs[ndrow1],&nSupRs );
#endif
#ifdef OPENBLAS
      dgemm_("N","C",&ndrow3,&ndrow1,&supWdts,one,srcL,&nSNRCur,
                               src,&nSNRCur,zero,contribs+ndrow1,&nSupRs );
#endif
#ifdef MYBLAS
      //TODO
#endif

     }
     //copying contrib to L
     for (int i = 0; i < ndrow1; ++i) {//Copy contribs to L
      int col=map[lR[Li_ptr_cSN+i+lb]];//col in the SN
      for (int j = i; j < nSupRs ; ++j) {
       int cRow= lR[Li_ptr_cSN+j+lb];//corresponding row in SN
       //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
       cur[col*nSupR+map[cRow]] -= contribs[i*nSupRs+j];
      }
     }
    }//end top

#ifdef MKL
    dpotrf("L",&supWdt,cur,&nSupR,&info);
#endif
#ifdef OPENBLAS
    dpotrf_("L",&supWdt,cur,&nSupR,&info);
#endif
#ifdef MYBLAS
    Cholesky_col(nSupR,supWdt,cur);
#endif

    int rowNo=nSupR-supWdt;
#ifdef MKL
    dtrsm("R", "L", "C", "N", &rowNo, &supWdt,one,
          cur,&nSupR,&cur[supWdt],&nSupR);
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


    //        }

   }//end lIter
   delete []xi;
   delete []contribs;
   delete []map;
  }//omp
 }//end Level
 return true;
}
#endif //PROJECT_PARALLEL_PB_CHOLESKY_WAVEFRONT_H
