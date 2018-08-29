#include <omp.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <chrono>
#include "Util.h"
#include "Triangular_CSC.h"
#include "Inspection_Level.h"
#include "Inspection_Block.h"
#include "Triangular_BCSC.h"

#undef CSC
#undef BCSC
#define RUNALL

int main(int argc, char *argv[])  {
 std::string fName = argv[1];
#ifndef RUNALL
 int numThreads = atoi(argv[2]);
    int chunk = atoi(argv[3]);
    omp_set_num_threads(numThreads);
#endif
 int *col, *row;
 double  *y, *val, *x;
 int n, nnz;
 std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;
 double duration4 = 0 ,duration3 = 0, duration2=0, duration1=0;
 if (!readMatrix_old(fName,n,nnz,col,row,val))
  return -1;
 x=new double[n]();

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

#ifdef BCSC
 //*************** BCSC
    int *levelbPtr, *levelbSet, blevels=0;
    rhsInit(n,col,row,val,x);
    int *col2sup = new int[n];
    int supNo=0, newNNZ=0, newRowSize=0;
    superNodeDetection(n,col,row,col2sup,supNo);
    int *sup2col = new int[supNo];
    int *newCol = new int[n+1];
    calcSize(n,col,newCol,col2sup,sup2col,supNo,newRowSize,newNNZ);
    //int average = averageSupNode(sup2col,supNo);
    int *newRow = new int[newRowSize+1];
    double *newVal = new double[newNNZ];
    int *rowP = new int[n+1];
    createFormat(n,col,row,val,nnz,newRow,newRowSize,newVal,rowP,newCol,
                 col2sup,sup2col,supNo);
    blevels = buildLevelSet_BCSC(n,nnz,col,rowP,newRow,supNo,
                                 sup2col,col2sup,levelbPtr,levelbSet);
    //*************** Serial Blocked
    start = std::chrono::system_clock::now();
    blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    duration2=elapsed_seconds.count();
    std::cout<<duration2<<",";

#if 0
    for (int j = 0; j < blevels; ++j) {
        for (int i = levelbPtr[j]; i < levelbPtr[j+1]; ++i) {
            std::cout<<levelbSet[i]<<",";
        }
        std::cout<<"\n";
    }
#endif
    //*************** Parallel Blocked
    rhsInit(n,col,row,val,x);
    start = std::chrono::system_clock::now();
    //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
    leveledBlockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,
                         x,blevels,levelbPtr,levelbSet,chunk);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    duration2=elapsed_seconds.count();
    std::cout<<fName<<","<<n<<duration2<<","<< blevels <<","<<newNNZ<<","
             <<supNo<<","<<nnz <<"\n";
#endif
#ifndef RUNALL
 //Testing
    int test=0;
    for (int i = 0; i < n; ++i) {
        //test+=x[i];
        if(1-x[i]<0.001)
            test++;
    }
    if(n-test>0.001)
        std::cout<<"Error margin is high:"<<n-test<<",";
    /*else
        std::cout<<"WELL DONE!\n";*/
#endif

#ifdef RUNALL
 int chunk=10, num_threads = atoi(argv[2]), num_chunk = 1;
 int chunks[20]={5,10,20,30,40,50,60,70,80,90,
                 100,150,200,250,300,350,400,450,500,550 };
 std::cout<<fName<<","<<n<<","<<num_threads<<","<<num_chunk<<",";

 //***************Serial, Dense RHS
 rhsInit(n,col,row,val,x);
 start = std::chrono::system_clock::now();
 lsolve(n,col,row,val,x);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 duration1=elapsed_seconds.count();
 std::cout<<duration1<<",";

 //***************Serial vectorized, Dense RHS
 rhsInit(n,col,row,val,x);
 start = std::chrono::system_clock::now();
 lsolveVectorize(n,col,row,val,x);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 duration1=elapsed_seconds.count();
 std::cout<<duration1<<",";

 //***************Pruned
 int rhsPercent=0.03*n;
 int *xi = new int[2*n]();
 double *Bx = new double[3];
 int *Bp = new int[2]; Bp[0]=0; Bp[1]=rhsPercent;
 int *Bi = new int[rhsPercent]; //Bi[0]=4000;Bi[1]=4001;
 for (int i = 0; i < n; ++i) {
  x[i] = 0;
 }
 for (int i = n-rhsPercent,cnt=0; i < n; ++i) {
  x[i]=1;
  Bi[cnt++]=i;
 }
 start = std::chrono::system_clock::now();
 lsolve_reach_dec(n,col,row,val,Bp,Bi,Bx,0,xi,x,0,duration1);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 duration1=elapsed_seconds.count();
 std::cout<<duration1<<",";

 //***************Serial vectorized, Sparse RHS
 rhsPercent=0.03*n;
 for (int i = 0; i < n; ++i) {
  x[i] = 0;
 }
 for (int i = n-rhsPercent,cnt=0; i < n; ++i) {
  x[i]=1;
 }
 start = std::chrono::system_clock::now();
 lsolveVectorize(n,col,row,val,x);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 duration1=elapsed_seconds.count();
 std::cout<<duration1<<",";

 // Making level set
 int *levelPtr, *levelSet, levels=0;
 levels= buildLevelSet_CSC(n, nnz, col, row, levelPtr, levelSet);

 //*************** BCSC
 int *levelbPtr, *levelbSet, blevels=0;
 rhsInit(n,col,row,val,x);
 int *col2sup = new int[n];
 int supNo=0, newNNZ=0, newRowSize=0;
 superNodeDetection(n,col,row,col2sup,supNo);
 int *sup2col = new int[supNo];
 int *newCol = new int[n+1];
 calcSize(n,col,newCol,col2sup,sup2col,supNo,newRowSize,newNNZ);
 //int average = averageSupNode(sup2col,supNo);
 int *newRow = new int[newRowSize+1];
 double *newVal = new double[newNNZ];
 int *rowP = new int[n+1];
 createFormat(n,col,row,val,nnz,newRow,newRowSize,newVal,rowP,newCol,
              col2sup,sup2col,supNo);
 blevels = buildLevelSet_BCSC(n,nnz,col,rowP,newRow,supNo,
                              sup2col,col2sup,levelbPtr,levelbSet);
 //*************** Serial Blocked
 start = std::chrono::system_clock::now();
 blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 duration2=elapsed_seconds.count();
 std::cout<<duration2<<",";
 for (int th = 1; th <= num_threads; ++th) {
  omp_set_num_threads(th);
  for (int ch = 0; ch < num_chunk; ++ch) {
   chunk=chunks[ch];

   //************** Serial Pruned Blocked
   double durationAnalysis=0;
   int *PBset = new int[supNo];
   int top =reach_sn(n,col,row,Bp,Bi,0,PBset,0,supNo,col2sup,durationAnalysis);
   start = std::chrono::system_clock::now();
   lsolve_reach_dec(n,col,row,val,Bp,Bi,Bx,0,xi,x,0,duration1);
   blockedPrunedLSolve(n,newCol,newRow,newVal,nnz,rowP,PBset,top,sup2col,supNo,x);
   end = std::chrono::system_clock::now();
   elapsed_seconds = end-start;
   duration1=elapsed_seconds.count();
   std::cout<<duration1<<",";
   delete []PBset;

   //****************Parallel CSC
   rhsInit(n,col,row,val,x);
   //if(th <= 5){
   start = std::chrono::system_clock::now();
   lsolvePar(n,col,row,val,x,levels,levelPtr,levelSet, chunk);
   end = std::chrono::system_clock::now();
   /*}else{
       end=start;
   }*/
   elapsed_seconds = end-start;
   duration3=elapsed_seconds.count();
   std::cout<<duration3<<",";

   //*************** Parallel Blocked
   rhsInit(n,col,row,val,x);
   start = std::chrono::system_clock::now();
   leveledBlockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,
                        supNo,x,blevels,levelbPtr,levelbSet,chunk);
   end = std::chrono::system_clock::now();
   elapsed_seconds = end-start;
   duration4=elapsed_seconds.count();
   std::cout<<duration4<<",";

  }
 }
 std::cout<<"\n";
#endif
}