//
// Created by kazem on 7/18/17.
//

#ifndef TRIANGOPENMP_INSPECTION_H
#define TRIANGOPENMP_INSPECTION_H

#include <assert.h>
/*
 *
 */
int buildLevelSet_CSC(size_t n, size_t nnz, int *Lp, int *Li, int *&levelPtr,
                      int *&levelSet){
 int begin=0,end=n-1;
 int curLevel=0, curLevelCol=0;
 levelPtr = new int[n+1]();
 levelSet = new int[n]();
 int *inDegree = new int[n]();
 bool *visited = new bool[n]();
 for (int i = 0; i < Lp[n]; ++i) {//O(nnz)
  inDegree[Li[i]]++;
 }
#if 0
 for (int k = 0; k < n; ++k) {
        std::cout<<inDegree[k]<<",";
    }
    std::cout<<"\n";
#endif
 while(begin <= end){
  for (int i = begin; i <= end; ++i) {//For level curLevel
   if (inDegree[i] == 1 && !visited[i]) {//if no incoming edge
    visited[i] = true;
    levelSet[curLevelCol] = i; //add it to current level
    curLevelCol++;//Adding to level-set
   }
  }
  curLevel++;//all nodes with zero indegree are processed.
  levelPtr[curLevel]=curLevelCol;
  while(inDegree[begin]==1)
   begin++;
  while(inDegree[end]==1 && begin <= end)
   end--;
  //Updating degrees after removing the nodes
  for (int l = levelPtr[curLevel-1]; l < levelPtr[curLevel]; ++l) {
   int cc=levelSet[l];
   for (int j = Lp[cc] + 1; j < Lp[cc + 1]; ++j) {
    inDegree[Li[j]]--;//removing corresponding edges
   }
  }
#if 0
  for (int k = 0; k < n; ++k) {
                std::cout<<inDegree[k]<<",";
            }
            std::cout<<"\n";
#endif

 }
 return curLevel;//return number of levels
}

/*
 *
 */

int buildLevelSet_BCSC(int n, int nnz, size_t *Lp,size_t *Li_ptr, int* Li,
                       size_t blockNo, const int* sup2col, const int* col2sup,
                       int* &levelPtr, size_t* &levelSet ){
 int begin=0,end=blockNo-1;
 int curLevel=0, curLevelCol=0, curCol, nxtCol, supWdt;

 levelPtr = new int[blockNo+1]();
 levelSet = new size_t[blockNo]();
 int *inDegree = new int[blockNo]();
 bool *visited = new bool[blockNo]();
 int *node2Level = new int[blockNo]();
 for (int i = 0; i < blockNo; ++i) {
  inDegree[i]=0; visited[0]= false;
 }
 for (int s = 1; s <= blockNo; ++s) { //In degree computation in assembly tree
  curCol = s!=0 ? sup2col[s-1]:0;
  nxtCol = sup2col[s];
  supWdt = nxtCol-curCol;
  for (int r = Li_ptr[curCol]+supWdt-1; r < Li_ptr[nxtCol]; ++r) {
   inDegree[col2sup[Li[r]]]++;
  }
 }
#if 0
 std::cout<<"\n";
    for (int k = 0; k < blockNo; ++k) {
        std::cout<<inDegree[k]<<",";
    }
    std::cout<<"\n";
#endif
 while(begin <= end){
  for (int i = begin; i <= end; ++i) {//For level curLevel
   if (inDegree[i] == 1 && !visited[i]) {//if no incoming edge
    visited[i] = true;
    assert(i>=0);
    levelSet[curLevelCol] = i; //add it to current level
    node2Level[i] = curLevel;
    curLevelCol++;//Adding to level-set
   }
  }
  curLevel++;//all nodes with 1 indegree are processed.
  assert(curLevelCol>=0);
  assert(curLevel <=blockNo);
  levelPtr[curLevel]=curLevelCol;
  while(inDegree[begin]==1)
   begin++;
  while(inDegree[end]==1 && begin <= end)
   end--;
  //Updating degrees after removing the nodes
  for (int l = levelPtr[curLevel-1]; l < levelPtr[curLevel]; ++l) {
   int cc=levelSet[l]+1;
   curCol = cc!=0 ? sup2col[cc-1]:0;
   nxtCol = sup2col[cc];
   supWdt = nxtCol-curCol;
   for (int j = Li_ptr[curCol]+supWdt; j<Li_ptr[nxtCol]; ++j) {
    inDegree[col2sup[Li[j]]]--;//removing corresponding edges
   }
  }


 }
#if 1
 for (int ii = 0; ii < blockNo; ++ii) {
  assert(visited[ii]== true);
  assert(inDegree[ii] == 1);
 }
 for (int k = 0; k < curLevel; ++k) {
  for (int ii = levelPtr[k]; ii < levelPtr[k+1]; ++ii) {
   int nnooddee = levelSet[ii];
   visited[nnooddee] = false;
  }
 }
 for (int ii = 0; ii < blockNo; ++ii) {
  assert(!visited[ii] );
 }
 for (int ii = 0; ii < curLevel; ++ii) {
  for (int jj = levelPtr[ii]; jj < levelPtr[ii+1]-1; ++jj) {
   int n1 = levelSet[jj];
   int n2 = levelSet[jj+1];
   assert(node2Level[n1] == node2Level[n2]);
  }
 }
#endif
 delete []inDegree;
 delete []visited;
 return curLevel;//return number of levels
}

/*
 *Builds the levelset from the BCSC format
 */
int buildLevelSet_BCSC_fix(size_t n, //FIXME
                       size_t *Lp,
                       size_t *LiP,
                       int *Li,
                       const int *sup2Col,
                       int *&levelPtr,
                       size_t *&levelSet){
 int begin=0,end=n-1;
 int curLevel=0, curLevelCol=0;
 levelPtr = new int[n+1]();
 levelSet = new size_t[n]();
 int *inDegree = new int[n]();
 bool *visited = new bool[n]();
 for (int i = 0; i < LiP[n]; ++i) {//O(nnz)
  inDegree[Li[i]]++;
 }
#if 0
 for (int k = 0; k < n; ++k) {
        std::cout<<inDegree[k]<<",";
    }
    std::cout<<"\n";
#endif
 while(begin <= end){
  for (int i = begin; i <= end; ++i) {//For level curLevel
   if (inDegree[i] == 1 && !visited[i]) {//if no incoming edge
    visited[i] = true;
    levelSet[curLevelCol] = i; //add it to current level
    curLevelCol++;//Adding to level-set
   }
  }
  curLevel++;//all nodes with zero indegree are processed.
  levelPtr[curLevel]=curLevelCol;
  while(inDegree[begin]==1)
   begin++;
  while(inDegree[end]==1 && begin <= end)
   end--;
  //Updating degrees after removing the nodes
  for (int l = levelPtr[curLevel-1]; l < levelPtr[curLevel]; ++l) {
   int cc=levelSet[l];
   for (int j = LiP[cc] + 1; j < LiP[cc + 1]; ++j) {
    inDegree[Li[j]]--;//removing corresponding edges
   }
  }
#if 0
  for (int k = 0; k < n; ++k) {
                std::cout<<inDegree[k]<<",";
            }
            std::cout<<"\n";
#endif

 }
#ifdef VERIFY
 bool *toCheck = new bool[n];
 for (int j = 0; j < n; ++j) {
  toCheck[j]=false;
 }
 for (size_t i = 0; i < n; ++i) {
  for (int j = levelPtr[i]; j < levelPtr[i+1]; ++j) {
   toCheck[levelSet[j]]=true;
  }
 }
 for (int j = 0; j < n; ++j) {
  assert(toCheck[j]);
 }
#endif
 return curLevel;//return number of levels
}
#endif //TRIANGOPENMP_INSPECTION_H
