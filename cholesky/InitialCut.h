//
// Created by kazem on 7/6/18.
//

#ifndef PROJECT_INITIALCUT_H
#define PROJECT_INITIALCUT_H

#include <iostream>
#include <Inspection_Level.h>
#include <vector>

int* getInitialCut_DAG(size_t *lC, size_t *Li_ptr, int* lR,
                       size_t n, const int *blk2col, const int *col2blk
){

 int *node2partition = new int[n];
 int *node2Level = new int[n];
 int *levelPtr = new int[n+1]();
 bool *isChanged = new bool[n]();
 bool *isMerged = new bool[n]();
 size_t *levelSet = new size_t[n]();
 std::vector<int> remainingNodes, remainingTmp;
 int clusterCnt=0;
 for (int i = 0; i < n; ++i) {
  node2partition[i]=-1;
 }
 int levelNo = buildLevelSet_BCSC(0,0,lC,Li_ptr,lR,
                                  n,blk2col,col2blk,levelPtr, levelSet);
 //COMPUTING NODE2lEVEL
 for (int i = 0; i < levelNo; ++i) {
  for (int j = levelPtr[i]; j < levelPtr[i+1]; ++j) {
   int node=levelSet[j];
   node2Level[node]=i;
  }
 }
 int *connectedCompPerLevel = new int[levelNo]();

 for (int i = levelPtr[0]; i < levelPtr[1]; ++i) {
  int node=levelSet[i];
  node2partition[node] = clusterCnt++;
 }
 connectedCompPerLevel[0]=clusterCnt;
 //connectedCompPerLevel[levelNo-1]=1;
 //TODO: assume that the max number of partitions are in leaf nodes.
 bool *clusterCount = new bool[clusterCnt]();
 for (int m = 0; m < levelNo; ++m) {
  for (int i = levelPtr[m]; i < levelPtr[m+1]; ++i) {
   std::cout<<levelSet[i]<<",";
  }
  std::cout<<"\n";
 }
 int supWdt=0;
 for (int i = 0; i < levelNo-1; ++i) {

  for (int j = levelPtr[i]; j < levelPtr[i+1]; ++j) {
   int curNode=levelSet[j];
   int curCol=blk2col[curNode];
   int nxtCol=blk2col[curNode+1];
   supWdt = nxtCol-curCol;
   int jj = Li_ptr[curCol]+supWdt;
   //looking into curNode neighbourhood
   for (; jj < Li_ptr[nxtCol]; ) {
    int k = col2blk[lR[jj]];
    supWdt = blk2col[k+1]-blk2col[k];
    jj+=supWdt;
    if(node2Level[k] == i+1){
     isMerged[curNode] = true;//We know the node is merged now.
     if(node2partition[k]<0){
      node2partition[k] = node2partition[curNode];
     }else{
      //Let's rename the cluster of curNode
      int curDestCluster = node2partition[k];
      int curSrcCluster = node2partition[curNode];
      if(!isChanged[curNode]){//First time merging
       isChanged[curNode]=true;
       node2partition[curNode] = curDestCluster;
       //node2partition[k]=curSrcCluster;
       //TODO should I update the other nodes in this level?
      }else{//Then this is the second renaming, has to update other nodes
       node2partition[curNode] = curDestCluster;
       //node2partition[k]=curSrcCluster;
       for (int l = levelPtr[i]; l < j; ++l) {
        //Make the remaining node consistent
        if(node2partition[l] == curSrcCluster){
         node2partition[l]=curDestCluster;
        }
       }
      }
     }
    }
   }// one is processed

   if(!isMerged[curNode]){
    remainingTmp.push_back(curNode);
   }
  }// Processing one level is finished
  //Let's process the remaining nodes from last levels
  int processedNodes=remainingNodes.size();
  for (int m = 0; remainingNodes.size()>0 && processedNodes>0;
       processedNodes--) {
   int curNode=remainingNodes[m];
   if(!isMerged[curNode]){
    int curCol=blk2col[curNode];
    int nxtCol=blk2col[curNode+1];
    supWdt = nxtCol-curCol;
    int jj = Li_ptr[curCol]+supWdt;
    //looking into curNode neighbourhood
    bool removed = false;
    for (; jj < Li_ptr[nxtCol]; ) {
     int k = col2blk[lR[jj]];
     supWdt = blk2col[k + 1] - blk2col[k];
     jj += supWdt;
     if (node2Level[k] == i + 1) {
      isMerged[curNode] = true;
      removed = true;
      if (node2partition[k] < 0) {
       node2partition[k] = node2partition[curNode];
      } else {
       //Let's rename the cluster of curNode
       int curDestCluster = node2partition[k];
       int curSrcCluster = node2partition[curNode];
       if (!isChanged[curNode]) {//First time merging
        isChanged[curNode] = true;
        node2partition[curNode] = curDestCluster;
        // node2partition[k] = curSrcCluster;
        //TODO should I update the other nodes in this level?
       } else {//Then this is the second renaming, has to update other nodes
        node2partition[curNode] = curDestCluster;
        // node2partition[k] = curSrcCluster;
        //Make the remaining node consistent
        for (int ll = levelPtr[i]; ll < levelPtr[i+1]; ++ll) {
         int l=levelSet[ll];
         if (node2partition[l] == curSrcCluster) {
          node2partition[l] = curDestCluster;
         }
        }
        //Make labels consistent in other remaining nodes
        for (int ll = 0; ll < remainingNodes.size(); ++ll) {
         int l = remainingNodes[ll];
         if (node2partition[l] == curSrcCluster) {
          node2partition[l] = curDestCluster;
         }
        }
       }
      }
     }
    }
    if(removed){//the node is processed
     remainingNodes.erase(remainingNodes.begin()+m);
    }else{//then, too early to remove this node, go to next one
     m++;
    }

   }else{//so, the node is just merged, no more processing, remove it!
    remainingNodes.erase(remainingNodes.begin()+m);
   }
  }
  remainingNodes.insert(remainingNodes.end(),remainingTmp.begin(),
                        remainingTmp.end());
  for (int kkk = levelPtr[i]; kkk < levelPtr[i+1]; ++kkk) {
   int kk=levelSet[kkk];
   int par = node2partition[kk];
   if(par>=0 && !clusterCount[par]){
    clusterCount[par]=true;
    connectedCompPerLevel[i+1]++;
   }
  }
  for (int kk = 0; kk < clusterCnt; ++kk) {
   clusterCount[kk]=false;
  }

  for (int i1 = 0; i1 < n; ++i1) {
   std::cout<<node2partition[i1]<<";";
  }
  std::cout<<"\n";
 }
 for (int k1 = 0; k1 < levelNo; ++k1) {
  std::cout<<connectedCompPerLevel[k1]<<";";
 }
 std::cout<<"\n";
 return connectedCompPerLevel;
}
#endif //PROJECT_INITIALCUT_H
