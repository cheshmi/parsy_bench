//
// Created by kazem on 7/12/18.
//

#ifndef PROJECT_INSPECTIONDAG_03_H
#define PROJECT_INSPECTIONDAG_03_H
#include "../common/TreeUtils.h"
#include "../common/SparseUtils.h"
#include "../triangularSolve/Inspection_Level.h"
#include "../common/BFS.h"



int getCoarseLevelSet_DAG_CSC03(size_t n,
                                int *lC,
                                int* lR,
                                int &finaLevelNo,
                                int* &finaLevelPtr,
                                int* &parLevelSet,
                                int &partNo,
                                int* &finalPartPtr,
                                int* &finalNodePtr,
                                int innerParts,
                                int minLevelDist,
                                int divRate,
                                double *nodeCost){
 int *node2partition = new int[n];
 double *outCost = new double[n];
 double *newOutCost = new double[n];
 int *node2Level = new int[n];
 int *levelPtr; //= new int[n+1]();
 bool *isChanged = new bool[n]();
 bool *visited = new bool[n]();
 int *isMarked = new int[n]();
 int *levelSet; //= new size_t[n]();
 int *xi = new int[2*n];
 int curNumOfPart=0;
 std::vector<int> remainingNodes, remainingTmp;
 std::vector<std::vector<int>> newLeveledParList,mergedLeveledParList;
 int clusterCnt=0;
 int originalHeight=0;
 finaLevelPtr = new int[n+1];
 finalPartPtr = new int[n]();
 finalNodePtr = new int[n];
 int *inDegree = new int[n];
 finaLevelPtr[0]=0;
 for (int i = 0; i < n; ++i) {
  node2partition[i]=-1;
  outCost[i]=0.0;
  newOutCost[i]=0.0;
  inDegree[i]=0;
 }
 int averageCC=0;
 //making levelset
 int levelNo = buildLevelSet_CSC(n,0,lC,lR,levelPtr,levelSet);
 //COMPUTING NODE2lEVEL
 for (int i = 0; i < levelNo; ++i) {
  for (int j = levelPtr[i]; j < levelPtr[i+1]; ++j) {
   int node=levelSet[j];
   node2Level[node]=i;
  }
 }
#if 0
 for (int i = 0; i < levelNo; ++i) {
  std::cout<<i<<"::";
  for (int j = levelPtr[i]; j < levelPtr[i+1]; ++j) {
   std::cout<<levelSet[j]<<";";
  }
  std::cout<<"\n";
 }
#endif
 //Filling indegree array with incoming edges
 for (int s = 0; s < n; ++s) {
  for (int r = lC[s];  r < lC[s+1]; ++r) {
   int cn=lR[r];
   inDegree[cn]++;
  }
 }
 //H-partitioning
 int *partition2Level = new int[levelNo+1]();
 std::vector<int> innerPartsSize;
 originalHeight=levelNo;
 std::vector<std::vector<int>> slackGroups(originalHeight+1);
 std::vector<std::vector<int>> slackedLevelSet(originalHeight+1);
 int lClusterCnt = heightPartitioning_DAG_Trng(levelNo, levelPtr,
                                      NULL, originalHeight, innerParts,
                                      minLevelDist, divRate,
                                      innerPartsSize, slackGroups,
                                      NULL, partition2Level);

 for (int l = 0; l < lClusterCnt; ++l) {//for each leveled partition
  int lbLevel = partition2Level[l]-1;
  int ubLevel = partition2Level[l+1];
  int dfsLevel = partition2Level[l];
  int curLeveledParCost=0;
  //Marking lower bound
  // FIXME: we might need to do all levels below for general DAG
  for (int j = levelPtr[lbLevel>0?lbLevel:0];
       j < levelPtr[lbLevel + 1]; ++j) {
   int curNode = levelSet[j];
   isMarked[curNode]=true;
  }
  //Marking upper bound
  for (int ii = ubLevel; ii < originalHeight ; ++ii) {
   for (int j = levelPtr[ii]; j < levelPtr[ii+1]; ++j) {
    int curNode = levelSet[j];
    isMarked[curNode]=true;
   }
  }
  //Iterating over non-visited leaf nodes to compute CCs
  //CC:connected component
  int stackStart=0, cc=0;
  std::vector<int> needAliased;
  bool *isUniq = new bool[n]();
  for (int k = levelPtr[dfsLevel]; k < levelPtr[dfsLevel+1]; ++k) {
   int curLeaf = levelSet[k];
   bool isCC=true;
   int minAliasedPar=INT32_MAX;
   if(!isMarked[curLeaf]){
    stackStart = dfs_CSC_CC(n, curLeaf,lC,lR,
                            isMarked,
                            n,xi,xi+n,needAliased,NULL);
    //Finding unique clusters from needAliased
    minAliasedPar=makeUnique(node2partition,needAliased,n,isUniq);

    isCC=needAliased.size()==0;
    if(!isCC){//There are some intersection between found CCs.
     for (int j = 0; j < needAliased.size(); ++j) {//the first is min
      int tn=node2partition[needAliased[j]];
      if(tn!=minAliasedPar){
       cc--;
       for (int i = 0; i < n; ++i) {
        //Replace all needAliased node with their min part number.
        if(node2partition[i]==tn){
         node2partition[i]=minAliasedPar;
        }
       }
      }
     }
     needAliased.erase(needAliased.begin(),needAliased.end());
     //Set the nodes in the current cluster
     for (int i = stackStart; i < n; ++i) {
      int node=xi[i];
      node2partition[node]=minAliasedPar;
      //Compute the cost of each CC
      outCost[minAliasedPar]+=nodeCost[node];
      //The cost of cur h-partition
      curLeveledParCost+=nodeCost[node];
      //reseting all nodes but leaf node
      //marke it with -1
      if(node2Level[node] != dfsLevel)
       isMarked[node]=-1;
     }
    }else{
     for (int i = stackStart; i < n; ++i) {
      int node=xi[i];
      node2partition[node]=cc;
      //Compute the cost of each CC
      outCost[cc]+=nodeCost[node];
      curLeveledParCost+=nodeCost[node];
      //reseting all nodes but leaf node
      //marke it with -1
      if(node2Level[node] != dfsLevel)
       isMarked[node]=-1;
     }
     cc++;// one more CC.
    }
   }
  }
  //Reset all marked node in the DAG
    //std::cout<<cc<<"\n";
  for (int j = levelPtr[lbLevel>0?lbLevel:0];
       j < levelPtr[lbLevel + 1]; ++j) {
   int curNode = levelSet[j];
   isMarked[curNode]= false;
  }
  //Marking upper bound
  for (int ii = ubLevel; ii < originalHeight ; ++ii) {
   for (int j = levelPtr[ii]; j < levelPtr[ii+1]; ++j) {
    int curNode = levelSet[j];
    isMarked[curNode]= false;
    visited[curNode] = true;//Make it ready for mod-BFS
   }
  }
  //Topological sort of each cc, the fastest way, FIXME: make it more
  //local
  std::vector<int> extraDim;
  for (int i = 0; i < cc; ++i) {
   newLeveledParList.push_back(extraDim);
  }
  modifiedBFS_CSC(n,lC,lR,inDegree,visited,
                  node2partition,levelPtr,levelSet,dfsLevel,
                  newLeveledParList);
  /*for (int ll = dfsLevel; ll < ubLevel; ++ll) {
   for (int ii = levelPtr[ll]; ii < levelPtr[ll+1]; ++ii) {
    int curNode=levelSet[ii];
    assert(node2partition[curNode]>=0);
    newLeveledParList[node2partition[curNode]].push_back(curNode);
   }
  }*/
  //Marking upper bound
  for (int ii = ubLevel; ii < originalHeight ; ++ii) {
   for (int j = levelPtr[ii]; j < levelPtr[ii+1]; ++j) {
    int curNode = levelSet[j];
    visited[curNode] = false;//Make it ready for mod-BFS
   }
  }
  //Bin packing and form W-partitions
  int levelParCostThresh=curLeveledParCost/innerParts;
  levelParCostThresh+=(0.1*levelParCostThresh);
  int outinnerParts=0;
  averageCC+=newLeveledParList.size();
  mergedLeveledParList.resize(innerPartsSize[l]);//FIXME
  if(newLeveledParList.size() > innerPartsSize[l]) {
   outinnerParts= worstFitBinPack(newLeveledParList, outCost,
                                  mergedLeveledParList, newOutCost,
                                  levelParCostThresh, innerPartsSize[l]);
   //assert(outinnerParts<=innerParts);
  }else {
   mergedLeveledParList.erase(mergedLeveledParList.begin(),mergedLeveledParList.end());
   mergedLeveledParList=newLeveledParList;
   outinnerParts=newLeveledParList.size();
#if 0
   if(outinnerParts>1) {
    for (int ii = 0; ii < outinnerParts; ++ii) {
     std::cout << outCost[ii] << ";";
    }
    for (int ii = outinnerParts; ii < innerParts; ++ii) {
     std::cout << "0;";
    }
   }
#endif
  }

  double curPartCost=0;
  finaLevelPtr[l+1] = finaLevelPtr[l]+outinnerParts;

  for (int i = 0; i < outinnerParts; ++i) {
   int curPartElem=0;
   curPartCost=0;
   for (int j = 0; j < mergedLeveledParList[i].size(); ++j) {
    curPartCost += nodeCost[mergedLeveledParList[i][j]];
    finalNodePtr[finalPartPtr[curNumOfPart]+curPartElem]=mergedLeveledParList[i][j];
    node2partition[mergedLeveledParList[i][j]] = curNumOfPart;
    curPartElem++;
   }
#if 0
   //std::cout<<"parts: "<<newLeveledParList.size()<<","<<curPartCost<<", ";
   std::cout<<curPartCost<<", ";
#endif
   finalPartPtr[curNumOfPart+1] = finalPartPtr[curNumOfPart] + curPartElem;
   curNumOfPart++;
  }

  //Cleaning the current sets.
  for (int i = 0; i < mergedLeveledParList.size(); ++i) {
   mergedLeveledParList[i].erase(mergedLeveledParList[i].begin(),mergedLeveledParList[i].end());
  }
  mergedLeveledParList.erase(mergedLeveledParList.begin(),mergedLeveledParList.end());
  for (int i = 0; i < newLeveledParList.size(); ++i) {
   newLeveledParList[i].erase(newLeveledParList[i].begin(),newLeveledParList[i].end());
  }
  newLeveledParList.erase(newLeveledParList.begin(),newLeveledParList.end());

 }

 finaLevelNo=lClusterCnt;
 if(true) {//Verification of the set.
  bool *checkExist = new bool[n];
  for (int i = 0; i < n; ++i) checkExist[i] = false;
  for (int i = 0; i < lClusterCnt; ++i) {
   for (int k = finaLevelPtr[i]; k < finaLevelPtr[i + 1]; ++k) {
    for (int j = finalPartPtr[k]; j < finalPartPtr[k + 1]; ++j) {
     assert(checkExist[finalNodePtr[j]] == false);
     checkExist[finalNodePtr[j]] = true;
    }
   }
  }
  for (int i = 0; i < n; ++i) {
   assert(checkExist[i] == true);
  }
  delete[] checkExist;
 }
 delete [] outCost;
 delete [] newOutCost;
 delete [] partition2Level;
 delete [] levelPtr;
 delete [] levelSet;
 delete [] node2partition;
 delete [] node2Level;
 delete [] isChanged;
 delete [] isMarked;
 delete []xi;
 delete []visited;
 delete []inDegree;

 return averageCC/lClusterCnt;

}
#endif //PROJECT_INSPECTIONDAG_03_H
