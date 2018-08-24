//
// Created by kazem on 8/21/17.
//

#include <iostream>
#include <vector>

#ifndef CHOLOPENMP_BUILDLEVELSET_H
#define CHOLOPENMP_BUILDLEVELSET_H






void copyChild2Partition(const int *inActPtr, const int *inActNode,
                         const int &outSize, const int *outChildPtr,
                         const int *outChildNo, int *outActPtr,
                         int *outActNode) {
 int cnt=0;
 for (int i = outChildPtr[outSize]; i < outChildPtr[outSize + 1]; ++i) {
  int curParChil = outChildNo[i];
  for (int j = inActPtr[curParChil]; j < inActPtr[curParChil+1]; ++j) {
   outActNode[outActPtr[outSize]+cnt]=inActNode[j];
   cnt++;
  }
 }
 outActPtr[outSize+1]=outActPtr[outSize]+cnt;
}

void copyChild2PartitionAll(const int *inChildPtr, const int *inChildNo,
                            const int &outSize, const int *outParChildPtr,
                            const int *outParChildNo, int* outParChildPtrTmp,
                            int *outParChildNoTmp, int *outChildPtr,
                         int *outChildNo) {
 for (int h = 0; h < outSize; ++h) {//for each node in par ETree
  int cnt=0;
  //for each partitioned child
  /*for (int i = outParChildPtr[h]; i < outParChildPtr[h + 1]; ++i) {
   int curParChil = outParChildNo[i];*/
   //for each node in actual node list
   for (int k = outParChildPtrTmp[h];k < outParChildPtrTmp[h+1]; ++k) {
    int currentChild = outParChildNoTmp[k];
    //Now take the actual node from the initial partition.
    for (int j = inChildPtr[currentChild]; j < inChildPtr[currentChild+1]; ++j){
     outChildNo[outChildPtr[h]+cnt]=inChildNo[j];
     cnt++;
    }
   }
   outChildPtr[h+1]=outChildPtr[h]+cnt;
  //}
 }//end for h
}



int coarsening(int inSize, int *inTree, const int *origETree, int *inCost,
               int *nChild,
               int *inParChildPtr, int *inParChildNo, //Partitioned children
               int *inChildPtr, int *inChildNo,//Actual nodes in partitions
               int *inNod2Par, int Threshold, int mergeFactor,
               /*Outputs*/
               int &outSize, int* outTree, int* outCost,
               int *outParChildPtr, int *outParChildNo,
               int* outChildPtr, int* outChildNo,
               int* outNode2Par){
 auto *visited = new bool[inSize]();
 int tmpN=0;
 outChildPtr[0]=0;
 outChildNo[0]=0;
 int *outParChildPtrTmp=new int[inSize](), *outParChildNoTmp=new int[inSize]();
 outParChildPtrTmp[0]=0;
 //std::fill_n(outCost,inSize,0);
 for (int i = 0; i < inSize; ++i)
  outCost[i]=0;
 //visited[0]=true;
 int curPartElem=0;
 outSize=0;
 for (int k = 0; k < inSize; ++k){//k is in partitioned node
  if(outCost[outSize] < Threshold){
   //Add current node if not visited before
   if(!visited[k]){
    if(nChild[k]>0){//If it has any unvisited children, add them first.
     for (int i = inParChildPtr[k]; i < inParChildPtr[k+1]; ++i) {
      int cn = inParChildNo[i];
      //int cnp = inNod2Par[cn];
      if(!visited[cn] && nChild[cn]==0){
       outParChildNoTmp[curPartElem++]=cn;
       outCost[outSize]+=inCost[cn];
       visited[cn]=true;
       outNode2Par[cn]=outSize;//FIXME
       //Making sure parent is in the same outParChildNoTmp with its children
       nChild[k]--;
      }
     }
    }
    //tmpN=inNod2Par[k];
    outParChildNoTmp[curPartElem++]=k;//Add the node itself
    outCost[outSize]+=inCost[k];
    visited[k]=true;
    outNode2Par[k]=outSize;//FIXME
    if(inTree[k]!=-1)nChild[inTree[k]]--;
    if(outCost[outSize]>Threshold){//close the outParChildNoTmp and start the new one
     outParChildPtrTmp[outSize+1]=curPartElem;
     //Copy actual children from input to the output children list
     /*copyChild2Partition(inParChildPtr, inParChildNo, outSize, outParChildPtrTmp, outParChildNoTmp,
                         outParChildPtrTmp, outParChildNoTmp);*/
     outSize++;
     continue;
    }
    int par=k;
    int mfactor = mergeFactor;
    //going up till making the outParChildNoTmp large enough
    while(outCost[outSize] < Threshold
          && inTree[par] != -1
          && mfactor>0){
     par= inTree[par];
     if(!visited[par]){
      mfactor--;
      if(nChild[par]>0){//If parent has any children, add them to outParChildNoTmp
       for (int i = inParChildPtr[par]; i < inParChildPtr[par+1]; ++i) {
        int cn = inParChildNo[i];//real node
        //int cnp = inNod2Par[cn];
        if(!visited[cn]) {
         if (nChild[cn] == 0) {
          outParChildNoTmp[curPartElem++] = cn;
          outCost[outSize] += inCost[cn];
          visited[cn] = true;
          outNode2Par[cn] = outSize;//FIXME
          nChild[par]--;
         } else {
          for (int ii = inParChildPtr[cn]; ii < inParChildPtr[cn + 1]; ++ii) {
           int cnc = inParChildNo[ii];
           //int cncp = inNod2Par[cnc];
           if (!visited[cnc] && nChild[cnc]==0){
            outParChildNoTmp[curPartElem++] = cnc;
            outCost[outSize] += inCost[cnc];
            visited[cnc] = true;
            outNode2Par[cnc] = outSize;//FIXME
            nChild[cn]--;//the father of cnc
           }
          }
          outParChildNoTmp[curPartElem++] = cn;//Add parent
          outCost[outSize] += inCost[cn];
          visited[cn] = true;
          outNode2Par[cn] = outSize;//FIXME
          nChild[par]--;//par is the child of cn
         }
        }//if node is not visited
       }//End for, checking parent of current node
      }//if parent has any children
      outParChildNoTmp[curPartElem++]=par;//Add parent
      outCost[outSize]+=inCost[par];
      visited[par]=true;
      outNode2Par[par]=outSize;//FIXME
      if(inTree[par]!=-1)nChild[inTree[par]]--;
     }
    }
    //Close the outParChildNoTmp anyway
    if(outParChildPtrTmp[outSize] != curPartElem){//If the partion is not empty
     outParChildPtrTmp[outSize+1]=curPartElem;
     //Copy actual children from input to the output children list
     /*copyChild2Partition(inParChildPtr, inParChildNo, outSize, outParChildPtrTmp, outParChildNoTmp,
                         outParChildPtrTmp, outParChildNoTmp);*/
     outSize++;//This is the maximum size of this outParChildNoTmp
    }
   }

  }else{//close this outParChildNoTmp and go to the next one
   outParChildPtrTmp[outSize+1]=curPartElem;
   //Copy actual children from input to the output children list
   /*copyChild2Partition(inParChildPtr, inParChildNo, outSize, outParChildPtrTmp, outParChildNoTmp,
                       outParChildPtrTmp, outParChildNoTmp);*/
   outSize++;
  }
 }

#if 0
 std::cout<<"partitions: "<<outSize<<"\n";
 for (int i1 = 0; i1 < outSize + 1; ++i1) {
  std::cout<<outParChildPtrTmp[i1]<<",";
 }
 std::cout<<"\n";
#endif

 //Creating partitioned ETree
 for (int i = 0; i < outSize; ++i) {
//  assert(outParChildPtrTmp[i+1]-1 < inSize);
  //the parent of the parent in the outParChildNoTmp i
  //int cp = outParChildNoTmp[outParChildPtrTmp[i+1]-1];
  int cp = outParChildNoTmp[outParChildPtrTmp[i+1]-1];
 // assert(outNode2Par[inTree[cp]] < inSize);
  outTree[i] = inTree[cp]>0 ? outNode2Par[inTree[cp]] : -1;
 }
#if 0
 std::cout<<"par ETRee: "<<outSize<<"\n";
 for (int i1 = 0; i1 < outSize ; ++i1) {
  std::cout<<outTree[i1]<<",";
 }
 std::cout<<"\n Actual Node partition\n";
 for (int j = 0; j < outSize; ++j) {
  for (int i = outParChildPtrTmp[j]; i < outParChildPtrTmp[j + 1]; ++i) {
   std::cout<<outParChildNoTmp[i]<<",";
  }
  std::cout<<"\n";
 }
 std::cout<<"\n";
#endif

 //Populate actual children using parChild
 //populateChildren(outSize,outTree,outParChildPtr,outParChildNo); //FIXME add nChild param
#if 0
 std::cout<<"\n partiton ETREE\n";
 for (int j = 0; j < outSize; ++j) {
  for (int i = outParChildPtr[j]; i < outParChildPtr[j + 1]; ++i) {
   std::cout<<outParChildNo[i]<<",";
  }
  std::cout<<"\n";
 }
 std::cout<<"\n";
#endif

 copyChild2PartitionAll(inChildPtr,inChildNo,
                        outSize,outParChildPtr,outParChildNo,
                        outParChildPtrTmp,outParChildNoTmp,
                        outChildPtr,outChildNo);
#if 0
 std::cout<<"\n Actual nodes overall:\n";
 for (int j = 0; j < outSize; ++j) {
  for (int i = outChildPtr[j]; i < outChildPtr[j + 1]; ++i) {
   std::cout<<outChildNo[i]<<",";
  }
  std::cout<<"\n";
 }
 std::cout<<"\n";
#endif
 delete []visited;
 return outSize;
}



#endif //CHOLOPENMP_BUILDLEVELSET_H
