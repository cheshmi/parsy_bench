//
// Created by kazem on 10/17/17.
//

#ifndef CHOLOPENMP_MERGEWITHSLACK_H
#define CHOLOPENMP_MERGEWITHSLACK_H
#include <vector>
#include <iostream>
#include "../common/TreeUtils.h"
#include "Partitioning.h"


int mergingWithTightSlack(int n, int *eTree, int *nodeCost, int *parCost, int *Node2Par, int *nChild,
                          std::vector<std::vector<int>> parList, int initiLevelPartition, int innerParts,
                          int minLevelsforMerge,
        //Outputs
                          int *outTree, std::vector<std::vector<int>> newParList,
                          int &finaLevelNo, int *finaLevelPtr,
                          int *finalPartPtr, int *finalNodePtr
){
    int partitionNo = parList.size();
    int overallCost=0;
    auto nodeHeight = new int[n]();
    auto slackNumber = new int[n]();
    std::vector<int> criticalLeaves;
    auto nChildArtTree = new int[partitionNo]();
    auto partitionMaxLevels = new int[partitionNo]();
    auto partitionMinLevels = new int[partitionNo]();
    auto initLevelPtr = new int[n+1]();
    auto initLevelSet = new int[n]();
    auto levelPtr = new int[n+1]();
    auto levelSet = new int[n]();
    auto node2Level = new int[n]();
    auto childPtrSubtree=new int[n+1](),
            childNoSubtree=new int[n](),
            nChildTmp=new int[n]();
    double fracThr=0.2;
    int curNumOfPart=0;
    finaLevelPtr[0]=0;
    std::vector<std::vector<int>> newLeveledParList, mergedLeveledParList;
    for (int i = 0; i < partitionNo; ++i) {
        partitionMinLevels[i]=INT_MAX;
    }
    for (int i = 0; i < n; ++i) {
        nChildTmp[i]=nChild[i];
    }
    //Find the number of connected parts
    int treeParts=0, actualArtiHeight=0;
    for (int j = 0; j < n; ++j) {
        if(eTree[j]<0)
            treeParts++;
    }
    //Finding node levels in original ETree
    //TODO also find the unweighted version of eTree, there might be more levels ignoring them.
    int originalHeight = getTreeHeight(n, eTree, nChildTmp);
    for (int i = 0; i < n; ++i) {
        nodeHeight[i]= getNodeDepth(i, n, eTree);
        if(nodeHeight[i]==originalHeight)
            criticalLeaves.push_back(i);
        assert(nodeHeight[i]>=0);
    }
    //Computing slack number for each node
    int cntSlack=0,costSlack=0;
    std::vector<std::vector<int>> slackGroups(originalHeight+1);
    //First compute node2level
    int levelNo = getLevelSet(n,eTree,levelPtr,levelSet);
    for (int i = 0; i < levelNo; ++i) {
        for (int j = levelPtr[i]; j < levelPtr[i + 1]; ++j) {
            node2Level[levelSet[j]]=i;
        }
    }
    int maxLevelWithSlack=0;
    for (int i = 0; i < n; ++i) {
        slackNumber[i] = originalHeight-nodeHeight[i]-node2Level[i];
        assert(slackNumber[i]<originalHeight);
        //Postpone the node i, to higher levels
        if(node2Level[i]+slackNumber[i]>maxLevelWithSlack)
            maxLevelWithSlack=node2Level[i]+slackNumber[i];
        slackGroups[node2Level[i]+slackNumber[i]].push_back(i);
        if(slackNumber[i]>22 && (node2Level[i]==0 || node2Level[i]==1)){
            cntSlack++;
            costSlack+=nodeCost[i];
        }
        //slackGroups[slackNumber[i]].push_back(i);
    }
#if 0
    for (int i = 0; i < levelNo; ++i) {
        std::cout<<"l# "<<i<<": ";
        for (int j = levelPtr[i]; j < levelPtr[i+1]; ++j) {
            std::cout<<levelSet[j]<<";";
        }
        std::cout<<"\n";
    }
    for (int i = 0; i < slackGroups.size(); ++i) {
        std::cout<<"s# "<<i<<": ";
        for (int j = 0; j < slackGroups[i].size(); ++j) {
            std::cout<<slackGroups[i][j]<<";";
        }
        std::cout<<"\n";
    }
#endif
#if 1
    std::cout<<"# of nodes with slack: "<<cntSlack<<"\n";
    std::cout<<" slack cost: "<<costSlack<<"\n";
#endif
    //Creating new level-set after Slack number applied.
    int curElement=0;
    levelPtr[0]=0;
    for (int i = 0; i < slackGroups.size(); ++i) {
        for (int j = 0; j < slackGroups[i].size(); ++j) {
            node2Level[slackGroups[i][j]]=i;
            levelSet[curElement++]=slackGroups[i][j];
        }
        levelPtr[i+1]=curElement;
    }
#if 0
    for (int i = 0; i < levelNo; ++i) {
        std::cout<<"lNew# "<<i<<": ";
        for (int j = levelPtr[i]; j < levelPtr[i+1]; ++j) {
            std::cout<<levelSet[j]<<";";
        }
        std::cout<<"\n";
    }
#endif
    //computing the cost of each level,
    auto levelCost = new int[levelNo]();
    auto levelParCost = new int[levelNo]();
    for (int i = 0; i < n; ++i) {
        assert(node2Level[i]<levelNo);
        levelCost[node2Level[i]] += nodeCost[i];
        overallCost+=nodeCost[i];
    }
    //Compute the average cost of each level
    int sumOfSlackedLevel=0;
    for (int i = maxLevelWithSlack; i > 0; --i) {
        sumOfSlackedLevel+=levelCost[i];
    }
    int costPerLevelPartition = overallCost/initiLevelPartition;
    //TODO: reordering the levels and shift all nodes with slack number'
    //Merging levels till get to 'costPerLevelPartition' number
    auto partition2Level = new int[levelNo+1]();
    partition2Level[0]=0;
    int lClusterCnt=0;
    for (int i = 0; i < levelNo; ) {
        int curCost=0;
        while (curCost<costPerLevelPartition && i<levelNo){
            curCost+=levelCost[i];
            i++;
        }
        if(i+minLevelsforMerge>levelNo){//last nodes
            levelParCost[lClusterCnt]=curCost;
            lClusterCnt++;
            partition2Level[lClusterCnt] = i ;
            break;
        }
        //if it is too large
        if(curCost>costPerLevelPartition+(fracThr*costPerLevelPartition) &&
           i-minLevelsforMerge>partition2Level[lClusterCnt] ){
            i--;
        }
        if(i-partition2Level[lClusterCnt]>minLevelsforMerge-1){//At least two levels
            levelParCost[lClusterCnt]=curCost;
            lClusterCnt++;
            partition2Level[lClusterCnt] = i ;
        }
    }
    partition2Level[lClusterCnt] = originalHeight+1;//FIXME why +1
#if 1
    for (int i = 0; i < levelNo; ++i) {
        std::cout<<levelCost[i]<<";";
    }
    std::cout<<"\n";
    for (int i = 0; i < lClusterCnt+1; ++i) {
        std::cout<<partition2Level[i]<<";";
    }
    std::cout<<"\n";
    std::cout<<cntSlack<<","<<costPerLevelPartition<<"\n";
#endif

    //TODO:  removing the nodes that have many edge cuts

    //TODO: partition and reorder the nodes within each level partition
    //TODO: First make the partition that has nodes with slack zero
    //Create a tree for each level partition
    auto tmpTree = new int[n];
    auto outCost = new int[n](),
            newOutCost = new int[n](),
            outNode2Par = new int[n]();
    int outSize=0;
    int ccc=0;
    for (int i = 0; i < n; ++i)
        tmpTree[i]=-2;
    for (int l = 0; l < lClusterCnt; ++l) {//for each leveled partition
        int curLeveledParCost=0;
        for (int i = partition2Level[l]; i < partition2Level[l+1]; ++i) {
            curLeveledParCost+=levelCost[i];
            for (int j = levelPtr[i]; j < levelPtr[i + 1]; ++j) {
                int curNode=levelSet[j];
                int par=eTree[curNode];
                /*if(curNode==25 || par==25)
                    std::cout<<"\n";*/
                if(node2Level[par] < partition2Level[l+1] &&
                   l+1 < lClusterCnt &&
                   node2Level[par] >= partition2Level[l]){//The parent is in the same partition
                    tmpTree[curNode]=eTree[curNode];//copy only those in the partition
                    ccc+=nodeCost[curNode];
                }else{
                    tmpTree[curNode]=-1; //The edge-cut happens here
                }
            }
        }
        //Here we have the ETree for the lth leveled partition

        //Now compute the number of children for the etree of current partition
        for (int i = 0; i < n; ++i) nChildTmp[i]=0;
        populateChildren(n,tmpTree,childPtrSubtree,childNoSubtree,nChildTmp);
#if 0
        for (int m = 0; m < n; ++m) {
            std::cout<<nChildTmp[m]<<";";
        }
        std::cout<<"\n";
#endif
        //partition it using the partitioning function call
        outSize=0;
        for (int k = 0; k < n; ++k) {
            outNode2Par[k]=0;outCost[k]=0;
        }
        treeParts=0;
        for (int j = 0; j < n; ++j) {
            if(tmpTree[j]==-1)
                treeParts++;
        }
#if 0
        for (int i = 0; i < n; ++i) {
            std::cout<<tmpTree[i]<<";";
        }
        std::cout<<"\n";
#endif
        int outinnerParts=0;
        int levelParCostThresh=curLeveledParCost/innerParts;
        levelParCostThresh+=(0.1*levelParCostThresh);
        partitioning(n,tmpTree,nodeCost,childPtrSubtree,childNoSubtree,
                     nChildTmp,curLeveledParCost,innerParts,
                     outSize,outCost,outNode2Par,newLeveledParList);
        int tmpp=0;
        for (int i = 0; i < newLeveledParList.size(); ++i) {
            for (int j = 0; j < newLeveledParList[i].size(); ++j) {
                tmpp += nodeCost[newLeveledParList[i][j]];
            }
        }
        assert(tmpp == curLeveledParCost);
        //the partitioning is found, now we need merge some of those if the number is larger than
        //inner parts. It should be typically bigger since the input is a forest.
        // merging the final partitioning

        mergedLeveledParList.resize(innerParts);
        if(newLeveledParList.size() > innerParts) {
            outinnerParts=mergeInnerPart(newLeveledParList,outCost,mergedLeveledParList,
                                         newOutCost,levelParCostThresh);
            assert(outinnerParts<=innerParts);
        }else{
            mergedLeveledParList.erase(mergedLeveledParList.begin(),mergedLeveledParList.end());
            mergedLeveledParList=newLeveledParList;
            outinnerParts=newLeveledParList.size();
        }


        //resetting temporary ETree
        for (int j = 0; j < n; ++j) {
            tmpTree[j]=-2;
        }
        //putting the inner partitions into the final level set
        // int *finaLevelPtr, int *finalPartPtr, int *finalNodePtr
        finaLevelPtr[l+1] = finaLevelPtr[l]+outinnerParts;
        for (int i = 0; i < outinnerParts; ++i) {
            int curPartElem=0;
            for (int j = 0; j < mergedLeveledParList[i].size(); ++j) {
                finalNodePtr[finalPartPtr[curNumOfPart]+curPartElem]=mergedLeveledParList[i][j];
                outNode2Par[mergedLeveledParList[i][j]] = curNumOfPart;
                curPartElem++;
            }
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
#if 0
        int cc=0;
        for (int m = 0; m < mergedLeveledParList.size(); ++m) {
            std::cout<<"Partition# "<<m<<": ";
            for (int n = 0; n < mergedLeveledParList[m].size(); ++n) {
                std::cout<<mergedLeveledParList[m][n]<<";";
                cc+=nodeCost[mergedLeveledParList[m][n]];
            }
            std::cout<<"\n";
        }
        std::cout<<"\n"<<cc<<"\n";
#endif
    }
    std::cout<<"parts: "<<curNumOfPart<<"\n";
#if OLD
    //Making partitioned ETree
    for (int i = 0; i < parList.size(); ++i) {
        int endEl = parList[i].size()-1;//Last element is the highest level
        int par = eTree[parList[i][endEl]];
        outTree[i]= par>0 ? Node2Par[par] : -1;
        assert(outTree[i]<partitionNo);
    }
    //Finding the number of children in artETree
    for (int i = 0; i < partitionNo; ++i)
        if(outTree[i]>=0)
            nChildArtTree[outTree[i]]++;

    //Finding the max and min level in each partition
    for (int i = 0; i < partitionNo; ++i) {
        for (int j = 0; j < parList[i].size(); ++j) {
            int tmp=nodeHeight[parList[i][j]];
            if(tmp > partitionMaxLevels[i]){
                partitionMaxLevels[i]=tmp;
            }
            if(tmp < partitionMinLevels[i]){
                partitionMinLevels[i]=tmp;
            }
        }
    }
    int artHeight = getTreeHeight(n, outTree, nChildArtTree, parCost);
    //Finding the height of the artificial ETree.
    for (int k = 0; k < partitionNo; ++k) {
        actualArtiHeight += (partitionMaxLevels[k]-partitionMinLevels[k]);
    }

    std::cout<<"The original tree height is: " <<originalHeight<<"\n";
    std::cout<<"The actual artificial tree height is:"<<actualArtiHeight<<"\n";
    std::cout<<"The artificial tree height is:"<<artHeight<<"\n";
    std::cout<< "Tree parts are: "<<treeParts<< "\n";
    for (int i = 0; i < partitionNo; ++i) {
        std::cout<<i<<" :"<<partitionMaxLevels[i]-partitionMinLevels[i]<<"; ";
    }
    std::cout<<"\n";
#endif
    finaLevelNo=lClusterCnt;
#if 1
    bool *checkExist = new bool[n];
    for (int i = 0; i < n; ++i) checkExist[i]=false;
    for (int i = 0; i < lClusterCnt; ++i) {
        for (int k = finaLevelPtr[i]; k < finaLevelPtr[i+1]; ++k) {
            for (int j = finalPartPtr[k]; j < finalPartPtr[k+1]; ++j) {
                assert(checkExist[finalNodePtr[j]]== false);
                int par=eTree[finalNodePtr[j]]>= 0 ? eTree[finalNodePtr[j]]:n-1;
                assert(outNode2Par[finalNodePtr[j]] <= outNode2Par[par]);
                checkExist[finalNodePtr[j]]=true;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        assert(checkExist[i]==true);
    }
    delete []checkExist;
#endif
    delete [] nodeHeight;
    delete [] slackNumber;
    delete []nChildArtTree;
    delete [] partitionMaxLevels;
    delete [] partitionMinLevels;
    delete []tmpTree;
    delete [] outCost;
    delete [] newOutCost;
    delete [] outNode2Par;
    delete [] levelPtr;
    delete [] levelSet;
    delete [] node2Level;
    delete []  childPtrSubtree;
    delete [] childNoSubtree;
    delete [] nChildTmp;
    delete [] levelCost;
    delete [] levelParCost;
    delete [] partition2Level;
    delete [] initLevelPtr;
    delete [] initLevelSet;
    return 1;

}

#endif //CHOLOPENMP_MERGEWITHSLACK_H
