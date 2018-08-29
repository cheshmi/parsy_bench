//
// Created by kazem on 7/23/17.
//

#ifndef TRIANGOPENMP_COARSE_INSPECTION_LEVEL_H
#define TRIANGOPENMP_COARSE_INSPECTION_LEVEL_H

/*
 *
 */
int buildCoarseLevelSet_CSC(int n, int nnz, int *Lp, int *Li, int *&levelPtr,
                      int *&levelSet){
    int begin=0,end=n-1;
    int curLevel=0, curLevelCol=0;
    levelPtr = new int[n+1]();
    levelSet = new int[n]();
    int *inDegree = new int[n]();
    bool *visited = new bool[n]();
    int *node2Level = new int[n]();
    for (int i = 0; i < nnz; ++i) {//O(nnz)
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
                node2Level[i] = curLevel; //map node to level
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
int coarsening(int n, int nnz, int *Lp, int *Li,
               int nLevels, int *levelPtr, int *levelSet, int *node2Level,
               int *coarseLevelPtr, int *treeLevelPtr, int *treeLevelSet){
    int levelArray0=0, levelArray1=0, levelArray2=0;
    int *mark = new int[n]();
    coarseLevelPtr[0]=0;
    for (int l = 0; l < nLevels; ++l) {
        for (int ls = levelPtr[l]; ls < levelPtr[l+1]; ++ls) {
            int curNode = levelSet[ls];
            treeLevelSet[levelArray2++];
            mark[curNode] = 1;//visited
            levelArray1=1;
            //for every edge of curNode
            for (int nxt = Lp[curNode]; nxt < Lp[curNode + 1]; ++nxt) {
                //see if it is existed in the next level
                if(node2Level[ Li[nxt] ] == l && !mark[Li[nxt]]){
                    treeLevelPtr[levelArray2++] = Li[nxt];
                    mark[Li[nxt]] = 1;
                }
            }
            //treeLevelPtr[]
        }

    }
    return 1;
}
#endif //TRIANGOPENMP_COARSE_INSPECTION_LEVEL_H
