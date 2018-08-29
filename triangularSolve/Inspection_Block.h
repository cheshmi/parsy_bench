//
// Created by kazem on 7/18/17.
//

#ifndef TRIANGOPENMP_INSPECTION_BCSC_H
#define TRIANGOPENMP_INSPECTION_BCSC_H
void superNodeDetection(int n, int *col, int *row,int *col2sup,int &supNo){
    int prev, cur;
    bool sim;
    supNo=0;
    col2sup[0]=0;
    for (int i = 1; i < n; ++i) {
        sim=true;
        for (prev = col[i-1], cur = col[i]; 1; ) {
            if((row[prev] == i-1 && prev < col[i]) || (row[prev] == i && prev < col[i])){
                //skip diagonal block of prev col
                ++prev;
                continue;
            }
            if(row[cur]==i){//skip diagonal block of cur col
                ++cur;
                continue;
            }
            if(prev - col[i] != cur - col[i+1]){ // the off-diagonals length
                sim=false;
                break;
            }else if(prev - col[i]==0)
                break;
            if(row[prev] != row[cur]){//now off-diagonals
                sim=false;
                break;
            }else{
                ++prev;++cur;
                if (prev >= col[i] || cur >= col[i+1])
                    break;
            }
        }
#if DEBUG >0
        if(i<0 || i>=n)
            printf("TTT \n");
#endif
        if(sim ){//col cur and nxt are similar
            col2sup[i]=supNo;
        }else{
            supNo++;
            col2sup[i]=supNo;
        }
    }
    supNo+= sim ;
}

void calcSize(int n, int *col, int *newCol, int *col2sup, int *sup2col, int supNo,
              int& newRowSize, int& newNNZ){
    newNNZ=0;newRowSize=0;
    newCol[0]=0;
    int curCol = 0, cnt=0, firstCol=0, tmpSize=0;
    for (int j = 0; j < supNo; ++j) {
        for (cnt=0; col2sup[curCol]==j && curCol<n; ++curCol, ++cnt);
        sup2col[j]=curCol;
        firstCol= j!=0 ? sup2col[j-1] : 0;
        cnt--;//To find the last col offset of supernode
        tmpSize=col[firstCol+cnt+1]-col[firstCol+cnt]+cnt;
        newRowSize+=tmpSize;//diagonal dense part + off diagonal sparse one
        for (int i = firstCol; i < curCol; ++i) {
            newNNZ+=tmpSize;
            newCol[i+1]=newNNZ;
        }
    }
#if DEBUG>0
    for (int i = 0; i < supNo; ++i) {//printing the result
        std::cout<<"Snode: "<<i<<" col: "<<sup2col[i]<<"\n";
    }
#endif
}
void createFormat(int n, int *col, int *row, double *val, int NNZ, int* &newRow, int newRowSize, double* &newVal,
                  int* &row_ptr, int* &newCol, int* &col2sup, int* &sup2col, int &supNo){
    int l=0, curOrig, wdth=0, c, origRow, firstCol;
    int firstRow;
    int offDiagofFirst;
    //Creating the new blocked format
    for (int s = 0; s < supNo; ++s) {
        //specifying the compressed row
        firstCol=s!=0 ? sup2col[s-1] : 0;
        firstRow = row[col[firstCol]];
        wdth=sup2col[s]-firstCol;
        row_ptr[firstCol]=l;
        offDiagofFirst=0;
        for (int r = 0; r < wdth; ++r, ++l, ++firstRow) {//row indices for the dense part
            newRow[l]=firstRow;
            if(row[col[firstCol]+offDiagofFirst]<=firstRow)//zero in dense part of first col
                offDiagofFirst++;
        }
        for (int r = col[firstCol]+offDiagofFirst; r < col[firstCol+1]; ++r, ++l) {
            newRow[l]=row[r];
        }
        for (int i = firstCol; i < sup2col[s]; ++i) {
            row_ptr[i] = row_ptr[firstCol];//copying the row number

            firstRow = row[col[firstCol]];
            for (c = newCol[i], curOrig=col[i]; c < newCol[i]+wdth; ++c, ++firstRow) {//Diagonal zero padding
                if(row[curOrig] == firstRow){
                    newVal[c] = val[curOrig];
                    curOrig++;
                }else{
                    newVal[c]=0;
                }
            }
            for (c = newCol[i]+wdth; c < newCol[i+1]; ++curOrig, ++c) {
                newVal[c]=val[curOrig];
            }
        }
    }
    row_ptr[n]=newRowSize;
}
void toSuperNode(int n, int *col, int *row, double *val, int NNZ, int* &newRow, double* &newVal,
                 int* &row_ptr, int* &newCol, int* &col2sup, int* &sup2col, int &supNo){
    int prev, cur;
    int newRowSize=0, newNNZ=0;
    bool sim;
    supNo=0;
    col2sup[0]=0;
    for (int i = 1; i < n; ++i) {
        sim=true;
        for (prev = col[i-1], cur = col[i]; 1; ) {
            if((row[prev] == i-1 && prev < col[i]) || (row[prev] == i && prev < col[i])){
                //skip diagonal block of prev col
                ++prev;
                continue;
            }
            if(row[cur]==i){//skip diagonal block of cur col
                ++cur;
                continue;
            }
            if(prev - col[i] != cur - col[i+1]){ // the off-diagonals length
                sim=false;
                break;
            }else if(prev - col[i]==0)
                break;
            if(row[prev] != row[cur]){//now off-diagonals
                sim=false;
                break;
            }else{
                ++prev;++cur;
                if (prev >= col[i] || cur >= col[i+1])
                    break;
            }
        }
        if(sim ){//col cur and nxt are similar
            col2sup[i]=supNo;
        }else{
            supNo++;
            col2sup[i]=supNo;
        }
    }
    supNo+= sim ;
    sup2col = new int[supNo];
    newCol = new int[n];
    newCol[0]=0;
    int curCol = 0, cnt=0, firstCol=0, tmpSize=0;
    for (int j = 0; j < supNo; ++j) {
        for (cnt=0; col2sup[curCol]==j; ++curCol, ++cnt);
        sup2col[j]=curCol;
        firstCol= j!=0 ? sup2col[j-1] : 0;
        cnt--;//To find the last col offset of supernode
        tmpSize=col[firstCol+cnt+1]-col[firstCol+cnt]+cnt;
        newRowSize+=tmpSize;//diagonal dense part + off diagonal sparse one
        for (int i = firstCol; i < curCol; ++i) {
            newNNZ+=tmpSize;
            newCol[i+1]=newNNZ;
        }
    }
#if DEBUG>0
    for (int i = 0; i < supNo; ++i) {//printing the result
        std::cout<<"Snode: "<<i<<" col: "<<sup2col[i]<<"\n";
    }
#endif
    newRow = new int[newRowSize];
    newVal = new double[newNNZ];
    row_ptr = new int[n+1];
    int l=0, curOrig, wdth=0, c, origRow;
    int firstRow;
    //Creating the new blocked format
    for (int s = 0; s < supNo; ++s) {
        //specifying the compressed row
        firstCol=s!=0 ? sup2col[s-1] : 0;
        firstRow = row[col[firstCol]];
        wdth=sup2col[s]-firstCol;
        row_ptr[firstCol]=l;
        for (int r = 0; r < wdth; ++r, ++l, ++firstRow) {//row indices for the dense part
            newRow[l]=firstRow;
        }
        for (int r = col[firstCol]+wdth; r < col[firstCol+1]; ++r, ++l) {
            newRow[l]=row[r];
        }
        for (int i = firstCol; i < sup2col[s]; ++i) {
            row_ptr[i] = row_ptr[firstCol];//copying the row number

            firstRow = row[col[firstCol]];
            for (c = newCol[i], curOrig=col[i]; c < newCol[i]+wdth; ++c, ++firstRow) {//Diagonal zero padding
                if(row[curOrig] == firstRow){
                    newVal[c] = val[curOrig];
                    curOrig++;
                }else{
                    newVal[c]=0;
                }
            }
            for (c = newCol[i]+wdth; c < newCol[i+1]; ++curOrig, ++c) {
                newVal[c]=val[curOrig];
            }
        }
    }
    row_ptr[n]=newRowSize;
#if DEBUG>0
    for (int i = 0; i < newRowSize; ++i) {//printing the result
        std::cout<<"#"<<i<<" idx: "<<newRow[i]<<"\n";
    }
    for (int i = 0; i < newNNZ; ++i) {//printing the result
        std::cout<<"#"<<i<<" col: "<<newVal[i]<<"\n";
    }
    for (int i = 0; i < n+1; ++i) {//printing the result
        std::cout<<"#"<<i<<" row: "<<row_ptr[i]<<"\n";
    }
#endif
}
#endif //TRIANGOPENMP_INSPECTION_BCSC_H
