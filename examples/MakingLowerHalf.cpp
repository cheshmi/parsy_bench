//
// Created by kazem on 7/16/18.
//

#include <fstream>
#include <iostream>
#include <sstream>

bool printLower(std::string fName){
 /*This function reads the input matrix from "fName" file and
  * allocate memory for matrix A, L and U.
  * - The input file is a coordinate version and e
  * ach row of the file shows (col, row, nnz)
  * - The matrices are zero-indexed
  */

 std::ifstream inFile;
 double tol=0.1;
 size_t n, NNZ;
 inFile.open(fName);
 std::string line,banner, mtx, crd, arith, sym;
 /*  File format:
  *    %%MatrixMarket matrix coordinate real general/symmetric/...
  *    % ...
  *    % (optional comments)
  *    % ...
  *    #rows    #non-zero
  *    Triplet in the rest of lines: row    col    value
  */
 std::getline(inFile,line);
 //std::cout<<line<<"\n";
 std::cout<<"%%MatrixMarket matrix coordinate real symmetric\n";
 for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
 std::istringstream iss(line);
 if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
  std::cout<<"Invalid header (first line does not contain 5 tokens)\n";
  return false;
 }

 if(banner.compare("%%matrixmarket")) {
  std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")\n";
  return false;
 }
 if(mtx.compare("matrix")) {
  std::cout<<"Not a matrix; this driver cannot handle that.\"\n";
  return false;
 }
 if(crd.compare("coordinate")) {
  std::cout<<"Not in coordinate format; this driver cannot handle that.\"\n";
  return false;
 }
 if(arith.compare("real")) {
  if(!arith.compare("complex")) {
   std::cout<<"Complex matrix; use zreadMM instead!\n";
   return false;
  }
  else if(!arith.compare("pattern")) {
   std::cout<<"Pattern matrix; values are needed!\n";
   return false;
  }
  else {
   std::cout<<"Unknown arithmetic\n";
   return false;
  }
 }
 while (!line.compare(0,1,"%"))
 {
  std::getline(inFile, line);
  //std::cout<<line<<"\n";
 }
 std::istringstream issDim(line);
 if (!(issDim >> n >> n >> NNZ)){
  std::cout<<"The matrix dimension is missing\n";
  return false;
 }
 if(n <= 0 || NNZ <= 0)
  return false;
 //Initializing the result vector
 int y, x, colCnt=0, nnzCnt=0;
 double value;
 std::cout<<n<< " "<<n<<" " << (NNZ-n)/2+n <<"\n";

 for (; nnzCnt<NNZ; ) {//Reading from file row by row
  inFile>>x;
  inFile>>y;
  inFile>>value;
  if(y > n)
   return false;
  if(x>=y){
   if(x==y){
    double valTmp=value>=0?value+tol:value-tol;
    std::cout<<x<<" "<<y<<" " <<valTmp<<"\n";

   }
   else
    std::cout<<x<<" "<<y<<" " <<value<<"\n";
  }
  nnzCnt++;
 }

 return true;
}


int main(int argc, char *argv[]) {


 //std::string f1 = "/home/kazem/UFDB/SymFull/cbuckle.mtx";
 //string fName1= "/home/kazem/UFDB/SymSparsity/cbuckle_sparsity_amd.dat";

 if (argc < 2)
  printf("input args are missing");
 std::string f1 = argv[1];
//    string fName1 = argv[2];
 int *col, *row;
 int *rowL;
 size_t *colL;
 double *valL;
 double *y, *val, *x;
 int maxSupWid, maxCol;
 size_t n, NNZ;

 if (!printLower(f1))
  return -1;
}