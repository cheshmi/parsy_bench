//
// Created by kazem on 3/4/18.
//

#ifndef PROJECT_MKL_UTILS_H
#define PROJECT_MKL_UTILS_H

#include <mkl.h>

#include <cstring> // memset
#include <cassert> // for assert
#include <omp.h>   // just for timer
#include <cstdlib> // for random generation
#include <cstdio>  // for printf

// Utils macro
#define Min(x,y) ((x)<(y)?(x):(y))
#define Max(x,y) ((x)>(y)?(x):(y))
#define Abs(x) ((x)>(0)?(x):-(x))

//////////////////////////////////////////////////////////////////////////
// COO part
//////////////////////////////////////////////////////////////////////////
struct COOArrays{
 MKL_INT m;      //< the dimension of the matrix
 MKL_INT nnz;    //< the number of nnz inside the matrix
 double *val;    //< the values (size = nnz)
 MKL_INT *rowind;//< the row indexes (size = nnz)
 MKL_INT *colind;//< the col indexes (size = nnz)

 /** simply set ptr to null */
 COOArrays(){
  val = NULL;
  rowind = NULL;
  colind = NULL;
 }

 /** delete ptr */
 ~COOArrays(){
  delete[] val;
  delete[] rowind;
  delete[] colind;
 }
};

/** see https://software.intel.com/fr-fr/node/520817#38F0A87C-7884-4A96-B83E-CEE88290580F */
void compute_COO(COOArrays& coo, double *x , double *y){
 char transa = 'N';
 // void mkl_cspblas_dcoogemv(char *transa, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
 mkl_cspblas_dcoogemv(&transa, &coo.m, coo.val, coo.rowind, coo.colind, &coo.nnz, x, y);
}

//////////////////////////////////////////////////////////////////////////
// CRS part
//////////////////////////////////////////////////////////////////////////

struct CRSArrays{
 MKL_INT m;  //< the dim of the matrix
 MKL_INT nnz;//< the number of nnz (== ia[m])
 double *a;  //< the values (of size NNZ)
 MKL_INT *ia;//< the usual rowptr (of size m+1)
 MKL_INT *ja;//< the colidx of each NNZ (of size nnz)

 CRSArrays(){
  a = NULL;
  ia = NULL;
  ja= NULL;
 }

 ~CRSArrays(){
/*  delete[] a;
  delete[] ia;
  delete[] ja;*/
 }
};

/** See https://software.intel.com/fr-fr/node/520849#449CA855-CE5B-4061-B003-70D078CA5E05 */
void COO_to_CRS(COOArrays& coo, CRSArrays* crs){
 MKL_INT job[6] = {1,//if job(1)=1, the matrix in the coordinate format is converted to the CRS format.
                   0,//If job(2)=0, zero-based indexing for the matrix in CRS format is used;
                   0,//If job(3)=0, zero-based indexing for the matrix in coordinate format is used;
                   0,
                   coo.nnz,//job(5)=nnz - sets number of the non-zero elements of the matrix A if job(1)=1.
                   0 //If job(6)=0, all arrays acsr, ja, ia are filled in for the output storage.
 };
 // Init crs
 crs->m = coo.m;
 crs->nnz = coo.nnz;
 crs->a = new double[crs->nnz];
 crs->ia = new MKL_INT[crs->m+1];
 crs->ja = new MKL_INT[crs->nnz];
 MKL_INT nnz = coo.nnz;
 MKL_INT info;
 mkl_dcsrcoo(job , &coo.m,
             crs->a , crs->ja , crs->ia , &nnz ,
             coo.val, coo.rowind , coo.colind , &info );
}

/** See https://software.intel.com/fr-fr/node/520815#D840F0E5-E41A-4E91-94D2-FEB320F93E91 */
void compute_CRS( CRSArrays& crs, double *x , double *y){
 char transa = 'N';
 // void mkl_cspblas_dcsrgemv (const char *transa , const MKL_INT *m , const double *a , const MKL_INT *ia , const MKL_INT *ja , const double *x , double *y );
 mkl_cspblas_dcsrgemv(&transa, &crs.m , crs.a , crs.ia , crs.ja , x , y);
}

//////////////////////////////////////////////////////////////////////////
// DIA part
//////////////////////////////////////////////////////////////////////////

struct DIAArrays{
 MKL_INT m;      //< the dimensio of the matrix
 MKL_INT nnz;    //< the number of nnz inside the matrix
 double *val;    //< the NNZ values - may include zeros (of size lval*ndiag)*/
 MKL_INT *idiag; //< distance from the diagonal (of size ndiag)
 MKL_INT lval;   //< leading where the diagonals are stored >= m,  which is the declared leading dimension in the calling (sub)programs
 MKL_INT ndiag;  //< number of diagonals that have at least one nnz

 DIAArrays(){
  val = NULL;
  idiag = NULL;
 }

 ~DIAArrays(){
  delete[] val;
  delete[] idiag;
 }
};

/** https://software.intel.com/fr-fr/node/520852#00B3CA58-E0E4-4ED9-B42B-BC6338AB461D */
void CRS_to_DIA(CRSArrays& crs, DIAArrays* dia){
 MKL_INT job[6] = {0,//If job(1)=0, the matrix in the CRS format is converted to the diagonal format;
                   0,//If job(2)=0, zero-based indexing for the matrix in CRS format is used;
                   1,//if job(3)=1, one-based indexing for the matrix in the diagonal format is used.
                   0,
                   0,
                   10//If job(6)=10, diagonals are selected internally, and acsr_rem, ja_rem, ia_rem are not filled in for the output storage.
 };
 dia->m = crs.m;
 dia->nnz = crs.nnz;
 dia->lval = dia->m;
 dia->ndiag = 0;
 // We need to count the number of diagonals with NNZ
 {
  unsigned* usedDiag = new unsigned[crs.m*2-1];
  memset(usedDiag, 0, sizeof(unsigned)*(crs.m*2-1));

  for(int idxRow = 0 ; idxRow < crs.m ; ++idxRow){
   for(int idxVal = crs.ia[idxRow] ; idxVal < crs.ia[idxRow+1] ; ++idxVal){
    const int idxCol = crs.ja[idxVal];
    const int diag = crs.m-idxRow+idxCol-1;
    assert(0 <= diag && diag < crs.m*2-1);
    if(usedDiag[diag] == 0){
     usedDiag[diag] = 1;
     dia->ndiag += 1;
    }
   }
  }

  delete[] usedDiag;
 }
 // Allocate the working arrays
 dia->val = new double[dia->ndiag*dia->lval];
 dia->idiag = new MKL_INT[dia->ndiag];

 // void mkl_dcsrdia (const MKL_INT *job , const MKL_INT *n , double *acsr , MKL_INT *ja ,
 // MKL_INT *ia , double *adia , const MKL_INT *ndiag , MKL_INT *distance , MKL_INT *idiag ,
 // double *acsr_rem , MKL_INT *ja_rem , MKL_INT *ia_rem , MKL_INT *info );

 MKL_INT info;
 mkl_dcsrdia(job , &crs.m ,crs.a , crs.ja, crs.ia ,
             dia->val , &dia->lval , dia->idiag , &dia->ndiag ,
             NULL , NULL , NULL , &info );
}

/** https://software.intel.com/fr-fr/node/520806#FCB5B469-8AA1-4CFB-88BE-E2F22E9E2AF0 */
void compute_DIA_1idx(DIAArrays& dia , double *x , double *y){
 char transa = 'N';
 // void mkl_ddiagemv (const char *transa , const MKL_INT *m , const double *val , const MKL_INT *lval , const MKL_INT *idiag , const MKL_INT *ndiag , const double *x , double *y );
 mkl_ddiagemv(&transa, &dia.m , dia.val , &dia.lval , dia.idiag , &dia.ndiag , x , y);
}

//////////////////////////////////////////////////////////////////////////
// BCRS part
//////////////////////////////////////////////////////////////////////////

struct BCRSArrays{
 MKL_INT m;
 MKL_INT nnz;
 MKL_INT nbBlocks;
 MKL_INT nbBlockRows;
 MKL_INT lb;/*size of blocks*/
 MKL_INT ldabsr;/*leading >= lb*lb*/
 double *a;/*values(m*lb*lb)*/
 MKL_INT *ia;/*i(m+1)*/
 MKL_INT *ja;/*j(m+1)*/
 MKL_INT allocatedBlocks;

 BCRSArrays(){
  a = NULL;
  ia = NULL;
  ja = NULL;
 }

 ~BCRSArrays(){
  delete[] a;
  delete[] ia;
  delete[] ja;
 }
};

/** https://software.intel.com/fr-fr/node/520850#3A22B45C-4604-4444-B6FE-205A5CD4E667 */
void CRS_to_BCRS(CRSArrays& crs, BCRSArrays* bcrs, const int blockSize){
 MKL_INT job[6] = {0,//If job(1)=0, the matrix in the CSR format is converted to the BSR format;
                   0,//If job(2)=0, zero-based indexing for the matrix in CSR format is used;
                   0,//If job(3)=0, zero-based indexing for the matrix in the BSR format is used;
                   0,
                   0,
                   1 //If job(6)>0, all output arrays absr, jab, and iab are filled in for the output storage.
 };
 bcrs->m = crs.m;
 bcrs->nnz = crs.nnz;
 bcrs->nbBlocks = 0;
 // We need to count the number of blocks (aligned!)
 {
  const MKL_INT maxBlockPerRow = (bcrs->m+blockSize-1)/blockSize;
  unsigned* usedBlocks = new unsigned[maxBlockPerRow];

  for(int idxRow = 0 ; idxRow < crs.m ; ++idxRow){
   if(idxRow%blockSize == 0){
    memset(usedBlocks, 0, sizeof(unsigned)*maxBlockPerRow);
   }
   for(int idxVal = crs.ia[idxRow] ; idxVal < crs.ia[idxRow+1] ; ++idxVal){
    const int idxCol = crs.ja[idxVal];
    if(usedBlocks[idxCol/blockSize] == 0){
     usedBlocks[idxCol/blockSize] = 1;
     bcrs->nbBlocks += 1;
    }
   }
  }

  delete[] usedBlocks;
 }
 bcrs->nbBlockRows = (bcrs->m+blockSize-1)/blockSize;
 bcrs->lb = blockSize;
 bcrs->ldabsr = bcrs->lb*bcrs->lb;
 bcrs->a = new double[bcrs->nbBlocks*bcrs->lb*bcrs->lb]();
 bcrs->ia = new MKL_INT[2*bcrs->nbBlockRows+1]();
 bcrs->ja = new MKL_INT[2*bcrs->nbBlocks]();

 // void mkl_dcsrbsr (const MKL_INT *job , const MKL_INT *m , const MKL_INT *mblk ,
 // const MKL_INT *ldabsr , double *acsr , MKL_INT *ja , MKL_INT *ia , double *absr ,
 // MKL_INT *jab , MKL_INT *iab , MKL_INT *info );
 MKL_INT info;
 mkl_dcsrbsr(job , &crs.m ,
             &bcrs->lb , &bcrs->ldabsr,
             crs.a , crs.ja , crs.ia ,
             bcrs->a, bcrs->ja , bcrs->ia, &info );
 assert(bcrs->ia[bcrs->nbBlockRows] == bcrs->nbBlocks);
 if(info!=0)
  printf("Unsuccessful conversion CSR to BSR;");
}

/** https://software.intel.com/fr-fr/node/520816#366F2854-A2C0-4661-8CE7-F478F8E6B613 */
void compute_BSR(BCRSArrays& bcsr, double *x , double *y ){
 char transa = 'N';
 // void mkl_cspblas_dbsrgemv (const char *transa , const MKL_INT *m , const MKL_INT *lb , const double *a , const MKL_INT *ia , const MKL_INT *ja , const double *x , double *y );
 mkl_cspblas_dbsrgemv(&transa, &bcsr.nbBlockRows , &bcsr.lb ,
                      bcsr.a , bcsr.ia , bcsr.ja , x , y);
}

#endif //PROJECT_MKL_UTILS_H
