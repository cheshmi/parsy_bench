//
// Created by kazem on 8/13/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <c++/7/string>
#include <fstream>
#include <sstream>
#include <iostream>




/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                             double *, int    *,    int *, int *,   int *, int *,
                             int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                                    double *, int *);




/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
bool readMatrix(std::string fName, int &n, int &NNZ, int* &col,//FIXME change col type to size_t
                int* &row, double* &val){
 /*This function reads the input matrix from "fName" file and
  * allocate memory for matrix A, L and U.
  * - The input file is a coordinate version and e
  * ach row of the file shows (col, row, nnz)
  * - The matrices are zero-indexed
  */

 std::ifstream inFile;
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
 }
 std::istringstream issDim(line);
 if (!(issDim >> n >> n >> NNZ)){
  std::cout<<"The matrix dimension is missing\n";
  return false;
 }
 if(n <= 0 || NNZ <= 0)
  return false;
 //col = new int[n + 1]();
 col = (int*) calloc((n+1), sizeof(int));
 // colL = new int[n + 1]; colU = new int[n + 1];
 //row = new int[NNZ];
 row = (int*) calloc((NNZ), sizeof(int));
 // rowL = new int[factorSize]; rowU = new int[factorSize];
 //val = new double[NNZ];
 val = (double*) calloc((NNZ), sizeof(double));
 // valL = new double[factorSize]; valU = new double[factorSize];
 if(!val || !col || !row)
  return false;
 //Initializing the result vector
 int y, x, colCnt=0, nnzCnt=0;
 double value;

 col[0]=0;
 for (int i = 0; nnzCnt<NNZ; ) {//Reading from file row by row
  inFile>>x;x--;
  inFile>>y;y--;//zero indexing
  inFile>>value;
  if(y > n)
   return false;
  if(y==i){
   val[nnzCnt]=value;
   row[nnzCnt]=x;
   colCnt++; nnzCnt++;
  }
  else{//New col
   col[i+1]=col[i]+colCnt;
   i++;//next iteration
   colCnt=1;
   val[nnzCnt]=value;
   row[nnzCnt]=x;
   nnzCnt++;
  }

 }
 col[n]= col[n - 1] + colCnt;//last col

 return true;
}

int main( int argc, char *argv[] )
{
 /* Matrix data. */
/*
 int    n = 8;
 int    ia[ 9] = { 0, 4, 7, 9, 11, 14, 16, 17, 18 };
 int    ja[18] = { 0,    2,       5, 6,
                   1, 2,    4,
                   2,             7,
                   3,       6,
                   4, 5, 6,
                   5,    7,
                   6,
                   7 };
 double  a[18] = { 7.0,      1.0,           2.0, 7.0,
                   -4.0, 8.0,           2.0,
                   1.0,                     5.0,
                   7.0,           9.0,
                   5.0, 1.0, 5.0,
                   0.0,      5.0,
                   11.0,
                   5.0 };




 int      nnz = ia[n];
*/
 std::string f1 = argv[1];
 int n, nnz;
 int *ia, *ja;
 double *a;

 int iter_no=5;

 readMatrix(f1, n, nnz, ia, ja, a);

 int      mtype = 2;        /* Real symmetric matrix */

 /* RHS and solution vectors. */
 //double   b[8], x[8];
 double *b, *x;
 x = (double*)calloc(n, sizeof(double));
 b = (double*)calloc(n, sizeof(double));
 int      nrhs = 1;          /* Number of right hand sides. */

 /* Internal solver memory pointer pt,                  */
 /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
 /* or void *pt[64] should be OK on both architectures  */
 void    *pt[64];

 /* Pardiso control parameters. */
 int      iparm[64];
 double   dparm[64];
 for (int m = 0; m < 64; ++m) {
  iparm[m] = 0;
  dparm[m] = 0.0;
 }
 int      maxfct, mnum, phase, error, msglvl, solver;

 /* Number of processors. */
 int      num_procs;

 /* Auxiliary variables. */
 char    *var;
 int      i, k;

 double   ddum;              /* Double dummy */
 int      idum;              /* Integer dummy. */


/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters.                                */
/* -------------------------------------------------------------------- */
 //int numThread=4;
 //int blasThreads=1;
 //omp_set_num_threads(numThread);



 error = 0;
 solver=0;/* use sparse direct solver */
 pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

 if (error != 0)
 {
  if (error == -10 )
   printf("No license file found \n");
  if (error == -11 )
   printf("License is expired \n");
  if (error == -12 )
   printf("Wrong username or hostname \n");
  return 1;
 }
 else
  printf("[PARDISO]: License check was successful ... \n");

 /* Numbers of processors, value of OMP_NUM_THREADS */
 var = getenv("OMP_NUM_THREADS");
 if(var != NULL)
  sscanf( var, "%d", &num_procs );
 else {
  printf("Set environment OMP_NUM_THREADS to 1");
  exit(1);
 }
 iparm[2]  = num_procs;
 iparm[51] = 1;
 //iparm[0]=0;

 maxfct = 1;		/* Maximum number of numerical factorizations.  */
 mnum   = 1;         /* Which factorization to use. */

 msglvl = 1;         /* Print statistical information  */
 error  = 0;         /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
 for (i = 0; i < n+1; i++) {
  ia[i] += 1;
 }
 for (i = 0; i < nnz; i++) {
  ja[i] += 1;
 }

 /* Set right hand side to i. */
 for (i = 0; i < n; i++) {
  b[i] = i;
 }

/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

 pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
 if (error != 0) {
  printf("\nERROR in consistency of matrix: %d", error);
  exit(1);
 }

/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

 pardiso_chkvec (&n, &nrhs, b, &error);
 if (error != 0) {
  printf("\nERROR  in right hand side: %d", error);
  exit(1);
 }

/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */

 pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
 if (error != 0) {
  printf("\nERROR right hand side: %d", error);
  exit(1);
 }

/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
 phase = 11;

 pardiso (pt, &maxfct, &mnum, &mtype, &phase,
          &n, a, ia, ja, &idum, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error, dparm);

 if (error != 0) {
  printf("\nERROR during symbolic factorization: %d", error);
  exit(1);
 }
 printf("\nReordering completed ... ");
 printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
 printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

/* -------------------------------------------------------------------- */
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */
 phase = 22;
 iparm[32] = 0; /* compute determinant */
 for (int j = 0; j < iter_no; ++j) {
  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
           &n, a, ia, ja, &idum, &nrhs,
           iparm, &msglvl, &ddum, &ddum, &error,  dparm);

 }


 if (error != 0) {
  printf("\nERROR during numerical factorization: %d", error);
  exit(2);
 }
 printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */
 phase = 33;

 iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

 pardiso (pt, &maxfct, &mnum, &mtype, &phase,
          &n, a, ia, ja, &idum, &nrhs,
          iparm, &msglvl, b, x, &error,  dparm);

 if (error != 0) {
  printf("\nERROR during solution: %d", error);
  exit(3);
 }

 printf("\nSolve completed ... ");
/* printf("\nThe solution of the system is: ");
 for (i = 0; i < n; i++) {
  printf("\n x [%d] = % f", i, x[i] );
 }
 printf ("\n\n");*/


/* -------------------------------------------------------------------- */
/* ... Inverse factorization.                                           */
/* -------------------------------------------------------------------- */

 if (solver == 0)
 {
  printf("\nCompute Diagonal Elements of the inverse of A ... \n");
  phase = -22;
  iparm[35]  = 1; /*  no not overwrite internal factor L */

  pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
           iparm, &msglvl, b, x, &error,  dparm);

  /* print diagonal elements */
  /*for (k = 0; k < n; k++)
  {
   int j = ia[k]-1;
   printf ("Diagonal element of A^{-1} = %d %d %32.24e\n", k, ja[j]-1, a[j]);
  }*/

 }


/* -------------------------------------------------------------------- */
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */
 for (i = 0; i < n+1; i++) {
  ia[i] -= 1;
 }
 for (i = 0; i < nnz; i++) {
  ja[i] -= 1;
 }

/* -------------------------------------------------------------------- */
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */
 phase = -1;                 /* Release internal memory. */

 pardiso (pt, &maxfct, &mnum, &mtype, &phase,
          &n, &ddum, ia, ja, &idum, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error,  dparm);

 free(ia);
 free(ja);
 free(a);
 free(x);
 free(b);

 return 0;
}
