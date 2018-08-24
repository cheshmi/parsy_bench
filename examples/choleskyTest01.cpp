#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <cholmod.h>
#include <cholmod_function.h>
#include <omp.h>

#include "Ordering.h"
#include "Inspection_Prune.h"
#include "Inspection_Block.h"
#include "Util.h"
#include "PB_Cholesky.h"
#include "LSparsity.h"

#include "mkl.h"

#include "parallel_PB_Cholesky_05.h"
#include "Cholesky_Perf_Comm_modeling.h"

#include <complex>
//#define ANALYZE
//#define ANALYZE_THEORY
#define CPUTIME (SuiteSparse_time ( ))
#undef DEBUG
using namespace std;


struct cholTime{
 double alltogether, parallel, rootNodes;
 std::vector<double> tArray;
 cholTime():alltogether(0),parallel(0),rootNodes(0){

 }
};

bool retCMP(cholTime a, cholTime b){
 return a.alltogether>b.alltogether;}

bool cmpPI(profilingInfo a, profilingInfo b){
 return a.theoreticalTime>b.theoreticalTime;}

int main(int argc, char *argv[]) {


 //std::string f1 = "/home/kazem/UFDB/SymFull/cbuckle.mtx";
 //string fName1= "/home/kazem/UFDB/SymSparsity/cbuckle_sparsity_amd.dat";

 if(argc<2)
  printf("input args are missing");
 std::string f1 = argv[1];
//    string fName1 = argv[2];
 int *col, *row;
 int  *rowL;
 size_t *colL;
 double *valL;
 double  *y, *val, *x;
 int  maxSupWid, maxCol;
 size_t n, NNZ;
 std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;
 double durationSym = 0 ,duration3 = 0, duration2=0, duration1=0;
 long totalIntraContribNNZs=0,totalInterCotribNNZs=0, numOfEdgeCuts=0;
 int numberOfIntraCore=0,numberOfInterCore=0;

 profilingInfo pi;
 std::vector<profilingInfo> piArray;

 if (!readMatrix(f1,n,NNZ,col,row,val))
  return -1;

/* if(wasteFile.fail())
  return -1;*/
 int numThread = atoi(argv[2]);
 int chunk = atoi(argv[3]);
 int costParam = atoi(argv[4]);//Inner parts
 int levelParam = atoi(argv[5]);// level distance
 int blasThreads = atoi(argv[6]);
 int finalSeqNode = atoi(argv[7]);
 size_t *inPerm;
 if(argc>8){
  std::string orderFileName = argv[8];
  inPerm = new size_t[n]();
  readOrdering(orderFileName,n,inPerm);
 }
 omp_set_num_threads(numThread);

 // MKL_Set_Num_Threads(1);
 MKL_Domain_Set_Num_Threads(blasThreads,MKL_DOMAIN_BLAS);
 //cout<<"---" <<MKL_Domain_Get_Max_Threads(MKL_DOMAIN_BLAS)<<"\n";
/*    chrono::time_point<std::chrono::system_clock> start, end;
    double durationAlltogether = 0, durationPruned=0, durationSym=0,
            ordering=0, durationBlock=0, durationAllSmall=0;
    chrono::duration<double> elapsed_seconds;*/

 int factorSize=0;
 double timing[4];//for time measurement

 ifstream spFile1;
//    spFile1.open(fName1);
 int *prunePtr, *pruneSet;
 int *levelPtr = NULL, *levelSet = NULL, *parPtr = NULL,
   *partition =NULL;
 int nLevels=0, nPar=0;

 //double timingChol[4]={.0,.0,.0,.0};//for time measurement
 double *timingChol = new double[4+numThread]();//for time measurement
 double orderingTime=0;
/* int nrelax[3] = {4,16,48};//TODO
 double zrelax[3] = {0.8,0.1,0.05};*/
 int nrelax[3] = {4,16,48};//TODO
 double zrelax[3] = {0.8,0.1,0.05};
 int status=0;
 double *contribs;
 int super_max = 164; //tunig parameter for the max size of supernodes TODO: find it in analysis
 int col_max = n;
// int *col2sup=new int[n]();
 //int *blockSet;
 //contribs = new double[super_max * col_max]();
 size_t *li_ptr; //= new size_t[n+1];
 //int *map = new int[n]();
 //colL = new size_t[n + 1]();
 CSC *Amat = new CSC;
 Amat->nzmax = NNZ; Amat->ncol=Amat->nrow=n;
 Amat->stype=-1;Amat->xtype=CHOLMOD_REAL;Amat->packed=TRUE;
 Amat->p = col; Amat->i = row; Amat->x=val; Amat->nz = NULL;
 Amat->sorted = TRUE;
 start = std::chrono::system_clock::now();
 BCSC *L = analyze_p2(1,Amat,NULL,NULL,nrelax,zrelax,
                      n,prunePtr,pruneSet,
                      nLevels, levelPtr,levelSet,
                      nPar, parPtr, partition,
                      costParam,levelParam,finalSeqNode,
                      status, maxSupWid, maxCol, orderingTime,
                      inPerm);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 durationSym=elapsed_seconds.count();
#if 0
 cout<<"\n";
    for (int j = 0; j < L->nsuper; ++j) {
        int colLength = L->pi[j+1]-L->pi[j];
        int supWid = L->super[j+1]-L->super[j];
        cout<<"Supernode: "<<j<<" Len and Wid are: "
            <<colLength<<","<<supWid<<"\n";
        for (int i = L->pi[j]; i < L->pi[j+1]; ++i) {
            cout<<L->s[i]<<",";
        }
        cout<<"\n";
    }
#endif
 //Some conversion for sympiler

 int colLength=0;

/* for (int j = 0; j < L->nsuper; ++j) {
  int curCol = L->super[j];
  int nxtCol = L->super[j+1];
  colLength = L->pi[j+1]-L->pi[j];
  //int supWid = L->super[j+1]-L->super[j];
  for (int i = curCol+1; i < nxtCol+1; ++i) {
   li_ptr[i-1] = L->pi[j];
   colL[i]= colL[curCol] + (i-curCol)*colLength;
  }
 }
 li_ptr[n] = L->pi[L->nsuper];
 colL[n] = colL[n-1] + colLength;*/
 li_ptr = L->i_ptr;
 colL = L->p;
 valL = new double[L->xsize]();
 //delete []L->x;
 //delete []L->px;
 //delete []L->p;
 delete []L->pi;
 delete []L->i;
 delete []L->ColCount;
#if 0

 for (int j = 0; j < 10; ++j) {
        int curCol = L->super[j];
        int nxtCol = L->super[j+1];
        for (int k = li_ptr[curCol]; k < li_ptr[nxtCol]; ++k) {
            cout<<L->s[k]<<",";

        }
        cout<<"\n";
    }
#endif

 CSC *A1 = ptranspose(Amat,2,L->Perm,NULL,0,status);
 CSC *A2 = ptranspose(A1,2,NULL,NULL,0,status);
#ifndef VERIFY
 delete []L->Perm;
#endif
// int *map = new int[n]();
 std::vector<cholTime> timeArray;
 cholTime ct;
 //enableColdCache(1200,wasteFile);
 int iterNo=5;
 for (int k = 0; k < iterNo; ++k) {

  for (int i = 0; i < L->xsize; ++i) {
   valL[i]=0.0;
  }
  for (int i = 0; i < numThread + 4; ++i) {
   timingChol[i]=0.0;
  }
  bool retval=false;
 start = std::chrono::system_clock::now();



  retval=cholesky_left_par_05(n,A2->p,A2->i,A2->x,L->p,L->s,L->i_ptr,valL,
                       L->super,L->nsuper, timingChol,
#ifndef PRUNE
                       L->sParent,A1->p, A1->i, L->col2Sup,
#else
                       prunePtr,pruneSet,
#endif
                       nLevels,levelPtr,levelSet,
                       nPar, parPtr, partition,
                       chunk, numThread, maxSupWid+1,maxCol+1);




 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 duration2=elapsed_seconds.count();
 if(!retval)
     return -1;
#ifdef ANALYZE_THEORY
  for (int i = 0; i < L->xsize; ++i) {
   valL[i]=0.0;
  }
    pi=cholesky_Perf_Comm_Modeling(n,A2->p,A2->i,A2->x,L->p,L->s,L->i_ptr,valL,
                       L->super,L->nsuper, timingChol,
#ifndef PRUNE
                       L->sParent,A1->p, A1->i, L->col2Sup,
#else
    prunePtr,pruneSet,
#endif
                       nLevels,levelPtr,levelSet,
                       nPar, parPtr, partition,
                       chunk, numThread, maxSupWid+1,maxCol+1);
#endif
 ct.alltogether=duration2;
  ct.parallel=timingChol[0];
  ct.rootNodes=timingChol[1];
  for (int tt = 4; tt < numThread+4; ++tt) {
   ct.tArray.push_back(timingChol[tt]);
  }
  timeArray.push_back(ct);
  piArray.push_back(pi);
 }
#ifdef ANALYZE
 double costPerThread=0;
 computeComCost(n,A2->p,A2->i,A2->x,L->p,L->s,L->i_ptr,valL,
                L->super,L->nsuper, timingChol,
                L->sParent,A1->p, A1->i, L->col2Sup,
                nLevels,levelPtr,levelSet,
                nPar, parPtr, partition,
                totalIntraContribNNZs,totalInterCotribNNZs,numberOfIntraCore,
                numberOfInterCore,numOfEdgeCuts,costPerThread);
#endif
/* std::sort(timeArray.begin(),timeArray.end(),
 //          [&](cholTime a, cholTime b)-> bool{ return a.alltogether>b.alltogether;}
           retCMP
 );
 std::sort(piArray.begin(),piArray.end(),cmpPI);*/
 int mid = iterNo%2==1 ? (iterNo/2) : iterNo/2;
 mid = iterNo==1 ? 0 : mid;
 cout<<f1<<","<<numThread<<","<<chunk<<","<<costParam<<","
     <<levelParam<<","<<blasThreads<<","<<finalSeqNode<<",";
 cout<<timeArray[mid].alltogether<<","<<timeArray[mid].parallel<<","
     <<timeArray[mid].rootNodes<<",";
 cout<<durationSym<<","<<orderingTime<<",";
#if 0
 for (int l = 0; l < iterNo; ++l) {
  cout<<timeArray[l].alltogether<<",";
 }

#endif

 #ifdef ANALYZE

 std::cout<<totalIntraContribNNZs<<";"<<totalInterCotribNNZs<<";"<<numberOfIntraCore<<
   ";"<<numberOfInterCore<<";"<<numOfEdgeCuts<<";*;";

#endif
#ifdef ANALYZE_THEORY
 std::cout<<piArray[mid].theoreticalTime<<";";
 for (int k1 = 0; k1 < piArray[mid].allInfo.size(); ++k1) {
  for (int k2 = 0; k2 < piArray[mid].allInfo[k1].size(); ++k2) {
   std::cout<<piArray[mid].allInfo[k1][k2]<<";";
  }
  std::cout<<";;";
 }
#endif


/* std::cout<<"BLAS kernel time: \n";
 double sumBlas=0;
 for (int i = 4; i < 4 + numThread; ++i) {
  sumBlas += timeArray[mid].tArray[i];
 }
  std::cout<<sumBlas<<"\n";*/
 #ifdef VERIFY
 MKL_Domain_Set_Num_Threads(4,MKL_DOMAIN_BLAS);
 double resid [4], t, ta, tf, ts [3], tot, bnorm, xnorm, anorm, rnorm, fl,
   anz, axbnorm, rnorm2, resid2, rcond ;
 FILE *f ;
 cholmod_sparse *Ac1 ;
 cholmod_dense *X = NULL, *B, *W, *R = NULL ;
 double one [2], zero [2], minusone [2], beta [2], xlnz ;
 cholmod_common Common, *cm ;
 cholmod_factor *L1 ;
 double *Bx, *Rx, *Xx ;
 int i, isize, xsize, ordering, xtype, s, ss, lnz ;
 int trial, method, L_is_super ;
 int ver [3] ;

 ts[0] = 0.;
 ts[1] = 0.;
 ts[2] = 0.;

 /* ---------------------------------------------------------------------- */
 /* get the file containing the input matrix */
 /* ---------------------------------------------------------------------- */

 if ((f = fopen (argv [1], "r")) == NULL)
 {
  return -1;
 }


 /* ---------------------------------------------------------------------- */
 /* start CHOLMOD and set parameters */
 /* ---------------------------------------------------------------------- */

 cm = &Common ;
 cholmod_start (cm) ;
 //CHOLMOD_FUNCTION_DEFAULTS ;     /* just for testing (not required) */

 /* use default parameter settings, except for the error handler.  This
  * demo program terminates if an error occurs (out of memory, not positive
  * definite, ...).  It makes the demo program simpler (no need to check
  * CHOLMOD error conditions).  This non-default parameter setting has no
  * effect on performance. */
 // cm->error_handler = my_handler ;

 /* Note that CHOLMOD will do a supernodal LL' or a simplicial LDL' by
  * default, automatically selecting the latter if flop/nnz(L) < 40. */

 cm->final_asis = TRUE;
 cm->nmethods=1;
 cm->method[0].ordering = CHOLMOD_GIVEN ;
 //cm->postorder = TRUE;
 cm->supernodal_switch = 0.00000001;

 cm->nrelax [0] = nrelax[0] ;
 cm->nrelax [1] = nrelax[1] ;
 cm->nrelax [2] = nrelax[2] ;

 /* ---------------------------------------------------------------------- */
 /* create basic scalars */
 /* ---------------------------------------------------------------------- */

 zero [0] = 0 ;
 zero [1] = 0 ;
 one [0] = 1 ;
 one [1] = 0 ;
 minusone [0] = -1 ;
 minusone [1] = 0 ;
 beta [0] = 1e-6 ;
 beta [1] = 0 ;

 /* ---------------------------------------------------------------------- */
 /* read in a matrix */
 /* ---------------------------------------------------------------------- */

//    printf ("\n---------------------------------- cholmod_demo:\n") ;
 cholmod_version (ver) ;
//    printf ("cholmod version %d.%d.%d\n", ver [0], ver [1], ver [2]) ;
 SuiteSparse_version (ver) ;
//    printf ("SuiteSparse version %d.%d.%d\n", ver [0], ver [1], ver [2]) ;
 Ac1 = cholmod_read_sparse (f, cm) ;


 xtype = Ac1->xtype ;
 anorm = 1 ;
#ifndef NMATRIXOPS
 anorm = cholmod_norm_sparse (Ac1, 0, cm) ;
//    printf ("norm (A,inf) = %g\n", anorm) ;
//    printf ("norm (A,1)   = %g\n", cholmod_norm_sparse (Ac1, 1, cm)) ;
#endif
//    cholmod_print_sparse (Ac1, "A", cm) ;

 if (Ac1->nrow > Ac1->ncol)
 {
  /* Transpose A so that A'A+beta*I will be factorized instead */
  cholmod_sparse *C = cholmod_transpose (Ac1, 2, cm) ;
  cholmod_free_sparse (&Ac1, cm) ;
  Ac1 = C ;
//        printf ("transposing input matrix\n") ;
 }

 /* ---------------------------------------------------------------------- */
 /* create an arbitrary right-hand-side */
 /* ---------------------------------------------------------------------- */

 n = Ac1->nrow ;
 B = cholmod_zeros (n, 1, xtype, cm) ;
 Bx = (double *)B->x ;

#if GHS
 {
	/* b = A*ones(n,1), used by Gould, Hu, and Scott in their experiments */
	cholmod_dense *X0 ;
	X0 = cholmod_ones (A->ncol, 1, xtype, cm) ;
	cholmod_sdmult (A, 0, one, zero, X0, B, cm) ;
	cholmod_free_dense (&X0, cm) ;
    }
#else
 if (xtype == CHOLMOD_REAL)
 {
  /* real case */
  for (i = 0 ; i < n ; i++)
  {
   double x = n ;
   Bx [i] = 1 + i / x ;
  }
 }
 else
 {
  /* complex case */
  for (i = 0 ; i < n ; i++)
  {
   double x = n ;
   Bx [2*i  ] = 1 + i / x ;		/* real part of B(i) */
   Bx [2*i+1] = (x/2 - i) / (3*x) ;	/* imag part of B(i) */
  }
 }
#endif

//    cholmod_print_dense (B, "B", cm) ;
 bnorm = 1 ;
#ifndef NMATRIXOPS
 bnorm = cholmod_norm_dense (B, 0, cm) ;	/* max norm */
//    printf ("bnorm %g\n", bnorm) ;
#endif

 /* ---------------------------------------------------------------------- */
 /* analyze and factorize */
 /* ---------------------------------------------------------------------- */

 //t = CPUTIME ;
 //L1 = cholmod_analyze (Ac1, cm) ;
 L1 = cholmod_analyze_p (Ac1, L->Perm,NULL,NULL, cm) ;
 //ta = CPUTIME - t ;
 ta = MAX (ta, 0) ;

 printf ("Analyze flop, %g , lnz %g ,", cm->fl, cm->lnz) ;

 if (Ac1->stype == 0)
 {
  //    printf ("Factorizing A*A'+beta*I\n") ;
  //t = CPUTIME ;
  cholmod_factorize_p (Ac1, beta, NULL, 0, L1, cm) ;
  //tf = CPUTIME - t ;
  tf = MAX (tf, 0) ;
 }
 else
 {
  //  enableColdCache(1200,wasteFile);
  //start = std::chrono::system_clock::now();
  t = CPUTIME ;
  cholmod_factorize (Ac1, L1, cm) ;
  tf = CPUTIME - t ;
  tf = MAX (tf, 0) ;
  /*end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  duration3=elapsed_seconds.count();*/
 }


/* anz = cm->anz ;
 for (i = 0 ; i < CHOLMOD_MAXMETHODS ; i++)
 {
  fl = cm->method [i].fl ;
  xlnz = cm->method [i].lnz ;
  cm->method [i].fl = -1 ;
  cm->method [i].lnz = -1 ;
  ordering = cm->method [i].ordering ;
  if (fl >= 0)
  {
   printf ("Ordering: ") ;
   if (ordering == CHOLMOD_POSTORDERED) printf ("postordered ") ;
   if (ordering == CHOLMOD_NATURAL)     printf ("natural ") ;
   if (ordering == CHOLMOD_GIVEN)           printf ("user    ") ;
   if (ordering == CHOLMOD_AMD)             printf ("AMD     ") ;
   if (ordering == CHOLMOD_METIS)           printf ("METIS   ") ;
   if (ordering == CHOLMOD_NESDIS)      printf ("NESDIS  ") ;
   if (xlnz > 0)
   {
    printf ("fl/lnz %10.1f", fl / xlnz) ;
   }
   if (anz > 0)
   {
    printf ("  lnz/anz %10.1f", xlnz / anz) ;
   }
   printf ("\n") ;
  }
 }*/

#if 1
 int *Lpi = static_cast<int*> (L1->pi);
 int *Lsuper = static_cast<int*> (L1->super);
 int *Ls = static_cast<int*> (L1->s);
 int *Lpx = static_cast<int*> (L1->px);
 int cnt=0;
 double *Lx = static_cast<double*> (L1->x);
 int *LPerm = static_cast<int*> (L1->Perm);
 // cout<<"\n";
 for (int i = 0; i < L1->n; ++i) {
  assert(L->Perm[i]==LPerm[i]);
 }
// ASSERT(L1->nsuper == L->nsuper);
 for (int j = 0; j < L1->nsuper; ++j) {
  int colLength = Lpi[j+1]-Lpi[j];
  int supWid = Lsuper[j+1]-Lsuper[j];
  //cout<<"Supernode: "<<j<<" Len and Wid are: "
  //    <<colLength<<","<<supWid<<"\n";
//        for (int k = Lpi[j]; k < Lpi[j+1]; ++k) {
//            cout<<Ls[k]<<",";
//        }
//        cout<<"\n";
  for (int i = Lpx[j]; i < Lpx[j+1]; ++i) {
   if(Lx[i] - valL[i] > 0.001){
    cnt++;
    // cout<<"SN: "<<j<<";";
   }
   //cout<<Lx[i]<<",";
  }
  //  cout<<"\n";
 }

#endif
 /* solve with Bset will change L from simplicial to supernodal */
 L_is_super = L->is_super ;

 //cholmod_print_factor (L1, "L", cm) ;
 int NNZL = 0;//Lpx[L1->nsuper];
 cholmod_free_factor (&L1, cm) ;
 cholmod_free_dense (&X, cm) ;

 /* ---------------------------------------------------------------------- */
 /* free matrices and finish CHOLMOD */
 /* ---------------------------------------------------------------------- */

 cholmod_free_sparse (&Ac1, cm) ;
 cholmod_free_dense (&B, cm) ;
 cholmod_finish (cm) ;
#endif
 allocateAC(Amat,0,0,0,FALSE);
 allocateAC(A1,0,0,0,FALSE);
 allocateAC(A2,0,0,0,FALSE);
 //allocateLC(L,FALSE);
 //TODO HACK
 delete []L->super;
 delete []L->sParent;
 delete []L->s;
 delete []L->col2Sup;
 delete []L->p;
 delete []L->i_ptr;

#ifdef VERIFY
 cout<<tf <<","<< NNZL <<",";
 if(cnt>0)
  cout<<"#"<<cnt<<";";
#endif

#if DEBUG > 0
 for (int i = n-10; i < n; ++i) {
            std::cout<<i<<":\n";
            for (int m = colL[i],cnt=0; m < colL[i+1]; ++m, ++cnt) {
                if(!std::isfinite(valL[m])) {
                    std::cout << "Error in col "<< i;
                    return -1;
                }
                if(rowL[li_ptr[i]+cnt] >= i )
                    std::cout<<valL[m]<<",";
            }
            std::cout<<"\n";
        }
        std::cout<<"\n";
#endif
// delete []col2sup;
#ifdef PRUNE
 delete []prunePtr; delete []pruneSet;
#endif
 if(levelPtr!=NULL)
  delete []levelPtr;
 if(levelPtr!=NULL)
  delete []levelSet;
 if(parPtr!=NULL )
  delete []parPtr;
 if(partition!=NULL)
  delete []partition;
 //delete []contribs;
 //delete []map;
 delete []valL;
 //delete []colL;
 //delete []li_ptr;
 delete []timingChol;

 return 0;
}
