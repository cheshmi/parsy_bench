//
// Created by kazem on 7/26/17.
//

#ifndef CHOLOPENMP_LSPARSITY_H
#define CHOLOPENMP_LSPARSITY_H


#include "../common/def.h"
#include "../common/PostOrder.h"
#include "../common/Etree.h"
#include "ColumnCount.h"
#include "Inspection_BlockC.h"
#include "amd.h"


#ifdef SCOTCH
#include "scotch.h"
#endif

//From CHOLMOD
#include "amd.h"
#include "metis.h"
#include "InspectionLevel.h"

#include "InspectionLevel_06.h"
#include "performanceModel.h"

int permute_matrices
  (
    /* ---- input ---- */
    CSC *A,	/* matrix to permute */
    int ordering,	/* ordering method used */
    int *Perm,		/* fill-reducing permutation */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int do_rowcolcounts,/* if TRUE, compute both S and F.  If FALSE, only
			 * S is needed for the symmetric case, and only F for
			 * the unsymmetric case */

    /* ---- output --- */
    CSC **A1_handle,	    /* see comments below for A1, A2, S, F */
    CSC **A2_handle,
    CSC **S_handle,
    CSC **F_handle,
    /* --------------- */
    //cholmod_common *Common
    int &status
  )
{
 CSC *A1, *A2, *S, *F ;

 *A1_handle = NULL ;
 *A2_handle = NULL ;
 *S_handle = NULL ;
 *F_handle = NULL ;
 A1 = NULL ;
 A2 = NULL ;

 if (ordering == CHOLMOD_NATURAL)
 {

  /* ------------------------------------------------------------------ */
  /* natural ordering of A */
  /* ------------------------------------------------------------------ */

  if (A->stype < 0)
  {
   /* symmetric lower case: A already in lower form, so S=A' */
   /* workspace: Iwork (nrow) */
   A2 = ptranspose(A, 0, NULL, NULL, 0, status) ;
   F = A ;
   S = A2 ;
  }
  else if (A->stype > 0)
  {
   /* symmetric upper case: F = pattern of triu (A)', S = A */
   /* workspace: Iwork (nrow) */
   if (do_rowcolcounts)
   {
    /* F not needed for symmetric case if do_rowcolcounts FALSE */
    A1 = ptranspose(A, 0, NULL, fset, fsize, status) ;
   }
   F = A1 ;
   S = A ;
  }
  else
  {
   /* unsymmetric case: F = pattern of A (:,f)',  S = A */
   /* workspace: Iwork (nrow if no fset, MAX(nrow,ncol) if fset) */
   A1 = ptranspose(A, 0, NULL, fset, fsize, status) ;
   F = A1 ;
   S = A ;
  }

 }
 else
 {

  /* ------------------------------------------------------------------ */
  /* A is permuted */
  /* ------------------------------------------------------------------ */

  if (A->stype < 0)
  {
   /* symmetric lower case: S = tril (A (p,p))' and F = S' */
   /* workspace: Iwork (2*nrow) */
   A2 = ptranspose(A, 0, Perm, NULL, 0, status) ;
   S = A2 ;
   /* workspace: Iwork (nrow) */
   if (do_rowcolcounts)
   {
    /* F not needed for symmetric case if do_rowcolcounts FALSE */
    A1 = ptranspose(A2, 0, NULL, NULL, 0, status) ;
   }
   F = A1 ;
  }
  else if (A->stype > 0)
  {
   /* symmetric upper case: F = triu (A (p,p))' and S = F' */
   /* workspace: Iwork (2*nrow) */
   A1 = ptranspose(A, 0, Perm, NULL, 0, status) ;
   F = A1 ;
   /* workspace: Iwork (nrow) */
   A2 = ptranspose(A1, 0, NULL, NULL, 0, status) ;
   S = A2 ;
  }
  else
  {
   /* unsymmetric case:     F = A (p,f)'         and S = F' */
   /* workspace: Iwork (nrow if no fset, MAX(nrow,ncol) if fset) */
   A1 = ptranspose(A, 0, Perm, fset, fsize, status) ;
   F = A1 ;
   if (do_rowcolcounts)
   {
    /* S not needed for unsymmetric case if do_rowcolcounts FALSE */
    /* workspace: Iwork (nrow) */
    A2 = ptranspose(A1, 0, NULL, NULL, 0, status) ;
   }
   S = A2 ;
  }
 }

 /* If any cholmod_*transpose fails, one or more matrices will be NULL */
 *A1_handle = A1 ;
 *A2_handle = A2 ;
 *S_handle = S ;
 *F_handle = F ;
 return (TRUE) ;
}

/* ========================================================================== */
/* === cholmod_analyze_ordering ============================================= */
/* ========================================================================== */

/* Given a matrix A and its fill-reducing permutation, compute the elimination
 * tree, its (non-weighted) postordering, and the number of nonzeros in each
 * column of L.  Also computes the flop count, the total nonzeros in L, and
 * the nonzeros in A (Common->fl, Common->lnz, and Common->anz).
 *
 * The column counts of L, flop count, and other statistics from
 * cholmod_rowcolcounts are not computed if ColCount is NULL.
 *
 * workspace: Iwork (2*nrow if symmetric, 2*nrow+ncol if unsymmetric),
 *	Flag (nrow), Head (nrow+1)
 */

int analyze_ordering
  (
    /* ---- input ---- */
    CSC *A,	/* matrix to analyze */
    int ordering,	/* ordering method used */
    int *Perm,		/* size n, fill-reducing permutation to analyze */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    int *Parent,	/* size n, elimination tree */
    int *Post,		/* size n, postordering of elimination tree */
    int *ColCount,	/* size n, nnz in each column of L */

    /* ---- workspace  */
    int *First,		/* size n workspace for cholmod_postorder */
    int *Level,		/* size n workspace for cholmod_postorder */
    /* --------------- */
    // cholmod_common *Common
    int &status
  )
{
 CSC *A1, *A2, *S, *F ;
 int n, ok, do_rowcolcounts ;
 int fl, aatfl, lnz;

 /* check inputs */
 //   RETURN_IF_NULL_COMMON (FALSE) ;
 //   RETURN_IF_NULL (A, FALSE) ;

 n = A->nrow ;

 do_rowcolcounts = (ColCount != NULL) ;

 /* permute A according to Perm and fset */
 ok = permute_matrices (A, ordering, Perm, fset, fsize, do_rowcolcounts,
                        &A1, &A2, &S, &F, status) ;

 /* find etree of S (symmetric upper/lower case) or F (unsym case) */
 /* workspace: symmmetric: Iwork (nrow), unsym: Iwork (nrow+ncol) */
 ok = ok && etreeC (A->stype ? S:F, Parent, status) ;
#if 0
 int cc=0;
 for (int i = 0; i < n; ++i) {
        //printf("%d :", Parent[i]);
  if(Parent[i]<0)
   cc++;
 }
 printf("\n");
#endif

 /* postorder the etree (required by cholmod_rowcolcounts) */
 /* workspace: Iwork (2*nrow) */
 ok = ok && postOrderC(Parent, n, NULL, Post, status) == n ;

 /* cholmod_postorder doesn't set Common->status if it returns < n */
 //status = (!ok && status ) ?
 // CHOLMOD_INVALID : Common->status ;

 /* analyze LL'=S or SS' or S(:,f)*S(:,f)' */
 /* workspace:
  *	if symmetric:   Flag (nrow), Iwork (2*nrow)
  *	if unsymmetric: Flag (nrow), Iwork (2*nrow+ncol), Head (nrow+1)
  */
 if (do_rowcolcounts)
 {
  ok = ok && rowcolcounts(A->stype ? F:S, fset, fsize, Parent,
                          Post, NULL, ColCount, First, Level,
                          fl,
                          aatfl,
                          lnz, status) ;
 }

 /* free temporary matrices and return result */
 // CHOLMOD(free_sparse) (&A1, Common) ;
 // CHOLMOD(free_sparse) (&A2, Common) ;
 allocateAC(A1,0,0,0,FALSE);
 delete A1;
 allocateAC(A2,0,0,0,FALSE);
 delete A2;
 return (ok) ;
}


/* ========================================================================== */
/* === cholmod_analyze_p2 =================================================== */
/* ========================================================================== */

/* Ordering and analysis for sparse Cholesky or sparse QR.  */

BCSC* analyze_p2
  (
    /* ---- input ---- */
    int for_whom,       /* FOR_SPQR     (0): for SPQR but not GPU-accelerated
                           FOR_CHOLESKY (1): for Cholesky (GPU or not)
                           FOR_SPQRGPU  (2): for SPQR with GPU acceleration */
    CSC *A,	/* matrix to order and analyze */
    int *UserPerm,	/* user-provided permutation, size A->nrow */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    int *nrelax,
    double *zrelax,
    size_t fsize,	/* size of fset */
    int* &prunePtr,
    int* &pruneSet,
    int &levelNo,
    int* &levelPtr,
    int* &levelSet,
    int &parNo,
    int* &parPtr,
    int* &partition,
    int costParam,
    int levelParam, //min level distance
    int finalSeqNodes, //number of level in final partition
    /* --------------- */
    int status,
    //out
    int &maxSupWid,
    int &maxCol,
    double &orderingTime,
    size_t *inPerm=NULL
  )
{
 double lnz_best ;
 int *First, *Level, *Work4n, *Cmember, *CParent, *ColCount, *Lperm, *Parent,
   *Post, *Perm, *Lparent, *Lcolcount ;
 BCSC *L ;
 int n, ncol, k, ordering, method, nmethods,  default_strategy, uncol,
   skip_analysis, skip_best ;
 int amd_backup ;
 size_t s ;
 int ok = TRUE ;
 std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;
 /* ---------------------------------------------------------------------- */
 /* check inputs */
 /* ---------------------------------------------------------------------- */

/*    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    Common->status = CHOLMOD_OK ;
    status = CHOLMOD_OK ;
    Common->selected = EMPTY ;
    Common->called_nd = FALSE ;*/

 /* ---------------------------------------------------------------------- */
 /* get inputs */
 /* ---------------------------------------------------------------------- */

 n = A->nrow ;
 ncol = A->ncol ;
 uncol = (A->stype == 0) ? (A->ncol) : 0 ;

 /* ---------------------------------------------------------------------- */
 /* set the default strategy */
 /* ---------------------------------------------------------------------- */

 lnz_best = (double) EMPTY ;
 skip_best = FALSE ;
 //  nmethods = MIN (Common->nmethods, CHOLMOD_MAXMETHODS) ;
 // nmethods = MAX (0, nmethods) ;



/*    default_strategy = (nmethods == 0) ;
    if (default_strategy)
    {
        *//* default strategy: try UserPerm, if given.  Try AMD for A, or AMD
         * to order A*A'.  Try METIS for the symmetric case only if AMD reports
             * a high degree of fill-in and flop count.  METIS is not tried if the
             * Partition Module isn't installed.   If Common->default_nesdis is
             * TRUE, then NESDIS is used as the 3rd ordering instead. *//*
        Common->method [0].ordering = CHOLMOD_GIVEN ;*//* skip if UserPerm NULL *//*
        Common->method [1].ordering = CHOLMOD_AMD ;
        Common->method [2].ordering =
                (Common->default_nesdis ? CHOLMOD_NESDIS : CHOLMOD_METIS) ;
        amd_backup = FALSE ;
#ifndef NPARTITION
        nmethods = 3 ;
#else
        nmethods = 2 ;
#endif
    }
    else
    {
        *//* If only METIS and NESDIS are selected, or if 2 or more methods are
         * being tried, then enable AMD backup *//*
        amd_backup = (nmethods > 1) || (nmethods == 1 &&
                                        (Common->method [0].ordering == CHOLMOD_METIS ||
                                         Common->method [0].ordering == CHOLMOD_NESDIS)) ;
    }*/

#ifdef NSUPERNODAL
 /* CHOLMOD Supernodal module not installed, just do simplicial analysis */
    Common->supernodal = CHOLMOD_SIMPLICIAL ;
#endif

 /* ---------------------------------------------------------------------- */
 /* allocate workspace */
 /* ---------------------------------------------------------------------- */

 /* Note: enough space needs to be allocated here so that routines called by
  * cholmod_analyze do not reallocate the space.
  */

 /* s = 6*n + uncol */
 /* s = CHOLMOD(mult_size_t) (n, 6, &ok) ;
  s = CHOLMOD(add_size_t) (s, uncol, &ok) ;
  if (!ok)
  {
//        ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
      return (NULL) ;
  }

  CHOLMOD(allocate_work) (n, s, 0, Common) ;
  if (Common->status < CHOLMOD_OK)
  {
      return (NULL) ;	    *//* out of memory *//*
    }*/
 //   ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

 /* ensure that subsequent routines, called by cholmod_analyze, do not
  * reallocate any workspace.  This is set back to FALSE in the
  * FREE_WORKSPACE_AND_RETURN macro, which is the only way this function
  * returns to its caller. */
 // Common->no_workspace_reallocate = TRUE ;

 /* Use the last 4*n int's in Iwork for Parent, First, Level, and Post, since
  * other CHOLMOD routines will use the first 2n+uncol space.  The ordering
  * routines (cholmod_amd, cholmod_colamd, cholmod_ccolamd, cholmod_metis)
  * are an exception.  They can use all 6n + ncol space, since the contents
  * of Parent, First, Level, and Post are not needed across calls to those
  * routines. */
 //s = 6*n + uncol;
 s=4*n;
 Work4n = new int[s]() ;
 // Work4n += 2*((size_t) n) + uncol ;
 Parent = Work4n ;
 First  = Work4n + n ;
 Level  = Work4n + 2*((size_t) n) ;
 Post   = Work4n + 3*((size_t) n) ;

 /* note that this assignment means that cholmod_nested_dissection,
  * cholmod_ccolamd, and cholmod_camd can use only the first 4n+uncol
  * space in Common->Iwork */
 Cmember = Post ;
 CParent = Level ;

 /* ---------------------------------------------------------------------- */
 /* allocate more workspace, and an empty simplicial symbolic factor */
 /* ---------------------------------------------------------------------- */

 //L = CHOLMOD(allocate_factor) (n, Common) ;//TODO
 L = new BCSC;
 L->Perm     = new int[n]();
 L->ColCount = new int[n]();

 Lparent  = new int[n]();//CHOLMOD(malloc) (n, sizeof (int), Common) ;
 //  Perm     = new int[n]();//CHOLMOD(malloc) (n, sizeof (int), Common) ;
 //  ColCount = new int[n]();//CHOLMOD(malloc) (n, sizeof (int), Common) ;
 /*if (Common->status < CHOLMOD_OK)
 {
     *//* out of memory *//*
        FREE_WORKSPACE_AND_RETURN ;
    }*/
 Lperm = L->Perm ;
 Lcolcount = L->ColCount ;
 //size_t *colL = new size_t[n](); //colL = L->p;
 //L->p = colL;
 //size_t *Li_ptr = L->i_ptr;
 //Common->anz = EMPTY ;

 /* ---------------------------------------------------------------------- */
 /* try all the requested ordering options and backup to AMD if needed */
 /* ---------------------------------------------------------------------- */

 /*Common->method [i].prune_dense2 = -1 ;	*//* COLAMD dense row control *//*
    320  	Common->method [i].aggressive = TRUE ;	*//* aggressive absorption *//*
    321  	Common->method [i].order_for_lu = FALSE ;*//* order for Cholesky not LU */
 start = std::chrono::system_clock::now();
#ifdef GIVEN
 //pastix_data_t **pastix_data;
 L->ordering = CHOLMOD_METIS;
 for (int l = 0; l < A->nrow; ++l) {
  Lperm[l] = inPerm[l];
 }

#elif SCOTCH
 L->ordering = CHOLMOD_METIS;
 CSC *ATrans;
 unsigned long nnzFull = A->nzmax*2;//Symmetric case
 ATrans = ptranspose(A, 0, NULL, NULL, 0, status) ;

 SCOTCH_Num baseVal=0;
 SCOTCH_Num          vertnbr;                    /* Number of vertices */
 SCOTCH_Graph        grafdat;                    /* Source graph       */
 SCOTCH_Ordering     ordedat;                    /* Graph ordering     */
 SCOTCH_Strat        stradat;                    /* Ordering strategy  */
 SCOTCH_Num          straval;
 //Making the graph for passing it to metis, it should have
 //both upper and lower parts
 //allocateAC(AFull,ncol,nnzFull,0,TRUE);
 SCOTCH_Num *permtab = new SCOTCH_Num[ncol]();                    /* Permutation array  */
 SCOTCH_Num *AFullp = new SCOTCH_Num[ncol+1]();
 SCOTCH_Num *AFulli = new SCOTCH_Num[nnzFull]();
 SCOTCH_Num ncolIDXT = ncol;

 AFullp[0]=0;
 for (int i = 0; i < ncol; ++i) {
  int nnzOfCurCol = ATrans->p[i+1]-ATrans->p[i]-1;
  nnzOfCurCol += A->p[i+1]-A->p[i]-1;
  AFullp[i+1] =(long int) AFullp[i]+nnzOfCurCol;
  //copying Upper part, ignoring diagonal since it is in L already
  int base=AFullp[i];
  for (int j = ATrans->p[i],k=0; j < ATrans->p[i+1]-1; ++j,++k) {
   AFulli[base+k] =(long int) ATrans->i[j];
  }
  //copying L part
  base+=ATrans->p[i+1]-ATrans->p[i]-1;
  for (int j = A->p[i]+1,k=0; j < A->p[i+1]; ++j, ++k) {
   AFulli[base+k] =(long int) A->i[j];
  }
 }
 if(SCOTCH_stratInit (&stradat))
  return NULL;
 char straStr[550];
 sprintf(straStr, SCOTCH_STRAT_DIRECT);
 if(SCOTCH_stratGraphOrder(&stradat,straStr))
  return NULL;
 nnzFull=AFullp[ncol];
 if(SCOTCH_graphBuild (&grafdat,baseVal,ncol,AFullp,NULL,NULL,NULL,
                   nnzFull,AFulli,NULL))
  return NULL;
 if(SCOTCH_graphCheck(&grafdat))
  return NULL;
 /*if(SCOTCH_graphOrderList(&grafdat, ncol, NULL, &stradat,
                          permtab, NULL, NULL, NULL, NULL)) { *//* Create ordering *//*
  return NULL;
 }*/
 /*if(SCOTCH_graphOrderInit    (&grafdat, &ordedat, permtab, NULL, NULL, NULL, NULL)) {  Create ordering
  return NULL;
 }
 if(SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat)) {  Perform ordering
  return NULL;
 }*/

 if(SCOTCH_graphOrderInit    (&grafdat, &ordedat, permtab, NULL, NULL, NULL, NULL)) {  // Create ordering
  return NULL;
 }
 if(SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat)) { // Perform ordering
  return NULL;
 }

 if (SCOTCH_graphOrderCheck (&grafdat, &ordedat) != 0)
  return NULL;
 for (int i = 0; i < ncol; ++i) {
  Lperm[i]=permtab[i];
  //std::cout<<Lperm[i];
 }
 SCOTCH_graphOrderExit (&grafdat, &ordedat);
 SCOTCH_stratExit      (&stradat);
 SCOTCH_graphExit      (&grafdat);

 allocateAC(ATrans,ATrans->nrow,ATrans->nzmax,ATrans->stype,false);
 delete [] AFullp;
 delete [] AFulli;
 delete [] permtab;

#elif METIS
 L->ordering = CHOLMOD_METIS;
 CSC *ATrans;
 unsigned long nnzFull = A->nzmax*2;//Symmetric case
 ATrans = ptranspose(A, 0, NULL, NULL, 0, status) ;
#if 0
 for (int i = 0; i < ncol; ++i) {
  for (int j = A->p[i]; j < A->p[i+1]; ++j) {
   std::cout<<A->i[j]<<";";
  }
  std::cout<<"\n";
 }
 std::cout<<"---\n";
 for (int i = 0; i < ncol; ++i) {
  for (int j = ATrans->p[i]; j < ATrans->p[i+1]; ++j) {
   std::cout<<ATrans->i[j]<<";";
  }
  std::cout<<"\n";
 }
 std::cout<<"==\n";
#endif

 //Making the graph for passing it to metis, it should have
 //both upper and lower parts
 //allocateAC(AFull,ncol,nnzFull,0,TRUE);
 idx_t options1[METIS_NOPTIONS];
 METIS_SetDefaultOptions(options1);

 idx_t *AFullp = new idx_t[ncol+1]();
 idx_t *AFulli = new idx_t[nnzFull]();
 idx_t ncolIDXT = ncol;
 idx_t *weigt = new idx_t[ncol];
 idx_t *LpermIDX = new idx_t[ncol];
 idx_t *ILpermIDX = new idx_t[ncol];
 for (int i = 0; i < ncol; ++i) {
  LpermIDX[i]=0;ILpermIDX[i]=0;weigt[i]=1;
 }
 AFullp[0]=0;
 for (int i = 0; i < ncol; ++i) {
  int nnzOfCurCol = ATrans->p[i+1]-ATrans->p[i]-1;
  nnzOfCurCol += A->p[i+1]-A->p[i]-1;
  AFullp[i+1] =(long int) AFullp[i]+nnzOfCurCol;
  //copying Upper part, ignoring diagonal since it is in L already
  int base=AFullp[i];
  for (int j = ATrans->p[i],k=0; j < ATrans->p[i+1]-1; ++j,++k) {
   AFulli[base+k] =(long int) ATrans->i[j];
  }
  //copying L part
  base+=ATrans->p[i+1]-ATrans->p[i]-1;
  for (int j = A->p[i]+1,k=0; j < A->p[i+1]; ++j, ++k) {
   AFulli[base+k] =(long int) A->i[j];
  }
 }
#if 0
 for (int i = 0; i < ncol; ++i) {
  for (int j = AFullp[i]; j < AFullp[i+1]; ++j) {
   std::cout<<AFulli[j]<<";";
  }
  std::cout<<"\n";
 }
 std::cout<<"\n";
#endif

 int retMet= METIS_NodeND(&ncolIDXT,AFullp,AFulli,NULL,options1,
              LpermIDX,ILpermIDX);
 assert(retMet==METIS_OK);
 if(retMet!=METIS_OK){
  std::cout<<" "<<retMet<<"\n";
  exit(10);
 }
 for (int i = 0; i < ncol; ++i) {
  Lperm[i]=LpermIDX[i];
  //std::cout<<Lperm[i];
 }
 allocateAC(ATrans,ATrans->nrow,ATrans->nzmax,ATrans->stype,false);
 METIS_Free(AFullp);
 METIS_Free(AFulli);
 METIS_Free(weigt);
 METIS_Free(LpermIDX) ;
 METIS_Free(ILpermIDX);
#else
 double info[20]={0};
 double Control[2];
 Control [0] = 10; //TODO check later //AMD_Dense
 Control [1] = TRUE; //AMD_AGGRESSIVE
 L->ordering = CHOLMOD_AMD;
 amd_order(ncol,A->p,A->i,Lperm,NULL,info);
#endif
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 orderingTime = elapsed_seconds.count();
//printf("ddddd %f ddd \n",orderingTime);
#ifdef VERIFY
 auto checkOrder = new bool[ncol]();
 for (int i = 0; i < ncol; ++i) checkOrder[i]=false;
 for (int i = 0; i < ncol; ++i) {
  checkOrder[Lperm[i]]=true;
 }
 for (int i = 0; i < ncol; ++i) {
  assert(checkOrder[i]==true);
 }
 delete checkOrder;
#endif
#if 0
 for (int i = 0; i < n; ++i) {
        printf("%d, ", Lperm[i]);
    }
    printf("\n");
#endif
//////// Can be added later
 skip_best=TRUE;

 /* ---------------------------------------------------------------------- */
 /* do the analysis for AMD, if skipped */
 /* ---------------------------------------------------------------------- */

/*    Common->fl  = Common->method [Common->selected].fl  ;
    Common->lnz = Common->method [Common->selected].lnz ;
    ASSERT (Common->lnz >= 0) ;*/

 if (skip_best)
 {
  if (!analyze_ordering (A, L->ordering, Lperm, fset, fsize,
                         Lparent, Post, Lcolcount, First,
                         Level, status))
  {
   /* out of memory, or method failed */
   //FREE_WORKSPACE_AND_RETURN ;
   return FALSE;
  }
 }

 /* ---------------------------------------------------------------------- */
 /* postorder the etree, weighted by the column counts */
 /* ---------------------------------------------------------------------- */

 //if (Common->postorder)
 if(true)//TODO
 {
  /* combine the fill-reducing ordering with the weighted postorder */
  /* workspace: Iwork (2*nrow) */
  if (postOrderC(Lparent, n, Lcolcount, Post,status) == n)
  {
   /* use First and Level as workspace [ */
   int *Wi = First, *InvPost = Level ;
   int newchild, oldchild, newparent, oldparent ;

   for (k = 0 ; k < n ; k++)
   {
    Wi [k] = Lperm [Post [k]] ;
   }
   for (k = 0 ; k < n ; k++)
   {
    Lperm [k] = Wi [k] ;
   }

   for (k = 0 ; k < n ; k++)
   {
    Wi [k] = Lcolcount [Post [k]] ;
   }
   for (k = 0 ; k < n ; k++)
   {
    Lcolcount [k] = Wi [k] ;
   }
   for (k = 0 ; k < n ; k++)
   {
    InvPost [Post [k]] = k ;
   }


   /* updated Lparent needed only for supernodal case */
   for (newchild = 0 ; newchild < n ; newchild++)
   {
    oldchild = Post [newchild] ;
    oldparent = Lparent [oldchild] ;
    newparent = (oldparent == EMPTY) ? EMPTY : InvPost [oldparent] ;
    Wi [newchild] = newparent ;
   }
   for (k = 0 ; k < n ; k++){
    Lparent [k] = Wi [k] ;
   }
   /* done using Iwork as workspace ] */

   /* L is now postordered, no longer in natural ordering */
   if (L->ordering == CHOLMOD_NATURAL)
   {
    L->ordering = CHOLMOD_POSTORDERED ;
   }
  }
 }


 /* ---------------------------------------------------------------------- */
 /* supernodal analysis, if requested or if selected automatically */
 /* ---------------------------------------------------------------------- */

#ifndef NSUPERNODAL
 /*  if (Common->supernodal > CHOLMOD_AUTO
       || (Common->supernodal == CHOLMOD_AUTO &&
           Common->lnz > 0 &&
           (Common->fl / Common->lnz) >= Common->supernodal_switch))*/
 int supernodal = TRUE;
 int *Sparent = new int[n]();
 if(supernodal)
 {
  CSC *S, *F, *A2, *A1 ;

  permute_matrices (A, L->ordering, Lperm, fset, fsize, TRUE,
                    &A1, &A2, &S, &F, status) ;

  /* workspace: Flag (nrow), Head (nrow), Iwork (5*nrow) */
  /*int nrelax[3] = {4,16,48};//TODO
  double zrelax[3] = {0.8,0.1,0.05};*/
  super_symbolic2(for_whom, S, F, Lparent, L,nrelax,zrelax,Sparent, status) ;
#ifdef PRUNE
  prunePtr = new int[L->nsuper+1]();
  pruneSet = new int[L->ssize];
#endif
  int *col2Sup = new int[F->nrow];

  // Computing col2sup
  maxSupWid=0;
  for (int i = 0; i < L->nsuper; ++i) {
   int k1 = L->super[i];
   int k2 = L->super[i+1];
   if(maxSupWid<(k2-k1))
    maxSupWid=(k2-k1);
   ASSERT(k1 <= F->nrow && k2 <= F->nrow);
   for (int j = k1; j < k2; ++j) {
    col2Sup[j] = i;
   }
  }
  //Computing the row pointer of L
  int colLength=0;
  maxCol=0;
  for (int j = 0; j < L->nsuper; ++j) {
   int curCol = L->super[j];
   int nxtCol = L->super[j+1];
   colLength = L->pi[j+1]-L->pi[j];
   if(colLength>maxCol)
    maxCol=colLength;
   //int supWid = L->super[j+1]-L->super[j];
   for (int i = curCol+1; i < nxtCol+1; ++i) {
    L->i_ptr[i-1] = L->pi[j];
    L->p[i]= L->p[curCol] + (i-curCol)*colLength;
   }
  }
  L->i_ptr[n] = L->pi[L->nsuper];
  L->p[n] = L->p[n-1] + colLength;


  L->col2Sup=col2Sup;
  L->sParent=Sparent;
  //Computing the node cost

  start = std::chrono::system_clock::now();
  double *nodeCost = new double[L->nsuper];
  int *xi = new int[2*F->ncol]();
  for (int s = 1; s <= L->nsuper; ++s) {
   int curCol = s != 0 ? L->super[s - 1] : 0;
   int nxtCol = L->super[s];
   assert(s-1 < L->nsuper && s-1 >= 0 );
   nodeCost[s-1]=computeCostperBlock(F->ncol, curCol, nxtCol, Lparent,
     F->p, F->i, col2Sup, L->super, L->s, L->i_ptr,xi);
  }
  delete []xi;
  //PRINT1 (("status %d\n", Common->status)) ;
#ifdef PRUNE
  getBlockedPruneSet(L->nsuper,S->p,S->i,col2Sup,Sparent,L->super,
                     prunePtr,pruneSet);
#endif

  //levelNo = getLevelSet(L->nsuper,Sparent,levelPtr,levelSet);

  getCoarseLevelSet_6(L->nsuper,Sparent,L->super,
                    levelNo, levelPtr, levelSet,//Partition level set
                    parNo, parPtr, partition,
                    costParam,levelParam,finalSeqNodes,nodeCost);



  assert(levelNo<=L->nsuper);

  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  double duration1 = elapsed_seconds.count();
  //printf("time lbc: %f ;",duration1);
  //CHOLMOD(free_sparse) (&ATrans, Common) ;
  //CHOLMOD(free_sparse) (&A2, Common) ;
  //delete []col2Sup;
  allocateAC(A1,0,0,0,FALSE);
  allocateAC(A2,0,0,0,FALSE);
  delete []nodeCost;
 }
#endif

 /* ---------------------------------------------------------------------- */
 /* free temporary matrices and workspace, and return result L */
 /* ---------------------------------------------------------------------- */

 //FREE_WORKSPACE_AND_RETURN ;
 //delete []Sparent;
 delete []Work4n;
 return L;
}
//#endif
#endif //CHOLOPENMP_LSPARSITY_H
