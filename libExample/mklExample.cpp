//
// Created by kazem on 3/10/18.
//

/*******************************************************************************
* Copyright 2013-2017 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
*   Content : Intel(R) MKL SpMV Format Prototype Package C native example
*
********************************************************************************
*/
/*
!
! Consider the matrix A (see 'Sparse Storage Formats for Sparse BLAS Level 2
! and Level 3 in the Intel(R) MKL Reference Manual')
!
!                 |   1   -1   7   -3   0    0   |
!                 |  -2    5   9   -2   0    0   |
!   A    =        |   0    0   4    6   4   -9   |
!                 |   0    0   2    7   1    2   |
!                 |   3    8   0    0   5   -2   |
!                 |   4    1   0    0   5    7   |
!
!  The matrix A is represented in a zero-based compressed sparse row (CSR) storage
!  scheme with three arrays (see 'Sparse Matrix Storage Schemes' in the
!  Intel(R) MKL Reference Manual) as follows:
!
!         values  =  ( 1 -1 -2 5 7 -3 9 -2 4 6 2 7 4 -9 1 2 3 8 4 1 5 -2 5 7 )
!         columns =  ( 0  1  1  2  0  2 )
!         rowIndex = ( 0  2  4  6 )
!
!*******************************************************************************
*/
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "mkl_spblas.h"

int main() {
    //*******************************************************************************
    //     Declaration and initialization of parameters for sparse representation of
    //     the matrix A in the block CSR format:
    //*******************************************************************************
#define M 3
#define NRHS 2
#define NNZ 24
#define NNZB 6
#define LB 2
    //*******************************************************************************
    //    Sparse representation of the matrix A
    //*******************************************************************************
    double csrVal[NNZ]    = {  1.0, -1.0, -2.0, 5.0, 7.0, -3.0, 9.0, -2.0,
                               4.0, 6.0, 2.0, 7.0, 4.0, -9.0, 1.0, 2.0,
                               3.0, 8.0, 4.0, 1.0, 5.0, -2.0, 5.0, 7.0, };
    MKL_INT    csrColInd[NNZB] = { 0, 1, 1, 2, 0, 2 };
    MKL_INT    csrRowPtr[M+1] = { 0, 2, 4, 6 };
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;
    // // Structure with sparse matrix stored in CSR format
    sparse_matrix_t       bsrA;
    //*******************************************************************************
    //    Declaration of local variables:
    //*******************************************************************************
    double x_m[M*LB*NRHS]  = { 1.0, 5.0, 3.0, 4.0, 2.0, 6.0, 2.0, 10.0, 6.0, 8.0, 4.0, 12.0};
    double y_m[M*LB*NRHS]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double x_v[M*LB]  = { 3.0, 2.0, 5.0, 4.0, 1.0, 6.0};
    double y_v[M*LB]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double tmp_v[M*LB]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double alpha = 1.0, beta = 0.0;
    MKL_INT    i;

    printf( "\n EXAMPLE PROGRAM FOR CSR format routines from package\n" );
    printf( "-------------------------------------------------------\n" );

    // Create handle with matrix stored in CSR format
    mkl_sparse_d_create_bsr ( &bsrA, SPARSE_INDEX_BASE_ZERO,
                                     SPARSE_LAYOUT_ROW_MAJOR,
                                     M,  // number of rows
                                     M,  // number of cols
                                     LB, // block size
                                     csrRowPtr,
                                     csrRowPtr+1,
                                     csrColInd,
                                     csrVal );

    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_optimize ( bsrA );

//  Task 1: Obtain matrix-matrix multiply (L+D)' *x_v --> y_v
//          and solve triangular system   (L+D)' *tmp_v = y_v
//          Array tmp_v must be equal to the array x_v
    printf( "                                  \n" );
    printf( "   Task 1:                        \n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mv \n" );
    printf( "   WITH TRIANGULAR SPARSE MATRIX  \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_NON_TRANSPOSE     \n" );
    printf( "   Input vector                   \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", x_v[i] );
    }

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descrA.mode = SPARSE_FILL_MODE_LOWER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;

    mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE,
                      alpha,
                      bsrA,
                      descrA,
                      x_v,
                      beta,
                      y_v );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", y_v[i] );
    }

    printf("   Solve triangular system   \n");
    printf("   with obtained             \n");
    printf("   right hand side           \n");

    mkl_sparse_d_trsv ( SPARSE_OPERATION_NON_TRANSPOSE,
                      alpha,
                      bsrA,
                      descrA,
                      y_v,
                      tmp_v );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_trsv \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", tmp_v[i] );
    }
    printf( "-------------------------------------------------------\n" );

//  Task 2: Obtain matrix-matrix multiply (U+I)' *x_v --> y_v
//          and solve triangular system   (U+I)' *tmp_v = y_v
//          Array tmp_v must be equal to the array x_v
    printf( "                                  \n" );
    printf( "   Task 2:                        \n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mv \n" );
    printf( "   WITH TRIANGULAR SPARSE MATRIX  \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_TRANSPOSE     \n" );
    printf( "   Input vector                   \n" );
    // Release matrix handle and deallocate matrix
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", x_v[i] );
    }

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_UNIT;

    mkl_sparse_d_mv ( SPARSE_OPERATION_TRANSPOSE,
                      alpha,
                      bsrA,
                      descrA,
                      x_v,
                      beta,
                      y_v );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", y_v[i] );
    }

    printf("   Solve triangular system   \n");
    printf("   with obtained             \n");
    printf("   right hand side           \n");

    mkl_sparse_d_trsv ( SPARSE_OPERATION_TRANSPOSE,
                      alpha,
                      bsrA,
                      descrA,
                      y_v,
                      tmp_v );
    printf( "   OUTPUT DATA FOR mkl_sparse_d_trsv \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f\n", tmp_v[i] );
    }
    printf( "-------------------------------------------------------\n" );

//  Task 3: Obtain matrix-matrix multiply A' *x_m --> y_m
//          A - zero-based indexing,
//          x_m - column major ordering
    printf( "                                  \n" );
    printf( "   Task 3:                        \n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mm \n" );
    printf( "   WITH GENERAL SPARSE MATRIX     \n" );
    printf( "   COLUMN MAJOR ORDERING for RHS  \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_TRANSPOSE     \n" );
    printf( "   Input vectors                  \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f, %7.1f\n", x_m[i], x_m[M*LB+i] );
    }

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    mkl_sparse_d_mm ( SPARSE_OPERATION_TRANSPOSE,
                      alpha,
                      bsrA,
                      descrA,
                      SPARSE_LAYOUT_COLUMN_MAJOR,
                      x_m,
                      NRHS,   // number of right-hand-side
                      M*LB,   // ldx
                      beta,
                      y_m,
                      M*LB ); // ldy

    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f, %7.1f\n", y_m[i], y_m[M*LB+i] );
    }
    printf( "-------------------------------------------------------\n" );

//  Task 4: Obtain matrix-matrix multiply A*x_m --> y_m
//          A - zero-based indexing,
//          x_m - row major ordering
    printf( "                                  \n" );
    printf( "   Task 4:                        \n" );
    printf( "   INPUT DATA FOR mkl_sparse_d_mm \n" );
    printf( "   WITH GENERAL SPARSE MATRIX     \n" );
    printf( "   ROW MAJOR ORDERING for RHS     \n" );
    printf( "   ALPHA = %4.1f  BETA = %4.1f    \n", alpha, beta );
    printf( "   SPARSE_OPERATION_TRANSPOSE     \n" );
    printf( "   Input vectors                  \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f, %7.1f\n", x_m[2*i], x_m[2*i+1] );
    }

    // Create matrix descriptor
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    mkl_sparse_d_mm ( SPARSE_OPERATION_TRANSPOSE,
                      alpha,
                      bsrA,
                      descrA,
                      SPARSE_LAYOUT_ROW_MAJOR,
                      x_m,
                      NRHS,   // number of right hand sides
                      NRHS,   // ldx
                      beta,
                      y_m,
                      NRHS ); // ldy

    printf( "   OUTPUT DATA FOR mkl_sparse_d_mv \n" );
    for ( i = 0; i < M*LB; i++ )
    {
        printf( "%7.1f, %7.1f\n", y_m[2*i], y_m[2*i+1] );
    }

    // Release matrix handle and deallocate matrix
    mkl_sparse_destroy ( bsrA );

    printf( "-------------------------------------------------------\n" );
    return 0;
}