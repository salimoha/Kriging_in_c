/* Calculating QR decompostion
 * using SuiteSparse library
 * For:
 *     Ax=b    A matrix but x and b are vector
 * and
 *     AX=B   A, X, B all are matrix
 *
 *
 * here we have
 *  matrix and vector in boost as an input and out put
 *  and calculation in cholmod sparse matrices
 *
 *  calculation time does not get affected since objective function is
 * very expensive but one can improve this funtion.
 *
 * IMPROTANT FUNCTION
 *
 * Author: Shahrouz Alimohammadi
 * Modified: August 2014
 * Email: shahrouz.alm@gmail.com
 */
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "SuiteSparseQR.hpp"

// =========================
// EVAL QR FACTOR FOR VECTOR
// =========================
using namespace boost::numeric::ublas;
void evalQRFactor(matrix<double> &inputMat, vector<double> &inputRHS, vector<double>& sol){
    // DECLARE

    // CHOLMOD
    cholmod_common Common, *cc ;
    cholmod_dense *Residual;
    cholmod_dense* X;
    double rnorm, one [2] = {1,0}, minusone [2] = {-1,0} ;
    //int mtype ;
    //double* solCoeff;

    // start CHOLMOD
    cc = &Common ;
    cholmod_l_start(cc);

    // GET THE SIZE OF THE MATRIX
    long nrow = inputMat.size1();
    long ncol = inputMat.size2();

    // CREATE A NULL CHOLMOD MATRIX and RHS
    cholmod_dense* A = (cholmod_dense*)cholmod_l_eye(nrow,ncol,CHOLMOD_REAL,cc);
    cholmod_dense* B = (cholmod_dense*)cholmod_l_ones(nrow,1,CHOLMOD_REAL,cc);

    // CONVERT BOOST MATRIX AND RHS INTO CHOLMOD MATRIX
    for(int loopA=0;loopA<nrow;loopA++){
        ((double*)B->x)[loopA] = inputRHS(loopA);
        for(int loopB=0;loopB<ncol;loopB++){
            ((double*)A->x)[loopA + loopB*A->d] = inputMat(loopA,loopB);
        }
    }
/*
    std::cout <<"-----------"<<std::endl;
    for(int loopA=0;loopA<nrow;loopA++){
        for(int loopB=0;loopB<ncol;loopB++){
            std::cout << ((double*)A->x)[loopA + loopB*A->d]<<" ";
        }
    }
*/	

    // PRINT DENSE MATRIX
//    cholmod_l_print_dense(A,"A",cc);

    // PRINT DENSE VECTOR
  //  cholmod_l_print_dense(B,"B",cc);

    // CONVERT CHOLMOD MATRIX TO SPARSE MATRIX
    cholmod_sparse* spA = cholmod_l_dense_to_sparse(A,1,cc);

    // PRINT SPARSE MATRIX
    //cholmod_l_print_sparse(spA,"A",cc);

    // X = A\B
    X = SuiteSparseQR<double>(spA,B,cc);

    // PRINT DENSE MATRIX
   // cholmod_l_print_dense(X,"X",cc);

    // rnorm = norm (B-A*X)
    Residual = cholmod_copy_dense (B, cc) ;
    cholmod_sdmult (spA, 0, minusone, one, X, Residual, cc) ;
    rnorm = cholmod_norm_dense (Residual, 2, cc) ;
   // printf ("2-norm of residual: %8.1e\n", rnorm) ;
   // printf ("rank %ld\n", cc->SPQR_istat [4]) ;

    // COPY THE SOLUTION OVER THE RHS
    sol.resize(ncol);
    if(X != NULL){
        for(int loopA=0;loopA<ncol;loopA++){
            sol(loopA) = ((double*)X->x)[loopA];
        }
    }else{
        for(int loopA=0;loopA<ncol;loopA++){
            sol(loopA) = 0.0;
        }
    }

    // free everything and finish CHOLMOD
    cholmod_free_sparse(&spA, cc) ;
    cholmod_free_dense(&Residual, cc) ;
    cholmod_free_dense(&A, cc) ;
    cholmod_free_dense(&B, cc) ;
    cholmod_finish(cc) ;

}

// =========================
// EVAL QR FACTOR FOR MATRIX
// =========================

void evalQRFactorM(matrix<double> &A, matrix<double> &B, matrix<double>& X){
for ( int k = 0 ; k < X.size2(); k++){

    VecD tmpSol(X.size1(),0);
    VecD tmpRHS(B.size1(),0);
    for ( int i=0; i < X.size1(); i++){
        //tmpSol(i)= X(i,k);
        tmpRHS(i)=B(i,k);
    }

   evalQRFactor(A, tmpRHS, tmpSol);

   for ( int i=0; i < X.size1(); i++)
        X(i,k)=tmpSol(i);
}
return;
}

