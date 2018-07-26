// regression function
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

void regpoly1(MatrixD S, MatrixD& f/*, MatrixD& df*/){

    int row = S.size1();
    int col = S.size2();

    for ( int i=0; i < row; i++){
        f(i,0) =1;
        for ( int j=0; j < col; j++)
            f(i,j+1) = S(i,j);
    }
    return;
}

void regpoly0(MatrixD S, MatrixD& f/*, MatrixD& df*/){

    int row = S.size1();
    int col = S.size2();

    for ( int i=0; i < row; i++)
        f(i) =1;
    return;
}


