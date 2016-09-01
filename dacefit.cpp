/* dacefit.cpp
 * DACEFIT Constrained non-linear least-squares fit of a given correlation
 *  model to the provided data set and regression model
 *
 *  based on:
 *              Sacks et al 1989
 *
 * Author: Shahrouz Alimohammadi
 * Modified: August 2014
 * Email: shahrouz.alm@gmail.com
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <fstream>
#include <cassert>
#include <math.h>
//#include "SuiteSparseQR.hpp"
//#include "cholesky.cpp"
//#include "normalize.cpp"
//#include "objfunc.cpp"
//#include "boxmin.cpp"


using namespace boost::numeric::ublas;



void dacefit(matrix<double>& S, VecD& Y, matrix<double>& Ssc, vector<double>& Ysc ,vector<double>& theta,
        void (*corr)(matrix<double>& , matrix<double>& ,vector<double> )
        ,void(*regr) (MatrixD , MatrixD& ),MatrixD C, VecD& gamma, VecD& Beta   ){

    // DACEFIT***********************
    /***********************initial points*****************************************/
    // calculating mean and variance of data
        // and normalize them which we have mean=0 and variance=1
      normalize( S,  Y, Ssc, Ysc);

  //  std:: cout << "****** normalize \n S , Y "<< S << Y << "\n =========\nSsc , Ysc   " << Ssc << Ysc << std::endl;

    double J;
    VecD sigma2(n,0);
    J = objfunc(S, Y,theta,sigma2, (*corr), (*regr), C, gamma, Beta);
    //J = objfunc(S, Y,p,theta,sigma2, corr_gauss, regpoly1);
  //  std:: cout << "******likeihood = "<< J << std::endl;
    /*******************************finding the theta_star****************************/
  double delta=20;

  // finding optimization parameters to generate the
  // kriging or dacefit model
 boxmin( S, Y, theta, sigma2, C, gamma, Beta, delta);
        std:: ofstream betaOut, gammaOut, sigma2Out;
    betaOut.open("beta.txt", std::ios::app);
    gammaOut.open("gamma.txt", std::ios::app );
    sigma2Out.open("sigma2.txt", std::ios::app);
    betaOut << Beta<< std::endl;
    sigma2Out << sigma2<< std::endl;
for ( int i =0 ; i < gamma.size(); i++)
 gammaOut << gamma(i)  << std::endl;
 gammaOut << std::endl;



    std:: cout << "******\n sigma= "<< sigma2 << std::endl;
    std:: cout << "******\n Beta= "<< Beta << std::endl;
    std:: cout << "******\n sigma= "<< sigma2 << std::endl;
//std::cout << "********\n theta =   " << theta << std::endl;



    /******************************* S = gridsamp(range, q) **********************/


    return;
}
