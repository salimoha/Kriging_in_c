/*Notes:
 *
 *  liklihood function
 * Calculating the likelihood function and objective function for boxmin.cpp
 * output:
 *         obj:  likelihood value
 *
 *  This function needs to be improved (This is called inverse problem) via
 * regularization algorithms such as:
 *
 * 2 gaussian processes -> booker method, Gtech method, a group in Isreal
 *
 * Or other method such as modified cholesky ( I like this method more. )
 *
 * Advantage of this code:
 *                      QR decomposition does not change the phase of the answer
 *                      and it helps to improve the ill-conditioning.
 *                      Modified Cholesky decompostion also helps to make the
 *                      problem well conditioned
 *
 *
 * Author: Shahrouz Alimohammadi
 * Modified: August 2014
 * Email: shahrouz.alm@gmail.com
*/
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cassert>
#include "cholesky.cpp"
#include "evalQRsolver.cpp"
#include "corr.cpp"
#include "regr.cpp"


double objfunc(MatrixD S, VecD Y, VecD& theta, VecD& sigma2,
		void (*corr)(matrix<double>& S, matrix<double>& R,vector<double> theta )
		, void(*regr) (MatrixD , MatrixD& ) ,MatrixD& C,VecD& gamma, VecD& Beta  ){

	// finding the regression for the points
	MatrixD F(m,n+1,0);
	(*regr)(S, F);

	for ( int i =0; i < n; i++)
		if (theta(i) < 0.0000000001){
			theta(i) =0.00000001;
			std::cout << " theta got close to 0.....!!!" << std::endl;
		}

	// constructing the correlation matrix
	// and cholesky factorization
	MatrixD R(m,m,0);
	(*corr)(S,R,theta); //calculating R
	//    std::cout << "\n \n ==========\n R ==========\n =  " << R <<std::endl;

	///**************///// needs to be improved///////////////////////////////////
    /**** here is the cholesky problem and the ill conditioning problem that we have***/

// modified cholesky with epsilon
    double epsilon = 0.00001;
	for ( int i=0 ; i < R.size1(); i++)
        R(i,i)=R(i,i) + epsilon;

// using cholesky decompostion via boost libarary
	size_t res = cholesky_decompose(R, C);
	//std::cout << "------------------"<< std::endl;
	//std::cout << res << std::endl;
	//std::cout << "R:\n" << R << std::endl;
	//std::cout << "C:\n" << C << std::endl;
	//std::cout << " \n C test:\n" << R - prod(C,trans(C) ) << std::endl;
	assert(!res);

	///////////////////////////////////////////////////////////////////////////////////
	for ( int i=0 ; i < C.size1(); i++)
		if ( C(i,i) < 0.000001 ){
			std::cout << "C has a problem in "<< i << "element" << std::endl;
			return 0;
		}
	VecD Yt(m,0);
	Y.resize(m);
	evalQRFactor( C, Y, Yt);

	VecD V =  Y - prod(C,Yt);

	std :: cout <<" \n **********test forthe Yt for SparseQR scheme:********\n V= \n "<< V << std::endl;
	double errY=0;
	for ( int i =0; i < V.size(); i++)
		errY+=V(i)*V(i);
	errY = sqrt(errY);

	std :: cout <<" \n **********test forthe Yt for SparseQR scheme:******** errY=  "<< errY << std::endl;

	// QG for F~
	MatrixD Ft(m , n + 1,0);
	F.resize(m,n+1);
	evalQRFactorM(C, F, Ft );
	//double errF = norm_M(F - prod(C,Ft));

	MatrixD M = F - prod(C,Ft);
	double d=0;
	for (int i=0; i < M.size1() ; i++)
		for ( int j =0; j< M.size2() ; j++)
			d += M(i,j)*M(i,j);
	d = sqrt(d);



//	std :: cout <<" \n ***************test for the Ft:*********** errF= "<< d << std::endl;
//	std::cout << "F=" << F << std::endl << "Ft=" << Ft << std::endl;
//	std::cout << "Y=" << Y << std::endl << "Yt=" << Yt << std::endl;

	std :: cout <<"==========\n" << std::endl;

//	std :: cout <<"==========> C: \n "<< C << std::endl;

	// beta = G\Q'Y~
	// backsubstitutoion
	// inputMat Ft, inputRhs Yt, sol Beta
	// SOLVE SYSTEM WITH QR DECOMPOSITION
	//VecD Beta(n+1, 0); // went to main
	evalQRFactor( Ft,Yt, Beta);
//	std:: cout << "Beta= " << Beta << std::endl;

	// sigma2

	//VecD sigma2(n,0);
	VecD FtB= prod(Ft, Beta);
//	std::cout << "FtB****************** \n " << FtB <<std::endl;
  //  std::cout << "size of FtB****************** \n " << FtB.size() <<std::endl;
    vector<double> rho(m);
    for ( int i=0; i < m; i++){
    //	sigma2(i) = 1.0/m*pow( Yt(i)-FtB(i),2 );
//        sigma2(i) = pow( Yt(i)-FtB(i),2 );
        rho(i) = pow( Yt(i)-FtB(i),2 );
	}
// std::cout << "rho****************** \n " << rho <<std::endl;

 // likelihood
 double obj=0;
 double sumsigma2=0;
 for (int i=0; i<m ; i++) {
     sumsigma2 += rho(i);
 }
 sumsigma2=sumsigma2 / m;

//	std:: cout << "variance=" << sigma2 << std::endl;
	// det R
	double lRl=1;

	for( int i=0; i < C.size1(); i++)
        lRl *= pow(C(i,i), 2.0/m);

	std::cout << "det(R) = " << lRl << std::endl;
    obj = sumsigma2*lRl;
//	std::cout << "\n \n +++++++++sigms2+++++++++= \n"<< sumsigma2 << std::endl;
	//   std::cout << "\n \n +++++++++det R+++++++++= \n"<< lRl << std::endl;
	std::cout << "likelihood= "<< obj << std::endl;
	// calculating gamma
	VecD RHS = Yt - FtB;
	MatrixD Ct = trans(C);
	evalQRFactor( Ct,RHS ,gamma );

    // obj : likelihood value and objective function for boxmin.cpp
	return obj;

}



