/* Kriging predictor function
 * one of the 2 important function
 * in kriging code.
 *
 * y = corrterm + regrterm;
 *
 * Author: Shahrouz Alimohammadi
 * Modified: Dec 2013
 * Email: shahrouz.alm@gmail.com
 */

#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "SuiteSparseQR.hpp"


//FIXME
// correlation function for a vector: r
void my_corr_x( VecD x, vector<double> theta, MatrixD S , vector<double> p, vector<double>& r){
	int m = r.size();
	int n = theta.size();
	double temp;
	for ( int j=0 ; j < m ; j++) {
		temp=0;
		for ( int l=0; l < n ; l++)
			temp += theta(l) * pow( S(j,l)-x(l), p(l) );
		r(j)=exp(-temp );
	}
	return;
}

// regresion function
void regpoly1_x(VecD x, VecD& f/*, MatrixD& df*/){
	
    int row = x.size();
     f(0) =1;
    for ( int i=0; i < row; i++){
        f(i+1) = x(i);
    }
    return;
}


/*****************predictor*******************/
double predictor(VecD x, MatrixD S, MatrixD Ssc, VecD Y, VecD Ysc, VecD theta,
			   MatrixD C, VecD gamma, VecD Beta){

	int N = S.size2();
	int M = S.size1();
    //N=n; M = m;


    for ( int i =0; i < x.size(); i++)
    x(i) = ( x(i) - Ssc(i,0) ) / (Ssc(i,1) );


	double y=0;
	// MatrixD x(1,m,0);
	VecD r(M,0);
  	VecD rt(M,0);
	VecD p(N,2); // gaussian correlation
    vector<double> f( N+1,0 );
	my_corr_x(x, theta, S ,  p,  r);
	regpoly1_x(x, f);


   // std::cout << "f =" << f << std::endl;
  //  std::cout << "r=" << r <<std::endl;

	// rt = C \ r
	evalQRFactor(C , r , rt);
 //   std::cout << "**********rt=" << rt <<std::endl;
	
	double corrterm=0;
	double regrterm=0;
    for( int i=0; i < N+1; i++)
        regrterm += f(i)*Beta(i);


	  corrterm = inner_prod(r, gamma);
// std::cout << " Z1 =  " << corrterm << std::endl;
/*
	for (int j=0; j < M; j++) {
        corrterm += r(j)*gamma(j);
	}
*/
	y = corrterm + regrterm; 

   y = Ysc(1)*y+Ysc(0);

     std::cout << "**********rgterm1=" <<  regrterm <<std::endl;
// std::cout << " r =  " << r << std::endl;
// std::cout << " gamma =  " << gamma << std::endl;
// std::cout << " Z1 =  " << corrterm << std::endl;

return y;
}








