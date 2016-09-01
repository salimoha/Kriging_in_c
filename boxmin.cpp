/* finding the optimum value for the likelihood fuction (objfun.cpp) and the minimizer theta_opt via: boxmin.cpp
 * boxmin.cpp is the equvalent implimentation of boxmin function in
 * MATLAB dacefit package. Here GPS method is used to find the theta_optimum; however,
 * in the MATLAB code Hooke Jeeves method is implimented.
 * Future work: One can improve the searching point from GPS to the other methods,
 * for example, Hooke and jeeves.
 *
 * Advantage of Hooke and Jeeves method: it is a self starting method.
 * Disadvantage of GPS: needs a initial mesh size for THETA (different from the grid mesh size)
 * Disscussion: As you can see the value of the local minimum is not criticaly important
 *, if you choose a big (not very big though!) you can get a decent value for theta_opt.
 * In this code, delta is initialized in dacefit.cpp function.
 *
 * IMPORTNAT NOTE: finding a right value for theta_opt is extreremly important.
 * If one uses Cholesky method or LU or Jones et. al. 1999 method; will end up
 * with different theta value and different interpolation.
 * Cholesky with Suitsparse library works the best in ill-conditioned cases.
 *
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
#include "cholesky.cpp"


using namespace boost::numeric::ublas;

// finding the theta optimum
// by using a derivative free optimization scheme
// via General Pattern Serach Algorithm (GPS)
void boxmin(MatrixD S, VecD Y, VecD& theta, VecD& sigma2,MatrixD&C, VecD& gamma, VecD& Beta, double& delta){
    int N = 2*n; //GPS strategy
  //  double J_opt;
    //double delta=10;
    std::cout << "enter the initial mesh size delta=  "<<std::endl;
//    std::cin >> delta;

    //delta = theta(0)*4; /// needs to be imporve******coefficient has a problem if I put 5 it will have problem at 5
    double cur_best_J;
    double J_prev;
   // double J;
    VecD theta_opt;
    double residual = 0.00001;
  //  int id;

// to check the value of the likelihood fucntion
 // the value is saved at J.txt file
  // one can remove this proccess ( all myfile*)
	  std::ofstream myfile;
  myfile.open ("J.txt");

    //VecD f(n+1,0);
    VecD t(n,10);
    VecD theta_old(n,0);
    // follow 2n neighbors
    MatrixD T(n,N,0); // neighbors

//	cur_best_J = fun_eval(theta); // for checking the algorith with other functions
    cur_best_J = objfunc(S, Y,theta,sigma2, corr_gauss, regpoly1, C, gamma, Beta);

    int iter=0;
    int max_iter = 1000;
    VecD f(N+1,0);
    while (delta > residual && iter < max_iter ){

    iter++;
    // N = 2*n
        //std::cout << " best theta is: " << std::endl;
            for ( int i=0; i < n; i++){
            for ( int j = 0; j < N; j++){
            //	if ( i==j ||  i + N/2 == j)
                //assert(theta(i) >= delta );
                    if ( i==j) 	T(i,j) = theta(i) + delta;
                    else if( i + N/2 == j )	 T(i,j) = theta(i) - delta;
                else
                    T(i,j) = theta(i);
              //  std::cout <<"   "<< T(i,j);
            }
      //  std::cout << std::endl;
    }
            J_prev = objfunc(S , Y, theta, sigma2, corr_gauss, regpoly1, C, gamma, Beta );
//	J_prev = fun_eval(theta);

    cur_best_J = J_prev;
    theta_old = theta;
    for ( int j = 0 ; j < N; j++){
        for ( int i =0; i < n; i++)
                t(i) = T(i,j); // t is the current point to evaluate
        f(j) = objfunc(S,Y, t, sigma2, corr_gauss, regpoly1, C, gamma, Beta);
//		f(j) = fun_eval(t); // for checking the algorith with other functions
                J_prev = cur_best_J;
                if ( cur_best_J - f(j) > 0 ){
                     cur_best_J = f(j);
                    theta = t;

                }

    }
    if (cur_best_J == J_prev  ){
        delta = delta / 2.0;
    std::cout << " 2: mesh refinment " << std::endl;
    }
    else 
	//	delta = 2*delta;
            std::cout << " 1: polling+++++++++ " << std::endl;
//	std::cout << " delta = " << delta << std:: endl;
//    std::cout << " best point theta= " << theta << std:: endl;
  //  std::cout << " current best J = " << cur_best_J << std::endl;

    myfile << cur_best_J<<"\t" << theta << std::endl; 

}

 myfile.close();


return;
}
