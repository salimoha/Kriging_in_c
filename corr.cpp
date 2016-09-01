/* Gaussian Correlation function in boost library
 * One can use other correlation functions
 * like: cubic, exponenstial, spline, mixed,
 *      the one is introduced in booker paper
 *
 *
 * Author: Shahrouz Alimohammadi
 * Modified: August 2014
 * Email: shahrouz.alm@gmail.com
 */
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric::ublas;

// guassian correlation p=2
void corr_gauss(matrix<double>& S, matrix<double>& R,
        vector<double> theta ){
double eps = 1E-20;
    VecD p(n,2); // 0 < pl <= 2
    long double s=0; // temporary variable
    for ( int i = 0; i < m; i++)
        for ( int j = 0; j < m; j++){
            s=0;
            for ( int l=0; l < n; l++){
                s = s + theta(l)*pow( S(i,l)-S(j,l), p(l) );
            }
            R(i,j) = exp(-s);
    if ( R(i,j) < eps)
        R(i,j) = 0.0;
        }
    return;
}
