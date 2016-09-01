/* Normalizing the initial data
 * such as matlab code
 * this process is very important.
 * making:
 *    mean = 0
 *    variance = 1
 *
 * one can improve this function
 *
 * Author: Shahrouz Alimohammadi
 * Modified: Dec 2013
 * Email: shahrouz.alm@gmail.com
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
//#include "normalize.h"

using namespace boost::numeric::ublas;


void normalize(matrix<double>& S, vector<double>& Y, matrix<double>& Ssc, vector<double>& Ysc ){


int n = S.size2();
int m = Y.size();
 
VecD mS(n,0);
//mean( S, mS);
double sum =0;
    for ( int j=0; j < n ; j++){
    sum=0;
        for ( int i=0; i < m ; i++)
            sum+=S(i,j);
    mS(j)=sum/S.size1();
    }
std::cout<<"mean S: mS= "<< mS << std::endl;
VecD var(n,0);
//variance(S, var);
 sum=0;
        for ( int j=0; j < n ; j++){
         sum=0;
            for ( int i =0; i < m ; i++){
            sum+=pow( (S(i,j)-mS(j)) , 2);
		}
    var(j)=sqrt(sum*1.0/(m-1));
        }
        std::cout<<"variance S: varS= "<< var<< std::endl;


// writing mS and var in Ssc
	for ( int i=0; i < n; i++){
	Ssc(i,0) = mS(i);
	Ssc(i,1) = var(i);
	}
	


//std::cout<<"1. S= "<< S << std::endl;



//normalize(S,Y);
for( int j = 0; j < n; j++)
    for ( int i = 0; i < m; i++)
        S(i,j) = (S(i,j) - mS(j)) / var(j);
//std::cout<<"2. S= "<< S << std::endl;
//mean(Y,mY)
double mY=0;
for( int i=0; i < m; i++)
mY+=Y(i);
mY=mY/Y.size();
//variance(Y,varY);
double varY=0;
for ( int i=0; i < m; i++)
   varY += pow( Y(i) - mY,2);
varY=sqrt(varY*1.0/(m-1));

if(Y.size()!=m) std:: cout << "errooor: Check the size of Y. It is not equal to m"<<std::endl;
//normalize(Y)

for ( int i=0; i < Y.size(); i++)
Y(i) = ( Y(i) - mY ) / varY;

// in order to check mean = 0 and variance = 1
sum=0;//S
for ( int j=0; j < n ; j++){
    sum=0;
        for ( int i=0; i < S.size1() ; i++)
            sum+=S(i,j);
    mS(j)=sum/S.size1();
    }
// writing the value of mean(Y) and variance(Y)

Ysc(0) = mY;
Ysc(1) = varY;


// checking
mY=0;//Y
for( int i=0; i < m; i++)
mY+=Y(i);
mY=mY/Y.size();

std::cout<<"mean S: mS="<<mS<< std::endl;
std::cout<<"mean Y: mY="<<mY<< std::endl;


//also to check variance = 1
 sum=0;//S
        for ( int j=0; j < n ; j++){
         sum=0;
            for ( int i =0; i < m ; i++)
            sum+=pow(S(i,j)-mS(j),2);
    var(j)=sqrt(sum/(m-1) );
        }
varY=0; //Y
for ( int i=0; i < m; i++)
   varY += pow( Y(i) - mY,2);
varY=sqrt(varY/(m-1) );

std::cout<<"variance S: var="<<var<< std::endl;
std::cout<<"variance Y: varY="<<varY<< std::endl;
std::cout<<"sample pionts S="<<S<< std::endl;
return;
}

