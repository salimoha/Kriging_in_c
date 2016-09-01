/* this is a script to test kriging code
 *
 *
                 Shahrouz Alimohammadi
                 March 2014
 */



#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
using namespace boost::numeric::ublas;


int n=1, m=17;
typedef vector<double> VecD;
typedef matrix<double> MatrixD;
#include "normalize.cpp"
#include "objfunc.cpp"
#include "boxmin.cpp"
#include "predictor.cpp"
#include "dacefit.cpp"


int main(){
/*
    std::ofstream m_file;
    m_file.open("m.txt");
    m_file << m << std::endl;
*/
std::ifstream m_file;
    m_file.open("m.txt");
    m_file >> m ;

	vector<double> beta(n+1);
	vector<double> gamma(m);
	matrix<double> S(m,n);
	matrix<double> C(m,m);
	vector<double> Y(m);
	vector<double> sigma2(n);
	vector<double> theta(n,2);

	double a;
//	std:: cout << " enter theta: \n ";
//    std::cin >> a;
 //   theta(0) =a;


/*	std::ifstream b_file("beta.txt"), g_file("gamma.txt"),
 S_file("/home/shahrouz/Desktop/test_C/SMF_STAND_ALONE_CPP/gitlib_test/kriging_jones/testMydacefit/S.txt"),
 C_file("C.txt"), Y_file("/home/shahrouz/Desktop/test_C/SMF_STAND_ALONE_CPP/gitlib_test/kriging_jones/testMydacefit/Y.txt");
*/ std::ifstream b_file("beta.txt"), g_file("gamma.txt"), S_file("S.txt"), C_file("C.txt"), Y_file("Y.txt");

	for(int i=0; i< S.size1(); i++)
		for(int j=0; j< S.size2(); j++)
			S_file >> S(i,j);

	for ( int i=0; i < Y.size(); i++)
		Y_file >> Y(i);


//	for(int i=0; i< C.size1(); i++)
//		for(int j=0; j< C.size2(); j++)
//			C_file >> C(i,j);

//	for(int i=0; i<beta.size(); i++)
//		b_file >> beta(i);

//	for(int i=0; i<gamma.size(); i++)
//		g_file >> gamma(i);




	std::cout << "S ="<<S << std::endl;
	std::cout << "Y ="<<Y << std::endl;
//	std::cout <<  "beta= " << beta << std::endl;
//	std::cout << "gamma= " << gamma << std::endl;
//	std::cout << "C= " << C << std::endl;


    matrix<double> Ssc(n,2);
    vector<double>Ysc(2);
    std::ofstream Jout, tout;
    Jout.open("j.txt");
    tout.open("theta.txt");
    double J;




dacefit( S, Y ,Ssc, Ysc, theta, corr_gauss , regpoly1, C, gamma, beta);

/*
 tout << theta(0) << std::endl;
                Jout << J <<std::endl;
	std:: ofstream betaOut, gammaOut, COut, sigma2Out;
    betaOut.open("beta.txt");
    gammaOut.open("gamma.txt");
    COut.open("C.txt");
    sigma2Out.open("sigma2.txt");


    betaOut << beta<< std::endl;
    gammaOut << gamma << std::endl;
    COut << C << std::endl;
    sigma2Out << sigma2 << std::endl;
*/
    vector<double> x(n,6);

    double yy = predictor(x,S, Ssc, Y, Ysc, theta,C, gamma, beta);

std::cout << " yy = " << yy <<std::endl;



std::ofstream Fout, Xout;
//Fout.open("compare_to_matlab/yy.txt");
//Xout.open("compare_to_matlab/xx.txt");
VecD xx(n,0);
Fout.open("yy.txt");
Xout.open("xx.txt");

while ( xx(0) < 1){
    yy = predictor(xx,S, Ssc, Y, Ysc, theta,C, gamma, beta);
    Xout << xx(0) <<std::endl;
    Fout << yy << std::endl;
xx(0) += 0.005;
}


std::cout << " theta = " << theta <<std::endl;
std::cout << " yy = " << yy <<std::endl;

    return 0;
}




