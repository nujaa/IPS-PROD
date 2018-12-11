#include <iostream>
#include <armadillo>
#include <fstream>
#include "poly.h"
#include "basis.h"


int main(){	
	//calcul des coefficients
	double inf = -0.1, sup = 0.1, nbX = 10;
	int n = 5, i;
	// arma::vec X = arma::linspace(inf, sup, nbX);
	arma::vec X = {-3.1, -2.3, -1.0, -0.3, 0.1, 4.3, 9.2, 13.7};
/*	// Poly* P = new Poly((uint)5, (uint)5, (uint)5, X);
	std::ofstream PsiFile("toPlot/Psis201.dat");
	// std::ofstream Cube("toPlot/Cube.dat");
	Poly* P = new Poly();
	P->calcHermite(1, X);
	P->calcLaguerre(1, 4, X);
	P->hermite(0);
	P->laguerre(0, 3);
	//calcul des solutions
	// P->getHermiteValues().print();
	// P->getLaguerreValues().print();
	PsiFile << P->getHermiteValues() << std::endl;*/
	// Cube << P->getLaguerreValues() << std::endl;
	Basis b = Basis(1.935801664793151, 2.829683956491218, 14, 1.3);
	arma::vec r = {-10.1, -8.4, -1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
	//for(i = 0; i < 15; i ++)
	// b.zPart(r, 15).print();
	printf("%lf", (double)1.0 / (double)(sqrt(2.829683956491218 * pow(2, 15) * sqrt(pi) * factorial(15))));
	arma::exp(-r % r / (2 * (2.829683956491218 * 2.829683956491218))).print();
}