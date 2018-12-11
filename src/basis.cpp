#include <stdio.h>
#include "basis.h"
#include "poly.h"

#define NZ (N + 2) * pow(Q, 2./3) + 1./2

Basis::Basis(double Br, double Bz, int N, double Q){
	br = Br, bz = Bz;
	mMax = 0;
	setmMax(N, Q);
	nMax = arma::vec(mMax);
	setNMax();
	n_zMax = arma::mat(mMax, nMax(0), arma::fill::zeros);
	setN_zMax(N, Q);
}

int Basis::n_ZMax(uint i, int N, double Q){
	return (NZ - (i * Q));
}

void Basis::setmMax(int N, double Q){
	while(n_ZMax(mMax+1, N, Q) >= 1.)
		mMax ++; //find mMax
}

void Basis::setNMax(){
	int m;
	for(m = 0; m < mMax; m ++){
		nMax(m) = ((mMax - m - 1) / 2) + 1;
	}
}

void Basis::setN_zMax(int N, double Q){	
	int i, j;
	for(i = 0; i < mMax; i ++){
		n_zMax(i, 0) = n_ZMax(i+1, N, Q);
		for(j = 1; j < nMax(0); j ++){ // on remplit les valeurs dans les autres colonnes qui ont la meme.
			if(i - 2 * j < 0)
				break;
			n_zMax(i - 2 * j, j) = n_zMax(i, 0);
		}
	}
}


arma::vec Basis::rPart(arma::vec r, int m, int n){
	Poly poly;
    poly.calcLaguerre(m, n, r % r / (br * br));
    double f = sqrt(factorial(n) / (pi * factorial(n + m))) / br;
    arma::vec e = arma::exp(- r % r / (2* (br * br)));
    arma::vec p = arma::pow(r / br, m);
    arma::vec L = poly.laguerre(m, n);
    return f * e % p % L;
} 
arma::vec Basis::zPart(arma::vec z, int nZ){
	Poly poly;
	printf("yo\n");
    poly.calcHermite(nZ + 1, z / bz);
    printf("yo\n");
    poly.hermite(nZ).print();
    long double f    = 1.0 / sqrt(bz * pow(2, nZ) * sqrt(pi) * factorial(nZ));
	printf("yo\n");
    arma::vec e = arma::exp(-z % z / (2 * (bz * bz)));
	e.print();
	fflush(stdout);

    arma::vec h = poly.hermite(nZ);
    printf("yo\n");
    return f * e % h;
}
/*int main(){
	Basis b = Basis(1.935801664793151, 2.829683956491218, 14, 1.3);
	arma::vec r = {-10.1, -8.4, -1.0, 0.0, 0.1, 4.3, 9.2, 13.7};
	b.zPart(r, 15).print();
}*/