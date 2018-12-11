#include <stdio.h>
#include "poly.h"

long factorial(unsigned long f)
{
	double res = 1;
	int i;
	for(i=1; (uint)i < f;i++){
		res = res * i;
	}
	return res;
}

Poly::Poly(){}

Poly::Poly(uint nHermite, uint nLaguerre, uint mLaguerre, arma::vec Z){
	this->calcHermite(nHermite, Z);
	this->calcLaguerre(nLaguerre, mLaguerre, Z);
}

//getters
uint Poly::getNH(){
	return nH;
} 

uint Poly::getNL() {
	return nL;
}

uint Poly::getML() {
	return mL;
}

arma::mat Poly::getHermiteValues() {
	return hermiteValues;
}

arma::cube Poly::getLaguerreValues() {
	return laguerreValues;
}
arma::vec Poly::getZValues(){
	return zValues;
}

//Setters
void Poly::setNH(uint n) {
	nH = n;
}

void Poly::setNL(uint n) {
	nL = n;
}

void Poly::setM(uint m) {
	mL = m;
}

//Hermite
void Poly::setHermiteValues(arma::mat values) {
	hermiteValues = values;
}
void Poly::setLaguerreValues(arma::cube values) {
	laguerreValues = values;
}

//Laguerre
void Poly::setZValues(arma::vec values) {
	zValues = values;
}

void Poly::calcHermite(uint n, arma::vec X) {
	arma::vec ones = arma::vec(X.n_elem, arma::fill::ones);
	arma::vec zVal = X * sqrt(omega0), Z = zVal % zVal * omega0 / 2;
	int nbX = zVal.n_elem, i;
	if(n<1 || nbX < 1){
		return;
	}
	resHermite = arma::mat(nbX, n, arma::fill::ones);
	if(n > 1)
		resHermite.col(1) = 2 * zVal;
	for(int i = 2; (unsigned)i < n; i ++){
		resHermite.col(i) = 2 * zVal % resHermite.col(i-1) - 2 * (i-1) * resHermite.col(i-2);
	}

	//calcul des hermiteValuesutions
	
	this->hermiteValues = arma::mat(nbX, n, arma::fill::ones);
	for(i = 0; (unsigned)i < n; i++){
		hermiteValues.col(i) =  
						(1 / sqrt(pow(2, i) * factorial(i))
						* pow(omega0 / pi, 1/4))
						* arma::exp(-Z)
						% resHermite.col(i);
	}
}

arma::vec Poly::hermite(uint n) {
	if(!hermiteValues.is_empty() && hermiteValues.n_rows > n)
		return resHermite.col(n);
	else
		return NULL;
}

void Poly::calcLaguerre(uint m, uint n, arma::vec zVal) {
	int nbX = zVal.n_elem; 
	double i, j;
	resLaguerre = arma::cube(nbX, m, n, arma::fill::ones);

	arma::vec courant = arma::vec(nbX);
	resLaguerre.slice(0) = arma::mat(nbX, m, arma::fill::ones); //set slice 0
	for(i=0;i<m;i++){											//set slice 1
		resLaguerre.slice(1).col(i) = arma::vec(nbX, arma::fill::ones) * (1 + i) - zVal;
	}
	

	// récurrence

	for(i=2;i<n;i++){
		for(j=0;j<m;j++){
			resLaguerre.slice(i).col(j) =
			(2 * arma::vec(nbX, arma::fill::ones) + 1 / i * ((j - 1) * arma::vec(nbX, arma::fill::ones) - zVal)) % resLaguerre.slice(i-1).col(j)
			- (1 + 1/i * (j - 1)) * resLaguerre.slice(i-2).col(j);
		}
	}
	laguerreValues = arma::cube(nbX, m, n, arma::fill::ones);
}

arma::vec Poly::laguerre(uint n, uint m) {
	if(!laguerreValues.is_empty() && laguerreValues.n_rows > n && laguerreValues.n_cols > m){
		return resLaguerre.slice(m).col(n);
	}
	else
		return NULL;
}
