#include <iostream>
#include <armadillo>
#include <fstream>


#define omega 1
#define mass 9.109e-31
#define hBar 1.0545718e-34	
#define pi 3.14159265358979323846
#define omega0 sqrt(mass * omega / hBar)

unsigned long factorial(unsigned long f)
{
	if(f < 0){
		exit(-1);
	}
    if ( f == 0 ) 
        return 1;
    return(f * factorial(f - 1));
}

int main(){
	//calcul des coefficients

	double inf = -5, sup = 5, nbX = 201;
	int n = 5, i;
	std::ofstream PsiFile("toPlot/Psis201.dat"), HermitFile("toPlot/hermit.dat");

	arma::mat res = arma::mat(nbX, n, arma::fill::ones);
	arma::colvec X = arma::linspace(inf, sup, nbX);
	res.col(1) = 2 * X;
	for(int i = 2; i < n; i ++){
		// printf("%d\n", i );
		res.col(i) = 2 * X % res.col(i-1) - 2 * (i-1) * res.col(i-2);
	}
	// X.print();
	HermitFile << res<< std::endl;

	//calcul des solutions

	arma::mat sol = arma::mat(nbX, n, arma::fill::ones);
	for(i = 0; i < n; i++){
		sol.col(i) =  sqrt(omega0) / sqrt(
					 	pow(2, i) *
					 	factorial(i)*
					 	sqrt(pi)) *
					 	arma::exp(-pow(X, 2)/2) %
					res.col(i);

	}
	PsiFile << sol << std::endl;
}