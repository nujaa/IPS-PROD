#include <armadillo>

class Basis{
    
    protected:
    public:
    	//members
    	int mMax;
        double br, bz;
    	arma::vec nMax;
    	arma::mat n_zMax;
    	// constructors
    	Basis(double Br, double Bz, int N, double Q);
    	//setters

    	//getters

    	//methods
        void setmMax(int N, double Q);
        void setNMax();
        void setN_zMax(int N, double Q);
        int n_ZMax(uint i, int N, double Q);

    	arma::vec rPart(arma::vec r, int m, int n); 
    	arma::vec zPart(arma::vec z, int nZ);
};