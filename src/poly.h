#include <iostream>
#include <armadillo>
#include <math.h>

    #define omega 1
    #define mass 9.109e-31
    #define hBar 1.0545718e-34  
    #define pi 3.14159265358979323846
    // #define omega0 sqrt(mass * omega / hBar)
    #define omega0 1
/**
 * \brief a Polynome is defined by its value in abscissa and its orderly.
 */

class Poly{
    
    protected:
        uint nH, nL, mL; // max indexes for hermit//Laguerre
        arma::mat resHermite,hermiteValues;
        arma::cube resLaguerre, laguerreValues;
        arma::vec zValues;
    public:
        Poly();
        Poly(uint nHermite, uint nLaguerre, uint mLaguerre, arma::vec Z);
        // unsigned long factorial(unsigned long f);

        //Getters
        uint getNH();
        uint getNL();
        uint getML();
        arma::mat getHermiteValues();
        arma::cube getLaguerreValues();
        arma::vec getZValues();

        //Setters
        void setNH(uint n);
        void setNL(uint n);
        void setM(uint m);
        void setHermiteValues(arma::mat values);
        void setLaguerreValues(arma::cube values);
        void setZValues(arma::vec values);

        /** 
         * \brief computes the values of n first Hermitian Polynoms at zVal values. 
         * \assigns this.polyValues
         * \assigns this.zValues
         * \assigns this.n
         * \param n the n first Hermitian will be computed.
         * \param zVal List of values to compute Hermitian
         */
        void calcHermite(uint n, arma::vec zVal);

        /**
        * \brief gives the values of the n-th polynom evaluated at this.zValues.
        * \warning this.calcHermite has to be called Before calling this.hermite().
        * \param n the n-th hermitian polynom's values will we be returned.
        */
        arma::vec hermite(uint n);

        /** 
         * \brief computes the values of n first Laguerre Polynoms at zVal values. 
         * \assigns this.polyValues
         * \assigns this.zValues
         * \assigns this.n
         * \assigns this.m
         * \param n the n first Laguerre will be computed.
         * \param zVal List of values to compute Laguerre
         */
        void calcLaguerre(uint n, uint m, arma::vec zVal);

        /**
        * \brief gives the the values of the n-th Laguerre polynom evaluated at this.zValues.
        * \warning this.calcLaguerre has to be called Before calling this.laguerre().
        * \param n the 1st index of laguerre polynomial values will we be returned.
        * \param m the 2nd index of laguerre polynomial values will we be returned.
        */  
        arma::vec laguerre(uint n, uint m);



};


long factorial(unsigned long f);