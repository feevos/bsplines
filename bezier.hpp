/*

    Copyright (c) F.I.Diakogiannis  2015


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with popmcmc++. If not, see <http://www.gnu.org/licenses/>.

*/



/* File Description:
	
 	Implementation of Bezier functions. 

	References:

	[1] The nurbs book - Piegl & Tiller 2nd edition. 

*/


#ifndef _bezier_
#define _bezier_

#include <iostream> 
#include <vector> 
#include <cmath>
#include <algorithm>

#include "macros.h"



using namespace std; 

class bezier{

	private:
		int degree; 								/*! Degree of Bezier polynomial */
		
		vector <double> Bin; 							/*! Value of Bernstein basis at u: 0 <= u <=1 */
		vector <double> coeffs; 						/*! Coefficients of Bezier functions */  

		double u_saved; 							/*! Temporary stored x value,  initialize to NAN. Usage: avoide unnecessary function calls */
	
		/*! Based on algorithm A1.3 - p.20 [1]. */  	
		void eval_Bernstein_basis(const int &n,  const double &u); 		/*! Evaluates non zero Bernstein basis of degree n at u \im [0,1] */

	public:
		bezier();
		bezier(vector<double> _coeffs);

		
		double get_Binu(const int &i, const double &u); 			/*! Value Bin(u) */
		// double get_DjBinu(const int &j, const int &i, const double &u); 	/*! Value d^j B_i(x)/ dx^j */ 
		int get_degree(); 							/*! Returns degree of Bernstein basis */
		vector<double> get_coeffs();
		double eval(const double &u);						/*! Evaluates value of Bezier function at u */
};

bezier:: bezier (vector<double> _coeffs) : coeffs(_coeffs) 
	{
	// Since 0 <= u <= 1, this indicates no evaluation has taken place yet. 
	u_saved = -1.;
	degree=coeffs.size()-1; 
}

/*! 
Computes all Bernstein polynomials of degree n, at value u: 0 <= u <=1, 
and stores them in vector Bin. 
*/
void bezier:: eval_Bernstein_basis(const int &n,  const double &u){
	/*
	// Sanity check: 
	if (u < 0.0 || u >1 ) {
	cerr<< "Gave parameter value outside domain of definition for Bernstein basis, aborting ..." <<endl; 
	cerr<< "Value u:= " << u << endl; 
	DEBUG (u);
	}
	*/
	// Initialize everything to zero 
	vector<double> Bin_temp(n+1,0.0);
	
	Bin_temp[0]=1.;
	double u1 = 1. - u;

	// Temporary variables. 
	double saved, temp; 
	for (int j=1; j<=n; ++j)
		{
		saved = 0.0;
		for (int k=0; k<j; ++k)
			{
			temp = Bin_temp[k];
			Bin_temp[k] = saved + u1*temp;
			saved = u*temp;
			}
		Bin_temp[j] = saved;
		}

	Bin = std::move(Bin_temp);
	//Bin = Bin_temp;

}

double bezier:: get_Binu(const int &i, const double &u){
	// Sanity check: 
	if (u < 0.0 || u >1 ) {
	cerr<< "Gave parameter value outside domain of definition for Bernstein basis, aborting ..." <<endl; 
	cerr<< "Value u:= " << u << endl; 
	DEBUG (u);
	}else if (i <0 || i > degree){
	cerr<< "Request for evaluation of Bernstein basis out of range: i: in [0,n], aborting ..." <<endl; 
	cerr<< "Value u:= " << u << endl; 
	DEBUG (u);
	}


	if (u == u_saved){
		return 
			Bin[i];
	} else {

	eval_Bernstein_basis(degree,u);
		return 
			Bin[i];
	}
}

int bezier:: get_degree()
	{
	return 
		degree;
}

vector<double> bezier::get_coeffs(){
	return 
		coeffs;
}
double bezier:: eval(const double &u){
	
	double tempsum = 0.0;
	for(int i=0; i < (int) coeffs.size(); ++i) // Change (int) coeffs.size() --> degree+1
		tempsum += coeffs[i] * get_Binu(i, u);

	return 
		tempsum;

}





#endif
