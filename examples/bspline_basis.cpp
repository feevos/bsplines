// Creating a set of b-spline basis and calling them. 

#include <bsplines/bsplines.hpp>
#include <iostream>

using namespace std; 


int main(){

	// Construct a vector of breakpoints.  See bsplines.hpp file for alternative definition. 
	vector<double> breakpts {0., .25, .5, .75, 1.};


	// Define order of b-spline
	int k=5; 

	// Define bspline basis: 
	bspline_basis mybasis(breakpts,k); 

	// Print all basis for x \in breakpts
	/* Derivative order of the b-spline basis. der_order <k
	for der_order = 0, we get the basis $B_i(x)$ 
	*/
	int der_order = 0;  
	for (double x=breakpts.front(); x <= breakpts.back(); x+=0.01)
		{ 
		cout << x << " "; 
		for (auto i=0; i < mybasis.get_nbasis(); ++i)
			cout << mybasis.get_DjBix(der_order,i,x)<< " "; // Zeroth order derivative = bspline basis. 
		cout << endl; 


		}





}
