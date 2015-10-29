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
	for (double x=breakpts.front(); x <= breakpts.back(); x+=0.01)
		{ 
		cout << x << " "; 
		for (auto j=0; j < mybasis.get_nbasis(); ++j)
			cout << mybasis.get_DjBix(0,j,x)<< " "; // Zeroth order derivative = bspline basis. 
		cout << endl; 


		}





}
