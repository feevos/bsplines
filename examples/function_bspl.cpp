// Demonstration of creating and manipulating a b-spline function of order k. 


#include <bsplines/f_bspline.hpp>
#include <iostream>
#include <random> 

using namespace std; 


int main(){

	// Construct a vector of breakpoints. 
	vector<double> breakpts {0., .25, .5, .75, 1.};

	// Define order of b-spline
	int k=5; 

	// Define bspline basis: 
	bspline_basis mybasis(breakpts,k); 

	// Create a set of (random) coefficients: 
        std::random_device rd;
        std::mt19937 gen;
	gen.seed(123456789); 
	//gen.seed(rd());		// True random seed: 

	uniform_real_distribution<> unif_real(-1.0,1.0);

	vector<double> ai; 
	
	for (auto j=0; j < mybasis.get_nbasis(); ++j)
		ai.push_back(unif_real(gen)); 

	// Construct function b-spline: 
	f_bspline myf(ai,mybasis); 

	// Print function and it's derivatives:  
	for (double x=breakpts.front(); x <= breakpts.back(); x+=0.01)
		{ 
		cout << x << " "; 

		// Output derivatives. 
		for (auto j=0; j < k; ++j)
			cout << myf.Djf_eval(j,x)<< " "; // Zeroth order derivative = function. 

		// Output integrals (by default 1 up to 3rd order): 
		for (auto j=1; j <= 3; ++j)
			cout << myf.Intf_eval(j,x)<< " ";  
		cout << endl; 

		}





}
