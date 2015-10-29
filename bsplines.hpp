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
	
This file contains definition of the bspline basis class and bspline functions (pp polynomials). 
Sofisticated algorithms are based mainly in the book "The NURBS Book" - Piegl, Tiller, 2nd edition. 


	References: 
		[1] The NURBS Book, Piegl - Tiller, 1996, 2nd edition 
		[2] Curves and surfaces for CAGD, G. Farin, 2002. 
*/
 

#ifndef _bspline_basis_
#define _bspline_basis_

#include <iostream> 
#include <utility>
#include <algorithm>

#include "macros.h"

using namespace std; 




class bspline_basis{

	private:
	//public:
		int k; 						/*! order of Bspline basis */
		int nbreak; 					/*! Dimension of breakpoints vector. */
		int nknots; 					/*! Dimension of knots vector. */
		int nbasis; 					/*! Number of basis functions */
		int nderivs; 					/*! Total number of derivatives evaluated. Allow (nderivs < k) ?  */

/*add*/		int nintegrals; 				/*! Total number of evaluated integrals. Default: 1 */		

		/* Single knot <--> multiplicity =1 */
		vector<pair<double,int> > knots_wth_mult; 	/*! Vector of knots, each entry is the pair (knot value, multiplicity)*/
		vector<double> breakpts;			/*! Represents strictly increasing values of knots, excluding all multiplicities */
		vector<double> knots; 				/*! Raw knot vector of BSpline basis, it may contain multiple entries due to multiplicity */
		
		bool greville_evaluated; 				/*! Logical operator, informs if Greville abscissa have been evaluated */
		vector<double> greville; 			/* Contains values of Greville abscissa according to definitions [2] */

		vector<double> Bix_nonzero; 			/*! size: (k,1). Stores  nonzero components of bspline basis */
		vector<double> Bix; 				/*! size: nbasis Stores  all  components of bspline basis.  Not necessary - remove? */
		
		vector<double> DjBix_nonzero; 			/*! size: (nderivs+1,k). Stores  the jth derivative of the ith bspline basis function at x, nonzero terms only */
		vector<double> DjBix; 				/*! Stores  the jth derivative of the ith bspline basis function at x */


		// Need custom copy and assignment constructor for this, I think. 
		int i_saved; // initialize to negative; 					/*! Saved index i of Bix evaluation, avoids repetetion of find_nonzero_interval */ 
		double x_saved; // initialize to NAN; 						/*! Temporary stored x value, used for avoiding unnecessary function calls */
		bool eval_Bix_flag; 	// Informs program if previous calls to get_Bix took place 
		bool eval_DjBix_flag;   // Informs program if previous calls to get_DjBix took place

		/* Experimental */
		bool eval_DjBix_lim_flag; 	// Informs program if previous calls to get_DjBix_lim took place


	public:
		/* Experimental */
		int find_knot_span_of_x_lim(const double &x, const char * lim); 	/*! Returns integer i: t_i < x <= t_{i+k}. Upon call it stores i in i_saved */		
		const char * low = "low";
		const char * high = "high";



		int find_knot_span_of_x(const double &x); 					/*! Returns integer i: t_i <= x < t_{i+k}. Upon call it stores i in i_saved */		
		int find_knot_index (const double &x );						/*! Returns integer i that corresponds to the first instance of x in knot vector, t_i==x */
		pair<int,int> find_nonzero_basis_at_x(const double &x); 			/*! Returns first, last index of nonzero basis B_i(x) at particular x. */
		pair<int,int> find_base_nonzero_interval(const double &x); 			/*! Returns first (i) , last (i+k) index of knots t_i at particular x. */


	private:
		/* !ESSENTIAL ROUTINES FOR EVALUATION! Add as optional argument another knot vector for use in evaluation of integrals */
		void eval_nonzero_basis(const int &i, const double &x); 	/*! Evaluates non zero basis functions at x */
		void eval_nonzero_Djbasis(const int &i, const double &x);	/*! Evaluates non zero basis functions and their derivatives at x */

		void eval_Bix(const int &i, const double &x);			/*! Evaluates all basis functions at x */

		int idx_nonzero(const int &i, const int &j);			/*! collective index for DjBix_nonzero 2D --> 1D */
		int idx_all(const int &i, const int &j); 			/*! collective index for DjBix 2D --> 1D */
		void eval_DjBix(const int &i, const double &x);				/*! Evaluates all basis functions and their derivatives at x */


	public:
		bspline_basis(){};
		virtual ~bspline_basis(){};
		/*! Constructor from custom-multiplicity knot vector */
		bspline_basis(const vector<pair<double,int> > &_knots_wth_mult, const int &_k);
		/*! Default clamped  knot vector constructor */
		bspline_basis(const vector<double> &_breakpts, const int &_k);
		//bspline_basis(const vector<double> &_knots, const int &_k);	/*! Need to add this constructor */

		void set_nderivs (const int & new_nderivs);			/*! Sets total number of derivatives to be evaluated. Default nderivs = min(3,k-1) */
		void set_nintegrals(const int & new_nintegrals); 		/*! Sets total number of Bix integrals. Default = 1 */

		
		int get_nbasis();						/*! Returns number of Bix basis functions */
		int get_nbreak();						/*! Returns dimension of breakpoints vector. */
		int get_order();						/*! Returns number of Bix basis functions */
		vector<pair<double,int> >  get_knots_wth_mult();		/*! Returns knot vector with multiplicities */
		vector<double> 	get_knots();					/*! Returns knot vector in vector<double> format, with possible multiple knots */
		vector<double> 	get_breakpts();					/*! Returns breakpoints vector */
		vector<double> get_abscissa();					/*! If necessary, evaluates and returns vector of Greville abscissa */


		void knots_TO_knots_wth_mult();					/*! Pass information from vector<double> to vector<pair<double,int> > */
		void knots_wth_mult_TO_knots();					/*! Pass information from vector<pair<double,int> > to vector<double> */
		void set_knots(vector<double> & knots_new);			/*! Update knot vector in vector<double> form  */ 
		void set_knots(vector<pair<double,int> > & knots_wth_mult_new);	/*! Update knot vector in vector<pair<double,int>> form  */ 
		


		/* Experimental */
		double get_DjBix_lim(const int &j, const int &i, const double &x, const char * lim); 	/*! Value d^j B_i(x)/ dx^j for upper or lower limit.  */

		/* Evaluation functions */
		double get_Bix(const int &i, const double &x); 			/*! Value B_i(x) */
		double get_DjBix(const int &j, const int &i, const double &x); 	/*! Value d^j B_i(x)/ dx^j */
		void eval_greville();						/*! Evaluates all Greville abscissa according to [2] */


		// TO BE ADDED, requires vectors of fixed coefficients to be stored for faster evaluation. 
		double get_IntjBix(const int &j, const int &i, const double &x); /*! Value of j nested integrals \int_0^x B_i(u) du */


};

/*!
 Returns integer i: t_i < x <= t_{i+k} OR t_i <= x < t_{i+k}.  This is used when one wants to calculate limiting values of B_ix functions and their derivatives between different knot points.
 It is necessary especially for cases of non continuous derivatives of b-spline basis.   
*/		
int bspline_basis::find_knot_span_of_x_lim(const double &x, const char * lim) 	
	{
	// Sanity check 
	if( x < knots.front() || x > knots.back()){
		DEBUG(x);
                cerr<< "Value x outside of knot interval, aborting ... " <<endl;
                throw;
        }

	//Case t_i <= x < t_{i+k}
	if (lim == low){
        /**
        In case x=t_{m}, special treatment, see p.68 the NURBS book. 
        */
        if( x== knots.back() ){
        i_saved = nknots - k -1; // Second last knot position.  t_{last -1 } 

        return
                i_saved;
        }


        // Upper bound is what I need to avoid index multiplicity. 
        auto lower = std::upper_bound(knots.begin(),knots.end(),x); // This method goes beyond last index if it is used. 
        i_saved = distance(knots.begin(),lower); // Corresponds to i+1 value, i.e. t_i <= x < t_{i+1}
        i_saved--; // now i_saved --> i. 

	
	return 
		i_saved; 

	
	//Case t_i < x <= t_{i+k}
	}else if (lim == high){
	// Offer 


        /**
        In case x=t_{m}, special treatment, see p.68 the NURBS book. 
        */
        if( x == knots.front() ){
        i_saved = k -1; // first breakpoint position. 

        return
                i_saved;
        }


        // Lower bound is what I need to avoid index multiplicity from above. 
        auto lower = std::lower_bound(knots.begin(),knots.end(),x); // This method goes beyond last index if it is used. 
        i_saved = distance(knots.begin(),lower); // Corresponds to i+1 value, i.e. t_i < x <= t_{i+1}
        i_saved--; // now i_saved --> i. 

	return 
		i_saved; 

	}
}


/*!
Returns integer i, such that: 
	
	t_i <= x < t_{i+k}


Special treatment, if x==knots.back(), where then:

	t_i < x <= t_{i+k}, i: nknots-k-2

*/
// TESTED, OK! 
int bspline_basis::find_knot_span_of_x(const double &x)
	{
	
	 if( x < knots.front() || x > knots.back()){
		DEBUG(x);
                cerr<< "Value x outside of knot interval, aborting ... " <<endl;
                throw;
        }


        /**
        In case x=t_{m}, special treatment, see p.68 the NURBS book. 
        */
        if( x== knots.back() ){
        i_saved = nknots - k -1; // Second last knot position.  t_{last -1 } 

        return
                i_saved;
        }


        // Upper bound is what I need to avoid index multiplicity. 
        auto lower = std::upper_bound(knots.begin(),knots.end(),x); // This method goes beyond last index if it is used. 
        i_saved = distance(knots.begin(),lower); // Corresponds to i+1 value, i.e. t_i <= x < t_{i+1}
        i_saved--; // now i_saved --> i. 

	
	return 
		i_saved; 
}

/*!
This function returns the index position of the FIRST knot element equal to the value x. 
If x is not a knot element, throws an error. 
*/
int bspline_basis::find_knot_index (const double &x )
	{
        auto it = std::find(knots.begin(),knots.end(),x);
        int idx = distance(knots.begin(),it);
        if (it == knots.end() )
        	{
                cout << "Knot value not contained in knot vector, aborting ... " << endl;
		DEBUG(idx);
                throw 0;
                }

        return 
		idx;


}




/*!
For an x in [t_i,t_{i+1}], the only nonzero basis functions are k basis functions, namely: {B_{i-k+1}(x), ... , B_i(x)}. 
This function evaluates pair<istart, iend> such that: 

	t_{iend} <=  x < t_{iend+k}
and 

	istart = iend-k+1
*/
// TEST IT
pair<int,int> bspline_basis::find_nonzero_basis_at_x(const double &x)
	{
	find_knot_span_of_x(x);
	pair<int,int> temp = make_pair(i_saved-k+1,i_saved);
		
	return 
		temp; 

}

/*!
For an x in [t_i,t_{i+1}], the corresponding basis B_i(x) is nonzero in [t_i,t_{i+k})	
This function evaluates pair<istart, iend> such that: 

	t_{istart} <=  x < t_{istart+k}
and 
	iend = istart + k
*/
// TEST IT
pair<int,int> bspline_basis::find_base_nonzero_interval(const double &x)
	{
	find_knot_span_of_x(x);
	pair<int,int> temp = make_pair(i_saved,i_saved + k);

	return temp; 
	
}




/*!
Input: value x for evaluation, index i: t_{i} <= x < t_{i+k}, except if x==knots.back();
On exit, the values of the nonzero bspline basis are stored in Bix_nonzero
*/
void bspline_basis::eval_nonzero_basis(const int &i, const double &x){

	// Calculate interval where i lies, i.t. knots[i] <= x < knots[i+1];
	//int i=find_nonzero_basis_at_x(x).second; 

	if(x < knots.front() || x > knots.back()){
	DEBUG(x);
	std::cerr << "Value x outside of knot interval, aborting ..." <<endl; 
	throw 0;
	}

	vector<double> N(k,0.0), left(k,0.0), right(k,0.0);
	N[0]=1.0; 

	double saved,temp;
	for(int j=1; j < k; ++j)
		{
		left[j]	= x-knots[i+1-j];
		right[j]= knots[i+j]-x;
		saved=0.0;
		for(int r=0; r<j; ++r)
			{
			temp 	= N[r]/(right[r+1]+left[j-r]); 
			N[r] 	= saved + right[r+1]*temp;
			saved 	= left[j-r]*temp;
			}
		N[j]=saved;
		}
	// Move the value to the class member Bix_nonzero
	Bix_nonzero=std::move(N);

}


/*!
Input: value x for evaluation, 
index i: t_{i} <= x < t_{i+k}, except if x==knots.back(), 
nderivs: Order up to which derivatives are evaluated (rarely we will need to evaluate all of them). 
Restriction: nderifs < k 

On exit, the values of the nonzero bspline basis as well as their derivatives are stored in Bix_nonzero vector. 
*/
void bspline_basis::eval_nonzero_Djbasis(const int &i, const double &x){
	
	if(nderivs >= k){
	cerr << "Request for derivative order equal or higher than bspline order, aborting ...." << endl; 
	throw 0;
	}

	// Sanity check. 
	if(x < knots.front() || x > knots.back()){
	DEBUG(x);
	std::cerr << "Value x outside of knot interval, aborting ..." <<endl; 
	throw 0;
	}

	// Collective index 2D-->1D
	auto idx_ndu = [&] (const int &ii, const int &jj)->int{
		return 	
			ii+k*jj;
	};

	
	auto idx_ders = [&] (const int &ii, const int &jj)->int{
		return 	
			ii+(nderivs+1)*jj;
	};

	
	auto idx_a = [&] (const int &ii, const int &jj)->int{
		return 	
			ii+2*jj;
	};

	/* local objects */
	vector<double> ndu(k*k,0.0), a(2*k,0.0) ;
	ndu[idx_ndu(0,0)]=1.0;

	vector<double>  left(k,0.0), right(k,0.0);

	/* loop variables */
	
	int j, r, m, j1, j2;
	double saved,temp;

	for( j=1; j<k; ++j)
		{
		left[j]		= 	x-knots[i+1-j];
		right[j]	=	knots[i+j]-x;
		saved=0.0;
		for(int r=0; r<j; ++r)
			{
			ndu[idx_ndu(j,r)]	=	right[r+1] + left[j-r];
			temp 			=	ndu[idx_ndu(r,j-1)] / ndu[idx_ndu(j,r)];

			ndu[idx_ndu(r,j)]	=	saved + right[r+1]*temp;
			saved			=	left[j-r]*temp;
			}
		
		ndu[idx_ndu(j,j)] = saved;
		}
	
	// Matrix that stores the derivatives. First column: stores values of B_{i,k}(x)
	vector<double> ders((nderivs+1)*k,0.0);
	/* Load the basis functions */
	for( j=0; j<k; j++ ) 
		ders[idx_ders(0,j)] = ndu[idx_ndu(j,k-1)];		
		
	/* This section computes the derivatives */
	for( r=0; r<k; r++) 
		{
		/* Alternate rows in array a */
		int s1=0, s2=1;
		a[idx_a(0,0)] = 1.0;
		/* loop to compute mth derivative */
		for( m=1; m<= nderivs; ++m)	
			{
			double d=0.0;
			int rk	= r-m;
			int pk	= (k-1) - m;		
			if (r >= m)
				{
				a[idx_a(s2,0)] = a[idx_a(s1,0)] / ndu[idx_ndu(pk+1,rk)];
				d = a[idx_a(s2,0)]*ndu[idx_ndu(rk,pk)];	
				}

			rk >= -1 	? 	j1 = 1 	: 	j1 = -rk;			
			r-1 <= pk 	? 	j2 = m-1:	j2 = (k-1)-r;

			for( j=j1; j<= j2; ++j)
				{
				a[idx_a(s2,j)] = (a[idx_a(s1,j)] - a[idx_a(s1,j-1)]) / ndu[idx_ndu(pk+1,rk+j)];
				d +=  a[idx_a(s2,j)] * ndu[idx_ndu(rk+j,pk)];
				}
			if (r <= pk)
				{
				a[idx_a(s2,m)] = - a[idx_a(s1,m-1)] / ndu[idx_ndu(pk+1,r)];
				d += a[idx_a(s2,m)] * ndu[idx_ndu(r,pk)];
				}
			ders[idx_ders(m,r)]=d;
			/* switch rows */
			j=s1;
			s1=s2;
			s2=j;
			}
		}

		/* Multiply through by the correct factors */
		r = k-1;
		for(m=1; m <= nderivs; ++m)
			{
			for (j=0; j < k; ++j)
				ders[idx_ders(m,j)] *= r;

			r *= ((k-1)-m);
			}
			
	DjBix_nonzero = std::move(ders);

}		


/*!
Wrapper function for eval_nonzero_basis. Passes nonzero basis values to vector<double> Bix. 
*/
void bspline_basis::eval_Bix (const int &ii, const double &x){
	
        //pair<int,int> i_start_end = find_nonzero_basis_at_x(x);
        pair<int,int> i_start_end = make_pair(ii-k+1,ii);
        
	// Evaluate all nonzero entries. for this index. 
        eval_nonzero_basis(i_start_end.second, x);

        // Initialize (to zeros) temporary vector of dimension nbasis 
        vector<double> Bix_temp(nbasis,0.0);
	        
        // Pass nonzero entries to temporary vector 
	for(int j= i_start_end.first; j <= i_start_end.second; ++j)
        	Bix_temp[j] = Bix_nonzero[j-i_start_end.first];

	// move temporary vector to Bix
	Bix=std::move(Bix_temp);
}


int bspline_basis::idx_nonzero(const int &i, const int &j)
	{

        int temp = i + (nderivs+1)*j;
        if (temp > (nderivs+1)*k){
	cerr << "Error in function: " << __PRETTY_FUNCTION__ << ", in file: " << __FILE__ << " --- line number: " << __LINE__ <<endl; 
	cerr << "Aborting ... " << endl; 
	
        throw 0;
	}
         
	return 
		temp;

}


int bspline_basis::idx_all (const int &i, const int &j)
	{
        int temp = i + (nderivs+1)*j;
        if (temp > (nderivs+1)*nbasis){
	cerr << "Error in function: " << __PRETTY_FUNCTION__ << ", in file: " << __FILE__ << " --- line number: " << __LINE__ <<endl; 
	DEBUG(i);
	DEBUG(j);
	DEBUG(temp);
	cerr << "Aborting ... " << endl; 
	
        throw 0;
	}
         
	return 
		temp;

}

// TESTED -  OK      
void bspline_basis::eval_DjBix (const int &ii, const double &x)
	{
        pair<int,int> i_start_end = make_pair(ii-k+1,ii); // range of indices for nonzero basis Bix. 

        // Evaluate all nonzero entries. for this index. 
        eval_nonzero_Djbasis(i_start_end.second, x);

        // Pass entries to temporary vector of dimension nbasis 
        vector<double> temp_DjBix( nbasis * (nderivs+1),0.0);
                for(int i=0; i<=nderivs; ++i)
                        for(int j= i_start_end.first; j <= i_start_end.second; ++j)
                                temp_DjBix[idx_all(i,j)] = DjBix_nonzero[idx_nonzero(i,j-i_start_end.first)];

        //DjBix=std::move(temp_DjBix);
        DjBix=temp_DjBix;
}

/*!
Constructs knot vector from a given strictly increasing breakpoints sequence. In this, the multiplicity of first and last knot points 
is equal to the order k of the bspline, such as, any bspline function (curve), passes from end points. 
*/
bspline_basis::bspline_basis(const vector<double>&_breakpts, const int &_k):  k(_k), breakpts(_breakpts){

	
	// Set default number of nderivs: 
	nderivs = k-1; 
	nbreak = breakpts.size();

	// Perform sanity check, that all breakpoint elements are UNIQUE
	vector<double> temp = _breakpts;
        auto last = std::unique(temp.begin(),temp.end());
        temp.erase(last,temp.end());

        if (temp.size() != nbreak )
		{
		cerr<< "Gave verctor of NON unique elements in breakpoints, aborting ..." << endl; 
		DEBUG(temp.size());
		throw 0;
		}


	
        knots_wth_mult.resize(nbreak);
	// Pass values of breakpoints to knots. 
        for(int i=1; i<nbreak-1; i++)
        	{
                knots_wth_mult[i].first=_breakpts[i];
                knots_wth_mult[i].second=1;
                }

        // Construct clamped 
        knots_wth_mult.front().first=_breakpts.front();
        knots_wth_mult.front().second = k;
        knots_wth_mult.back().first=_breakpts.back();
        knots_wth_mult.back().second= k;
	
	knots.resize(nbreak+2*(k-1));

	for(int i=0; i<k; ++i){
		knots[i]=breakpts.front();
		knots[nbreak+2*(k-1)-1-i] = breakpts.back();
	}
	for(int i=0; i<nbreak; ++i)
		knots[i+k-1]=breakpts[i];
	nknots=knots.size();
	nbasis = nknots - k; // Holds for arbirtrary multipliciy. 

	/* Sanity check - passed. 
	// Total number of multiplicity affected basis.  
	double tempsum=0.;
	for(int i=0; i<nbreak; i++)	
		tempsum += knots_wth_mult[i].second-1.;
	*/

	// initialize i_saved, x_saved:
	i_saved = -10;
	x_saved = NAN;
	// Initialize flags for matrix evaluation to false. 
	eval_Bix_flag = false; 	
	eval_DjBix_flag = false;   
	eval_DjBix_lim_flag = false; 


	// Calculate and store Greville:  or not ... ? 
	greville.resize(nbasis);
	greville_evaluated = false;


	// More sanity checks. 
	if(breakpts[0] ==  breakpts[1] || breakpts[nbreak-2] == breakpts[nbreak-1]){
	cerr<< "Multiplicity higher than k in first/last knots is not supported, aborting ... "<<endl;
	cerr<< "Error in function: " << __PRETTY_FUNCTION__ << ", in file: " << __FILE__ << " -- line number: " << __LINE__<<endl;
	throw 0;
	} 


	if ( knots_wth_mult[0].second != k ){
	cerr<< "Multiplicity != k in first/last knots is not supported, aborting ... "<<endl;
	cerr<< "Error in function: " << __PRETTY_FUNCTION__ << ", in file: " << __FILE__ << " -- line number: " << __LINE__<<endl;
	throw 0;
	} 

	bool sorted = is_sorted(breakpts.begin(),breakpts.end(), [](double &x, double &y) { if (x < y){ return true; }else if (y <= x){return false;}});
	if (sorted==false){
	DEBUG(0.0);
	cerr<<"breakpoints vector is not strictly increasing,  aborting ..."<<endl;
	throw 0;
	}

	sorted = std::is_sorted(knots.begin(),knots.end());
	if (sorted==false){
	cerr<<"Knot vector is not sorted, error in function: " << __PRETTY_FUNCTION__ << ", in file: " << __FILE__ << " -- line number: ";
	cerr << __LINE__ <<endl; 
	cerr << "Aborting ..."<<endl;
	throw 0;
	}

}


/*!
bspline basis constructor for a custom knot vector. 
*/
bspline_basis::bspline_basis(const vector<pair<double,int> >&_knots_wth_mult, const int &_k):  k(_k), knots_wth_mult(_knots_wth_mult){


	// Set default number of nderivs: 
	nderivs = k-1; 
	// Construct knots vector: 
	nbreak = knots_wth_mult.size();
	
	// Construct breakpoints and knot vector. 
	breakpts.resize(nbreak);
	knots.reserve(nbreak*k);
	for(int i=0; i<nbreak; ++i)
		{
		if (knots_wth_mult[i].second > k){
		DEBUG(knots_wth_mult[i].second);
		std::cerr<< "Multiplicity of " << i << " knot value higher than order k, aborting... "<<endl;
		throw 0;
		}
		breakpts[i] = knots_wth_mult[i].first;
		knots.push_back(knots_wth_mult[i].first);
		for(int j=1; j < knots_wth_mult[i].second; j++)
			{
			knots.push_back(knots_wth_mult[i].first);
			}	
	}
	nknots=knots.size();
	// Now that you know # of nknots, store number of basis: 
	nbasis = nknots - k; // Holds for arbirtrary multipliciy. 

	/* Sanity check - passed. 
	// Total number of multiplicity affected basis.  
	double tempsum=0.;
	for(int i=0; i<nbreak; i++)	
		tempsum += knots_wth_mult[i].second-1.;
	*/

	// initialize i_saved, x_saved:
	i_saved = -10;
	x_saved = NAN;

	// Initialize flags for matrix evaluation to false. 
	eval_Bix_flag = false; 	
	eval_DjBix_flag = false;   
	eval_DjBix_lim_flag = false; 


	// Calculate and store Greville:  or not ... ? 
	greville.resize(nbasis);
	greville_evaluated = false;


	// More sanity checks. 
	if(breakpts[0] ==  breakpts[1] || breakpts[nbreak-2] == breakpts[nbreak-1]){
	cerr<< "Multiplicity higher than k in first/last knots is not supported, aborting ... "<<endl;
	cerr<< "Error in function: " << __PRETTY_FUNCTION__ << ", in file: " << __FILE__ << " -- line number: " << __LINE__<<endl;
	throw 0;
	} 


	if ( knots_wth_mult[0].second != k ){
	cerr<< "Multiplicity != k in first/last knots is not supported, aborting ... "<<endl;
	cerr<< "Error in function: " << __PRETTY_FUNCTION__ << ", in file: " << __FILE__ << " -- line number: " << __LINE__<<endl;
	throw 0;
	} 

	bool sorted = is_sorted(breakpts.begin(),breakpts.end(), [](double &x, double &y) { if (x < y){ return true; }else if (y <= x){return false;}});
	if (sorted==false){
	DEBUG(0.0);
	cerr<<"breakpoints vector is not strictly increasing,  aborting ..."<<endl;
	throw 0;
	}

	sorted = std::is_sorted(knots.begin(),knots.end());
	if (sorted==false){
	cerr<<"Knot vector is not sorted, error in function: " << __PRETTY_FUNCTION__ << ", in file: " << __FILE__ << " -- line number: ";
	cerr << __LINE__ <<endl; 
	cerr << "Aborting ..."<<endl;
	throw 0;
	}

}

/*! 
Sets total number of derivatives to be evaluated. Default nderivs = 3 
*/
void bspline_basis::set_nderivs (const int & new_nderivs)
	{
	nderivs = new_nderivs;
}			

/*! 
Sets total number of Bix integrals. Default = 1 
*/
void bspline_basis::set_nintegrals(const int & new_nintegrals)
	{	
	nintegrals = new_nintegrals; 
}

int bspline_basis::get_nbasis()
	{
	return 	
		nbasis;
}


int bspline_basis::get_nbreak()
	{
	return 	
		nbreak;
}





int bspline_basis::get_order()
	{
	return 	
		k;
}
/*! Returns knot vector with multiplicities */
vector<pair<double,int> >  bspline_basis::get_knots_wth_mult()			
	{
	return knots_wth_mult;
}
/*! Returns knot vector in vector<double> format, with possible multiple knots */
vector<double> 	bspline_basis::get_knots()					
	{
	return knots; 
}
/*! Returns breakpoints vector */
vector<double> 	bspline_basis::get_breakpts()					
	{
	return breakpts;
}



void bspline_basis::knots_TO_knots_wth_mult(){					/*! Pass information from vector<double> to vector<pair<double,int> > */

	vector<pair<double,int> > _knots_wth_mult;
	_knots_wth_mult.reserve(nknots);

        int counter=1;
        double temp=knots[0];
        for(int i=1; i<nknots; i++)
        	{
                if( knots[i] == temp){
                counter++;
                }else{
                _knots_wth_mult.emplace_back(temp,counter);
                temp = knots[i];
                counter=1;
                }

                }
         // Pass last knot value.        
         _knots_wth_mult.emplace_back(temp,counter);
	// Pass values to knots_wth_mult composite form 
	knots_wth_mult = _knots_wth_mult;

}

void bspline_basis::knots_wth_mult_TO_knots(){

	vector<double> knots_new;
	// Reserve as if each knot value had multiplicity equal to order. 
	int nreserve = nbreak * k; 
	knots_new.reserve (nreserve);
        for(int i=0; i<nbreak; ++i)
                {
                knots_new.push_back(knots_wth_mult[i].first);
                for(int j=1; j < knots_wth_mult[i].second; j++)
                        {
                        knots_new.push_back(knots_wth_mult[i].first);
                        }
                }



       knots=knots_new;


}

/*! 
Update knot vector, in vector<double> form. 
*/
void bspline_basis::set_knots(vector<double> & knots_new)			/*! Update knot vector */ 
	{
	// Pass new knot vector
	knots=knots_new;
	// Update class info
	nknots=knots.size();
	nbasis = nknots - k; // Holds for arbirtrary multipliciy. 

	// update knots_wth_mult
	knots_TO_knots_wth_mult();
	nbreak = knots_wth_mult.size();

	breakpts.resize(nbreak);

	for(int i=0; i<nbreak; i++)
		breakpts[i] = knots_wth_mult[i].first;

	greville_evaluated = false; 
}


/*! 
Update knot vector, in vector<double> form. 
*/
void bspline_basis::set_knots(vector<pair<double,int> > & knots_wth_mult_new)			/*! Update knot vector */ 
	{
	// Pass new knot vector
	knots_wth_mult=knots_wth_mult_new;
	// Update class info
	nbreak=knots_wth_mult.size();
	knots_wth_mult_TO_knots();
	nknots=knots.size();
	nbasis = nknots - k; // Holds for arbirtrary multipliciy. 

	// update knots_wth_mult
	knots_wth_mult_TO_knots();
	// Update breakpoints. 
	breakpts.resize(nbreak);
	for(int i=0; i<nbreak; i++)
		breakpts[i] = knots_wth_mult[i].first;

	
	greville_evaluated = false; 
}



 /*! Returns value B_i(x) */
double bspline_basis::get_Bix(const int &i, const double &x)
	{

	if (i<0 || i> nbasis){
	DEBUG(i);
	std::cerr<< "Index of Bix out of range, aborting ..." << endl; 
	throw 0;
	}

	if(x==x_saved && eval_Bix_flag==true){
		eval_DjBix_flag = false;
		eval_DjBix_lim_flag = false;
		return 
			Bix[i];
	}else{
		
		eval_DjBix_flag = false;
		eval_DjBix_lim_flag = false;
		// a. Find knot span of x: 
		find_knot_span_of_x(x); 
		// b. Evaluate all nonzero Bix and pass them to Bix: 
		eval_Bix(i_saved,x);
	
		// Used for optimization - avoids unnecessary cals to function eval_Bix
		x_saved=x;	// Store x for subsequent evaluations.
		eval_Bix_flag=true;
		return 
			Bix[i];

	}

}



/*! 
Returns value d^j B_i(x)/ dx^j. This is a wrapper function for eval_DjBix in order to avoid unnecessary calculations 

*/
double bspline_basis::get_DjBix_lim(const int &j, const int &i, const double &x,const char * lim)
	{

	if (i<0 || i> nbasis){
	DEBUG(i);
	std::cerr<< "Index of Bix out of range, aborting ..." << endl; 
	throw 0;
	} 

	if ( j> nderivs) {
	DEBUG(j);
	std::cerr<< "Requested derivative order greater than evaluated nderivs. Set variable nderivs to a higher value and recalculate, aborting ..." << endl; 
	throw 0;
	}
	if(x==x_saved && eval_DjBix_lim_flag==true){
		eval_Bix_flag = false;
		eval_DjBix_flag = false;
		return 
			DjBix[idx_all(j,i)];
	}else{
		eval_Bix_flag = false;
		eval_DjBix_flag = false;
		// a. Find knot span of x: 
		find_knot_span_of_x_lim(x,lim); 
		// b. Evaluate all nonzero Bix and pass them to Bix: 
		eval_DjBix(i_saved,x);
		
		// I will use these for optimization, don't know how yet ... 
		x_saved=x;		// Store x for subsequent evaluations.
		//i_saved=i;		// Store knot span i for possible subsequent evaluations.
		eval_DjBix_lim_flag = true;
		
		return 
			DjBix[idx_all(j,i)];
	}

}


/*! 
Returns value d^j B_i(x)/ dx^j. This is a wrapper function for eval_DjBix in order to avoid unnecessary calculations 

*/
double bspline_basis::get_DjBix(const int &j, const int &i, const double &x)
	{
	if (i<0 || i> nbasis){
	DEBUG(i);
	std::cerr<< "Index of Bix out of range, aborting ..." << endl; 
	throw 0;
	} 

	if ( j> nderivs) {
	DEBUG(j);
	std::cerr<< "Requested derivative order greater than evaluated nderivs. Set variable nderivs to a higher value and recalculate, aborting ..." << endl; 
	throw 0;
	}
	if(x==x_saved && eval_DjBix_flag==true){
		eval_Bix_flag = false;
		eval_DjBix_lim_flag = false;
		return 
			DjBix[idx_all(j,i)];
	}else{
		eval_Bix_flag = false;
		eval_DjBix_lim_flag = false;
		// a. Find knot span of x: 
		find_knot_span_of_x(x); 
		// b. Evaluate all nonzero Bix and pass them to Bix: 
		eval_DjBix(i_saved,x);
		
		// I will use these for optimization, don't know how yet ... 
		x_saved=x;		// Store x for subsequent evaluations.
		//i_saved=i;		// Store knot span i for possible subsequent evaluations.
		eval_DjBix_flag = true;
		
		return 
			DjBix[idx_all(j,i)];
	}
}


void bspline_basis::eval_greville()
	{
	greville.front() = knots.front();
	greville.back()  = knots.back();
	
	for(int i=1; i<nbasis-1; ++i)
		{
		double tempsum=0.0;
		for (int j=0; j < k-1; ++j)
			tempsum += knots[(i+1)+j];
		greville[i] = tempsum/(double(k-1));
		
		}
	greville_evaluated = true;

}

/*!
Returns vector containing Greville abscissa (x coordinate of control coefficiens a_i for a bspline function f(x) = a^i B_i(x) ). 
*/
vector<double>  bspline_basis::get_abscissa()
	{
	if (greville_evaluated){
		return 
			greville;
	}else{
		eval_greville();
		return 
			greville;
	}
}






#endif
