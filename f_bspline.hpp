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
        
This file contains definition for a bspline representation of a function f(x):= a^i B_i,k}(x).  
Algorithms are based mainly in the book "The NURBS Book" - Piegl, Tiller, 2nd edition. 


        References: 
                [1] The NURBS Book, Piegl - Tiller, 1996, 2nd edition 
                [2] Curves and surfaces for CAGD, G. Farin, 2002. 
*/


#ifndef _f_bspline_
#define _f_bspline_

#include "macros.h"
#include "bsplines.hpp"

using namespace std; 



// This gives the value of function f(x) in interval [rmin,...,rmax,rt] that is expanded in a bspline basis, according to: 
//
//      f(x) = \sum_i a_i * B[i,p](x)  
// 
// User supplies rmin,rmax,rt, a_i vector and optionally order k (=4 default) and gamma (=3.0 default) 
class f_bspline: public bspline_basis{
        private:
                vector<double> coeffs; // Vector of coefficients that are used in the basis of expansion. 
                int Nintegrations=3; // How many integrations will be requested. 
                vector<bspline_basis> basis_int; // Vector that holds the bspline basis of the integration. 
                vector<vector<double> > coeffs_int; // coefficients for each integral.                                                                                                         
       		bool base_mod; 		/*! Informs for all possible modifications of bspline_basis. If true, integrations of Bix need recalculated */

        public:
		void eval_integrals(const bool &_base_mod);	/*! Up on call repeats evaluation for coefficients and integrals basis */
		// Default constructor to be used inside other classes definition                                                                                                              
               	f_bspline(){};        
               	f_bspline(const vector<double> &_coeffs, const vector<pair<double,int> > &_knots_wth_mult, const int &_k);  
               	f_bspline(const vector<double> &_coeffs, const vector<double> &_breakpts, const int &_k);                                                       
               	f_bspline(const vector<double> &_coeffs, bspline_basis mybasis);                                                                                                      

		void set_base_mod(bool new_base_mod);

// -------------Move to individual class - bspline_manipulator ? or something like this -------------------------------------------------------------------------------------
		void refine_knot(vector<pair<double,int> > &X_wth_mult);	/*! Inserts new elements X into knot vector, and calculates new coefficients. coeffs and bspline_basis are updated */
		void refine_knot(vector<double> &_X);  				/*! Inserts new elements X into knot vector, and calculates new coefficients. coeffs and bspline_basis are updated */
		// NEEDS MORE WORK on TOL value/distance function. 
		void remove_knot(const int &r, const int &num); 		/*! Remove knot t_r (t_r != t_{r+1}), num times (try to remove it, with TOL = 1e.-10 */
		void f_Bez_decomp(int &nbezier, vector<vector<double>> &Qbez);	/*! Decomposes function to the number of Bezier segments nb, and coefficients Qbez */	
// -------------Move to individual class - bspline_manipulator ? or something like this -------------------------------------------------------------------------------------


               	void 	set_coeffs(vector<double> &_coeffs); 			/*! Sets a new vector of coeffs, without needing to re-define knot vector etc. */
               	vector<double> 	get_coeffs(); 					/*! Returns coefficients of f_bspline */ 
               	double 	f_eval(const double &x );
               	double 	Djf_eval(const int &j,  const double &x ); 		/*! Returns D^j (f) the jth derivative of f. j < order k */

		bspline_basis get_basis_int (const int &j); 
		vector<double>  get_coeffs_int (const int &j); 
               	double 	Intf_eval(const int &j,  const double &x ); 		/*! Returns \int_0^x ...\int_0^x f(x)dx , j integrations in total. */


};



void f_bspline:: eval_integrals(const bool &_base_mod)	/*! Up on call repeats evaluation for coefficients and integrals basis */
	{

	if (_base_mod == false)
		return ; // If not necessary to evaluate integrals, do nothing 
	
                /* ALGORITHM:  
                        1. Definte total number of objects:  -- allready done so in Nintegrations. 
                        2. Give vectors of class correct size - according to 1 above. 
                        3. Construct vector of bspline_basis. 
                        4. Construct first integral "manually" in order to avoid unnecessary copy of *this
                        5. Construct rest of integrals with a for loop. 

                */
                // STEP 2. 
                basis_int.resize(Nintegrations);
                coeffs_int.resize(Nintegrations);
                int ncoeffs=get_nbasis(); // Coefficients of integration are +1 for each nested integration. 
		vector<pair<double,int> > knots_wth_mult_ref =  get_knots_wth_mult();

                for(int i=0; i<Nintegrations; ++i)
                        {
                        vector<double> temp(ncoeffs+1+i,0.0);
                        coeffs_int[i] = std::move(temp);  // First element will be instantatiated to zeros  // CHECK THIS!
                        }
                // STEP 3. 
                //basis_int[0]=*this; // BASE bspline_basis - special treatment :)
                int kk = get_order();
                for(int i=0; i<Nintegrations; ++i) /* constrct rest of basis */
                        {
			
			// First fix  
			knots_wth_mult_ref.front().second = (kk+1)+i;
			knots_wth_mult_ref.back().second  = (kk+1)+i;
                        bspline_basis temp(knots_wth_mult_ref,(kk+1)+i); //Temporary objects, gets destroyed after each iteration. 
                        basis_int[i]= std::move(temp);
                        }

        	// STEP 4. 
                auto factor_first = [&](int &j)->double
                        {
                        return  // This is the basis BEFORE integration. 
                                (get_knots()[j+kk] - get_knots()[j])/double(kk);
                        };
                // Evaluate coefficients to be used for the interpretation of \int f_bspline 
                for(int i=0; i < get_nbasis(); ++i)
                        {
                        double tempsum=0.0;
                        for(int _idx=0; _idx <=i; ++_idx)
                                tempsum+= factor_first(_idx)*coeffs[_idx];

                        coeffs_int[0][i+1]=std::move(tempsum);
                        }


                // STEP 5. 
                // First index: integration level, second index: knot vector index. 
                auto factor = [&](int &i, int &j)->double
                        {
                        if(i<0){
                                cout<<"You fucked it mate, i must be >= 1, aborting ..." << endl;
                                throw(0);
                        }
                        return  // This is the basis BEFORE integration. 
                              (basis_int[i].get_knots()[j+(kk+i+1)] - basis_int[i].get_knots()[j] )/double(kk+i+1);
                        };


                for(int i=1; i<Nintegrations; ++i)
                                for(int ii=0; ii < basis_int[i-1].get_nbasis(); ++ii)
                                        {
                                        double tempsum=0.0;
                                        for(int _idx=0; _idx <=ii; ++_idx){
                                                int temp_idx=i-1;
                                                tempsum+= factor(temp_idx,_idx)*coeffs_int[i-1][_idx];
                                                }
                                        coeffs_int[i][ii+1]=tempsum;
                                        }

		
		base_mod = false;

}


f_bspline::f_bspline(const vector<double> &_coeffs, const vector<pair<double,int> > &_knots_wth_mult, const int &_k)
                :bspline_basis(_knots_wth_mult, _k), coeffs(_coeffs)
	{

		base_mod = true;
                if( (int) coeffs.size() != get_nbasis())
                        {
			DEBUG(0.0);                        
			std::cerr<< "Error in f_bspline constructor, missmatch of coefficients size with number of basis Bix, aborting ..." <<endl; 
			throw 0;
                        }



}

f_bspline::f_bspline(const vector<double> &_coeffs,const  vector<double> &_breakpts, const int &_k)
                :bspline_basis(_breakpts, _k), coeffs(_coeffs)
                {

		base_mod = true;
                if( (int) coeffs.size() != get_nbasis())
                        {
			DEBUG(0.0);                        
			std::cerr<< "Error in f_bspline constructor, missmatch of coefficients size with number of basis Bix, aborting ..." <<endl; 
			throw 0;
                        }

}



f_bspline::f_bspline(const vector<double> &_coeffs, bspline_basis mybasis)
                : bspline_basis(mybasis), coeffs(_coeffs)
                {	
		
		base_mod = true;
                if((int) coeffs.size() != get_nbasis())
                        {
			DEBUG(0.0);                        
			std::cerr<< "Error in f_bspline constructor, missmatch of coefficients size with number of basis Bix, aborting ..." <<endl; 
			throw 0;
                        }

}




void f_bspline:: set_base_mod(bool new_base_mod){
	base_mod = new_base_mod; 
}


/*!
Insert elements X into knot vector, thus create a new knot vector. Coefficients coeffs are updated. Old results are discarded. 
*/
void f_bspline::refine_knot(vector<pair<double,int> > &X_wth_mult)
	{ 
	
	// Transorm X_wth_mult --> X :: vector<double>
	vector<double> X;
	int nbreak = X_wth_mult.size();
        int nreserve = nbreak * get_order(); // Reserve as if each knot value had multiplicity equal to order. 
	X.reserve(nreserve);


        for(int i=0; i<nbreak; ++i)
                {
                X.push_back(X_wth_mult[i].first);
                for(int j=1; j < X_wth_mult[i].second; j++)
                        {
                        X.push_back(X_wth_mult[i].first);
                        }
                }



	/*	
	A. Create function that makes the change from knots_wth_mult --> knots, and versa. 
	*/

	vector<double> coeffs_new(coeffs.size() + X.size()); 
	vector<double> knots_new(get_knots().size()+X.size()); 
	sort(X.begin(),X.end());



	int r = X.size()-1;
	int n = get_nbasis()-1;
	int p = get_order()-1;
	int m = n + p + 1;

	int a = find_knot_span_of_x(X.front());
	int b = find_knot_span_of_x(X.back());

	b=b+1;

	for(int j=0; j<= a-p; ++j)
		coeffs_new[j] = coeffs[j];
	for(int j=b-1; j<=n; ++j)
		coeffs_new[j+r+1]=coeffs[j];
	for(int j=0; j<=a; ++j)
		knots_new[j] = get_knots()[j];
	for(int j=b+p; j<=m; ++j)
		knots_new[j+r+1] = get_knots()[j];

	int i=b+p-1;
	int kk=b+p+r;

	for(int j=r; j >= 0; j--)
		{
		while(X[j] <= get_knots()[i] && i>a){
		coeffs_new[kk-p-1] = coeffs[i-p-1];
		knots_new[kk] = get_knots()[i];
		kk = kk-1;
		i=i-1;
		};

		coeffs_new[kk-p-1] = coeffs_new[kk-p];

		for(int l=1; l <= p; l++)
			{
			int index = kk-p+l;
			double alfa = knots_new[kk+l] - X[j];
			if (alfa==0.0){
				coeffs_new[index-1] = coeffs_new[index];
			}else{
				alfa = alfa/ (knots_new[kk+l]-get_knots()[i-p+l]);
				coeffs_new[index-1] = alfa * coeffs_new[index-1] + (1.0-alfa)*coeffs_new[index];
			}
			}
		
		knots_new[kk]=X[j];
		kk--;
		}

	

	set_knots(knots_new);
	set_coeffs(coeffs_new);
	
	base_mod = true;
}







/*!
Insert elements X into knot vector, thus create a new knot vector. Coefficients coeffs are updated. Old results are discarded. 
*/
void f_bspline::refine_knot(vector<double> &X)
	{ 

	/*	
	A. Create function that makes the change from knots_wth_mult --> knots, and versa. 
	*/

	vector<double> coeffs_new(coeffs.size() + X.size());
	vector<double> knots_new(get_knots().size()+X.size()); 
	sort(X.begin(),X.end());



	int r = X.size()-1;
	int n = get_nbasis()-1;
	int p = get_order()-1;
	int m = n + p + 1;

	int a = find_knot_span_of_x(X.front());
	int b = find_knot_span_of_x(X.back());
	b++; // This makes sure that x_i < b for all x_i \in X 

	for(int j=0; j<= a-p; ++j)
		coeffs_new[j] = coeffs[j];
	for(int j=b-1; j<=n; ++j)
		coeffs_new[j+r+1]=coeffs[j];
	for(int j=0; j<=a; ++j)
		knots_new[j] = get_knots()[j];
	for(int j=b+p; j<=m; ++j)
		knots_new[j+r+1] = get_knots()[j];

	int i=b+p-1;
	int kk=b+p+r;

	for(int j=r; j >= 0; j--)
		{
		while(X[j] <= get_knots()[i] && i>a){
		coeffs_new[kk-p-1] = coeffs[i-p-1];
		knots_new[kk] = get_knots()[i];
/* CHECK */	//knots_new_wth_mult[kk] = knots_wth_mult[i];	
		kk--;
		i--;
		};

		coeffs_new[kk-p-1] = coeffs_new[kk-p];

		for(int l=1; l <= p; l++)
			{
			int index = kk-p+l;
			double alfa = knots_new[kk+l] - X[j];
			if (alfa==0.0){
				coeffs_new[index-1] = coeffs_new[index];
			}else{
				alfa = alfa/ (knots_new[kk+l]-get_knots()[i-p+l]);
				coeffs_new[index-1] = alfa * coeffs_new[index-1] + (1.0-alfa)*coeffs_new[index];
			}
			}
		
		knots_new[kk]=X[j];
		kk--;
		}

	set_knots(knots_new);
	set_coeffs(coeffs_new);

	
	base_mod = true;

}



/*
// Tiller 92 Algorithm. 
void f_bspline::remove_knot(const int &s, const int &num) 	// Remove knot t_r, num times (try to remove it, with TOL = 0.0 
	{
	double TOL=0.0; // Allow removal of only exact type - write better. 
	// This function calculates the distance/difference between original and evaluated control points, or I don't know what else it does. . . 
	auto Distance = [] (double x, double y){
		//return 0.0;
		double temp=x-y;
		return 
			temp*temp; 
	};

	vector<double> knots_local = get_knots();
	



// ****************************************
	int mult = 2; // Initial multiplicity of knot index idx. 
	double u = knots_local[s];
	
	int n = get_nbasis()-1;
	int p = get_order()-1;
	int ord = get_order();
	int m = n + p + 1; // Total number of knots - 1;  


	int fout  = (2*s-mult-p)/2; // First control point out 
	int last  = s - mult; 
	int first = s - p;
	
	// loop variables 
	vector<double> temp(2*p+1);
	int i,j,ii,jj,t,off;

	int remflag;

	double alfi, alfj;

	for(t=0; t<num; first--, last++,t++)
		{
		off = first - 1;

		temp[0] = coeffs[off];
		temp[last+1-off] = coeffs[last+1];

		i=first;
		j=last;

		ii=1;
		jj=last-off;		

		remflag=0;
		while (j-i>t)
			{
			alfi = (u - knots_local[i]) / (knots_local[i+ord+t] - knots_local[i]);
			alfj = (u-knots_local[j-t]) / (knots_local[j+ord] - knots_local[j-t]);

			temp[ii] = (coeffs[i++] - (1.-alfi)*temp[ii-1])/alfi;
			temp[jj] = (coeffs[j--] - alfj*temp[jj+1]) / (1.0-alfj);

			ii++;
			jj--;

			} 

		if(j-i<t){
			if(Distance(temp[ii-1],temp[jj+1])  <=TOL)
				remflag = 1;
			}else{
			alfi = (u - knots_local[i]) / (knots_local[i+ord+t]-knots_local[i]);
			if (Distance(coeffs[i], alfi*temp[ii+t+1]+(1.-alfi)*temp[ii-1] ) <= TOL )
				remflag = 1;

			}
	
		if (remflag == 0){
				break;	
			}else{ // Succesfull removal, save new control points

			i=first;
			j=last;
			while (j-i>t)
				{
				coeffs[i] = temp[i-off];
				coeffs[j] = temp[j-off];
				i++;
				j--;
				}
			}
		} // end of for (t=0 ... ) loop. 

		if (t==0) // If no knot removal occured, get out of here 
			return ;  
		int kk;
		for(kk=s+1; kk<m; kk++)
			knots_local[kk-t] = knots_local[kk];
		
		j=fout;
		i=j;
		for(kk=1; kk<t; kk++)
			if (kk%2 == 1){
				i++;
			}else{
				j--;
			}  

		for(kk=i+1; kk <= n; j++, kk++)
			coeffs[j] = coeffs[kk];


// ************************************************
	for(int pop=0; pop<t; pop++)
		{
		coeffs.pop_back();
		knots_local.pop_back();
		}

	set_knots(knots_local);
	set_coeffs(coeffs);

cout<< "coefficients AFTER removal: " << endl; 
for(auto it: coeffs )
	cout<< it << " ";
cout<<endl; 

cout<< "Knots AFTER removal: " << endl; 
for(auto it: knots_local )
	cout<< it << " ";
cout<<endl; 
exit(1);

}
*/






/*!
Remove knot t_r != t_{r+1} that appears with multiplicity s in the knot vector. 
*/
void f_bspline::remove_knot(const int &r, const int &num) 	//! Remove knot t_r, num times (try to remove it, with TOL = 0.0 
	{

	// Tollerance factor, 
	double TOL=1.e-10; // Allow removal of only exact type - write better. 
	// This function calculates the distance/difference between original and evaluated control points, or I don't know what else it does. . . 
	auto Distance = [] (double x, double y){
		//return 0.0;
		double temp=x-y;
		return 
			abs(temp); 
	};

	vector<double> knots_local = get_knots();
	vector<double> breakpts_local = get_breakpts();
	// local copy of evaluation variable. 
	double u = knots_local[r];
	// find corresponding index of knot_wth_mult to get the multiplicity. 
	auto it = find(breakpts_local.begin(),breakpts_local.end(),u);	
	int index = distance (breakpts_local.begin(),it);	
	if (it == breakpts_local.end())
		{
		for(auto &it: breakpts_local )
			cout<< it << " ";
		cout << endl; 
		DEBUG(u);	
		}

	int  s = get_knots_wth_mult()[index].second; 

	int n = get_nbasis()-1;
	int p = get_order()-1;
	int ord = get_order();
	int m = n + p + 1;

	
	vector<double> temp(2*p+1,0.0);

	// Is this a correct result? 
	int fout  = (2*r-s-p)/2;
	int last  = r-s;
	int first = r-p;

	// loop variables 
	int i, ii, j, jj, remflag, t, off;
	double alfi, alfj;


	for(t=0; t<num; t++)
		{
		off = first - 1; // difference in index between temp and coeffs
		temp[0] = coeffs[off];
		temp[last+1-off] = coeffs[last+1];
		i=first;
		j=last;
		ii=1;
		jj=last-off;
		remflag=0;
		while(j-i > t){ // Compute new control points for one removal step. 
		
			alfi = (u - knots_local[i]) / ( knots_local[i+ord+t] -knots_local[i]);
			alfj = (u - knots_local[j-t]) / (knots_local[j+ord]-knots_local[j-t]);

			temp[ii] = (coeffs[i] - (1.0-alfi)*temp[ii-1])/alfi;
			temp[jj] = (coeffs[j] - alfj*temp[jj+1]) / (1. - alfj);


			i++;
			ii++;
			j--;
			jj--;
			}

		if(j-i <t){
			if(Distance(temp[ii-1],temp[jj+1])<= TOL ){
				remflag=1;
				}

		}else { 
			alfi = (u - knots_local[i])/ (knots_local[i+ord+t] - knots_local[i]);
			if ( Distance(coeffs[i], alfi*temp[ii+t+1] + (1.-alfi)*temp[ii-1] ) <= TOL ){
					remflag = 1;
					}

			}
		if(remflag==0){
			break;
		}else{

		i=first;
		j=last;
		while( j-i > t)
			{
			coeffs[i] = temp[i-off];
			coeffs[j] = temp[j-off];
			i++;
			j--;
			}
		}
		first--;
		last++;
		}
		// This tells the algorithm to exit, if no removal took place 
		//{

		if(t==0){
			base_mod=true;
			return;
		}


		for(int kk=r+1; kk <=m; kk++)
			knots_local[kk-t] = knots_local[kk];	

		j=fout;
		i=j;
		
		for(int kk=1; kk<t; kk++)
			if(kk%2 ==1){
			i++;
			}else{
			j--;
			}
		
		for(int kk=i+1; kk<=n; kk++)
			{
			coeffs[j] = coeffs[kk];
			j++;
			}

		//return;
	for(int pop=0; pop<t; pop++)
		{
		coeffs.pop_back();
		knots_local.pop_back();
		}

	set_knots(knots_local);
	set_coeffs(coeffs);

	base_mod=true;
}


/*!
Decomposes an f_bspine function to it's bezier segments. This is important when one wishes to perform operations between different functions, e.g. multiply f(x) * g(x). 

*/
void f_bspline::f_Bez_decomp(int &nbezier,vector<vector<double> > &Qbez)	/*! Decomposes function to the number of Bezier segments nb, and coefficients Qbez */	
	{

	// Decomposing into bezier segments: Each Bezier segment will have the same degree as the bspline. The number of Bezier segments are the number of intervals between breakpoints. 
	int dim1 = get_breakpts().size()-1;
	int dim2 = get_order();
	vector<vector<double> > Qbez_local(dim1, vector<double>(dim2,0.0)); 
	

	int n=get_nbasis()-1;
	int p=get_order()-1;

	int m=n+p+1;

	int a=p;
	int b=p+1;

	vector<double> knots_local = get_knots();

	// loop variables 
	int i,j,nb,mult,r,save,s,kk;
	double numer,alpha; 

	// Knot sure for the total dimentions of alpha yet. 
	vector<double> alphas(p); 

	nb=0;

	for( i=0; i<=p ; i++)
		Qbez_local[nb][i] = coeffs[i];

	while ( b < m)
		{
		i = b;
		while (b < m && knots_local[b+1] == knots_local[b])
			{ 
			b++;
			}
		// Finds multiplicitely of knot (I think!).
		mult = b-i+1;

		if (mult < p )
			{
			numer = knots_local[b] - knots_local[a];
			// Compute and store alphas: 
			for (j=p; j > mult; j--)
				alphas[j-mult-1] = numer / (knots_local[a+j] - knots_local[a]);

			r = p - mult;  // Insert knot r  times; 

			for (j=1; j <=r; j++)
				{
				save = r - j;
				s = mult + j;
				for (kk=p; kk>=s; kk--)
					{
					alpha = alphas [kk-s];
					Qbez_local[nb][kk] = alpha * Qbez_local[nb][kk]+ (1.0- alpha) * Qbez_local[nb][kk-1];
					}
				if (b < m) 
					Qbez_local[nb+1][save] =  Qbez_local[nb][p]; // next segment 
				}
			}


		nb++;
		if (b < m)
			{ // initialize the next segment 
			for (i=p-mult; i<=p; i++)
				Qbez_local[nb][i] = coeffs[b-p+i];
			a=b;
			b++;
			}
		} // END OF while 
	nbezier = std::move(nb); 
	Qbez = std::move(Qbez_local);

}

void f_bspline::set_coeffs(vector<double> &_coeffs) // Sets a new vector of coeffs, without needing to re-define knot vector etc. 
        {
        coeffs=_coeffs;
        if( (int) coeffs.size() != get_nbasis()){
        cout<< "Error in routine f_bspline, ncoeffs is not consistent with Bspline basis definition points ! "<<endl;
	cout<< "coeffs.size(): " << coeffs.size() << ", nbasis: " << get_nbasis() << endl; 
	DEBUG(0.0);                        
        throw 0;
        }
}



vector<double> 	f_bspline::get_coeffs(){ 				/*! Returns coefficients of f_bspline */ 
	return 
		coeffs;


}
double f_bspline::f_eval(const double &x )
        {
        if((int) coeffs.size() != get_nbasis()){
	DEBUG(coeffs.size());                        
        cout<< "Error in routine f_bspline, ncoeffs is not consistent with Bspline basis definition points ! "<<endl;
        throw 0;
        }

	pair<int,int> temp = find_nonzero_basis_at_x(x); 

        double ttempsum(0.0);
        for(int i = temp.first; i <= temp.second; i++)
                ttempsum += coeffs[i]*get_Bix(i,x);
        return  ttempsum;
}

// Gives derivative of order j (e.g. j=1 first derivative Df/dx). 
double f_bspline::Djf_eval(const int &j, const double &x )
        {
        if((int) coeffs.size() != get_nbasis()){
	DEBUG(coeffs.size());                        
        cout<< "Error in routine Djf_eval, ncoeffs is not consistent with number of Bspline basis functions, aborting ...  "<<endl;
        throw 0;
        }

	pair<int,int> temp = find_nonzero_basis_at_x(x); 

        double tempsum(0.0);
        for(int i = temp.first; i <= temp.second; i++)
                tempsum += coeffs[i]*get_DjBix(j,i,x);
        return  tempsum;


}


bspline_basis f_bspline::get_basis_int (const int &j){
	
	
	// If integrals are all ready evaluated then this function does nothing :)
	eval_integrals(base_mod);

	return
		basis_int[j];

}



vector<double>  f_bspline::get_coeffs_int (const int &j){
	

	
	// If integrals are all ready evaluated then this function does nothing :)
	eval_integrals(base_mod);

	return 
		coeffs_int[j];


}



double f_bspline::Intf_eval(const int &j,  const double &x ) // Returns \int_0^x ...\int_0^x f(x)dx , j integrations in total. 
        {
        if(j<1 || j > Nintegrations ){
                cout<< "Request for integral evaluation j:= "<< j << " outside of pre-specified, [1,"<<Nintegrations<<"]."<<endl;
                cout<< "If you want more integrations, change variable Nintegrations in class f_bspline. Aborting ..." << endl;
                throw(0);
        }
	
        if((int) coeffs.size() != get_nbasis()){
	DEBUG(coeffs.size());                        
        cout<< "Error in routine Djf_eval, ncoeffs is not consistent with number of Bspline basis functions, aborting ...  "<<endl;
        throw 0;
        }

	// If integrals are all ready evaluated then this function does nothing :)
	eval_integrals(base_mod);

	pair<int,int> temp = basis_int[j-1].find_nonzero_basis_at_x(x); 

        double tempsum(0.0);
        for(int i = temp.first; i <= temp.second; i++)
                tempsum += coeffs_int[j-1][i]*basis_int[j-1].get_Bix(i,x);
        return  tempsum;


}











#endif
