#pragma once
#include "stdafx.h"

Polynomial::Polynomial(double* coef, int order) {
		variable = 'x';
		int largestExp = 0;
		coefficient = new double[order];
		for (int i = 0; i < order; i++) {
			coefficient[i] = coef[i];
			if (coef[i] != 0) {
				largestExp = i;
			}
		}
		terms = largestExp + 1;
	}

Polynomial::Polynomial(long double* coef, int order) {
		variable = 'x';
		int largestExp = 0;
		coefficient = new double[order];
		for (int i = 0; i < order; i++) {
			coefficient[i] = coef[i];
			if (coef[i] != 0) {
				largestExp = i;
			}
		}
		terms = largestExp + 1;
	}

Polynomial::Polynomial(std::vector<double> coef) {
		variable = 'x';
		int largestExp = 0;
		coefficient = new double[coef.size()];
		for (int i = 0; i < coef.size(); i++) {
			coefficient[i] = coef[i];
			if (coef[i] != 0) {
				largestExp = i;
			}
		}
		terms = largestExp + 1;
	}

Polynomial::Polynomial(std::vector<double> coef, char x) {
		variable = x;
		int largestExp = 0;
		coefficient = new double[coef.size()];
		for (int i = 0; i < coef.size(); i++) {
			coefficient[i] = coef[i];
			if (coef[i] != 0) {
				largestExp = i;
			}
		}
		terms = largestExp + 1;
	}

Polynomial::Polynomial(double* coef, int order, char var) {
		variable = var;
		int largestExp = 0;
		coefficient = new double[order];
		for (int i = 0; i < order; i++) {
			coefficient[i] = coef[i];
			if (coef[i] != 0) {
				largestExp = i;
			}
		}
		terms = largestExp + 1;
	}

Polynomial::Polynomial() {}
Polynomial::~Polynomial() {}

	void Polynomial::randomize() {
		srand(time(0) + clock());
		coefficient = new double[terms];
		for (int i = 0; i < terms; i++) {
			double n = ((double)rand() / 10000) * pow(-1, rand() % 2);
			coefficient[i] = n;
		}
	}

	double Polynomial::largestCoefficient() {
		double biggest = 0;
		for (int i = 0; i < terms; ++i) {
			if (fabsf(coefficient[i]) > biggest) { biggest = abs(coefficient[i]); }
		}
		return biggest;
	}

	double Polynomial::evaluate(double n) {//slower, more exact method of evaluation
		double answer = 0;
		for (int i = 0; i < terms; i++) {
			if (coefficient[i] != 0) { answer += (coefficient[i] * pow(n, i)); }
		}
		return answer;
	}

	ComplexNumber Polynomial::evaluate(ComplexNumber z) {
		ComplexNumber ans(0, 0);
		for (int i = 0; i < terms; i++) {

			ans = ans.add(ans, ans.multiply(ans.exponent(z, i), coefficient[i]));
		}
		return ans;
	}

	double Polynomial::HornerEvaluate(double x) {//fastest evaluation
		double result = 0.0;
		for (int idx = terms - 1; idx >= 0; idx--) {
			result = fma(result, x, coefficient[idx]);
		}
		return result;
	}

	ComplexNumber Polynomial::HornerEvaluateComplex(ComplexNumber z) {
		std::vector<double> v = toSTLVector(coefficient, terms);
		ComplexNumber s(0, 0);
		for (std::vector<double>::const_reverse_iterator i = v.rbegin(); i != v.rend(); i++) {
			s = z.add(z.multiply(s, z), *i);
		}
		return s;
	}

	double Polynomial::NewtonsMethod(double x) {//given some root, this method will refine it until it is much more accurate
		double x1;
		double f, dfdx;
		const double epsilon = 0.0000001; //tolerance
		const int n = 15; //# steps
		Polynomial deriv = derivative();//derivative of polynomial

		x1 = x;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = HornerEvaluate(x); //Defines the function 
			dfdx = deriv.HornerEvaluate(x); //Derivative of the function
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (abs(x - x1) < epsilon) { //checks if guesses are within a certain tolerance value 
			return x1;
		}
		return 0;
	}

	double Polynomial::evaluateY(double y) {	//solve for the value of x which yields a particular f(x) -- sort of like the inverse function of x, f^-1(x)
		for (double i = 0; i < 600; i += 0.001) {//IN PROGRESS
			double answer = 0;
			for (int j = 0; j < terms; j++) {
				answer += (coefficient[j] * pow(i, j));
			}
			if (answer == y) {
				return i;
			}
		}
	}

	Polynomial Polynomial::add(Polynomial a, Polynomial b) {
		int new_length = a.terms;
		int shorter_length = b.terms;
		if (b.terms > a.terms) { new_length = b.terms; shorter_length = a.terms; }
		double* coef1 = new double[new_length];
		for (int i = 0; i < new_length; i++) {
			if (i < shorter_length) { coef1[i] = a.coefficient[i] + b.coefficient[i]; }
			else {
				if (b.terms > a.terms) { coef1[i] = b.coefficient[i]; }
				else { coef1[i] = a.coefficient[i]; }
			}
		}
		return Polynomial(coef1, new_length);
	}

	Polynomial Polynomial::add(Polynomial a, double  b) {
		double* coef1 = coefficient;
		coef1[0] = a.coefficient[0] + b;
		return Polynomial(coef1, terms);
	}

	Polynomial Polynomial::add(double b, Polynomial a) {
		double* coef1 = coefficient;
		coef1[0] = a.coefficient[0] + b;
		return Polynomial(coef1, terms);
	}

	Polynomial Polynomial::subtract(Polynomial a, Polynomial b) {
		int new_length = a.terms;
		int shorter_length = b.terms;
		if (b.terms > a.terms) { new_length = b.terms; shorter_length = a.terms; }
		double* coef1 = new double[new_length];
		for (int i = 0; i < new_length; i++) {
			if (i < shorter_length) { coef1[i] = a.coefficient[i] - b.coefficient[i]; }
			else {
				if (b.terms > a.terms) { coef1[i] = b.coefficient[i]; }
				else { coef1[i] = a.coefficient[i]; }
			}
		}
		return Polynomial(coef1, new_length);
	}

	Polynomial Polynomial::subtract(Polynomial a, double  b) {
		double* coef1 = coefficient;
		coef1[0] = a.coefficient[0] - b;
		return Polynomial(coef1, terms);
	}

	Polynomial Polynomial::subtract(double b, Polynomial a) {
		double* coef1 = coefficient;
		coef1[0] = a.coefficient[0] - b;
		return Polynomial(coef1, terms);
	}

	Polynomial Polynomial::multiply(Polynomial a, Polynomial b) {
		int new_length = a.terms + b.terms;
		double* coef1 = new double[new_length];
		for (int i = 0; i < new_length; i++) {
			coef1[i] = 0;
		}
		for (int i = 0; i < a.terms; i++) {
			for (int j = 0; j < b.terms; j++) {
				int index = i + j;
				coef1[index] += a.coefficient[i] * b.coefficient[j];
			}
		}
		return Polynomial(coef1, new_length);
	}
	Polynomial Polynomial::multiply(Polynomial a, double b) {
		double* coef1 = new double[a.terms];
		for (int i = 0; i < a.terms; i++) {
			coef1[i] = a.coefficient[i] * b;
		}
		return Polynomial(coef1, a.terms);
	}
	Polynomial Polynomial::multiply(double b, Polynomial a) {
		double* coef1 = new double[a.terms];
		for (int i = 0; i < a.terms; i++) {
			coef1[i] = a.coefficient[i] * b;
		}
		return Polynomial(coef1, a.terms);
	}

	Polynomial Polynomial::divide(Polynomial p, Polynomial q) {//factors out a binomial by synthetic division
		double n = q.coefficient[0];
		if (q.terms > 2) { return Polynomial(); }
		double* vec = new double[p.terms - 1];
		vec[terms - 2] = p.coefficient[terms - 1];
		for (int i = p.terms - 2; i > 0; --i) {
			vec[i - 1] = p.coefficient[i] - (vec[i] * n);
		}
		if (p.coefficient[0] - (vec[0] * n) != 0) { return Polynomial(); }
		return Polynomial(vec, terms - 1);
	}

	Polynomial Polynomial::divide(Polynomial p2, double q) {
		Polynomial p = p2;
		for (int i = 0; i < p.terms; ++i) { p.coefficient[i] *= (1 / q); }
		return p;
	}

	Polynomial Polynomial::exponentiate(Polynomial p, int b) {
		if (b > 1) {
			for (int i = 1; i < b; ++i) {
				p = p.multiply(p, *this);
			}
		}
		return p;
	}

	//Overloaded operators:
	Polynomial& Polynomial::operator*=(Polynomial& rhs) { *this = multiply(*this, rhs); return *this; }
	Polynomial& Polynomial::operator*=(double x) { *this = multiply(*this, x); return *this; }
	Polynomial& Polynomial::operator/=(Polynomial& rhs) { *this = divide(*this, rhs); return *this; }
	Polynomial& Polynomial::operator/=(double x) { return divide(*this, x); return *this; }
	Polynomial& Polynomial::operator+=(Polynomial& rhs) { *this = add(*this, rhs); return *this; }
	Polynomial& Polynomial::operator+=(double x) { *this = add(*this, x); return *this; }
	Polynomial& Polynomial::operator-=(Polynomial& rhs) { *this = subtract(*this, rhs); return *this; }
	Polynomial& Polynomial::operator-=(double x) { *this = subtract(*this, x); return *this; }
	Polynomial& Polynomial::operator*(Polynomial& rhs) { return multiply(*this, rhs); }
	Polynomial& Polynomial::operator*(double x) { Polynomial p = multiply(*this, x); return p; }
	Polynomial& Polynomial::operator/(Polynomial& rhs) { return divide(*this, rhs); }
	Polynomial& Polynomial::operator/(double x) { Polynomial p = divide(*this, x); return p; }
	Polynomial& Polynomial::operator+(Polynomial& rhs) { return add(*this, rhs); }
	Polynomial& Polynomial::operator+(double x) { return add(*this, x); }
	Polynomial& Polynomial::operator-(Polynomial& rhs) { return subtract(*this, rhs); }
	Polynomial& Polynomial::operator-(double x) { return subtract(*this, x); }
	std::wostream& operator<<(std::wostream& os, Polynomial& rhs) {
		os << rhs.toString();
		return os;
	}

	bool Polynomial::isIntegerValued() {
		for (int i = 0; i < terms; ++i) {
			if (coefficient[i] != floor(coefficient[i])) { return false; }
		}
		return true;
	}

	bool Polynomial::isMonic() {
		if (terms>0 && coefficient[terms - 1] != 1) { return false; }
		return true;
	}

	Polynomial Polynomial::makeMonic() {
		if (coefficient[terms - 1] != 1) {
			double a = 1 / coefficient[terms - 1];
			for (int i = 0; i < terms; ++i) { coefficient[i] *= a; }
		}
		return *this;
	}

	Polynomial Polynomial::factorOutBinomial(Polynomial p) {//factors out a binomial by synthetic division where binomial p = (x + ...)
		double n = p.coefficient[0] - .1; //get -n from the input binomial (x + n)
		if (p.terms > 2) { return Polynomial(); }//can't factor out larger polynomials


												  //synthetic division w/ adjustment for imprecise n
												  //================================================
		double* vec = NULL;
		int counter = 0;
		double remainder = 99999;
		double delta = 0.001;
		int sign = 1;
		while (std::abs(remainder) > 0.001 && counter < 500) {
			vec = new double[terms - 1];
			vec[terms - 2] = coefficient[terms - 1];//drop down first coefficient of synth div
			for (int i = terms - 2; i >= 0; --i) {
				if (i == 0) { remainder = coefficient[i] + (vec[i] * n); }
				else { vec[i - 1] = coefficient[i] + (vec[i] * n); }
			}
			++counter;
			if (counter <= 1) { sign = sgn(remainder); }
			if (counter > 1) { if (sign != sgn(remainder)) { n -= delta; delta /= 10; } }
			if (std::abs(remainder) < 0.00001) {
				//Polynomial test = Polynomial(vec, length - 1).factorOutBinomial(n);
				//if (test.length != 0) { return test; }//check to see if this is a repeated root
				return Polynomial(vec, terms - 1);
			}
			if (std::abs(remainder) >= 0.00001) { n += delta; }
		}
		return Polynomial(vec, terms - 1);
		//return Polynomial();
	}

	Polynomial Polynomial::factorOutBinomial(double n) {//factors out a binomial by synthetic division
											//if (n == 0) { return Polynomial(); }
		std::vector<double> cof;
		cof.push_back(n);
		cof.push_back(1);
		Polynomial p(cof);
		return factorOutBinomial(cof);
	}

	double Polynomial::syntheticDivision(Polynomial p) {//factors out a binomial by synthetic division where binomial p = (x + ...) and returns the remainder
		double n = p.coefficient[0]; //get -n from the input binomial (x + n)
		if (p.terms > 2) { return 0; }//can't factor out larger polynomials


									   //synthetic division w/ adjustment for imprecise n
									   //================================================
		double* vec = NULL;
		int counter = 0;
		double remainder;
		double delta = 0.001;
		int sign = 1;
		vec = new double[terms - 1];
		vec[terms - 2] = coefficient[terms - 1];//drop down first coefficient of synth div
		for (int i = terms - 2; i >= 0; --i) {
			if (i == 0) { remainder = coefficient[i] + (vec[i] * n); }
			else { vec[i - 1] = coefficient[i] + (vec[i] * n); }
		}
		return remainder;
	}

	double Polynomial::syntheticDivision(double n) {//factors out a binomial by synthetic division, return remainder
		std::vector<double> cof;
		cof.push_back(n);
		cof.push_back(1);
		Polynomial p(cof);
		return syntheticDivision(cof);
	}

	ComplexNumber Polynomial::complexSyntheticDivision(ComplexNumber z) {//factors out a binomial by synthetic division where binomial p = (x + ...) and returns the remainder
		ComplexNumber* vec = NULL;
		int counter = 0;
		ComplexNumber remainder;

		vec = new ComplexNumber[terms - 1];
		vec[terms - 2] = ComplexNumber(coefficient[terms - 1], 0);//drop down first coefficient of synth div
		vec[terms - 2].display();//debug
		for (int i = terms - 2; i >= 0; --i) {
			if (i == 0) { remainder = z.add(ComplexNumber(coefficient[i], 0), ComplexNumber(vec[i].multiply(z, vec[i]))); }
			else {
				vec[i - 1] = z.add(coefficient[i], vec[i].multiply(vec[i], z));
				vec[i - 1].display();//debug
			}
		}
		remainder.display();//debug
		std::cin.ignore();//debug
		vec = NULL;
		return remainder;
	}

	bool Polynomial::checkLowBound(Polynomial p) {//factors out a binomial by synthetic division where binomial p = (x + ...) and returns the true is signs alternate
		double n = p.coefficient[0]; //get n from the input binomial (x + n)
		if (p.terms > 2) { return 0; }//can't factor out larger polynomials

									   //synthetic division
									   //==================
		double* vec = NULL;
		int counter = 0;
		double remainder;
		double delta = 0.001;
		int signflips = 0;
		int sign = 0;

		vec = new double[terms - 1];
		vec[terms - 2] = coefficient[terms - 1];//drop down first coefficient of synth div
		sign = sgn(vec[terms - 2]);
		//std::wcout << vec[terms - 2] << L"\n";//debug
		for (int i = terms - 2; i >= 0; --i) {
			if (i == 0) { remainder = coefficient[i] + (vec[i] * n); }
			else {
				vec[i - 1] = coefficient[i] + (vec[i] * n);
				//std::wcout << vec[i - 1] << L"\n";//debug
				int sn = sgn(vec[i - 1]);
				if (sign != sn) {
					++signflips;
					sign = sn;
				}
			}
		}
		if (sgn(remainder) != sign) { ++sign; }
		//std::wcout << signflips << L" ---> " << terms - 2 << L"\n";//debug
		//std::cin.ignore();//debug
		if (signflips >= terms - 2) { return true; }
		return false;
	}

	bool Polynomial::checkLowBound(double n) {//factors out a binomial by synthetic division, return remainder
		std::vector<double> cof;
		cof.push_back(n);
		cof.push_back(1);
		Polynomial p(cof);
		return checkLowBound(cof);
	}

	Polynomial Polynomial::derivative() {
		if (terms == 0 || terms == 1) { return Polynomial(); }
		double* coef1 = new double[terms - 1];
		for (int i = 0; i < terms - 1; i++) {
			coef1[i] = coefficient[i + 1] * (i + 1);
		}
		Polynomial p(coef1, terms - 1);
		coef1 = NULL;
		return p;
	}

	Polynomial Polynomial::integral() {
		double* coef1 = new double[terms + 1];
		coef1[0] = 0;
		for (int i = 1; i < terms + 1; i++) {
			if (coefficient[i - 1] != 0) { coef1[i] = ((double)1 / i)*coefficient[i - 1]; }
			else { coef1[i] = 0; }
		}
		return Polynomial(coef1, terms + 1);
	}

	Polynomial Polynomial::CauchyProduct(Polynomial A, Polynomial B) {	//NEEDS TO BE CHECKED
		if (A.terms != B.terms) {
			std::wcout << L"ERROR!  Mismatched array lengths.";
			return Polynomial();
		}

		int new_length = A.terms;
		std::vector<double> coef;

		for (int i = 0; i < new_length; i++) {
			coef.push_back(0);
		}

		for (int i = 0; i < new_length; i++) {
			int index = i + (new_length - i);
			coef[index] = A.coefficient[i] * B.coefficient[new_length - i];
		}
		Polynomial answer(coef);
		return answer;
	}

	std::vector<double> Polynomial::getRootBounds() {
		Polynomial temp(coefficient, terms);
		if (temp.isMonic() == false) { temp.makeMonic(); }
		double M = temp.largestCoefficient() + 1;//by Cauchy's theorem
		std::vector<double> answer;
		int lb = 0;
		int hb = 0;

		//Now, test for smaller bounds by using synthetic division:

		//check lower bound
		for (double i = -M; i < M; ++i) {
			if (checkLowBound(i) == true) { lb = i; }
			if (checkLowBound(i) == false) { i = M; }//once you've passed beyond the bound, break
		}
		if (lb != 0) { answer.push_back(lb); }
		else { answer.push_back(-M); }

		//check upper bound
		for (double i = 0; i < M; i++) {
			if (std::abs(syntheticDivision(i)) < 0.01) { hb = i; i = M; }
		}
		if (hb != 0) { answer.push_back(hb); }
		else { answer.push_back(M); }

		return answer;
	}

	double Polynomial::root(Polynomial p, double lower_bound, double upper_bound) {//bounded search for a root
		//given a set of bounds for a polynomial, this method will find the first root
		if (p.terms <= 1) { return 0; }
		if (p.terms == 2) { return (-p.coefficient[0] / p.coefficient[1]); }
		if (p.terms == 3) {
			double part1 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part1 >= 0) { return ((-p.coefficient[1] + sqrt(part1)) / (2 * p.coefficient[2])); }
			else { return 0; }
		}

		//else, p is of order 3 or higher....
		double TOL = 0.0001;    // tolerance
		double MAX_ITER = 1000; // maximum number of iterations
		double a = lower_bound;
		double b = upper_bound;
		double fa = p.HornerEvaluate(a);   // calculated now to save function calls
		double fb = p.HornerEvaluate(b);   // calculated now to save function calls
		double fs = 0;      // initialize 

		if (!(fa * fb < 0)) { return 0; }

		if (fabsf(fa) < fabsf(b)) { // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
			std::swap(a, b);
			std::swap(fa, fb);
		}

		double c = a;           // c now equals the largest magnitude of the lower and upper bounds
		double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
		bool mflag = true;      // boolean flag used to evaluate if statement later on
		double s = 0;           // Our Root that will be returned
		double d = 0;           // Only used if mflag is unset (mflag == false)

		for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
		{
			// stop if converged on root or error is less than tolerance
			if (std::abs(b - a) < TOL) {
				goto NarrowResult;
			} // end if

			if (fa != fc && fb != fc) {
				// use inverse quadratic interopolation
				s = (a * fb * fc / ((fa - fb) * (fa - fc)))
					+ (b * fa * fc / ((fb - fa) * (fb - fc)))
					+ (c * fa * fb / ((fc - fa) * (fc - fb)));
			}
			else { s = b - fb * (b - a) / (fb - fa); }//secant method

		  /*	Crazy condition statement!:
		  -------------------------------------------------------
		  (condition 1) s is not between  (3a+b)/4  and b or
		  (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
		  (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
		  (condition 4) (mflag is set and |b−c| < |TOL|) or
		  (condition 5) (mflag is false and |c−d| < |TOL|) */
			if (((s < (3 * a + b) * 0.25) || (s > b)) ||
				(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
				(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
				(mflag && (std::abs(b - c) < TOL)) ||
				(!mflag && (std::abs(c - d) < TOL)))
			{
				// bisection method
				s = (a + b)*0.5;
				mflag = true;
			}
			else { mflag = false; }

			fs = p.HornerEvaluate(s);  // calculate fs
			d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
			c = b;      // set c equal to upper bound
			fc = fb;    // set f(c) = f(b)

			if (fa * fs < 0) {// fa and fs have opposite signs
				b = s;
				fb = fs;    // set f(b) = f(s)
			}
			else {
				a = s;
				fa = fs;    // set f(a) = f(s)
			}

			if (std::abs(fa) < std::abs(fb)) {// if magnitude of fa is less than magnitude of fb
				std::swap(a, b);     // swap a and b
				std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
			}
		}
		return 0;

	NarrowResult:
		//use root as guess for Newton's method now
		double x, x1; //user input value
		double f, dfdx; //function and its derivative 
		const double epsilon = 0.0000001; //tolerance
		const int n = 15; //# steps
		Polynomial deriv = p.derivative();//derivative of polynomial

		x1 = s;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.HornerEvaluate(x); //Defines the function 
			dfdx = deriv.HornerEvaluate(x); //Derivative of the function
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (abs(x - x1) < epsilon) { //checks if guesses are within a certain tolerance value 
			return x1;
		}

		return 0;
	}

	std::vector<double> Polynomial::realRoots(Polynomial p) {
		std::vector<double> candidates;
		if (p.terms <= 1) { return candidates; }
		if (p.terms == 2) {
			candidates.push_back(-p.coefficient[0] / p.coefficient[1]);
			return candidates;
		}
		if (p.terms == 3) {
			double part2 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part2 >= 0) {
				double part1 = (-p.coefficient[1]) / (2 * p.coefficient[2]);
				part2 = sqrt(part2) / (2 * p.coefficient[2]);
				candidates.push_back(part1 + part2);
				candidates.push_back(part1 - part2);
			}
			return candidates;
		}
		//if the polynomial is longer then 2nd order...
		std::vector<double> bounds = p.getRootBounds();
		double lowbound = bounds[0];
		double highbound = bounds[1];

		//Since 0 is used as a null value for roots in the rootfinding algorithms,
		//we must determine whether or not it is a root outside of these algorithms.
		if (p.HornerEvaluate(0) == 0) { candidates.push_back(0); }

		//now sweep through the bounds of the polynomial to look for roots using the
		//mixed Brent/Newton search method in findRoot()
		Polynomial temp = p;
		double lastroot = 0;
		double delta = 0.01;
		if (getLength(p.largestCoefficient()) > 3) { delta *= 10 * getLength(p.largestCoefficient()); }
		while (lowbound < highbound) {
			double x = root(temp, lowbound, lowbound + delta);
			if (x == 0) { lowbound += delta; }//using 0 as a null return value here
			else {
				if (x != lastroot) { //prevent same root from continually showing up
					candidates.push_back(x);
					lastroot = x;
					lowbound = x + delta;
					if (x == floor(x)) {
						temp = temp.factorOutBinomial(x); //perform a shift to reduce the polynomial
						bounds = temp.getRootBounds();
						double lowbound = bounds[0];
						double highbound = bounds[1];
					}//if x is an integer, you can factor it out of the polynomial
				}
				else { lowbound += delta; }
			}
		}
		return candidates;
	}

	std::vector<double> Polynomial::realRoots() {
		Polynomial p(coefficient, terms - 1);
		return realRoots(p);
	}

	std::vector<double> Polynomial::getComplexRootBounds() {
		Polynomial temp(coefficient, terms);
		if (temp.isMonic() == false) { temp.makeMonic(); }
		double M = temp.largestCoefficient() + 1;//by Cauchy's theorem
		std::vector<double> answer;
		double Hbound = 0.5 * ((std::abs(coefficient[terms - 2]) - 1) + sqrt(pow(coefficient[terms - 2] - 1, 2) + (4 * M)));//Sun and Hsieh bound
		answer.push_back(-Hbound);
		answer.push_back(Hbound);

		std::vector<double> co;
		co.push_back(1);
		co.push_back(2 - std::abs(coefficient[terms - 2]));
		co.push_back(1 - std::abs(coefficient[terms - 2]) - std::abs(coefficient[terms - 3]));
		co.push_back(-M);
		Polynomial P(co);
		co = realRoots(P);
		//displayarr(co);//debug
		//std::cin.ignore();//debug
		if (co.size() > 1) {
			answer[0] = co[0];
			answer[1] = co[1];
		}
		return answer;
	}

	std::wstring Polynomial::toString() {
		std::wstringstream s;

		int lastnonzero = 0;
		for (int i = terms - 1; i >= 0; i--) { if (coefficient[i] != 0) { lastnonzero = i; } }

		for (int i = terms - 1; i >= 0; i--) {
			if (coefficient[i] != 0) {
				if (coefficient[i] == floor(coefficient[i])) { polynomialPrecision = 0; }
				if (coefficient[i] != floor(coefficient[i])) { polynomialPrecision = precision; }

				if (i == 0) {
					if (coefficient[i] >= 0) {
						s << L" + " << std::setprecision(polynomialPrecision) << coefficient[i];
					}
					if (coefficient[i] < 0) {
						std::wstring str = L" - ";
						if (i == terms - 1) { str = L"-"; }
						s << str << std::setprecision(polynomialPrecision) << std::abs(coefficient[i]);
					}
				}
				if (i == 1) {
					if (std::abs(coefficient[i]) != 1) {
						if (coefficient[i] >= 0) {
							s << L" + " << std::setprecision(polynomialPrecision) << coefficient[i] << variable;
						}
						if (coefficient[i] < 0) {
							std::wstring str = L" - ";
							if (i == terms - 1) { str = L"-"; }
							s << str << std::setprecision(polynomialPrecision) << std::abs(coefficient[i]) << variable;
						}
					}
					if (coefficient[i] == 1) {
						s << L" + " << variable;
					}
					if (coefficient[i] == -1) {
						std::wstring str = L" - ";
						if (i == terms - 1) { str = L" - "; }
						s << str << variable;
					}
				}
				if (i > 1) {
					if (coefficient[i] != 1 && coefficient[i] != -1) {

						if (coefficient[i] >= 0) {
							s << L" + " << std::setprecision(polynomialPrecision) << coefficient[i] << variable << L"^" << i;
						}
						if (coefficient[i] < 0) {
							std::wstring str = L" - ";
							if (i == terms - 1) { str = L"-"; }
							s << str << std::setprecision(polynomialPrecision) << std::abs(coefficient[i]) << variable << L"^" << i;
						}

					}
					if (coefficient[i] == 1) {
						s << L" + " << variable << L"^" << i;
					}
					if (coefficient[i] == -1) {
						std::wstring str = L" - ";
						if (i == terms - 1) { str = L"-"; }
						s << str << variable << L"^" << i;
					}
				}
			}
		}
		std::wstring temp = s.str();
		while (isNumber(temp[0]) == false && temp[0] != variable && temp[0] != '-') { temp = temp.substr(1); }
		return temp;
	}

	void Polynomial::printToFile(std::wstring filename) {
		std::wofstream outf(filename, std::ios_base::app);
		outf << toString() << L"\n";
		outf.close();
	}

	void Polynomial::display() { std::wcout << toString() << std::endl; }
//END Polynomial CLASS=========================================================

	Polynomial BellPolynomial(int n) {
		if (n < 0) { return Polynomial(); }
		std::vector<double> B0;//B0 = 1
		B0.push_back(1);
		if (n == 0) { return Polynomial(B0); }
		std::vector<double> B1;//B1 = x
		B1.push_back(0);
		B1.push_back(1);
		if (n == 1) { return Polynomial(B1); }
		std::vector<double> B2;//B2 = x^2 + x
		B2.push_back(0);
		B2.push_back(1);
		B2.push_back(1);
		if (n == 2) { return Polynomial(B2); }

		//if n>2, determine polynomial by recursive method:
		Polynomial Bnext;
		Polynomial Bprev = B2;

		std::vector<double> temp;
		temp.push_back(0);
		temp.push_back(1);
		Polynomial mult(temp);//multiplicand is p(x) = x

		for (int i = 3; i <= n; ++i) {
			// loop:  Bn+1 = x [Bn + Bn']
			Bnext = Bprev.multiply(mult, Bprev.add(Bprev, Bprev.derivative()));
			Bprev = Bnext;
		}
		return Bnext;
	}

	Polynomial BernoulliPolynomial(int n) {
		std::vector<double> vec;
		for (int i = 0; i <= n; ++i) {
			double p1 = combination(n, i)*BernoulliNumber(-(n - i));//coefficients are not exactly precise, but numerically "close enough"
			vec.push_back(p1);
		}
		return Polynomial(vec);
	}

	Polynomial BesselPolynomial(int n) {
		std::vector<double> vec;
		for (int i = 0; i <= n; ++i) {
			double p1 = factorial(n + i) / (factorial(n - i)*factorial(i)*pow(2, i));
			vec.push_back(p1);
		}
		return Polynomial(vec);
	}

	int complementaryBellNumber(int n) {
		return BellPolynomial(n).evaluate(-1);
	}

	/*Polynomial Matrix::characteristicPolynomial() {
		std::vector<double> ans = characteristicMatrix().row(0);
		std::vector<double> answer = initZeros(answer, columns + 1);
		for (int i = 0; i < columns; ++i) { answer[columns - 1 - i] = -ans[i]; }
		answer[columns] = 1;
		return Polynomial(answer, columns + 1);
	}*/

	Polynomial ChebyshevPolynomial1stKind(int n) {
		if (n < 0) { return Polynomial(); }
		std::vector<double> T0;
		T0.push_back(1);
		if (n == 0) { return Polynomial(T0); }
		std::vector<double> T1;
		T1.push_back(0);
		T1.push_back(1);
		if (n == 1) { return Polynomial(T1); }
		std::vector<double> T2;
		T2.push_back(-1);
		T2.push_back(0);
		T2.push_back(2);
		if (n == 2) { return Polynomial(T2); }

		//if n>2, determine polynomial by recursive method:
		Polynomial Tnext;
		Polynomial Tprev = T2;
		Polynomial Tprevprev = T1;

		std::vector<double> temp;
		temp.push_back(0);
		temp.push_back(2);
		Polynomial mult(temp);//multiplicand is p(x) = 2x

		for (int i = 3; i <= n; ++i) {
			// loop:  Tn+1 = (2x*Tn) - Tn-1
			Tnext = Tprev.subtract(Tprev.multiply(mult, Tprev), Tprevprev);
			Tprevprev = Tprev;
			Tprev = Tnext;
		}
		return Tnext;
	}

	Polynomial ChebyshevPolynomial2ndKind(int n) {
		if (n < 0) { return Polynomial(); }
		std::vector<double> U0; // U0 = 1
		U0.push_back(1);
		if (n == 0) { return Polynomial(U0); }
		std::vector<double> U1; //U1 = 2x
		U1.push_back(0);
		U1.push_back(2);
		if (n == 1) { return Polynomial(U1); }
		std::vector<double> U2;	//U2 = 4x^2 -1
		U2.push_back(-1);
		U2.push_back(0);
		U2.push_back(4);
		if (n == 2) { return Polynomial(U2); }

		//if n>2, deUermine polynomial by recursive meUhod:
		Polynomial Unext;
		Polynomial Uprev = U2;
		Polynomial Uprevprev = U1;

		std::vector<double> temp;
		temp.push_back(0);
		temp.push_back(2);
		Polynomial mult(temp);//mulUiplicand is p(x) = 2x

		for (int i = 3; i <= n; ++i) {
			// loop:  Un+1 = (2x*Un) - Un-1
			Unext = Uprev.subtract(Uprev.multiply(mult, Uprev), Uprevprev);
			Uprevprev = Uprev;
			Uprev = Unext;
		}
		return Unext;
	}

	Matrix companionMatrix(Polynomial p) {//Froebenius form
		p = p.makeMonic();
		int N = p.terms - 1;
		std::vector<double> elem = initZeros(elem, N*N);
		int count = 0;
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				if (j != i - 1 && j < N - 1) { elem[i*N + j] = 0; }
				if (j == i - 1) { elem[i*N + j] = 1; }
				if (j == N - 1) { elem[i*N + j] = -p.coefficient[count]; count++; }
			}
		}
		return Matrix(N, N, elem);
	}

	Matrix companionMatrix2(Polynomial p) {//alternate form, top row has coefficients
		p = p.makeMonic();
		int N = p.terms - 1;
		std::vector<double> elem = initZeros(elem, N*N);
		for (int i = 0; i < N*N; ++i) { elem[i] = 0; }//initialize as all 0's
		int count = 0;
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				if (j == i - 1) { elem[i*N + j] = 1; }
				if (i == 0) { elem[i*N + j] = -p.coefficient[p.terms - 2 - count]; count++; }
			}
		}
		return Matrix(N, N, elem);
	}

	Quaternion directionQuaternion(Vector v) {
		Vector v2(4);
		for (int i = 0; i < 3; ++i) { v2.vec.push_back(v.vec[i]); }
		v2.vec.push_back(0);
		return v2.toQuaternion();
	}

	Polynomial EulerPolynomial(int m) {
		std::vector<double> temp;
		temp.push_back(1);
		if (m == 0) { return Polynomial(temp); }//E0 = 1
		temp.clear();
		temp.push_back(-0.5);
		temp.push_back(1);
		if (m == 1) { return Polynomial(temp); }//E1 = x-0.5
		temp.clear();

		Polynomial Euler;
		for (int n = 0; n <= m; ++n) {
			double temp1 = pow(2, n);
			temp1 = 1 / temp1;
			Polynomial p;
			for (int k = 0; k <= n; ++k) {
				double p2 = pow(-1, k)*combination(n, k);
				temp.push_back(k);//multiplicand, p(x) = (x+k)
				temp.push_back(1);
				Polynomial ptemp(temp);
				temp.clear();
				ptemp = ptemp.exponentiate(ptemp, m);
				ptemp = ptemp.multiply(ptemp, p2);
				if (k == 0) { p = ptemp; }
				else { p = p.add(p, ptemp); }
			}
			if (n == 0) { Euler = p.multiply(p, temp1); }
			else { Euler = Euler.add(Euler, p.multiply(p, temp1)); }
		}
		return Euler;
	}

	double EulerNumber(int n) {
		if (n % 2 == 1) { return 0; }//if n is odd, En=0;
		return pow(2, n)*EulerPolynomial(n).evaluate(0.5);
	}

	double findRoot(Polynomial p) {
		//given a set of bounds for a polynomial, this method will find the first root
		if (p.terms <= 1) { return 0; }
		if (p.terms == 2) { return (-p.coefficient[0] / p.coefficient[1]); }
		if (p.terms == 3) {
			double part1 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part1 >= 0) { return ((-p.coefficient[1] + sqrt(part1)) / (2 * p.coefficient[2])); }
			else { return 0; }
		}

		//else, p is of order 3 or higher....
		double bound = 1 + (fabsf(p.largestCoefficient()));
		double TOL = 0.0001;    // tolerance
		double MAX_ITER = 1000; // maximum number of iterations
		double a = -bound; //lower_bound;
		double b = bound; //upper_bound;
		double fa = p.HornerEvaluate(a);   // calculated now to save function calls
		double fb = p.HornerEvaluate(b);   // calculated now to save function calls
		double fs = 0;      // initialize 

		if (!(fa * fb < 0)) { return 0; }

		if (fabsf(fa) < fabsf(b)) { // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
			std::swap(a, b);
			std::swap(fa, fb);
		}

		double c = a;           // c now equals the largest magnitude of the lower and upper bounds
		double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
		bool mflag = true;      // boolean flag used to evaluate if statement later on
		double s = 0;           // Our Root that will be returned
		double d = 0;           // Only used if mflag is unset (mflag == false)

		for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
		{
			// stop if converged on root or error is less than tolerance
			if (std::abs(b - a) < TOL) {
				goto NarrowResult;
			} // end if

			if (fa != fc && fb != fc) {
				// use inverse quadratic interopolation
				s = (a * fb * fc / ((fa - fb) * (fa - fc)))
					+ (b * fa * fc / ((fb - fa) * (fb - fc)))
					+ (c * fa * fb / ((fc - fa) * (fc - fb)));
			}
			else { s = b - fb * (b - a) / (fb - fa); }//secant method

													  /*	  (condition 1) s is not between  (3a+b)/4  and b or
													  (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
													  (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
													  (condition 4) (mflag is set and |b−c| < |TOL|) or
													  (condition 5) (mflag is false and |c−d| < |TOL|)
													  */
			if (((s < (3 * a + b) * 0.25) || (s > b)) ||
				(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
				(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
				(mflag && (std::abs(b - c) < TOL)) ||
				(!mflag && (std::abs(c - d) < TOL)))
			{
				// bisection method
				s = (a + b)*0.5;
				mflag = true;
			}
			else { mflag = false; }

			fs = p.HornerEvaluate(s);  // calculate fs
			d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
			c = b;      // set c equal to upper bound
			fc = fb;    // set f(c) = f(b)

			if (fa * fs < 0) {// fa and fs have opposite signs
				b = s;
				fb = fs;    // set f(b) = f(s)
			}
			else {
				a = s;
				fa = fs;    // set f(a) = f(s)
			}

			if (std::abs(fa) < std::abs(fb)) {// if magnitude of fa is less than magnitude of fb
				std::swap(a, b);     // swap a and b
				std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
			}
		}
		return 0;

	NarrowResult:
		//use root as guess for Newton's method now
		double x, x1; //user input value
		double f, dfdx; //function and its derivative 
		const double epsilon = 0.0000001; //tolerance
		const int n = 15; //# steps
		Polynomial deriv = p.derivative();//derivative of polynomial

		x1 = s;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.HornerEvaluate(x); //Defines the function 
			dfdx = deriv.HornerEvaluate(x); //Derivative of the function
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (abs(x - x1) < epsilon) { //checks if guesses are within a certain tolerance value 
			return x1;
		}

		return 0;
	}

	double findRoot(Polynomial p, double lower_bound, double upper_bound) {//bounded search for a root
																		   //given a set of bounds for a polynomial, this method will find the first root
		if (p.terms <= 1) { return 0; }
		if (p.terms == 2) { return (-p.coefficient[0] / p.coefficient[1]); }
		if (p.terms == 3) {
			double part1 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part1 >= 0) { return ((-p.coefficient[1] + sqrt(part1)) / (2 * p.coefficient[2])); }
			else { return 0; }
		}

		//else, p is of order 3 or higher....
		double TOL = 0.0001;    // tolerance
		double MAX_ITER = 1000; // maximum number of iterations
		double a = lower_bound;
		double b = upper_bound;
		double fa = p.HornerEvaluate(a);   // calculated now to save function calls
		double fb = p.HornerEvaluate(b);   // calculated now to save function calls
		double fs = 0;      // initialize 

		if (!(fa * fb < 0)) { return 0; }

		if (fabsf(fa) < fabsf(b)) { // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
			std::swap(a, b);
			std::swap(fa, fb);
		}

		double c = a;           // c now equals the largest magnitude of the lower and upper bounds
		double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
		bool mflag = true;      // boolean flag used to evaluate if statement later on
		double s = 0;           // Our Root that will be returned
		double d = 0;           // Only used if mflag is unset (mflag == false)

		for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
		{
			// stop if converged on root or error is less than tolerance
			if (std::abs(b - a) < TOL) {
				goto NarrowResult;
			} // end if

			if (fa != fc && fb != fc) {
				// use inverse quadratic interopolation
				s = (a * fb * fc / ((fa - fb) * (fa - fc)))
					+ (b * fa * fc / ((fb - fa) * (fb - fc)))
					+ (c * fa * fb / ((fc - fa) * (fc - fb)));
			}
			else { s = b - fb * (b - a) / (fb - fa); }//secant method

													  /*
													  Crazy condition statement!:
													  -------------------------------------------------------
													  (condition 1) s is not between  (3a+b)/4  and b or
													  (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
													  (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
													  (condition 4) (mflag is set and |b−c| < |TOL|) or
													  (condition 5) (mflag is false and |c−d| < |TOL|)
													  */
			if (((s < (3 * a + b) * 0.25) || (s > b)) ||
				(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
				(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
				(mflag && (std::abs(b - c) < TOL)) ||
				(!mflag && (std::abs(c - d) < TOL)))
			{
				// bisection method
				s = (a + b)*0.5;
				mflag = true;
			}
			else { mflag = false; }

			fs = p.HornerEvaluate(s);  // calculate fs
			d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
			c = b;      // set c equal to upper bound
			fc = fb;    // set f(c) = f(b)

			if (fa * fs < 0) {// fa and fs have opposite signs
				b = s;
				fb = fs;    // set f(b) = f(s)
			}
			else {
				a = s;
				fa = fs;    // set f(a) = f(s)
			}

			if (std::abs(fa) < std::abs(fb)) {// if magnitude of fa is less than magnitude of fb
				std::swap(a, b);     // swap a and b
				std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
			}
		}
		return 0;

	NarrowResult:
		//use root as guess for Newton's method now
		double x, x1; //user input value
		double f, dfdx; //function and its derivative 
		const double epsilon = 0.0000001; //tolerance
		const int n = 15; //# steps
		Polynomial deriv = p.derivative();//derivative of polynomial

		x1 = s;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.HornerEvaluate(x); //Defines the function 
			dfdx = deriv.HornerEvaluate(x); //Derivative of the function
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (abs(x - x1) < epsilon) { //checks if guesses are within a certain tolerance value 
			return x1;
		}

		return 0;
	}

	std::vector<double> findRealRoots(Polynomial p) {
		std::vector<double> candidates;
		if (p.terms <= 1) { return candidates; }
		if (p.terms == 2) {
			candidates.push_back(-p.coefficient[0] / p.coefficient[1]);
			return candidates;
		}
		if (p.terms == 3) {
			double part2 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part2 >= 0) {
				double part1 = (-p.coefficient[1]) / (2 * p.coefficient[2]);
				part2 = sqrt(part2) / (2 * p.coefficient[2]);
				candidates.push_back(part1 + part2);
				candidates.push_back(part1 - part2);
			}
			return candidates;
		}
		//if the polynomial is longer then 2nd order...
		std::vector<double> bounds = p.getRootBounds();
		double lowbound = bounds[0];
		double highbound = bounds[1];

		//Since 0 is used as a null value for roots in the rootfinding algorithms,
		//we must determine whether or not it is a root outside of these algorithms.
		if (p.HornerEvaluate(0) == 0) { candidates.push_back(0); }

		//now sweep through the bounds of the polynomial to look for roots using the
		//mixed Brent/Newton search method in findRoot()
		Polynomial temp = p;
		double lastroot = 0;
		double delta = 0.01;
		if (getLength(p.largestCoefficient()) > 3) { delta *= 10 * getLength(p.largestCoefficient()); }
		while (lowbound < highbound) {
			double x = findRoot(temp, lowbound, lowbound + delta);
			if (x == 0) { lowbound += delta; }//using 0 as a null return value here
			else {
				if (x != lastroot) { //prevent same root from continually showing up
					candidates.push_back(x);
					lastroot = x;
					lowbound = x + delta;
					if (x == floor(x)) {
						temp = temp.factorOutBinomial(x); //perform a shift to reduce the polynomial
						bounds = temp.getRootBounds();
						double lowbound = bounds[0];
						double highbound = bounds[1];
					}//if x is an integer, you can factor it out of the polynomial
				}
				else { lowbound += delta; }
			}
		}
		return candidates;
	}

	Polynomial HermitePolynomial(int n) {
		if (n < 0) { return Polynomial(); }
		std::vector<double> H0; // H0 = 1
		H0.push_back(1);
		if (n == 0) { return Polynomial(H0); }
		std::vector<double> H1; //H1 = 2x
		H1.push_back(0);
		H1.push_back(2);
		if (n == 1) { return Polynomial(H1); }
		std::vector<double> H2;	//H2 = 4x^2 -2
		H2.push_back(-2);
		H2.push_back(0);
		H2.push_back(4);
		if (n == 2) { return Polynomial(H2); }

		//if n>2, determine polynomial by recursive meUhod:
		Polynomial Hnext;
		Polynomial Hprev = H2;
		Polynomial Hprevprev = H1;

		std::vector<double> temp;
		temp.push_back(0);
		temp.push_back(2);
		Polynomial mult(temp);//multiplicand is p(x) = 2x

		for (int i = 3; i <= n; ++i) {
			// loop:  Hn+1 = (2x*Hn) - (2n*Hn-1)
			Hnext = Hprev.subtract(Hprev.multiply(mult, Hprev), Hprevprev.multiply(Hprevprev, 2 * (i - 1)));
			Hprevprev = Hprev;
			Hprev = Hnext;
		}
		return Hnext;
	}

	double HermiteNumber(int n) {
		if (n < 0) { n = abs(n); }
		return HermitePolynomial(n).evaluate(0);
	}

	Polynomial LaguerrePolynomial(int n) {
		if (n < 0) { return Polynomial(); }
		std::vector<double> coef;
		for (int i = 0; i <= n; ++i) {
			coef.push_back((pow(-1, i) / factorial(i)) * combination(n, i));
		}
		return Polynomial(coef);
	}

	Polynomial LagrangeInterpolatingPolynomial(std::vector<double> x, std::vector<double> y) {
		if (x.size() != y.size()) { return Polynomial(); }
		Polynomial p;
		for (int i = 0; i < x.size(); ++i) {
			double denom = 1;
			Polynomial num;
			bool first = true;
			for (int j = 0; j < x.size(); ++j) {
				if (j != i) {
					denom *= (x[i] - x[j]);
					if (first == true) {
						//initialize first Polynomial -- (x - x_j)
						std::vector<double> coef;
						coef.push_back(-x[j]);
						coef.push_back(1);
						num = Polynomial(coef);
						first = false;
					}
					else {
						std::vector<double> coef;
						coef.push_back(-x[j]);
						coef.push_back(1);
						num = num.multiply(num, Polynomial(coef));
					}
				}
			}
			denom = 1.0 / denom;
			denom *= y[i];
			num = num.multiply(denom, num);
			if (i == 0) { p = num; }
			else { p = p.add(p, num); }

		}
		return p;
	}

	Polynomial LegendrePolynomial(int n) {
		if (n <= 0) { return Polynomial(); }
		if (n == 10) {
			double coef1[11] = {//special case -- n=10, used for Gauss Quadrature for numerical integration
				-63.0 / 256,
				0,
				3465 / 256.0,
				0,
				-30030 / 256.0,
				0,
				90090 / 256.0,
				0,
				-109395 / 256.0,
				0,
				46189 / 256.0
			};
			return Polynomial(coef1, n + 1);
		}

		double* coef1 = new double[n + 1];
		for (int k = 0; k < n + 1; ++k) { coef1[k] = 0; }//initialize as all '0' values

		double p1 = pow(2, n);
		p1 = 1 / p1;
		for (int k = 0; k <= floor(n / 2); ++k) {
			coef1[n - (2 * k)] = p1 * (pow(-1, k) * combination(n, k) * combination((2 * n) - (2 * k), n));
		}
		return Polynomial(coef1, n + 1);
	}

	double* LegendrePolynomialRootTable(int n) {
		if (n == 2) { double v[2] = { -0.577350269189626 , 0.577350269189626 }; return v; }
		if (n == 3) { double v[3] = { 0, -0.774596669241483 , 0.774596669241483 }; return v; }
		if (n == 7) {
			double v[7] = {
				-0.949107912342759,
				-0.741531185599394,
				-0.405845151377397,
				0,
				0.405845151377397,
				0.741531185599394,
				0.949107912342759
			}; return v;
		}
		if (n == 8) {
			double v[8] = {
				-0.960289856497536,
				-0.796666477413627,
				-0.525532409916329,
				-0.183434642495650,
				0.183434642495650,
				0.525532409916329,
				0.796666477413627,
				0.960289856497536
			}; return v;
		}
		if (n == 9) {
			double v[9] = {
				-0.968160239507626,
				-0.836031107326636,
				-0.613371432700590,
				-0.324253423403809,
				0,
				0.324253423403809,
				0.613371432700590,
				0.836031107326636,
				0.968160239507626
			}; return v;
		}
		if (n == 10) {
			double v[10] = {
				-0.973906528517172,
				-0.865063366688985,
				-0.679409568299024,
				-0.433395394129247,
				-0.148874338981631,
				0.148874338981631,
				0.433395394129247,
				0.679409568299024,
				0.865063366688985,
				0.973906528517172
			}; return v;
		}
		if (n == 14) {
			double v[14] = {
				-0.108054948707344,
				0.108054948707344,
				-0.319112368927890,
				0.319112368927890,
				-0.515248636358154,
				0.515248636358154,
				-0.687292904811685,
				0.687292904811685,
				-0.827201315069765,
				0.827201315069765,
				-0.928434883663574,
				0.928434883663574,
				-0.986283808696812,
				0.986283808696812
			}; return v;
		}
		if (n == 16) {
			double v[16] = {
				-0.095012509837637,
				0.095012509837637,
				-0.281603550779259,
				0.281603550779259,
				-0.458016777657227,
				0.458016777657227,
				-0.617876244402644,
				0.617876244402644,
				-0.755404408355003,
				0.755404408355003,
				-0.865631202387832,
				0.865631202387832,
				-0.944575023073233,
				0.944575023073233,
				-0.989400934991650,
				0.989400934991650 }; return v;
		}
		return NULL;
	}

	std::vector<double> ChristoffelNumbers(int n) {
		std::vector<double> answer;
		if (n <= 0) { answer.push_back(0); return answer; }
		if (n == 1) { answer.push_back(2); return answer; }
		if (n == 2) { answer.push_back(1); answer.push_back(1); return answer; }
		Polynomial p = LegendrePolynomial(n);
		std::vector<double> roots = p.realRoots(p);
		p = p.derivative();
		for (int i = 0; i < roots.size(); ++i) {
			answer.push_back(2 / ((1 - pow(roots[i], 2))*(pow(p.evaluate(roots[i]), 2))));
		}
		return answer;
	}

	std::vector<double> ChristoffelNumbers(std::vector<double> roots, int n) {
		std::vector<double> answer;
		if (n <= 0) { answer.push_back(0); return answer; }
		if (n == 1) { answer.push_back(2); return answer; }
		if (n == 2) { answer.push_back(1); answer.push_back(1); return answer; }
		Polynomial p = LegendrePolynomial(n);
		p = p.derivative();
		for (int i = 0; i < n; ++i) {
			answer.push_back(2 / ((1 - pow(roots[i], 2))*(pow(p.evaluate(roots[i]), 2))));
		}
		return answer;
	}

	std::vector<double> ChristoffelNumbers(double* roots, int n) {
		std::vector<double> answer;
		if (n <= 0) { answer.push_back(0); return answer; }
		if (n == 1) { answer.push_back(2); return answer; }
		if (n == 2) { answer.push_back(1); answer.push_back(1); return answer; }
		Polynomial p = LegendrePolynomial(n);
		p = p.derivative();
		for (int i = 0; i < n; ++i) {
			double tp = roots[i];
			tp *= roots[i];
			tp = 1 - tp;
			tp *= p.evaluate(roots[i]);
			tp *= p.evaluate(roots[i]);
			tp = 2.0 / tp;
			answer.push_back(tp);

			//answer.push_back(2 / ((1 - pow(roots[i], 2))*(pow(p.evaluate(roots[i]), 2))));
		}
		return answer;
	}

	ComplexNumber NewtonsMethodComplex(Polynomial p, ComplexNumber z) {
		//use root as guess for Newton's method now
		ComplexNumber x, x1; //user input value
		ComplexNumber f, dfdx; //function and its derivative 
		const double epsilon = 0.00000001; //tolerance
		const int n = 1000; //# steps
		Polynomial deriv = p.derivative();//derivative of polynomial

		x1 = z;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.HornerEvaluateComplex(x); //Defines the function 
			dfdx = deriv.HornerEvaluateComplex(x); //Derivative of the function
			x1 = x.subtract(x, f.divide(f, dfdx)); //Newton Method formula = x1 = x - f/dfdx
		}

		if (x.Im != 0 && x.subtract(x, x1).Re < epsilon && x.subtract(x, x1).Im < epsilon) { //checks if guesses are within a certain tolerance value AND complex (not real)
			return x1;
		}
		return ComplexNumber();//if not close enough, return 0 value
	}

	Quaternion positionQuaternion(Vector v) {
		Vector v2(4);
		for (int i = 0; i < 3; ++i) { v2.vec.push_back(v.vec[i]); }
		v2.vec.push_back(1);
		return v2.toQuaternion();
	}

	Polynomial reverseBesselPolynomial(int n) {
		std::vector<double> vec;
		for (int i = 0; i <= n; ++i) {
			double p1 = factorial(n + i) / (factorial(n - i)*factorial(i)*pow(2, i));
			vec.push_back(p1);
		}
		return Polynomial(reverse(vec));
	}

	Polynomial RogersSzegoPolynomial(double q, int n) {
		std::vector<double> coef;
		long double qqn = qPochhammerFactorial(q, q, n);
		for (int i = 0; i < n; ++i) {
			long double qqi = qPochhammerFactorial(q, q, i);
			long double qqni = qPochhammerFactorial(q, q, n - i);
			coef.push_back(qqn / (qqi*qqni));
		}
		return Polynomial(coef);
	}

	std::vector<ComplexNumber> roots(Polynomial p) {
		std::vector<double> realroots = findRealRoots(p);
		std::vector<ComplexNumber> complexRoots;
		for (int i = 0; i < realroots.size(); ++i) {
			//std::wcout << realRoots[i] << L"\n";//debug
			complexRoots.push_back(ComplexNumber(realroots[i], 0));
		}

		std::vector<double> bounds = p.getComplexRootBounds();
		double lowbound = bounds[0];
		double highbound = bounds[1];
		double lowboundC = bounds[0];

		//std::wcout << L"root bounds = [" << lowbound << L"," << highbound << L"]\n";//debug
		for (int i = 0; i < realroots.size(); ++i) {
			if (std::abs(p.terms - (int)realroots.size()) <= 3) {
				//Use the shift method only if the real roots reduce it down to an easily calculable size (order 2 or less),
				//otherwise the results will be inaccurate
				if (p.terms <= 3) { i = realroots.size(); break; }
				if (p.terms > 3) { p = p.factorOutBinomial(realroots[i]); }
			}
		}
		p.display();//debug
		if (p.terms < 3) { return complexRoots; }

		if (p.terms == 3) {
			double part2 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part2 < 0) {
				part2 = part2 / (2 * p.coefficient[2]);
				double part1 = (-p.coefficient[1] / (2 * p.coefficient[2]));
				complexRoots.push_back(ComplexNumber(part1, part2));
				complexRoots.push_back(ComplexNumber(part1, -part2));
			}
			return complexRoots;
		}

		if (p.terms > 3) {


			//start finding roots with Durand-Kerner method
			int iterations = 10000;
			ComplexNumber z = ComplexNumber(lowbound, lowboundC);
			int size = sizeof(z);
			std::vector<ComplexNumber> R;
			for (int i = 0; i < p.terms; i++) { R.push_back(z.exponent(z, i)); }

			for (int i = 0; i < iterations; i++) {
				for (int j = 0; j < p.terms; j++) {
					ComplexNumber B = p.HornerEvaluateComplex(R[j]);
					for (int k = 0; k < p.terms; k++) {
						if (k != j) { B /= (R[j] - R[k]); }
					}
					R[j] -= B;
				}
			}


			for (int i = 0; i < R.size(); ++i) { //now filter out all the bad results
				if (R[i].Im != 0) { //only complex roots accepted
					R[i] = NewtonsMethodComplex(p, R[i]);
					ComplexNumber temp = p.HornerEvaluateComplex(R[i]);
					double a = std::abs(temp.Re);
					double b = std::abs(temp.Im);
					if (a < 0.1 && b < 0.1) {
						a = std::abs(R[i].Re);
						b = std::abs(R[i].Im);
						bool isSimilar = false;
						for (int i = 0; i < complexRoots.size(); ++i) {
							double x = std::abs(complexRoots[i].Re);
							double y = std::abs(complexRoots[i].Im);
							if (std::abs(a - x) < 0.001 && std::abs(b - y) < 0.001) { isSimilar = true; }
						}
						if (isSimilar == false) { //if this is indeed a new root, save the root and its conjugate
							complexRoots.push_back(R[i]);
							//R[i].display();
							complexRoots.push_back(ComplexNumber(R[i].Re, R[i].Im*-1));
						}
					}
				}
			}
			R.clear();
		}
		return complexRoots;
	}

	Polynomial StieltjesPolynomial(double n) {
		//ADD CODE HERE
		return Polynomial();
	}

	/*Quaternion Matrix::toQuaternion() {//for use in homogeneous coordinate system for 3D graphics
		if (columns == rows == 4) {//if using a 4x4 matrix representation of a quaternion
			std::vector<double> d = diagonal();
			return Quaternion(d);
		}

		if (columns == 1 && rows == 4) {//if using a single vector
			std::vector<double> c = column(0);
			return Quaternion(c);
		}
		return Quaternion();//else return NULL
	}*/

	Quaternion transformVector(Quaternion T, Quaternion R, Quaternion S, Quaternion Vec) {
		//order is Transformed Matrix = T*(R*(S*Vector)))
		Matrix Tr = T.toTranslationMatrix();
		Matrix Sc = S.toScalingMatrix();
		Matrix Ro = R.toRotationMatrix();
		Matrix V = Vec.toMatrix();
		Matrix M = modelMatrix(Tr, Ro, Sc);
		M *= V;
		return M.toQuaternion();
	}
