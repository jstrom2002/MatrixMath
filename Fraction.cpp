#pragma once
#include "stdafx.h"
#include "discrete math.h"

Fraction::Fraction() { num = 0; denom = 0; }
Fraction::~Fraction() {};

Fraction::Fraction(int n, int m) {
	if (m == 0 || n == 0) { num = 0; denom = 0; }
	else {
		num = n;
		denom = m;
		if (n < 0 && m < 0) { 
			num = abs(num); 
			denom = abs(denom); 
		}
	}
}

Fraction::Fraction(long long int n, long long int m) {
	if (m == 0 || n == 0) { num = 0; denom = 0; }
	else {
		num = n;
		denom = m;
		if (n < 0 && m < 0) {
			num = abs(num);
			denom = abs(denom);
		}
	}
}

Fraction::Fraction(double n, double m) {
	if (m == 0 || n == 0) { num = 0; denom = 0; }
	else {
		num = (int)n;
		denom = (int)m;
		if (n < 0 && m < 0) {
			num = abs(num);
			denom = abs(denom);
		}
	}
}

Fraction::Fraction(int n) {
	num = n;
	denom = 1;
}

Fraction::Fraction(long long int n) {
	num = n;
	denom = 1;
}

Fraction::Fraction(double x) {
	if (floor(x) == x) { num = (int)x; denom = 1; return; }
	bool isNegative = false;
	double n = x;
	if (n < 0) {
		n = abs(n);
		isNegative = true;
	}
	int num = 0;
	int denom = 1;
	double start = floor(x);
	//testing integer division loop
	for (double i = start; i < 99999; ++i) {
		for (double j = 1; j < 99999; ++j) {
			if (abs(n - (i / j)) < 0.0000001) {  //this is how close you want to approximate the Fraction
				num = (int)i;
				denom = (int)j;
				i = j = 999999999999;
			}
		}
	}
	if (isNegative == true) {
		num *= -1;
	}
}

	   //FRACTION METHODS
	   //================
	   Fraction Fraction::simplify() {
		   int d = GCD(num, denom);
		   if (d <= abs(num) && d <= abs(denom)) { 
			   return Fraction(num/d, denom/d); 
		   }
		   return Fraction(num, denom);
	   }

	   Fraction Fraction::add(Fraction fr1, Fraction fr2) {
		   int n1, n2, d1, d2;
		   d1 = fr1.denom;
		   d2 = fr2.denom;
		   n1 = fr1.num;
		   n2 = fr2.num;
		   d1 *= d2;
		   n1 *= d2;           //multiply both by the other's denom
		   d2 *= fr1.denom;
		   n2 *= fr1.denom;
		   n1 += n2;
		   Fraction answer(n1, d1);
		   answer.simplify();
		   return answer;
	   }

	   Fraction Fraction::add(Fraction fr1, double fr2) {
		   Fraction F(fr2);
		   return add(fr1, F);
	   }

	   Fraction Fraction::subtract(Fraction fr1, Fraction fr2) {
		   if (fr1.num == 0 || fr1.denom == 0) { 
			   return Fraction(fr2.num*-1, fr2.denom); 
		   }
		   if (fr2.num == 0 || fr2.denom == 0) {
			   return fr1;
		   }
		   double nm = ((double)fr2.denom*fr1.num) - (fr1.denom*fr2.num);
		   Fraction answer(nm, (double)fr1.denom*fr2.denom);
		   answer.simplify();
		   return answer;
	   }

	   Fraction Fraction::subtract(Fraction fr1, double fr2) {
		   Fraction F(fr2,1.0);
		   return subtract(fr1, F);
	   }

	   Fraction Fraction::multiply(Fraction fr1, Fraction fr2) {
		   Fraction answer(fr1.num*fr2.num, fr1.denom*fr2.denom);
		   answer.simplify();
		   return answer;
	   }

	   Fraction Fraction::multiply(Fraction fr1, double fr2) {
		   return Fraction((long long int)(fr1.num*fr2), fr1.denom);
	   }

	   Fraction Fraction::multiply(Fraction fr1, int fr2) {
		   return Fraction(fr1.num*fr2, fr1.denom);
	   }

	   Fraction Fraction::divide(Fraction fr1, Fraction fr2) {
		   // a/b / c/d  =  a/b * d/c 
		   Fraction answer(fr1.num*fr2.denom, fr1.denom*fr2.denom);
		   answer.simplify();
		   return answer;
	   }
	   Fraction Fraction::divide(Fraction fr1, double fr2) {
		   Fraction F(fr2);
		   return divide(fr1, F);
	   }

	   Fraction Fraction::exponent(double d) { return Fraction(pow(num, d), pow(denom, d)); }
	   Fraction Fraction::exponent(Fraction d) { return Fraction(pow(num, d.num / d.denom), pow(denom, d.num / d.denom)); }

	   //OVERLOADED OPERATORS
	   //====================
	   Fraction& Fraction::operator*=(Fraction& rhs) { *this = multiply(*this, rhs); return *this; }
	   Fraction& Fraction::operator*=(double x) { *this = multiply(*this, x); return *this; }
	   Fraction& Fraction::operator/=(Fraction& rhs) { *this = divide(*this, rhs); return *this; }
	   Fraction& Fraction::operator/=(double x) { *this = divide(*this, x); return *this; }
	   Fraction& Fraction::operator+=(Fraction& rhs) { *this = add(*this, rhs); return *this; }
	   Fraction& Fraction::operator+=(double x) { *this = add(*this, x); return *this; }
	   Fraction& Fraction::operator-=(Fraction& rhs) { *this = subtract(*this, rhs); return *this; }
	   Fraction& Fraction::operator-=(double x) { *this = subtract(*this, x); return *this; }

	   Fraction& Fraction::operator*(Fraction& rhs) { return multiply(*this, rhs); }
	   Fraction& Fraction::operator*(double x) { Fraction p = multiply(*this, x); return p; }
	   Fraction& Fraction::operator/(Fraction& rhs) { return divide(*this, rhs); }
	   Fraction& Fraction::operator/(double x) { Fraction p = divide(*this, x); return p; }
	   Fraction& Fraction::operator+(Fraction& rhs) { return add(*this, rhs); }
	   Fraction& Fraction::operator+(double x) { return add(*this, x); }
	   Fraction& Fraction::operator-(Fraction& rhs) { return subtract(*this, rhs); }
	   Fraction& Fraction::operator-(double x) { return subtract(*this, x); }
	   Fraction& Fraction::operator^(Fraction& rhs) { return exponent(rhs); }
	   Fraction& Fraction::operator^(double x) { return exponent(x); }

	   std::wstring Fraction::toString() {
		   std::wstring answer = L"";
		   std::wostringstream strs;
		   strs << num;
		   answer.append(strs.str());
		   answer.append(L"/");
		   std::wostringstream strs2;
		   strs2 << denom;
		   answer.append(strs2.str());
		   strs.clear();
		   strs2.clear();
		   return answer;
	   }

	   void Fraction::display() { std::wcout << toString() << std::endl; }

	   void Fraction::printToFile(std::wstring filename) {
		   std::wofstream outf(filename, std::ios_base::app);
		   outf << toString() << L"\n";
		   outf.close();
	   }