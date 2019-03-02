#pragma once
#include "stdafx.h"

DualNumber::DualNumber() {}
DualNumber::DualNumber(double x) {
	a = x;
	b = 0;
	int c = getLength(a);
	dualNumberPrecision = c;
}
DualNumber::DualNumber(double x, double y) {
	a = x;
	b = y;
	int c = getLength(a);
	if (getLength(b) > c) { c = getLength(b); }
	dualNumberPrecision = c;
}
DualNumber::~DualNumber() {}

	   DualNumber DualNumber::add(DualNumber A, DualNumber B) { return DualNumber(A.a + B.a, A.b + B.b); }
	   DualNumber DualNumber::add(DualNumber A, double B) { return DualNumber(A.a + B, A.b + B); }
	   DualNumber DualNumber::subtract(DualNumber A, DualNumber B) { return DualNumber(A.a - B.a, A.b - B.b); }
	   DualNumber DualNumber::subtract(DualNumber A, double B) { return DualNumber(A.a - B, A.b - B); }
	   DualNumber DualNumber::multiply(DualNumber A, DualNumber B) { return DualNumber(A.a*B.a, (A.a*B.b) + (A.b*B.a)); }
	   DualNumber DualNumber::multiply(DualNumber A, double B) { return DualNumber(A.a*B, A.b*B); }
	   DualNumber DualNumber::divide(DualNumber A, DualNumber B) {
		   if (B.a == 0) { return DualNumber(); }
		   return DualNumber(A.a / B.a, (A.b / B.a) - ((A.a*B.b) / (B.a*B.a)));
	   }
	   DualNumber DualNumber::divide(DualNumber A, double B) { if (B != 0) { return DualNumber(A.a*B, A.b / B); } }

	   DualNumber& DualNumber::operator+=(const DualNumber& rhs) {
		   a += rhs.a;
		   b += rhs.b;
		   return *this;
	   }

	   DualNumber& DualNumber::operator+=(double x) {
		   a += x;
		   return *this;
	   }

	   const DualNumber DualNumber::operator+(const DualNumber& rhs) const {
		   DualNumber sum = *this;
		   sum.a += rhs.a;
		   sum.b += rhs.b;
		   return sum;
	   }
	   std::wostream& operator<<(std::wostream& os, const DualNumber& rhs) {
		   os << rhs.a << L"+" << rhs.b << L"ε";
		   return os;
	   }

	   void DualNumber::randomize() {
		   srand(time(0) + clock());
		   double a2 = (rand() / 1000000.00) * pow(-1, rand() % 2);
		   a = a2;
		   a2 = (rand() / 1000000.00) * pow(-1, rand() % 2);
		   b = a2;
	   }

	   DualNumber DualNumber::inverse() { return DualNumber(1.0 / a, b / (a*a)); }
	   DualNumber DualNumber::conjugate() { return DualNumber(a, -b); }

	   DualNumber DualNumber::exponent(int a) {
		   DualNumber z = *this;
		   if (a < 0) {
			   z = inverse();
			   a = abs(a);
		   }
		   return DualNumber(pow(z.a, a), a*pow(z.a, a - 1)*z.b);
		   if (a == 0) { return DualNumber(1, 0); }
	   }

	   DualNumber DualNumber::exponent(DualNumber z, double a) {
		   DualNumber c = z;
		   if (a < 0) {
			   c = c.inverse();
			   a = abs(a);
		   }
		   if (a > 0) {
			   for (int i = 0; i < a - 1; i++) {
				   c = c.multiply(c, z);
			   }
			   return c;
		   }
		   if (a == 0) { return DualNumber(1, 0); }
	   }
	   DualNumber DualNumber::exponent(DualNumber d1, DualNumber d2) { return DualNumber(d1.a, d1.b + (d1.a*d1.a*d2.b) - (d1.a*d2.b)); }

	   DualNumber& DualNumber::operator*=(DualNumber& rhs) { *this = multiply(*this, rhs); return *this; }
	   DualNumber& DualNumber::operator*=(double x) { *this = multiply(*this, x); return *this; }
	   DualNumber& DualNumber::operator/=(DualNumber& rhs) { *this = divide(*this, rhs); return *this; }
	   DualNumber& DualNumber::operator/=(double x) { *this = divide(*this, x); return *this; }
	   DualNumber& DualNumber::operator+=(DualNumber& rhs) { *this = add(*this, rhs); return *this; }
	   //DualNumber& DualNumber::operator+=(double x) { *this = add(*this, x); return *this; }
	   DualNumber& DualNumber::operator-=(DualNumber& rhs) { *this = subtract(*this, rhs); return *this; }
	   DualNumber& DualNumber::operator-=(double x) { *this = subtract(*this, x); return *this; }

	   DualNumber& DualNumber::operator*(DualNumber& rhs) { return multiply(*this, rhs); }
	   DualNumber& DualNumber::operator*(double x) { DualNumber p = multiply(*this, x); return p; }
	   DualNumber& DualNumber::operator/(DualNumber& rhs) { return divide(*this, rhs); }
	   DualNumber& DualNumber::operator/(double x) { DualNumber p = divide(*this, x); return p; }
	   DualNumber& DualNumber::operator+(DualNumber& rhs) { return add(*this, rhs); }
	   DualNumber& DualNumber::operator+(double x) { return add(*this, x); }
	   DualNumber& DualNumber::operator-(DualNumber& rhs) { return subtract(*this, rhs); }
	   DualNumber& DualNumber::operator-(double x) { return subtract(*this, x); }
	   DualNumber& DualNumber::operator^(int x) { return exponent(*this, x); }
	   DualNumber& DualNumber::operator^(double x) { return exponent(*this, x); }
	   DualNumber& DualNumber::operator^(DualNumber x) { return exponent(*this, x); }

	   std::wstring DualNumber::toString() {
		   std::wstring answer;
		   std::wostringstream strs;
		   dualNumberPrecision = precision;
		   if (a == floor(a)) { dualNumberPrecision = 0; }
		   strs << std::fixed << std::setprecision(dualNumberPrecision) << a;
		   answer.append(strs.str());
		   dualNumberPrecision = precision;
		   if (b == floor(b)) { dualNumberPrecision = 0; }
		   if (b != 0) {
			   answer.append(L" + ");
			   if (b != 1) {
				   std::wostringstream strs2;
				   strs2 << std::fixed << std::setprecision(dualNumberPrecision) << b;
				   answer.append(strs2.str());
			   }
			   answer.append(L"ε");
		   }
		   return answer;
	   }

	   void DualNumber::display() {
		   std::wstring temp = L"ε";
		   std::wcout << a << L"+" << b << temp << std::endl;
	   }
//==============END DUAL NUMBER CLASS FUNCTIONS======================


  //Dual Number valued functions
  //============================
DualNumber exp(DualNumber d) { return DualNumber(exp(d.a), exp(d.a)*(1 + d.b)); }
DualNumber log(DualNumber d) { return DualNumber(log(d.a), d.b / d.a); }

DualNumber sin(DualNumber d) {
	double re = 1;
	double du = 0;
	for (int i = 1; i < 200; ++i) {
		double divisor = factorial(2 * i)*pow(-1, i);
		re += pow(d.a, 2 * i) / divisor;
		du += (pow(d.a, 2 * i - 1) * d.b) / divisor;
	}
	return DualNumber(re, du);
}
DualNumber cos(DualNumber d) {
	double re = 0;
	double du = 0;
	for (int i = 1; i < 200; ++i) {
		double divisor = factorial(2 * i - 1)*pow(-1, i - 1);
		re += pow(d.a, 2 * i - 1) / divisor;
		du += (pow(d.a, 2 * i - 2) * d.b) / divisor;
	}
	return DualNumber(re, du);
}
DualNumber tan(DualNumber d) { return sin(d) / cos(d); }
DualNumber sec(DualNumber d) { return cos(d).inverse(); }
DualNumber csc(DualNumber d) { return sin(d).inverse(); }
DualNumber cot(DualNumber d) { return cos(d) / sin(d); }
DualNumber sinh(DualNumber d) { return DualNumber(0.5 * (exp(d.a) - (1.0 / exp(d.a))), 0.5*d.b*(exp(d.a) + (1.0 / exp(d.a)))); }
DualNumber cosh(DualNumber d) { return DualNumber(0.5 * (exp(d.a) + (1.0 / exp(d.a))), 0.5*d.b*(exp(d.a) - (1.0 / exp(d.a)))); }
DualNumber tanh(DualNumber d) { return sinh(d) / cosh(d); }
DualNumber sech(DualNumber d) { return cosh(d).inverse(); }
DualNumber csch(DualNumber d) { return sinh(d).inverse(); }
DualNumber coth(DualNumber d) { return cosh(d) / sinh(d); }
