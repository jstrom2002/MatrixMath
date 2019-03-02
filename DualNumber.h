#ifndef DUALNUMBER_H
#define DUALNUMBER_H
#pragma once
#include "stdafx.h"

class DualNumber {
	/*Dual numbers are similar to imaginary numbers (i.e. 'a+bd' where d ~ i), but with nilpotent extension d^2 = 0.
	A dual number input to a function will return the derivative in its dual or 'd' component.	*/
public:
	double a;//first, real value
	double b;//dual number
	unsigned int dualNumberPrecision;
	DualNumber();
	DualNumber(double x);
	DualNumber(double x, double y);
	~DualNumber();
	DualNumber add(DualNumber A, DualNumber B);
	DualNumber add(DualNumber A, double B);
	DualNumber subtract(DualNumber A, DualNumber B);
	DualNumber subtract(DualNumber A, double B);
	DualNumber multiply(DualNumber A, DualNumber B);
	DualNumber multiply(DualNumber A, double B);
	DualNumber divide(DualNumber A, DualNumber B);
	DualNumber divide(DualNumber A, double B);
	DualNumber& operator+=(const DualNumber& rhs);
	DualNumber& DualNumber::operator+=(double x);
	const DualNumber DualNumber::operator+(const DualNumber& rhs) const;
//	std::wostream& operator<<(std::wostream& os, const DualNumber& rhs);
	void randomize();
	DualNumber inverse();
	DualNumber conjugate();
	DualNumber exponent(DualNumber z, double a);
	DualNumber exponent(DualNumber d1, DualNumber d2);
	DualNumber exponent(int a);
	DualNumber& operator*=(DualNumber& rhs);
	DualNumber& operator*=(double x);
	DualNumber& operator/=(DualNumber& rhs);
	DualNumber& operator/=(double x);
	DualNumber& operator+=(DualNumber& rhs);
	//DualNumber& operator+=(double x);
	DualNumber& operator-=(DualNumber& rhs);
	DualNumber& operator-=(double x);
	DualNumber& operator*(DualNumber& rhs);
	DualNumber& operator*(double x);
	DualNumber& operator/(DualNumber& rhs);
	DualNumber& operator/(double x);
	DualNumber& operator+(DualNumber& rhs);
	DualNumber& operator+(double x);
	DualNumber& operator-(DualNumber& rhs);
	DualNumber& operator-(double x);
	DualNumber& operator^(int x);
	DualNumber& operator^(double x);
	DualNumber& operator^(DualNumber x);
	std::wstring toString();
	void display();
};//==============END DUAL NUMBER CLASS======================
DualNumber exp(DualNumber d);
DualNumber log(DualNumber d);
DualNumber sin(DualNumber d);
DualNumber cos(DualNumber d);
DualNumber tan(DualNumber d);
DualNumber sec(DualNumber d);
DualNumber csc(DualNumber d);
DualNumber cot(DualNumber d);
DualNumber sinh(DualNumber d);
DualNumber cosh(DualNumber d);
DualNumber tanh(DualNumber d);
DualNumber sech(DualNumber d);
DualNumber csch(DualNumber d);
DualNumber coth(DualNumber d);
#endif