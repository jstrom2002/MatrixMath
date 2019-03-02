#ifndef FRACTION_H
#define FRACTION_H
#pragma once
#include "stdafx.h"
#include "discrete math.h"

class Fraction {
public:	
	long long int num;
	long long int denom;
	Fraction();
	~Fraction();
	Fraction(int n, int m);
	Fraction(long long int n, long long int m);
	Fraction(double n, double m);
	Fraction(int n);
	Fraction(long long int n);
	Fraction(double x);
	Fraction simplify();
	Fraction add(Fraction fr1, Fraction fr2);
	Fraction add(Fraction fr1, double fr2);
	Fraction subtract(Fraction fr1, Fraction fr2);
	Fraction subtract(Fraction fr1, double fr2);
	Fraction multiply(Fraction fr1, Fraction fr2);
	Fraction multiply(Fraction fr1, double fr2);
	Fraction multiply(Fraction fr1, int fr2);
	Fraction divide(Fraction fr1, Fraction fr2);
	Fraction divide(Fraction fr1, double fr2);
	Fraction exponent(double d);
	Fraction exponent(Fraction d);
	Fraction& operator*=(Fraction& rhs);
	Fraction& operator*=(double x);
	Fraction& operator/=(Fraction& rhs);
	   Fraction& operator/=(double x);
	   Fraction& operator+=(Fraction& rhs);
	   Fraction& operator+=(double x);
	   Fraction& operator-=(Fraction& rhs);
	   Fraction& operator-=(double x);
	   Fraction& operator*(Fraction& rhs);
	   Fraction& operator*(double x);
	   Fraction& operator/(Fraction& rhs);
	   Fraction& operator/(double x);
	   Fraction& operator+(Fraction& rhs);
	   Fraction& operator+(double x);
	   Fraction& operator-(Fraction& rhs);
	   Fraction& operator-(double x);
	   Fraction& operator^(Fraction& rhs);
	   Fraction& operator^(double x);
	   std::wstring toString();
	   void display();
	   void printToFile(std::wstring filename);
};
#endif