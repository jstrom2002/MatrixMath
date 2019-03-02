#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#pragma once
#include "stdafx.h"

class Polynomial {	//if declared with array = {1,2,3,...} will be 1 + 2x + 3x^2 + ...
public:
	double* coefficient;
	char variable;		//symbol for Polynomial i.e. 't^2 + 2'
	int terms;			//length is defined as the greatest exponent in the Polynomial.
	unsigned int polynomialPrecision;
	Polynomial(double* coef, int order);
	Polynomial(long double* coef, int order);
	Polynomial(std::vector<double> coef);
	Polynomial(std::vector<double> coef, char x);
	Polynomial(double* coef, int order, char var);
	Polynomial();
	~Polynomial();
	void randomize();
	double largestCoefficient();
	double evaluate(double n);
	ComplexNumber evaluate(ComplexNumber z);
	double HornerEvaluate(double x);
	ComplexNumber HornerEvaluateComplex(ComplexNumber z);
	double NewtonsMethod(double x);
	double evaluateY(double y);
	Polynomial add(Polynomial a, Polynomial b);
	Polynomial add(Polynomial a, double  b);
	Polynomial add(double b, Polynomial a);
	Polynomial subtract(Polynomial a, Polynomial b);
	Polynomial subtract(Polynomial a, double  b);
	Polynomial subtract(double b, Polynomial a);
	Polynomial multiply(Polynomial a, Polynomial b);
	Polynomial multiply(Polynomial a, double b);
	Polynomial multiply(double b, Polynomial a);
	Polynomial divide(Polynomial p, Polynomial q);
	Polynomial divide(Polynomial p2, double q);
	Polynomial exponentiate(Polynomial p, int b);
	Polynomial& operator*=(Polynomial& rhs);
	Polynomial& operator*=(double x);
	Polynomial& operator/=(Polynomial& rhs);
	Polynomial& operator/=(double x);
	Polynomial& operator+=(Polynomial& rhs);
	Polynomial& operator+=(double x);
	Polynomial& operator-=(Polynomial& rhs);
	Polynomial& operator-=(double x);
	Polynomial& operator*(Polynomial& rhs);
	Polynomial& operator*(double x);
	Polynomial& operator/(Polynomial& rhs);
	Polynomial& operator/(double x);
	Polynomial& operator+(Polynomial& rhs);
	Polynomial& operator+(double x);
	Polynomial& operator-(Polynomial& rhs);
	Polynomial& operator-(double x);
	friend std::wostream& operator<<(std::wostream& os, Polynomial& rhs);
	bool isIntegerValued();
	bool isMonic();
	Polynomial makeMonic();
	Polynomial factorOutBinomial(Polynomial p);
	Polynomial factorOutBinomial(double n);
	double syntheticDivision(Polynomial p);
	double syntheticDivision(double n);
	ComplexNumber complexSyntheticDivision(ComplexNumber z);
	bool checkLowBound(double n);
	bool checkLowBound(Polynomial p);
	Polynomial derivative();
	Polynomial integral();
	Polynomial CauchyProduct(Polynomial A, Polynomial B);
	std::vector<double> getRootBounds();
	double root(Polynomial p, double lower_bound, double upper_bound);
	std::vector<double> realRoots(Polynomial p);
	std::vector<double> realRoots();
	std::vector<double> getComplexRootBounds();
	std::wstring toString();
	void printToFile(std::wstring filename);
	void display();
};//END Polynomial CLASS=========================================================

//function prototypes
Polynomial BellPolynomial(int n);
Polynomial BernoulliPolynomial(int n);
Polynomial BesselPolynomial(int n);
int complementaryBellNumber(int n);
Polynomial ChebyshevPolynomial1stKind(int n);
Polynomial ChebyshevPolynomial2ndKind(int n);
Matrix companionMatrix(Polynomial p);
Matrix companionMatrix2(Polynomial p);
Quaternion directionQuaternion(Vector v);
Polynomial EulerPolynomial(int m);
double EulerNumber(int n);
double findRoot(Polynomial p);
double findRoot(Polynomial p, double lower_bound, double upper_bound);
std::vector<double> findRealRoots(Polynomial p);
Polynomial HermitePolynomial(int n);
double HermiteNumber(int n);
Polynomial LaguerrePolynomial(int n);
Polynomial LagrangeInterpolatingPolynomial(std::vector<double> x, std::vector<double> y);
Polynomial LegendrePolynomial(int n);
double* LegendrePolynomialRootTable(int n);
std::vector<double> ChristoffelNumbers(int n);
std::vector<double> ChristoffelNumbers(std::vector<double> roots, int n);
std::vector<double> ChristoffelNumbers(double* roots, int n);
ComplexNumber NewtonsMethodComplex(Polynomial p, ComplexNumber z);
Quaternion positionQuaternion(Vector v);
Polynomial reverseBesselPolynomial(int n);
Polynomial RogersSzegoPolynomial(double q, int n);
std::vector<ComplexNumber> roots(Polynomial p);
Polynomial StieltjesPolynomial(double n);
Quaternion transformVector(Quaternion T, Quaternion R, Quaternion S, Quaternion Vec);
#endif