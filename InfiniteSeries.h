#ifndef INFINITESERIES_H
#define INFINITESERIES_H
#pragma once
#include "stdafx.h"

class InfiniteSeries {
public:
	Function f;
	bool isProductSeries;
	InfiniteSeries();
	~InfiniteSeries();
	InfiniteSeries(Function f1);
	InfiniteSeries(std::vector<wchar_t> var, std::wstring funct);
	InfiniteSeries add(InfiniteSeries a, InfiniteSeries b);
	InfiniteSeries add(InfiniteSeries a, double b);
	InfiniteSeries subtract(InfiniteSeries a, InfiniteSeries b);
	InfiniteSeries subtract(InfiniteSeries a, double b);
	InfiniteSeries multiply(InfiniteSeries a, InfiniteSeries b);
	InfiniteSeries multiply(InfiniteSeries a, double b);
	InfiniteSeries divide(InfiniteSeries a, InfiniteSeries b);
	InfiniteSeries divide(InfiniteSeries a, double b);
	InfiniteSeries inverse(InfiniteSeries a);
	InfiniteSeries exponent(InfiniteSeries a, InfiniteSeries b);
	InfiniteSeries exponent(InfiniteSeries a, double b);
	InfiniteSeries& operator*=(InfiniteSeries& rhs);
	InfiniteSeries& operator*=(double& rhs);
	InfiniteSeries& operator/=(InfiniteSeries& rhs);
	InfiniteSeries& operator/=(double& rhs);
	InfiniteSeries& operator+=(InfiniteSeries& rhs);
	InfiniteSeries& operator+=(double& rhs);
	InfiniteSeries& operator-=(InfiniteSeries& rhs);
	InfiniteSeries& operator-=(double& rhs);
	InfiniteSeries& operator*(InfiniteSeries& rhs);
	InfiniteSeries& operator*(double& rhs);
	InfiniteSeries& operator/(InfiniteSeries& rhs);
	InfiniteSeries& operator/(double& rhs);
	InfiniteSeries& operator+(InfiniteSeries& rhs);
	InfiniteSeries& operator+(double& rhs);
	InfiniteSeries& operator-(InfiniteSeries& rhs);
	InfiniteSeries& operator-(double& rhs);
	InfiniteSeries& operator^(InfiniteSeries x);
	InfiniteSeries& operator^(double x);
	double evaluate(int lo, int hi);
	double evaluate(std::vector<double> vals, int lo, int hi);
	double evaluateLoToInf(double lo, double hi);
	void setProductSeries();
	std::wstring testConvergence(unsigned int n);
	std::wstring InfiniteSeries::toString();
	void display();
};

class ComplexInfiniteSeries {
public:
	ComplexFunction cf;
	bool isProductSeries;
	ComplexInfiniteSeries();
	~ComplexInfiniteSeries();
	ComplexInfiniteSeries(ComplexFunction f1);
	ComplexInfiniteSeries(std::vector<wchar_t> var, std::wstring funct);
	ComplexInfiniteSeries add(ComplexInfiniteSeries a, ComplexInfiniteSeries b);
	ComplexInfiniteSeries add(ComplexInfiniteSeries a, double b);
	ComplexInfiniteSeries subtract(ComplexInfiniteSeries a, ComplexInfiniteSeries b);
	ComplexInfiniteSeries subtract(ComplexInfiniteSeries a, double b);
	ComplexInfiniteSeries multiply(ComplexInfiniteSeries a, ComplexInfiniteSeries b);
	ComplexInfiniteSeries multiply(ComplexInfiniteSeries a, double b);
	ComplexInfiniteSeries divide(ComplexInfiniteSeries a, ComplexInfiniteSeries b);
	ComplexInfiniteSeries divide(ComplexInfiniteSeries a, double b);
	ComplexInfiniteSeries inverse(ComplexInfiniteSeries a);
	ComplexInfiniteSeries exponent(ComplexInfiniteSeries a, ComplexInfiniteSeries b);
	ComplexInfiniteSeries exponent(ComplexInfiniteSeries a, double b);
	ComplexInfiniteSeries& operator*=(ComplexInfiniteSeries& rhs);
	ComplexInfiniteSeries& operator*=(double& rhs);
	ComplexInfiniteSeries& operator/=(ComplexInfiniteSeries& rhs);
	ComplexInfiniteSeries& operator/=(double& rhs);
	ComplexInfiniteSeries& operator+=(ComplexInfiniteSeries& rhs);
	ComplexInfiniteSeries& operator+=(double& rhs);
	ComplexInfiniteSeries& operator-=(ComplexInfiniteSeries& rhs);
	ComplexInfiniteSeries& operator-=(double& rhs);
	ComplexInfiniteSeries& operator*(ComplexInfiniteSeries& rhs);
	ComplexInfiniteSeries& operator*(double& rhs);
	ComplexInfiniteSeries& operator/(ComplexInfiniteSeries& rhs);
	ComplexInfiniteSeries& operator/(double& rhs);
	ComplexInfiniteSeries& operator+(ComplexInfiniteSeries& rhs);
	ComplexInfiniteSeries& operator+(double& rhs);
	ComplexInfiniteSeries& operator-(ComplexInfiniteSeries& rhs);
	ComplexInfiniteSeries& operator-(double& rhs);
	ComplexInfiniteSeries& operator^(ComplexInfiniteSeries x);
	ComplexInfiniteSeries& operator^(double x);
	ComplexNumber evaluate(double lo, double hi);
	ComplexNumber evaluate(ComplexNumber vals2, int lo, int hi);
	ComplexNumber evaluate(std::vector<ComplexNumber> vals, int lo, int hi);
	ComplexNumber evaluateWithCoefficients(std::vector<ComplexNumber> coef, double lo, double hi);
	ComplexNumber evaluateLoToInf(std::vector<ComplexNumber> vals, double lo, double hi);
	void setProductSeries();
	bool testConvergence();
	std::wstring toString();
	void display();
};
#endif