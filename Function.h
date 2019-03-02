#ifndef FUNCTION_H
#define FUNCTION_H
#pragma once
#include "stdafx.h"

class Function {
public:
	std::vector <wchar_t> variables;	//strictly contains variable chars i.e. it will contain x,y if f(x,y) = x^2 + y^2
	std::wstring function;			//contains just the rhs of the the function (using the example above, it would be 'x^2 + y^2'
	Function();
	~Function();
	Function(std::wstring str);
	Function(std::wstring vars, std::wstring funct);
	Function(std::vector <wchar_t> vars, std::wstring functs);
	Function(std::vector<std::wstring> vec);
	void clear();
	std::vector<wchar_t> varsToChars(std::wstring vars);
	std::vector<wchar_t> combineVariables(std::vector<wchar_t> v1, std::vector<wchar_t> v2);
	bool variablesMatch(Function a, Function b);
	Function abs(Function a);
	Function add(Function a, Function b);
	Function subtract(Function a, Function b);
	Function multiply(Function a, Function b);
	Function divide(Function a, Function b);
	Function inverse(Function a);
	Function exponent(Function a, Function b);
	Function& operator*=(Function& rhs);
	Function& operator/=(Function& rhs);
	Function& operator+=(Function& rhs);
	Function& operator-=(Function& rhs);
	Function& operator*(Function& rhs);
	Function& operator/(Function& rhs);
	Function& operator+(Function& rhs);
	Function& operator-(Function& rhs);
	friend std::wostream& operator<<(std::wostream& os, Function& rhs);
	std::wstring replaceVariable(std::wstring str, wchar_t replace, std::wstring replacer);
	Function compose(Function a, Function b);
	double evalFunction(std::wstring str);
	double evaluate(std::vector<double> vals);
	double evaluate(double x);
	double Function::findARoot(Function f, double lower_bound, double upper_bound);
	double Function::findARoot(Function p, double lower_bound, double upper_bound, std::vector<double> vals, int pos);//for multivariate functions
	std::vector<double> findRoots();
	std::vector<double> findRoots(double lo, double hi);
	std::vector<double> findRoots(double lo, double hi, std::vector<double> vals, int pos);//for multivariate functions
	bool hasDiscontinuities(std::vector<double> testVals);//for real-valued univariate functions
	bool hasDiscontinuities(Matrix testVals);//for real-valued multivariate functions
	std::wstring variablesToString();
	std::wstring toString();
	void display();
};//================================================================

class VectorValuedFunction : public Function {
public:
	std::vector<wchar_t> variables;
	std::vector<Function> f;
	VectorValuedFunction();
	~VectorValuedFunction();
	VectorValuedFunction(std::vector<wchar_t> vars, std::vector<Function> f1);
	VectorValuedFunction(std::vector<Function> f1);
	VectorValuedFunction::VectorValuedFunction(std::wstring f1);
	void clear();
	int size();
	std::vector<double> evaluate(std::vector<double> vals);
	std::vector<double> evaluate(double val);
	std::wstring variablesToString();
	std::wstring toString();
	void display();
};//================================================================

class ComplexFunction : public Function {
public:
	Function f;
	ComplexFunction();
	~ComplexFunction();
	ComplexFunction(std::wstring str);
	ComplexFunction(std::wstring vars, std::wstring funct);
	ComplexFunction(std::vector <wchar_t> vars, std::wstring funct);
	ComplexFunction(Function funct);
	void clear();
	ComplexFunction abs(ComplexFunction a);
	ComplexFunction add(ComplexFunction a, ComplexFunction b);
	ComplexFunction subtract(ComplexFunction a, ComplexFunction b);
	ComplexFunction multiply(ComplexFunction a, ComplexFunction b);
	ComplexFunction divide(ComplexFunction a, ComplexFunction b);
	ComplexFunction exponent(ComplexFunction a, ComplexFunction b);
	ComplexFunction& operator*=(ComplexFunction& rhs);
	ComplexFunction& operator/=(ComplexFunction& rhs);
	ComplexFunction& operator+=(ComplexFunction& rhs);
	ComplexFunction& operator-=(ComplexFunction& rhs);
	ComplexFunction& operator*(ComplexFunction& rhs);
	ComplexFunction& operator/(ComplexFunction& rhs);
	ComplexFunction& operator+(ComplexFunction& rhs);
	ComplexFunction& operator-(ComplexFunction& rhs);
	friend std::wostream& operator<<(std::wostream& os, ComplexFunction& rhs);
	ComplexNumber evalComplexFunction(std::wstring str);
	ComplexNumber evaluateComplex(std::vector<ComplexNumber> vals);
	ComplexNumber evaluateComplex(ComplexNumber z);
	bool hasDiscontinuities(std::vector<ComplexNumber> testVals);//for complex-valued univariate functions
	bool hasDiscontinuities(ComplexMatrix testVals);//for complex-valued multivariate functions
	std::wstring variablesToString();
	std::wstring toString();
	void display();
};//================================================================
class ParametricFunction {
	//a parametric function is a function comprised of many functions all unified with a set of parameters 
public:
	std::vector <wchar_t> variables;	//strictly contains variable chars i.e. it will contain x,y if F(x,y) = x^2i + y^2j
	std::wstring function;			//contains just the rhs of the the function (using the example above, it would be 'x^2i + y^2j'
	ParametricFunction();
	~ParametricFunction();
	ParametricFunction(std::wstring str);
	ParametricFunction(std::vector<wchar_t> params, std::wstring functs);
	ParametricFunction(std::wstring vars, std::wstring funct);
	std::wstring variablesToString();
	std::wstring toString();
	void display();
};//========================================================

	//Helper functions for string parsing:
	int findFunction(std::wstring in);
	bool hasFunction(std::wstring in);
	bool hasOperator(std::wstring in);
	int findComplexNumber(std::wstring str);
	int findNearestDouble(std::wstring str);
	std::wstring getNearestComplexNumber(std::wstring str);
	std::wstring removeComplexNumbers(std::wstring str);
	std::wstring removeDoubles(std::wstring str);
	int countComplexNumbers(std::wstring str);
	int countOperators(std::wstring str);
	bool hasOperatorNotComplexNumber(std::wstring in);
	int findOperator(std::wstring in);
	std::wstring convertDoubleToComplex(std::wstring funct);
	bool isOperableExpressionComplex(std::wstring tempstr);
	bool isOperableExpression(std::wstring tempstr);
	std::wstring processInnerTermsComplex(std::wstring str);
	std::wstring processInnerTerms(std::wstring str);
	std::wstring processTwoTermOperationComplex(std::wstring funct, std::wstring tempstr);
	std::wstring processTwoTermOperation(std::wstring funct, std::wstring tempstr);
	ComplexNumber calculateArithmeticComplex(std::wstring funct);
	double calculateArithmetic(std::wstring funct);
	bool testForFunctions(std::wstring in);
#endif