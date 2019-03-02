/*#pragma once
#include "stdafx.h"
#include "Function.h"

class FractionalFunction {
public:
	Function numerator;
	Function denominator;

	FractionalFunction(Function p, Function q) {
		numerator = p;
		denominator = q;
	}

	//class methods
	//=============
	double evaluate(double x) {
		return (numerator.evaluate(x) / denominator.evaluate(x));
	}
	double evaluate(std::vector<double> vals) {
		
	}

	std::wstring toString() {
		return (numerator.toString() + "/" + denominator.toString());
	}

	void display() { std::wcout << toString() << std::endl; }
};
*/