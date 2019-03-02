#pragma once
#include "stdafx.h"
#include "polynomial.h"

class PolynomialFraction {
public:
	Polynomial numerator;
	Polynomial denominator;

	PolynomialFraction() {}
	~PolynomialFraction() {}
	PolynomialFraction(Polynomial p, Polynomial q) {
		numerator = p;
		denominator = q;
	}

	//class methods
	//=============



	std::wstring toString() {return std::wstring(numerator.toString() + L"/" + denominator.toString());}

	void display() { std::wcout << toString() << std::endl; }
};