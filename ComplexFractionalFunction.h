/*#pragma once
#include "stdafx.h"
#include "ComplexFunction.h"

class ComplexFractionalFunction {
public:
	ComplexFunction numerator;
	ComplexFunction denomenator;


	//class methods
	//=============


	std::wstring toString() {
		return (numerator.toString() + L"/" + denomenator.toString());
	}

	void display() { std::wcout << toString() << std::endl; }
};
*/