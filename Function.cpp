#pragma once
#include "stdafx.h"

std::wstring functionList[78] = {//list of all functions recognized by the 'evaluate()' function of the Function class
	L"abs(", L"acos(", L"acosh(", L"arg(", L"asin(", L"asinh(", L"atan(", L"atan2(", L"atanh(", L"BesselFunctionFirstKind(",
	L"BesselFunctionSecondKind(", L"beta", L"cbrt(", L"Ci(", L"conj(", L"ceil(", L"combination(", L"cos(", L"cosh(", 
	L"cot(", L"csc(", L"digamma(", L"doubleFactorial(", L"Ei(", L"erf(", L"erfc", L"exp(", L"exp2(", 
	L"factorial(", L"fallingFactorial(", L"floor(", L"gamma(", L"hypergeometricFunction(", L"hypot(", 
	L"if(", L"imag(", L"int(", L"inverseErfc", L"LambertW(", L"Li(", L"log(", L"log2(", L"log10(", L"logGamma(",
	L"logit(", L"logisticFunction(",
	L"max(", L"min(", L"permutation(", L"prime(", L"sumPrimes(", L"productPrimes(", L"primeCountingFunction(",
	L"polar(", L"pow(", L"polygamma(", L"real(", L"regularModifiedBesselFunctionFirstKind(",
	L"rect(", L"risingFactorial(", L"sec(", L"Si(", L"sin(", L"sinc(",
	L"sinh(", L"sphericalBesselFunctionFirstKind(",
	L"squarewave(", L"sgn(", L"sqrt(", L"step(", L"StirlingNumberFirstKind(", L"StirlingNumberSecondKind(",
	L"tan(", L"tanh(", L"tri(", L"trunc(", L"WrightOmegaFunction(",
	L"zeta("
};

unsigned int functionCount = 78;

Function::Function() {}
Function::~Function() {}

Function::Function(std::wstring str) {		//input of format f(x,y,w,z...) = x*3 + exp(-3*y*z) - ...
		if (str.find(L"=") != std::wstring::npos) {	//case:  string is 'f(x) = 2*x + exp(-y)'
			while (str[0] != '(') { str = str.substr(1); }//clip off everything leading up to the variable definition
			str = str.substr(1);
			for (int i = 0; i < str.find(L"=") + 1; ++i) {
				if (isLetter(str[i]) == true) { variables.push_back(str[i]); }
			}
			function = str.substr(str.find(L"=") + 1);
		}

		else { function = str; }					//case:  string has no variables, i.e. '23-12+exp(3)'
	}

Function::Function(std::wstring vars, std::wstring funct) {
		for (int i = 0; i < vars.size(); ++i) {
			if (isLetter(vars[i]) == true) { variables.push_back(vars[i]); }
		}
		function = funct;
	}

Function::Function(std::vector <wchar_t> vars, std::wstring functs) {
		function = functs;
		variables = vars;
	}

	Function::Function(std::vector<std::wstring> vec) {
		function = vec[0];
		variables = varsToChars(vec[1]);
	}

	//class methods
	//=============
	void Function::clear() {
		variables.clear();
		function.clear();
	}

	std::vector<wchar_t> Function::varsToChars(std::wstring vars) {
		std::vector<wchar_t> variables;
		for (int i = 0; i < vars.size(); ++i) {
			if (isLetter(vars[i]) == true && vars[i] != ',') { variables.push_back(vars[i]); }
		}
		return variables;
	}

	std::vector<wchar_t> Function::combineVariables(std::vector<wchar_t> v1, std::vector<wchar_t> v2) {
		for (int i = 0; i < v2.size(); ++i) {
			bool skip = 0;
			for (int j = 0; j < v1.size(); ++j) {//if any char in v2 is the same as a char in v1, do not append it to v1
				if (v2[i] == v1[j]) {
					skip = 1;
				}
			}
			if (skip == 0) {//if the char is not a duplicate, append it to v1
				v1.push_back(v2[i]);
			}
		}
		return v1;
	}

	bool Function::variablesMatch(Function a, Function b) {
		if (a.variables.size() != b.variables.size()) { return false; }
		for (int i = 0; i < a.variables.size(); ++i) {
			bool findVar = 0;
			for (int j = 0; j < b.variables.size(); ++j) {
				if (a.variables[i] == b.variables[j]) { findVar = 1; }
			}
			if (findVar == 0) { return false; }
		}
		//if all else fails, return true:
		return true;
	}

	//ARITHMETIC METHODS
	//===================
	Function Function::abs(Function a) {
		std::wstring func = L"abs(";
		func.append(a.function);
		func.append(L")");
		return Function(a.variables, func);
	}
	Function Function::add(Function a, Function b) {
		a.function.append(L"+(");
		a.function.append(b.function);
		a.function.append(L")");
		return Function(combineVariables(a.variables, b.variables), a.function);
	}
	Function Function::subtract(Function a, Function b) {
		a.function.append(L"-(");
		a.function.append(b.function);
		a.function.append(L")");
		return Function(combineVariables(a.variables, b.variables), a.function);
	}
	Function Function::multiply(Function a, Function b) {
		a.function.append(L"*(");
		a.function.append(b.function);
		a.function.append(L")");
		return Function(combineVariables(a.variables, b.variables), a.function);
	}
	Function Function::divide(Function a, Function b) {
		a.function.append(L"/(");
		a.function.append(b.function);
		a.function.append(L")");
		return Function(combineVariables(a.variables, b.variables), a.function);
	}

	Function Function::inverse(Function a) {
		std::wstring func = L"1/(";
		func.append(a.function);
		func.append(L")");
		return Function(a.variables, func);
	}
	Function Function::exponent(Function a, Function b)
	{
		a.function.append(L"^(");
		a.function.append(b.function);
		a.function.append(L")");
		return Function(combineVariables(a.variables, b.variables), a.function);
	}
	//OVERLOADED OPERATORS
	//====================
	Function& Function::operator*=(Function& rhs) { *this = multiply(*this, rhs); return *this; }
	Function& Function::operator/=(Function& rhs) { *this = divide(*this, rhs); return *this; }
	Function& Function::operator+=(Function& rhs) { *this = add(*this, rhs); return *this; }
	Function& Function::operator-=(Function& rhs) { *this = subtract(*this, rhs); return *this; }

	Function& Function::operator*(Function& rhs) { return multiply(*this, rhs); }
	Function& Function::operator/(Function& rhs) { return divide(*this, rhs); }
	Function& Function::operator+(Function& rhs) { return add(*this, rhs); }
	Function& Function::operator-(Function& rhs) { return subtract(*this, rhs); }

	std::wostream& operator<<(std::wostream& os, Function& rhs) {
		os << rhs.toString();
		return os;
	}

	std::wstring Function::replaceVariable(std::wstring str, wchar_t replace, std::wstring replacer) {
		for (int i = 0; i < str.length(); ++i) {//run replace loop ONCE only
			std::wstring sub = str.substr(i, str.length());
			int pos = sub.find(replace) + i;
			if (sub.find(replace) == std::wstring::npos || pos < 0 || sub.find(replace) >= str.length() - 1) {
				//break if done replacing
				i = str.length() + 1;
				pos = std::wstring::npos;//set position to NULL for string type
			}
			if (
				(pos > 0 && pos < str.length() && isLetter(sub[sub.find(replace) + 1]) == false && isLetter(sub[sub.find(replace) - 1]) == false) ||
				(pos == 0 && isLetter(sub[sub.find(replace) + 1]) == false) ||
				(pos == str.length() && isLetter(sub[sub.find(replace) - 1]) == false)
				)
			{
				str.erase(pos, 1);
				str.insert(pos, replacer);
				i = pos + replacer.length();
			}
		}
		return str;
	}

	Function Function::compose(Function a, Function b) {
		//variables of both functions must match
		if (variablesMatch(a, b) == false) { return Function(); }
		std::wstring replacer = L"(";
		replacer.append(b.function);
		replacer.append(L")");
		for (int i = 0; i < a.variables.size(); ++i) {
			a.function = replaceVariable(a.function, a.variables[i], replacer);
		}
		return a;
	}

	double Function::evalFunction(std::wstring str) {//convert string to double if the string is a function of the kind that is parsable
		if (str.length() < 2) { return 0; }//all function names have more than 2 letters in a row

		if (str.find(L"abs(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"abs(") + 1));
			return std::abs(tmp);
		}
		if (str.find(L"acos(") != std::wstring::npos && str.find(L"acosh(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"acos(") + 1));
			return acos(tmp);
		}
		if (str.find(L"acosh(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"acosh(") + 1));
			return acosh(tmp);
		}
		/*	if (str.find(L"arg(L") != std::wstring::npos) {
		double tmp = ParseInputNumber(str.substr(str.find(L"arg(L") + 1));
		return acos(tmp);
		}*/
		if (str.find(L"asin(") != std::wstring::npos && str.find(L"asinh(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"asin(") + 1));
			return asin(tmp);
		}
		if (str.find(L"asinh(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"asinh(") + 1));
			return asinh(tmp);
		}
		if (str.find(L"atan(") != std::wstring::npos && str.find(L"atan2(") == std::wstring::npos && str.find(L"atanh(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"atan(") + 1));
			return atan(tmp);
		}
		if (str.find(L"atan2(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return atan2(tmp[1], tmp[0]);
		}
		if (str.find(L"atanh(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"atanh(") + 1));
			return atanh(tmp);
		}
		if (
			str.find(L"BesselFunctionFirstKind(") != std::wstring::npos
			&& str.find(L"sphericalBesselFunctionFirstKind") == std::wstring::npos
			&& str.find(L"regularModifiedBesselFunctionFirstKind(") == std::wstring::npos
			) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return BesselFunction1stKind(tmp[0], tmp[1]);
		}
		if (str.find(L"BesselFunctionSecondKind(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return BesselFunction2ndKind(tmp[0], tmp[1]);
		}
		if (str.find(L"beta(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return betaFunction(tmp);
		}
		/*	if (str.find(L"cbrt(L") != std::wstring::npos) {
		double tmp = ParseInputNumber(str.substr(str.find(L"cbrt(L") + 1));
		return acos(tmp);
		}*/
		if (str.find(L"combination(L") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"combination(L") + 1));
			if (tmp.size() == 2) { return combination(tmp[0], tmp[1]); }
			else { return 0; }
		}
		/*if (str.find(L"conj(L") != std::wstring::npos) {
		double tmp = ParseInputNumber(str.substr(str.find(L"conj(L") + 1));
		return conj(tmp);
		}*/
		if (str.find(L"ceil(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"ceil(") + 1));
			return ceil(tmp);
		}
		if (str.find(L"Ci(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return Ci(tmp);
		}
		if (str.find(L"cos(") != std::wstring::npos && str.find(L"acos(") == std::wstring::npos && str.find(L"cosh(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"cos(") + 1));
			return cos(tmp);
		}
		if (str.find(L"cosh(") != std::wstring::npos && str.find(L"acosh(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"cosh(L") + 1));
			return cosh(tmp);
		}
		if (str.find(L"cot(") != std::wstring::npos && str.find(L"coth(") == std::wstring::npos && str.find(L"acot(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"cot(") + 1));
			return cot(tmp);
		}
		if (str.find(L"csc(") != std::wstring::npos && str.find(L"acsc(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"csc(") + 1));
			return (csc(tmp));
		}
		if (str.find(L"digamma(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return digamma(tmp);
		}
		if (str.find(L"doubleFactorial(") != std::wstring::npos && str.find(L"doubleFactorial(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return doubleFactorial(tmp);
		}
		if (str.find(L"Ei(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return Ei(tmp);
		}
		if (str.find(L"erf(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return errorFunction(tmp);
		}
		if (str.find(L"erfc(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"erfc(") + 1));
			return errorFunctionComplement(tmp);
		}
		if (str.find(L"exp(") != std::wstring::npos && str.find(L"exp2(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"exp(") + 1));
			return exp(tmp);
		}
		if (str.find(L"exp2(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"exp2(") + 1));
			return exp2(tmp);
		}
		if (str.find(L"factorial(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"factorial(") + 1));
			return factorial(tmp);
		}
		if (str.find(L"fallingFactorial(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return fallingFactorial(tmp[0], tmp[1]);
		}
		if (str.find(L"floor(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"floor(") + 1));
			return floor(tmp);
		}
		if (str.find(L"gamma(") != std::wstring::npos &&
			str.find(L"digamma(") == std::wstring::npos &&
			str.find(L"polygamma(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"gamma(") + 1));
			return gamma(tmp);
		}
		if (str.find(L"hypergeometricFunction(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return hypergeometricFunction(tmp);
		}
		if (str.find(L"hypot(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return sqrt(pow(tmp[0], 2) + pow(tmp[1], 2));
		}
		/*if (str.find(L"if(") != std::wstring::npos) {
		double tmp = ParseInputNumber(str.substr(str.find(L"acos(L") + 1));
		return acos(tmp);
		}*/
		/*	if (str.find(L"imag(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"imag(L") + 1));
		return tmp.Im;
		}*/
		if (str.find(L"int(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"int(") + 1));
			return (int)(tmp);
		}
		if (str.find(L"inverseErfc(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"inverseErfc(") + 1));
			return inverseErfc(tmp);
		}
		if (str.find(L"LambertW(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return LambertW(tmp);
		}
		if (str.find(L"Li(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return Li(tmp);
		}
		if (str.find(L"log(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			if (tmp.size() == 1) { return log(tmp[0]); }//standard log base e
			if (tmp.size() == 2) { return logarithm(tmp[0],tmp[1]); }//log w/different base
			return 0;
		}
		if (str.find(L"log2(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"log2(") + 1));
			return log2(tmp);
		}
		if (str.find(L"log10(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"log10(") + 1));
			return log10(tmp);
		}
		if (str.find(L"logGamma(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return logGamma(tmp);
		}
		if (str.find(L"logisticFunction(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"logisticFunction(") + 1));
			return logisticFunction(tmp);
		}
		if (str.find(L"logit(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"logit(") + 1));
			return logit(tmp);
		}
		/*if (str.find(L"max(") != std::wstring::npos) {
		double tmp = ParseInputNumber(str.substr(str.find(L"acos(L") + 1));
		return acos(tmp);
		}
		if (str.find(L"min(") != std::wstring::npos) {
		double tmp = ParseInputNumber(str.substr(str.find(L"acos(L") + 1));
		return acos(tmp);
		}
*/		if (str.find(L"permutation(L") != std::wstring::npos) {
		std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(L") + 1));
		if (tmp.size() == 2) { return permutation(tmp[0], tmp[1]); }
		else { return 0; }
		}
/*		if (str.find(L"polar(") != std::wstring::npos) {
		double tmp = ParseInputNumber(str.substr(str.find(L"acos(L") + 1));
		return acos(tmp);
		}*/
		if (str.find(L"polygamma(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return polygamma(tmp[0], tmp[1]);
		}
		if (str.find(L"pow(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return pow(tmp[0], tmp[1]);
		}
		if (str.find(L"prime(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return prime(tmp);
		}
		if (str.find(L"sumPrimes(") != std::wstring::npos) {
			int tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return sumPrimes(tmp);
		}
		if (str.find(L"productPrimes(") != std::wstring::npos) {
			int tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return productPrimes(tmp);
		}
		if (str.find(L"primeCountingFunction(") != std::wstring::npos) {
			int tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return primeCountingFunction(tmp);
		}
		/*	if (str.find(L"real(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acos(L") + 1));
		return tmp.Re;
		}*/
		if (str.find(L"rect(") != std::wstring::npos) {
			std::wstring s = str.substr(str.find(L"rect(L"));
			s = s.substr(str.find(L"(L") + 1, str.find(L")") - 1);
			double tmp = ParseInputNumber(s);
			return rect(tmp);
		}
		if (str.find(L"regularModifiedBesselFunctionFirstKind(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return regularModifiedBesselFunction1stKind(tmp[0], tmp[1]);
		}
		if (str.find(L"risingFactorial(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return risingFactorial(tmp[0], tmp[1]);
		}
		if (str.find(L"sec(") != std::wstring::npos && str.find(L"asec(") == std::wstring::npos && str.find(L"sech(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"sec(") + 1));
			return sec(tmp);
		}
		if (str.find(L"Si(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return Si(tmp);
		}
		if (str.find(L"sin(") != std::wstring::npos && str.find(L"asin(") == std::wstring::npos && str.find(L"sinh(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"sin(") + 1));
			return sin(tmp);
		}
		if (str.find(L"sinc(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return sinc(tmp);
		}
		if (str.find(L"sinh(") != std::wstring::npos && str.find(L"asinh(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"sinh(") + 1));
			return sinh(tmp);
		}
		if (str.find(L"sgn(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"sgn(") + 1));
			return sgn(tmp);
		}
		if (str.find(L"sphericalBesselFunctionFirstKind(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return sphericalBesselFunction1stKind(tmp[0], tmp[1]);
		}
		if (str.find(L"sqrt(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"sqrt(") + 1));
			return sqrt(tmp);
		}
		if (str.find(L"squarewave(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"squarewave(") + 1));
			return squarewave(tmp);
		}
		if (str.find(L"step(") != std::wstring::npos) {//Heaviside step function
			std::wstring s = str.substr(str.find(L"step("));
			s = s.substr(str.find(L"(") + 1, str.find(L")") - 1);
			double tmp = ParseInputNumber(s);
			return step(tmp);
		}
		if (str.find(L"StirlingNumberFirstKind(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return StirlingNumber1stKind(tmp[0], tmp[1]);
		}
		if (str.find(L"StirlingNumberSecondKind(") != std::wstring::npos) {
			std::vector<double> tmp = ParseNInputNumbers(str.substr(str.find(L"(") + 1));
			return StirlingNumber2ndKind(tmp[0], tmp[1]);
		}
		if (str.find(L"tan(") != std::wstring::npos && str.find(L"tanh(") == std::wstring::npos && str.find(L"atan(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"tan(") + 1));
			return tan(tmp);
		}
		if (str.find(L"tanh(") != std::wstring::npos && str.find(L"atanh(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"tanh(") + 1));
			return tanh(tmp);
		}
		if (str.find(L"tri(") != std::wstring::npos) {
			std::wstring s = str.substr(str.find(L"tri("));
			s = s.substr(str.find(L"(") + 1, str.find(L")") - 1);
			double tmp = ParseInputNumber(s);
			return tri(tmp);
		}
		/*if (str.find(L"trunc(") != std::wstring::npos) {
		double tmp = ParseInputNumber(str.substr(str.find(L"acos(") + 1));
		return acos(tmp);
		}*/
		if (str.find(L"WrightOmegaFunction(") != std::wstring::npos && str.find(L"WrightOmegaFunction(") == std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"WrightOmegaFunction(") + 1));
			return tanh(tmp);
		}
		if (str.find(L"zeta(") != std::wstring::npos) {
			double tmp = ParseInputNumber(str.substr(str.find(L"(") + 1));
			return RiemannZetaFunction(tmp);
		}
		return 0;	//else, if none of the other functions match, return NULL
	}

	int findFunction(std::wstring in) {//test a string to see if it has a function in it
									   //search defined functions
		for (int i = 0; i < functionCount; ++i) {
			if (in.find(functionList[i]) != std::wstring::npos) {
				return i;//return index of function in array of pre-defined functions
			}
		}
		return -1;	//if nothing is thrown, return -1 by default
	}

	bool hasFunction(std::wstring in) {//test a string to see if it has a function in it
									   //defined functions
		for (int i = 0; i < functionCount; ++i) {
			if (in.find(functionList[i]) != std::wstring::npos) {
				return true;
			}
		}
		return false;	//if nothing is thrown, return false by default
	}

	bool hasOperator(std::wstring in) {//test a string to see if it has an operator in it (does not include parenthesis)
		if (
			//defined operators
			in.find(L"+") != std::wstring::npos ||
			in.find(L"-") != std::wstring::npos ||
			in.find(L"*") != std::wstring::npos ||
			in.find(L"/") != std::wstring::npos ||
			in.find(L"^") != std::wstring::npos ||
			in.find(L"%") != std::wstring::npos
			) {
			return true;
		}

		//defined functions
		for (int i = 0; i < functionCount; ++i) {
			if (in.find(functionList[i]) != std::wstring::npos) {
				return true;
			}
		}
		return false;	//if nothing is thrown, return false by default
	}

	int findComplexNumber(std::wstring str) {//return index if ComplexNumber is found, -1 if not
		for (int i = 0; i <= str.length(); ++i) {
			if (i == 0 && str[i] == L'i' && isLetter(str[i + 1]) == false) { return 0; }
			if (i > 0 && str[i] == L'i' && (isNumber(str[i - 1]) == true || str[i - 1] == L'+' || str[i - 1] == L'-')) {
				bool halfway = false;//indicates when you are past the '+' or '-' sign joining the real and imaginary parts of the string
				for (int j = i - 1; j >= 0; --j) {
					if (j == 1 && isNumber(str[j]) == true && isNumber(str[j - 1]) == false && halfway == true) { return 1; }
					if (j == 1 && str[0] == L'-' && halfway == true) { return 0; }
					if (j == 0 && halfway == true) { return 0; }
					if (halfway == true && isNumber(str[j]) == false){
						if (str[j] == L'-') { return j; }
						if (str[j] != L'-' && str[j] != L'.') {
							return j + 1;
						}
					}
					if ((str[j] == L'+' || str[j] == L'-') && halfway == false) { halfway = true; }
				}
			}
		}
		return std::wstring::npos;//if all else fails return null position value
	}

	std::wstring getNearestComplexNumber(std::wstring str) {
		int i = findComplexNumber(str);
		int i2 = str.substr(i).find(L"i") + 1;
		if (str.substr(i).find(L"i") < 0) { return NULL; }
		else {
			std::wstring cn = str.substr(i, i2);
			while (cn[0] != L'-' && isNumber(cn[0]) == false) { cn = cn.substr(1); }
			while (cn[cn.length() - 1] != L'i') { cn.pop_back(); }
			return cn;
		}
		return NULL;
	}

	std::wstring getNearestDouble(std::wstring str) {
		int i = findNearestDouble(str);
		if (str.size() == 1 && i == 0) {
			std::wstring s1 = L" ";
			s1[0] = (str[0]);
			return s1;
		}
		int i2 = i;
		for (int j = i+1; j < str.size(); ++j) {
			if (isNumber(str[j]) == false && str[j] != L'.') {
				i2 = j - 1;
				if (i == 0 && j != 0 && j + 1 < str.size()) { ++i2; }
				j = str.size() + 1;
			}
			if (j == str.size() - 1 && isNumber(str[j])==true) { i2 = str.size(); j += 2; }
		}
		if (i<0 || i2<i) { return NULL; }
		else {
			std::wstring cn = str.substr(i, i2-i+1);
			while (cn[0] != L'-' && isNumber(cn[0]) == false) { cn = cn.substr(1); }
			return cn;
		}
		return NULL;
	}

	std::wstring removeDoubles(std::wstring str) {
		if (findNearestDouble(str) < 0) { return str; }
		while (findNearestDouble(str) >= 0) {
			std::wstring cpnum = getNearestDouble(str);
			str.erase(findNearestDouble(str), cpnum.length());
		}
		return str;
	}

	std::wstring removeComplexNumbers(std::wstring str) {
		if (findComplexNumber(str) < 0) { return str; }
		while (findComplexNumber(str) >= 0) {
			std::wstring cpnum = getNearestComplexNumber(str);
			str.erase(findComplexNumber(str), cpnum.length());
		}
		return str;
	}

	int countDoubles(std::wstring str) {//returns number of ComplexNumbers present in a string
		if (str == L"" || str.size() == 0) { return 0; }//handle null string
		int count = 0;
		while (findNearestDouble(str) >= 0) {
			std::wstring cpnum = getNearestDouble(str);
			size_t end = cpnum.size();
			size_t start = findNearestDouble(str);
			str.erase(findNearestDouble(str), cpnum.size());
			++count;
		}
		return count;
	}

	int countComplexNumbers(std::wstring str) {//returns number of ComplexNumbers present in a string
		if (str == L"" || str.size() == 0) { return 0; }//handle null string
		int count = 0;
		while (findComplexNumber(str) >= 0) {
			std::wstring cpnum = getNearestComplexNumber(str);
			str.erase(findComplexNumber(str), cpnum.size());
			++count;
		}
		return count;
	}

	int countOperators(std::wstring str) {//returns number of operators present in a string
		if (findComplexNumber(str)>0) { str = removeComplexNumbers(str); }
		str = removeDoubles(str);
		if (str == L"" || str.size() == 0) { return 0; }//handle null string
		int count = 0;
		for (int i = 0; i < str.length(); ++i) {
			if (isOperatorNotParen(str[i]) == true) { ++count; }
		}
		return count;
	}

	bool hasOperatorNotComplexNumber(std::wstring in) {//test a string to see if it has an operator in it (does not include parenthesis)
		if (in.size() == 0) { return false; }//handle null string input

		in = removeComplexNumbers(in);

		if (//find defined operators
			in.find(L"*") != std::wstring::npos ||
			in.find(L"/") != std::wstring::npos ||
			in.find(L"^") != std::wstring::npos ||
			in.find(L"%") != std::wstring::npos
			) {
			return true;
		}

		//defined functions
		for (int i = 0; i < functionCount; ++i) {
			if (in.find(functionList[i]) != std::wstring::npos) {
				return true;
			}
		}

		//finally, check to see if there is an operator not attached to a complexNumber
		if (in.find(L"+") != std::wstring::npos) {//check '+' operator, rejecting any complex number of form '1+2i'
			int i = in.find(L'+');
			for (i + 1; i < in.length(); ++i) {
				if (isNumber(in[i]) == false && in[i]!=L',') {
					if (in[i] == L'i') {
						in = in.substr(i);
						i = in.find(L'+');
					}
					else {
						return true;
					}
				}
			}
		}

		if (in.find(L"-") != std::wstring::npos) {
			bool halfway = false;
			int i = in.find(L'-');
			for (i + 1; i < in.length(); ++i) {
				if (halfway == true && i == in.length() - 1) { return true; }
				if (isNumber(in[i]) == false && in[i] != L'.' && in[i]!=L',') {
					if (halfway == false && (in[i] == L'-' || in[i] == L'+')) { halfway = true; }
					if (in[i] == L'i') {
						in = in.substr(i + 1);
						i = in.find(L'+');
						int i2 = in.find(L'-');
						if (i2 >= 0 && i2 < i) { i = i2; }
						halfway = false;
					}
				}
			}
		}

		return false;	//if nothing is thrown, return false by default
	}

	int findOperator(std::wstring in) {//gives position of nearest operator (does not include parenthesis)
		int a = std::wstring::npos;
		if (in.find(L"+") >= 0) { a = ((int)in.find(L"+") < a) ? in.find(L"+") : a; }
		if (in.find(L"-") >= 0) { a = ((int)in.find(L"-") < a) ? in.find(L"-") : a; }
		if (in.find(L"*") >= 0) { a = ((int)in.find(L"*") < a) ? in.find(L"*") : a; }
		if (in.find(L"/") >= 0) { a = ((int)in.find(L"/") < a) ? in.find(L"/") : a; }
		if (in.find(L"^") >= 0) { a = ((int)in.find(L"^") < a) ? in.find(L"^") : a; }
		if (in.find(L"%") >= 0) { a = ((int)in.find(L"%") < a) ? in.find(L"%") : a; }
		return a;
	}

	int findNearestDouble(std::wstring str) {
		int start = 0;
		for (int i = 0; i < str.length(); ++i) {
			if (isNumber(str[i]) == true) {
				bool isNegativeNumber = false;
				if (i > 0 && str[i - 1] == L'-') {
					int j = i - 1;
					while (j >= 0 && str[j] != L'(' && isNumber(str[j]) == false) {
						if ( (j>0 && str[j - 1] == L'(') || j==0 ) { isNegativeNumber = true; j = -1; }
						--j;
					}
				}
				start = i;
				if (isNegativeNumber == true) { start = i - 1; }
				return start;
			}
		}
		return std::wstring::npos;
	}

	std::wstring convertDoubleToComplexNumber(std::wstring str, int start, int end) {
		std::wstring replacer = L"(";
		replacer.append(ComplexNumber(ParseInputNumber(str.substr(start, end - start))).toStringBothParts());
		replacer.append(L")");
		str.erase(start, end - start);
		str.insert(start, replacer);
		return str;
	}

	std::wstring convertDoubleToComplex(std::wstring funct) {
		//first, convert any double values in the string to ComplexNumbers
		bool foundNum = false;
		int start = 0;
		int end = 0;
		for (int i = 0; i < funct.length(); ++i) {
			if (funct[i] == L'i') {	//if it is a complex number, skip
				start = end = i;
				foundNum = false;
			}
			if (i == funct.length() - 1 && foundNum == true && isNumber(funct[i])==true) {
				end = i;
				funct = convertDoubleToComplexNumber(funct, start, end);
				start = end = i;
				foundNum = false;
			}
			if (isNumber(funct[i]) == false && start != end && funct[i] != L'.' && findComplexNumber(funct.substr(i - 1)) != 0) {//find end of double
				if (foundNum == true) {
					end = i;
					funct = convertDoubleToComplexNumber(funct, start, end);
				}
				start = end = i;
				foundNum = false;
			}
			if (isNumber(funct[i]) == true) {//if the char is a number, include it in the range
				bool isNegativeNumber = false;
				if (i > 0 && funct[i - 1] == L'-') {
					int j = i-1;
					while (j > 0 && funct[j] != L'(' && isNumber(funct[j]) == false) {
						if (funct[j - 1] == L'(') { isNegativeNumber = true; j = -1; }
						--j;
					}
				}

				if (start == end) {					
					start = i;
					if (isNegativeNumber == true) {
						start = i - 1;
					}
					foundNum = true;
				}
				end = i + 1;
				if (i == funct.length() - 1 && foundNum == true) { 
					funct = convertDoubleToComplexNumber(funct, start, end); 
					foundNum = false;
					i += end - start;
				}
			}
		}
		return funct;
	}

	std::wstring processTwoTermOperation(std::wstring funct, std::wstring tempstr) {
		double x = calculateArithmetic(tempstr);
		std::wstring replc = (L"(");
		replc.append(to_stringPrecision(x));
		replc.append(L")");
		funct = replaceString(funct, tempstr, replc);
		return funct;
	}

	std::wstring processTwoTermOperationComplex(std::wstring funct, std::wstring tempstr) {
		ComplexNumber x = calculateArithmeticComplex(tempstr);
		std::wstring replc = (L"(");
		replc.append(x.toStringBothParts());
		replc.append(L")");
		funct = replaceString(funct, tempstr, replc);
		return funct;
	}

	bool isOperableExpression(std::wstring tempstr) {
		if (countDoubles(tempstr) == 2 && countOperators(tempstr) == 1) { return true; }
		return false;
	}

	bool isOperableExpressionComplex(std::wstring tempstr) {
		if (countComplexNumbers(tempstr) == 2 && countOperators(tempstr) == 1){return true;}
		return false;
	}

	std::wstring processInnerTerms(std::wstring str) {
		if (str.find(L'(') == std::wstring::npos) { return str; }
		bool hasOperableTerms = true;
		while (hasOperableTerms == true) {
			std::vector<std::wstring> exps = loadExpressionsIntoArray(str);
			for (int i = 0; i < exps.size(); ++i) {
				if (countDoubles(exps[i]) >=2 && countOperators(exps[i])>=1) {
					str = processTwoTermOperation(str, exps[i]);
					i = exps.size() + 1;
				}
				if (i == exps.size() - 1) { hasOperableTerms = false; }
			}
		}
		return str;
	}

	std::wstring processInnerTermsComplex(std::wstring str) {
		bool hasOperableTerms = true;
		while (hasOperableTerms == true) {
			std::vector<std::wstring> exps = loadExpressionsIntoArray(str);
			for (int i = 0; i < exps.size(); ++i) {
				if (countComplexNumbers(exps[i]) >= 2 && countOperators(exps[i]) >= 1) {
					str = processTwoTermOperationComplex(str, exps[i]);
					i = exps.size() + 1;
				}
				if (i == exps.size() - 1) { hasOperableTerms = false; }
			}
		}
		return str;
	}

	double Function::evaluate(std::vector<double> vals) {
		/*	NOTE:
		Suitable arithmetic operators are +,-,*,/,^.  Suitable functions are found in evalFunction method.
		This function will properly evaluate a function of the form f(x,y) = x^2 without error.		*/

		/*=================
		|DECLARE VARIABLES|
		=================*/
		//make sure when using to_string precision or other functions which truncate values, we do as little truncation as possible
		int prcs = precision;
		precision = 30;

		std::wstring funct = function;//local copy of function string, necessary as we will erase/replace parts of the string as we go

		funct = removeSpaces(funct);//remove all spaces in function string
		funct = replaceConstants(funct);//replace all system constants in string
		//========================================
		
		/*================================
		|REPLACE VARIABLES WITH CONSTANTS|
		================================*/
		//replace all variables with their double value from double* vals
		if (vals.size() > 0) {
			for (int i = 0; i < variables.size(); ++i) {
				std::wstring val = to_stringPrecision(vals[i]);//double values for variable made into a string

			   //check string for spots to replace char, then perform replace function
				for (int j = 0; j < funct.size(); ++j) {
					if (j == 0 && funct[j] == variables[i] && isLetter(funct[j + 1]) == false) {
						//std::wcout << L"replacing char #" << j << L"...\n";//debug
						funct.erase(j, 1);
						funct.insert(j, val);
						j = 0;
					}
					if (j>0 && j<funct.size() - 1 && funct[j] == variables[i] && isLetter(funct[j + 1]) == false && isLetter(funct[j - 1]) == false) {
						//std::wcout << L"replacing char #" << j << L"...\n";//debug
						funct.erase(j, 1);
						funct.insert(j, val);
					}
					if (j == funct.size() - 1 && funct[j] == variables[i] && isLetter(funct[j - 1]) == false) {
						//std::wcout << L"replacing char #" << j << L"...\n";//debug
						funct.erase(j, 1);
						funct.insert(j, val);
					}
				}
			}
		}//================================================================

		 //remove any double '-' signs that might have resulted accidentally
		for (int i = 0; i < funct.size() - 1; ++i) {
			if (funct[i] == L'-' && funct[i + 1] == L'-') {
				funct[i] = L' ';
				funct[i + 1] = L'+';
				++i;
			}
		}
		funct = removeSpaces(funct);

		/*================
		|HANDLE FUNCTIONS|
		================*/
		bool searchForFuncts = hasFunction(funct);
		while (searchForFuncts == true) {
			bool foundOne = false;
			for (int i = 0; i < funct.size() - 1; ++i) {
			topOfTheFunctsLoop:
				if (isLetter(funct[i]) == true && isLetter(funct[i + 1]) == true) {		//find functions, which have more than one letter in a row
					std::wstring tempstr = funct.substr(i);
					tempstr = tempstr.substr(0, findMatchingParen(tempstr, tempstr.find(L"("))+1);		//save function (i.e. 'exp(12.2)') as a substr to evaluate
					if (hasFunction(tempstr.substr(tempstr.find(L"("))) == true) {//keep searching for the innermost paren without a function inside
						i += tempstr.find(L"(")-1;//skip past the function (i.e. i+=4 if outer function is 'cos(' for example)
						goto topOfTheFunctsLoop;
					}
					std::wstring toReplace = tempstr;

					//make sure there's nothing within the parenthesis that needs to be simplified (i.e. 'sin(12-2)')
					int posit = tempstr.find(L'(');
					if (hasOperator(tempstr.substr(posit+1, findMatchingParen(tempstr,posit)+ 1)) == true && countDoubles(tempstr)>1) {
						int opers = 0;
						for (int k = 0; k < tempstr.length() - 2; ++k) {
							if (isOperator(tempstr[k]) == true) {
								if (tempstr[k] == '-') {
									if (k > 0) {
										if (isNumber(tempstr[k - 1]) == true && isNumber(tempstr[k + 1]) == true) {
											++opers;
										}
									}
								}
								else { ++opers; }
							}
						}

						if (opers>0) {
							size_t a = tempstr.find(L"(") + 1;
							size_t b = tempstr.find(L")") - 1;
							std::wstring chk = tempstr.substr(a, b);
							chk = replaceChar(chk, L')', L' ');
							chk = replaceChar(chk, L'(', L' ');
							chk = removeSpaces(chk);
							Function f2(chk);
							double p = f2.evaluate(0.0);
							std::wstring replc = to_stringPrecision(p);
							tempstr = replaceString(tempstr, chk, replc);
						}
					}

					//now evaluate the function, get the double value, and replace the whole thing in the original string (i.e. 'sin(0)' --> '0')
					double tmpdbl = evalFunction(tempstr);
					funct = replaceString(funct, toReplace, to_stringPrecision(tmpdbl));
					foundOne = true;
				}
			}
			if (foundOne == false) { searchForFuncts = false; }
		}//====================================================
		if (countDoubles(funct) == 1) { 
			precision = prcs;
			return stringToDouble(getNearestDouble(funct)); 
		}
		funct = processInnerTerms(funct); //handle parens
		double answer = calculateArithmetic(funct);//calculate arithmetic of the last few values
		precision = prcs; //finally, restore precision and return a value
		return answer; //note: if unsuccessful, will return 0
	}//END EVALUATE FUNCTION================================================

	double Function::evaluate(double x) {
		std::vector<double> vals;
		vals.push_back(x);
		double n = evaluate(vals);
		return n;
	}

	double Function::findARoot(Function p, double lower_bound, double upper_bound, std::vector<double> vals, int pos) {
		//else, p is of order 3 or higher....
		double TOL = 0.0001;    // tolerance
		double MAX_ITER = 1000; // maximum number of iterations
		double a = lower_bound;
		double b = upper_bound;
		vals.insert(vals.begin()+pos, a);
		double fa = p.evaluate(vals);   // calculated now to save function calls
		vals.erase(vals.begin() + pos);
		vals.insert(vals.begin() + pos, b);
		double fb = p.evaluate(vals);   // calculated now to save function calls
		vals.erase(vals.begin() + pos);
		double fs = 0;      // initialize 

		if (!(fa * fb < 0)) { return 0; }

		if (fabsf(fa) < fabsf(b)) { // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
			std::swap(a, b);
			std::swap(fa, fb);
		}

		double c = a;           // c now equals the largest magnitude of the lower and upper bounds
		double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
		bool mflag = true;      // boolean flag used to evaluate if statement later on
		double s = 0;           // Our Root that will be returned
		double d = 0;           // Only used if mflag is unset (mflag == false)

		for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
		{
			// stop if converged on root or error is less than tolerance
			if (std::abs(b - a) < TOL) {
				goto NarrowResult;
			} // end if

			if (fa != fc && fb != fc) {
				// use inverse quadratic interopolation
				s = (a * fb * fc / ((fa - fb) * (fa - fc)))
					+ (b * fa * fc / ((fb - fa) * (fb - fc)))
					+ (c * fa * fb / ((fc - fa) * (fc - fb)));
			}
			else { s = b - fb * (b - a) / (fb - fa); }//secant method

			 /*
			Crazy condition statement!:
			 --------------------------------------------------
			(condition 1) s is not between  (3a+b)/4  and b or
			(condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
			(condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
			(condition 4) (mflag is set and |b−c| < |TOL|) or
			(condition 5) (mflag is false and |c−d| < |TOL|)
			*/
			if (((s < (3 * a + b) * 0.25) || (s > b)) ||
				(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
				(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
				(mflag && (std::abs(b - c) < TOL)) ||
				(!mflag && (std::abs(c - d) < TOL)))
			{
				// bisection method
				s = (a + b)*0.5;
				mflag = true;
			}
			else { mflag = false; }

			vals.insert(vals.begin() + pos, s);
			fs = p.evaluate(vals);  // calculate fs
			vals.erase(vals.begin()+pos);
			d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
			c = b;      // set c equal to upper bound
			fc = fb;    // set f(c) = f(b)

			if (fa * fs < 0) {// fa and fs have opposite signs
				b = s;
				fb = fs;    // set f(b) = f(s)
			}
			else {
				a = s;
				fa = fs;    // set f(a) = f(s)
			}

			if (std::abs(fa) < std::abs(fb)) {// if magnitude of fa is less than magnitude of fb
				std::swap(a, b);     // swap a and b
				std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
			}
		}
		return 0;

	NarrowResult:
		//use root as guess for Newton's method now
		double x, x1; //user input value
		double f, dfdx; //function and its derivative 
		const double epsil = 0.0000001; //tolerance
		const int n = 15; //# steps

		x1 = s;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			vals.insert(vals.begin() + pos, x);
			f = p.evaluate(vals); //Defines the function 
			dfdx = nthPartialDerivative(1, pos, vals, p); //derivative of Function
			vals.erase(vals.begin()+pos);
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (std::abs(x - x1) < epsil) { //checks if guesses are within a certain tolerance value 
			return x1;
		}
		return 0;
	}

	double Function::findARoot(Function p, double lower_bound, double upper_bound) {
		//else, p is of order 3 or higher....
		double TOL = 0.0001;    // tolerance
		double MAX_ITER = 1000; // maximum number of iterations
		double a = lower_bound;
		double b = upper_bound;
		double fa = p.evaluate(a);   // calculated now to save function calls
		double fb = p.evaluate(b);   // calculated now to save function calls
		double fs = 0;      // initialize 

		if (!(fa * fb < 0)) { return 0; }

		if (fabsf(fa) < fabsf(b)) { // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
			std::swap(a, b);
			std::swap(fa, fb);
		}

		double c = a;           // c now equals the largest magnitude of the lower and upper bounds
		double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
		bool mflag = true;      // boolean flag used to evaluate if statement later on
		double s = 0;           // Our Root that will be returned
		double d = 0;           // Only used if mflag is unset (mflag == false)

		for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
		{
			// stop if converged on root or error is less than tolerance
			if (std::abs(b - a) < TOL) {
				goto NarrowResult;
			} // end if

			if (fa != fc && fb != fc) {
				// use inverse quadratic interopolation
				s = (a * fb * fc / ((fa - fb) * (fa - fc)))
					+ (b * fa * fc / ((fb - fa) * (fb - fc)))
					+ (c * fa * fb / ((fc - fa) * (fc - fb)));
			}
			else { s = b - fb * (b - a) / (fb - fa); }//secant method

													  /*
													  Crazy condition statement!:
													  --------------------------------------------------
													  (condition 1) s is not between  (3a+b)/4  and b or
													  (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
													  (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
													  (condition 4) (mflag is set and |b−c| < |TOL|) or
													  (condition 5) (mflag is false and |c−d| < |TOL|)
													  */
			if (((s < (3 * a + b) * 0.25) || (s > b)) ||
				(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
				(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
				(mflag && (std::abs(b - c) < TOL)) ||
				(!mflag && (std::abs(c - d) < TOL)))
			{
				// bisection method
				s = (a + b)*0.5;
				mflag = true;
			}
			else { mflag = false; }

			fs = p.evaluate(s);  // calculate fs
			d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
			c = b;      // set c equal to upper bound
			fc = fb;    // set f(c) = f(b)

			if (fa * fs < 0) {// fa and fs have opposite signs
				b = s;
				fb = fs;    // set f(b) = f(s)
			}
			else {
				a = s;
				fa = fs;    // set f(a) = f(s)
			}

			if (std::abs(fa) < std::abs(fb)) {// if magnitude of fa is less than magnitude of fb
				std::swap(a, b);     // swap a and b
				std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
			}
		}
		return 0;

	NarrowResult:
		//use root as guess for Newton's method now
		double x, x1; //user input value
		double f, dfdx; //function and its derivative 
		const double epsil = 0.0000001; //tolerance
		const int n = 15; //# steps

		x1 = s;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.evaluate(x); //Defines the function 
			dfdx = nthDerivative1D(1, x, p); //derivative of Function
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (std::abs(x - x1) < epsil) { //checks if guesses are within a certain tolerance value 
			return x1;
		}
		return 0;
	}

	std::vector<double> Function::findRoots() { return findRoots(-20,20); }

	std::vector<double> Function::findRoots(double lo, double hi) {
		std::vector<double> candidates;
		double lowbound = lo;
		double highbound = hi;

		//Since 0 is used as a null value for roots in the rootfinding algorithms,
		//we must determine whether or not it is a root outside of these algorithms.
		if (evaluate(0) == 0) { candidates.push_back(0); }

		//now sweep through the bounds of the polynomial to look for roots using the
		//mixed Brent/Newton search method in findRoot()
		double lastroot = 0;
		double delta = 0.01;
		while (lowbound < highbound) {
			double x = findARoot(*this, lowbound, lowbound + delta);
			if (x == 0) { lowbound += delta; }//using 0 as a null return value here
			else {
				if (x != lastroot) { //prevent same root from continually showing up
					candidates.push_back(x);
					lastroot = x;
					lowbound = x + delta;
					if (x == floor(x)) {
						//temp = temp.factorOutBinomial(x); //perform a shift to reduce the polynomial
						//bounds = temp.getRootBounds();
						//double lowbound = bounds[0];
						//double highbound = bounds[1];
					}//if x is an integer, you can factor it out of the polynomial
				}
				else { lowbound += delta; }
			}
		}
		return candidates;
	}

	std::vector<double> Function::findRoots(double lo, double hi, std::vector<double> vals, int pos) {
		std::vector<double> candidates;
		double lowbound = lo;
		double highbound = hi;

		//Since 0 is used as a null value for roots in the rootfinding algorithms,
		//we must determine whether or not it is a root outside of these algorithms.
		vals.insert(vals.begin() + pos, 0);
		double vl = evaluate(vals);
		if (vl == 0) { candidates.push_back(0); }
		vals.erase(vals.begin() + pos);

		//now sweep through the bounds of the polynomial to look for roots using the
		//mixed Brent/Newton search method in findRoot()
		double lastroot = 0;
		double delta = 0.01;
		while (lowbound < highbound) {
			double x = findARoot(*this, lowbound, lowbound + delta,vals,pos);
			if (x == 0) { lowbound += delta; }//using 0 as a null return value here
			else {
				if (x != lastroot) { //prevent same root from continually showing up
					candidates.push_back(x);
					lastroot = x;
					lowbound = x + delta;
					if (x == floor(x)) {
						//temp = temp.factorOutBinomial(x); //perform a shift to reduce the polynomial
						//bounds = temp.getRootBounds();
						//double lowbound = bounds[0];
						//double highbound = bounds[1];
					}//if x is an integer, you can factor it out of the polynomial
				}
				else { lowbound += delta; }
			}
		}
		return candidates;
	}

	bool Function::hasDiscontinuities(std::vector<double> testVals) {//looks for possible discontinuities by examining the lh derivative vs the rh deriviative
		//for univariate functions
		double epsil = 0.00001;
		for (int i = 0; i < testVals.size(); ++i) {
			double left = evaluate(testVals[i]+epsil);
			double right = evaluate(testVals[i]-epsil);
			if ((0.5*(left - right)) > epsil) { return true; }
		}
		return false;
	}

	bool Function::hasDiscontinuities(Matrix testVals) {//for real-valued multivariate functions
		double epsil = 0.00001;
		for (int i = 0; i < testVals.size(); ++i) {
			double left = evaluate(testVals.element[i] + epsil);
			double right = evaluate(testVals.element[i] - epsil);
			if ((0.5*(left - right)) > epsil) { return true; }
		}
		return false;
	}

	std::wstring Function::variablesToString() { //returns string of variables in format 'x,y,z,w,p,...'
		std::wstring vars = L"";
		if (variables.size() == 1) { vars.push_back(variables[0]); return vars; }
		if (variables.size() > 1) {
			for (int i = 0; i < variables.size() - 1; ++i) {
				vars.push_back(variables[i]);
				vars.push_back(',');
			}
			vars.push_back(variables[variables.size() - 1]);
			return vars;
		}
		return vars;
	}

	std::wstring Function::toString() { return function; }	//returns function in format 'x * exp(y+2)...'
	void Function::display() { std::wcout << toString() << L"\n"; }
//================================================================	   	 

	VectorValuedFunction::VectorValuedFunction() {}
	VectorValuedFunction::~VectorValuedFunction() {}

	VectorValuedFunction::VectorValuedFunction(std::vector<wchar_t> vars, std::vector<Function> f1) {
		sort(vars.begin(), vars.end());
		vars.erase(unique(vars.begin(), vars.end()), vars.end());
		variables = vars;
		f = f1;
	}

	VectorValuedFunction::VectorValuedFunction(std::wstring f1) {
		VectorValuedFunction vvf(ParseNFunctions(f1));
		variables = vvf.variables;
		f = vvf.f;
	}

	VectorValuedFunction::VectorValuedFunction(std::vector<Function> f1) {
		std::vector<wchar_t> varbs;
		f = f1;
		for (int i = 0; i < f1.size(); ++i) {
			varbs = f1[i].combineVariables(f1[i].variables, varbs);
		}
		sort(varbs.begin(), varbs.end());
		varbs.erase(unique(varbs.begin(), varbs.end()), varbs.end());
		variables = varbs;
	}

	//class methods
	//=============
	void VectorValuedFunction::clear() { f.clear(); variables.clear(); }
	int VectorValuedFunction::size() { return f.size(); }

	std::vector<double> VectorValuedFunction::evaluate(std::vector<double> vals) {
		std::vector<double> answer;
		for (int i = 0; i < size(); ++i) {
			answer.push_back(f[i].evaluate(vals));
		}
		return answer;
	}

	std::vector<double> VectorValuedFunction::evaluate(double val) {
		std::vector<double> vals;
		vals.push_back(val);
		return evaluate(vals);
	}

	std::wstring VectorValuedFunction::variablesToString() { //returns string of variables in format 'x,y,z,w,p,...'
		std::wstring vars = L"";
		if (variables.size() == 1) { vars.push_back(variables[0]); return vars; }
		if (variables.size() > 1) {
			for (int i = 0; i < variables.size() - 1; ++i) {
				vars.push_back(variables[i]);
				vars.push_back(',');
			}
			vars.push_back(variables[variables.size() - 1]);			
			return vars;
		}
		return vars;
	}

	std::wstring VectorValuedFunction::toString() {
		std::wstring str = L"<";
		if (size() > 0) {
			for (int i = 0; i < size() - 1; ++i) {
				str.append(f[i].toString());
				str.append(L",");
			}
		}
		str.append(f[f.size()-1].toString());
		str.append(L">");
		return str;
	}

	void VectorValuedFunction::display() { std::wcout << toString() << L"\n"; }
//================================================================

	ComplexFunction::ComplexFunction() {}
	ComplexFunction::~ComplexFunction() {}

	ComplexFunction::ComplexFunction(std::wstring str) {		//input of format f(x,y,w,z...) = x*3 + exp(-3*y*z) - ...
		std::wstring function;
		std::vector<wchar_t> variables;
		if (str.find(L"=") != std::wstring::npos) {	//case:  string is 'f(x) = 2*x + exp(-y)'
			while (str[0] != '(') { str = str.substr(1); }//clip off everything leading up to the variable definition
			str = str.substr(1);
			for (int i = 0; i < str.find(L"=") + 1; ++i) {
				if (isLetter(str[i]) == true) { variables.push_back(str[i]); }
			}
			function = str.substr(str.find(L"=") + 1);
		}
		else { function = str; }					//case:  string has no variables, i.e. '23-12+exp(3)'
		f = Function(variables, function);
	}

	ComplexFunction::ComplexFunction(std::wstring vars, std::wstring funct) {
		f = Function(vars, funct);
	}

	ComplexFunction::ComplexFunction(std::vector <wchar_t> vars, std::wstring funct) {
		f = Function(vars, funct);
	}

	ComplexFunction::ComplexFunction(Function funct) { f = funct; }

	void ComplexFunction::clear() { f.clear(); }

	//ARITHMETIC METHODS
	//===================
	ComplexFunction ComplexFunction::abs(ComplexFunction a) { return ComplexFunction(f.abs(a.f)); }
	ComplexFunction ComplexFunction::add(ComplexFunction a, ComplexFunction b) { return ComplexFunction(a.add(a.f, b.f)); }
	ComplexFunction ComplexFunction::subtract(ComplexFunction a, ComplexFunction b) { return ComplexFunction(a.subtract(a.f, b.f)); }
	ComplexFunction ComplexFunction::multiply(ComplexFunction a, ComplexFunction b) { return ComplexFunction(a.multiply(a.f, b.f)); }
	ComplexFunction ComplexFunction::divide(ComplexFunction a, ComplexFunction b) { return ComplexFunction(a.divide(a.f, b.f)); }
	ComplexFunction ComplexFunction::exponent(ComplexFunction a, ComplexFunction b) { return ComplexFunction(a.exponent(a.f, b.f)); }

	//OVERLOADED OPERATORS
	//====================
	ComplexFunction& ComplexFunction::operator*=(ComplexFunction& rhs) { *this = multiply(*this, rhs); return *this; }
	ComplexFunction& ComplexFunction::operator/=(ComplexFunction& rhs) { *this = divide(*this, rhs); return *this; }
	ComplexFunction& ComplexFunction::operator+=(ComplexFunction& rhs) { *this = add(*this, rhs); return *this; }
	ComplexFunction& ComplexFunction::operator-=(ComplexFunction& rhs) { *this = subtract(*this, rhs); return *this; }

	ComplexFunction& ComplexFunction::operator*(ComplexFunction& rhs) { return multiply(*this, rhs); }
	ComplexFunction& ComplexFunction::operator/(ComplexFunction& rhs) { return divide(*this, rhs); }
	ComplexFunction& ComplexFunction::operator+(ComplexFunction& rhs) { return add(*this, rhs); }
	ComplexFunction& ComplexFunction::operator-(ComplexFunction& rhs) { return subtract(*this, rhs); }

	std::wostream& operator<<(std::wostream& os, ComplexFunction& rhs) {
		os << rhs.toString();
		return os;
	}

	ComplexNumber ComplexFunction::evalComplexFunction(std::wstring str) {//convert string to double if the string is a function of the kind that is parsable
		if (str.length() < 2) { return ComplexNumber(0); }//all function names have more than 2 letters in a row

		if (str.find(L"abs(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"abs(") + 1));
			return tmp.modulus();
		}
		if (str.find(L"acos(") != std::wstring::npos && str.find(L"acosh(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acos(") + 1));
			return acos(tmp);
		}
		if (str.find(L"acosh(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acosh(") + 1));
			return acosh(tmp);
		}
		if (str.find(L"arg(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return tmp.arg();
		}
		if (str.find(L"asin(") != std::wstring::npos && str.find(L"asinh(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"asin(") + 1));
			return asin(tmp);
		}
		if (str.find(L"asinh(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"asinh(") + 1));
			return asinh(tmp);
		}
		if (str.find(L"atan(") != std::wstring::npos && str.find(L"atan2(") == std::wstring::npos && str.find(L"atanh(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"atan(") + 1));
			return atan(tmp);
		}
		/*		if (str.find(L"atan2(") != std::wstring::npos) {
		std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"atan2(") + 1));
		return atan2(tmp[1], tmp[0]);
		}
		*/		if (str.find(L"atanh(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"atanh(") + 1));
			return atanh(tmp);
		}
		if (str.find(L"BesselFunctionFirstKind(") != std::wstring::npos
			&& str.find(L"sphericalBesselFunctionFirstKind(") == std::wstring::npos
			&& str.find(L"regularModifiedBesselFunctionFirstKind") == std::wstring::npos
			) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"(") + 1));
			if (tmp[0].Im==0) {
				return BesselFunction1stKind(tmp[0].Re,tmp[1]);
			}
			else {
				return BesselFunction1stKind(tmp[0], tmp[1]);
			}
		}
		if (str.find(L"BesselFunctionSecondKind(") != std::wstring::npos	) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"(") + 1));
			if (tmp[0].Im == 0) {
				return BesselFunction2ndKind(tmp[0].Re, tmp[1]);
			}
			else {
				return BesselFunction2ndKind(tmp[0], tmp[1]);
			}
		}
		if (str.find(L"beta(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"beta(") + 1));
			return betaFunction(tmp);
		}
		/*	if (str.find(L"cbrt(L") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"cbrt(L") + 1));
		return acos(tmp);
		}*/
		if (str.find(L"conj(L") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"conj(L") + 1));
			return tmp.conjugate();
		}
		if (str.find(L"ceil(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"ceil(") + 1));
			return ceil(tmp);
		}
		if (str.find(L"Ci(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return Ci(tmp);
		}
		if (str.find(L"cos(") != std::wstring::npos && str.find(L"acos(") == std::wstring::npos && str.find(L"cosh(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"cos(") + 1));
			return cos(tmp);
		}
		if (str.find(L"cosh(") != std::wstring::npos && str.find(L"acosh(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"cosh(L") + 1));
			return cosh(tmp);
		}
		if (str.find(L"cot(") != std::wstring::npos && str.find(L"coth(") == std::wstring::npos && str.find(L"acot(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"cot(") + 1));
			return cot(tmp);
		}
		if (str.find(L"csc(") != std::wstring::npos && str.find(L"acsc(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"csc(") + 1));
			return (csc(tmp));
		}
		if (str.find(L"digamma(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return digamma(tmp);
		}
/*		if (str.find(L"doubleFactorial(") != std::wstring::npos && str.find(L"doubleFactorial(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return doubleFactorial(tmp);
		}
*/		if (str.find(L"Ei(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return Ei(tmp);
		}
/*		if (str.find(L"erf(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"erf(") + 1));
		return erf(tmp);
		}
		*//*		if (str.find(L"erfc(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"erfc(") + 1));
		return erfc(tmp);
		}
		*/if (str.find(L"exp(") != std::wstring::npos && str.find(L"exp2(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"exp(") + 1));
			return exp(tmp);
		}
		/*		if (str.find(L"exp2(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"exp2(") + 1));
		return exp2(tmp);
		}*/
/*		if (str.find(L"factorial(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"factorial(") + 1));
		return factorial(tmp);
		}
*/		if (str.find(L"fallingFactorial(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"(") + 1));
			return fallingFactorial(tmp[0], tmp[1].Re);
		}
		if (str.find(L"floor(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"floor(") + 1));
			return floor(tmp);
		}
		if (str.find(L"gamma(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"gamma(") + 1));
			return gamma(tmp);
		}
		if (str.find(L"hypergeometricFunction(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"(") + 1));
			return hypergeometricFunction(tmp);
		}
		if (str.find(L"hypot(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"hypot(") + 1));
			return sqrt(pow(tmp[0], 2) + pow(tmp[1], 2));
		}
		/*if (str.find(L"if(L") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acos(L") + 1));
		return acos(tmp);
		}*/
		if (str.find(L"imag(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"imag(L") + 1));
			return tmp.Im;
		}
		/*		if (str.find(L"int(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"int(") + 1));
		return (int)(tmp);
		}
*//*		if (str.find(L"inverseErfc(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"inverseErfc(") + 1));
			return inverseErfc(tmp);
		}
*/		if (str.find(L"LambertW(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"LambertW(") + 1));
			return LambertW(tmp);
		}
		if (str.find(L"Li(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return Li(tmp);
		}
		if (str.find(L"log(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"log(") + 1));
			return log(tmp);
		}
		if (str.find(L"log2(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"log2(") + 1));
			return log2(tmp);
		}
		if (str.find(L"log10(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"log10(") + 1));
			return log10(tmp);
		}
		if (str.find(L"logGamma(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return logGamma(tmp);
		}
/*		if (str.find(L"logit(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"logit(") + 1));
			return logit(tmp);
		}
*//*	if (str.find(L"logisticFunction(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"logisticFunction(") + 1));
			return logisticFunction(tmp);
		}
*//*	if (str.find(L"max(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acos(L") + 1));
		return acos(tmp);
		}
		if (str.find(L"min(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acos(L") + 1));
		return acos(tmp);
		}
		*/		if (str.find(L"polar(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acos(L") + 1));
			return tmp.polar();
		}
		if (str.find(L"pow(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"pow(") + 1));
			return pow(tmp[0], tmp[1]);
		}
		if (str.find(L"polygamma(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"(") + 1));
			return polygamma(tmp[0], tmp[1].Re);
		}
/*		if (str.find(L"prime(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return prime(tmp);
		}
*/		if (str.find(L"real(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acos(L") + 1));
			return tmp.Re;
		}
		/*		if (str.find(L"rect(") != std::wstring::npos) {
		std::wstring s = str.substr(str.find(L"rect(L"));
		s = s.substr(str.find(L"(L") + 1, str.find(L")") - 1);
		ComplexNumber tmp = ParseComplexNumber(s);
		return rect(tmp);
		}
		*/		
		if (str.find(L"regularModifiedBesselFunctionFirstKind(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"(") + 1));
			if (tmp[0].Im == 0) {
				return regularModifiedBesselFunction1stKind(tmp[0].Re, tmp[1]);
			}
			else {
				return regularModifiedBesselFunction1stKind(tmp[0], tmp[1]);
			}
		}
		if (str.find(L"risingFactorial(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"(") + 1));
			return risingFactorial(tmp[0], tmp[1].Re);
		}
		if (str.find(L"sec(") != std::wstring::npos && str.find(L"asec(") == std::wstring::npos && str.find(L"sech(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"sec(") + 1));
			return sec(tmp);
		}
		if (str.find(L"Si(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return Si(tmp);
		}
		if (str.find(L"sin(") != std::wstring::npos && str.find(L"asin(") == std::wstring::npos && str.find(L"sinh(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"sin(") + 1));
			return sin(tmp);
		}
		if (str.find(L"sinc(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return sinc(tmp);
		}
		if (str.find(L"sinh(") != std::wstring::npos && str.find(L"asinh(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"sinh(") + 1));
			return sinh(tmp);
		}
		/*		if (str.find(L"sgn(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"sgn(") + 1));
		return sgn(tmp);
		}
		*/
		if (str.find(L"sphericalBesselFunctionFirstKind(") != std::wstring::npos) {
			std::vector<ComplexNumber> tmp = ParseNComplexNumbers(str.substr(str.find(L"(") + 1));
			if (tmp[0].Im == 0) {
				return sphericalBesselFunction1stKind(tmp[0].Re, tmp[1]);
			}
			else {
				return sphericalBesselFunction1stKind(tmp[0], tmp[1]);
			}
		}
		if (str.find(L"sqrt(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"sqrt(") + 1));
			return sqrt(tmp);
		}
		/*		if (str.find(L"squarewave(") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"squarewave(") + 1));
		return squarewave(tmp);
		}
		if (str.find(L"step(") != std::wstring::npos) {//Heaviside step function
		std::wstring s = str.substr(str.find(L"step("));
		s = s.substr(str.find(L"(") + 1, str.find(L")") - 1);
		ComplexNumber tmp = ParseComplexNumber(s);
		return step(tmp);
		}
		*/		if (str.find(L"tan(") != std::wstring::npos && str.find(L"tanh(") == std::wstring::npos && str.find(L"atan(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"tan(") + 1));
			return tan(tmp);
		}
		if (str.find(L"tanh(") != std::wstring::npos && str.find(L"atanh(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"tanh(") + 1));
			return tanh(tmp);
		}
		/*		if (str.find(L"tri(") != std::wstring::npos) {
		std::wstring s = str.substr(str.find(L"tri("));
		s = s.substr(str.find(L"(") + 1, str.find(L")") - 1);
		ComplexNumber tmp = ParseComplexNumber(s);
		return tri(tmp);
		}
		*//*	if (str.find(L"trunc(L") != std::wstring::npos) {
		ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"acos(") + 1));
		return acos(tmp);
		}*/
		if (str.find(L"WrightOmegaFunction(") != std::wstring::npos && str.find(L"WrightOmegaFunction(") == std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"WrightOmegaFunction(") + 1));
			return tanh(tmp);
		}
		if (str.find(L"zeta(") != std::wstring::npos) {
			ComplexNumber tmp = ParseComplexNumber(str.substr(str.find(L"(") + 1));
			return RiemannZetaFunction(tmp);
		}
		return 0;	//else, if none of the other functions match, return NULL
	}

	ComplexNumber ComplexFunction::evaluateComplex(std::vector<ComplexNumber> vals) {
		/*	NOTE:
		Suitable arithmetic operators are +,-,*,/,^.  Suitable functions are found in evalFunction method.
		This function will properly evaluate a function of the form f(x,y) = x^2 without error.		*/

		/*=================
		|DECLARE VARIABLES|
		=================*/
		//make sure when using to_string precision or other functions which truncate values, we do as little truncation as possible
		int prcs = precision;
		precision = 30;
		std::wstring funct = f.function;//local copy of function string, necessary as we will erase/replace parts of the string as we go
		std::vector<ComplexNumber> dbls;	 //copy of every double value
		std::vector<wchar_t> ops;		 //copy of every operator in left-to-right order
		funct = removeSpaces(funct);//remove all spaces in function string
		funct = replaceConstants(funct);//replace all system constants in string
		//========================================

		/*================================
		|REPLACE VARIABLES WITH CONSTANTS|
		================================*/
		while (vals.size() < f.variables.size()) { vals.push_back(ComplexNumber(0)); }//make sure there are enough values for every variable
		//replace all variables with their ComplexNumber value from double* vals
		if (vals.size() != NULL) {
			for (int i = 0; i < f.variables.size(); ++i) {
				std::wstring val = L"(";
				val.append(vals[i].toStringBothParts());//double values for variable made into a string
				val.append(L")");

				//check string for spots to replace char, then perform replace function
				for (int j = 0; j < funct.size(); ++j) {
					if (j == 0 && funct[j] == f.variables[i] && isLetter(funct[j + 1]) == false) {
						funct.erase(j, 1);
						funct.insert(j, val);
						j = 0;
					}
					if (j>0 && j<funct.size() - 1 && funct[j] == f.variables[i] && isLetter(funct[j + 1]) == false && isLetter(funct[j - 1]) == false) {
						funct.erase(j, 1);
						funct.insert(j, val);
					}
					if (j == funct.size() - 1 && funct[j] == f.variables[i] && isLetter(funct[j - 1]) == false) {
						funct.erase(j, 1);
						funct.insert(j, val);
					}
				}

			}
		}//================================================================
		funct = convertDoubleToComplex(funct);

		/*================
		|HANDLE FUNCTIONS|
		================*/
		bool searchForFuncts = hasFunction(funct);
		while (searchForFuncts == true) {
			bool foundOne = false;
			for (int i = 0; i < funct.size() - 1; ++i) {
			topOfTheFunctsLoop:
				if (isLetter(funct[i]) == true && isLetter(funct[i + 1]) == true) {		//find functions, which have more than one letter in a row
					std::wstring tempstr = funct.substr(i);
					size_t start = tempstr.find(L'(');
					size_t end = findMatchingParen(tempstr, start);
					tempstr = tempstr.substr(0, end+1);//save function (i.e. 'exp(12.2)') as a substr to evaluate
					if (hasFunction(tempstr.substr(tempstr.find(L"("))) == true) {//keep searching for the innermost paren without a function inside
						i += tempstr.find(L"(");//skip past the function (i.e. i+=4 if outer function is 'cos(' for example)
						goto topOfTheFunctsLoop;
					}
					std::wstring toReplace = tempstr;

					//make sure there's nothing within the parenthesis that needs to be simplified (i.e. 'sin(12-2)')
					if (hasOperatorNotComplexNumber(tempstr) == true && 
						countComplexNumbers(tempstr)>1 && 
						tempstr.find(L',')==std::wstring::npos
						) {
						size_t a = tempstr.find(L"(");
						size_t b = findMatchingParen(tempstr, a);
						std::wstring chk = tempstr.substr(a, b - a+1);
						ComplexFunction f2(chk);
						ComplexNumber x = f2.evaluateComplex(vals);
						std::wstring replc = x.toStringBothParts();
						if (replc.find(L"(") == std::wstring::npos) {
							std::wstring copystr = L"(";
							copystr.append(replc);
							replc = copystr;
							copystr.clear();
						}
						if (replc.find(L")") == std::wstring::npos) {replc.append(L")");}
						tempstr = replaceString(tempstr, chk, replc);
					}

					//now evaluate the function, get the value, and replace the whole thing in the original string (i.e. 'sin(2+3i)' --> '(9.5+4.3i)')
					ComplexNumber tmp1 = evalComplexFunction(tempstr);
					std::wstring replacer = L"(";
					replacer.append(tmp1.toStringBothParts());
					replacer.append(L")");
					funct = replaceString(funct, toReplace, replacer);
					foundOne = true;
					i = 0;
				}
			}
			if (foundOne == false) { searchForFuncts = false; }
		}//====================================================

		/*================
		|  HANDLE PARENS |
		================*/
		funct = processInnerTermsComplex(funct);

		/*=================
		|HANDLE ARITHMETIC|
		=================*/
		ComplexNumber answer = calculateArithmeticComplex(funct);
		if (answer.Re != answer.Re || answer.Im != answer.Im) { 
			answer = ComplexNumber();
		}
		//finally, restore precision and return a value
		precision = prcs;
		return answer; //note: if unsuccessful, will return 0
	}//END EVALUATE FUNCTION================================================

	ComplexNumber ComplexFunction::evaluateComplex(ComplexNumber z) {
		std::vector<ComplexNumber> vals;
		vals.push_back(z);
		ComplexNumber n = evaluateComplex(vals);
		vals.clear();
		return n;
	}

	ComplexNumber calculateArithmeticComplex(std::wstring funct) {
		int start = 0;
		int end = 0;
		std::vector<ComplexNumber> dbls;
		std::vector<wchar_t> ops;

		//load all ComplexNumbers into the array
		while (findComplexNumber(funct) >= 0) {
			std::wstring cpnum = getNearestComplexNumber(funct);
			dbls.push_back(ComplexNumber(ParseComplexNumber(cpnum)));
			funct.erase(findComplexNumber(funct), cpnum.size());
		}
		if (dbls.size() == 1) { return dbls[0]; }

		//collect leftover operators for arithmetic evaluation
		while (hasOperator(funct) == true) {
			for (int i = 0; i < funct.length(); ++i) {
				if (isOperatorNotParen(funct[i]) == true) {
					ops.push_back(funct[i]);
					funct.erase(i, 1);
				}
			}
		}
		funct.clear();

		if (ops.size() == 0) { return ParseComplexNumber(funct); }//if operator vector is empty, there are no operators, only a value to return

		ComplexNumber answer(0);
		int pass = 0;

		//evaluate in PEMDMAS order, 6 passes through the vector, eliminating evaluated terms and operators as you go
		while (dbls.size()>1) {

			for (int i = 0; i < ops.size(); ++i) {
				//note: format = 'a+b'
				ComplexNumber a = dbls[i];
				ComplexNumber b = dbls[i + 1];
				if (i >= 0 && ops[i] == L'+' && pass == 4) {
					answer = a + b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == L'-' && pass == 5) {
					answer = a - b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == L'*' && pass == 1) {
					answer = a * b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == L'/' && pass == 2) {
					if (b == ComplexNumber(0)) { 
						return (ComplexNumber());
					}else {
						answer = a / b;
					}
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == L'^' && pass == 0) {
					answer = pow(a, b);
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == L'%'&& pass == 3) {
					answer = (int)a.modulus() % (int)b.modulus();
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
			}
			++pass;
			pass %= 6;//keep going around until only one dbl remains
		}
		return answer; //note: if unsuccessful, will return 0
	}

	double calculateArithmetic(std::wstring funct){
	std::vector<double> dbls;	 //copy of every double value
	std::vector<wchar_t> ops;	 //copy of every operator in left-to-right order

	//collect leftover operators for arithmetic evaluation
	//while (hasOperator(funct) == true) {
		for (int i = 0; i < funct.length(); ++i) {
			if (
				(isOperatorNotParen(funct[i]) == true && funct[i]!=L'-') ||
				(
					funct[i]==L'-' && 
					( 
						(i==0 && isNumber(funct[i+1])==true) || 
						(i>0 && funct[i-1]!=L',' && funct[i-1]!=L')' && isNumber(funct[i-1])==false)
					)
				)
				) {
				ops.push_back(funct[i]);
				funct.insert(funct.begin()+i, L',');
				funct.erase(i+1, 1);
			}
		}
	//}
	//load all Doubles into the array
	/*while (findNearestDouble(funct) >= 0) {
		std::wstring cpnum = getNearestDouble(funct);
		dbls.push_back(ParseInputNumber(cpnum));
		funct.erase(findNearestDouble(funct), cpnum.size());
	}*/
	funct = replaceString(funct, L"(", L" ");
	funct = replaceString(funct, L")", L" ");
	funct = removeSpaces(funct);
	funct = funct.append(L")");
	dbls = ParseNInputNumbers(funct);
	if (dbls.size() == 1) { return dbls[0]; }


	funct.clear();

	if (ops.size() == 0) { return ParseInputNumber(funct); }//if operator vector is empty, there are no operators, only a value to return
	double answer = 0;
	int pass = 0;

	//evaluate in PEMDMAS order, 6 passes through the vector, eliminating evaluated terms and operators as you go
	while (dbls.size()>1) {

		for (int i = 0; i < ops.size(); ++i) {
			//note: format = 'a+b'

			double a = dbls[i];
			double b = dbls[i + 1];
			if (i >= 0 && ops[i] == L'+' && pass == 4) {
				answer = a + b;
				dbls[i + 1] = answer;
				dbls.erase(dbls.begin() + i);
				ops.erase(ops.begin() + i);
				--i;
			}
			if (i >= 0 && ops[i] == L'-' && pass == 5) {
				answer = a - b;
				dbls[i + 1] = answer;
				dbls.erase(dbls.begin() + i);
				ops.erase(ops.begin() + i);
				--i;
			}
			if (i >= 0 && ops[i] == L'*' && pass == 1) {
				answer = a * b;
				dbls[i + 1] = answer;
				dbls.erase(dbls.begin() + i);
				ops.erase(ops.begin() + i);
				--i;
			}
			if (i >= 0 && ops[i] == L'/' && pass == 2) {
				if (b == 0) {
					return NAN;//division by 0 is impossible
				}
				else {
					answer = a / b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
			}
			if (i >= 0 && ops[i] == L'^' && pass == 0) {
				answer = pow(a, b);
				dbls[i + 1] = answer;
				dbls.erase(dbls.begin() + i);
				ops.erase(ops.begin() + i);
				--i;
			}
			if (i >= 0 && ops[i] == L'%'&& pass == 3) {
				answer = (int)a % (int)b;
				dbls[i + 1] = answer;
				dbls.erase(dbls.begin() + i);
				ops.erase(ops.begin() + i);
				--i;
			}
		}
		++pass;
		pass %= 6;//keep going around until only one dbl remains
	}
	return answer;
	}

	bool ComplexFunction::hasDiscontinuities(std::vector<ComplexNumber> testVals) {//for complex-valued univariate functions
		double epsil = 0.0000001;
		ComplexNumber vals[8] = {
			ComplexNumber(epsil,epsil), ComplexNumber(0,epsil), 
			ComplexNumber(epsil,0),     ComplexNumber(0,-epsil), 
			ComplexNumber(-epsil,0),    ComplexNumber(epsil,-epsil), 
			ComplexNumber(-epsil,epsil),ComplexNumber(-epsil,-epsil)
		};

		for (int i = 0; i < testVals.size(); ++i) {
			for (int k = 0; k < 9; ++k) {
				ComplexNumber left = evaluateComplex(testVals[i] + vals[k]);
				ComplexNumber right = evaluateComplex(testVals[i] - vals[k]);
				if ( 
					isNan(left) ||
					isNan(right) ||
					(0.5*(left - right).modulus()) > epsil
					) 
				{ return true; }
			}
		}
		return false;
	}

	bool ComplexFunction::hasDiscontinuities(ComplexMatrix testVals) {
		double epsil = 0.0000001;
		ComplexNumber vals[9] = {
			ComplexNumber(0,0),
			ComplexNumber(epsil,epsil), ComplexNumber(0,epsil),
			ComplexNumber(epsil,0),     ComplexNumber(0,-epsil),
			ComplexNumber(-epsil,0),    ComplexNumber(epsil,-epsil),
			ComplexNumber(-epsil,epsil),ComplexNumber(-epsil,-epsil)
		};

		for (int i = 0; i < testVals.rows; ++i) {
			for (int j = 0; j < testVals.columns; ++j) {
				for (int k = 0; k < 9; ++k) {
					ComplexNumber left = evaluateComplex(testVals.element[i*testVals.rows+j] + vals[k]);
					ComplexNumber right = evaluateComplex(testVals.element[i*testVals.rows+j] - vals[k]);
					if (
						isNan(left) ||
						isNan(right) ||
						(0.5*(left - right).modulus()) > epsil
						)
					{
						return true;
					}
				}
			}
		}
		return false;
	}

	std::wstring ComplexFunction::variablesToString() { return f.variablesToString(); }
	std::wstring ComplexFunction::toString() { return f.toString(); }
	void ComplexFunction::display() { std::wcout << toString() << std::endl; }
//================================================================
	ParametricFunction::ParametricFunction() {}
	ParametricFunction::~ParametricFunction() {}

	ParametricFunction::ParametricFunction(std::wstring str) {
			if (str.find(L"=") != std::wstring::npos) {	//case:  string is 'f(x) = 2*x + exp(-y)'
				while (str[0] != L'(') { str = str.substr(1); }//clip off everything leading up to the variable definition
				str = str.substr(1);
				for (int i = 0; i < str.find(L"=") + 1; ++i) {
					if (isLetter(str[i]) == true) { variables.push_back(str[i]); }
				}
				function = str.substr(str.find(L"=") + 1);
			}
			else { function = str; }					//case:  string has no variables, i.e. '23-12+exp(3)'
		}

	ParametricFunction::ParametricFunction(std::vector<wchar_t> params, std::wstring functs) {
			if (variables.size() == function.size()) {
				function = functs;
				variables = params;
			}
		}

	ParametricFunction::ParametricFunction(std::wstring vars, std::wstring funct) {
			for (int i = 0; i < vars.size(); ++i) {
				if (isLetter(vars[i]) == true) { variables.push_back(vars[i]); }
			}
			function = funct;
		}

		//class methods
		//=============
		std::wstring ParametricFunction::variablesToString() { //returns string of variables in format 'x,y,z,w,p,...'
			std::wstring vars = L"";
			if (variables.size() == 1) { vars.push_back(variables[0]); return vars; }
			if (variables.size() > 1) {
				for (int i = 0; i < variables.size() - 1; ++i) {
					vars.push_back(variables[i]);
					vars.push_back(',');
				}
				vars.push_back(variables[variables.size() - 1]);
				return vars;
			}
			return vars;
		}

		std::wstring ParametricFunction::toString() { return function; }	//returns function in format 'x * exp(y+2)...'
		void ParametricFunction::display() { std::wcout << toString() << L"\n"; }
		//=======================================================================

		bool testForFunctions(std::wstring in) {
			if (hasFunction(in) == true) { return true; }
			if (
				in.find(L"Vector") != std::wstring::npos ||
				in.find(L"Matrix") != std::wstring::npos ||
				in.find(L"Quaternion") != std::wstring::npos ||
				in.find(L"ComplexNumber") != std::wstring::npos ||
				in.find(L"ComplexMatrix") != std::wstring::npos ||
				in.find(L"CM[") != std::wstring::npos ||
				in.find(L"S[") != std::wstring::npos ||
				in.find(L"ODE") != std::wstring::npos ||
				in.find(L"PDE") != std::wstring::npos ||
				in.find(L"Octonian") != std::wstring::npos ||
				in.find(L"DualNumber") != std::wstring::npos ||
				in.find(L"Group") != std::wstring::npos ||
				in.find(L"Tensor") != std::wstring::npos ||
				in.find(L"Octonian") != std::wstring::npos ||
				in.find(L"InfiniteSeries") != std::wstring::npos) {
				return false;
			}
			return false;
		}