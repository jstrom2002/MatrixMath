#pragma once
#include "stdafx.h"

InfiniteSeries::InfiniteSeries() {}
InfiniteSeries::~InfiniteSeries() {}
InfiniteSeries::InfiniteSeries(Function f1) {
	f = f1;
}
InfiniteSeries::InfiniteSeries(std::vector<wchar_t> var, std::wstring funct) {
	f = Function(var, funct);
	isProductSeries = false;
}


//ARITHMETIC METHODS
//===================
InfiniteSeries InfiniteSeries::add(InfiniteSeries a, InfiniteSeries b) {
	a.f.function.append(L"+");
	a.f.function.append(b.f.function);
	return InfiniteSeries(f.combineVariables(a.f.variables, b.f.variables), a.f.function);
}
InfiniteSeries InfiniteSeries::add(InfiniteSeries a, double b)
{
	a.f.function.append(L"+");
	a.f.function.append(to_stringPrecision(b));
	return InfiniteSeries(a.f.variables, a.f.function);
}
InfiniteSeries InfiniteSeries::subtract(InfiniteSeries a, InfiniteSeries b) {
	a.f.function.append(L"-");
	a.f.function.append(b.f.function);
	return InfiniteSeries(f.combineVariables(a.f.variables, b.f.variables), a.f.function);
}
InfiniteSeries InfiniteSeries::subtract(InfiniteSeries a, double b)
{
	a.f.function.append(L"-");
	a.f.function.append(to_stringPrecision(b));
	return InfiniteSeries(a.f.variables, a.f.function);
}
InfiniteSeries InfiniteSeries::multiply(InfiniteSeries a, InfiniteSeries b) {
	a.f.function.append(L"*");
	a.f.function.append(b.f.function);
	return InfiniteSeries(f.combineVariables(a.f.variables, b.f.variables), a.f.function);
}
InfiniteSeries InfiniteSeries::multiply(InfiniteSeries a, double b)
{
	a.f.function.append(L"*");
	a.f.function.append(to_stringPrecision(b));
	return InfiniteSeries(a.f.variables, a.f.function);
}
InfiniteSeries InfiniteSeries::divide(InfiniteSeries a, InfiniteSeries b) {
	a.f.function.append(L"/");
	a.f.function.append(b.f.function);
	return InfiniteSeries(f.combineVariables(a.f.variables, b.f.variables), a.f.function);
}
InfiniteSeries InfiniteSeries::divide(InfiniteSeries a, double b)
{
	a.f.function.append(L"+");
	a.f.function.append(to_stringPrecision(b));
	return InfiniteSeries(a.f.variables, a.f.function);
}

InfiniteSeries InfiniteSeries::inverse(InfiniteSeries a) {
	std::wstring func = L"1/";
	func.append(a.f.function);
	return InfiniteSeries(a.f.variables, func);
}

InfiniteSeries InfiniteSeries::exponent(InfiniteSeries a, InfiniteSeries b)
{
	a.f.function.append(L"^");
	a.f.function.append(b.f.function);
	return InfiniteSeries(f.combineVariables(a.f.variables, b.f.variables), a.f.function);
}
InfiniteSeries InfiniteSeries::exponent(InfiniteSeries a, double b)
{
	a.f.function.append(L"^");
	a.f.function.append(to_stringPrecision(b));
	return InfiniteSeries(a.f.variables, a.f.function);
}

//OVERLOADED OPERATORS
//====================
InfiniteSeries& InfiniteSeries::operator*=(InfiniteSeries& rhs) { *this = multiply(*this, rhs); return *this; }
InfiniteSeries& InfiniteSeries::operator*=(double& rhs) { *this = multiply(*this, rhs); return *this; }
InfiniteSeries& InfiniteSeries::operator/=(InfiniteSeries& rhs) { *this = divide(*this, rhs); return *this; }
InfiniteSeries& InfiniteSeries::operator/=(double& rhs) { *this = divide(*this, rhs); return *this; }
InfiniteSeries& InfiniteSeries::operator+=(InfiniteSeries& rhs) { *this = add(*this, rhs); return *this; }
InfiniteSeries& InfiniteSeries::operator+=(double& rhs) { *this = add(*this, rhs); return *this; }
InfiniteSeries& InfiniteSeries::operator-=(InfiniteSeries& rhs) { *this = subtract(*this, rhs); return *this; }
InfiniteSeries& InfiniteSeries::operator-=(double& rhs) { *this = subtract(*this, rhs); return *this; }

InfiniteSeries& InfiniteSeries::operator*(InfiniteSeries& rhs) { return multiply(*this, rhs); }
InfiniteSeries& InfiniteSeries::operator*(double& rhs) { return multiply(*this, rhs); }
InfiniteSeries& InfiniteSeries::operator/(InfiniteSeries& rhs) { return divide(*this, rhs); }
InfiniteSeries& InfiniteSeries::operator/(double& rhs) { return divide(*this, rhs); }
InfiniteSeries& InfiniteSeries::operator+(InfiniteSeries& rhs) { return add(*this, rhs); }
InfiniteSeries& InfiniteSeries::operator+(double& rhs) { return add(*this, rhs); }
InfiniteSeries& InfiniteSeries::operator-(InfiniteSeries& rhs) { return subtract(*this, rhs); }
InfiniteSeries& InfiniteSeries::operator-(double& rhs) { return subtract(*this, rhs); }
InfiniteSeries& InfiniteSeries::operator^(InfiniteSeries x) { return exponent(*this, x); }
InfiniteSeries& InfiniteSeries::operator^(double x) { return exponent(*this, x); }

double InfiniteSeries::evaluate(int lo, int hi) {//evaluate series from some integer to another integer (i.e. from 0 to 15)
	if (isProductSeries == true) {
		double ans = 1;
		for (int i = lo; i <= hi; ++i) {
			ans *= f.evaluate(i);
		}
		return ans;
	}
	else {
		double ans = 0;
		for (int i = lo; i <= hi; ++i) {
			ans += f.evaluate(i);
		}
		return ans;
	}
}

double InfiniteSeries::evaluate(std::vector<double> vals, int lo, int hi) {//evaluate series from some integer to another integer (i.e. from 0 to 15)
	if (isProductSeries == true) {
		double ans = 1;
		for (int i = lo; i <= hi; ++i) {
			vals.push_back(i);
			ans *= f.evaluate(vals);
			vals.pop_back();
		}
		return ans;
	}
	else {
		double ans = 0;
		for (int i = lo; i <= hi; ++i) {
			vals.push_back(i);
			ans += f.evaluate(vals);
			vals.pop_back();
		}
		return ans;
	}
}

double InfiniteSeries::evaluateLoToInf(double lo, double hi) {
	/*evaluate from some lower bound to an upper bound of infinity.
	the second parameter is the farthest bound to calculate to
	before integrating to calculate the remainder. */
	double ans = 0;
	for (int i = lo; i <= hi; ++i) {
		ans += f.evaluate(i);
	}

	/*calculate remainder of series as integral. This is done by integrating
	from 0 to inf and subtract from 0 to lo, thus integrating from lo to inf */
	ans += integrateGaussLaguerreQuadrature(f) - integrateGaussLegendreQuadrature(0, lo, f);
	return ans;
}

void InfiniteSeries::setProductSeries() { isProductSeries = true; }

std::wstring InfiniteSeries::testConvergence(unsigned int n) {
	//ratio test
	std::wstring answer = L"";
	double test = abs(evaluate(n + 1, n + 1) / evaluate(n, n));
	if (test < 1) { answer = L"true"; }
	if (test > 1) { answer = L"false"; }
	answer = L"inconclusive";


	return answer;
}

std::wstring InfiniteSeries::toString() {
	return f.toString();
}

void InfiniteSeries::display() { std::wcout << toString() << std::endl; }
//===================================================================


ComplexInfiniteSeries::ComplexInfiniteSeries() {}
ComplexInfiniteSeries::~ComplexInfiniteSeries() {}
ComplexInfiniteSeries::ComplexInfiniteSeries(ComplexFunction f1) {
		cf = f1;
		isProductSeries = false;
	}
ComplexInfiniteSeries::ComplexInfiniteSeries(std::vector<wchar_t> var, std::wstring funct) {
		cf = ComplexFunction(var, funct);
		isProductSeries = false;
	}

	//ARITHMETIC METHODS
	//===================
	ComplexInfiniteSeries ComplexInfiniteSeries::add(ComplexInfiniteSeries a, ComplexInfiniteSeries b) {
		a.cf.f.function.append(L"+");
		a.cf.f.function.append(b.cf.f.function);
		return ComplexInfiniteSeries(cf.f.combineVariables(a.cf.f.variables, b.cf.f.variables), a.cf.f.function);
	}
	ComplexInfiniteSeries ComplexInfiniteSeries::add(ComplexInfiniteSeries a, double b)
	{
		a.cf.f.function.append(L"+");
		a.cf.f.function.append(to_stringPrecision(b));
		return ComplexInfiniteSeries(a.cf.f.variables, a.cf.f.function);
	}
	ComplexInfiniteSeries ComplexInfiniteSeries::subtract(ComplexInfiniteSeries a, ComplexInfiniteSeries b) {
		a.cf.f.function.append(L"-");
		a.cf.f.function.append(b.cf.f.function);
		return ComplexInfiniteSeries(cf.f.combineVariables(a.cf.f.variables, b.cf.f.variables), a.cf.f.function);
	}
	ComplexInfiniteSeries ComplexInfiniteSeries::subtract(ComplexInfiniteSeries a, double b)
	{
		a.cf.f.function.append(L"-");
		a.cf.f.function.append(to_stringPrecision(b));
		return ComplexInfiniteSeries(a.cf.f.variables, a.cf.f.function);
	}
	ComplexInfiniteSeries ComplexInfiniteSeries::multiply(ComplexInfiniteSeries a, ComplexInfiniteSeries b) {
		a.cf.f.function.append(L"*");
		a.cf.f.function.append(b.cf.f.function);
		return ComplexInfiniteSeries(cf.f.combineVariables(a.cf.f.variables, b.cf.f.variables), a.cf.f.function);
	}
	ComplexInfiniteSeries ComplexInfiniteSeries::multiply(ComplexInfiniteSeries a, double b)
	{
		a.cf.f.function.append(L"*");
		a.cf.f.function.append(to_stringPrecision(b));
		return ComplexInfiniteSeries(a.cf.f.variables, a.cf.f.function);
	}
	ComplexInfiniteSeries ComplexInfiniteSeries::divide(ComplexInfiniteSeries a, ComplexInfiniteSeries b) {
		a.cf.f.function.append(L"/");
		a.cf.f.function.append(b.cf.f.function);
		return ComplexInfiniteSeries(cf.f.combineVariables(a.cf.f.variables, b.cf.f.variables), a.cf.f.function);
	}
	ComplexInfiniteSeries ComplexInfiniteSeries::divide(ComplexInfiniteSeries a, double b)
	{
		a.cf.f.function.append(L"+");
		a.cf.f.function.append(to_stringPrecision(b));
		return ComplexInfiniteSeries(a.cf.f.variables, a.cf.f.function);
	}

	ComplexInfiniteSeries ComplexInfiniteSeries::inverse(ComplexInfiniteSeries a) {
		std::wstring func = L"1/";
		func.append(a.cf.f.function);
		return ComplexInfiniteSeries(a.cf.f.variables, func);
	}

	ComplexInfiniteSeries ComplexInfiniteSeries::exponent(ComplexInfiniteSeries a, ComplexInfiniteSeries b)
	{
		a.cf.f.function.append(L"^");
		a.cf.f.function.append(b.cf.f.function);
		return ComplexInfiniteSeries(cf.f.combineVariables(a.cf.f.variables, b.cf.f.variables), a.cf.f.function);
	}
	ComplexInfiniteSeries ComplexInfiniteSeries::exponent(ComplexInfiniteSeries a, double b)
	{
		a.cf.f.function.append(L"^");
		a.cf.f.function.append(to_stringPrecision(b));
		return ComplexInfiniteSeries(a.cf.f.variables, a.cf.f.function);
	}

	//OVERLOADED OPERATORS
	//====================
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator*=(ComplexInfiniteSeries& rhs) { *this = multiply(*this, rhs); return *this; }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator*=(double& rhs) { *this = multiply(*this, rhs); return *this; }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator/=(ComplexInfiniteSeries& rhs) { *this = divide(*this, rhs); return *this; }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator/=(double& rhs) { *this = divide(*this, rhs); return *this; }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator+=(ComplexInfiniteSeries& rhs) { *this = add(*this, rhs); return *this; }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator+=(double& rhs) { *this = add(*this, rhs); return *this; }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator-=(ComplexInfiniteSeries& rhs) { *this = subtract(*this, rhs); return *this; }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator-=(double& rhs) { *this = subtract(*this, rhs); return *this; }

	ComplexInfiniteSeries& ComplexInfiniteSeries::operator*(ComplexInfiniteSeries& rhs) { return multiply(*this, rhs); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator*(double& rhs) { return multiply(*this, rhs); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator/(ComplexInfiniteSeries& rhs) { return divide(*this, rhs); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator/(double& rhs) { return divide(*this, rhs); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator+(ComplexInfiniteSeries& rhs) { return add(*this, rhs); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator+(double& rhs) { return add(*this, rhs); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator-(ComplexInfiniteSeries& rhs) { return subtract(*this, rhs); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator-(double& rhs) { return subtract(*this, rhs); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator^(ComplexInfiniteSeries x) { return exponent(*this, x); }
	ComplexInfiniteSeries& ComplexInfiniteSeries::operator^(double x) { return exponent(*this, x); }

	ComplexNumber ComplexInfiniteSeries::evaluate(double lo, double hi) {
		ComplexNumber ans = 0;
		for (int i = lo; i <= hi; ++i) {
			ans += cf.evaluateComplex(ComplexNumber(i));
		}
		return ans;
	}

	ComplexNumber ComplexInfiniteSeries::evaluate(ComplexNumber vals2, int lo, int hi) {//evaluate series from some integer to another integer (i.e. from 0 to 15)
		std::vector<ComplexNumber> vals;
		vals.push_back(vals2);
		if (isProductSeries == true) {
			ComplexNumber ans = 1;
			for (double i = lo; i <= hi; ++i) {
				ComplexNumber i2(i);
				vals.push_back(i2);
				ans *= cf.evaluateComplex(vals);
				vals.pop_back();
			}
			return ans;
		}
		else {
			ComplexNumber ans = 0;
			for (double i = lo; i <= hi; ++i) {
				ComplexNumber i2(i);
				vals.push_back(i2);
				ans += cf.evaluateComplex(vals);
				vals.pop_back();
			}
			return ans;
		}
	}

	ComplexNumber ComplexInfiniteSeries::evaluate(std::vector<ComplexNumber> vals, int lo, int hi) {//evaluate series from some integer to another integer (i.e. from 0 to 15)
		if (isProductSeries == true) {
			ComplexNumber ans = 1;
			for (int i = lo; i <= hi; ++i) {
				ComplexNumber i2(i);
				vals.push_back(i2);
				ans *= cf.evaluateComplex(vals);
				vals.pop_back();
			}
			return ans;
		}
		else {
			ComplexNumber ans = 0;
			for (int i = lo; i <= hi; ++i) {
				ComplexNumber i2(i);
				vals.push_back(i2);
				ans = ans.add(ans,cf.evaluateComplex(vals));
				vals.pop_back();
			}
			return ans;
		}
	}

	ComplexNumber ComplexInfiniteSeries::evaluateWithCoefficients(std::vector<ComplexNumber> coef, double lo, double hi) {
		ComplexNumber ans = 0;
		for (int i = lo; i <= hi; ++i) {
			ans += coef[i].multiply(coef[i], cf.evaluateComplex(ComplexNumber(i)));
		}
		return ans;
	}

	ComplexNumber ComplexInfiniteSeries::evaluateLoToInf(std::vector<ComplexNumber> vals, double lo, double hi) {
		/*evaluate from some lower bound to an upper bound of infinity.
		the second parameter is the farthest bound to calculate to
		before integrating to calculate the remainder. */
		ComplexNumber ans;
		if (isProductSeries == true) {
			ans = ComplexNumber(1);
			for (int i = lo; i <= hi; ++i) {
				ComplexNumber i2(i);
				vals.push_back(i2);
				ans *= cf.evaluateComplex(vals);
				vals.pop_back();
			}
		}else {
			ans = ComplexNumber(0);
			for (int i = lo; i <= hi; ++i) {
				ComplexNumber i2(i);
				vals.push_back(i2);
				ans += cf.evaluateComplex(vals);
				vals.pop_back();
			}
		}

		/*calculate remainder of series as integral. This is done by integrating
		from 0 to inf and subtract from 0 to lo, thus integrating from lo to inf */
		ans += integrateGaussLaguerreQuadrature(cf.f) - integrateGaussLegendreQuadrature(0, lo, cf.f);
		return ans;
	}

	void ComplexInfiniteSeries::setProductSeries() { isProductSeries = true; }

	bool ComplexInfiniteSeries::testConvergence() {
		//ADD LATER!
		return false;
	}

	std::wstring ComplexInfiniteSeries::toString() {
		return cf.f.toString();
	}

	void ComplexInfiniteSeries::display() { std::wcout << toString() << std::endl; }
//================================================================================