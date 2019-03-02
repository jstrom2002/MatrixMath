#ifndef NUMBERTHEORY_H
#define NUMBERTHEORY_H
#pragma once

double csc(double x);
double sec(double x);
double cot(double x);
double sech(double x);
double csch(double x);
double coth(double x);
double asech(double x);
double acsch(double x);
double acoth(double x);
int AckermannsFunction(int m, int n);
double BernoulliNumber(int n);
double BernoulliNumber(int m, int n);
float CarmackFastInverseSqrt(float x);
int CollatzConjecture(int n);
long long unsigned int doubleFactorial(int n);
double factorial(double n);
double FibonacciNumber(int n);
bool isSquarefree(int n);
double logarithm(double n, double base);
int permutation(int n, int r);
int combination(int n, int r);
int KroneckerDelta(double i, double j);
bool isPowerOfTwo(int x);
int prime(int n);
int sumPrimes(int n);
int productPrimes(int n);
int primeCountingFunction(int n);
int MoebiusFunction(int n);
double fallingFactorial(double n, double r);
ComplexNumber fallingFactorial(ComplexNumber n, double r);
double risingFactorial(double n, double r);
ComplexNumber risingFactorial(ComplexNumber n, double r);
double hypergeometricFunction(std::vector<double> vec);
ComplexNumber hypergeometricFunction(std::vector<double> vec, ComplexNumber z2);
ComplexNumber hypergeometricFunction(std::vector<ComplexNumber> vec);
double betaFunction(double x, double y);
double betaFunction(std::vector<double> x);
ComplexNumber betaFunction(ComplexNumber x, ComplexNumber y);
ComplexNumber betaFunction(std::vector<ComplexNumber> x);
double betaFunction(double x, double a, double b);
double incompleteBetaFunction(double x, double a, double b);
double completeEllipticIntegral1stKind(double k);
double completeEllipticIntegral2ndKind(double k);
double completeEllipticIntegral3rdKind(double k, double v);
double incompleteEllipticIntegral1stKind(double k, double v);
double incompleteEllipticIntegral2ndKind(double k, double v);
double incompleteEllipticIntegral3rdKind(double k, double v, double psi);
long double factorialApproximation(long double x);
long double gammaApproximation(long double x);
double cylindricalNeumannFunction(double v, double x);
double sphericalNeumannFunction(double n, double x);
double sphericalAssociatedLegendreFunction(double n, double x, double theta);
double lowerIncompleteGamma(double s, double x);
double upperIncompleteGamma(double s, double x);
double regularizedGamma(double s, double x);
double gamma(double x);
ComplexNumber gamma(ComplexNumber z);
ComplexNumber gamma(ComplexNumber z, double epsiln);
double logGamma(double z);
ComplexNumber logGamma(ComplexNumber z);
double polygamma(double z, double order);
ComplexNumber polygamma(ComplexNumber z, double order);
double digamma(double z);
ComplexNumber digamma(ComplexNumber z);
int pentagonalNumber(int n);
double percentError(double theoretical, double experimental);
long long unsigned int convertToBinary(int a);
double LambertW(double x);
ComplexNumber LambertW(ComplexNumber z);
double WrightOmegaFunction(double x);
ComplexNumber WrightOmegaFunction(ComplexNumber z);
double RiemannZetaFunction(double x);
ComplexNumber RiemannZetaFunction(ComplexNumber x);
long long unsigned int BellNumber(double n);
int StirlingNumber1stKind(int n, int k);
int UnsignedStirlingNumber1stKind(int n, int k);
double StirlingNumber2ndKind(double n, double k);
int EulersTotientFunction(int n);
double        BesselFunction1stKind(double x, double order);
ComplexNumber BesselFunction1stKind(double x, ComplexNumber z);
ComplexNumber BesselFunction1stKind(ComplexNumber w, ComplexNumber z);
double        BesselFunction2ndKind(double x, double order);
ComplexNumber BesselFunction2ndKind(double x, ComplexNumber z2);
ComplexNumber BesselFunction2ndKind(ComplexNumber w, ComplexNumber z);
long double   sphericalBesselFunction1stKind(double x, double order);
ComplexNumber sphericalBesselFunction1stKind(double x, ComplexNumber z);
ComplexNumber sphericalBesselFunction1stKind(ComplexNumber w, ComplexNumber z);
long double   regularModifiedBesselFunction1stKind(double x, double order);
ComplexNumber regularModifiedBesselFunction1stKind(double x, ComplexNumber z);
ComplexNumber regularModifiedBesselFunction1stKind(ComplexNumber w, ComplexNumber z);
long double   irregularModifiedBesselFunction1stKind(double x, double order);
ComplexNumber irregularModifiedBesselFunction1stKind(double x, ComplexNumber z);
ComplexNumber irregularModifiedBesselFunction1stKind(ComplexNumber w, ComplexNumber z);
double qPochhammerFactorial(double a, double q, int n);
double Ci(double z);
ComplexNumber Ci(ComplexNumber z);
double Cin(double z);
ComplexNumber Cin(ComplexNumber z);
double Chi(double x);
ComplexNumber Chi(ComplexNumber z2);
double Ei(double z);
ComplexNumber Ei(ComplexNumber z);
double Li(double z);
ComplexNumber Li(ComplexNumber z);
double Si(double z);
ComplexNumber Si(ComplexNumber z);
double Shi(double x);
ComplexNumber Shi(ComplexNumber z2);
#endif