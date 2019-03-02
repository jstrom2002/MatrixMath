#ifndef ANALYSIS_H
#define ANALYSIS_H
#pragma once
#include "stdafx.h"

double nthDerivative1D(double D, double x, Function f);
ComplexNumber nthDerivative1D(double D, ComplexNumber x, ComplexFunction f);
double nthDerivative1DEpsilon(double D, double x, double e, Function f);
double nthPartialDerivative(double D, int wrt, std::vector<double> x, Function f);
std::vector<double> GaussLegendreTable();
std::vector<double> GaussHermiteTable();
std::vector<long double> GaussLaguerreTable();
std::vector<double> GaussLaguerreTableWithAlpha();
std::vector<double> GaussKronodTable();
std::vector<double> GaussKronodTableCorrespondingGaussTable();
double integrateSimpsonsRule(double a, double b, Function f);
double integrateGaussLegendreQuadrature(double a, double b, Function f);
double integrateGaussLegendreQuadrature(double a, double b, Function f, std::vector<double> vals, int pos);
double integrateGaussHermiteQuadrature(Function f);
double integrateGaussLaguerreQuadrature(Function f);
double integrateGaussLaguerreQuadrature(Function f, std::vector<double> vals, int pos);
double integrateGaussKronodQuadrature(double a, double b, Function f);
double integrateGaussKronodQuadrature(double a, double b, Function f, std::vector<double> vals, int pos);
std::vector<double> gradient(std::vector<double> X, Function f);
std::vector<double> gradient(std::vector<double> X, VectorValuedFunction vf);
double divergence(std::vector<double> X, Function f);
double divergence(std::vector<double> X, VectorValuedFunction vf);
std::vector<double> directionalDerivative(std::vector<double> X, Function f, std::vector<double> v);
std::vector<double> directionalDerivative(std::vector<double> X, VectorValuedFunction vf, std::vector<double> v);
/*
double compose(double(*f)(double), double(*g)(double), double x);
double* composeNthDVectorValuedFunction(double*(*g)(double[]), double*(*f)(double[]), double X[]);
double* composeNthDParametricVectorValuedFunction(double*(*g)(double[], double), double*(*f)(double[], double), double X[], double t);
double arcLength(double a, double b, double(*funct)(double));
double arcLengthMonteCarlo(double a, double b, double(*funct)(double));
double ODEFirstOrder(double t, double IVP, double(*funct)(double, double));
double ODEFirstOrderImprovedEuler(double t, double IVP, double(*funct)(double, double));
double ODEFirstOrderRK(double t, double IVP, double(*funct)(double, double));
double ODESecondOrder(double t, double IVP_one, double IVP_two, double(*funct)(double, double, double));
double ODESecondOrderRK(double t, double IVP_one, double IVP_two, double(*funct)(double, double, double));
double convolve(double(*funct1)(double), double(*funct2)(double), double t);
*/
#endif