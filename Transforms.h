#ifndef TRANSFORMS_H
#define TRANSFORMS_H
#pragma once
//#include "stdafx.h"

double FourierTransformAngularFrequency(double w, Function f);
double FourierTransform(double w, Function f);
double inverseFourierTransform(double w, Function f);
double LaplaceTransform(double s, Function f);
double twoSidedLaplaceTransform(double s, Function f);
double MellinTransform(double s, Function f);
double twoSidedMellinTransform(double s, Function f);
double convolve1D(double t, Function f1, Function f2);
std::vector<ComplexNumber> discreteFourierTransform(std::vector<ComplexNumber> data);
std::vector<ComplexNumber> discreteFourierTransform(std::vector<double> data);
std::vector<ComplexNumber> discreteFourierTransform(std::vector<double> x, Function f);
std::vector<ComplexNumber> discreteFourierTransform(std::vector<ComplexNumber> x, ComplexFunction f);
ComplexMatrix discreteFourierTransform(Matrix A);
ComplexMatrix discreteFourierTransform(ComplexMatrix A);
std::vector<ComplexNumber> inverseDiscreteFourierTransform(std::vector<ComplexNumber> data);
std::vector<ComplexNumber> inverseDiscreteFourierTransform(std::vector<double> data);
ComplexMatrix inverseDiscreteFourierTransform(Matrix A);
ComplexMatrix inverseDiscreteFourierTransform(ComplexMatrix A);
std::vector<ComplexNumber> fastFourierTransform(std::vector<ComplexNumber> data);
void fft(ComplexNumber *data, ComplexNumber *temp, int n, int ndata);
Matrix discreteCosineTransform(std::vector<double> data);
Matrix discreteCosineTransform(Matrix data);
Matrix inverseDiscreteCosineTransform(std::vector<double> X);
Matrix inverseDiscreteCosineTransform(Matrix data);
#endif