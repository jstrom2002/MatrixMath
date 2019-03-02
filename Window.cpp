#pragma once
#include "stdafx.h"

std::vector<ComplexNumber> BartlettWindow(std::vector<ComplexNumber> data) {
	int N = data.size();
	for (int i = 0; i < data.size(); i++) {
		std::complex<double> data2 = data[i].toComplex();
		data2 *= 1 - abs( (i - (0.5*(N-1)))/(0.5*(N+1)) );
		data[i] = ComplexNumber(data2);
	}
	return data;
}

std::vector<ComplexNumber> GaussianWindow(std::vector<ComplexNumber> data) {
	int N = data.size();
	std::complex<double> sigma = mean(data).toComplex();
	for (int i = 0; i < data.size(); i++) {
		std::complex<double> data2 = data[i].toComplex();
		data2 *= exp(-0.5*((i - (0.5*(N - 1))) / sigma));
		data[i] = ComplexNumber(data2);
	}
	return data;
}

std::vector<ComplexNumber> HammingWindow(std::vector<ComplexNumber> data) {
	for (int i = 0; i < data.size(); i++) {
		std::complex<double> data2 = data[i].toComplex();
		data2 *= HammingAlpha + (1.0 - HammingAlpha) * (cos((2 * PI*i) / (data.size() - 1)));
		data[i] = ComplexNumber(data2);
	}
	return data;
}

std::vector<ComplexNumber> WelchWindow(std::vector<ComplexNumber> data) {
	int N = data.size();
	for (int i = 0; i < data.size(); i++) {
		std::complex<double> data2 = data[i].toComplex();
		data2 *= 1 - pow(abs((i - (0.5*(N - 1))) / (0.5*(N - 1))), 2);
		data[i] = ComplexNumber(data2);
	}
	return data;
}