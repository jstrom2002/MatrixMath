#ifndef TENSOR_H
#define TENSOR_H
#pragma once
#include "stdafx.h"

class Tensor {
public:
	std::vector<double> element;  //store all tensor elements as row-major indexed vector
	bool isCovariant;			  //if true, covariant.  if false, contravariant.
	std::vector<int> index;		  //list of all dimension sizes for the tensor, i.e. <1,3,2> = 1x3x2 tensor (a 3D column vector with a depth of 2).
	int TensorPrecision;
	int largestElement;
	Tensor(std::vector<int> ind, bool cov);
	Tensor(std::vector<double> el, std::vector<int> ind, bool cov);
	Tensor(std::vector<double> v);
	Tensor();
	~Tensor();
	double get(int n);	
	double get(std::vector<int> rw);
	Matrix getSubmatrix(int ind);
	int rank();
	void randomize();
	Tensor tensorProduct(Tensor A, Tensor B);
	std::wstring toString();
	void display();
};
#endif