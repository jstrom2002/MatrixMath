#pragma once
#include "stdafx.h"

Tensor::Tensor(std::vector<int> ind, bool cov) {
	int size = 1;
	for (int i = 0; i < ind.size(); ++i) { size *= ind[i]; }
	std::vector<double> vec(size,0);
	element = vec;
	index = ind;
	isCovariant = cov;
	TensorPrecision = precision;
}

Tensor::Tensor(std::vector<double> el, std::vector<int> ind, bool cov) {
		element = el;
		index = ind;
		isCovariant = cov;
		TensorPrecision = precision;
}

Tensor::Tensor(std::vector<double> v) {
		element = v;
		index.push_back(v.size());
		TensorPrecision = precision;
		isCovariant = true;
}

Tensor::Tensor() {}
Tensor::~Tensor() { element.clear(); }


	double Tensor::get(int n) { return element[n]; }

	double Tensor::get(std::vector<int> rw) {//analogous to matrix index [i*columns + j], except there are many indices to consider
		int indx = 0;
		for (int i = 0; i < rw.size() - 1; ++i) {
			indx += (rw[i] * index[i]);
		}
		indx += (rw[rw.size() - 1]);
		return element[indx];
	}

	Matrix Tensor::getSubmatrix(int ind) {
		std::vector<double> elements;
		if (index.size() == 3) {
			for (int i = 0; i < index[0]; ++i) {
				for (int j = 0; j < index[1]; ++j) {
					elements.push_back(get(i*index[1] + j + ind*index[2]));
				}
			}
			return Matrix(index[0], index[1], elements);
		}
		return Matrix();
	}

	int Tensor::rank() { return index.size(); }

	void Tensor::randomize() {
		for (int i = 0; i < element.size(); ++i) {
			srand(time(0) + clock());
			double a = rand()*pow(-1, rand() % 2);
			int sz = rand() % 3;
			while (getLength(a) > sz) { a /= 10; }
			element[i] = a;
		}
	}

	Tensor Tensor::tensorProduct(Tensor A, Tensor B) {
		//finish this later
		return Tensor();
	}

	std::wstring Tensor::toString() {//Uses Einstein notational convention for displaying elements (i.e. no '+' char)
		std::wostringstream s;
		s << L"(";
		for (int i = 0; i < element.size(); ++i) {
			s <<
				std::setprecision(TensorPrecision) <<
				get(i) << L"e";
			if (isCovariant == false) { s << L"_"; }//convention of Ricci calculus
			else { s << L"^"; }
			s << std::to_wstring(i + 1);
			if (i < element.size() - 1) { s << L", "; }
		}
		s << L")";
		std::wstring temp = s.str();
		s.clear();
		return temp;
	}

	void Tensor::display() { std::wcout << toString() << L"\n"; }