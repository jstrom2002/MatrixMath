#pragma once
#include "stdafx.h"
#include "Function.h"

class FunctionMatrix {
public:
	std::vector<std::wstring> element;
	int rows;
	int columns;
	int largestElement;
	int FunctionMatrixPrecision;
	
	FunctionMatrix() {}
	~FunctionMatrix() {}
	FunctionMatrix(int r, int c, std::vector<std::wstring> f) {
		element = f;
		rows = r;
		columns = c;
		largestElement = 0;
		for (int i = 0; i < f.size(); ++i) {
			if ( ((int)f[i].length()) > largestElement) { largestElement = f[i].length(); }
		}
	}

	//class methods
	//=============
	std::wstring get(int i, int j) { return element[i*columns + j]; }

	std::wstring toString() {
		std::wostringstream temp;
		int defaultLength = 2 + largestElement;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				temp << std::setw(defaultLength) << std::left << get(i,j);
			}
			temp << L"\n";
		}
		std::wstring ret = temp.str();
		temp.clear();
		temp.flush();
		return ret;
	}

	void display() { std::wcout << toString() << L"\n"; }
};