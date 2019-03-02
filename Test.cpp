#pragma once
#include "stdafx.h"

void MatrixTest() {
	Matrix* Mat = NULL;
	int counter,n,m;
	counter = 0;
	n = m = 3;

	//TEST MATRICES
	double temp[9] = {//det = 4, should be diagonalizable
		 2,-1, 0,
		-1, 2,-1,
		 0,-1,2
	};
	double temp2[9] = {//det = -2, should be diagonalizable
			 1, 1, 1,
			 2, 3, 5,
			 4, 6, 8
		};
	double temp3[9] = {//det = 0, not diagonalizable, rank = 2
			 1, 2, 3,
			 4, 5, 6,
			 7, 8, 9
	};
	
	while (counter < 3) {
		if (counter == 0) { 
			Mat = new Matrix(n, m, toSTLVector(temp, n*m)); 
		}
		if (counter == 1) { 
			Mat = new Matrix(n, m, toSTLVector(temp2, n*m)); 
		}
		if (counter == 2) { 
			Mat = new Matrix(n, m, toSTLVector(temp3, n*m)); 
		}

		Matrix* temp;
				
		output.append(L"Main matrix:\r\n");//display main matrix
		output.append((*Mat).toString());
		output.append(L"\r\n");

		output.append(L"Gaussian Elimination:\r\n");
		temp = new Matrix((*Mat).GaussianElimination());
		output.append((*temp).toString());
		output.append(L"\r\n");

		output.append(L"LU decomp\r\nL:\r\n");
		temp = new Matrix((*Mat).lowerTriangularize());
		output.append((*temp).toString());
		output.append(L"\r\n");

		output.append(L"U:\r\n");
		temp = new Matrix((*Mat).upperTriangularize());
		output.append((*temp).toString());
		output.append(L"\r\n");

		output.append(L"determinant = ");
		double d = (*Mat).detExact();
		output.append(to_stringPrecision(d));
		output.append(L"\r\n");

		output.append(L"L*U:\r\n");//make sure L*U = original matrix
		temp = new Matrix((*Mat).multiply((*Mat).lowerTriangularize(), (*Mat).upperTriangularize()));
		output.append((*temp).toString());
		output.append(L"\r\n");

		output.append(L"==============\r\n");
		++counter;
	}
}