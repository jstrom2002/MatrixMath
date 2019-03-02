#pragma once
#include "stdafx.h"

	Matrix::Matrix() { rows = columns = matrixPrecision = 0; }//default constructor
	
	Matrix::~Matrix() { element.clear(); }//destructor

	Matrix::Matrix(int n) {
		rows = n;
		columns = n;
		modulus = 0;
		for (int i = 0; i < n*n; ++i) {
			element.push_back(0);
		}
		matrixPrecision = 4;
		largestElement = 0;
	}

	Matrix::Matrix(int row, int col) {//constructor -- matrix of 0's
		modulus= 0;
		int n = row;
		int m = col;
		std::vector<double> elm(row*col);
		element = elm;
		rows = row; columns = col; matrixPrecision = 4; largestElement = 0;
	}

	Matrix::Matrix(int row, int col, std::vector<double> element2) {//constructor
		int n = row;
		int m = col;
		if (col > row) { n = m; }
		if (col < row) { m = n; }
		element.clear();
		element.resize(row*col,0);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				double n2 = element2[i*col + j];
				if (modulus != 0 && modulus != NULL) { int b = n2; b %= modulus; n2 = b; }
				element[i*col + j] = n2;
				if (getLength(n2) > largestElement) { largestElement = getLength(n2); }
			}
		}
		rows = row; columns = col; matrixPrecision = 4;
	}

	Matrix::Matrix(int row, int col, std::vector<unsigned char> element2) {
		int n = row;
		int m = col;
		modulus = 255;
		if (col > row) { n = m; }
		if (col < row) { m = n; }
		element.clear();
		element.resize(row*col, 0);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				double n2 = element2[i*col + j];
				if (modulus != 0 && modulus != NULL) { int b = n2; b %= modulus; n2 = b; }
				element[i*col + j] = n2;
				if (getLength(n2) > largestElement) { largestElement = getLength(n2); }
			}
		}
		rows = row; columns = col; matrixPrecision = 3;
	}


	Matrix::Matrix(std::vector<std::vector<double>> elm) {
		int n = elm.size();
		int m = elm[0].size();
		element.clear();
		element.resize(elm.size()*elm[0].size(), 0);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				double n2 = elm[i][j];
				if (modulus != 0 && modulus != NULL) { int b = n2; b %= modulus; n2 = b; }
				element[i*m + j] = n2;
				if (getLength(n2) > largestElement) { largestElement = getLength(n2); }
			}
		}
		rows = n; columns = m; matrixPrecision = 4;
	}

	std::wostream& operator<<(std::wostream& os, Matrix rhs) {
		os << rhs.toString();
		return os;
	}


  //Matrix METHODS
  //================
double Matrix::get(int i, int j) { return element[(i*columns) + j]; }
void Matrix::set(int i, int j, double x) { element[(i*columns) + j] = x; }
int Matrix::size() { return rows*columns; }

void Matrix::setModulus(int n) {
	modulus = n;
	if (n == 0) { return; }
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			set(i, j, (int)get(i, j) % n);
		}
	}
}

double Matrix::biggestValue() {
	double x = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (std::abs(get(i, j)) > x) { x = get(i, j); }
		}
	}
	return x;
}

void Matrix::identity() {//turns the matrix into an identity matrix
	for (int i = 0; i<rows; ++i) {
		for (int j = 0; j<columns; ++j) {
			if (i == j) { element[i*columns + j] = 1; }
			else { element[i*columns + j] = 0; }
		}
	}
	largestElement = 0;
	matrixPrecision = 0;
}

Matrix Matrix::submatrix(int i, int j, int i2, int j2) {
	int size1 = (((i2 + 1)*(j2 + 1)) - ((i + 1)*(j + 1)));
	std::vector<double> vals = initZeros(vals, size1);
	int counter = 0;
	for (int a = i; a <= i2; ++a) {
		for (int b = j; b <= j2; ++b) {
			vals[counter] = get(a, b);
			++counter;
		}
	}
	return Matrix(i2 - i + 1, j2 - j + 1, vals);
}

std::vector<double> Matrix::diagonal() {
	std::vector<double> vals;
	for (int i = 0; i < rows; ++i) {
		vals.push_back(get(i, i));
	}
	return vals;
}

Matrix Matrix::multiply(Matrix A, Matrix B) {
//	if (A.modulus != B.modulus) { return Matrix(); }
	int size = A.rows*B.columns;
	std::vector<double> answerN(size);
	if (B.columns == 1) {//if B is a column vector, use a faster method for calculation
		for (int i = 0; i < A.rows; ++i) {
			for (int j = 0; j < A.columns; ++j) {
				answerN[i] += A.get(i, j) * B.get(j, 0);
			}
		}
		return Matrix(A.rows, B.columns, answerN);
	}

	for (int i = 0; i < A.rows; ++i) {
		for (int j = 0; j < B.columns; ++j) {
			for (int k = 0; k < A.columns; ++k) {
				if (modulus == 0 || modulus == 2) { answerN[(i*B.columns) + j] += A.get(i, k) * B.get(k, j); }
			}
		}
	}

	if (modulus == 2) {//if modulus == 2, multiplication is a Boolean product
		Matrix temp(A.rows, A.columns, answerN);
		for (int i = 0; i < size; ++i) {
			if (temp.element[i] != 0) { temp.element[i] = 1; }
		}
		return temp;
	}

	return Matrix(A.rows, B.columns, answerN);
}

Matrix Matrix::BooleanProduct(Matrix A, Matrix B) {
	if (A.modulus != 2 && B.modulus != 2) { return Matrix(); }
	return A.multiply(A, B);//Boolean Product
}

Matrix Matrix::multiply(Matrix A, double B) {
	std::vector<double> answerN(size());
	for (int i = 0; i < A.rows; ++i) {
		for (int j = 0; j < A.columns; ++j) {
			if (modulus == 0) { answerN[(i*A.columns) + j] = A.get(i, j) * B; }
			if (modulus == 2) {
				double n = A.get(i, j) * B;
				int p = n;
				p %= 2;
				answerN[(i*A.columns) + j] = p;
			}
		}
	}
	return Matrix(A.rows, A.columns, answerN);
}

Matrix Matrix::add(Matrix A, Matrix B) {
	if (A.size() != B.size()) { return Matrix(); }
	int size = A.rows*A.columns;
	std::vector<double> answerN = initZeros(answerN, size);
	for (int i = 0; i < A.rows; ++i) {
		for (int j = 0; j < B.columns; ++j) {
			if (modulus == 0) { answerN[(i*A.columns) + j] = A.get(i, j) + B.get(i, j); }
			if (modulus == 2) { answerN[(i*A.columns) + j] = XOR(A.get(i, j), B.get(i, j)); }//boolean matrices treat addition as a XOR operation
		}
	}
	return Matrix(B.rows, A.columns, answerN);
}

Matrix Matrix::add(Matrix A, double B) {
	int size = A.rows*A.columns;
	std::vector<double> answerN = initZeros(answerN, size);
	for (int i = 0; i < A.rows; ++i) {
		for (int j = 0; j < A.columns; ++j) {
			{
				if (modulus == 0) { answerN[(i*A.columns) + j] = A.get(i, j) + B; }
				if (modulus == 2) {
					double n = A.get(i, j) + B;
					int p = n;
					p %= 2;
					answerN[(i*A.columns) + j] = p;
				}
			}
		}
	}
	return Matrix(A.rows, A.columns, answerN);
}

Matrix Matrix::subtract(Matrix A, Matrix B) {
	if (A.size() != B.size()) { return Matrix(); }
	int size = A.rows*A.columns;
	std::vector<double> answerN = initZeros(answerN, size);
	for (int i = 0; i < size; ++i) { answerN[i] = 0; }

	for (int i = 0; i < A.rows; ++i) {
		for (int j = 0; j < B.columns; ++j) {
			{
				if (modulus == 0) { answerN[(i*A.columns) + j] += A.get(i, j) - B.get(i, j); }
			}
		}
	}
	return Matrix(B.rows, A.columns, answerN);
}

Matrix Matrix::subtract(Matrix A, double B) {
	int size = A.rows*A.columns;
	std::vector<double> answerN = initZeros(answerN, size);
	for (int i = 0; i < A.rows; ++i) {
		for (int j = 0; j < A.columns; ++j) {
			{
				answerN[(i*A.columns) + j] = A.get(i, j) - B;
			}
		}
	}
	return Matrix(A.rows, A.columns, answerN);
}

Matrix Matrix::inverseExact() {
	if (rows != columns) { return Matrix(); }
	Matrix m(rows, columns, element);
	m = m.concatenate(m, identityMatrix(rows,columns));
	m = m.GaussianElimination();
	while (m.columns > columns) {m.removeColumn(0);}
	return m;
}


Matrix Matrix::inverse() {
	/*INPUT: n×m matrix A.
	OUTPUT: n×m matrix in reduced row echelon form.
	1.Set j←1
	2.For each row i from 1 to n do
	a. While column j h as all zero elements, set j <- j+1. If j>m return A.
	b. If element a_ij is zero, then interchange row i with a row x > i that has a_xj!=0.
	c. Divide each element of row i by a_ij, thus making the pivot aij equal to one.
	d. For each row k from 1 to n, with k != i, subtract row i multiplied by a_kj from row k.
	3. Return transformed matrix A.		*/

	if (modulus == 2) {
		Matrix A = *this;
		Matrix AT = *this;
		AT.transpose();
		Matrix I(rows, columns);
		I.identity();
		A = multiply(A, AT);
		A = subtract(A, I);
		if (A.sum() == 0) { return AT; }
		else { return Matrix(); }//if A*AT = I, then this matrix has an inverse (which is AT).
	}

	if (det() == 0) { return Matrix(); }
	Matrix A(rows, columns, element);
	Matrix Aug(rows, columns);
	Aug.identity();

	while (A.zeroInDiag() == true) {
		for (int i = 0; i < rows; ++i) {
			if (A.get(i, i) == 0) {
				for (int j = 0; j < columns; ++j) {
					if (A.get(i, j) != 0 && A.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create a non-zero diagonal
						A.swapRows(i, j);
						Aug.swapRows(i, j);
						j = columns;
					}
				}
			}
		}
	}

	for (int j = 0; j<columns; ++j) {//go down and eliminate the lower triangle of values
		double pivot = A.get(j, j);
		for (int i = 0; i < rows; ++i) {
			if (i != j && A.get(i, j) != 0) {
				std::vector<double> temp = A.row(j);
				std::vector<double> temp2 = Aug.row(j);

				if (abs(A.get(i, j)) == abs(pivot)) {
					if (A.get(i, j) == pivot) {
						arrmultiply(temp, -1.0);
						arrmultiply(temp2, -1.0);
						goto skip;
					}
					if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) { goto skip; }
				}
				if (abs(A.get(i, j)) != abs(pivot)) {
					double a = -A.get(i, j) / pivot;
					arrmultiply(temp, a);
					arrmultiply(temp2, a);
				}

			skip:
				A.addToRow(i, temp);
				Aug.addToRow(i, temp2);

			}
		}
	}

	for (int j = columns - 1; j >= 0; --j) {//go the other way to remove the upper triangle of non-diagonal values
		double pivot = A.get(j, j);
		for (int i = rows - 1; i >= 0; --i) {
			if (i != j && A.get(i, j) != 0) {
				std::vector<double> temp = A.row(j);
				std::vector<double> temp2 = Aug.row(j);

				if (abs(A.get(i, j)) == abs(pivot)) {
					if (A.get(i, j) == pivot) {
						arrmultiply(temp, -1);
						arrmultiply(temp2, -1);
						goto skip2;
					}
					if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) { goto skip2; }
				}
				if (abs(A.get(i, j)) != abs(pivot)) {
					double a = -A.get(i, j) / pivot;
					arrmultiply(temp, a);
					arrmultiply(temp2, a);
				}

			skip2:
				A.addToRow(i, temp);
				Aug.addToRow(i, temp2);
			}
		}
	}

	for (int i = 0; i < rows; ++i) {//make diagonals into all 1's
		if (A.get(i, i) != 0 && A.get(i, i) != 1) {
			double a = ((double)1 / (A.get(i, i)));
			A.multiplyRow(i, a);
			Aug.multiplyRow(i, a);
		}
	}
	return Matrix(rows, columns, Aug.element);
}

Matrix Matrix::exponent(Matrix A, int b) {
	if (b != floor(b)) { return Matrix(); }//temporary }//fractional exponent handling
	if (b < 0) {
		A = A.inverse();
		b = abs(b);
	}
	//std::wcout << L"\nMADE IT HERE\n";
	if (b == 0) {
		Matrix A2(A.rows, A.columns);
		A2.identity();
		return A2;
	}
	for (int i = 1; i < b; ++i) { A = A.multiply(A, A); }
	return A;
}

//Overloaded operators:
Matrix Matrix::operator*=(Matrix rhs) { *this = multiply(*this, rhs); return*this; }
Matrix Matrix::operator*=(double x) { *this = multiply(*this, x); return *this; }
Matrix Matrix::operator*(Matrix rhs) { return multiply(*this, rhs); }
Matrix Matrix::operator*(double x) { return multiply(*this, x); }

Matrix Matrix::operator+=(Matrix rhs) { *this = add(*this, rhs); return*this; }
Matrix Matrix::operator+=(double x) { *this = add(*this, x); return *this; }
Matrix Matrix::operator+(Matrix rhs) { return add(*this, rhs); }
Matrix Matrix::operator+(double x) { return add(*this, x); }

Matrix Matrix::operator-=(Matrix rhs) { *this = subtract(*this, rhs); return*this; }
Matrix Matrix::operator-=(double x) { *this = subtract(*this, x); return *this; }
Matrix Matrix::operator-(Matrix rhs) { return subtract(*this, rhs); }
Matrix Matrix::operator-(double x) { return subtract(*this, x); }

Matrix Matrix::operator^(double x) { return exponent(*this, x); }

bool Matrix::operator==(Matrix B) {
	if (rows != B.rows || columns != B.columns) { return false; }
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (get(i, j) != B.get(i, j)) { return false; }
		}
	}
	return true;
}
//========================

Matrix Matrix::outerProduct(Matrix A, Matrix B) {//tensore product that creates for every element of Matrix A, a submatrix of B * (A_i_j)
	if (A.columns != B.rows) { return Matrix(); }
	std::vector<double> ans(A.size()*B.size());
	int rw = A.rows*B.rows;
	int cl = A.columns*B.columns;
	
	for (int i = 0; i < rw; ++i) {
		for (int j = 0; j < cl; ++j) {
			ans[i*cl + j] = B.element[((i%B.rows)*B.columns) + (j%B.columns)] * A.element[floor(i / B.rows)*A.columns + floor(j / B.columns)];
		}
	}
	return Matrix(rw, cl, ans);
}

Matrix Matrix::HadamardProduct(Matrix A, Matrix B) {//multiplys two matrices by index so that C(i,j) = A(i,j)*B(i,j)
	if (A.rows != B.rows || A.columns != B.columns) { return Matrix(); }
	std::vector<double> elm(A.size());
	for (int i = 0; i <A.rows; ++i) {
		for (int j = 0; j < A.columns; ++j) {
			elm[i*A.columns + j] = A.get(i, j)*B.get(i, j);
		}
	}
	return Matrix(A.rows, A.columns, elm);
}

std::vector<double> Matrix::row(int rw) {
	if (rw > rows || rw < 0) { std::vector<double> n; return n; }
	std::vector<double> x = initZeros(x, columns);
	for (int k = 0; k < columns; ++k)
	{
		x[k] = get(rw, k);
	}
	return x;
}

std::vector<double> Matrix::column(int cl) {//returns a row vector
	if (cl > columns || cl<0) { std::vector<double> n; return n; }
	std::vector<double> x = initZeros(x, rows);
	for (int k = 0; k < rows; ++k) {
		x[k] = get(k, cl);
	}
	return x;
}

void Matrix::setRowNonZeroValues(int rw, double x) {
	for (int k = 0; k < columns; ++k) {
		if (abs(get(rw,k)) > SYSTEM_EPSILON) {
			set(rw, k, x);
		}
	}
}

void Matrix::setRow(int rw, std::vector<double> n) {
	for (int k = 0; k < columns; ++k) { set(rw, k, n[k]); }
}

void Matrix::setColumnNonZeroValues(int col, double x) {
	for (int k = 0; k < rows; ++k) { 
		if (abs(get(k, col)) > SYSTEM_EPSILON) {
			set(k, col, x);
		}
	}
}

void Matrix::setColumn(int col, std::vector<double> n) {
	for (int k = 0; k < rows; ++k) { set(k, col, n[k]); }
}

void Matrix::multiplyRow(int rw, double n) {
	for (int k = 0; k < columns; ++k)
	{
		element[(rw*columns) + k] *= n;
	}
}

void Matrix::multiplyRow(int rw, std::vector<double> n) {
	for (int k = 0; k < columns; ++k)
	{
		element[(rw*columns) + k] *= n[k];
	}
}

void Matrix::addToRow(int rw, double n) {
	for (int k = 0; k < columns; ++k)
	{
		element[(rw*columns) + k] += n;
	}
}

void Matrix::addToRow(int rw, std::vector<double> n) {
	for (int k = 0; k < columns; ++k)
	{
		element[(rw*columns) + k] += n[k];
	}
}

void Matrix::addToColumn(int cl, double n) {
	for (int k = 0; k < columns; ++k)
	{
		element[(k*columns) + cl] += n;
	}
}

void Matrix::addToColumn(int cl, std::vector<double> n) {
	for (int k = 0; k < columns; ++k)
	{
		element[(k*columns) + cl] += n[k];
	}
}

std::vector<double> Matrix::findRow(std::vector<double> x) {
	int index = 0;
	bool isRow = false;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (x[j] != get(i, j)) { isRow = false; j = columns; }
			if (j = columns - 1) { isRow = true; index = i; }
		}
	}
	return row(index);
}

std::vector<double> Matrix::findColumn(std::vector<double> x) {
	int index = -1;
	bool isColumn = false;
	for (int j = 0; j < columns; ++j) {
		for (int i = 0; i < rows; ++i) {
			if (x[j] != get(i, j)) { isColumn = false; i = rows; }
			if (i = rows - 1) { isColumn = true; index = j; }
		}
	}
	if (index > -1) { return column(index); }
	else { std::vector<double> n; return n; }
}

void Matrix::reverseRow(int n) {
	std::vector<double> vec = row(n);
	vec = reverse(vec);
	for (int j = 0; j < columns; ++j) {
		set(n, j, vec[j]);
	}
}

void Matrix::reverseColumn(int n) {//IN PROGRESS
	std::vector<double> vec = column(n);
	vec = reverse(vec);
	for (int j = 0; j < rows; ++j) {
		set(j, n, vec[j]);
	}
}

void Matrix::swapRows(int rw1, int rw2) {
	std::vector<double> r1 = row(rw1);
	std::vector<double> r2 = row(rw2);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (i == rw1) { element[i*columns + j] = r2[j]; }
			if (i == rw2) { element[i*columns + j] = r1[j]; }
			if (i > rw1 && i > rw2) { return; }
		}
	}
}

void Matrix::swapColumns(int col1, int col2) {//IN PROGRESS
	std::vector<double> c1 = column(col1);
	std::vector<double> c2 = column(col2);
	for (int j = 0; j < columns; ++j) {
		for (int i = 0; i < rows; ++i) {
			if (j == col1) { element[i*columns + j] = c2[i]; }
			if (j == col2) { element[i*columns + j] = c1[i]; }
			if (j > col1 && j > col2) { return; }
		}
	}
}

void Matrix::randomize() {//creates a random matrix of double values
	srand(time(0) + clock());
	largestElement = 0;
	for (int i = 0; i<rows; ++i) {
		for (int j = 0; j<columns; ++j) {
			double a = ((double)rand() / (pow(10, (rand() % 5) + 3))) *  pow(-1, rand() % 2);
			while (a == 0) { ((double)rand() / (pow(10, (rand() % 5) + 3))) *  pow(-1, rand() % 2); }
			element[i*columns + j] = a;
			if (getLength(element[i*columns + j]) > largestElement) { largestElement = getLength(element[i*columns + j]); }
		}
	}
	matrixPrecision = 4;
}

void Matrix::randomizeInteger() {//creates a random matrix of integer values
	srand(time(0) + clock());
	largestElement = 0;
	for (int i = 0; i<rows; ++i) {
		for (int j = 0; j<columns; ++j) {
			double a = (int)((double)rand() / (pow(10, (rand() % 4) + 2))) * pow(-1, rand() % 2);
			element[i*columns + j] = a;
			if (getLength(a) > largestElement) { largestElement = getLength(a); }
		}
	}
	matrixPrecision = 0;
}

void Matrix::randomizeBoolean() {//creates a random matrix of boolean values
	srand(time(0) + clock());
	largestElement = 0;
	for (int i = 0; i<rows; ++i) {
		for (int j = 0; j<columns; ++j) {
			double a = 1 * (rand() % 2);
			element[i*columns + j] = a;
			if (getLength(a) > largestElement) { largestElement = getLength(a); }
		}
	}
	matrixPrecision = 0;
	modulus = 2;
}

void Matrix::randomizeSymmetric() {
	srand(time(0) + clock());
	largestElement = 0;
	for (int i = 0; i<rows; ++i) {
		for (int j = 0; j <= i; ++j) {
			element[i*columns + j] = ((double)rand() / (30000)) *  pow(-1, rand() % 2);
			element[j*columns + i] = element[i*columns + j];
			if (getLength(element[i*columns + j]) > largestElement) { largestElement = getLength(element[i*columns + j]); }
		}
	}
	matrixPrecision = 4;
}

void Matrix::sortRowsSmallToLarge() {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			//NEEDS WORK
		}
	}
}

Matrix Matrix::expandToUpperLeft(int r, int c) {//takes original matrix and creates a new matrix with the original in the lower right corner.
	Matrix newEl(r, c);
	if (r > rows && c > columns) {
		newEl.identity();
		for (int i = r - rows; i < r; ++i) {
			for (int j = c - columns; j < c; ++j) {
				newEl.set(i, j, get((i - (r - rows)), (j - (c - columns))));
			}
		}
	}
	return Matrix(r, c, newEl.element);
}

Matrix Matrix::expandToLowerRight(int r, int c) {//takes original matrix and creates a new matrix with the original in the upper left corner.
	Matrix newEl(r, c);
	if (r > rows && c > columns) {
		newEl.identity();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				newEl.set(i, j, get(i, j));
			}
		}
	}
	return Matrix(r, c, newEl.element);
}

Matrix Matrix::extendRows(int r) {
	Matrix newEl(rows + r, columns);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			newEl.set(i, j, get(i, j));
		}
	}
	return Matrix(rows + r, columns, newEl.element);
}

Matrix Matrix::extendColumns(int c) {
	Matrix newEl(rows, columns + c);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			newEl.set(i, j, get(i, j));
		}
	}
	return Matrix(rows, columns + c, newEl.element);
}

Matrix Matrix::concatenate(Matrix A, Matrix B) {//creates a larger matrix from two matrices, where each matrix is next to the other one
	int r = A.rows;
	if (B.rows > A.rows) { r = B.rows; }
	int c = A.columns + B.columns;
	Matrix C(r, c);
	for (int i = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j) {
			if (j < A.columns && i<A.rows) { C.set(i, j, A.get(i, j)); }
			if (j >= A.columns && i<B.rows) { C.set(i, j, B.get(i, j - A.columns)); }
		}
	}

	return C;
}

Matrix Matrix::autoCorrelation(Matrix A) {return A.convolve(A, A);}//auto-correlation is the convolution of a signal with its own conjugate, which is the function itself in a Real valued domain

Matrix Matrix::crossCorrelation(Matrix A, Matrix B) {
	int kCenterX = A.columns / 2;
	int kCenterY = A.rows / 2;

	int mm, nn, ii, jj;

	std::vector<double> vals(B.element.size());

	for (int i = 0; i < B.rows; ++i) {
		for (int j = 0; j < B.columns; ++j) {
			for (int m = 0; m < A.rows; ++m) {
				for (int n = 0; n < A.columns; ++n) {

					// index of input signal, used for checking boundary
					ii = i + (kCenterY - m);
					jj = j + (kCenterX - n);

					// ignore input samples which are out of bound
					if (ii >= 0 && ii < B.rows && jj >= 0 && jj < B.columns) { vals[i*B.columns + j] += B.get(ii, jj) * A.get(m, n); }
				}
			}
		}
	}
	return Matrix(B.rows, B.columns, vals);
}


Matrix Matrix::convolve(Matrix A, Matrix B) {//covolve Matrix kernel A w/ Matrix B
	int kCenterX = A.columns / 2;
	int kCenterY = A.rows / 2;

	int mm, nn, ii, jj;

	std::vector<double> vals(B.element.size());

	for (int i = 0; i < B.rows; ++i){
		for (int j = 0; j < B.columns; ++j)   {
			for (int m = 0; m < A.rows; ++m){
				mm = A.rows - 1 - m;      // row index of flipped kernel
				for (int n = 0; n < A.columns; ++n){
					nn = A.columns - 1 - n;  // column index of flipped kernel

					// index of input signal, used for checking boundary
					ii = i + (kCenterY - mm);
					jj = j + (kCenterX - nn);

					// ignore input samples which are out of bound
					if (ii >= 0 && ii < B.rows && jj >= 0 && jj < B.columns) { vals[i*B.columns + j] += B.get(ii,jj) * A.get(mm,nn); }
				}
			}
		}
	}
	return Matrix(B.rows, B.columns, vals);
}

//BOOLEAN FUNCTIONS
//=================
bool Matrix::isOrderedSmallToLarge() {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			//NEEDS WORK
		}
	}
	return true;
}

bool Matrix::isIntegerMatrix() {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (floor(abs(get(i, j))) != abs(get(i, j))) { return false; }
		}
	}
	matrixPrecision = 0;
	return true;
}

bool Matrix::isSymmetric() {
	if (*this == transpose()) { return true; }
	return false;
}

bool Matrix::isBoolean() {
	if (modulus == 2) { return true; }
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (get(i, j) != 0 && get(i, j) != 1) { return false; }
		}
	}
	return true;
}

bool Matrix::isDiagonalized() {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (i != j) {
				if (abs(get(i, j)) > 0.001) { return false; }//only 4 digits of the mantissa are required for matrixPrecision
			}
		}
	}
	return true;
}

bool Matrix::isOrthonormal() {
	for (int i = 0; i < columns; ++i) {
		if (columnNorm(i) != 1) { return false; }
	}
	if (det() == 0) {
		std::wcout << det() << std::endl;
		return false;
	}
	return true;
}

bool Matrix::isLinearlyIndependent() {
	if (rows < columns) { return true; }
	if (det() == 0) { return false; }
	if (GramMatrix().det() == 0) { return false; }
	return true;
}

std::vector<double> Matrix::meanVector() {//returns a vector of all of the column means
	std::vector<double> vec;
	for (int i = 0; i < columns; ++i) {
		vec.push_back(columnMean(i));
	}
	return vec;
}

Matrix Matrix::varianceMatrix() {
	Matrix M(columns, columns);
	for (int i = 0; i < columns; ++i) {
		for (int j = 0; j < rows; ++j) {
			//p[i*rows + j] = get(j, i);
		}
	}
	return Matrix();
}

Matrix Matrix::populationCovarianceMatrix() {
	Matrix M(columns, columns);
	for (int i = 0; i < columns; ++i) {
		for (int j = 0; j < columns; ++j) {
			std::vector<double> X1 = arrsubtract(column(i), columnMean(i));
			std::vector<double> X2 = arrsubtract(column(j), columnMean(j));
			M.set(i, j, dotProduct(X1, X2) / (rows));
		}
	}
	return M;
}

Matrix Matrix::sampleCovarianceMatrix() {
	Matrix M(columns, columns);
	for (int i = 0; i < columns; ++i) {
		for (int j = 0; j < columns; ++j) {
			std::vector<double> X1 = arrsubtract(column(i), columnMean(i));
			std::vector<double> X2 = arrsubtract(column(j), columnMean(j));
			M.set(i, j, dotProduct(X1, X2) / ((rows)-1));
		}
	}
	return M;
}


Matrix Matrix::transpose() {
	int n = rows;
	int m = columns;
	std::vector<double> p(n*m);
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			p[i*rows + j] = get(j, i);
		}
	}
	return Matrix(m, n, p);
}

Matrix Matrix::GramMatrix() {//turn row vector or nxn matrix into its Gram matrix
	Matrix T = transpose();
	if (columns == 1) { return multiply(*this, T); }
	else { return multiply(T, *this); }
}

Matrix Matrix::GramMatrix(Matrix M) {//turn 2 vectors into a Gram matrix
	Matrix Trans = M.transpose();
	if (M.columns == 1) { return multiply(M, Trans); }
	else { return multiply(Trans, M); }
}

Matrix Matrix::Vandermonde() {//create Vandermonde matrix from vector
	int p = columns;
	if (columns == 1) { p = rows; }
	std::vector<double> temp = initZeros(temp, p*p);
	for (int i = 0; i < p; ++i) {
		for (int j = 0; j < p; ++j) {
			temp[(i*columns) + j] = pow(element[i], j);
		}
	}
	return Matrix(p, p, temp);
}

double Matrix::mean() {
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			answer += get(i, j);
		}
	}
	return answer / size();
}
double Matrix::sum() { double answer = 0;  for (int b = 0; b < rows*columns; ++b) { answer += element[b]; } return answer; }
double Matrix::sumSquares() { double answer = 0;  for (int b = 0; b < rows*columns; ++b) { answer += element[b] * element[b]; } return answer; }

double Matrix::geometricMean() {
	double answer = 1;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			answer *= get(i, j);
		}
	}
	return pow(answer, (1.0 / size()));
}

double Matrix::trace() {//the trace of a matrix is the sum of its diagonals
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += get(i, i);
	}
	return answer;
}

std::wstring Matrix::toString() {
	std::wostringstream s;
	int defaultLength = 4;
	bool tf = isIntegerMatrix();
	if (tf == true) { --defaultLength; }
	if (tf == false) { ++defaultLength; }
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			s <<
				std::fixed <<
				std::left <<
				std::setw(defaultLength + (matrixPrecision)+(largestElement)) <<
				std::setprecision(matrixPrecision) <<
				element[i*columns + j];
		}
		s << L"\r\n";
	}
	std::wstring temp = s.str();
	s.clear();
	return temp;
}

void Matrix::printToFile(std::wstring filename) {
	std::wofstream outf(filename, std::ios_base::app);
	outf << toString() << L"\r\n";
	outf.close();
}

void Matrix::display() {
	std::wcout << toString() << std::endl;
}

double Matrix::norm() {//the Euclidean or Schur norm of the matrix, equal to the square root of the sum of all squared elements.
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			answer += pow(get(i, j), 2);
		}
	}
	return sqrt(answer);
}

double Matrix::FrobeniusNorm() { //the Froebenius norm is equal to the sqrt of the trace of the matrix multiplied by its conjugate transpose.
	Matrix AT = *this;
	AT.transpose();
	return sqrt(multiply(*this, AT).trace());
}

double Matrix::pNorm(double p) {
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			answer += pow(get(i, j), p);
		}
	}
	return sqrt(answer);
}

double Matrix::sumRow(int r) {
	double answer = 0;
	for (int i = 0; i < columns; ++i) {
		answer += get(r, i);
	}
	return answer;
}

double Matrix::sumColumn(int r) {
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += get(i, r);
	}
	return answer;
}

double Matrix::sumAll() {
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; j++) {
			answer += get(i, j);
		}
	}
	return answer;
}

double Matrix::columnMean(int c) {
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += get(i, c);
	}
	return (answer / rows);
}

double Matrix::rowMean(int r) {
	double answer = 0;
	for (int i = 0; i < columns; ++i) {
		answer += get(r, i);
	}
	return (answer / columns);
}

double Matrix::rowNorm(int r) {
	double answer = 0;
	for (int i = 0; i < columns; ++i) {
		answer += pow(get(r, i), 2);
	}
	return sqrt(answer);
}

double Matrix::columnNorm(int c) {
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += pow(get(i, c), 2);
	}
	return sqrt(answer);
}

double Matrix::columnSumSquares(int c) {
	double colMean = columnMean(c);
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += pow(get(i, c) - colMean, 2);
	}
	return answer;
}

double Matrix::rowSumSquares(int r) {
	double rMean = rowMean(r);
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += pow(get(r, i) - rMean, 2);
	}
	return answer;
}

double Matrix::columnVariance(int c) {
	double colMean = columnMean(c);
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += pow(get(i, c) - colMean, 2);
	}
	return answer / (rows - 1);
}

double Matrix::columnCovariance(int c1, int c2) {
	int n = rows;
	double answer = 0;
	std::vector<double> x = column(c1);
	std::vector<double> y = column(c2);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			answer += 0.5 * (x[i] - x[j]) * (y[i] - y[j]);
		}
	}
	return answer / (n*n);
}

double Matrix::columnDeviation(int c) {
	double colMean = columnMean(c);
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += pow(get(i, c) - colMean, 2);
	}
	return sqrt(answer / (rows - 1));
}

double Matrix::rowCovariance(int r1, int r2) {
	int n = columns;
	double answer = 0;
	std::vector<double> x = row(r1);
	std::vector<double> y = row(r2);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			answer += 0.5 * (x[i] - x[j]) * (y[i] - y[j]);
		}
	}
	return answer / (n*n);
}

double Matrix::rowVariance(int r) {
	double rMean = rowMean(r);
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += pow(get(r, i) - rMean, 2);
	}
	return answer / (columns - 1);
}

double Matrix::rowDeviation(int r) {
	double rMean = rowMean(r);
	double answer = 0;
	for (int i = 0; i < rows; ++i) {
		answer += pow(get(r, i) - rMean, 2);
	}
	return sqrt(answer / (columns - 1));
}

std::vector<double> Matrix::normalizedRow(int r) {
	std::vector<double> rw = row(r);
	rw = arrmultiply(rw, rowNorm(r));
	return rw;
}

std::vector<double> Matrix::normalizedColumn(int c) {
	std::vector<double> col = row(c);
	col = arrmultiply(col, columnNorm(c));
	return col;
}

void Matrix::normalizeColumn(int c) {
	double temp = columnNorm(c);
	temp = 1 / temp;
	for (int i = 0; i < rows; ++i) {
		element[i*columns + c] *= temp;
	}
}

void Matrix::normalizeRow(int r) {
	double temp = rowNorm(r);
	temp = 1 / temp;
	for (int i = 0; i < columns; ++i) {
		element[r*columns + i] *= temp;
	}
}

void Matrix::normalize() {
	double d = norm();
	if (d != 0) {
		d = (1 / d);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; j++) {
				element[i*columns + j] *= d;
			}
		}
	}
}

double Matrix::sampleStandardDeviation() {
	double answer = 0;
	double Mean = mean();
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			answer += pow(get(i, j) - Mean, 2);
		}
	}
	return sqrt(answer / ((rows*columns) - 1));
}

double Matrix::populationStandardDeviation() {
	double answer = 0;
	double Mean = mean();
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			answer += pow(get(i, j) - Mean, 2);
		}
	}
	return sqrt(answer / (rows*columns));
}

double Matrix::ChiSquareTestStatistic() {
	double answer = 0;
	double df = (rows - 1) * (columns - 1);

	double sm = sumAll();
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			double Ei = (sumRow(i)*sumColumn(j) / sm);
			answer += pow(get(i, j) - Ei, 2) / Ei;
		}
	}
	return answer;
}

double Matrix::ChiSquareDegreesFreedom() { return  (rows - 1) * (columns - 1); }

void Matrix::removeRow(int n) {
	std::vector<double> newElements = initZeros(newElements, (rows-1)*columns);
	int counter = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (i != n) { newElements[counter] = get(i, j); ++counter; }
		}
	}
	element = newElements;
	--rows;
}

void Matrix::removeColumn(int n) {
	std::vector<double> newElements = initZeros(newElements, rows*(columns-1));
	int counter = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (j != n) { newElements[counter] = get(i, j); ++counter; }
		}
	}
	element = newElements;
	--columns;
}

void Matrix::trim() {//removes row and column vectors of 0's from inner, non-zero, matrix
	int firstNonZeroRow = 0;
	int firstNonZeroColumn = 0;

	//remove inner column/row 0 vectors
	for (int i = rows - 1; i >= 0; --i) {
		if (sumRow(i) == 0) { removeRow(i); }
		if (sumRow(i) != 0) { firstNonZeroRow = i; }
	}
	for (int i = columns - 1; i >= 0; --i) {
		if (sumColumn(i) == 0) { removeColumn(i); }
		if (sumColumn(i) != 0) { firstNonZeroColumn = i; }
	}

	//remove outer column/row 0 vectors;
	for (int i = 0; i < firstNonZeroRow; ++i) {
		if (sumRow(i) == 0) { removeRow(i); }
	}
	for (int i = 0; i < firstNonZeroColumn; ++i) {
		if (sumColumn(i) == 0) { removeColumn(i); }
	}
}

Matrix Matrix::directSum(Matrix A, Matrix B) {
	Matrix T(A.rows + B.rows, A.columns + B.columns);
	T.identity();
	for (int i = 0; i < T.rows; ++i) {
		for (int j = 0; j < T.columns; ++j) {
			if (i < A.rows && j < A.columns) { T.set(i, j, A.get(i, j)); }
			if (i >= A.rows && j >= A.columns) { T.set(i, j, B.get(i - A.rows, j - A.columns)); }
		}
	}
	return T;
}

Matrix Matrix::tensorProduct(Matrix A, Matrix B) {
	Matrix T(A.rows * B.rows, A.columns * B.columns);
	T.identity();
	int Rbound = (A.rows - 1);
	int Cbound = (A.columns - 1);
	int Rcounter = 0;
	int Ccounter = 0;
	for (int i = 0; i < T.rows; ++i) {
		for (int j = 0; j < T.columns; ++j) {
			T.set(i, j, A.get(Rcounter, Ccounter) * B.get(i%B.rows, j%B.columns));
			if (j % A.columns == Cbound && j != 0) { ++Ccounter; }
		}
		Ccounter = 0;
		if (i % A.rows == Rbound && i != 0) { ++Rcounter; }
	}
	return T;
}

bool Matrix::canSwapOutZeroDiagonals() {
	for (int i = 0; i < columns; ++i) {
		if (sumColumn(i) == 0) { return false; }
	}
	return true;
}

Matrix Matrix::swapOutZeroDiagonals() {
	Matrix M(rows, columns, element);
	while (M.zeroInDiag() == true) {
		for (int i = 0; i < rows; ++i) {
			if (M.get(i, i) == 0) {
				for (int j = 0; j < columns; ++j) {
					if (M.get(i, j) != 0 && M.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create M non-zero diagonal
						M.swapRows(i, j); j = columns;
					}
				}
			}
		}
	}
	return M;
}

double Matrix::columnAbsMax(int n) {
	double m = get(0,n);
	if (rows <= 1) { return m; }
	for (int i = 0; i < rows; ++i) {
		double temp = (get(i, n));
		if (abs(temp) > m) { m = temp; }
	}
	return m;
}

double Matrix::columnMax(int n) {
	double m = get(0,n);
	if (rows <= 1) { return m; }
	for (int i = 0; i < rows; ++i) {
		double temp = get(i, n);
		if (temp > m) { m = temp; }
	}
	return m;
}

double Matrix::columnAbsMin(int n) {
	double m = get(0, n);
	if (rows <= 1) { return m; }
	for (int i = 1; i < rows; ++i) {
		double temp = (get(i, n));
		if (abs(temp) < m) { m = temp; }
	}
	return m;
}

double Matrix::columnMin(int n) {
	double m = get(0, n);
	if (rows <= 1) { return m; }
	for (int i = 1; i < rows; ++i) {
		double temp = get(i, n);
		if (temp < m) { m = temp; }
	}
	return m;
}

double Matrix::rowAbsMax(int n) {
	double m = get(n,0);
	if (columns <= 1) { return m; }
	for (int i = 0; i < columns; ++i) {
		double temp = (get(n,i));
		if (abs(temp) > m) { m = temp; }
	}
	return m;
}

double Matrix::rowMax(int n) {
	double m = get(n,0);
	if (columns <= 1) { return m; }
	for (int i = 0; i < columns; ++i) {
		double temp = get(n,i);
		if (temp > m) { m = temp; }
	}
	return m;
}

double Matrix::rowAbsMin(int n) {
	double m = get(n, 0);
	if (columns <= 1) { return m; }
	for (int i = 0; i < columns; ++i) {
		double temp = (get(n,i));
		if (abs(temp) < m) { m = temp; }
	}
	return m;
}

double Matrix::rowMin(int n) {
	double m = get(n, 0);
	if (columns <= 1) { return m; }
	for (int i = 0; i < columns; ++i) {
		double temp = get(n,i);
		if (temp < m) { m = temp; }
	}
	return m;
}

int Matrix::rowNonzeroValues(int rw) {
	int n = 0; 
	for (int i = 0; i < columns; ++i) {
		if (abs(get(rw,i)) > SYSTEM_EPSILON) {
			++n;
		}
	}
	return n;
}

int Matrix::columnNonzeroValues(int cl) {
	int n = 0;
	for (int i = 0; i < rows; ++i) {
		if (abs(get(i,cl)) > SYSTEM_EPSILON) {
			++n;
		}
	}
	return n;
}

int Matrix::getPivot(int rw) {
	int val = 0;
	while (val<columns){
		if (abs(get(rw, val)) > SYSTEM_EPSILON) {
			return val;
		}
		++val; 
	}
	return -1;
}

int Matrix::getReversePivot(int rw) {
	int val = columns-1;
	while (val >= 0) {
		if (abs(get(rw, val)) > SYSTEM_EPSILON) {
			return val;
		}
		--val;
	}
	return -1;
}

Matrix Matrix::GaussianElimination() {
	int i, j, k;	
	double piv;
	Matrix m(rows, columns,element);
	//Gauss-Jordan elimination
	for (i = 0; i < rows; i++) {
		if (zeroInDiag() && canSwapOutZeroDiagonals()) { m = m.swapOutZeroDiagonals(); }//remove zeroes in diagonal if possible
		int ind = m.getPivot(i);
		if (ind >= 0) {
			piv = m.get(i, ind);
			//normalize pivot row so pivot = 1
			if (piv != 1 && piv != 0 && abs(piv) > SYSTEM_EPSILON) {
				for (int l = i; l < rows; ++l) {
					m.set(i, l, m.get(i, l) / piv);
				}
				m.set(i, i, 1);
			}
			//proceed down the column to make each non-pivot value = 0
			for (j = 0; j < rows; j++) {
				if (m.get(j, i) != 0) {
					if (j != i) {
						piv = m.get(j, i);
						for (k = i + 1; k < columns; k++) {
							m.set(j, k, m.get(j, k) - (m.get(i, k) * piv));
						}
						m.set(j, i, 0);
					}
				}
			}
		}
	}
	return m;
}

Matrix Matrix::lowerTriangularize() {
	int i, j, k;
	double piv;
	Matrix m(rows, columns, element);
	//Gauss-Jordan elimination
	for (i = rows-1; i >= 0; --i) {
		if (zeroInDiag() && canSwapOutZeroDiagonals()) { m = m.swapOutZeroDiagonals(); }//remove zeroes in diagonal if possible
		int ind = m.getReversePivot(i);
		if (ind >= 0) {
			piv = m.get(i, ind);
			//normalize pivot row so pivot = 1
			if (piv != 1 && piv != 0 && abs(piv) > SYSTEM_EPSILON) {
				for (int l = 0; l < i; ++l) {
					m.set(i, l, m.get(i, l) / piv);
				}
				m.set(i, i, 1);
			}
			//proceed up the column to make each non-pivot value = 0
			for (j = i-1; j >= 0; --j) {
				if (m.get(j, i) != 0) {
					if (j != i) {
						piv = m.get(j, i);//m.get(i,i);
						for (k = 0; k<i; ++k) {
							m.set(j, k, m.get(j, k) - (m.get(i, k) * piv));
						}
						m.set(j, i, 0);
					}
				}
			}
		}
	}
	return m;
}

Matrix Matrix::upperTriangularize() {
	int i, j, k;
	double piv;
	Matrix m(rows, columns, element);
	//Gauss-Jordan elimination
	for (i = 0; i < rows; i++) {
		if (zeroInDiag() && canSwapOutZeroDiagonals()) { m = m.swapOutZeroDiagonals(); }//remove zeroes in diagonal if possible
		int ind = m.getPivot(i);
		if (ind >= 0) {
			piv = m.get(i, ind);
			//normalize pivot row so pivot = 1
			if (piv != 1 && piv != 0 && abs(piv) > SYSTEM_EPSILON) {
				for (int l = i+1; l < rows; ++l) {
					m.set(i, l, m.get(i, l) / piv);
				}
				m.set(i, i, 1);
			}
			//proceed down the column to make each non-pivot value = 0
			for (j = i+1; j < rows; j++) {
				if (m.get(j, i) != 0) {
					if (j != i) {
						piv = m.get(j, i);//m.get(i,i);
						for (k = i + 1; k < columns; k++) {
							m.set(j, k, m.get(j, k) - (m.get(i, k) * piv));
						}
						m.set(j, i, 0);
					}
				}
			}
		}
	}
	return m;
}

int Matrix::rank() {
	int rowval = 0;
	for (int i = 0; i < rows; ++i) {
			double elem = sumRow(i);
			if (abs(elem) > SYSTEM_EPSILON) {
				++rowval;
			}
	}
	return rowval;
}

int Matrix::nullity() {return columns-rank();}

void Matrix::removeZeroColumns(){
	for (int i = 0; i < columns; ++i) {
		if (sumColumn(i) == 0) { removeColumn(i); --i; }
	}
}

void Matrix::removeZeroRows(){
	for (int i = 0; i < rows; ++i) {
		if (sumRow(i) == 0) { 
			removeRow(i); --i; 
		}
	}
}

std::vector<int> Matrix::pivotColumns() {
	std::vector<int> piv;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			double val = get(i, j);
			if (val != 0) {
				piv.push_back(j);
				j = columns + 1;
			}
		}
	}
	return piv;
}

Matrix Matrix::nullSpace(){
	Matrix RRE = GaussianElimination();
	int kernelDimension = RRE.rank();//the rank of the matrix is equivalent to the dimension of the kernel
	//if (columns == kernelDimension) { return Matrix(1, 1); }

	//for a RRE matrix, the non-pivot columns are the free variables
	std::vector<int> piv = RRE.pivotColumns();
	RRE.removeZeroRows();

	std::vector<std::vector<double>> ans;
	int counter = 0;
	while (ans.size() < kernelDimension && counter<columns*10) {
		for (int i = 0; i < piv.size(); ++i) {
			Matrix Mcopy(RRE);//set up matrix

			//make one of the pivot columns = 1, the others = 0. rotate which one is 1 over each iteration to find all distinct solution vectors
			for (int p = 0; p < piv.size(); ++p) {
				if (p != i) { Mcopy.setColumnNonZeroValues(piv[p], 0); }
			}
			Mcopy.setColumnNonZeroValues(piv[i], 1);

			//run the solver with for Ax=b where b is the 0 vector
			ans.push_back(Mcopy.solve(std::vector<double>(columns)));
			++counter;
		}
	}
	Matrix ansMat(ans);
	return ansMat;
}

Matrix Matrix::upperHessenbergForm() {
	if (rows != columns) { return Matrix(); }
	Matrix M(rows, columns, element);

	//elimination with 3 cases:  
	for (int j = 0; j<columns - 1; ++j) {
		double pivot = M.get(j + 1, j);
		for (int i = j + 2; i < rows; ++i) {
			std::vector<double> temp = M.row(j + 1);
			if (M.get(i, j) != 0) {
				if (abs(M.get(i, j)) == abs(pivot)) {
					if (M.get(i, j) == pivot) { arrmultiply(temp, -1); goto skipUH; }
					if ((M.get(i, j) == (-1 * pivot)) || ((-1 * M.get(i, j)) == pivot)) { goto skipUH; }
				}
				if (abs(M.get(i, j)) != abs(pivot)) { arrmultiply(temp, (-M.get(i, j) / pivot)); }

			skipUH:
				M.addToRow(i, temp);
			}
		}
	}
	return Matrix(rows, columns, M.element);
}

Matrix Matrix::lowerHessenbergForm() {
	if (rows != columns) { return Matrix(); }
	Matrix A(rows, columns, element);

	for (int j = columns - 1; j >= 0; --j) {//go the other way to remove the upper triangle of non-diagonal values
		double pivot = A.get(j, j);
		for (int i = j - 1; i >= 0; --i) {
			if (i != j && A.get(i, j) != 0) {
				std::vector<double> temp = A.row(j);

				if (abs(A.get(i, j)) == abs(pivot)) {
					if (A.get(i, j) == pivot) {
						arrmultiply(temp, -1);
						goto skipLH;
					}
					if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) { goto skipLH; }
				}
				if (abs(A.get(i, j)) != abs(pivot)) {
					double a = -A.get(i, j) / pivot;
					arrmultiply(temp, a);
				}

			skipLH:
				A.addToRow(i, temp);
			}
		}
	}

	return Matrix(rows, columns, A.element);
}

Matrix Matrix::characteristicMatrix() {//numerical method
	Matrix next(rows, columns, element);
	for (int k = 1; k < rows; ++k) {
		Matrix M(columns, rows);
		M.identity();
		for (int p = 1; p<columns + 1; ++p) {
			if (p == (rows - k)) {
				M.element[((rows - k - 1)*rows) + (p - 1)] = (1 / (next.element[((rows - k)*columns) + (rows - k - 1)]));
			}
			else { M.element[((rows - k - 1)*rows) + (p - 1)] = -next.element[((rows - k)*columns) + (p - 1)] / (next.element[((columns - k)*columns) + (columns - k - 1)]); }
		}
		Matrix Minv(columns, rows);
		Minv.identity();
		for (int p = 1; p < columns + 1; ++p) { Minv.element[(columns - k - 1) * columns + (p - 1)] = next.element[(columns - k)*columns + (p - 1)]; }
		next = next.multiply(Minv, next.multiply(next, M));
	}

	for (int i = 1; i<next.rows; ++i) {//blank out anything that's not in the first row;
		for (int j = 0; j<next.columns; ++j) {
			if (i == j) { next.element[i*columns + j] = 1; }
			else { next.element[i*columns + j] = 0; }
		}
	}
	return Matrix(rows, columns, next.element);
}


bool Matrix::zeroInDiag() {//test if there are any 0's in the diagonal of a matrix
	for (int i = 0; i < rows; ++i) {
		if (get(i, i) == 0) { return true; }
	}
	return false;
}

double Matrix::det() {//determinant by exact (ie non-numerical) method, achieved by taking the product of diagonals in the reduced row echelon form (via Gaussian elimination)
	if (rows != columns) { return 0; }
	if (rows == 1) {
		return get(0, 0);
	}

	//MATRICES OF 2nd ORDER
	if (rows == 2) {
		return (get(0, 0)*get(1, 1) - get(0, 1)*get(1, 0));
	}

	//NTH ORDER -- use elimination
	Matrix M(rows, columns, element);

	int signFlip = 0;
	while (M.zeroInDiag() == true) {
		for (int i = 0; i < rows; ++i) {
			if (M.get(i, i) == 0) {
				for (int j = 0; j < columns; ++j) {
					if (M.get(i, j) != 0 && M.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create M non-zero diagonal
						M.swapRows(i, j); j = columns;
						++signFlip;	//each time you swap rows, you flip the sign of the determinant
					}
				}
			}
		}
	}

	//elimination with 3 cases:  
	for (int j = 0; j<columns; ++j) {
		double pivot = M.get(j, j);
		for (int i = j + 1; i < rows; ++i) {
			std::vector<double> temp = M.row(j);
			if (M.get(i, j) != 0) {
				if (abs(M.get(i, j)) == abs(pivot)) {
					if (M.get(i, j) == pivot) { arrmultiply(temp, -1); goto skip; }
					if ((M.get(i, j) == (-1 * pivot)) || ((-1 * M.get(i, j)) == pivot)) { goto skip; }
				}
				if (abs(M.get(i, j)) != abs(pivot)) { arrmultiply(temp, (-M.get(i, j) / pivot)); }

			skip:
				M.addToRow(i, temp);
			}
		}
	}

	double answer = 1;
	for (int i = 0; i < columns; ++i) {
		answer *= M.get(i, i);//multiply diagonals to get determinant
	}
	answer *= pow(-1, signFlip % 2);//fix flipped determinant
	return answer;
}//END DETERMINANT FUNCTION=======================================

double Matrix::detExact() {
	Matrix m = upperTriangularize();
	double det = 1;	
	for(int i=0; i<m.rows; ++i){
		det *= m.get(i, i);
	}
	return det;
}

Matrix Matrix::reducedRowEchelonForm() {
	Matrix A(rows, columns, element);

	while (A.zeroInDiag() == true) {
		for (int i = 0; i < rows; ++i) {
			if (A.get(i, i) == 0) {
				for (int j = 0; j < columns; ++j) {
					if (A.get(i, j) != 0 && A.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create a non-zero diagonal
						A.swapRows(i, j);
						j = columns;
					}
				}
			}
		}
	}

	for (int j = 0; j<rows; ++j) {//go down and eliminate the lower triangle of values
		double pivot = A.get(j, j);
		for (int i = 0; i < rows; ++i) {
			if (i != j && A.get(i, j) != 0) {
				std::vector<double> temp = A.row(j);

				if (abs(A.get(i, j)) == abs(pivot)) {
					if (A.get(i, j) == pivot) {
						arrmultiply(temp, -1);
						goto skip;
					}
					if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) { goto skip; }
				}
				if (abs(A.get(i, j)) != abs(pivot)) {
					double a = -A.get(i, j) / pivot;
					arrmultiply(temp, a);
				}

			skip:
				A.addToRow(i, temp);

			}
		}
	}


	for (int i = 0; i < rows; ++i) {//make diagonals into all 1's
		if (A.get(i, i) != 0 && A.get(i, i) != 1) {
			double a = ((double)1 / (A.get(i, i)));
			A.multiplyRow(i, a);
		}
	}
	return Matrix(rows, columns, A.element);
}

Matrix Matrix::addRow(std::vector<double> vec) {//adds vector to right hand side by default
	Matrix M(rows + 1, columns);
	int rw = rows;
	for (int i = 0; i < M.rows; ++i) {
		for (int j = 0; j < M.columns; ++j) {
			if (i == rw) {
				M.set(i, j, vec[j]);
			}
			else {
				if (i < rw) {
					M.set(i, j, get(i, j));
				}
				if (i > rw) {
					M.set(i, j, get(i - 1, j));
				}
			}
		}
	}
	return M;
}

Matrix Matrix::addRow(int rw, std::vector<double> vec) {
	Matrix M(rows + 1, columns);
	for (int i = 0; i < M.rows; ++i) {
		for (int j = 0; j < M.columns; ++j) {
			if (i == rw) {
				M.set(i, j, vec[j]);
			}
			else { 
				if (i < rw) {
					M.set(i, j, get(i,j));
				}
				if (i > rw) {
					M.set(i, j, get(i-1, j));
				}
			}
		}
	}
	return M;
}

Matrix Matrix::addColumn(std::vector<double> vec) {//adds to bottom of matrix by default
	Matrix M(rows, columns + 1);
	int cl = columns;
	for (int i = 0; i < M.rows; ++i) {
		for (int j = 0; j < M.columns; ++j) {
			if (j == cl) {
				M.set(i, j, vec[i]);
			}
			else {
				if (j < cl) {
					M.set(i, j, get(i, j));
				}
				if (j > cl) {
					M.set(i, j, get(i, j-1));
				}
			}
		}
	}
	return M;
}

Matrix Matrix::addColumn(int cl, std::vector<double> vec) {
	Matrix M(rows, columns + 1);
	for (int i = 0; i < M.rows; ++i) {
		for (int j = 0; j < M.columns; ++j) {
			if (j == cl) {
				M.set(i, j, vec[i]);
			}
			else {
				if (j < cl) {
					M.set(i, j, get(i, j));
				}
				if (j > cl) {
					M.set(i, j, get(i, j-1));
				}
			}
		}
	}
	return M;
}

bool Matrix::isInconsistent(std::vector<double> b) {//by Rouche-Capelli Theorem,  any system of linear equations (underdetermined or otherwise) is inconsistent if the rank of the augmented matrix is greater than the rank of the coefficient matrix
	if (rank() <= addColumn(columns-1, b).rank()) { return true; }
	return false;
}

bool Matrix::isUnderDetermined() {return (rows < columns);}
bool Matrix::isOverDetermined() { return (rows >= columns);}

bool Matrix::isSingular() {

	//FILL IN
	return false;
}

bool Matrix::solutionCheck(std::vector<double> x, std::vector<double> b) {//tests if Ax=b is a solution
	for (int i = 0; i < rows; ++i) {
		double val = 0;
		for (int j = 0; j < columns; ++j) {
			val += get(i, j)*x[j];
		}
		if (abs(val) - abs(b[i]) > SYSTEM_EPSILON) { return false; }
	}
	return true;
}

std::vector<double> Matrix::solve(std::vector<double> b){//solve for Ax=b, where x will be the returned vector of variable values
	if (isUnderDetermined() == true) {

	}
	if (isOverDetermined() == true) {

	}

	//if m==n
	int n = rows;
	std::vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = get(i, n - 1) / get(i, i);
		for (int k = i - 1; k >= 0; k--) {
			set(k, n - 1, get(k, n - 1) - get(k, i) * x[i]);
		}
	}
	return x;

}

std::vector<Matrix> Matrix::SVD() {
	/*Takes an mxn matrix a and decomposes it into udv, where u,v are
	left and right orthogonal transformation matrices, and d is a
	diagonal matrix of singular values.*/
	std::vector<Matrix> Mat;
	std::vector<double> w(columns);
	std::vector<double> eigens1, eigens2;
	Matrix v(rows,columns);
	Matrix M(rows,columns,element); 
	Matrix MT = M.transpose();

	//search for U
	Matrix MMT = M.multiply(M, MT);
	Mat.push_back(MMT.eigenvectors());

	//search for V
	Matrix MTM = M.multiply(MT, M);
	Mat.push_back(MTM.eigenvectors());

	//search for S
	eigens1 = MMT.eigenvaluesRealExact();
	eigens2 = MTM.eigenvaluesRealExact();
	eigens1.insert(eigens1.begin(), eigens2.begin(), eigens2.end());
	eigens2.clear();
	eigens1 = removeNullValues(eigens1);
	Mat.push_back(Vector(eigens1).toDiagonalMatrix());
	return(Mat);
}

std::vector<Matrix> Matrix::LU() {
	std::vector<Matrix> answer;
	Matrix A(rows, columns, element);
	Matrix L(rows, columns);
	Matrix U(rows, columns);

	int i, j, k;
	double sum = 0;

	for (i = 0; i < rows; i++) {
		U.set(i, i, 1);
	}

	for (j = 0; j < columns; j++) {
		for (i = j; i < rows; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L.get(i, k) * U.get(k, j);
			}
			L.set(i, j, A.get(i, j) - sum);
		}

		for (i = j; i < rows; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L.get(j, k) * U.get(k, i);
			}
			if (L.get(j, j) == 0) {
				//printf("det(L) may be close to 0!\n Can't divide by 0...\n");
				answer[0] = Matrix();
				answer[1] = Matrix();
				return answer;
			}
			U.set(j, i, ((A.get(j, i) - sum) / L.get(j, j)));
		}
	}
	answer.push_back(L);
	answer.push_back(U);
	return answer;
}

Matrix Matrix::L() {
	Matrix A(rows, columns, element);
	Matrix L(rows, columns);
	Matrix U(rows, columns);

	int i, j, k;
	double sum = 0;

	for (i = 0; i < rows; i++) {
		U.set(i, i, 1);
	}

	for (j = 0; j < columns; j++) {
		for (i = j; i < rows; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L.get(i, k) * U.get(k, j);
			}
			L.set(i, j, A.get(i, j) - sum);
		}

		for (i = j; i < rows; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L.get(j, k) * U.get(k, i);
			}
			if (L.get(j, j) == 0) {
				printf("det(L) close to 0!\n Can't divide by 0...\n");
				return Matrix();
			}
			U.set(j, i, ((A.get(j, i) - sum) / L.get(j, j)));
		}
	}
	return L;
}

Matrix Matrix::U() {
	Matrix A(rows, columns, element);
	Matrix L(rows, columns);
	Matrix U(rows, columns);

	int i, j, k;
	double sum = 0;

	for (i = 0; i < rows; i++) {
		U.set(i, i, 1);
	}

	for (j = 0; j < columns; j++) {
		for (i = j; i < rows; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L.get(i, k) * U.get(k, j);
			}
			L.set(i, j, A.get(i, j) - sum);
		}

		for (i = j; i < rows; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L.get(j, k) * U.get(k, i);
			}
			if (L.get(j, j) == 0) {
				printf("det(L) close to 0!\n Can't divide by 0...\n");
				return Matrix();
			}
			U.set(j, i, ((A.get(j, i) - sum) / L.get(j, j)));
		}
	}
	return U;
}

std::vector<Matrix> Matrix::QR() {//QR decomposition by Gram-Schmidt process
	std::vector<Matrix> A;
	Matrix Q(rows, columns);
	Matrix R(rows, columns);

	//set up Q matrix
	//===============
	Q.setColumn(0, normedArray(column(0)));
	for (int i = 0; i < columns; ++i) {
		std::vector<double> c = column(i);
		std::vector<double> a = c;
		for (int j = i; j > 0; --j) {
			std::vector<double> e = Q.column(j - 1);
			double dp = dotProduct(c, e);
			std::vector<double> tp = arrmultiply(e, dp);
			a = arrsubtract(a, tp);
		}
		a = normedArray(a);
		Q.setColumn(i, a);
	}

	//set up R matrix
	//===============
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (j >= i) {
				std::vector<double> a = column(j);
				std::vector<double> b = Q.column(i);
				double temp = dotProduct(a, b);
				R.set(i, j, temp);
			}
		}
	}
	A.push_back(Q);
	A.push_back(R);
	return A;
}

Matrix Matrix::Q() {//QR decomposition by Gram-Schmidt process
	Matrix Q(rows, columns);
	Matrix R(rows, columns);

	//set up Q matrix
	//===============
	Q.setColumn(0, normedArray(column(0)));
	for (int i = 1; i < rows - 1; ++i) {
		std::vector<double> c = column(i);
		for (int j = i; j > 0; --j) {
			std::vector<double> e = Q.column(j - 1);
			double dp = dotProduct(c, e);
			std::vector<double> tp = arrmultiply(e, dp);
			c = arrsubtract(c, tp);
		}
		c = normedArray(c);
		Q.setColumn(i, c);
	}

	//set up R matrix
	//===============
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (j >= i) {
				std::vector<double> a = column(j);
				std::vector<double> b = Q.column(i);
				double temp = dotProduct(a, b);
				R.set(i, j, temp);
			}
		}
	}
	return Matrix(rows, columns, Q.element);
}

Matrix Matrix::R() {//QR decomposition by Gram-Schmidt process
	Matrix Q(rows, columns);
	Matrix R(rows, columns);

	//set up Q matrix
	//===============
	Q.setColumn(0, normedArray(column(0)));
	for (int i = 1; i < rows - 1; ++i) {
		std::vector<double> c = column(i);
		for (int j = i; j > 0; --j) {
			std::vector<double> e = Q.column(j - 1);
			double dp = dotProduct(c, e);
			std::vector<double> tp = arrmultiply(e, dp);
			c = arrsubtract(c, tp);
		}
		c = normedArray(c);
		Q.setColumn(i, c);
	}

	//set up R matrix
	//===============
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (j >= i) {
				std::vector<double> a = column(j);
				std::vector<double> b = Q.column(i);
				double temp = dotProduct(a, b);
				R.set(i, j, temp);
			}
		}
	}
	return Matrix(rows, columns, R.element);
}

Matrix Matrix::inverseByQR() {
	if (det() == 0) { return Matrix(); }
	std::vector<Matrix> qr = QR();
	Matrix A = multiply(qr[1], qr[0].transpose());
	return Matrix(rows, columns, A.element);
}

Matrix Matrix::GramSchmidt() {//orthogonalization by Gram-Schmidt process using Gaussian Elimination

	if (det() == 0) { return Matrix(); }
	Matrix Aug(rows, columns, element);
	Matrix AT(rows, columns, element);
	AT = AT.transpose();
	Matrix A = multiply(*this, AT);


	while (A.zeroInDiag() == true) {
		for (int i = 0; i < rows; ++i) {
			if (A.get(i, i) == 0) {
				for (int j = 0; j < columns; ++j) {
					if (A.get(i, j) != 0 && A.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create a non-zero diagonal
						A.swapRows(i, j);
						Aug.swapRows(i, j);
						j = columns;
					}
				}
			}
		}
	}

	for (int j = 0; j<columns; ++j) {//go down and eliminate the lower triangle of values
		double pivot = A.get(j, j);
		for (int i = j + 1; i < rows; ++i) {
			if (i != j && A.get(i, j) != 0) {
				std::vector<double> temp = A.row(j);
				std::vector<double> temp2 = Aug.row(j);

				if (abs(A.get(i, j)) == abs(pivot)) {
					if (A.get(i, j) == pivot) {
						arrmultiply(temp, -1);
						arrmultiply(temp2, -1);
						goto skip;
					}
					if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) { goto skip; }
				}
				if (abs(A.get(i, j)) != abs(pivot)) {
					double a = -A.get(i, j) / pivot;
					arrmultiply(temp, a);
					arrmultiply(temp2, a);
				}

			skip:
				A.addToRow(i, temp);
				Aug.addToRow(i, temp2);

			}
		}
	}

	for (int i = 0; i < rows; ++i) {//make diagonals into all 1's
		if (A.get(i, i) != 0 && A.get(i, i) != 1) {
			double a = ((double)1 / (A.get(i, i)));
			A.multiplyRow(i, a);
			Aug.multiplyRow(i, a);
		}
	}
	return Matrix(rows, columns, Aug.element);
}

std::vector<int> Matrix::JacobiIndexing() {//returns indices i,j for largest element in matrix quickly for use in JacobiTransformation method
	std::vector<int> answer;
	answer.push_back(-1);
	answer.push_back(-1);
	double largest = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			if (i != j && abs(get(i, j)) > largest) {
				largest = abs(get(i, j));
				answer[0] = i;
				answer[1] = j;
			}
		}
	}
	return answer;
}

Matrix Matrix::GivensRotationMatrix(double a1, double a2, int r1, int r2) {//eliminate a2 with a1, a1 is in row r1, a2 is in row r2
	Matrix el(rows, rows);
	el.identity();

	double r = sqrt((a1*a1) + (a2*a2));
	double s = a1 / r;
	double negs = -s;
	double cs = a2 / r;

	el.element[r1*rows + r1] = cs;
	el.element[r2*rows + r2] = cs;
	el.element[r1*rows + r2] = s;
	el.element[r2*rows + r1] = negs;
	return Matrix(rows, rows, el.element);
}

Matrix Matrix::JacobiRotationMatrix(int p, int q, double c, double s) {//used in Jacobi Transformation -- will diagonalize a matrix and produce an eigenvector matrix
	Matrix el(rows, rows);
	el.identity();

	el.element[p*rows + p] = c;
	el.element[q*rows + q] = c;
	el.element[p*rows + q] = s;
	el.element[q*rows + p] = -s;
	return Matrix(rows, rows, el.element);
}

std::vector<Matrix> Matrix::JacobiTransformation() {//outputs an array with both diagonalized matrix and eigenvectors
													//if (rows != columns || isSymmetric()==false) { return NULL;	}
	Matrix D(rows, columns, element);  //matrix to rotate/diagonalize
	Matrix Eigenvectors(rows, columns); //matrix of Eigenvectors produced by multiplying each successive Jacobi Rotation matrix
	Eigenvectors.identity();

	int counter = 0;//counts the number of iterations
	long unsigned int limit = 10000;
	while (D.isDiagonalized() == false && counter<limit) {
		std::vector<int> elem = D.JacobiIndexing();
		if (elem[0] > D.rows || elem[1] > D.columns || elem[0] < 0 || elem[1] < 0 || elem[0] == elem[1]) { break; }
		int p = elem[0];
		int q = elem[1];
		double angle = (D.get(q, q) - D.get(p, p)) / (2 * D.get(p, q));
		int sign = 1;
		if (angle < 0) { sign = -1; }
		double t = sign / (abs(angle) + sqrt((angle*angle) + 1));
		double c = 1 / sqrt((t*t) + 1);
		double s = t*c;
		Matrix rot = D.JacobiRotationMatrix(p, q, c, s);
		Matrix rotT = rot.transpose();
		D = D.multiply(rotT, D);
		D *= rot;
		Eigenvectors = Eigenvectors.multiply(Eigenvectors, rot);
		++counter;
	}
	std::vector<Matrix> answer;
	answer.push_back(D);
	answer.push_back(Eigenvectors);
	return answer;
}

Matrix Matrix::eigenvectors() {
	if (isSymmetric() == true) {
		std::vector<Matrix> EV = JacobiTransformation();
		return EV[1];
	}
	else {
		std::vector<Matrix> qr = QR();
		Matrix Eigenvalues = qr[0];
		Matrix A2 = multiply(qr[1], qr[0]);

		int iterations = pow(rows*columns, 3);
		if (iterations > 500) { iterations = 500; }
		for (int i = 1; i < iterations; ++i) {//then, QR decomposition. Use Q as multiplicand where Eigenvectors = Q1*Q2*...Qn
			qr = A2.QR();
			if (qr[0].rows == 0 && qr[1].rows == 0) { break; }
			Eigenvalues = Eigenvalues.multiply(Eigenvalues, qr[0]);
			A2 = A2.multiply(qr[1], qr[0]);
			qr.clear();
		}
		return Eigenvalues;
	}
	return Matrix();
}

Matrix Matrix::Cholesky() {	//returns matrix L decomposition of matrix where A = LL^T
	if (rows <= 1 || columns <= 1) { return Matrix(); }
	Matrix L(rows, columns);
	L.set(0, 0, sqrt(get(0, 0)));

	for (int j = 0; j < columns; ++j) {
		for (int i = j; i < rows; ++i) {
			if (i == j) {
				double temp = get(j, j);
				double temp2 = 0;
				for (int k = 0; k < j; ++k) {
					temp2 += pow(L.get(j, k), 2);
				}
				L.set(j, j, sqrt(temp - temp2));
			}

			if (i > j) {
				double temp = (1 / L.get(j, j));
				double temp2 = get(i, j);
				double temp3 = 0;
				for (int k = 0; k < j; ++k) {
					temp3 += (L.get(i, k)*L.get(j, k));
				}
				if (temp != 0) { L.set(i, j, temp*(temp2 - temp3)); }
			}

		}
	}
	return Matrix(rows, columns, L.element);
}

std::vector<Matrix> Matrix::LDL() {	//returns matrix L and D decomposition of matrix where A = LDL^T
	std::vector<Matrix> answer;
	if (rows <= 1 || columns <= 1) { return answer; }
	Matrix L(rows, columns);
	Matrix D(rows, columns);
	L.set(0, 0, sqrt(get(0, 0)));

	for (int j = 0; j < columns; ++j) {
		for (int i = j; i < rows; ++i) {
			if (i == j) {
				L.set(j, j, 1);
				double temp = 0;
				for (int k = 0; k < j; ++k) {
					temp += pow(L.get(j, k), 2)*D.get(k, k);
				}
				D.set(j, j, get(j, j) - temp);
			}

			if (i > j) {
				double temp = (1 / D.get(j, j));
				double temp2 = get(i, j);
				double temp3 = 0;
				for (int k = 0; k < j; ++k) {
					temp3 += (L.get(i, k)*L.get(j, k)*D.get(k, k));
				}
				L.set(i, j, temp*(temp2 - temp3));
			}
		}
	}

	answer.push_back(L);
	answer.push_back(D);
	return answer;
}

Matrix Matrix::dominantEigenvector() {//calculated by power method
	Matrix A(rows,columns,element);
	Matrix rd(A.rows, 1);
	rd.set(0, 0, 1);//set first value of matrix to 1
	Matrix init = rd;
	for (int i = 0; i < 500; ++i) {
		rd = A.multiply(A, rd);
		rd = rd.multiply(rd, 1 / rd.element[A.rows - 1]);
	}
	return rd;//dominant eigenvector
}

double Matrix::dominantEigenvalue() {//calculated by power method
	Matrix A = *this;
	Matrix rd(A.rows, 1);
	rd.set(0, 0, 1);//set first value of matrix to 1
	Matrix init = rd;
	for (int i = 0; i < 500; ++i) {
		rd = A.multiply(A, rd);
		rd = rd.multiply(rd, 1 / rd.element[A.rows - 1]);
	}
	double egv = A.multiply(A.multiply(A, rd), rd).sumColumn(0);
	egv /= rd.multiply(rd, rd).sumColumn(0);
	return egv;//dominant eigenvalue
}

bool Matrix::isPositiveDefinite() {//if all eigenvalues are positive, the matrix is positive definite
	std::vector<double> vec = eigenvaluesRealExact();
	for (int i = 0; i < vec.size(); ++i) {
		if (vec[i] <= 0) { return false; }
	}
	return true;
}

bool Matrix::isPositiveSemidefinite() {//if all eigenvalues are positive or 0, the matrix is positive semi-definite
	std::vector<double> vec = eigenvaluesRealExact();
	for (int i = 0; i < vec.size(); ++i) {
		if (vec[i] < 0) { return false; }
	}
	return true;
}

std::vector<double> Matrix::eigenvaluesByGaussianElimination() {
	addColumn(std::vector<double>(rows));//make augmented matrix w/ column vector of 0's for solution
	Matrix GE = GaussianElimination();
	return GE.characteristicPolynomial().realRoots();
}

std::vector<double> Matrix::eigenvaluesRealExact() {//finds real eigenvalues by method of QR algorithm with Hessenberg intermediate step
	Polynomial p = characteristicPolynomial();
	return p.realRoots();
}

std::vector<ComplexNumber> Matrix::eigenvaluesRealAndComplexExact() {
	//Finds real eigenvalues by method of QR algorithm with Hessenberg intermediate step,
	//then either factors those roots to get the last few complex roots or uses the 
	//numerical Durand-Kerner method of root-finding
	std::wcout << L"processing...this may take awhile for large matrices...\n";
	Polynomial p = characteristicPolynomial();
	std::vector<double> realRoots = findRealRoots(p);
	std::vector<ComplexNumber> complexRoots;
	for (int i = 0; i < realRoots.size(); ++i) {
		//std::wcout << realRoots[i] << L"\n";//debug
		complexRoots.push_back(ComplexNumber(realRoots[i], 0));
	}

	std::vector<double> bounds = p.getComplexRootBounds();
	double lowbound = bounds[0];
	double highbound = bounds[1];
	double lowboundC = bounds[0];

	//std::wcout << L"root bounds = [" << lowbound << L"," << highbound << L"]\n";//debug
	for (int i = 0; i < realRoots.size(); ++i) {
		if (std::abs(p.terms - (int)realRoots.size()) <= 3) {
			//Use the shift method only if the real roots reduce it down to an easily calculable size (order 2 or less),
			//otherwise the results will be inaccurate
			if (p.terms <= 3) { i = realRoots.size(); break; }
			if (p.terms > 3) { p = p.factorOutBinomial(realRoots[i]); }
		}
	}
	p.display();//debug
	if (p.terms < 3) { return complexRoots; }

	if (p.terms == 3) {
		double part2 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
		if (part2 < 0) {
			part2 = part2 / (2 * p.coefficient[2]);
			double part1 = (-p.coefficient[1] / (2 * p.coefficient[2]));
			complexRoots.push_back(ComplexNumber(part1, part2));
			complexRoots.push_back(ComplexNumber(part1, -part2));
		}
		return complexRoots;
	}

	if (p.terms > 3) {


		//start finding roots with Durand-Kerner method
		int iterations = 10000;
		ComplexNumber z = ComplexNumber(lowbound, lowboundC);
		int size = sizeof(z);
		std::vector<ComplexNumber> R;
		for (int i = 0; i < p.terms; i++) { R.push_back(z.exponent(z, i)); }

		for (int i = 0; i < iterations; i++) {
			for (int j = 0; j < p.terms; j++) {
				ComplexNumber B = p.HornerEvaluateComplex(R[j]);
				for (int k = 0; k < p.terms; k++) {
					if (k != j) { B /= (R[j] - R[k]); }
				}
				R[j] -= B;
			}
		}


		for (int i = 0; i < R.size(); ++i) { //now filter out all the bad results
			if (R[i].Im != 0) { //only complex roots accepted
				R[i] = NewtonsMethodComplex(p, R[i]);
				ComplexNumber temp = p.HornerEvaluateComplex(R[i]);
				double a = std::abs(temp.Re);
				double b = std::abs(temp.Im);
				if (a < 0.1 && b < 0.1) {
					a = std::abs(R[i].Re);
					b = std::abs(R[i].Im);
					bool isSimilar = false;
					for (int i = 0; i < complexRoots.size(); ++i) {
						double x = std::abs(complexRoots[i].Re);
						double y = std::abs(complexRoots[i].Im);
						if (std::abs(a - x) < 0.001 && std::abs(b - y) < 0.001) { isSimilar = true; }
					}
					if (isSimilar == false) { //if this is indeed a new root, save the root and its conjugate
						complexRoots.push_back(R[i]);
						//R[i].display();
						complexRoots.push_back(ComplexNumber(R[i].Re, R[i].Im*-1));
					}
				}
			}
		}
		R.clear();
	}
	return complexRoots;
}

std::vector<double> Matrix::eigenvaluesNumerical() {//calculated by power method w/shifts
	Matrix A = *this;
	Matrix v;
	std::vector<double> egv;
	for (int i = 0; i < A.rows; ++i) {	//A' = A - (eigenvalue)*U*UT
		v = A.dominantEigenvector();
		double eigenval = A.multiply(A.multiply(A, v), v).sumColumn(0);
		double val = v.multiply(v, v).sumColumn(0);
		if (val != 0) { eigenval /= val; }
		else { eigenval = 0; }
		egv.push_back(eigenval);
		v.normalize();

		Matrix U = v.multiply(v, v.transpose());

		//U.display();
		U = U.multiply(U, -1 * eigenval);
		A = A.add(A, U);
	}
	return egv;
}

Matrix Matrix::eigenvalueMatrix() {
	Matrix A = *this;
	std::vector<Matrix> QR;

	int iterations = pow((A.rows*A.columns), 2);
	for (int i = 0; i < iterations; ++i) {//first, LU decomposition
		QR = A.LU();
		if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
		A = A.multiply(QR[1], QR[0]);
	}

	for (int i = 0; i < iterations; ++i) {//then, QR decomposition
		QR = A.QR();
		if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
		A = A.multiply(QR[1], QR[0]);
	}
	for (int i = 0; i < iterations; ++i) {//then, QR decomposition
		QR = A.HouseholderQR();
		if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
		A = A.multiply(QR[1], QR[0]);
	}

	return A;
}

Matrix Matrix::eigenvalueMatrix(Matrix A, int loops) {//recursive version of eigenvalueMatrix calculation
	if (loops <= 1) { return A; }
	--loops;
	std::vector<Matrix> QR;

	int iterations = (A.rows*A.columns);
	for (int i = 0; i < iterations; ++i) {//first, LU decomposition
		QR = A.LU();
		if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
		A = A.multiply(QR[1], QR[0]);
	}

	for (int i = 0; i < iterations; ++i) {//then, QR decomposition
		QR = A.QR();
		if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
		A = A.multiply(QR[1], QR[0]);
	}
	for (int i = 0; i < iterations; ++i) {//then, QR decomposition
		QR = A.HouseholderQR();
		if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
		A = A.multiply(QR[1], QR[0]);
	}

	return eigenvalueMatrix(A, loops);
}

void Matrix::power(double n) {//exponentiates all elements in the matrix by some power
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			set(i, j, pow(get(i, j), n));
		}
	}
}

Matrix Matrix::leastSquares(Matrix X, Matrix Y) {//B=(X^T*X)^-1 * X^T * Y  
	Matrix XT = X.transpose();
	return Y.multiply(X.multiply(X.multiply(XT, X).inverse(), XT), Y);
}

Matrix Matrix::leastSquaresMatrix() {	//B=(X^T*X)^-1 * X^T * Y 
										/*note:  the matrix must be formatted in MINITAB format, such that the first column is the Y values, the rest are the
										column vectors of the X_i variable values. */

	Matrix M(rows, columns, element);//copy parent matrix
	if (columns != 2 && rows == 2) { M = M.transpose(); }//put matrix into column vector form

	Matrix X(M.rows, M.columns, element);	//isolate X variable values
	for (int i = 0; i < X.rows; ++i) {
		X.set(i, 0, 1);
	}

	Matrix Y(M.rows, 1, M.column(0));	//create seperate vector of y values

	return multiply(multiply(X.transpose(), X).inverse(), multiply(X.transpose(), Y));
}

Function Matrix::leastSquares() {	//B=(X^T*X)^-1 * X^T * Y 
									/*note:  the matrix must be formatted in MINITAB format, such that the first column is the Y values, the rest are the
									column vectors of the X_i variable values. See https://onlinecourses.science.psu.edu/stat501/node/292 for
									examples to prove this technique. There is test code also available at:  https://www.wessa.net/rwasp_multipleregression.wasp */

	Matrix M(rows, columns, element);//copy parent matrix
	if (columns != 2 && rows == 2) { M = M.transpose(); }//put matrix into column vector form

	Matrix X(M.rows, M.columns, element);	//isolate X variable values
	for (int i = 0; i < X.rows; ++i) {
		X.set(i, 0, 1);
	}

	Matrix Y(M.rows, 1, M.column(0));	//create seperate vector of y values

	Matrix coef = multiply(multiply(X.transpose(), X).inverse(), multiply(X.transpose(), Y));
	std::vector<double> cf;
	std::wstring funct = L"f(";
	for (int i = 0; i < coef.rows; ++i) {
		cf.push_back(coef.get(i, 0));
		funct += getLetterByAlphabeticalOrderString(i) + L",";
	}
	funct.pop_back();
	funct += L") = ";

	for (int i = 0; i < cf.size(); ++i) {
		funct += to_stringPrecision(cf[i]);
		if (i > 0) { funct += L"*" + getLetterByAlphabeticalOrderString(i - 1); }
		funct += L" + ";
	}
	//Remove final "+ " chars from string
	funct.pop_back();
	funct.pop_back();

	return Function(funct);
}

std::vector<double> Matrix::residuals() {
	//get least squares regression function
	Function f = leastSquares();

	//get y-hat values
	std::vector<double> yhat;
	for (int i = 0; i < rows; ++i) {
		std::vector<double> rw = row(i);
		std::vector<double> vars;
		for (int j = 1; j < columns; ++j) {
			vars.push_back(rw[j]);
		}
		yhat.push_back(f.evaluate(vars));
	}

	//calculate resid. = y-yhat
	std::vector<double> residuals;
	for (int i = 0; i < yhat.size(); ++i) {
		residuals.push_back(get(i, 0) - yhat[i]);//in MINITAB format, column 0 is the y-value column
	}

	return residuals;
}

double Matrix::Rsquared() {// aka Pearsson's "R"^2
						   //for explaination, see https://people.richland.edu/james/ictcm/2004/multiple.html

	Function f = leastSquares();//get leastSquares multivariate function

								//get y-hat values
	std::vector<double> yhat;
	for (int i = 0; i < rows; ++i) {
		std::vector<double> rw = row(i);
		std::vector<double> vars;
		for (int j = 1; j < columns; ++j) {
			vars.push_back(rw[j]);
		}
		yhat.push_back(f.evaluate(vars));
	}

	//get residuals and residual standard deviation
	std::vector<double> resi = residuals();
	double residualSD = arraySampleStandardDeviation(resi);

	double ymean = columnMean(0);	//ymean

									//fill arrays of error calculations
	std::vector<double> explainedError;
	std::vector<double> unexplainedError;
	std::vector<double> totalError;
	for (int i = 0; i < resi.size(); ++i) {
		explainedError.push_back(pow(yhat[i] - ymean, 2));		// sum(yhat - yman)^2
		unexplainedError.push_back(pow(get(i, 0) - yhat[i], 2));	// sum(y - yhat)^2
		totalError.push_back(pow(get(i, 0) - ymean, 2));			// sum(y - ymean)^2
	}

	//degress of freedom calculations
	//===============================
	double dfRegression = columns - 1;	//#predictor varibles aka 'k'
	double dfTotal = rows - 1;//total degrees freedom
	double dfResiduals = dfTotal - dfRegression;//# of residuals

												//sum of squares calculations
												//===========================
	double SSregression = arrsum(explainedError);
	double SSresiduals = arrsum(unexplainedError);
	double SStotal = arrsum(totalError);

	//mean square (MS) calculations
	//=============================
	double MSregression = SSregression / dfRegression;
	double MSresiduals = SSresiduals / dfResiduals;
	double MStotal = SStotal / dfTotal;

	//important statistics
	//====================
	double sqrtErrorVariance = sqrt(MSresiduals);	//equal to standard deviation of distribution of estimate
	double FTestStatistic = MSregression / MSresiduals;
	double Rsquared = (SStotal - SSresiduals) / SStotal;
	double adjustedRsquared = (MStotal - MSresiduals) / MStotal;

	return Rsquared;
}

double Matrix::adjustedRsquared() {// aka Pearsson's "R"^2
								   //for explaination, see https://people.richland.edu/james/ictcm/2004/multiple.html

	Function f = leastSquares();//get leastSquares multivariate function

								//get y-hat values
	std::vector<double> yhat;
	for (int i = 0; i < rows; ++i) {
		std::vector<double> rw = row(i);
		std::vector<double> vars;
		for (int j = 1; j < columns; ++j) {
			vars.push_back(rw[j]);
		}
		yhat.push_back(f.evaluate(vars));
	}

	//get residuals and residual standard deviation
	std::vector<double> resi = residuals();
	double residualSD = arraySampleStandardDeviation(resi);

	double ymean = columnMean(0);	//ymean

									//fill arrays of error calculations
	std::vector<double> explainedError;
	std::vector<double> unexplainedError;
	std::vector<double> totalError;
	for (int i = 0; i < resi.size(); ++i) {
		explainedError.push_back(pow(yhat[i] - ymean, 2));			// sum(yhat - yman)^2
		unexplainedError.push_back(pow(get(i, 0) - yhat[i], 2));	// sum(y - yhat)^2
		totalError.push_back(pow(get(i, 0) - ymean, 2));			// sum(y - ymean)^2
	}

	//degress of freedom calculations
	//===============================
	double dfRegression = columns - 1;	//#predictor varibles aka 'k'
	double dfTotal = rows - 1;//total degrees freedom
	double dfResiduals = dfTotal - dfRegression;//# of residuals

												//sum of squares calculations
												//===========================
	double SSregression = arrsum(explainedError);
	double SSresiduals = arrsum(unexplainedError);
	double SStotal = arrsum(totalError);

	//mean square (MS) calculations
	//=============================
	double MSregression = SSregression / dfRegression;
	double MSresiduals = SSresiduals / dfResiduals;
	double MStotal = SStotal / dfTotal;

	//important statistics
	//====================
	double sqrtErrorVariance = sqrt(MSresiduals);	//equal to standard deviation of distribution of estimate
	double FTestStatistic = MSregression / MSresiduals;
	double Rsquared = (SStotal - SSresiduals) / SStotal;
	double adjustedRsquared = (MStotal - MSresiduals) / MStotal;

	return adjustedRsquared;
}

double Matrix::FTestStatistic() {
	//for explaination, see https://people.richland.edu/james/ictcm/2004/multiple.html

	Function f = leastSquares();//get leastSquares multivariate function

								//get y-hat values
	std::vector<double> yhat;
	for (int i = 0; i < rows; ++i) {
		std::vector<double> rw = row(i);
		std::vector<double> vars;
		for (int j = 1; j < columns; ++j) {
			vars.push_back(rw[j]);
		}
		yhat.push_back(f.evaluate(vars));
	}

	//get residuals and residual standard deviation
	std::vector<double> resi = residuals();
	double residualSD = arraySampleStandardDeviation(resi);

	double ymean = columnMean(0);	//ymean

									//fill arrays of error calculations
	std::vector<double> explainedError;
	std::vector<double> unexplainedError;
	std::vector<double> totalError;
	for (int i = 0; i < resi.size(); ++i) {
		explainedError.push_back(pow(yhat[i] - ymean, 2));			// sum(yhat - yman)^2
		unexplainedError.push_back(pow(get(i, 0) - yhat[i], 2));	// sum(y - yhat)^2
		totalError.push_back(pow(get(i, 0) - ymean, 2));			// sum(y - ymean)^2
	}

	//degress of freedom calculations
	//===============================
	double dfRegression = columns - 1;	//#predictor varibles aka 'k'
	double dfTotal = rows - 1;//total degrees freedom
	double dfResiduals = dfTotal - dfRegression;//# of residuals

	//sum of squares calculations
	//===========================
	double SSregression = arrsum(explainedError);
	double SSresiduals = arrsum(unexplainedError);
	double SStotal = arrsum(totalError);

	//mean square (MS) calculations
	//=============================
	double MSregression = SSregression / dfRegression;
	double MSresiduals = SSresiduals / dfResiduals;
	double MStotal = SStotal / dfTotal;

	//important statistics
	//====================
	double sqrtErrorVariance = sqrt(MSresiduals);	//equal to standard deviation of distribution of estimate
	double FTestStatistic = MSregression / MSresiduals;
	double Rsquared = (SStotal - SSresiduals) / SStotal;
	double adjustedRsquared = (MStotal - MSresiduals) / MStotal;

	return FTestStatistic;
}

std::vector<Matrix> Matrix::thinQR() {
	std::vector<Matrix> answer = QR();
	answer[0];
	answer[1].trim();
	return answer;
}

std::vector<Matrix> Matrix::HouseholderQR() {//given some matrix, this produces its Householder QR factorization
	Matrix A = *this;
	Matrix Q(rows, columns);
	Q.identity();
	for (int j = 0; j < columns - 1; ++j) {
		std::vector<double> temp = A.column(j);
		std::vector<double> vec;
		for (int i = j; i < rows; ++i) { vec[i - j] = temp[i]; }
		double temp1 = sgn(vec[0])*arrNorm(vec);
		vec[0] += temp1;
		vec = normedArray(vec);//normalize here
		Matrix Vec(rows - j, 1, vec);
		Matrix H = Vec.GramMatrix(Vec);
		H *= 2;
		Matrix I(rows - j, columns - j);
		I.identity();
		H = I - H;
		if (j > 0) { H = H.expandToUpperLeft(rows, columns); }
		A = H*A;
		Q *= H;
	}
	std::vector<Matrix> QR;
	QR.push_back(Q);
	QR.push_back(A);
	return QR;
}

std::vector<Matrix> Matrix::thinHouseholderQR() {
	std::vector<Matrix> answer = HouseholderQR();
	answer[0];
	answer[1].trim();
	return answer;
}

Matrix Matrix::pseudoinverse() {//calculates the Moore-Penrose inverse: A* = R^-1 * Q^T.  Note:  A*A = I.
	std::vector<Matrix> qr = thinQR();
	return multiply(qr[1].inverse(), qr[0].transpose());
}

Polynomial Matrix::characteristicPolynomial() {
	std::vector<double> ans = characteristicMatrix().row(0);
	std::vector<double> answer = initZeros(answer, columns + 1);
	for (int i = 0; i < columns; ++i) { answer[columns - 1 - i] = -ans[i]; }
	answer[columns] = 1;
	return Polynomial(answer, columns + 1);
}

ComplexMatrix Matrix::toComplexMatrix() {
	std::vector<ComplexNumber> c;
	for (int i = 0; i < element.size(); ++i) {
		c.push_back(element[i]);
	}
	return ComplexMatrix(rows,columns,c);
}

std::vector<std::vector<double>> Matrix::toVectorArray() {
	std::vector<std::vector<double>> arr;
	for (int i = 0; i < rows; ++i) {
		std::vector<double> ans;
		for (int j = 0; j < columns; ++j) {
			ans.push_back(get(i, j));
		}
		arr.push_back(ans);
	}
	return arr;
}

Quaternion Matrix::toQuaternion() {//for use in homogeneous coordinate system for 3D graphics
	if (columns == rows == 4) {//if using a 4x4 matrix representation of a quaternion
		std::vector<double> d = diagonal();
		return Quaternion(d);
	}

	if (columns == 1 && rows == 4) {//if using a single vector
		std::vector<double> c = column(0);
		return Quaternion(c);
	}
	return Quaternion();//else return NULL
}//=============================================


SparseMatrix::SparseMatrix() {
	rows = columns = size = savedItems = matrixPrecision = 0;
}
SparseMatrix::~SparseMatrix() {};	//destructor

SparseMatrix::SparseMatrix(int row, int col) {//constructor -- matrix of 0's
	rows = row; columns = col; size = row*col; matrixPrecision = 4; largestElement = savedItems = 0;
}

SparseMatrix::SparseMatrix(int row, int col, int items, double* element2) {//constructor
	savedItems = items;
	largestElement = modulus = 0;
	for (int i = 0; i<items; i += 3) {
		element.push_back(element2[i]);
		element.push_back(element2[i + 1]);
		element.push_back(element2[i + 2]);
		if (getLength(element[i]) > largestElement) { largestElement = getLength(element[i]); }
	}
	rows = row; columns = col; size = row*col; matrixPrecision = 4;
}

SparseMatrix::SparseMatrix(int row, int col, std::vector<double> element2) {//constructor
	savedItems = element2.size() / 3;
	largestElement = modulus = 0;
	for (int i = 0; i<element2.size() / 3; i += 3) {
		element.push_back(element2[i]);
		element.push_back(element2[i + 1]);
		element.push_back(element2[i + 2]);
		if (getLength(element[i]) > largestElement) { largestElement = getLength(element[i]); }
	}
	rows = row; columns = col; size = row*col; matrixPrecision = 4;
}

		//Sparse Matrix methods
		//=====================
		double SparseMatrix::get(int i, int j) {
			for (int i = 0; i < savedItems; i += 3) {
				if (element[i + 1] == i && element[i + 2] == j) { return element[i]; }
			}
			return 0;
		}

		void SparseMatrix::set(int i, int j, double x) {
			for (int i = 0; i < savedItems; i += 3) {
				if (element[i + 1] == i && element[i + 2] == j) { element[i] = x; return; }
			}
			element.push_back(x); 
			element.push_back(i); 
			element.push_back(j);
			savedItems = element.size() / 3;
		}

		void SparseMatrix::setModulus(int n) {
			modulus = n;
			for (int i = 0; i < savedItems; ++i) {
				int b = element[i];
				b %= n;
				element[i] = b;
			}
		}

		double SparseMatrix::sparsity() {
			return (savedItems / size);
		}

		bool SparseMatrix::isIntegerMatrix() {
			for (int i = 0; i < savedItems; i += 3) {
				if (element[i] != floor(element[i])) { return false; }
			}
			return true;
		}

		Matrix SparseMatrix::toMatrix() {
			std::wstring str = toString();
			return ParseMatrix(str);
		}

		std::wstring SparseMatrix::toString() {
			std::wostringstream s;
			int defaultLength = 4;
			bool tf = isIntegerMatrix();
			if (tf == true) { --defaultLength; }
			if (tf == false) { ++defaultLength; }
			s << L"(" << rows << L" x " << columns << L")\n";
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < columns; ++j) {
					s <<
						std::fixed <<
						std::left <<
						std::setw(defaultLength + (matrixPrecision)+(largestElement)) <<
						std::setprecision(matrixPrecision) <<
						get(i, j);
				}
				s << L"\n";
			}
			std::wstring temp = s.str();
			s.clear();
			return temp;
		}

		void SparseMatrix::display() {std::wcout << toString() << std::endl;}
		//======================================================================
		ComplexMatrix::ComplexMatrix() {
			rows = columns = size = complexMatrixPrecision = 0;
		}
		ComplexMatrix::~ComplexMatrix() { element.clear(); }	//destructor

		ComplexMatrix::ComplexMatrix(int n) {
			for (int i = 0; i < n*n; ++i) {
				element.push_back(ComplexNumber(0));
			}
			rows = n;
			columns = n;
			size = n * n;
			complexMatrixPrecision = 4;
			largestElementR = largestElementC = 0;
		}

		ComplexMatrix::ComplexMatrix(int row, int col) {//constructor -- matrix of 0's
			for (int i = 0; i < row*col; ++i) {
				element.push_back(ComplexNumber(0));
			}
			rows = row; columns = col; size = row*col; complexMatrixPrecision = 4; largestElementR = largestElementC = 0;
		}

		ComplexMatrix::ComplexMatrix(int row, int col, double* element2, double* celement2) {//constructor
			largestElementR = largestElementC = 0;
			for (int i = 0; i < row; ++i) {
				for (int j = 0; j < col; ++j) {
					double n = element2[i*col + j];
					double n2 = celement2[i*col + j];
					element.push_back(ComplexNumber(element2[i*row + j], celement2[i*row + j]));
					if (getLength(n) > largestElementR) { largestElementR = getLength(n); }
					if (getLength(n2) > largestElementC) { largestElementC = getLength(n2); }
				}
			}
			rows = row; columns = col; size = row*col; complexMatrixPrecision = 4;
		}

		ComplexMatrix::ComplexMatrix(int row, int col, std::vector<double> element2, std::vector<double> celement2) {//constructor
			largestElementR = largestElementC = 0;
			for (int i = 0; i < row; ++i) {
				for (int j = 0; j < col; ++j) {
					double n = element2[i*col + j];
					double n2 = celement2[i*col + j];
					element.push_back(ComplexNumber(element2[i*row + j], celement2[i*row + j]));
					if (getLength(n) > largestElementR) { largestElementR = getLength(n); }
					if (getLength(n2) > largestElementC) { largestElementC = getLength(n2); }
				}
			}
			rows = row; columns = col; size = row*col; complexMatrixPrecision = 4;
		}

		ComplexMatrix::ComplexMatrix(int row, int col, std::vector<ComplexNumber> elements2) {//constructor
			element = elements2;
			rows = row;
			columns = col;
			size = row * col;
			complexMatrixPrecision = 4;
			largestElementR = largestElementC = 0;
			for (int i = 0; i < row; ++i) {
				for (int j = 0; j < col; ++j) {
					double n = elements2[i*col + j].Re;
					double n2 = elements2[i*col + j].Im;
					if (getLength(n) > largestElementR) { largestElementR = getLength(n); }
					if (getLength(n2) > largestElementC) { largestElementC = getLength(n2); }
				}
			}
		}

		ComplexMatrix::ComplexMatrix(int row, int col, std::vector<std::complex<double>> elements2) {//constructor
			rows = row;
			columns = col;
			size = row * col;
			complexMatrixPrecision = 4;
			largestElementR = largestElementC = 0;
			element = std::vector<ComplexNumber>(row*col, ComplexNumber(0));
			std::fill(element.begin(), element.end(), ComplexNumber(0));
			for (int i = 0; i < row; ++i) {
				for (int j = 0; j < col; ++j) {
					element[i*col + j] = ComplexNumber(elements2[i]);
					double n = real(elements2[i*col + j]);
					double n2 = imag(elements2[i*col + j]);
					if (getLength(n) > largestElementR) { largestElementR = getLength(n); }
					if (getLength(n2) > largestElementC) { largestElementC = getLength(n2); }
				}
			}
		}

				void ComplexMatrix::randomize() {//creates a random matrix of double values
					srand(time(0) + clock());
					largestElementR = 0;
					largestElementC = 0;
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							element[i*columns + j].Re = ((double)rand() / (30000)) *  pow(-1, rand() % 2);
							if (getLength(element[i*columns + j].Re) > largestElementR) { largestElementR = getLength(element[i*columns + j].Re); }

							element[i*columns + j].Im = ((double)rand() / (30000)) *  pow(-1, rand() % 2);
							if (getLength(element[i*columns + j].Im) > largestElementC) { largestElementC = getLength(element[i*columns + j].Im); }
						}
					}
					complexMatrixPrecision = 4;
				}

				void ComplexMatrix::randomizeInteger() {//creates a random matrix of integer values
					srand(time(0) + clock());
					largestElementR = 0;
					largestElementC = 0;
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							double a = (int)((double)rand() / (pow(10, (rand() % 4) + 2))) * pow(-1, rand() % 2);
							element[i*columns + j].Re = a;
							if (getLength(a) > largestElementR) { largestElementR = getLength(a); }

							a = (int)((double)rand() / (pow(10, (rand() % 4) + 2))) * pow(-1, rand() % 2);
							element[i*columns + j].Im = a;
							if (getLength(a) > largestElementC) { largestElementC = getLength(a); }
						}
					}
					complexMatrixPrecision = 0;
				}

				void ComplexMatrix::randomizeBoolean() {//creates a random matrix of double values
					srand(time(0) + clock());
					largestElementR = 1;
					largestElementC = 1;
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							element[i*columns + j].Re = rand() % 2;
							element[i*columns + j].Im = rand() % 2;
						}
					}
					complexMatrixPrecision = 0;
				}

				//class methods
				//==============
				double ComplexMatrix::getReal(int i, int j) { return element[(i*columns) + j].Re; }
				void ComplexMatrix::setReal(int i, int j, double x) { element[(i*columns) + j].Re = x; }
				double ComplexMatrix::getImaginary(int i, int j) { return element[(i*columns) + j].Im; }
				void ComplexMatrix::setImaginary(int i, int j, double x) { element[(i*columns) + j].Im = x; }
				ComplexNumber ComplexMatrix::get(int i, int j) { return ComplexNumber(element[(i*columns) + j]); }
				void ComplexMatrix::set(int i, int j, ComplexNumber x) { element[(i*columns) + j] = x; }

				void ComplexMatrix::identity() {//turns the matrix into an identity matrix
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							if (i == j) { element[i*columns + j] = ComplexNumber(1); }
							else { element[i*columns + j] = ComplexNumber(0); }
						}
					}
					largestElementR = largestElementC = 0;
					complexMatrixPrecision = 0;
				}

				ComplexMatrix ComplexMatrix::submatrix(int i, int j, int i2, int j2) {
					double* vals = new double[(i2*j2) - (i*j)];
					double* valsC = new double[(i2*j2) - (i*j)];
					int counter = 0;
					for (int a = i; a < i2; ++a) {
						for (int b = j; b < j2; ++b) {
							vals[counter] = getReal(a, b);
							valsC[counter] = getImaginary(a, b);
							counter++;
						}
					}
					return ComplexMatrix(i2 - i, j2 - j, vals, valsC);
				}

				ComplexMatrix ComplexMatrix::expandToUpperLeft(int r, int c) {//takes original matrix and creates a new matrix with the original in the lower right corner.
					ComplexMatrix newEl(r, c);
					if (r > rows && c > columns) {
						newEl.identity();
						for (int i = r - rows; i < r; ++i) {
							for (int j = c - columns; j < c; ++j) {
								newEl.setReal(i, j, getReal((i - (r - rows)), (j - (c - columns))));
								newEl.setImaginary(i, j, getImaginary((i - (r - rows)), (j - (c - columns))));
							}
						}
					}
					return ComplexMatrix(r, c, newEl.element);
				}

				ComplexMatrix ComplexMatrix::expandToLowerRight(int r, int c) {//takes original matrix and creates a new matrix with the original in the upper left corner.
					ComplexMatrix newEl(r, c);
					if (r > rows && c > columns) {
						newEl.identity();
						for (int i = 0; i < rows; ++i) {
							for (int j = 0; j < columns; ++j) {
								newEl.setReal(i, j, getReal(i, j));
								newEl.setImaginary(i, j, getImaginary(i, j));
							}
						}
					}
					return ComplexMatrix(r, c, newEl.element);
				}

				ComplexMatrix ComplexMatrix::extendRows(int r) {
					ComplexMatrix newEl(rows+r, columns);
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							newEl.setReal(i, j, getReal(i, j));
							newEl.setImaginary(i, j, getImaginary(i, j));
						}
					}
					return ComplexMatrix(rows+r, columns, newEl.element);
				}

				ComplexMatrix ComplexMatrix::extendColumns(int c) {
					ComplexMatrix newEl(rows, columns+c);
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							newEl.setReal(i, j, getReal(i, j));
							newEl.setImaginary(i, j, getImaginary(i, j));
						}
					}
					return ComplexMatrix(rows, columns+c, newEl.element);
				}

				ComplexMatrix ComplexMatrix::concatenate(ComplexMatrix A, ComplexMatrix B) {//creates a larger matrix from two matrices, where each matrix is next to the other one
					int r = A.rows;
					if (B.rows > A.rows) { r = B.rows; }
					int c = A.columns + B.columns;
					ComplexMatrix C(r, c);

					for (int i = 0; i < r; ++i) {
						for (int j = 0; j < c; ++j) {
							if (j < A.columns && i < A.rows) {
								C.setReal(i, j, A.getReal(i, j));
								C.setImaginary(i, j, A.getImaginary(i, j));
							}
							if (j >= A.columns && i < B.rows) {
								C.setReal(i, j, B.getReal(i, j - A.columns));
								C.setImaginary(i, j, B.getImaginary(i, j - A.columns));
							}
						}
					}
					return ComplexMatrix(r, c, C.element);
				}

				//ARITHMETIC METHODS
				//===================
				ComplexMatrix ComplexMatrix::multiply(ComplexMatrix A, ComplexMatrix B) {
					if (A.columns != B.rows) { return ComplexMatrix(); }
					std::vector<std::complex<double>> answer(A.rows*B.columns, std::complex<double>(0));
					std::fill(answer.begin(), answer.end(), std::complex<double>(0)); 

					if (B.columns == 1) {//if B is a column vector, use a faster method for calculation
						for (int i = 0; i < A.rows; ++i) {
							for (int j = 0; j < A.columns; ++j) {
								answer[i] += A.get(i, j).toComplex() * B.get(j, 0).toComplex();
							}
						}
						return ComplexMatrix(A.rows, B.columns, answer);
					}

					for (int i = 0; i < A.rows; ++i) {
						for (int j = 0; j < B.columns; ++j) {
							for (int k = 0; k < A.columns; ++k)
							{
								answer[(i*B.columns) + j] += (A.get(i, k).toComplex() * B.get(k, j).toComplex());
							}
						}
					}
					return ComplexMatrix(A.rows, B.columns, answer);
				}

				ComplexMatrix ComplexMatrix::multiply(ComplexMatrix A, double B) {
					std::vector<double> answerR(A.size,0), answerC(A.size,0);
					std::fill(answerR.begin(), answerR.end(), 0);
					std::fill(answerC.begin(), answerC.end(), 0);
					for (int i = 0; i < A.rows; ++i) {
						for (int j = 0; j < A.columns; ++j) {
							answerR[(i*A.columns) + j] = A.getReal(i, j) * B;
							answerC[(i*A.columns) + j] = A.getImaginary(i, j) * B;
						}
					}
					return ComplexMatrix(A.rows, A.columns, answerR, answerC);
				}

				ComplexMatrix ComplexMatrix::add(ComplexMatrix A, ComplexMatrix B) {
					if (A.size != B.size) { return ComplexMatrix(); }
					int size = A.rows*A.columns;
					double* answerN = new double[size];
					double* answerC = new double[size];
					for (int i = 0; i < A.rows; ++i) {
						for (int j = 0; j < B.columns; ++j) {
							for (int k = 0; k < A.columns; ++k)
							{
								answerN[(i*A.columns) + j] = A.getReal(i, k) + B.getReal(k, j);
								answerC[(i*A.columns) + j] = A.getImaginary(i, k) + B.getImaginary(k, j);

							}
						}
					}
					return ComplexMatrix(B.rows, A.columns, answerN, answerC);
				}

				ComplexMatrix ComplexMatrix::add(ComplexMatrix A, double B) {
					int size = A.rows*A.columns;
					double* answerN = new double[size];
					double* answerC = new double[size];
					for (int i = 0; i < A.rows; ++i) {
						for (int j = 0; j < A.columns; ++j) {
							{
								answerN[(i*A.columns) + j] = A.getReal(i, j) + B;
								answerC[(i*A.columns) + j] = A.getImaginary(i, j);
							}
						}
					}
					return ComplexMatrix(A.rows, A.columns, answerN, answerC);
				}

				ComplexMatrix ComplexMatrix::subtract(ComplexMatrix A, ComplexMatrix B) {
					if (A.size != B.size) { return ComplexMatrix(); }
					int size = A.rows*A.columns;
					double* answerN = new double[size];
					double* answerC = new double[size];
					for (int i = 0; i < size; ++i) { answerN[i] = 0; }

					for (int i = 0; i < A.rows; ++i) {
						for (int j = 0; j < B.columns; ++j) {
							{
								answerN[(i*A.columns) + j] += A.getReal(i, j) - B.getReal(i, j);
								answerC[(i*A.columns) + j] += A.getImaginary(i, j) - B.getImaginary(i, j);
							}
						}
					}
					return ComplexMatrix(B.rows, A.columns, answerN, answerC);
				}

				ComplexMatrix ComplexMatrix::subtract(ComplexMatrix A, double B) {
					int size = A.rows*A.columns;
					double* answerN = new double[size];
					double* answerC = new double[size];
					for (int i = 0; i < A.rows; ++i) {
						for (int j = 0; j < A.columns; ++j) {
							{
								answerN[(i*A.columns) + j] = A.getReal(i, j) - B;
								answerC[(i*A.columns) + j] = A.getImaginary(i, j);
							}
						}
					}
					return ComplexMatrix(A.rows, A.columns, answerN, answerC);
				}

				ComplexMatrix ComplexMatrix::exponent(ComplexMatrix A, int b) {
					int size = A.rows*A.columns;
					double* answerN = new double[size];
					double* answerC = new double[size];
					for (int i = 0; i < A.rows; ++i) {
						for (int j = 0; j < A.columns; ++j) {
							double n = 1;
							double n2 = 1;
							for (int m = 0; m < b; ++m) {
								n *= A.getReal(i, j);
								n2 *= A.getImaginary(i, j);
							}
							answerN[(i*A.columns) + j] = n;
							answerC[(i*A.columns) + j] = n;
						}
					}
					return ComplexMatrix(A.rows, A.columns, answerN, answerC);
				}

				//Overloaded operators:
				ComplexMatrix& ComplexMatrix::operator*=(ComplexMatrix& rhs) { *this = multiply(*this, rhs); return*this; }
				ComplexMatrix& ComplexMatrix::operator*=(double x) { *this = multiply(*this, x); return *this; }
				ComplexMatrix& ComplexMatrix::operator*(ComplexMatrix& rhs) { return multiply(*this, rhs); }
				ComplexMatrix& ComplexMatrix::operator*(double x) { return multiply(*this, x); }

				ComplexMatrix& ComplexMatrix::operator+=(ComplexMatrix& rhs) { *this = add(*this, rhs); return*this; }
				ComplexMatrix& ComplexMatrix::operator+=(double x) { *this = add(*this, x); return *this; }
				ComplexMatrix& ComplexMatrix::operator+(ComplexMatrix& rhs) { return add(*this, rhs); }
				ComplexMatrix& ComplexMatrix::operator+(double x) { return add(*this, x); }

				ComplexMatrix& ComplexMatrix::operator-=(ComplexMatrix& rhs) { *this = subtract(*this, rhs); return*this; }
				ComplexMatrix& ComplexMatrix::operator-=(double x) { *this = subtract(*this, x); return *this; }
				ComplexMatrix& ComplexMatrix::operator-(ComplexMatrix& rhs) { return subtract(*this, rhs); }
				ComplexMatrix& ComplexMatrix::operator-(double x) { return subtract(*this, x); }

				ComplexMatrix& ComplexMatrix::operator^(double x) { return exponent(*this, x); }

				bool ComplexMatrix::operator==(ComplexMatrix& B) {
					if (rows != B.rows || columns != B.columns) { return false; }
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							if (getReal(i, j) != B.getReal(i, j)) { return false; }
							if (getImaginary(i, j) != B.getImaginary(i, j)) { return false; }
						}
					}
					return true;
				}

				std::wostream& operator<<(std::wostream& os, ComplexMatrix& rhs) {
					os << rhs.toString();
					return os;
				}
				//========================


				//BOOLEAN FUNCTIONS
				//=================
				bool ComplexMatrix::isIntegerMatrix() {
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							if (floor(abs(getReal(i, j))) != abs(getReal(i, j))) { return false; }
							if (floor(abs(getImaginary(i, j))) != abs(getImaginary(i, j))) { return false; }
						}
					}
					complexMatrixPrecision = 0;
					return true;
				}

				bool ComplexMatrix::isSymmetric() {
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							if (getReal(i, j) != getReal(j, i)) { return false; }
							if (getImaginary(i, j) != getImaginary(j, i)) { return false; }
						}
					}
					return true;
				}

				bool ComplexMatrix::isBoolean() {
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							if (getReal(i, j) != 0 && getReal(i, j) != 1) { return false; }
							if (getImaginary(i, j) != 0 && getImaginary(i, j) != 1) { return false; }
						}
					}
					return true;
				}

				bool ComplexMatrix::isDiagonalized() {
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							if (i != j) {
								if (abs(getReal(i, j)) > 0.00001) { return false; }//only 4 digits of the mantissa are required for complexMatrixPrecision
								if (abs(getImaginary(i, j)) > 0.00001) { return false; }//only 4 digits of the mantissa are required for complexMatrixPrecision
							}
						}
					}
					return true;
				}

				bool ComplexMatrix::isHermitian() {
					if (*this == conjugateTranspose()) { return true; }
					return false;
				}

				/*	bool ComplexMatrix::isOrthoNormal() {
				for (int i = 0; i < columns; ++i) {
				if (columnNorm(i) != 1) { return false; }
				}
				if (det() == 0) {
				std::wcout << det() << std::endl;
				return false;
				}
				return true;
				}

				bool ComplexMatrix::isLinearlyIndependent() {
				if (det() == 0) { return false; }
				return true;
				}*/
				//================================

				ComplexMatrix ComplexMatrix::transpose() {
					double* p = new double[columns*rows];
					double* q = new double[columns*rows];
					if (columns > rows) { p = new double[columns*columns]; q = new double[columns*columns]; }
					if (columns < rows) { p = new double[rows*rows]; q = new double[rows*rows]; }
					for (int i = 0; i < columns; ++i) {
						for (int j = 0; j < rows; ++j) {
							p[i*rows + j] = getReal(j, i);
							q[i*rows + j] = getImaginary(j, i);
						}
					}
					return ComplexMatrix(columns, rows, p, q);
				}

				ComplexMatrix ComplexMatrix::conjugateTranspose() {
					double* p = new double[columns*rows];
					double* q = new double[columns*rows];
					if (columns > rows) { p = new double[columns*columns]; q = new double[columns*columns]; }
					if (columns < rows) { p = new double[rows*rows]; q = new double[rows*rows]; }
					for (int i = 0; i < columns; ++i) {
						for (int j = 0; j < rows; ++j) {
							p[i*rows + j] = getReal(j, i);
							q[i*rows + j] = getImaginary(j, i);
						}
					}
					return ComplexMatrix(columns, rows, q, p);
				}

				ComplexMatrix ComplexMatrix::outerProduct(ComplexMatrix A, ComplexMatrix B) {//tensore product that creates for every element of Matrix A, a submatrix of B * (A_i_j)
					if (A.columns != B.rows) { return ComplexMatrix(); }
					std::vector<ComplexNumber> ans;
					for (int i = 0; i < A.size*B.size; ++i) {
						ans.push_back(ComplexNumber(0));
					}
					int rw = A.rows*B.rows;
					int cl = A.columns*B.columns;

					for (int i = 0; i < rw; ++i) {
						for (int j = 0; j < cl; ++j) {
							ans[i*cl + j] = B.element[((i%B.rows)*B.columns) + (j%B.columns)] * A.element[floor(i / B.rows)*A.columns + floor(j / B.columns)];
						}
					}
					return ComplexMatrix(rw, cl, ans);
				}

				std::wstring ComplexMatrix::toString() {
					std::wostringstream s;
					int defaultLength = 2;
					bool tf = isIntegerMatrix();
					if (tf == true) { defaultLength -= 2; }
					if (tf == false) { defaultLength++; }
					for (int i = 0; i < rows; ++i) {
						for (int j = 0; j < columns; ++j) {
							int adj = 0;
							if (element[i*columns + j].Re < 0) { ++adj; }
							if (element[i*columns + j].Im < 0) { ++adj; }
							std::wostringstream temp;
							temp << std::fixed << std::setprecision(complexMatrixPrecision) << element[i*columns + j].Re;
							if (element[i*columns + j].Im >= 0) {
								temp << L'+';
							}
							temp << std::fixed << std::setprecision(complexMatrixPrecision) << element[i*columns + j].Im;
							temp << L"i";
							std::wstring tempStr = temp.str();
							temp.clear();

							s <<
								std::fixed <<
								std::left <<
								std::setw(defaultLength + (complexMatrixPrecision * 2) + (largestElementR)+largestElementC + 6) <<
								std::setprecision(complexMatrixPrecision) <<
								tempStr;
						}
						s << L"\r\n";
					}
					std::wstring temp = s.str();
					s.clear();
					return temp;
				}

				void ComplexMatrix::printToFile(std::wstring filename) {
					std::wofstream outf(filename, std::ios_base::app);
					outf << toString() << L"\r\n";
					outf.close();
				}

				void ComplexMatrix::display() { std::wcout << toString() << std::endl; }
	//================END ComplexMatrix CLASS
			
				DualNumberMatrix::DualNumberMatrix() {}
				DualNumberMatrix::~DualNumberMatrix() { element.clear(); }

				DualNumberMatrix::DualNumberMatrix(int row, int col) {
						largestElementR = largestElementD = 1;
						for (int i = 0; i < row; ++i) {
							for (int j = 0; j < col; ++j) {
								element.push_back(DualNumber(0, 0));
							}
						}
						rows = row; columns = col; DualNumberMatrixPrecision = 1;
					}

				DualNumberMatrix::DualNumberMatrix(int row, int col, std::vector<DualNumber> element2) {
						largestElementR = largestElementD = 1;
						for (int i = 0; i < row; ++i) {
							for (int j = 0; j < col; ++j) {
								DualNumber n = element2[i*col + j];
								element.push_back(n);
								if (getLength(n.a) > largestElementR) { largestElementR = getLength(n.a); }
								if (getLength(n.b) > largestElement) { largestElementD = getLength(n.b); }
							}
						}
						rows = row; columns = col; DualNumberMatrixPrecision = 4;
					}

					DualNumber DualNumberMatrix::get(int n) { return element[n]; }
					DualNumber DualNumberMatrix::get(int i, int j) { return element[i*columns + j]; }
					double DualNumberMatrix::getReal(int n) { return element[n].a; }
					double DualNumberMatrix::getReal(int i, int j) { return element[(i*columns) + j].a; }
					double DualNumberMatrix::getDual(int n) { return element[n].b; }
					double DualNumberMatrix::getDual(int i, int j) { return element[(i*columns) + j].b; }
					void DualNumberMatrix::set(int n, DualNumber d) { element[n] = d; }
					void DualNumberMatrix::set(int i, int j, DualNumber d) { element[i*columns + j] = d; }


					bool DualNumberMatrix::isDualNull() {//test if dual matrix has only real values
						for (int i = 0; i < rows; ++i) {
							for (int j = 0; j < columns; ++j) {
								if (floor(abs(getReal(i, j))) != abs(getReal(i, j))) { return false; }
								if (floor(abs(getDual(i, j))) != abs(getDual(i, j))) { return false; }
							}
						}
						DualNumberMatrixPrecision = 0;
						return true;
					}

					//class methods
					//=============
					void DualNumberMatrix::identity() {//turns the matrix into an identity matrix
						for (int i = 0; i<rows; ++i) {
							for (int j = 0; j<columns; ++j) {
								if (i == j) { element[i*columns + j].a = 1; element[i*columns + j].b = 0; }
								else { element[i*columns + j].a = 0; element[i*columns + j].b = 0; }
							}
						}
						largestElementR = largestElementD = 0;
						DualNumberMatrixPrecision = 0;
					}

					int DualNumberMatrix::size() { return rows*columns; }

					DualNumberMatrix DualNumberMatrix::submatrix(int i, int j, int i2, int j2) {
						double* vals = new double[(i2*j2) - (i*j)];
						double* valsC = new double[(i2*j2) - (i*j)];
						int counter = 0;
						for (int a = i; a < i2; ++a) {
							for (int b = j; b < j2; ++b) {
								vals[counter] = getReal(a, b);
								valsC[counter] = getDual(a, b);
								counter++;
							}
						}
						std::vector<DualNumber> d;
						for (int i = 0; i < counter; ++i) { d.push_back(DualNumber(vals[i], valsC[i])); }
						return DualNumberMatrix(i2 - i, j2 - j, d);
					}

					DualNumberMatrix DualNumberMatrix::transpose() {
						std::vector<DualNumber> p;
						for (int i = 0; i < rows*columns; ++i) { p.push_back(DualNumber(0, 0)); }
						for (int i = 0; i < columns; ++i) {
							for (int j = 0; j < rows; ++j) {
								p[i*rows + j] = get(j, i);
							}
						}
						return DualNumberMatrix(columns, rows, p);
					}

					DualNumberMatrix DualNumberMatrix::conjugateTranspose() {
						std::vector<DualNumber> p;
						for (int i = 0; i < rows*columns; ++i) { p.push_back(DualNumber(0, 0)); }
						for (int i = 0; i < columns; ++i) {
							for (int j = 0; j < rows; ++j) {
								p[i*rows + j] = get(j, i).conjugate();
							}
						}
						return DualNumberMatrix(columns, rows, p);
					}

					//ARITHMETIC METHODS
					//===================
					DualNumberMatrix DualNumberMatrix::multiply(DualNumberMatrix A, DualNumberMatrix B) {
						if (A.columns != B.rows) { return DualNumberMatrix(); }
						int size = A.rows*B.columns;
						std::vector<DualNumber> answer;
						for (int i = 0; i < size; ++i) { answer.push_back(DualNumber(0, 0)); }

						for (int i = 0; i < A.rows; ++i) {
							for (int j = 0; j < B.columns; ++j) {
								for (int k = 0; k < A.columns; ++k) {
									answer[(i*B.columns) + j] += A.get(i, k).multiply(A.get(i, k), B.get(k, j));
								}
							}
						}
						return DualNumberMatrix(A.rows, B.columns, answer);
					}

					DualNumberMatrix DualNumberMatrix::multiply(DualNumberMatrix A, double B) {
						int size = A.rows*A.columns;
						std::vector<DualNumber> answer;
						for (int i = 0; i < size; ++i) { answer.push_back(DualNumber(0, 0)); }

						for (int i = 0; i < A.rows; ++i) {
							for (int j = 0; j < A.columns; ++j) {
								for (int k = 0; k < A.columns; ++k) {
									answer[(i*A.columns) + j] += A.get(i, k).multiply(A.get(i, k), B);
								}
							}
						}
						return DualNumberMatrix(A.rows, A.columns, answer);
					}

					DualNumberMatrix DualNumberMatrix::add(DualNumberMatrix A, DualNumberMatrix B) {
						if (A.size() != B.size()) { return DualNumberMatrix(); }
						int size = A.rows*A.columns;
						std::vector<DualNumber> answer;
						for (int i = 0; i < size; ++i) { answer.push_back(DualNumber(0, 0)); }

						for (int i = 0; i < A.rows; ++i) {
							for (int j = 0; j < B.columns; ++j) {
								answer[(i*A.columns) + j] = A.get(i, j).add(A.get(i, j), B.get(i, j));
							}
						}
						return DualNumberMatrix(B.rows, A.columns, answer);
					}

					DualNumberMatrix DualNumberMatrix::add(DualNumberMatrix A, double B) {
						int size = A.rows*A.columns;
						std::vector<DualNumber> answer;
						for (int i = 0; i < size; ++i) { answer.push_back(DualNumber(0, 0)); }

						for (int i = 0; i < A.rows; ++i) {
							for (int j = 0; j < A.columns; ++j) {
								answer[(i*A.columns) + j] = A.get(i, j).add(A.get(i, j), B);
							}
						}
						return DualNumberMatrix(A.rows, A.columns, answer);
					}

					DualNumberMatrix DualNumberMatrix::subtract(DualNumberMatrix A, DualNumberMatrix B) {
						if (A.size() != B.size()) { return DualNumberMatrix(); }
						int size = A.rows*A.columns;
						std::vector<DualNumber> answer;
						for (int i = 0; i < size; ++i) { answer.push_back(DualNumber(0, 0)); }

						for (int i = 0; i < A.rows; ++i) {
							for (int j = 0; j < B.columns; ++j) {
								answer[(i*A.columns) + j] = A.get(i, j).subtract(A.get(i, j), B.get(i, j));
							}
						}
						return DualNumberMatrix(B.rows, A.columns, answer);
					}

					DualNumberMatrix DualNumberMatrix::subtract(DualNumberMatrix A, double B) {
						int size = A.rows*A.columns;
						std::vector<DualNumber> answer;
						for (int i = 0; i < size; ++i) { answer.push_back(DualNumber(0, 0)); }

						for (int i = 0; i < A.rows; ++i) {
							for (int j = 0; j < A.columns; ++j) {
								answer[(i*A.columns) + j] = A.get(i, j).subtract(A.get(i, j), B);
							}
						}
						return DualNumberMatrix(A.rows, A.columns, answer);
					}

					DualNumberMatrix DualNumberMatrix::exponent(DualNumberMatrix A, int b) {
						int size = A.rows*A.columns;
						std::vector<DualNumber> answer;
						for (int i = 0; i < size; ++i) { answer.push_back(DualNumber(0, 0)); }

						for (int i = 0; i < A.rows; ++i) {
							for (int j = 0; j < A.columns; ++j) {
								DualNumber d = A.get(i, j);
								for (int m = 1; m < b; ++m) { d = d.multiply(d, A.get(i, j)); }
								answer[(i*A.columns) + j] = d;
							}
						}
						return DualNumberMatrix(A.rows, A.columns, answer);
					}

					//Overloaded operators:
					DualNumberMatrix DualNumberMatrix::operator*=(DualNumberMatrix rhs) { *this = multiply(*this, rhs); return*this; }
					DualNumberMatrix DualNumberMatrix::operator*=(double x) { *this = multiply(*this, x); return *this; }
					DualNumberMatrix DualNumberMatrix::operator*(DualNumberMatrix rhs) { return multiply(*this, rhs); }
					DualNumberMatrix DualNumberMatrix::operator*(double x) { return multiply(*this, x); }

					DualNumberMatrix DualNumberMatrix::operator+=(DualNumberMatrix rhs) { *this = add(*this, rhs); return*this; }
					DualNumberMatrix DualNumberMatrix::operator+=(double x) { *this = add(*this, x); return *this; }
					DualNumberMatrix DualNumberMatrix::operator+(DualNumberMatrix rhs) { return add(*this, rhs); }
					DualNumberMatrix DualNumberMatrix::operator+(double x) { return add(*this, x); }

					DualNumberMatrix DualNumberMatrix::operator-=(DualNumberMatrix rhs) { *this = subtract(*this, rhs); return*this; }
					DualNumberMatrix DualNumberMatrix::operator-=(double x) { *this = subtract(*this, x); return *this; }
					DualNumberMatrix DualNumberMatrix::operator-(DualNumberMatrix rhs) { return subtract(*this, rhs); }
					DualNumberMatrix DualNumberMatrix::operator-(double x) { return subtract(*this, x); }

					DualNumberMatrix DualNumberMatrix::operator^(double x) { return exponent(*this, x); }

					bool DualNumberMatrix::operator==(DualNumberMatrix& B) {
						if (rows != B.rows || columns != B.columns) { return false; }
						for (int i = 0; i < rows; ++i) {
							for (int j = 0; j < columns; ++j) {
								if (getReal(i, j) != B.getReal(i, j)) { return false; }
								if (getDual(i, j) != B.getDual(i, j)) { return false; }
							}
						}
						return true;
					}

					std::wostream& operator<<(std::wostream& os, DualNumberMatrix& rhs) {
						os << rhs.toString();
						return os;
					}


					std::wstring DualNumberMatrix::toString() {
						std::wostringstream s;
						int defaultLength = 2;
						bool tf = isDualNull();
						if (tf == true) { defaultLength -= 2; }
						if (tf == false) { defaultLength++; }
						for (int i = 0; i < rows; ++i) {
							for (int j = 0; j < columns; ++j) {
								int adj = 0;
								if (getReal(i*columns + j) < 0) { ++adj; }
								if (getDual(i*columns + j) < 0) { ++adj; }
								std::wostringstream temp;
								temp << std::fixed << std::setprecision(DualNumberMatrixPrecision) << getReal(i*columns + j);
								if (getDual(i*columns + j) >= 0) {
									temp << L'+';
								}
								temp << std::fixed << std::setprecision(DualNumberMatrixPrecision) << getDual(i*columns + j);
								temp << L"ε";
								std::wstring tempStr = temp.str();
								temp.clear();

								s <<
									std::fixed <<
									std::left <<
									std::setw(defaultLength + (DualNumberMatrixPrecision * 2) + (largestElementR)+largestElementD + 6) <<
									std::setprecision(DualNumberMatrixPrecision) <<
									tempStr;
							}
							s << L"\n";
						}
						std::wstring temp = s.str();
						s.clear();
						return temp;
					}

					void DualNumberMatrix::printToFile(std::wstring filename) {
						std::wofstream outf(filename, std::ios_base::app);
						outf << toString() << L"\n";
						outf.close();
					}

					void DualNumberMatrix::display() { std::wcout << toString() << std::endl; }
	//==============================================
			

	Matrix directionMatrix(Vector v) {
		int sz = v.size() + 1;
		Matrix A(sz, sz);
		A.identity();
		for (int i = 0; i < v.size(); ++i) { A.set(i, i, v.vec[i]); }
		return A;
	}

	Matrix identityMatrix(int n) { return identityMatrix(n, n); }
	Matrix identityMatrix(int n, int m) {
		std::vector<double> vals(n*m);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				vals[i*m + j] = 0;
				if (i == j) { vals[i*m + j] = 1; }
			}
		}
		return Matrix(n,m,vals);
	}

	Matrix modelMatrix(Matrix T, Matrix R, Matrix S) {
		//model matrix = TranslationMatrix * (RotationMatrix * ScalingMatrix)
		Matrix Rot = R;
		Matrix Trans = T;
		Rot *= S;
		Trans *= Rot;
		return T;
	}

	Matrix positionMatrix(Vector v) {
		int sz = v.size() + 1;
		Matrix A(sz, sz);
		for (int i = 0; i < v.size(); ++i) { A.set(i, i, v.vec[i]); }
		return A;
	}

	Matrix rotationMatrix(std::vector<double> v) {
		//Input rotation vector should be of form <angle,x,y,z> for 3D rotation, or <angle,x,y> for 2D.
		//If x, y, or z == 1, then it is a rotation around that axis.  Only angle and one coordinate should be non-zero
		if (v.size() < 3) { return Matrix(); }
		if (v.size() == 3) {//2D rotation matrix
			Matrix M(2, 2);
			M.set(0, 0, cos(v[0]));
			M.set(0, 1, -sin(v[0]));
			M.set(1, 0, sin(v[0]));
			M.set(1, 1, cos(v[0]));
			return M;
		}
		if (v.size() == 4) {//3D rotation matrix (contained in a 4x4 homogeneous coordinate matrix)
			if (v[1] == 1) {//x-axis rotation
				Matrix M(4, 4);
				M.identity();
				M.set(1, 1, cos(v[0]));
				M.set(1, 2, -sin(v[0]));
				M.set(2, 1, sin(v[0]));
				M.set(2, 2, cos(v[0]));
				return M;
			}
			if (v[2] == 1) {//y-axis rotation
				Matrix M(4, 4);
				M.identity();
				M.set(0, 0, cos(v[0]));
				M.set(0, 2, -sin(v[0]));
				M.set(2, 0, sin(v[0]));
				M.set(2, 2, cos(v[0]));
				return M;
			}
			if (v[3] == 1) {//z-axis rotation
				Matrix M(4, 4);
				M.identity();
				M.set(0, 0, cos(v[0]));
				M.set(0, 1, -sin(v[0]));
				M.set(1, 0, sin(v[0]));
				M.set(1, 1, cos(v[0]));
				return M;
			}
		}
		return Matrix();
	}

	Matrix rotationMatrix(Quaternion v) {
		if (v.element[0] == 1) {//x-axis rotation
			Matrix M(4, 4);
			M.identity();
			M.set(1, 1, cos(v.element[0]));
			M.set(1, 2, -sin(v.element[0]));
			M.set(2, 1, sin(v.element[0]));
			M.set(2, 2, cos(v.element[0]));
			return M;
		}
		if (v.element[1] == 1) {//y-axis rotation
			Matrix M(4, 4);
			M.identity();
			M.set(0, 0, cos(v.element[0]));
			M.set(0, 2, -sin(v.element[0]));
			M.set(2, 0, sin(v.element[0]));
			M.set(2, 2, cos(v.element[0]));
			return M;
		}
		if (v.element[2] == 1) {//z-axis rotation
			Matrix M(4, 4);
			M.identity();
			M.set(0, 0, cos(v.element[0]));
			M.set(0, 1, -sin(v.element[0]));
			M.set(1, 0, sin(v.element[0]));
			M.set(1, 1, cos(v.element[0]));
			return M;
		}
	}

	Matrix scalingMatrix(std::vector<double> v) {
		int sz = v.size() + 1;
		Matrix A(sz, sz);
		A.identity();
		for (int i = 0; i < v.size(); ++i) { A.set(i, i, v[i]); }
		return A;
	}

	Matrix translationMatrix(std::vector<double> v) {
		if (v.size() <= 1) { return Matrix(); }
		if (v.size() == 2) {}
		if (v.size() == 3) {
			Matrix A(4, 4);
			A.identity();
			for (int i = 0; i < v.size(); ++i) { A.set(i, A.columns - 1, v[i]); }
			return A;
		}
		return Matrix();
	}

	Matrix Vandermonde(Polynomial p) {
		//create a Vandermonde matrix from the vector of coefficients from a_n-1 to a0,
		//which is suitable for finding the complex eigenvalues of the polynomial
		Matrix CM = companionMatrix(p);
		std::vector<double> co = toSTLVector(p.coefficient, p.terms);
		Matrix V((p.terms - 1), p.terms - 1);
		for (int i = 0; i < p.terms - 1; ++i) {
			for (int j = 0; j < p.terms - 1; ++j) {
				if (i == 0) { V.set(i, j, 1); }
				else { V.set(i, j, pow(co[i], j + 1)); }
			}
		}
		return V;
	}