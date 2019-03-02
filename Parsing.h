#ifndef PARSING_H
#define PARSING_H
#pragma once
#include "stdafx.h"

std::vector<std::wstring> tokenizeByComma(std::wstring str);
std::vector<std::wstring> tokenizeByChar(std::wstring str, wchar_t token);
double ParseInputNumber(std::string in);
double ParseInputNumber(std::wstring in);
std::vector<double> ParseNInputNumbers(std::string str);
std::vector<double> ParseNInputNumbers(std::wstring str);
ComplexNumber ParseComplexNumber(std::wstring in);
std::vector<ComplexNumber> ParseNComplexNumbers(std::wstring str);
Matrix ParseMatrix(std::wstring inp);
Matrix ParseMatrixFromFile(std::wstring input);
void ParseMatrixLoadFile(std::wstring filename);
ComplexMatrix ParseComplexMatrix(std::wstring input);
DualNumberMatrix ParseDualNumberMatrix(std::wstring input);
Tensor ParseTensor(std::wstring in);
Tensor NEWParseTensor(std::wstring in);
ComplexPolynomial ParseComplexPolynomial(std::wstring input);
DualNumber ParseDualNumber(std::wstring in);
std::vector<int> ParseGroup(std::wstring input);
ODE ParseODE(std::wstring in);
PDE ParsePDE(std::wstring in);
std::vector<double> ParseVector(std::wstring input);
Quaternion ParseQuaternion(std::wstring input);
DualQuaternion ParseDualQuaternion(std::wstring input);
SparseMatrix ParseSparseMatrix(std::wstring inp);
SparseMatrix ParseSparseMatrixFromFile(std::wstring input);
Polynomial ParsePolynomial(std::wstring input);
std::vector<std::wstring> ParseFunction(std::wstring in);
FunctionMatrix ParseFunctionMatrix(std::wstring inp);
std::vector<Function> ParseNFunctions(std::wstring str);
std::vector<wchar_t> varsToChars(std::wstring vars);
#endif