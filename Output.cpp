#pragma once
#include "stdafx.h"

std::wstring printOutput() {
	std::wstring op;

	//LIST FUNCTIONS FIRST
	//====================
	if (FunctionContainer.size() > 0) {
		for (int i = 0; i < FunctionContainer.size(); ++i) {    //list functions
			op.append(L"f");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L"(");
			op.append(FunctionContainer[i].variablesToString());
			op.append(L") = ");
			op.append(FunctionContainer[i].toString());
			op.append(L"\r\n");
			if (i == FunctionContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}
	if (ComplexFunctionContainer.size() > 0) {
		for (int i = 0; i < ComplexFunctionContainer.size(); ++i) {    //list complex functions
			op.append(L"cf");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L"(");
			op.append(ComplexFunctionContainer[i].variablesToString());
			op.append(L") = ");
			op.append(ComplexFunctionContainer[i].toString());
			op.append(L"\r\n");
			if (i == ComplexFunctionContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (ParametricFunctionContainer.size() > 0) {
		for (int i = 0; i < ParametricFunctionContainer.size(); ++i) {    //list Parametric functions
			op.append(L"pf");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L"(");
			op.append(ParametricFunctionContainer[i].variablesToString());
			op.append(L") = ");
			op.append(ParametricFunctionContainer[i].toString());
			op.append(L"\r\n");
			if (i == ParametricFunctionContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (VectorValuedFunctionContainer.size() > 0) {
		for (int i = 0; i < VectorValuedFunctionContainer.size(); ++i) {    //list vector-valued functions
			op.append(L"F");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L"(");
			op.append(VectorValuedFunctionContainer[i].variablesToString());
			op.append(L") = ");
			op.append(VectorValuedFunctionContainer[i].toString());
			op.append(L"\r\n");
			if (i == VectorValuedFunctionContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}
	if (InfiniteSeriesContainer.size() > 0) {
		for (int i = 0; i < InfiniteSeriesContainer.size(); ++i) {    //list infinite series
			op.append(L"IS");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			if (InfiniteSeriesContainer[i].isProductSeries == true) { op.append(L" = ∏ "); }
			else { op.append(L" = Σ "); }
			op.append(InfiniteSeriesContainer[i].toString());
			op.append(L"\r\n");
			if (i == InfiniteSeriesContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}
	if (ComplexInfiniteSeriesContainer.size() > 0) {
		for (int i = 0; i < ComplexInfiniteSeriesContainer.size(); ++i) {    //list complex infinite series
			op.append(L"CS");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L" = Σ ");
			op.append(ComplexInfiniteSeriesContainer[i].toString());
			op.append(L"\r\n");
			if (i == ComplexInfiniteSeriesContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	//LIST NUMBERS
	//============
	if (ComplexNumberContainer.size() > 0) {
		for (int i = 0; i < ComplexNumberContainer.size(); ++i) { //list complex numbers
			op.append(L"z");
			std::wostringstream strs;
			strs << i + 1;
			op.append(strs.str());
			op.append(L" = ");
			op.append(ComplexNumberContainer[i].toString());
			op.append(L"\r\n");
			if (i == ComplexNumberContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (DualNumberContainer.size() > 0) {
		for (int i = 0; i < DualNumberContainer.size(); ++i) { //list dual numbers
			op.append(L"d");
			std::wostringstream strs;
			strs << i + 1;
			op.append(strs.str());
			op.append(L" = ");
			op.append(DualNumberContainer[i].toString());
			op.append(L"\r\n");
			if (i == DualNumberContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (GroupContainer.size() > 0) {
		for (int j = 0; j < GroupContainer.size(); j++) {
			op.append(L"G");
			std::wostringstream strs2;
			strs2 << j + 1;
			op.append(strs2.str());
			op.append(L" = ");
			op.append(GroupContainer[j].toString());
			op.append(L"\r\n");
			if (j == GroupContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (OctonianContainer.size() > 0) {
		for (int j = 0; j < OctonianContainer.size(); j++) {
			op.append(L"o");
			std::wostringstream strs2;
			strs2 << j + 1;
			op.append(strs2.str());
			op.append(L" = ");
			op.append(OctonianContainer[j].toString());
			op.append(L"\r\n");
			if (j == OctonianContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (ODEContainer.size() > 0) {
		for (int j = 0; j < ODEContainer.size(); j++) {
			op.append(L"ODE");
			std::wostringstream strs2;
			strs2 << j + 1;
			op.append(strs2.str());
			op.append(L": ");
			op.append(ODEContainer[j].lhs);
			op.append(L" = ");
			op.append(ODEContainer[j].rhs);
			op.append(L"\r\n");
			if (j == ODEContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (PDEContainer.size() > 0) {
		for (int j = 0; j < PDEContainer.size(); j++) {
			op.append(L"PDE");
			std::wostringstream strs2;
			strs2 << j + 1;
			op.append(strs2.str());
			op.append(L": ");
			op.append(PDEContainer[j].lhs);
			op.append(L" = ");
			op.append(PDEContainer[j].rhs);
			op.append(L"\r\n");
			if (j == PDEContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (PolynomialContainer.size() > 0) {
		for (int i = 0; i < PolynomialContainer.size(); ++i) {    //list polynomials
			op.append(L"p");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L"(");
			std::wostringstream strs;
			strs << PolynomialContainer[i].variable;
			op.append(strs.str());
			op.append(L") = ");
			op.append(PolynomialContainer[i].toString());
			op.append(L"\r\n");
			if (i == PolynomialContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (QuaternionContainer.size() > 0) {
		for (int j = 0; j < QuaternionContainer.size(); j++) {
			op.append(L"q");
			std::wostringstream strs2;
			strs2 << j + 1;
			op.append(strs2.str());
			op.append(L" = ");
			op.append(QuaternionContainer[j].toString());
			op.append(L"\r\n");
			if (j == QuaternionContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (DualQuaternionContainer.size() > 0) {
		for (int j = 0; j < DualQuaternionContainer.size(); j++) {
			op.append(L"dq");
			std::wostringstream strs2;
			strs2 << j + 1;
			op.append(strs2.str());
			op.append(L" = ");
			op.append(DualQuaternionContainer[j].toString());
			op.append(L"\r\n");
			if (j == DualQuaternionContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (VectorContainer.size() > 0) {
		for (int j = 0; j < VectorContainer.size(); j++) {
			op.append(L"v");
			std::wostringstream strs2;
			strs2 << j + 1;
			op.append(strs2.str());
			op.append(L" = ");
			op.append(VectorContainer[j].toString());
			op.append(L"\r\n");
			if (j == VectorContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	//PRINT MATRICES LAST
	//===================
	if (SparseMatrixContainer.size() > 0) {
		for (int i = 0; i < SparseMatrixContainer.size(); ++i) {    //list sparse matrices
			op.append(L"SM");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L" = \r\n");
			op.append(SparseMatrixContainer[i].toString());
			op.append(L"\r\n");
		}
	}

	if (TensorContainer.size() > 0) {
		for (int j = 0; j < TensorContainer.size(); j++) {
			op.append(L"T");
			std::wostringstream strs2;
			strs2 << j + 1;
			op.append(strs2.str());
			op.append(L" = ");
			op.append(TensorContainer[j].toString());
			op.append(L"\r\n");
			if (j == TensorContainer.size() - 1) { op.append(L"\r\n"); }
		}
	}

	if (FunctionMatrixContainer.size() > 0) {
		for (int i = 0; i < FunctionMatrixContainer.size(); ++i) {    //list matrices
			op.append(L"FM");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L" = \r\n");
			op.append(FunctionMatrixContainer[i].toString());
			op.append(L"\r\n");
		}
	}

	if (ComplexMatrixContainer.size() > 0) {
		for (int i = 0; i < ComplexMatrixContainer.size(); ++i) {    //list complex matrices
			op.append(L"CM");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L" = \r\n");
			op.append(ComplexMatrixContainer[i].toString());
			op.append(L"\r\n");
		}
	}

	if (DualNumberMatrixContainer.size() > 0) {
		for (int i = 0; i < DualNumberMatrixContainer.size(); ++i) {    //list dual matrices
			op.append(L"DM");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L" = \r\n");
			op.append(DualNumberMatrixContainer[i].toString());
			op.append(L"\r\n");
		}
	}

	if (MatrixContainer.size() > 0) {
		for (int i = 0; i < MatrixContainer.size(); ++i) {    //list matrices
			op.append(L"M");
			std::wostringstream strs2;
			strs2 << i + 1;
			op.append(strs2.str());
			op.append(L" = \r\n");
			op.append(MatrixContainer[i].toString());
			op.append(L"\r\n");
		}
	}

	return op;
}

void ClearContainers() {

	if (ComplexMatrixContainer.size() > 0) {
		for (int i = 0; i < ComplexMatrixContainer.size(); ++i) {
			ComplexMatrixContainer[i].element.clear();
		}
		ComplexMatrixContainer.clear();
	}

	if (ComplexNumberContainer.size() > 0) {
		for (int i = 0; i < ComplexNumberContainer.size(); ++i) {
			ComplexNumberContainer[i].Re = ComplexNumberContainer[i].Im = NULL;
		}
		ComplexNumberContainer.clear();
	}

	if (FunctionContainer.size() > 0) {
		for (int i = 0; i < FunctionContainer.size(); ++i) {
			FunctionContainer[i].clear();
		}
		FunctionContainer.clear();
	}

	if (VectorValuedFunctionContainer.size() > 0) {
		for (int i = 0; i < VectorValuedFunctionContainer.size(); ++i) {
			VectorValuedFunctionContainer[i].clear();
		}
		VectorValuedFunctionContainer.clear();
	}

	if (GroupContainer.size() > 0) {
		for (int i = 0; i < GroupContainer.size(); ++i) {
			GroupContainer[i].element.clear();
		}
		GroupContainer.clear();
	}

	if (MatrixContainer.size() > 0) {
		MatrixContainer.clear();
	}
	if (ODEContainer.size() > 0) {
		for (int i = 0; i < ODEContainer.size(); ++i) {
			ODEContainer[i].lhs = L"";
			ODEContainer[i].rhs = L"";
		}
		ODEContainer.clear();
	}

	if (PDEContainer.size() > 0) {
		for (int i = 0; i < PDEContainer.size(); ++i) {
			PDEContainer[i].lhs = L"";
			PDEContainer[i].rhs = L"";
		}
		PDEContainer.clear();
	}

	if (PolynomialContainer.size() > 0) {
		for (int i = 0; i < PolynomialContainer.size(); ++i) {
			PolynomialContainer[i].coefficient = NULL;
		}
		PolynomialContainer.clear();
	}

	if (QuaternionContainer.size() > 0) {
		for (int i = 0; i < QuaternionContainer.size(); ++i) {
			QuaternionContainer[i].element.clear();
		}
		QuaternionContainer.clear();
	}

	if (SparseMatrixContainer.size() > 0) {
		for (int i = 0; i < SparseMatrixContainer.size(); ++i) {
			SparseMatrixContainer[i].element.clear();
		}
		SparseMatrixContainer.clear();
	}

	if (TensorContainer.size() > 0) {
		for (int j = 0; j < TensorContainer.size(); j++) {
			TensorContainer[j].element.clear();
		}
		TensorContainer.clear();
	}

	if (VectorContainer.size() > 0) {
		for (int j = 0; j < VectorContainer.size(); j++) {
			VectorContainer[j].clear();
		}
		VectorContainer.clear();
	}

}



void loadData() {
	std::wstring line;
	std::wifstream myfile(filename);
	if (myfile) {
		while (getline(myfile, line)) {
			//Functions==================
			if (line[0] == L'f') {
				std::wstring newline;
				newline = L"f";
				newline.append(line.begin() + 2, line.end());
				std::vector<std::wstring> a = ParseFunction(newline);
				Function a1(a[0], a[1]);
				FunctionContainer.push_back(a1);
			}
			//Complex Functions==========
			if (line[0] == L'c' && line[1] == L'f') {
				std::wstring newline;
				newline = L"cf";
				newline.append(line.begin() + 2, line.end());
				std::vector<std::wstring> a = ParseFunction(newline);
				ComplexFunction z(a[0], a[1]);
				ComplexFunctionContainer.push_back(z);
			}
			//Vector-Valued Functions====
			if (line[0] == L'F') {
				std::wstring newline;
				newline = L"F";
				newline.append(line.begin() + 2, line.end());
				VectorValuedFunction z(ParseNFunctions(newline));
				VectorValuedFunctionContainer.push_back(z);
			}
			//Parametric Functions=======
			if (line[0] == L'p' && line[1] == L'f') {
				std::wstring newline;
				newline = L"pf";
				newline.append(line.begin() + 2, line.end());
				std::vector<std::wstring> a = ParseFunction(newline);
				ParametricFunction p(a[0], a[1]);
				ParametricFunctionContainer.push_back(p);
			}
			//Polynomials==========
			if (line[0] == L'p') {
				std::wstring newline;
				newline = L"p";
				newline.append(line.begin() + 2, line.end());
				Polynomial a1 = ParsePolynomial(newline);
				if (a1.terms >= 1) {
					PolynomialContainer.push_back(a1);
				}
			}
			//Vectors=================
			if (line[0] == L'v') {
				std::wstring newline;
				newline.append(line.begin() + 5, line.end());
				std::vector<double> vec2 = ParseVector(newline);
				Vector vec3(vec2);
				VectorContainer.push_back(vec3);
			}
			//Groups=================
			if (line[0] == L'G') {
				std::wstring newline;
				newline.append(line.begin() + 5, line.end());
				std::vector<int> vec2 = ParseGroup(newline);
				Group vec3(vec2);
				GroupContainer.push_back(vec3);
			}
			//Quaternion=================
			if (line[0] == L'q') {
				std::wstring newline;
				newline.append(line.begin() + 5, line.end());
				std::vector<double> vec2 = ParseVector(newline);
				Quaternion vec3(vec2);
				QuaternionContainer.push_back(vec3);
			}
			//Complex Numbers================
			if (line[0] == L'z') {
				std::wstring newline;
				newline = L"z =";
				newline.append(line.begin() + 5, line.end());
				ComplexNumber z = ParseComplexNumber(newline);
				ComplexNumberContainer.push_back(z);
			}
			//ODEs=========================
			if (line[0] == L'O' && line[1] == L'D' && line[2] == L'E') {
				ODE ode1 = ParseODE(line);
				ODEContainer.push_back(ode1);
			}
			//PDEs=========================
			if (line[0] == L'P' && line[1] == L'D' && line[2] == L'E') {
				PDE pde1 = ParsePDE(line);
				PDEContainer.push_back(pde1);
			}
			//Octonian=================
			if (line[0] == L'O' && line[1] != 'D') {//keep in mind there might be some overlap with 'ODE' declarations
				std::wstring newline;
				newline.append(line.begin() + 5, line.end());
				std::vector<double> vec2 = ParseVector(newline);
				Octonian vec3(vec2);
				OctonianContainer.push_back(vec3);
			}
			//Groups=================
			if (line[0] == L'G') {
				std::wstring newline;
				newline.append(line.begin() + 5, line.end());
				std::vector<int> vec2 = ParseGroup(newline);
				Group vec3(vec2);
				GroupContainer.push_back(vec3);
			}
			//Tensors=================
			if (line[0] == L'T') {
				std::wstring newline = line;
				Tensor t = ParseTensor(newline);
				TensorContainer.push_back(t);
			}

			//Matrices================
		topOfLoop:
			if (line[0] == L'M') {
				bool isBoolean = false;
				std::vector<std::wstring> finalStr;
				std::wstring nextLine;
				getline(myfile, nextLine);
				while (!myfile.eof() && nextLine.find(L"M") == std::wstring::npos) {   //while file isn't at end, and there's not a M in the next line
					if (line.find(L"Boolean") != std::wstring::npos) { isBoolean = true; }
					finalStr.push_back(nextLine);
					getline(myfile, nextLine);
				}
				if (nextLine.find(L"M") == std::wstring::npos && nextLine != L"") {
					finalStr.push_back(nextLine);
				}


				std::wstring f;
				//now format the string -- insert a ',' after the end of each number, then remove all the spaces in the string
				for (int k = 0; k < finalStr.size(); ++k) {
					bool firstSpace = true;
					for (int i = 0;i < finalStr[k].length();++i) {      //now format the entire string
						if (finalStr[k][i] == L' ' && firstSpace == true) { //the first space between matrix entries ought to be turned into a ','
							finalStr[k].erase(finalStr[k].begin() + i);
							finalStr[k].insert(i, L",");
							firstSpace = false;
						}
						if (isNumber(finalStr[k][i]) == true && finalStr[k][i] != L'.') {    //remove any spaces after the 1st space
							firstSpace = true;
						}
					}
					if (finalStr[k] != L"") {
						finalStr[k] = removeSpaces(finalStr[k]);
						while (finalStr[k][finalStr[k].length() - 1] == L' ') { finalStr[k].pop_back(); }
						if (finalStr[k][finalStr[k].length() - 1] == L',') { finalStr[k].pop_back(); }
						while (finalStr[k][finalStr[k].length() - 1] == L' ') { finalStr[k].pop_back(); }
						finalStr[k].append(L";");
						f.append(finalStr[k]);
						finalStr[k].clear();
					}
				}
				f.insert(0, L"[");
				f.append(L"]");
				Matrix M = ParseMatrix(f);
				if (isBoolean == true) { M.modulus = 2; }
				MatrixContainer.push_back(M);
				f.clear();
				finalStr.clear();
				line = nextLine;
				isBoolean = false;
				goto topOfLoop;
			}


			//Sparse Matrices================
		topOfLoopSP:
			if (line[0] == L'S' && line[1] == L'M') {
				std::vector<std::wstring> finalStr;
				std::wstring nextLine;
				getline(myfile, nextLine);
				while (!myfile.eof() && nextLine.find(L"SM") == std::wstring::npos) {   //while file isn't at end, and there's not a M in the next line
					finalStr.push_back(nextLine);
					getline(myfile, nextLine);
				}
				if (nextLine.find(L"SM") == std::wstring::npos && nextLine != L"") { finalStr.push_back(nextLine); }


				std::wstring f;
				//now format the string -- insert a ',' after the end of each number, then remove all the spaces in the string
				for (int k = 0; k < finalStr.size(); ++k) {
					bool firstSpace = true;
					for (int i = 0;i < finalStr[k].length();++i) {      //now format the entire string
						if (finalStr[k][i] == L' ' && firstSpace == true) { //the first space between matrix entries ought to be turned into a ','
							finalStr[k].erase(finalStr[k].begin() + i);
							finalStr[k].insert(i, L",");
							firstSpace = false;
						}
						if (isNumber(finalStr[k][i]) == true && finalStr[k][i] != '.') {    //remove any spaces after the 1st space
							firstSpace = true;
						}
					}
					if (finalStr[k] != L"") {
						finalStr[k] = removeSpaces(finalStr[k]);
						while (finalStr[k][finalStr[k].length() - 1] == L' ') { finalStr[k].pop_back(); }
						if (finalStr[k][finalStr[k].length() - 1] == L',') { finalStr[k].pop_back(); }
						while (finalStr[k][finalStr[k].length() - 1] == L' ') { finalStr[k].pop_back(); }
						finalStr[k].append(L";");
						f.append(finalStr[k]);
						finalStr[k].clear();
					}
				}
				f.insert(0, L"[");
				f.append(L"]");
				Matrix M = ParseMatrix(f);
				MatrixContainer.push_back(M);
				f.clear();
				finalStr.clear();
				line = nextLine;
				goto topOfLoopSP;
			}





		}
		myfile.close();
	}

	else {};
}


void clearSavedData() {//deletes all data saved in the .txt file
	std::wofstream outf(filename);
	outf.open(filename, std::wofstream::out | std::wofstream::trunc);
	outf.close();
}

void saveData() {
	clearSavedData();       //empty out file first to make sure data isn't double written
	std::wofstream outf2(filename, std::ios_base::app);//ios_base::app will always ONLY append data

	if (FunctionContainer.size() > 0) {
		for (int i = 0; i < FunctionContainer.size(); ++i) {    //list functions
			outf2 << L"f" << i + 1 << L"(" << FunctionContainer[i].variablesToString() << L") = ";
			outf2 << FunctionContainer[i].toString();
			outf2 << L"\n";
		}
		outf2 << L"\n";
	}

	if (ComplexFunctionContainer.size() > 0) {
		for (int i = 0; i < ComplexFunctionContainer.size(); ++i) {    //list complex-valued functions
			outf2 << L"cf" << i + 1 << L"(" << ComplexFunctionContainer[i].variablesToString() << L") = ";
			outf2 << ComplexFunctionContainer[i].toString();
			outf2 << L"\n";
		}
		outf2 << L"\n";
	}

	if (ParametricFunctionContainer.size() > 0) {
		for (int i = 0; i < ParametricFunctionContainer.size(); ++i) {    //list complex-valued functions
			outf2 << L"pf" << i + 1 << L"(" << ParametricFunctionContainer[i].variablesToString() << L") = ";
			outf2 << ParametricFunctionContainer[i].toString();
			outf2 << L"\n";
		}
		outf2 << L"\n";
	}

	if (VectorValuedFunctionContainer.size() > 0) {
		for (int i = 0; i < VectorValuedFunctionContainer.size(); ++i) {    //list VectorValued functions
			outf2 << L"F" << i + 1 << L"(" << VectorValuedFunctionContainer[i].variablesToString() << L") = ";
			outf2 << VectorValuedFunctionContainer[i].toString();
			outf2 << L"\n";
		}
		outf2 << L"\n";
	}

	if (PolynomialContainer.size() > 0) {
		for (int i = 0; i < PolynomialContainer.size(); ++i) {    //list polynomials
			if (PolynomialContainer[i].terms >= 1) {
				outf2 << L"p" << i + 1 << L"(" << PolynomialContainer[i].variable << L") = ";
				outf2 << PolynomialContainer[i].toString();
				outf2 << L"\n";
			}
		}
		outf2 << L"\n";
	}
	if (VectorContainer.size() > 0) {
		for (int j = 0; j < VectorContainer.size(); j++) {
			outf2 << L"v" << j + 1 << L" = ";
			outf2 << VectorContainer[j].toString() << L"\n";
		}
		outf2 << L"\n";
	}
	if (GroupContainer.size() > 0) {
		for (int j = 0; j < GroupContainer.size(); j++) {
			outf2 << L"G" << j + 1 << L" = ";
			outf2 << GroupContainer[j].toString() << L"\n";
		}
		outf2 << L"\n";
	}
	if (QuaternionContainer.size() > 0) {
		for (int j = 0; j < QuaternionContainer.size(); j++) {
			outf2 << L"q" << j + 1 << L" = ";
			outf2 << QuaternionContainer[j].toString() << L"\n";
		}
		outf2 << L"\n";
	}
	if (ComplexNumberContainer.size() > 0) {
		for (int i = 0; i < ComplexNumberContainer.size(); ++i) { //list complex numbers
			outf2 << L"z" << i + 1 << L" = " << ComplexNumberContainer[i].toString() << L"\n";
		}
		outf2 << L"\n";
	}
	if (ODEContainer.size() > 0) {
		for (int i = 0; i < ODEContainer.size(); ++i) { //list ODEs
			outf2 << L"ODE" << i + 1 << L": " << ODEContainer[i].lhs << L" = " << ODEContainer[i].rhs;
			outf2 << L"\n";
		}
		outf2 << L"\n";
	}

	if (PDEContainer.size() > 0) {
		for (int i = 0; i < PDEContainer.size(); ++i) { //list PDEs
			outf2 << L"PDE" << i + 1 << L": " << PDEContainer[i].lhs << L" = " << PDEContainer[i].rhs;
			outf2 << L"\n";
		}
		outf2 << L"\n";
	}

	if (TensorContainer.size() > 0) {
		for (int i = 0; i < TensorContainer.size(); ++i) { //list tensors
			outf2 << L"T" << i + 1 << L"(";
			for (int j = 0; j < TensorContainer[i].index.size() - 1; ++j) { outf2 << std::to_wstring(TensorContainer[i].index[j]) << L","; }
			outf2 << std::to_wstring(TensorContainer[i].index[TensorContainer[i].index.size() - 1]);
			outf2 << L") = " << TensorContainer[i].toString() << L"\n";
		}
		outf2 << L"\n";
	}

	if (ComplexMatrixContainer.size() > 0) {
		for (int i = 0; i < ComplexMatrixContainer.size(); ++i) { //list complex matrices
			outf2 << L"CM" << i + 1 << L" = \n" << ComplexMatrixContainer[i].toString() << L"\n";
		}
		outf2 << L"\n";
	}
	if (MatrixContainer.size() > 0) {
		for (int i = 0; i < MatrixContainer.size(); ++i) { //list matrices
			outf2 << L"M" << i + 1 << L" = ";
			if (MatrixContainer[i].modulus == 2) { outf2 << L"Boolean "; }
			outf2 << L"\n" << MatrixContainer[i].toString() << L"\n";
		}
	}

	if (SparseMatrixContainer.size() > 0) {
		for (int i = 0; i < SparseMatrixContainer.size(); ++i) { //list sparse matrices
			outf2 << L"SM" << i + 1 << L" = \n" << SparseMatrixContainer[i].toString() << L"\n";
		}
	}
	//===================================
	outf2.close();
	output.append(L"Saved data printed to file '");
	output.append(filename);
	output.append(L"'.\r\n\r\n");
}