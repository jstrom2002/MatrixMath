#pragma once
#include "stdafx.h"
#include "MainParser.h"

std::wstring input;  //holds string from input text box
std::wstring output; //holds display to be shown to user

std::vector<ComplexFunction> ComplexFunctionContainer;
std::vector<ComplexInfiniteSeries> ComplexInfiniteSeriesContainer;
std::vector<ComplexMatrix> ComplexMatrixContainer;
std::vector<ComplexNumber> ComplexNumberContainer;
std::vector<ComplexPolynomial> ComplexPolynomialContainer;
std::vector<ComplexPolynomialFraction> ComplexPolynomialFractionContainer;
std::vector<DirectedGraph> DirectedGraphContainer;
std::vector<DualNumber> DualNumberContainer;
std::vector<DualNumberMatrix> DualNumberMatrixContainer;
std::vector<DualNumberPolynomial> DualNumberPolynomialContainer;
std::vector<DualQuaternion> DualQuaternionContainer;
std::vector<Function> FunctionContainer;
std::vector<FunctionMatrix> FunctionMatrixContainer;
std::vector<Group> GroupContainer;
std::vector<InfiniteSeries> InfiniteSeriesContainer;
std::vector<Matrix> MatrixContainer;
std::vector<Octonian> OctonianContainer;
std::vector<ODE> ODEContainer;
std::vector<ParametricFunction> ParametricFunctionContainer;
std::vector<PDE> PDEContainer;
std::vector<Polynomial> PolynomialContainer;
std::vector<PolynomialFraction> PolynomialFractionContainer;
std::vector<Quaternion> QuaternionContainer;
std::vector<SparseMatrix> SparseMatrixContainer;
std::vector<Tensor> TensorContainer;
std::vector<UndirectedGraph> UndirectedGraphContainer;
std::vector<Vector> VectorContainer;
std::vector<VectorValuedFunction> VectorValuedFunctionContainer;

void MainParser(std::wstring in) {
	//TESTING AREA
	//in = L"((1+0i)/((1+0i)+(1+0i)))*(2.71828+0i)";
/*	std::vector<ComplexNumber> a;
	a.push_back(ComplexNumber(11));
	a.push_back(ComplexNumber(22));
	a.push_back(ComplexNumber(23));
	a.push_back(ComplexNumber(-2));
	a.push_back(ComplexNumber(0.2));
	a.push_back(ComplexNumber(-0.1));
	a = fastFourierTransform(a);	*/
	
/*	double p[16] = {
		1,2,2,0,
		0,1,3,1,
		0,1,2,1,
		1,2,2,-1
	};
	double pn[36] = {
		1.0,0.8,1.0,1.0,0.8,1.0,
		1.0,0.5,0.3,0.0,0.5,1.0,
		1.0,0.3,0.2,0.0,0.3,1.0,
		1.0,0.2,0.0,0.0,0.2,1.0,
		1.0,0.3,0.2,0.0,0.3,1.0,
		1.0,0.8,1.0,1.0,0.8,1.0
	};

	std::vector<double> p2 = toSTLVector(p, 16);
	//std::vector<double> p2 = toSTLVector(pn, 36);
	Matrix A(4, 4, p2);
	//Matrix A(6, 6, p2);
	Matrix Mat = discreteCosineTransform(A);
	output.append(A.toString());
	output.append(L"\r\n");
	output.append(Mat.toString());
	output.append(L"\r\n\r\n");
	Matrix Mat2 = inverseDiscreteCosineTransform(Mat);
	output.append(Mat2.toString());
	output.append(L"\r\n\r\n");*/

	//IMAGE bmp(L"cat.bmp");
/*	Matrix Mat = testMatrix();
	IMAGE bmp(Mat);
	Matrix Mat2 = GaussianKernel2D(5,1);
	output.append(Mat2.toString());
	bmp.transpose();
	bmp.saveImage(L"test2.bmp");
	BMPtoText(L"test2.bmp");	
*/	

/*	std::vector<double> v1,v2;
	v1.push_back(2);
	v1.push_back(2);
	v1.push_back(1);
	v1.push_back(1);
	v2.push_back(-3);
	v2.push_back(-3);
	v2.push_back(-2);
	v2.push_back(-2);
	Line l(v1,v2);
	output.append(l.toString());
	output.append(L"\r\n");

	v1.clear(); v2.clear();
	v1.push_back(-2);
	v1.push_back(2);
	v2.push_back(2);
	v2.push_back(0);
	Line l2(v1, v2);
	output.append(l2.toString());
	output.append(L"\r\n");

	double temp[16] = {
		1, -2, 2, 4,
		2, -4, 5, 9,
		3, -6, 8, 14,
		5, -10, 12, 22
	};
	Matrix Mat(4, 4, toSTLVector(temp,16));
//	std::vector<double> vec = Mat.solve();
//	Matrix Mat2 = Matrix(vec.size(),1,vec);
	output.append(Mat.nullSpace().toString());
//	output.append(Mat2.toString());	
*/


/*	double temp[6] = {
		3,2,2,
		2,3,-2
	};

	Matrix Mat1(2, 3, toSTLVector(temp,6));
	std::vector<Matrix> Mats = Mat1.SVD();
	output.append(Mats[0].toString());
	output.append(L"\r\n");
	output.append(Mats[1].toString());
	output.append(L"\r\n");
	output.append(Mats[2].toString());//solution: http://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm
*/

	MatrixTest();


/*	output.append(L"\r\n");
	Matrix Mat2 = Mat.multiply(MatInv, Mat);
	output.append(Mat2.toString());
*/
	//std::vector<double> intersection = getIntersection(l, l2); //should be 1.09,0.45
	//output.append(toString(intersection));

	//Matrix Mat = bmp.getBlackAndWhiteMatrix();
	//TensorContainer.push_back(loadImageColor(L"test.bmp"));
	//Matrix Mat = TensorContainer[0].getSubmatrix(0);
	//output.append(Mat.toString());
	//ComplexMatrix fftMat = discreteFourierTransform(Mat);
	//output.append(fftMat.toString());
	output.append(L"\r\n\r\n");
	return;
	//============
	
	if (in == L"") { return; }//don't process empty strings
	in = removeSpaces(in);//remove all spaces in input to eliminate ambiguity

	 //COMMANDS RELATED TO ENVIROMENT VARIABLES
	 //========================================
	if (in == L"clear") { ClearContainers(); return; }
	if (in == L"clearFile") { clearSavedData(); return; }
	if (in == L"exit" || in == L"quit") { ClearContainers(); running = false; return; }
	if (in.find(L"filename(") != std::wstring::npos) {
		filename.clear();
		filename.append(in.substr(9));
		if (filename[filename.length() - 1] == ')') { filename.pop_back(); }
		if (in.find(L".txt") == std::wstring::npos) {
			filename.append(L".txt");
		}
		output.append(L"Save/load file name has been set to ");
		output.append(filename);
		output.append(L"\r\n\r\n");
		return;
	}
	if (in == L"fullscreen") {
		if (fullscreen == false) {
			fullscreen = true;
			changeScreen = true;
			output.append(L"fullscreen: ON\r\n\r\n");
			return;
		}
		if (fullscreen == true) {
			fullscreen = false;
			changeScreen = true;
			output.append(L"fullscreen: OFF\r\n\r\n");
			return;
		}
		return;
	}
	if (in == L"load") { 
		loadData(); 
		if (first == false) { output.append(L"Data loaded.\r\n\r\n"); } 
		return; 
	}
	if (in.find(L"loadImage(")) {		
		std::wstring filename = (in.substr(in.find(L"loadImage(")));
		if (filename[filename.length() - 1] == L')') { filename.pop_back(); }
		IMAGE img(filename);
		Tensor M = img.toTensor();
		TensorContainer.push_back(M);
		return;
	}
	if (in.find(L"precision(") != std::wstring::npos) {
		std::wstring tmp = in.substr(10);
		if (tmp[tmp.length() - 1] == L')') { tmp.pop_back(); }
		precision = stringToDouble(tmp);

		if (MatrixContainer.size() > 0) {//apply precision to all non-integer matrices
			for (int i = 0; i < MatrixContainer.size(); ++i) {
				if (MatrixContainer[i].isIntegerMatrix() == false) { MatrixContainer[i].matrixPrecision = precision; }
			}
		}
		if (ComplexMatrixContainer.size() > 0) {//apply precision to all non-integer complex matrices
			for (int i = 0; i < ComplexMatrixContainer.size(); ++i) {
				if (ComplexMatrixContainer[i].isIntegerMatrix() == false) { ComplexMatrixContainer[i].complexMatrixPrecision = precision; }
			}
		}

		if (DualNumberMatrixContainer.size() > 0) {//apply precision to all non-integer DualNumber matrices
			for (int i = 0; i < DualNumberMatrixContainer.size(); ++i) {
				if (DualNumberMatrixContainer[i].isDualNull() == false) { DualNumberMatrixContainer[i].DualNumberMatrixPrecision = precision; }
			}
		}

		output.append(L"System precision has been set to ");
		output.append(tmp);
		output.append(L"\r\n\r\n");
		return;
	}
	if (in == L"productKey") { output.append(L"Product key: " + productKey + L"\r\n\r\n"); return; }
	if (in == L"save") { saveData(); return; }
	if (in == L"serialnumber") { output.append(L"Serial number: " + serialnumber + L"\r\n\r\n"); return; }
	if (in.find(L"setEpsilon(")!=std::wstring::npos) {
		double x = ParseInputNumber(in.substr(in.find(L'(')));
		EPSILON = x;
		output.append(L"System epsilon value has been set to ");
		int prc = precision;	
		precision = getMantissaLength(x);
		output.append(to_stringPrecision(x));
		precision = prc;
		output.append(L"\r\n\r\n");
		return; 
	}
//	if (in == L"test" || in == L"test()") { test(); return; }
	if (in == L"windowed") { fullscreen = false; changeScreen = true; return; }
	//END ENVIRONMENT RELATED COMMANDS SECTION==================================

	//Complex Matrices================
	if (in[0] == 'C' && in[1] == 'M' && in[2] == '[') {
		ComplexMatrix CM = ParseComplexMatrix(in.substr(1));
		ComplexMatrixContainer.push_back(CM);
		return;
	}

	//Complex Functions========
	if (in[0] == 'c' && in[1] == 'f' && in[2] == '(') {
		ComplexFunction z(in);
		ComplexFunctionContainer.push_back(z);
		return;
	}

	//Complex Infinite Series================
	if (in[0] == 'C' && in[1] == 'I' && in[2] == 'S' && in[3] == ':') {
		ComplexFunction f(in.substr(4));
		ComplexInfiniteSeries M(f);
		ComplexInfiniteSeriesContainer.push_back(M);
		return;
	}

	//Complex Numbers================
	if (in[0] == 'z' && in[1] == '=') {
		ComplexNumber z = ParseComplexNumber(in);
		ComplexNumberContainer.push_back(z);
		return;
	}

	//Dual Numbers===============
	if (in[0] == 'd' && in[1] == '=') {
		DualNumber d = ParseDualNumber(in);
		DualNumberContainer.push_back(d);
		return;
	}

	//Dual Matrices================
	if (in[0] == 'D' && in[1] == 'M' && in[2] == '[') {
		DualNumberMatrix CM = ParseDualNumberMatrix(in.substr(1));
		DualNumberMatrixContainer.push_back(CM);
		return;
	}

	//Functions==========
	if (in[0] == 'f' && in[1] == '(') {
		Function f(in);
		FunctionContainer.push_back(f);
		return;
	}

	//Function Matrices================
	if (in[0] == 'F' && in[1] == 'M' && in[2] == '[') {
		FunctionMatrix M = ParseFunctionMatrix(in);
		FunctionMatrixContainer.push_back(M);
		return;
	}

	//Infinite Series================
	if (in[0] == 'I' && in[1] == 'S' && in[2] == ':') {
		Function f(in.substr(3));
		InfiniteSeries M(f);
		InfiniteSeriesContainer.push_back(M);
		return;
	}

	//Parametric Functions==========
	if (in[0] == 'p' && in[1] == 'f' && in[2] == '(') {
		ParametricFunction f(in);
		ParametricFunctionContainer.push_back(f);
		return;
	}

	//Polynomials==========
	if (in[0] == 'p' && in[1] == '(') {
		Polynomial a1 = ParsePolynomial(in);
		if (a1.terms >= 1) {
			PolynomialContainer.push_back(a1);
		}
		return;
	}

	//Complex Polynomials==========
	if (in[0] == 'c' && in[1] == 'p' && in[2] == '(') {
		ComplexPolynomial a1 = ParseComplexPolynomial(in);
		if (a1.length >= 1) {
			ComplexPolynomialContainer.push_back(a1);
		}
		return;
	}

	//Vectors=================
	if (in[0] == '<') {
		std::vector<double> vec2 = ParseVector(in);
		VectorContainer.push_back(Vector(vec2));
		return;
	}

	//Vector-Valued Functions==========
	if (in[0] == 'F' && in[1] == '(') {
		VectorValuedFunction vvf(ParseNFunctions(in));
		VectorValuedFunctionContainer.push_back(vvf);
		return;
	}

	//Quaternions=================
	if (in[0] == 'q' && in[1] == '<') {
		Quaternion vec2 = ParseQuaternion(in);
		QuaternionContainer.push_back(vec2);
		return;
	}

	//Dual Quaternions=================
	if (in[0] == 'd' && in[1] == 'q' && in[2] == '<') {
		DualQuaternion vec2 = ParseDualQuaternion(in);
		DualQuaternionContainer.push_back(vec2);
		return;
	}

	//Matrices================
	if (in[0] == '[' || (in[0] == 'M' && in[1] == '[')) {
		Matrix M = ParseMatrix(in);
		MatrixContainer.push_back(M);
		return;
	}

	//Sparse Matrices================
	if (in[0] == 'S' && in[1] == 'M' && in[2] == '[') {
		SparseMatrix M = ParseSparseMatrix(in);
		SparseMatrixContainer.push_back(M);
		return;
	}

	//ODEs=========================
	if (in[0] == 'O' && in[1] == 'D' && in[2] == 'E' && in[3] == ':') {
		ODE ode1 = ParseODE(in);
		ODEContainer.push_back(ode1);
		return;
	}

	//PDEs=========================
	if (in[0] == 'P' && in[1] == 'D' && in[2] == 'E' && in[3] == ':') {
		PDE pde1 = ParsePDE(in);
		PDEContainer.push_back(pde1);
		return;
	}

	//Tensors==========================
	if (in[0] == 'T' && in[1] == '[') {
		Tensor t = ParseTensor(in);
		TensorContainer.push_back(t);
		return;
	}

	//UndirectedGraphs=======================
	/*	if (in[0] == 'U' && in[1] == 'G') {
	std::wstring tempstr = in.substr(2);

	if (in[2] == '[') {                 //case -- input matrix declaration
	matrix2D M2 = ParseMatrix2D(in);
	Eigen::MatrixXd M = convertMatrix2DtoEigenMatrix(M2);
	graphs.push_back(UndirectedSimpleGraph(M));
	}

	if (in[2] == 'M') {//associate graph with a saved matrix
	double tempIn = ParseInputNumber(tempstr);
	//MessageBox(NULL, to_stringPrecision(tempIn).c_str(), NULL, NULL); //debug code
	graphs.push_back(UndirectedSimpleGraph(A[tempIn - 1]));
	}
	}*/


	//==============================================================================================

	//FUNCTIONS
	//=========
	//constructor
	if (in[0] == 'f' && isNumber(in[1]) == true) {    //reference to saved functions
		int position = 0;
		for (int i = 1; i < in.length(); ++i) {
			if (isNumber(in[i]) == true) { ++position; }
			if (isNumber(in[i]) == false) { i = in.length(); }
		}
		std::wstring a = in.substr(1, position);
		int tempIn = ParseInputNumber(a);//get index of function array
		--tempIn;

		//Class functions for "Function"
		//==============================
		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(7));
			FunctionContainer.push_back(FunctionContainer[tempIn].add(FunctionContainer[temp[0] - 1], FunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			FunctionContainer.push_back(FunctionContainer[tempIn]);
			return;
		}
		if (in.find(L".compose(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(11));
			FunctionContainer.push_back(FunctionContainer[tempIn].compose(FunctionContainer[temp[0] - 1], FunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".convolve(") != std::wstring::npos) {
			std::wstring tempStr = in.substr(12);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<double> temp = ParseNInputNumbers(tempStr);
			double result = 0;
			//case: 1D convolution
			if (FunctionContainer[temp[1] - 1].variables.size() == 1 && FunctionContainer[temp[2] - 1].variables.size() == 1) {
				result = convolve1D(temp[0], FunctionContainer[temp[1] - 1], FunctionContainer[temp[2] - 1]);
			}
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".convolve(");
			output.append(tempStr);
			output.append(L") = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			FunctionContainer[tempIn].clear();
			FunctionContainer.erase(FunctionContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".derivative(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

			if (tempStr.find(L",") != std::wstring::npos) {//case: multiple input values
				std::vector<double> temp = ParseNInputNumbers(tempStr);

				if (temp.size() == 2) {//case: one variable, multiple derivatives -- first value is the number of derivatives
					double result = nthDerivative1D(temp[0], temp[1], FunctionContainer[tempIn]);
					output.append(L"(d/d");
					output.append(charToString(FunctionContainer[tempIn].variables[0]));
					output.append(L")");
					output.append(L"^");
					output.append(to_stringPrecision(temp[0]));
					output.append(L" f");
					output.append(to_stringPrecision(tempIn + 1));
					output.append(L"(");
					output.append(to_stringPrecision(temp[1]));
					output.append(L") = ");
					output.append(to_stringPrecision(result));
					output.append(L"\r\n\r\n");
				}

				if (temp.size() > 2 && FunctionContainer[tempIn].variables.size()>1) {//case: multiple variables, multiple derivatives -- first value is number of derivatives, second is which variable to differentiate with respect to	
					std::vector<double> newTemp;
					for (int i = 2; i < temp.size();++i) {
						newTemp.push_back(temp[i]);
					}   //convert STL vector to pointer array
					double result = nthPartialDerivative(temp[0], temp[1], newTemp, FunctionContainer[tempIn]);
					output.append(L"(a/a");
					output.append(charToString(FunctionContainer[tempIn].variables[0]));
					output.append(L")");
					if (temp[0] != 1) {
						output.append(L"^");
						output.append(to_stringPrecision(temp[0]));
					}
					output.append(L" f");
					output.append(to_stringPrecision(tempIn + 1));
					output.append(L"(");
					tempStr = tempStr.substr(tempStr.find(L",") + 1);
					tempStr = tempStr.substr(tempStr.find(L",") + 1);
					output.append(tempStr);
					output.append(L") = ");
					output.append(to_stringPrecision(result));
					output.append(L"\r\n\r\n");
				}
			}
			else {//case: univariate function, one derivative
				double temp = ParseInputNumber(tempStr);
				std::vector<double> newTemp;
				newTemp.push_back(temp);
				double result = nthDerivative1D(1, newTemp[0], FunctionContainer[tempIn]);
				output.append(L"d/d");
				output.append(charToString(FunctionContainer[tempIn].variables[0]));
				output.append(L" ");
				output.append(L"f");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L"(");
				output.append(tempStr);
				output.append(L") = ");
				output.append(to_stringPrecision(result));
				output.append(L"\r\n\r\n");
			}
			return;
		}

		if (in.find(L".divergence(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

			std::vector<double> temp = ParseNInputNumbers(tempStr);
			if (temp.size() >= 1 && FunctionContainer[tempIn].variables.size() >= 1) {
				double result = divergence(temp, FunctionContainer[tempIn]);
				output.append(L"f");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".divergence(");
				output.append(tempStr);
				output.append(L") = ");
				output.append(to_stringPrecision(result));
				output.append(L"\r\n\r\n");
			}
			return;
		}
		if (in.find(L".divide(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(10));
			FunctionContainer.push_back(FunctionContainer[tempIn].divide(FunctionContainer[temp[0] - 1], FunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".evaluate(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			if (tempStr.find(L",") != std::wstring::npos) {//case: multiple variables
				std::vector<double> temp = ParseNInputNumbers(tempStr);
				double result = FunctionContainer[tempIn].evaluate(temp);
				output.append(L"f");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".evaluate(");
				output.append(tempStr);
				output.append(L") = ");
				output.append(to_stringPrecision(result));
				output.append(L"\r\n\r\n");
			}
			else {//case: univariate function
				double temp = ParseInputNumber(tempStr);
				if (tempIn >= FunctionContainer.size()) { return; }//check to make sure the function index exists
				double result = FunctionContainer[tempIn].evaluate(temp);
				output.append(L"f");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".evaluate(");
				output.append(tempStr);
				output.append(L") = ");
				output.append(to_stringPrecision(result));
				output.append(L"\r\n\r\n");
			}
			return;
		}
		if (in.find(L".exponent(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(10));
			FunctionContainer.push_back(FunctionContainer[tempIn].exponent(FunctionContainer[temp[0] - 1], FunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".findRoots()") != std::wstring::npos) {
			std::vector<double> roots = FunctionContainer[tempIn].findRoots();
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".findRoots() = ");
			if (roots.size() == 0) {
				output.append(L"(no roots found)");
			}else {
				for (int i = 0; i < roots.size(); ++i) {
					output.append(to_stringPrecision(roots[i]));
					if (i < roots.size() - 1) { output.append(L","); }
				}
			}
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".findRoots(") != std::wstring::npos) {
			std::vector<double> vals = ParseNInputNumbers(in.substr(13));
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".findRoots(");
			output.append(in.substr(13));
			output.append(L" = ");
			if(vals.size()==2) {//find roots of univariate function using lo/hi bounds as input
				std::vector<double> roots = FunctionContainer[tempIn].findRoots(vals[0], vals[1]);
				if (roots.size() == 0) {
					output.append(L"(no roots found)");
				}
				for (int i = 0; i < roots.size(); ++i) {
					output.append(to_stringPrecision(roots[i]));
					if (i < roots.size() - 1) { output.append(L","); }
				}
			}
			if (vals.size() > 2) {//find roots of multivariate function. input is (lo bound, hi bound, values for variables, position of variable in array to find root with respect to)
				double lo = vals[0];
				double hi = vals[1];
				double pos = vals[vals.size() - 1];
				vals.erase(vals.begin());
				vals.erase(vals.begin());
				vals.erase(vals.begin()+(vals.size()-1));
				std::vector<double> roots = FunctionContainer[tempIn].findRoots(lo, hi, vals, pos);
				if (roots.size() == 0) {
					output.append(L"(no roots found)");
				}
				for (int i = 0; i < roots.size(); ++i) {
					output.append(to_stringPrecision(roots[i]));
					if (i < roots.size() - 1) { output.append(L","); }
				}
			}
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".FourierTransform(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<double> temp = ParseNInputNumbers(tempStr);
			double result = FourierTransform(temp[0], FunctionContainer[tempIn]);
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".FourierTransform(");
			output.append(tempStr);
			output.append(L") = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");
		}

		if (in.find(L".gradient(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

			if (tempStr.find(L",") != std::wstring::npos) {//case: multiple input values
				std::vector<double> temp = ParseNInputNumbers(tempStr);

				if (temp.size() > 1 && FunctionContainer[tempIn].variables.size() > 1) {//case: multiple variables
					std::vector<double> result = gradient(temp, FunctionContainer[tempIn]);
					output.append(L"f");
					output.append(to_stringPrecision(tempIn + 1));
					output.append(L".gradient(");
					output.append(tempStr);
					output.append(L") = <");
					output.append(arrToString(result));
					output.append(L">\r\n\r\n");
				}
			}
			else {//case: univariate function
				std::vector<double> temp = ParseNInputNumbers(tempStr);
				std::vector<double> result = gradient(temp, FunctionContainer[tempIn]);
				output.append(L"f");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L"(");
				output.append(tempStr);
				output.append(L") = <");
				output.append(arrToString(result));
				output.append(L">\r\n\r\n");
			}
			return;
		}
		if (in.find(L".integral(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

			if (tempStr.find(L",") != std::wstring::npos) {//case: multiple input values
				if (tempStr.find(L"inf") != std::wstring::npos) {//case: improper integral, bounds contain an infinite value
					std::wstring strA = tempStr.substr(0, tempStr.find(L","));
					tempStr = tempStr.substr(tempStr.find(L",") + 1);
					std::wstring strB = tempStr.substr(0, tempStr.find(L")"));
					strA = trim(strA);
					strB = trim(strB);

					if (strB.find(L"inf") != std::wstring::npos) {
						if (strA.find(L"-inf") != std::wstring::npos) {
							double result = integrateGaussHermiteQuadrature(FunctionContainer[tempIn]);
							output.append(L"f");
							output.append(to_stringPrecision(tempIn + 1));
							output.append(L".integral(" + strA + L", " + strB + L") = ");
							output.append(to_stringPrecision(result));
							output.append(L"\r\n\r\n");
						}
						if (strA.find(L"0") != std::wstring::npos) {
							double result = integrateGaussLaguerreQuadrature(FunctionContainer[tempIn]);
							output.append(L"f");
							output.append(to_stringPrecision(tempIn + 1));
							output.append(L".integral(" + strA + L", " + strB + L") = ");
							output.append(to_stringPrecision(result));
							output.append(L"\r\n\r\n");
						}
						if (strA != L"0") {
							double result = integrateGaussLaguerreQuadrature(FunctionContainer[tempIn]);
							double result2 = integrateGaussLegendreQuadrature(0, stringToDouble(strA), FunctionContainer[tempIn]);
							result -= result2;
							output.append(L"f");
							output.append(to_stringPrecision(tempIn + 1));
							output.append(L".integral(" + strA + L", " + strB + L") = ");
							output.append(to_stringPrecision(result));
							output.append(L"\r\n\r\n");
						}
					}
				}

				//case: definite integral of a finite, bounded function
				std::vector<double> temp = ParseNInputNumbers(tempStr);
				if (temp.size() == 2) {//case: one variable definite integral from [a,b]
					double result = integrateGaussLegendreQuadrature(temp[0], temp[1], FunctionContainer[tempIn]);
					output.append(L"f");
					output.append(to_stringPrecision(tempIn + 1));
					output.append(L".integral(");
					output.append(to_stringPrecision(temp[0]));
					output.append(L", ");
					output.append(to_stringPrecision(temp[1]));
					output.append(L") = ");
					output.append(to_stringPrecision(result));
					output.append(L"\r\n\r\n");
				}
			}
			return;
		}
		if (in.find(L".inverse()") != std::wstring::npos) {
			int temp = ParseInputNumber(in.substr(11));
			FunctionContainer.push_back(FunctionContainer[tempIn].inverse(FunctionContainer[tempIn]));
			return;
		}
		if (in.find(L".inverseFourierTransform(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<double> temp = ParseNInputNumbers(tempStr);
			double result = inverseFourierTransform(temp[0], FunctionContainer[tempIn]);
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".inverseFourierTransform(");
			output.append(tempStr);
			output.append(L") = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");
		}

		if (in.find(L".LaplaceTransform(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<double> temp = ParseNInputNumbers(tempStr);
			double result = LaplaceTransform(temp[0], FunctionContainer[tempIn]);
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".LaplaceTransform(");
			output.append(tempStr);
			output.append(L") = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");
		}

		if (in.find(L".MellinTransform(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<double> temp = ParseNInputNumbers(tempStr);
			double result = MellinTransform(temp[0], FunctionContainer[tempIn]);
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".MellinTransform(");
			output.append(tempStr);
			output.append(L") = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");
		}
		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(11));
			FunctionContainer.push_back(FunctionContainer[tempIn].multiply(FunctionContainer[temp[0] - 1], FunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(11));
			FunctionContainer.push_back(FunctionContainer[tempIn].subtract(FunctionContainer[temp[0] - 1], FunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".toComplexFunction()") != std::wstring::npos) {
			ComplexFunctionContainer.push_back(ComplexFunction(FunctionContainer[tempIn]));
			return;
		}
		if (in.find(L".toInfiniteSeries()") != std::wstring::npos) {
			InfiniteSeriesContainer.push_back(InfiniteSeries(FunctionContainer[tempIn]));
			return;
		}
		if (in.find(L".twoSidedLaplaceTransform(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<double> temp = ParseNInputNumbers(tempStr);
			double result = twoSidedLaplaceTransform(temp[0], FunctionContainer[tempIn]);
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".twoSidedLaplaceTransform(");
			output.append(tempStr);
			output.append(L") = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");
		}
		if (in.find(L".twoSidedMellinTransform(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<double> temp = ParseNInputNumbers(tempStr);
			double result = twoSidedMellinTransform(temp[0], FunctionContainer[tempIn]);
			output.append(L"f");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".twoSidedMellinTransform(");
			output.append(tempStr);
			output.append(L") = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");
		}

		if (in.find(L".variables(") != std::wstring::npos) {
			std::wstring tempStr = in.substr(in.find(L"(") + 1);
			while (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<wchar_t> vars;
			for (int i = 0; i < tempStr.size(); ++i) {
				if (isLetter(tempStr[i]) == true) { vars.push_back(tempStr[i]); }
			}
			FunctionContainer[tempIn].variables = vars;
			return;
		}
	}//END FUNCTION HANDLING===============================================================


	 //COMPLEX FUNCTIONS
	 //==================
	 //constructor
	if (in[0] == 'c' && in[1] == 'f' && isNumber(in[2]) == true) {    //reference to saved complex functions
		int position = 0;
		for (int i = 2; i < in.length(); ++i) {
			if (isNumber(in[i]) == true) { ++position; }
			if (isNumber(in[i]) == false) { i = in.length(); }
		}
		std::wstring a = in.substr(2, position);
		int tempIn = ParseInputNumber(a);//get index of function array
		--tempIn;

		//Class functions for "ComplexFunction"
		//==============================
		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(7));
			ComplexFunctionContainer.push_back(ComplexFunctionContainer[tempIn].add(ComplexFunctionContainer[temp[0] - 1], ComplexFunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			ComplexFunctionContainer.push_back(ComplexFunctionContainer[tempIn]);
			return;
		}
		if (in.find(L".compose(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			ComplexFunctionContainer.push_back(ComplexFunction(ComplexFunctionContainer[tempIn].f.compose(ComplexFunctionContainer[temp[0] - 1].f, ComplexFunctionContainer[temp[1] - 1].f)));
			return;
		}
		/*		if (in.find(L".convolve(") != std::wstring::npos) {
		std::wstring tempStr = in.substr(12);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		double result = 0;
		//case: 1D convolution
		if (ComplexFunctionContainer[temp[1] - 1].variables.size() == 1 && ComplexFunctionContainer[temp[2] - 1].variables.size() == 1) {
		result = convolve1D(temp[0], ComplexFunctionContainer[temp[1] - 1], ComplexFunctionContainer[temp[2] - 1]);
		}
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".convolve(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		return;
		}
		*/		if (in.find(L".delete()") != std::wstring::npos) {
			ComplexFunctionContainer[tempIn].clear();
			ComplexFunctionContainer.erase(ComplexFunctionContainer.begin() + tempIn);
			return;
		}
		/*		if (in.find(L".derivative(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

		if (tempStr.find(L",") != std::wstring::npos) {//case: multiple input values
		std::vector<double> temp = ParseNInputNumbers(tempStr);

		if (temp.size() == 2) {//case: one variable, multiple derivatives -- first value is the number of derivatives
		double result = nthDerivative1D(temp[0], temp[1], ComplexFunctionContainer[tempIn]);
		output.append(L"(d/d");
		output.append(charToString(ComplexFunctionContainer[tempIn].variables[0]));
		output.append(L")");
		output.append(L"^");
		output.append(to_stringPrecision(temp[0]));
		output.append(L" f");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L"(");
		output.append(to_stringPrecision(temp[1]));
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}

		if (temp.size() > 2 && ComplexFunctionContainer[tempIn].variables.size()>1) {//case: multiple variables, multiple derivatives -- first value is number of derivatives, second is which variable to differentiate with respect to
		std::vector<double> newTemp;
		for (int i = 2; i < temp.size();++i) {
		newTemp.push_back(temp[i]);
		}   //convert STL vector to pointer array
		double result = nthPartialDerivative(temp[0], temp[1], newTemp, ComplexFunctionContainer[tempIn]);
		output.append(L"(a/a");
		output.append(charToString(ComplexFunctionContainer[tempIn].variables[0]));
		output.append(L")");
		if (temp[0] != 1) {
		output.append(L"^");
		output.append(to_stringPrecision(temp[0]));
		}
		output.append(L" f");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L"(");
		tempStr = tempStr.substr(tempStr.find(L",") + 1);
		tempStr = tempStr.substr(tempStr.find(L",") + 1);
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		}
		else {//case: univariate function, one derivative
		double temp = ParseInputNumber(tempStr);
		std::vector<double> newTemp;
		newTemp.push_back(temp);
		double result = nthDerivative1D(1, newTemp[0], ComplexFunctionContainer[tempIn]);
		output.append(L"d/d");
		output.append(charToString(ComplexFunctionContainer[tempIn].variables[0]));
		output.append(L" ");
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L"(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		return;
		}

		if (in.find(L".divergence(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

		std::vector<double> temp = ParseNInputNumbers(tempStr);
		if (temp.size() >= 1 && ComplexFunctionContainer[tempIn].variables.size() >= 1) {
		double result = divergence(temp, ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".divergence(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		return;
		}
		*/		if (in.find(L".divide(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(10));
			ComplexFunctionContainer.push_back(ComplexFunctionContainer[tempIn].divide(ComplexFunctionContainer[temp[0] - 1], ComplexFunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".evaluate(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			if (tempStr.find(L",") != std::wstring::npos) {//case: multiple variables
				std::vector<ComplexNumber> temp = ParseNComplexNumbers(tempStr);
				ComplexNumber result = ComplexFunctionContainer[tempIn].evaluateComplex(temp);
				output.append(L"cf");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".evaluate(");
				output.append(tempStr);
				output.append(L") = ");
				output.append(result.toString());
				output.append(L"\r\n\r\n");
			}
			else {//case: univariate function
				ComplexNumber temp = ParseComplexNumber(tempStr);
				if (tempIn >= ComplexFunctionContainer.size()) { return; }//check to make sure the function index exists
				ComplexNumber result = ComplexFunctionContainer[tempIn].evaluateComplex(temp);
				output.append(L"cf");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".evaluate(");
				output.append(tempStr);
				output.append(L") = ");
				output.append(result.toString());
				output.append(L"\r\n\r\n");
			}
			return;
		}
		if (in.find(L".exponent(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(10));
			ComplexFunctionContainer.push_back(ComplexFunctionContainer[tempIn].exponent(ComplexFunctionContainer[temp[0] - 1], ComplexFunctionContainer[temp[1] - 1]));
			return;
		}
		/*		if (in.find(L".FourierTransform(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		double result = FourierTransform(temp[0], ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".FourierTransform(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}

		if (in.find(L".gradient(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

		if (tempStr.find(L",") != std::wstring::npos) {//case: multiple input values
		std::vector<double> temp = ParseNInputNumbers(tempStr);

		if (temp.size() > 1 && ComplexFunctionContainer[tempIn].variables.size() > 1) {//case: multiple variables
		std::vector<double> result = gradient(temp, ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".gradient(");
		output.append(tempStr);
		output.append(L") = <");
		output.append(arrToString(result));
		output.append(L">\r\n\r\n");
		}
		}
		else {//case: univariate function
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		std::vector<double> result = gradient(temp, ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L"(");
		output.append(tempStr);
		output.append(L") = <");
		output.append(arrToString(result));
		output.append(L">\r\n\r\n");
		}
		return;
		}
		if (in.find(L".integral(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

		if (tempStr.find(L",") != std::wstring::npos) {//case: multiple input values
		if (tempStr.find(L"inf") != std::wstring::npos) {//case: improper integral, bounds contain an infinite value
		std::wstring strA = tempStr.substr(0, tempStr.find(L","));
		tempStr = tempStr.substr(tempStr.find(L",") + 1);
		std::wstring strB = tempStr.substr(0, tempStr.find(L")"));
		strA = trim(strA);
		strB = trim(strB);

		if (strB.find(L"inf") != std::wstring::npos) {
		if (strA.find(L"-inf") != std::wstring::npos) {
		double result = integrateGaussHermiteQuadrature(ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".integral(" + strA + L", " + strB + L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		if (strA.find(L"0") != std::wstring::npos) {
		double result = integrateGaussLaguerreQuadrature(ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".integral(" + strA + L", " + strB + L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		if (strA != L"0") {
		double result = integrateGaussLaguerreQuadrature(ComplexFunctionContainer[tempIn]);
		double result2 = integrateGaussLegendreQuadrature(0, stringToDouble(strA), ComplexFunctionContainer[tempIn]);
		result -= result2;
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".integral(" + strA + L", " + strB + L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		}
		}

		//case: definite integral of a finite, bounded function
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		if (temp.size() == 2) {//case: one variable definite integral from [a,b]
		double result = integrateGaussLegendreQuadrature(temp[0], temp[1], ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".integral(");
		output.append(to_stringPrecision(temp[0]));
		output.append(L", ");
		output.append(to_stringPrecision(temp[1]));
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		}
		return;
		}
		if (in.find(L".inverse()") != std::wstring::npos) {
		int temp = ParseInputNumber(in.substr(11));
		ComplexFunctionContainer.push_back(ComplexFunctionContainer[tempIn].inverse(ComplexFunctionContainer[tempIn]));
		return;
		}
		if (in.find(L".inverseFourierTransform(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		double result = inverseFourierTransform(temp[0], ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".inverseFourierTransform(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}

		if (in.find(L".LaplaceTransform(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		double result = LaplaceTransform(temp[0], ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".LaplaceTransform(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}

		if (in.find(L".MellinTransform(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		double result = MellinTransform(temp[0], ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".MellinTransform(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		*/		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(11));
			ComplexFunctionContainer.push_back(ComplexFunctionContainer[tempIn].multiply(ComplexFunctionContainer[temp[0] - 1], ComplexFunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(11));
			ComplexFunctionContainer.push_back(ComplexFunctionContainer[tempIn].subtract(ComplexFunctionContainer[temp[0] - 1], ComplexFunctionContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".toFunction()") != std::wstring::npos) {
			FunctionContainer.push_back(Function(ComplexFunctionContainer[tempIn]));
			return;
		}
		if (in.find(L".toInfiniteSeries()") != std::wstring::npos) {
			InfiniteSeriesContainer.push_back(InfiniteSeries(ComplexFunctionContainer[tempIn]));
			return;
		}
		/*		if (in.find(L".twoSidedLaplaceTransform(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		double result = twoSidedLaplaceTransform(temp[0], ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".twoSidedLaplaceTransform(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		if (in.find(L".twoSidedMellinTransform(") != std::wstring::npos) {
		int substringPart = 0;
		for (int i = 0; i < in.length();++i) {
		if (in[i] == '(') {
		substringPart = i + 1;
		}
		}
		std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
		if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
		std::vector<double> temp = ParseNInputNumbers(tempStr);
		double result = twoSidedMellinTransform(temp[0], ComplexFunctionContainer[tempIn]);
		output.append(L"cf");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".twoSidedMellinTransform(");
		output.append(tempStr);
		output.append(L") = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		}
		*/
		if (in.find(L".variables(") != std::wstring::npos) {
			std::wstring tempStr = in.substr(in.find(L"(") + 1);
			while (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<wchar_t> vars;
			for (int i = 0; i < tempStr.size(); ++i) {
				if (isLetter(tempStr[i]) == true) { vars.push_back(tempStr[i]); }
			}
			ComplexFunctionContainer[tempIn].f.variables = vars;
			return;
		}
	}//END COMPLEXFUNCTION HANDLING===============================================================



	 //VECTOR-VALUED FUNCTIONS
	 //=======================
	 //constructor
	if (in[0] == 'F' && isNumber(in[1]) == true) { //reference to saved functions
		int position = 0;
		for (int i = 1; i < in.length(); ++i) {
			if (isNumber(in[i]) == true) { ++position; }
			if (isNumber(in[i]) == false) { i = in.length(); }
		}
		std::wstring a = in.substr(1, position);
		int tempIn = ParseInputNumber(a);//get index of function array
		--tempIn;

		//Class functions for Vector-Valued Functions
		//============================================
		if (in.find(L".copy()") != std::wstring::npos) {
			VectorValuedFunctionContainer.push_back(VectorValuedFunctionContainer[tempIn]);
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			VectorValuedFunctionContainer[tempIn].clear();
			VectorValuedFunctionContainer.erase(VectorValuedFunctionContainer.begin() + tempIn);
			return;
		}

		if (in.find(L".divergence(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

			std::vector<double> temp = ParseNInputNumbers(tempStr);
			if (temp.size() >= 1 && VectorValuedFunctionContainer[tempIn].variables.size() >= 1) {
				double result = divergence(temp, VectorValuedFunctionContainer[tempIn]);
				output.append(L"F");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".divergence(");
				output.append(tempStr);
				output.append(L") = ");
				output.append(to_stringPrecision(result));
				output.append(L"\r\n\r\n");
			}
			return;
		}

		if (in.find(L".evaluate(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

			if (tempStr.find(L",") != std::wstring::npos) {//case: multiple variables
				std::vector<double> temp = ParseNInputNumbers(tempStr);
				std::vector<double> result = VectorValuedFunctionContainer[tempIn].evaluate(temp);
				output.append(L"F");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".evaluate(");
				output.append(tempStr);
				output.append(L") = <");
				output.append(arrToString(result));
				output.append(L">\r\n\r\n");
			}
			else {//case: univariate function
				double temp = ParseInputNumber(tempStr);
				std::vector<double> newTemp;
				newTemp.push_back(temp);
				std::vector<double> result = VectorValuedFunctionContainer[tempIn].evaluate(newTemp);
				output.append(L"F");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".evaluate(");
				output.append(tempStr);
				output.append(L") = <");
				output.append(arrToString(result));
				output.append(L">\r\n\r\n");
			}
			return;
		}
		if (in.find(L".gradient(") != std::wstring::npos) {
			int substringPart = 0;
			for (int i = 0; i < in.length();++i) {
				if (in[i] == '(') {
					substringPart = i + 1;
				}
			}
			std::wstring tempStr = in.substr(substringPart);  //get a substring for whatever values there might be to evaluate.
			if (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }

			if (tempStr.find(L",") != std::wstring::npos) {//case: multiple input values
				std::vector<double> temp = ParseNInputNumbers(tempStr);

				if (temp.size() > 2 && FunctionContainer[tempIn].variables.size() > 1) {//case: multiple variables
					std::vector<double> result = gradient(temp, VectorValuedFunctionContainer[tempIn]);
					output.append(L"F");
					output.append(to_stringPrecision(tempIn + 1));
					output.append(L".gradient(");
					output.append(tempStr);
					output.append(L") = <");
					output.append(arrToString(result));
					output.append(L">\r\n\r\n");
				}
			}
			else {//case: univariate function
				std::vector<double> temp = ParseNInputNumbers(tempStr);
				std::vector<double> result = gradient(temp, VectorValuedFunctionContainer[tempIn]);
				output.append(L"F");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L"(");
				output.append(tempStr);
				output.append(L") = <");
				output.append(to_stringPrecision(result[0]));
				output.append(L">\r\n\r\n");
			}
			return;
		}
		if (in.find(L".size()") != std::wstring::npos) {
			int temp = VectorValuedFunctionContainer[tempIn].size();
			output.append(L"F");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".size() = ");
			output.append(std::to_wstring(temp));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".variables(") != std::wstring::npos) {
			std::wstring tempStr = in.substr(in.find(L"(") + 1);
			while (tempStr.find(L")") != std::wstring::npos) { tempStr.pop_back(); }
			std::vector<wchar_t> vars;
			for (int i = 0; i < tempStr.size(); ++i) {
				if (isLetter(tempStr[i]) == true) { vars.push_back(tempStr[i]); }
			}
			VectorValuedFunctionContainer[tempIn].variables = vars;
			return;
		}

	}//END VECTOR-VALUED FUNCTION HANDLING=================================================


	 //POLYNOMIALS
	 //===========
	 //constructor
	if (in.find(L"Polynomial(") != std::wstring::npos && in.substr(0, 11) == (L"Polynomial(")) {
		std::wstring tempstr = in.substr(11);
		while (tempstr.find(L")") != std::wstring::npos) { tempstr.pop_back(); }
		std::vector<double> temp = ParseNInputNumbers(tempstr);
		std::reverse(temp.begin(), temp.end());
		Polynomial p(temp);
		PolynomialContainer.push_back(p);
		temp.clear();
		return;
	}

	//Polynomial class functions
	//==========================
	if (in[0] == L'p' && isNumber(in[1]) == true) {    //reference to saved Polynomial
		double tempIn = ParseInputNumber(in.substr(1, 2));
		tempIn--;

		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(7));
			PolynomialContainer.push_back(PolynomialContainer[tempIn].subtract(PolynomialContainer[temp[0] - 1], PolynomialContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".checkLowBound(") != std::wstring::npos) {
			double temp = ParseInputNumber(in.substr(17));
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(x) divided by ");
			output.append(to_stringPrecision(temp));
			output.append(L" is a low bound");
			bool tf = PolynomialContainer[tempIn].checkLowBound(temp);
			if (tf == true) { output.append(L" = true."); }
			if (tf == false) { output.append(L" = false."); }
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			PolynomialContainer.push_back(PolynomialContainer[tempIn]);
			return;
		}
		if (in.find(L".companionMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(companionMatrix(PolynomialContainer[tempIn]));
			return;
		}
		if (in.find(L".complexSyntheticDivision(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(28));
			ComplexNumber z(temp[0], temp[1]);
			ComplexNumber e = PolynomialContainer[tempIn].complexSyntheticDivision(z);
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(x) divided by ");
			output.append(z.toString());
			output.append(L" = ");
			output.append(e.toString());
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			PolynomialContainer[tempIn].coefficient = NULL;
			PolynomialContainer.erase(PolynomialContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".derivative()") != std::wstring::npos) {
			Polynomial p = PolynomialContainer[tempIn].derivative();
			PolynomialContainer.push_back(p);
			output.append(L"derivative of p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L" is ");
			output.append(p.toString());
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".evaluate(") != std::wstring::npos) {
			double temp = ParseInputNumber(in.substr(12));
			double answer = PolynomialContainer[tempIn].evaluate(temp);
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(");
			output.append(to_stringPrecision(temp));
			output.append(L") = ");
			output.append(to_stringPrecision(answer));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".evaluateComplex(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(19));
			ComplexNumber z(temp[0], temp[1]);
			ComplexNumber e = PolynomialContainer[tempIn].evaluate(z);
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(");
			output.append(z.toString());
			output.append(L") = ");
			output.append(e.toString());
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".exponentiate(") != std::wstring::npos) {
			int temp = ParseInputNumber(in.substr(16));
			PolynomialContainer.push_back(PolynomialContainer[tempIn].exponentiate(PolynomialContainer[tempIn], temp));
			return;
		}
		if (in.find(L".factorOutBinomial(") != std::wstring::npos) {
			double temp = ParseInputNumber(in.substr(21));
			Polynomial p = PolynomialContainer[tempIn].factorOutBinomial(temp);
			if (p.terms > 0) { PolynomialContainer.push_back(p); }
			else { output.append(L"Could not factor out binomial, only certain values \nlike integers can factor neatly\r\n\r\n"); }
			return;
		}
		if (in.find(L".findRoot(") != std::wstring::npos) {
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(x) = 0 at ");
			output.append(to_stringPrecision(findRoot(PolynomialContainer[tempIn])));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".findRoots(") != std::wstring::npos) {
			std::vector<double> answer = findRealRoots(PolynomialContainer[tempIn]);
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(x) = 0 at ");
			Vector V(answer);
			output.append(V.toString());
			if (V.size() > 0) { VectorContainer.push_back(V); }
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".findRealRoots(") != std::wstring::npos) {
			std::vector<double> answer = findRealRoots(PolynomialContainer[tempIn]);
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(x) = 0 at ");
			Vector V(answer);
			output.append(V.toString());
			if (V.size() > 0) { VectorContainer.push_back(V); }
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".getRootBounds(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(17));
			std::vector<double> bounds = PolynomialContainer[tempIn].getRootBounds();
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L" has roots bounded between [");
			output.append(to_stringPrecision(bounds[0]));
			output.append(L" , ");
			output.append(to_stringPrecision(bounds[1]));
			output.append(L"]\r\n\r\n");
			return;
		}
		if (in.find(L".HornerEvaluate(") != std::wstring::npos) {
			double temp = ParseInputNumber(in.substr(18));
			double answer = PolynomialContainer[tempIn].HornerEvaluate(temp);
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(");
			output.append(to_stringPrecision(temp));
			output.append(L") = ");
			output.append(to_stringPrecision(answer));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".HornerEvaluateComplex(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(25));
			ComplexNumber z(temp[0], temp[1]);
			ComplexNumber answer = PolynomialContainer[tempIn].HornerEvaluateComplex(z);
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(");
			output.append(z.toString());
			output.append(L") = ");
			output.append(answer.toString());
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".integral()") != std::wstring::npos) {
			Polynomial p = PolynomialContainer[tempIn].integral();
			PolynomialContainer.push_back(p);
			output.append(L"integeral of p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L" is ");
			output.append(p.toString());
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".length") != std::wstring::npos) {
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(x) has ");
			output.append(to_stringPrecision(PolynomialContainer[tempIn].terms));
			output.append(L" terms\r\n\r\n");
			return;
		}
		if (in.find(L".makeMonic()") != std::wstring::npos) {
			PolynomialContainer[tempIn].makeMonic();
			return;
		}
		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			PolynomialContainer.push_back(PolynomialContainer[tempIn].multiply(PolynomialContainer[temp[0] - 1], PolynomialContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".randomize()") != std::wstring::npos) {
			PolynomialContainer[tempIn].randomize();
			return;
		}
		if (in.find(L".roots()") != std::wstring::npos) {
			std::vector<ComplexNumber> eigs = roots(PolynomialContainer[tempIn]);
			if (eigs.size() > 0) {
				output.append(L"p");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".roots() = \n");
				for (int i = 0; i < eigs.size(); ++i) {
					output.append(eigs[i].toString());
					if (i < eigs.size() - 1) { output.append(L", "); }
					ComplexNumberContainer.push_back(eigs[i]);
				}
				output.append(L"\r\n\r\n");
			}
			else { output.append(L"roots of the polynomial could not be found.\n"); }
			return;
		}
		if (in.find(L".syntheticDivision(") != std::wstring::npos) {
			double temp = ParseInputNumber(in.substr(20));
			output.append(L"p");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"(x) divided by ");
			output.append(to_stringPrecision(temp));
			output.append(L" has a remainder of ");
			output.append(to_stringPrecision(PolynomialContainer[tempIn].syntheticDivision(temp)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			PolynomialContainer.push_back(PolynomialContainer[tempIn].subtract(PolynomialContainer[temp[0] - 1], PolynomialContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".VandermondeMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(Vandermonde(PolynomialContainer[tempIn]));
			return;
		}


	}//END POLYNOMIAL FUNCTIONS======================================================================================


	 //ComplexNumbers
	 //==============
	 //constructor
	if (in.find(L"ComplexNumber(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(14);
		std::vector<double> temp = ParseNInputNumbers(tempstr);
		ComplexNumber z(temp[0], temp[1]);
		ComplexNumberContainer.push_back(z);
		temp.clear();
		return;
	}

	//assorted class functions of ComplexNumber
	//=========================================
	if (in[0] == 'z' && isNumber(in[1]) == true) {    //reference to saved ComplexNumber
		double tempIn = ParseInputNumber(in.substr(1, 2));
		tempIn--;
		//------------------------------------
		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(7));
			ComplexNumberContainer.push_back(ComplexNumberContainer[tempIn].add(ComplexNumberContainer[temp[0] - 1], ComplexNumberContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			ComplexNumberContainer.push_back(ComplexNumberContainer[tempIn]);
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			ComplexNumberContainer.erase(ComplexNumberContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".exponentiate(") != std::wstring::npos) {
			int temp = ParseInputNumber(in.substr(16));
			ComplexNumberContainer.push_back(ComplexNumberContainer[tempIn].exponent(ComplexNumberContainer[tempIn], temp));
			return;
		}
		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			ComplexNumberContainer.push_back(ComplexNumberContainer[tempIn].multiply(ComplexNumberContainer[temp[0] - 1], ComplexNumberContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".randomize()") != std::wstring::npos) {
			ComplexNumberContainer[tempIn].randomize();
			return;
		}
		if (in.find(L".roots(") != std::wstring::npos) {
			double temp = ParseInputNumber(in.substr(9));
			std::vector<ComplexNumber> c = roots(ComplexNumberContainer[tempIn], temp);
			for (int i = 0; i < c.size(); ++i) {
				ComplexNumberContainer.push_back(c[i]);
			}
			c.clear();
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			ComplexNumberContainer.push_back(ComplexNumberContainer[tempIn].subtract(ComplexNumberContainer[temp[0] - 1], ComplexNumberContainer[temp[1] - 1]));
			return;
		}


	}//END COMPLEXNUMBER FUNCTIONS======================================================================

	 //DualNumber
	 //==============
	 //constructor
	if (in.find(L"DualNumber(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(10);
		std::vector<double> temp = ParseNInputNumbers(tempstr);
		DualNumber d(temp[0], temp[1]);
		DualNumberContainer.push_back(d);
		temp.clear();
		return;
	}

	//assorted class functions of DualNumber
	//=========================================
	if (in[0] == 'd' && isNumber(in[1]) == true) {    //reference to saved DualNumber
		double tempIn = ParseInputNumber(in.substr(1, 2));
		tempIn--;
		//------------------------------------
		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(7));
			DualNumberContainer.push_back(DualNumberContainer[tempIn].subtract(DualNumberContainer[temp[0] - 1], DualNumberContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".conjugate(") != std::wstring::npos) {
			DualNumberContainer.push_back(DualNumberContainer[tempIn].conjugate());
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			DualNumberContainer.push_back(DualNumberContainer[tempIn]);
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			DualNumberContainer.erase(DualNumberContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".exponentiate(") != std::wstring::npos) {
			int temp = ParseInputNumber(in.substr(16));
			DualNumberContainer.push_back(DualNumberContainer[tempIn].exponent(DualNumberContainer[tempIn], temp));
			return;
		}
		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			DualNumberContainer.push_back(DualNumberContainer[tempIn].multiply(DualNumberContainer[temp[0] - 1], DualNumberContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".randomize()") != std::wstring::npos) {
			DualNumberContainer[tempIn].randomize();
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			DualNumberContainer.push_back(DualNumberContainer[tempIn].subtract(DualNumberContainer[temp[0] - 1], DualNumberContainer[temp[1] - 1]));
			return;
		}


	}//END DualNUMBER FUNCTIONS======================================================================

	 //VECTORS
	 //=======
	 //constructor
	if (in.find(L"Vector(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(6);
		std::vector<double> temp = ParseNInputNumbers(tempstr);
		Vector t(temp);
		VectorContainer.push_back(t);
		return;
	}

	//class functions
	//===============
	if (in[0] == 'v' && isNumber(in[1]) == true) {    //reference to saved Vector
		double tempIn = ParseInputNumber(in.substr(1, 2));
		tempIn--;

		if (in.find(L".AndersonDarlingTest(") != std::wstring::npos) {
			double answer = (Vector(AndersonDarlingTest(VectorContainer[tempIn]))).vec[0];
			output.append(in + L" = " + to_stringPrecision(answer) + L"\r\n\r\n");
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			VectorContainer.push_back(VectorContainer[tempIn]);
			return;
		}
		if (in.find(L".convertCartesianToPolar()") != std::wstring::npos) {
			VectorContainer[tempIn].convertCartesianToPolar();
			return;
		}
		if (in.find(L".convertCylindricalToRectangular()") != std::wstring::npos) {
			VectorContainer[tempIn].convertCylindricalToRectangular();
			return;
		}
		if (in.find(L".convertCylindricalToRectangular()") != std::wstring::npos) {
			VectorContainer[tempIn].convertCylindricalToRectangular();
			return;
		}
		if (in.find(L".convertPolarToCartesian()") != std::wstring::npos) {
			VectorContainer[tempIn].convertPolarToCartesian();
			return;
		}
		if (in.find(L".convertRectangularToCylindrical()") != std::wstring::npos) {
			VectorContainer[tempIn].convertRectangularToCylindrical();
			return;
		}
		if (in.find(L".convertRectangularToSpherical()") != std::wstring::npos) {
			VectorContainer[tempIn].convertRectangularToSpherical();
			return;
		}
		if (in.find(L".convertSphericalToRectangular()") != std::wstring::npos) {
			VectorContainer[tempIn].convertSphericalToRectangular();
			return;
		}
		/*if (in.find(L".covariance(") != std::wstring::npos) {
		std::vector<double> temp1 = ParseNInputNumbers(in.substr(14));
		double dp = VectorContainer[tempIn].covariance(VectorContainer[temp1[0] - 1], VectorContainer[temp1[1] - 1]);
		output.append(L"v");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".covariance(");
		output.append(to_stringPrecision(temp1[0]));
		output.append(L",");
		output.append(to_stringPrecision(temp1[1]));
		output.append(L") = ");
		output.append(to_stringPrecision(dp));
		output.append(L"\r\n\r\n");
		return;
		}*/
		if (in.find(L".correlationCoefficient(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(25));
			double dp = correlationCoefficient(VectorContainer[temp1[0] - 1], VectorContainer[temp1[1] - 1]);
			output.append(L"R(v");
			output.append(to_stringPrecision(temp1[0]));
			output.append(L", v");
			output.append(to_stringPrecision(temp1[1]));
			output.append(L") = ");
			output.append(to_stringPrecision(dp));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".coefficientOfDetermination(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(30));
			double dp = coefficientOfDetermination(VectorContainer[temp1[0] - 1], VectorContainer[temp1[1] - 1]);
			output.append(L"R^2(v");
			output.append(to_stringPrecision(temp1[0]));
			output.append(L", v");
			output.append(to_stringPrecision(temp1[1]));
			output.append(L") = ");
			output.append(to_stringPrecision(dp));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".createRotationQuaternion()") != std::wstring::npos) {
			QuaternionContainer.push_back(VectorContainer[tempIn].createRotationQuaternion());
			return;
		}
		if (in.find(L".createRotationQuaternion(") != std::wstring::npos) {
			int v = ParseInputNumber(in.substr(28));
			QuaternionContainer.push_back(VectorContainer[tempIn].createRotationQuaternion(VectorContainer[v - 1].vec));
			return;
		}
		if (in.find(L".crossProduct(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(15));
			VectorContainer.push_back(VectorContainer[tempIn].crossProduct(VectorContainer[temp1[0] - 1], VectorContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".delete(") != std::wstring::npos) {
			VectorContainer[tempIn].vec.clear();
			VectorContainer.erase(VectorContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".directSum(") != std::wstring::npos) {
			std::vector<int> temp1 = STLDoubleToInt(ParseNInputNumbers(in.substr(13)));
			VectorContainer.push_back(VectorContainer[tempIn].directSum(VectorContainer[temp1[0] - 1], VectorContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".discreteFourierTransform(") != std::wstring::npos) {
			std::vector<int> temp1 = STLDoubleToInt(ParseNInputNumbers(in.substr(13)));
			VectorContainer.push_back(VectorContainer[tempIn].directSum(VectorContainer[temp1[0] - 1], VectorContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".dotProduct(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(14));
			double dp = VectorContainer[tempIn].dotProduct(VectorContainer[temp1[0] - 1], VectorContainer[temp1[1] - 1]);
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".dotProduct(");
			output.append(to_stringPrecision(temp1[0]));
			output.append(L",");
			output.append(to_stringPrecision(temp1[1]));
			output.append(L") = ");
			output.append(to_stringPrecision(dp));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".geometricMean()") != std::wstring::npos) {
			double mean = VectorContainer[tempIn].geometricMean();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".geometricMean() = ");
			output.append(to_stringPrecision(mean));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".getAngle(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(12));
			double angle = VectorContainer[tempIn].getAngle(temp1[0], temp1[1]);
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".mean() = ");
			output.append(to_stringPrecision(angle));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".mean()") != std::wstring::npos) {
			double mean = VectorContainer[tempIn].mean();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".mean() = ");
			output.append(to_stringPrecision(mean));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".norm()") != std::wstring::npos) {
			double norm = VectorContainer[tempIn].norm();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".norm() = ");
			output.append(to_stringPrecision(norm));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".pNorm(") != std::wstring::npos) {
			double p = ParseInputNumber(in.substr(9));
			double pNorm = VectorContainer[tempIn].pNorm(p);
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".pNorm() = ");
			output.append(to_stringPrecision(pNorm));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".populationStandardDeviation()") != std::wstring::npos) {
			double stddev = VectorContainer[tempIn].populationStandardDeviation();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".populationStandardDeviation() = ");
			output.append(to_stringPrecision(stddev));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".populationVariance()") != std::wstring::npos) {
			double var = VectorContainer[tempIn].populationVariance();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".populationVariance() = ");
			output.append(to_stringPrecision(var));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".pvalues()") != std::wstring::npos) {
			Vector v = pvalues(VectorContainer[tempIn]);
			if (v.size() > 0) { VectorContainer.push_back(v); }
			return;
		}
		if (in.find(L".randomize()") != std::wstring::npos) {
			VectorContainer[tempIn].randomize();
			return;
		}
		if (in.find(L".randomizeBoolean()") != std::wstring::npos) {
			VectorContainer[tempIn].randomizeBoolean();
			return;
		}
		if (in.find(L".randomizeInteger()") != std::wstring::npos) {
			VectorContainer[tempIn].randomizeInteger();
			return;
		}
		if (in.find(L".sampleStandardDeviation()") != std::wstring::npos) {
			double stddev = VectorContainer[tempIn].sampleStandardDeviation();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sampleStandardDeviation() = ");
			output.append(to_stringPrecision(stddev));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".sampleVariance()") != std::wstring::npos) {
			double var = VectorContainer[tempIn].sampleVariance();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sampleVariance() = ");
			output.append(to_stringPrecision(var));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".scalarTripleProduct(") != std::wstring::npos) {
			std::vector<int> temp = STLDoubleToInt(ParseNInputNumbers(in.substr(23)));
			double spd = VectorContainer[tempIn].scalarTripleProduct(
				VectorContainer[temp[0]],
				VectorContainer[temp[1]],
				VectorContainer[temp[2]]
			);
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".scalarTripleProduct(v");
			output.append(std::to_wstring(temp[0]));
			output.append(L", v");
			output.append(std::to_wstring(temp[1]));
			output.append(L", v");
			output.append(std::to_wstring(temp[2]));
			output.append(L") = ");
			output.append(to_stringPrecision(spd));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".scalingMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(scalingMatrix(VectorContainer[tempIn].vec));
			return;
		}
		if (in.find(L".set(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(7));
			if (temp1.size() >= 2) { VectorContainer[tempIn].vec[temp1[0]] = temp1[1]; }
			return;
		}
		if (in.find(L".size()") != std::wstring::npos) {
			double sz = VectorContainer[tempIn].size();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".size() = ");
			output.append(to_stringPrecision(sz));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".standardErrorOfMean()") != std::wstring::npos) {
			double sz = VectorContainer[tempIn].standardErrorOfMean();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".size() = ");
			output.append(to_stringPrecision(sz));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".sum()") != std::wstring::npos) {
			double sum = VectorContainer[tempIn].sum();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sum() = ");
			output.append(to_stringPrecision(sum));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".sumSquares()") != std::wstring::npos) {
			double sum = VectorContainer[tempIn].sumSquares();
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sumSquares() = ");
			output.append(to_stringPrecision(sum));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".toMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(VectorContainer[tempIn].toMatrix());
			return;
		}
		if (in.find(L".toMatrix(") != std::wstring::npos) {
			std::vector<double> vec = ParseNInputNumbers(in.substr(12));
			MatrixContainer.push_back(VectorContainer[tempIn].toMatrix(vec[0], vec[1]));
			return;
		}
		if (in.find(L".toQuaternion()") != std::wstring::npos) {
			QuaternionContainer.push_back(VectorContainer[tempIn].toQuaternion());
			return;
		}
		if (in.find(L".tensorProduct(") != std::wstring::npos) {
			std::vector<int> temp1 = STLDoubleToInt(ParseNInputNumbers(in.substr(17)));
			VectorContainer.push_back(VectorContainer[tempIn].tensorProduct(VectorContainer[temp1[0] - 1], VectorContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".translationMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(translationMatrix(VectorContainer[tempIn].vec));
			return;
		}
		if (in.find(L".vectorScalarProjection(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(26));
			double proj = VectorContainer[tempIn].vectorScalarProjection(temp1[0], temp1[1]);
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".vectorScalarProjection(");
			output.append(to_stringPrecision(temp1[0]));
			output.append(L",");
			output.append(to_stringPrecision(temp1[1]));
			output.append(L")= ");
			output.append(to_stringPrecision(proj));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".vectorTripleProduct(") != std::wstring::npos) {
			std::vector<int> temp = STLDoubleToInt(ParseNInputNumbers(in.substr(23)));
			Vector v = VectorContainer[tempIn].vectorTripleProduct(
				VectorContainer[temp[0]],
				VectorContainer[temp[1]],
				VectorContainer[temp[2]]
			);
			VectorContainer.push_back(v);
			output.append(L"v");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".scalarTripleProduct(v");
			output.append(std::to_wstring(temp[0]));
			output.append(L", v");
			output.append(std::to_wstring(temp[1]));
			output.append(L", v");
			output.append(std::to_wstring(temp[2]));
			output.append(L") = ");
			output.append(v.toString());
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".unitNormal(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(14));
			VectorContainer.push_back(VectorContainer[tempIn].unitNormal(temp1[0], temp1[1]));
			return;
		}
		if (in.find(L".zscores()") != std::wstring::npos) {
			Vector v = zscores(VectorContainer[tempIn]);
			if (v.size() > 0) { VectorContainer.push_back(v); }
			return;
		}
	}


	//MATRICES
	//========
	//constructor
	if (in.find(L"Matrix(") != std::wstring::npos && in.find(L"ComplexMatrix(") == std::wstring::npos && in[0] == 'M' && in[1] == 'a' && isNumber(in[7]) == true) {
		std::wstring tempstr = in.substr(7);
		if (tempstr.find(L",") != std::wstring::npos) {   //if there are both a row and column
			std::vector<double> temp = ParseNInputNumbers(tempstr);
			Matrix t(temp[0], temp[1]);
			t.identity();
			MatrixContainer.push_back(t);
		}
		else {
			double temp = ParseInputNumber(tempstr);
			Matrix t(temp, temp);
			t.identity();
			MatrixContainer.push_back(t);
		}
		return;
	}
	//Assorted class functions of matrices
	//====================================
	if (in[0] == 'M' && isNumber(in[1]) == true) {    //reference to saved Matrix
		double tempIn = ParseInputNumber(in.substr(1, 2));
		tempIn--;

		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(7));
			MatrixContainer.push_back(MatrixContainer[tempIn].add(MatrixContainer[temp1[0] - 1], MatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".addScalar(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(13));
			MatrixContainer.push_back(MatrixContainer[tempIn].add(MatrixContainer[temp1[0] - 1], temp1[1]));
			return;
		}
		if (in.find(L".adjustedRsquared()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".adjustedRsquared() = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].adjustedRsquared()));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".AndersonDarlingTest(") != std::wstring::npos) {
			VectorContainer.push_back(Vector(AndersonDarlingTest(MatrixContainer[tempIn])));
			return;
		}
		if (in.find(L".BooleanProduct(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(18));
			MatrixContainer.push_back(MatrixContainer[tempIn].BooleanProduct(MatrixContainer[temp1[0] - 1], MatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".characteristicMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].characteristicMatrix());
			return;
		}
		if (in.find(L".characteristicPolynomial()") != std::wstring::npos) {
			PolynomialContainer.push_back(MatrixContainer[tempIn].characteristicPolynomial());
			return;
		}
		if (in.find(L".chiSquareTest()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".ChiSquareTest() = ");
			output.append(to_stringPrecision(ChiSquareTest(MatrixContainer[tempIn])));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".chiSquareDegreesFreedom()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".chiSquareDegreesFreedom() = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].ChiSquareDegreesFreedom()));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".chiSquareTestStatistic()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".ChiSquareTestStatistic() = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].ChiSquareTestStatistic()));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".Cholesky()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].Cholesky());
			return;
		}
		if (in.find(L".column(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(10));
			VectorContainer.push_back(Vector(MatrixContainer[tempIn].column(temp1)));
			return;
		}
		if (in.find(L".columns()") != std::wstring::npos) {
			output.append(in + L" = " + std::to_wstring(MatrixContainer[tempIn].columns) + L"\r\n\r\n");
			return;
		}
		if (in.find(L".columnCovariance(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(21));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".columnCovariance(");
			output.append(to_stringPrecision(temp1[0]));
			output.append(L",");
			output.append(to_stringPrecision(temp1[1]));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].columnCovariance(temp1[0], temp1[1])));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".columnDeviation(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(20));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".columnDeviation(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].columnDeviation(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".columnMean(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(13));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".columnMean(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].columnMean(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".columnNorm(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(10));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".columnNorm(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].columnNorm(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".columnSumSquares(") != std::wstring::npos) {
			double temp1 = ParseInputNumber(in.substr(17));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".columnSumSquares(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].columnSumSquares(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".columnVariance(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(19));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".columnVariance(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].columnVariance(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn]);
			return;
		}
		if (in.find(L".correlationMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(correlationMatrix(MatrixContainer[tempIn]));
			return;
		}
		if (in.find(L".covarianceMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(covarianceMatrix(MatrixContainer[tempIn]));
			return;
		}
		if (in.find(L".DAgostinoTest(") != std::wstring::npos) {
			double critVal = ParseInputNumber(in.substr(17));
			VectorContainer.push_back(Vector(DAgostinoTest(MatrixContainer[tempIn], critVal)));
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			if (tempIn >= 0 && tempIn < MatrixContainer.size()) {
				MatrixContainer[tempIn].element.clear();
				MatrixContainer.erase(MatrixContainer.begin() + tempIn);
			}
			return;
		}
		if (in.find(L".det()") != std::wstring::npos || in.find(L".determinant()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".det() = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].det()));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".diagonal(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(12));
			VectorContainer.push_back(Vector(MatrixContainer[tempIn].diagonal()));
			return;
		}
		if (in.find(L".directSum(") != std::wstring::npos) {
			std::vector<int> temp1 = STLDoubleToInt(ParseNInputNumbers(in.substr(13)));
			MatrixContainer.push_back(MatrixContainer[tempIn].directSum(MatrixContainer[temp1[0] - 1], MatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".eigenvaluesRealAndComplex()") != std::wstring::npos) {
			std::vector<ComplexNumber> eigs = MatrixContainer[tempIn].eigenvaluesRealAndComplexExact();
			if (eigs.size() > 0) {
				output.append(L"M");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".eigenvalues() = \n");
				for (int i = 0; i < eigs.size(); ++i) {
					output.append(eigs[i].toString());
					if (i < eigs.size() - 1) { output.append(L", "); }
					ComplexNumberContainer.push_back(eigs[i]);
				}
				output.append(L"\r\n\r\n");
			}
			else { output.append(L"Eigenvalues of the matrix could not be found.\n"); }
			return;
		}
		if (in.find(L".eigenvaluesRealExact()") != std::wstring::npos) {
			std::vector<double> eigs = MatrixContainer[tempIn].eigenvaluesRealExact();
			if (eigs.size() > 0) {
				output.append(L"M");
				output.append(to_stringPrecision(tempIn + 1));
				output.append(L".eigenvalues() = \nreal: ");
				for (int i = 0; i < eigs.size(); ++i) {
					output.append(to_stringPrecision(eigs[i]));
					if (i < eigs.size() - 1) { output.append(L","); }
				}
				output.append(L"\r\n\r\n");
				VectorContainer.push_back(Vector(eigs));
			}
			else { output.append(L"real eigenvalues of matrix could not be found.\nEither they are all complex or \nthe matrix is too large or there was some error.\r\n\r\n"); }
			return;
		}
		if (in.find(L".eigenvectors()") != std::wstring::npos) {
			Matrix egv = MatrixContainer[tempIn].eigenvectors();
			if (egv.size() > 0) { MatrixContainer.push_back(egv); }
			else { output.append(L"Could not find eigenvectors of this matrix."); }
			return;
		}
		if (in.find(L".eigenvaluesNumerical()") != std::wstring::npos) {
			Matrix cpy = MatrixContainer[tempIn];
			std::vector<double> evm = MatrixContainer[tempIn].eigenvaluesNumerical();
			std::vector<double> eigs;
			if (evm.size() > 0) {
				Polynomial p = cpy.characteristicPolynomial();
				for (int i = 0; i < evm.size(); ++i) {
					double rt = evm[i];
					rt = findRoot(p, rt - 1, rt + 1);
					if (std::abs(p.HornerEvaluate(rt) < 0.1) && rt != 0) { eigs.push_back(rt); }
				}
				if (eigs.size() > 0) {
					output.append(L"M");
					output.append(to_stringPrecision(tempIn + 1));
					output.append(L".eigenvalues() = \nreal: ");
					for (int i = 0; i < eigs.size(); ++i) {
						output.append(to_stringPrecision(eigs[i]));
						if (i < eigs.size() - 1) { output.append(L","); }
					}
					output.append(L"\r\n\r\n");
					VectorContainer.push_back(Vector(eigs));
				}
			}
			return;
		}
		if (in.find(L".eigenvalueMatrixRecursive(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(30));
			if (tempIn < 0 || tempIn >= MatrixContainer.size()) { return; }
			MatrixContainer.push_back(MatrixContainer[tempIn].eigenvalueMatrix(MatrixContainer[temp1[0] - 1], temp1[1]));
			return;
		}
		if (in.find(L".exponentiate(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(16));
			MatrixContainer.push_back(MatrixContainer[tempIn].exponent(MatrixContainer[temp1[0] - 1], temp1[1]));
			return;
		}
		if (in.find(L".FrobeniusNorm()") != std::wstring::npos) {
			double FrobeniusNorm = MatrixContainer[tempIn].FrobeniusNorm();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".FrobeniusNorm() = ");
			output.append(to_stringPrecision(FrobeniusNorm));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".FTestStatistic()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".FTestStatistic() = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].FTestStatistic()));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".FTest()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".FTest() = ");
			output.append(to_stringPrecision(FDistributionApproximation(MatrixContainer[tempIn].FTestStatistic(), MatrixContainer[tempIn].columns, MatrixContainer[tempIn].size() - MatrixContainer[tempIn].columns)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".geometricMean()") != std::wstring::npos) {
			double mean = MatrixContainer[tempIn].geometricMean();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".geometricMean() = ");
			output.append(to_stringPrecision(mean));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".GramSchmidt()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].GramSchmidt());
			return;
		}
		if (in.find(L".GramMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].GramMatrix());
			return;
		}
		if (in.find(L".GivensRotationMatrix(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(22));
			MatrixContainer.push_back(MatrixContainer[tempIn].GivensRotationMatrix(temp1[0], temp1[1], temp1[2], temp1[3]));
			return;
		}
		if (in.find(L".get(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(6));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".get(");
			output.append(to_stringPrecision(temp1[0]));
			output.append(L",");
			output.append(to_stringPrecision(temp1[1]));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].get(temp1[0], temp1[1])));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".HouseholderQR()") != std::wstring::npos) {
			std::vector<Matrix> qr = MatrixContainer[tempIn].HouseholderQR();
			MatrixContainer.push_back(qr[0]);
			MatrixContainer.push_back(qr[1]);
			qr.clear();
			return;
		}
		if (in.find(L".identity()") != std::wstring::npos) {
			MatrixContainer[tempIn].identity();
			return;
		}
		if (in.find(L".inverse()") != std::wstring::npos) {
			Matrix inv = MatrixContainer[tempIn].inverse();
			if (inv.size()>0) { MatrixContainer.push_back(inv); }
			return;
		}
		if (in.find(L".inverseByQR()") != std::wstring::npos) {
			Matrix inv = MatrixContainer[tempIn].inverseByQR();
			if (inv.size() > 0) { MatrixContainer.push_back(inv); }
			return;
		}
		if (in.find(L".isBoolean()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isBoolean() = ");
			bool tf = MatrixContainer[tempIn].isBoolean();
			std::wstring trf = L"";
			if (tf == false) { trf = L"false"; }
			if (tf == true) { trf = L"true"; }
			output.append(trf);
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".isDiagonalized()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isBoisDiagonalizedolean() = ");
			bool tf = MatrixContainer[tempIn].isDiagonalized();
			std::wstring trf = L"";
			if (tf == false) { trf = L"false"; }
			if (tf == true) { trf = L"true"; }
			output.append(trf);
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".isIntegerMatrix()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isIntegerMatrix() = ");
			bool tf = MatrixContainer[tempIn].isIntegerMatrix();
			std::wstring trf = L"";
			if (tf == false) { trf = L"false"; }
			if (tf == true) { trf = L"true"; }
			output.append(trf);
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".isLinearlyIndependent()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isLinearlyIndependent() = ");
			bool tf = MatrixContainer[tempIn].isLinearlyIndependent();
			std::wstring trf = L"";
			if (tf == false) { trf = L"false"; }
			if (tf == true) { trf = L"true"; }
			output.append(trf);
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".isOrthonormal()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isOrthonormal() = ");
			bool tf = MatrixContainer[tempIn].isOrthonormal();
			std::wstring trf = L"";
			if (tf == false) { trf = L"false"; }
			if (tf == true) { trf = L"true"; }
			output.append(trf);
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".isSymmetric()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isSymmetric() = ");
			bool tf = MatrixContainer[tempIn].isSymmetric();
			std::wstring trf = L"";
			if (tf == false) { trf = L"false"; }
			if (tf == true) { trf = L"true"; }
			output.append(trf);
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".JacobiRotationMatrix(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(22));
			MatrixContainer.push_back(MatrixContainer[tempIn].JacobiRotationMatrix(temp1[0], temp1[1], temp1[2], temp1[3]));
			return;
		}
		if (in.find(L".JacobiTransformation()") != std::wstring::npos) {
			std::vector<Matrix> jt = MatrixContainer[tempIn].JacobiTransformation();
			if (jt.size() != NULL) {
				MatrixContainer.push_back(jt[0]);
				MatrixContainer.push_back(jt[1]);
				jt.clear();
			}
			return;
		}
		if (in.find(L".JarqueBeraKurtosis()") != std::wstring::npos) {
			VectorContainer.push_back(Vector(JarqueBeraKurtosis(MatrixContainer[tempIn])));
			return;
		}
		if (in.find(L".JarqueBeraSkewness()") != std::wstring::npos) {
			VectorContainer.push_back(Vector(JarqueBeraSkewness(MatrixContainer[tempIn])));
			return;
		}
		if (in.find(L".JarqueBeraTest(") != std::wstring::npos) {
			double critVal = ParseInputNumber(in.substr(18));
			VectorContainer.push_back(Vector(JarqueBeraTest(MatrixContainer[tempIn], critVal)));
			return;
		}
		if (in.find(L".L()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].L());
			return;
		}
		if (in.find(L".LDL()") != std::wstring::npos) {
			std::vector<Matrix> ldl = MatrixContainer[tempIn].LDL();
			MatrixContainer.push_back(ldl[0]);
			MatrixContainer.push_back(ldl[1]);
			ldl.clear();
			return;
		}
		if (in.find(L".leastSquares()") != std::wstring::npos) {
			FunctionContainer.push_back(MatrixContainer[tempIn].leastSquares());
			return;
		}
		if (in.find(L".leastSquaresMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].leastSquaresMatrix());
			return;
		}
		if (in.find(L".leastSquares(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(15));
			Matrix M = MatrixContainer[tempIn].leastSquares(MatrixContainer[temp1[0] - 1], MatrixContainer[temp1[1] - 1]);
			if (M.size() > 0) { MatrixContainer.push_back(M); }
			return;
		}
		if (in.find(L".lowerHessenbergForm()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].lowerHessenbergForm());
			return;
		}
		if (in.find(L".LU()") != std::wstring::npos) {
			std::vector<Matrix> lu = MatrixContainer[tempIn].LU();
			MatrixContainer.push_back(lu[0]);
			MatrixContainer.push_back(lu[1]);
			lu.clear();
			return;
		}
		if (in.find(L".mean()") != std::wstring::npos) {
			double mean = MatrixContainer[tempIn].mean();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".mean() = ");
			output.append(to_stringPrecision(mean));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".meanVector()") != std::wstring::npos) {
			Vector v(MatrixContainer[tempIn].meanVector());
			if (v.size()>0) { VectorContainer.push_back(v); }
			return;
		}
		if (in.find(L".modulus()") != std::wstring::npos) {
			int mod = MatrixContainer[tempIn].modulus;
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".modulus() = ");
			output.append(to_stringPrecision(mod));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(12));
			MatrixContainer.push_back(MatrixContainer[tempIn].multiply(MatrixContainer[temp1[0] - 1], MatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".multiplyScalar(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(18));
			MatrixContainer.push_back(MatrixContainer[tempIn].multiply(MatrixContainer[temp1[0] - 1], temp1[1]));
			return;
		}
		if (in.find(L".norm()") != std::wstring::npos) {
			double norm = MatrixContainer[tempIn].norm();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".norm() = ");
			output.append(to_stringPrecision(norm));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".normalize()") != std::wstring::npos) {
			MatrixContainer[tempIn].normalize();
			return;
		}
		if (in.find(L".pNorm(") != std::wstring::npos) {
			double p = ParseInputNumber(in.substr(9));
			double pNorm = MatrixContainer[tempIn].pNorm(p);
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".pNorm() = ");
			output.append(to_stringPrecision(pNorm));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".populationCovarianceMatrix()") != std::wstring::npos) {
			Matrix cov = MatrixContainer[tempIn].populationCovarianceMatrix();
			if (cov.size()>0) { MatrixContainer.push_back(cov); }
			return;
		}
		if (in.find(L".populationStandardDeviation()") != std::wstring::npos) {
			double pop = MatrixContainer[tempIn].populationStandardDeviation();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".populationStandardDeviation() = ");
			output.append(to_stringPrecision(pop));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".power(") != std::wstring::npos) {
			double n = ParseInputNumber(in.substr(9));
			MatrixContainer[tempIn].power(n);
			return;
		}
		if (in.find(L".pseudoinverse()") != std::wstring::npos) {
			Matrix pinv = MatrixContainer[tempIn].pseudoinverse();
			if (pinv.size() > 0) { MatrixContainer.push_back(pinv); pinv.element.clear(); }
			return;
		}
		if (in.find(L".QR()") != std::wstring::npos) {
			std::vector<Matrix> qr = MatrixContainer[tempIn].QR();
			MatrixContainer.push_back(qr[0]);
			MatrixContainer.push_back(qr[1]);
			qr.clear();
			return;
		}
		if (in.find(L".randomize()") != std::wstring::npos) {
			MatrixContainer[tempIn].randomize();
			return;
		}
		if (in.find(L".randomizeBoolean()") != std::wstring::npos) {
			MatrixContainer[tempIn].randomizeBoolean();
			return;
		}
		if (in.find(L".randomizeInteger()") != std::wstring::npos) {
			MatrixContainer[tempIn].randomizeInteger();
			return;
		}
		if (in.find(L".randomizeSymmetric()") != std::wstring::npos) {
			MatrixContainer[tempIn].randomizeSymmetric();
			return;
		}
		if (in.find(L".reducedRowEchelonForm()") != std::wstring::npos) {
			MatrixContainer[tempIn].reducedRowEchelonForm();
			return;
		}
		if (in.find(L".removeColumn(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(16));
			MatrixContainer[tempIn].removeColumn(temp1);
			return;
		}
		if (in.find(L".removeRow(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(13));
			MatrixContainer[tempIn].removeRow(temp1);
			return;
		}
		if (in.find(L".residuals()") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(13));
			VectorContainer.push_back(Vector(MatrixContainer[tempIn].residuals()));
			return;
		}
		if (in.find(L".row(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(7));
			VectorContainer.push_back(Vector(MatrixContainer[tempIn].row(temp1)));
			return;
		}
		if (in.find(L".rows()") != std::wstring::npos) {
			output.append(in + L" = " + std::to_wstring(MatrixContainer[tempIn].rows) + L"\r\n\r\n");
			return;
		}
		if (in.find(L".rowDeviation(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(17));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".rowDeviation(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].rowDeviation(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".rowMean(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(10));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".rowMean(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].rowMean(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".rowNorm(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(10));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".rowNorm(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].rowNorm(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".rowSumSquares(") != std::wstring::npos) {
			double temp1 = ParseInputNumber(in.substr(16));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".rowSumSquares(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].rowSumSquares(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".rowVariance(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(16));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".rowVariance(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].rowVariance(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".Rsquared()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".Rsquared() = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].Rsquared()));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".sampleCovarianceMatrix()") != std::wstring::npos) {
			Matrix cov = MatrixContainer[tempIn].sampleCovarianceMatrix();
			if (cov.size()>0) { MatrixContainer.push_back(cov); }
			return;
		}
		if (in.find(L".sampleStandardDeviation()") != std::wstring::npos) {
			double std = MatrixContainer[tempIn].sampleStandardDeviation();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sampleStandardDeviation() = ");
			output.append(to_stringPrecision(std));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".set(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(7));
			if (temp1.size() >= 3) { MatrixContainer[tempIn].set(temp1[0], temp1[1], temp1[2]); }
			return;
		}
		if (in.find(L".setModulus(") != std::wstring::npos) {
			int temp = ParseInputNumber(in.substr(14));
			MatrixContainer[tempIn].setModulus(temp);
			output.append(L"Matrix M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L"'s modulus is set to ");
			output.append(to_stringPrecision(temp));
			output.append(L"\r\n\r\n");

			return;
		}
		if (in.find(L".size()") != std::wstring::npos) {
			double sz = MatrixContainer[tempIn].size();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".size() = ");
			output.append(to_stringPrecision(sz));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".submatrix(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(13));
			MatrixContainer.push_back(MatrixContainer[tempIn].submatrix(temp1[0], temp1[1], temp1[2], temp1[3]));
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(12));
			MatrixContainer.push_back(MatrixContainer[tempIn].subtract(MatrixContainer[temp1[0] - 1], MatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".subtractScalar(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(18));
			MatrixContainer.push_back(MatrixContainer[tempIn].subtract(MatrixContainer[temp1[0] - 1], temp1[1]));
			return;
		}
		if (in.find(L".sum()") != std::wstring::npos) {
			double sum = MatrixContainer[tempIn].sum();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sum() = ");
			output.append(to_stringPrecision(sum));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".sumColumn(") != std::wstring::npos) {
			double temp1 = ParseInputNumber(in.substr(13));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sumColumn(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].sumColumn(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".sumRow(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(10));
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sumRow(");
			output.append(to_stringPrecision(temp1));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].sumRow(temp1)));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".sumSquares()") != std::wstring::npos) {
			double sumSquares = MatrixContainer[tempIn].sumSquares();
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".sumSquares() = ");
			output.append(to_stringPrecision(sumSquares));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".swapColumns(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(15));
			MatrixContainer[tempIn].swapColumns(temp1[0], temp1[1]);
			return;
		}
		if (in.find(L".swapRows(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(12));
			MatrixContainer[tempIn].swapRows(temp1[0], temp1[1]);
			return;
		}
		if (in.find(L".tensorProduct(") != std::wstring::npos) {
			std::vector<int> temp1 = STLDoubleToInt(ParseNInputNumbers(in.substr(17)));
			MatrixContainer.push_back(MatrixContainer[tempIn].tensorProduct(MatrixContainer[temp1[0] - 1], MatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".thinQR()") != std::wstring::npos) {
			std::vector<Matrix> qr = MatrixContainer[tempIn].thinQR();
			MatrixContainer.push_back(qr[0]);
			MatrixContainer.push_back(qr[1]);
			qr.clear();
			return;
		}
		if (in.find(L".thinHouseholderQR()") != std::wstring::npos) {
			std::vector<Matrix> qr = MatrixContainer[tempIn].thinHouseholderQR();
			MatrixContainer.push_back(qr[0]);
			MatrixContainer.push_back(qr[1]);
			qr.clear();
			return;
		}
		if (in.find(L".transpose()") != std::wstring::npos) {
			MatrixContainer[tempIn] = MatrixContainer[tempIn].transpose();
			return;
		}
		if (in.find(L".trim()") != std::wstring::npos) {
			MatrixContainer[tempIn].trim();
			return;
		}
		if (in.find(L".U()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].U());
			return;
		}
		if (in.find(L".upperHessenbergForm()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].upperHessenbergForm());
			return;
		}
		if (in.find(L".VandermondeMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(MatrixContainer[tempIn].Vandermonde());
			return;
		}
	}

	//Tensors
	//=======
	//constructor
	if (in.find(L"Tensor(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(7);
		std::vector<double> temp = ParseNInputNumbers(tempstr);
		Tensor t(temp);
		TensorContainer.push_back(t);
		t.element.clear();
		temp.clear();
		return;
	}

	//assorted class functions of Tensor
	//=========================================
	if (in[0] == 'T' && isNumber(in[1]) == true) {    //reference to saved Tensor
		double tempIn = ParseInputNumber(in.substr(1, 2));
		tempIn--;
		//------------------------------------

		if (in.find(L".copy()") != std::wstring::npos) {
			TensorContainer.push_back(TensorContainer[tempIn]);
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			TensorContainer.erase(TensorContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".get(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(6));
			output.append(L"T");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".get(");
			for (int i = 0; i < temp1.size() - 1; ++i) {
				output.append(to_stringPrecision(temp1[i]));
				output.append(L",");
			}
			output.append(to_stringPrecision(temp1[temp1.size() - 1]));
			output.append(to_stringPrecision(temp1[1]));
			output.append(L") = ");
			output.append(to_stringPrecision(MatrixContainer[tempIn].get(temp1[0], temp1[1])));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".indices()") != std::wstring::npos) {
			output.append(L"T");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".indices() = ");
			for (int i = 0; i < TensorContainer[tempIn].index.size() - 1; ++i) {
				output.append(std::to_wstring(TensorContainer[tempIn].index[i]));
				output.append(L", ");
			}
			output.append(std::to_wstring(TensorContainer[tempIn].index[TensorContainer[tempIn].index.size() - 1]));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".isCovariant(") != std::wstring::npos) {
			int v = ParseInputNumber(in.substr(15));
			if (v == 1) { TensorContainer[tempIn].isCovariant = true; }
			if (v == 0) { TensorContainer[tempIn].isCovariant = false; }
			return;
		}
		if (in.find(L".randomize()") != std::wstring::npos) {
			TensorContainer[tempIn].randomize();
			return;
		}
		if (in.find(L".setDimensions(") != std::wstring::npos) {
			std::vector<int> v = STLDoubleToInt(ParseNInputNumbers(in.substr(17)));
			TensorContainer[tempIn].index = v;
			return;
		}


	}//END TENSOR FUNCTIONS======================================================================


	 //InfiniteSeries
	 //==============
	 //constructor
	if (in.find(L"InfiniteSeries(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(16);
		Function temp(ParseFunction(tempstr));
		InfiniteSeries s(temp);
		InfiniteSeriesContainer.push_back(s);
		temp.clear();
		return;
	}

	//assorted class functions of InfiniteSeries
	//=========================================
	if (in[0] == 'I' && in[1] == 'S' && isNumber(in[2]) == true) {    //reference to saved InfiniteSeries
		double tempIn = ParseInputNumber(in.substr(1, 3));
		tempIn--;
		//------------------------------------
		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(7));
			InfiniteSeriesContainer.push_back(InfiniteSeriesContainer[tempIn].subtract(InfiniteSeriesContainer[temp[0] - 1], InfiniteSeriesContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			InfiniteSeriesContainer.push_back(InfiniteSeriesContainer[tempIn]);
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			InfiniteSeriesContainer.erase(InfiniteSeriesContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".evaluate(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(13));
			output.append(L"IS");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".evaluate(");
			output.append(in.substr(13));
			output.append(L" = ");
			output.append(to_stringPrecision(InfiniteSeriesContainer[tempIn].evaluate(temp[0], temp[1])));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".exponentiate(") != std::wstring::npos) {
			int temp = ParseInputNumber(in.substr(16));
			InfiniteSeriesContainer.push_back(InfiniteSeriesContainer[tempIn].exponent(InfiniteSeriesContainer[tempIn], temp));
			return;
		}
		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			InfiniteSeriesContainer.push_back(InfiniteSeriesContainer[tempIn].multiply(InfiniteSeriesContainer[temp[0] - 1], InfiniteSeriesContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			InfiniteSeriesContainer.push_back(InfiniteSeriesContainer[tempIn].subtract(InfiniteSeriesContainer[temp[0] - 1], InfiniteSeriesContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".testConvergence(") != std::wstring::npos) {
			unsigned int temp = ParseInputNumber(in.substr(20));
			output.append(L"IS");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".testConvergence(");
			output.append(in.substr(20));
			output.append(L" = ");
			output.append(InfiniteSeriesContainer[tempIn].testConvergence(temp));
			output.append(L"\r\n\r\n");
			return;
		}


	}//END INFINITESERIES FUNCTIONS======================================================================

	 //InfiniteSeries
	 //==============
	 //constructor
	if (in.find(L"InfiniteSeries(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(16);
		Function temp(ParseFunction(tempstr));
		InfiniteSeries s(temp);
		InfiniteSeriesContainer.push_back(s);
		temp.clear();
		return;
	}

	//assorted class functions of ComplexInfiniteSeries
	//=========================================
	if (in[0] == 'C' && in[1] == 'I' && in[2] == 'S' && isNumber(in[3]) == true) {    //reference to saved ComplexInfiniteSeries
		double tempIn = ParseInputNumber(in.substr(1, 4));
		tempIn--;
		//------------------------------------
		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(7));
			ComplexInfiniteSeriesContainer.push_back(ComplexInfiniteSeriesContainer[tempIn].subtract(ComplexInfiniteSeriesContainer[temp[0] - 1], ComplexInfiniteSeriesContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			ComplexInfiniteSeriesContainer.push_back(ComplexInfiniteSeriesContainer[tempIn]);
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			ComplexInfiniteSeriesContainer.erase(ComplexInfiniteSeriesContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".evaluate(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(13));
			output.append(L"IS");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".evaluate(");
			output.append(in.substr(13));
			output.append(L") = ");
			output.append(to_stringPrecision(InfiniteSeriesContainer[tempIn].evaluate(temp[0], temp[1])));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".exponentiate(") != std::wstring::npos) {
			int temp = ParseInputNumber(in.substr(16));
			ComplexInfiniteSeriesContainer.push_back(ComplexInfiniteSeriesContainer[tempIn].exponent(ComplexInfiniteSeriesContainer[tempIn], temp));
			return;
		}
		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			ComplexInfiniteSeriesContainer.push_back(ComplexInfiniteSeriesContainer[tempIn].multiply(ComplexInfiniteSeriesContainer[temp[0] - 1], ComplexInfiniteSeriesContainer[temp[1] - 1]));
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(in.substr(12));
			ComplexInfiniteSeriesContainer.push_back(ComplexInfiniteSeriesContainer[tempIn].subtract(ComplexInfiniteSeriesContainer[temp[0] - 1], ComplexInfiniteSeriesContainer[temp[1] - 1]));
			return;
		}


	}//END COMPLEXINFINITESERIES FUNCTIONS======================================================================

	 //ODEs
	 //=====
	if (input.find(L"ODE") != std::wstring::npos) {	//get reference to saved ODEs
		int position = 0;
		for (int i = 0; i < input.length(); i++) {
			if (input[i] == '.') {
				position = i;
				i = input.length();
			}
		}
		std::wstring a = input.substr(0, position);
		double tempIn = ParseInputNumber(a);//get index of function array

											// Class functions for ODEs
											//========================
		if (in.find(L".copy()") != std::wstring::npos) {
			ODEContainer.push_back(ODEContainer[tempIn]);
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			if (tempIn >= 0 && tempIn < MatrixContainer.size()) {
				ODEContainer[tempIn].lhs = L"";
				ODEContainer[tempIn].rhs = L"";
				ODEContainer.erase(ODEContainer.begin() + tempIn);
			}
			return;
		}
		if (input.find(L".evaluate(") != std::wstring::npos) { //format is ODE1.evaluate(2xy + 23wz, 1.3, -2.023, 0, ...etc)
			std::wstring funct = L"";
			int pos = input.find(L"(");
			int posComma = input.find(L",");		//find first comma in string
			funct = input.substr(pos + 1, posComma - pos - 1);		//get the function that is used to replace 'y" in the 
			std::wstring evalNumb2 = input.substr(posComma + 1);
			std::vector<double> evalNumb = ParseNInputNumbers(evalNumb2);	//get arguments for input 
			std::vector<double> vecArgs = initZeros(vecArgs, evalNumb.size());
			for (int i = 0; i < evalNumb.size(); i++) {
				vecArgs[i] = evalNumb[i];		//vecArgs holds the final array of input numbers
			}

			std::wstring tempstr = ODEContainer[tempIn - 1].lhs;
			std::wstring tempstr2 = ODEContainer[tempIn - 1].rhs;
			for (int i = 0; i < tempstr.length();i++) {
				if (tempstr[i] == 'y') {
					tempstr.erase(tempstr.begin() + i, tempstr.begin() + i + 1);
					tempstr.insert(i, funct);
				}
			}
			for (int i = 0; i < tempstr2.length();i++) {
				if (tempstr2[i] == 'y') {
					tempstr2.erase(tempstr.begin() + i, tempstr2.begin() + i + 1);
					tempstr2.insert(i, funct);
				}
			}
			ODEContainer[tempIn - 1].lhs = tempstr;		//replace all the 'y' values in the original ODE equations
			ODEContainer[tempIn - 1].rhs = tempstr2;
			std::vector<std::wstring> f1 = ParseFunction(ODEContainer[tempIn - 1].lhs);
			std::vector<std::wstring> f2 = ParseFunction(ODEContainer[tempIn - 1].rhs);
			Function fp1(f1[0], f1[1]);
			Function fp2(f2[0], f2[1]);

			double result = fp1.evaluate(vecArgs);
			result += fp2.evaluate(vecArgs);
			output.append(L"ODE");
			output.append(to_stringPrecision(tempIn));
			output.append(L".evaluate(");
			output.append(input.substr(pos + 1));
			output.append(L" = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");

			vecArgs.clear();
			evalNumb.clear();
		}

	}//END ODE FUNCTIONS==================================================


	 //COMPLEX MATRICES
	 //================
	 //constructor
	if (in.find(L"ComplexMatrix(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(14);
		if (tempstr.find(L",") != std::wstring::npos) {   //if there are both a row and column
			std::vector<double> temp = ParseNInputNumbers(tempstr);
			ComplexMatrix t(temp[0], temp[1]);
			t.identity();
			ComplexMatrixContainer.push_back(t);
		}
		else {
			double temp = ParseInputNumber(tempstr);
			ComplexMatrix t(temp, temp);
			t.identity();
			ComplexMatrixContainer.push_back(t);
		}
		return;
	}
	//Assorted class functions of complex matrices
	//============================================
	if (in[0] == 'C' && in[1] == 'M' && isNumber(in[2]) == true) {    //reference to saved ComplexMatrix
		double tempIn = ParseInputNumber(in.substr(2, 2));
		tempIn--;

		/*	if (in.find(L".ChiSquareTestStatistic()") != std::wstring::npos) {
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".ChiSquareTestStatistic() = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].ChiSquareTestStatistic()));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".Cholesky()") != std::wstring::npos) {
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].Cholesky());
		return;
		}
		if (in.find(L".column(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(10));
		VectorContainer.push_back(Vector(ComplexMatrixContainer[tempIn].column(temp1), ComplexMatrixContainer[tempIn].rows));
		return;
		}
		if (in.find(L".columnCovariance(") != std::wstring::npos) {
		std::vector<double> temp1 = ParseNInputNumbers(in.substr(21));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".columnCovariance(");
		output.append(to_stringPrecision(temp1[0]));
		output.append(L",");
		output.append(to_stringPrecision(temp1[1]));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].columnCovariance(temp1[0], temp1[1])));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".columnDeviation(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(20));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".columnDeviation(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].columnDeviation(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".columnMean(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(13));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".columnMean(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].columnMean(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".columnNorm(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(10));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".columnNorm(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].columnNorm(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".columnSumSquares(") != std::wstring::npos) {
		double temp1 = ParseInputNumber(in.substr(17));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".columnSumSquares(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].columnSumSquares(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".columnVariance(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(19));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".columnVariance(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].columnVariance(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		*/	if (in.find(L".conjugateTranspose()") != std::wstring::npos) {
			ComplexMatrixContainer[tempIn] = ComplexMatrixContainer[tempIn].conjugateTranspose();
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn]);
			return;
		}
		if (in.find(L".delete(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(10));
			ComplexMatrixContainer[tempIn].element.clear();
			ComplexMatrixContainer.erase(ComplexMatrixContainer.begin() + temp1);
			return;
		}
		/*	if (in.find(L".det()") != std::wstring::npos) {
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".det() = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].det()));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".eigenvalues()") != std::wstring::npos) {
		std::vector<double> eigs = eigenvalues(ComplexMatrixContainer[tempIn]);
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".eigenvalues() = ");
		for (int i = 0; i < eigs.size(); ++i) {
		output.append(to_stringPrecision(eigs[i]));
		if (i < eigs.size() - 1) { output.append(L","); }
		}
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".FrobeniusNorm()") != std::wstring::npos) {
		double FrobeniusNorm = ComplexMatrixContainer[tempIn].FrobeniusNorm();
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".FrobeniusNorm() = ");
		output.append(to_stringPrecision(FrobeniusNorm));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".FTestStatistic()") != std::wstring::npos) {
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".FTestStatistic() = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].FTestStatistic()));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".FTest()") != std::wstring::npos) {
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".FTest() = ");
		output.append(to_stringPrecision(FDistributionApproximation(ComplexMatrixContainer[tempIn].FTestStatistic(), ComplexMatrixContainer[tempIn].columns, ComplexMatrixContainer[tempIn].size() - ComplexMatrixContainer[tempIn].columns)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".GramSchmidt()") != std::wstring::npos) {
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].GramSchmidt());
		return;
		}
		if (in.find(L".GivensRotationComplexMatrix(") != std::wstring::npos) {
		std::vector<double> temp1 = ParseNInputNumbers(in.substr(22));
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].GivensRotationComplexMatrix(temp1[0], temp1[1], temp1[2], temp1[3]));
		return;
		}
		*/	if (in.find(L".getComplexNumber(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(21));
			ComplexNumberContainer.push_back(ComplexMatrixContainer[tempIn].get(temp1[0], temp1[1]));
			return;
		}
		/*	if (in.find(L".HouseholderQR()") != std::wstring::npos) {
		Complexstd::vector<Matrix> qr = ComplexMatrixContainer[tempIn].HouseholderQR();
		ComplexMatrixContainer.push_back(qr[0]);
		ComplexMatrixContainer.push_back(qr[1]);
		qr = NULL;
		return;
		}
		if (in.find(L".inverse()") != std::wstring::npos) {
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].inverse());
		return;
		}
		if (in.find(L".inverseByQR()") != std::wstring::npos) {
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].inverseByQR());
		return;
		}
		*/	if (in.find(L".isBoolean()") != std::wstring::npos) {
			output.append(L"CM");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isBoolean() = ");
			output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].isBoolean()));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".isDiagonalized()") != std::wstring::npos) {
			output.append(L"CM");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isBoisDiagonalizedolean() = ");
			output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].isDiagonalized()));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".isIntegerMatrix()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isIntegerMatrix() = ");
			output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].isIntegerMatrix()));
			output.append(L"\r\n\r\n");
			return;
		}
		/*	if (in.find(L".isOrthonormal()") != std::wstring::npos) {
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".isOrthonormal() = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].isOrthonormal()));
		output.append(L"\r\n\r\n");
		return;
		}
		*/	if (in.find(L".isSymmetric()") != std::wstring::npos) {
			output.append(L"M");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".isSymmetric() = ");
			output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].isSymmetric()));
			output.append(L"\r\n\r\n");
			return;
		}
		/*	if (in.find(L".JacobiRotationComplexMatrix(") != std::wstring::npos) {
		std::vector<double> temp1 = ParseNInputNumbers(in.substr(22));
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].JacobiRotationComplexMatrix(temp1[0], temp1[1], temp1[2], temp1[3]));
		return;
		}
		if (in.find(L".JacobiTransformation()") != std::wstring::npos) {
		Complexstd::vector<Matrix> jt = ComplexMatrixContainer[tempIn].JacobiTransformation();
		ComplexMatrixContainer.push_back(jt[0]);
		ComplexMatrixContainer.push_back(jt[1]);
		jt = NULL;
		return;
		}
		if (in.find(L".L()") != std::wstring::npos) {
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].L());
		return;
		}
		if (in.find(L".LDL()") != std::wstring::npos) {
		Complexstd::vector<Matrix> ldl = ComplexMatrixContainer[tempIn].LDL();
		ComplexMatrixContainer.push_back(ldl[0]);
		ComplexMatrixContainer.push_back(ldl[1]);
		ldl = NULL;
		return;
		}
		if (in.find(L".leastSquares(") != std::wstring::npos) {
		std::vector<double> temp1 = ParseNInputNumbers(in.substr(15));
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].leastSquares(ComplexMatrixContainer[temp1[0] - 1], ComplexMatrixContainer[temp1[1] - 1]));
		return;
		}
		if (in.find(L".LU()") != std::wstring::npos) {
		Complexstd::vector<Matrix> lu = ComplexMatrixContainer[tempIn].LU();
		ComplexMatrixContainer.push_back(lu[0]);
		ComplexMatrixContainer.push_back(lu[1]);
		lu = NULL;
		return;
		}
		if (in.find(L".mean()") != std::wstring::npos) {
		double mean = ComplexMatrixContainer[tempIn].mean();
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".mean() = ");
		output.append(to_stringPrecision(mean));
		output.append(L"\r\n\r\n");
		return;
		}
		*/	if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(13));
			ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].multiply(ComplexMatrixContainer[temp1[0] - 1], ComplexMatrixContainer[temp1[1] - 1]));
			return;
		}
		/*	if (in.find(L".norm()") != std::wstring::npos) {
		double norm = ComplexMatrixContainer[tempIn].norm();
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".norm() = ");
		output.append(to_stringPrecision(norm));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".normalize()") != std::wstring::npos) {
		ComplexMatrixContainer[tempIn].normalize();
		return;
		}
		if (in.find(L".populationStandardDeviation()") != std::wstring::npos) {
		double pop = ComplexMatrixContainer[tempIn].populationStandardDeviation();
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".populationStandardDeviation() = ");
		output.append(to_stringPrecision(pop));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".pseudoinverse()") != std::wstring::npos) {
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].inverse());
		return;
		}
		if (in.find(L".QR()") != std::wstring::npos) {
		Complexstd::vector<Matrix> qr = ComplexMatrixContainer[tempIn].QR();
		ComplexMatrixContainer.push_back(qr[0]);
		ComplexMatrixContainer.push_back(qr[1]);
		qr = NULL;
		return;
		}
		*/	if (in.find(L".randomize()") != std::wstring::npos) {
			ComplexMatrixContainer[tempIn].randomize();
			return;
		}
		if (in.find(L".randomizeBoolean()") != std::wstring::npos) {
			ComplexMatrixContainer[tempIn].randomizeBoolean();
			return;
		}
		if (in.find(L".randomizeInteger()") != std::wstring::npos) {
			ComplexMatrixContainer[tempIn].randomizeInteger();
			return;
		}
		/*	if (in.find(L".randomizeSymmetric()") != std::wstring::npos) {
		ComplexMatrixContainer[tempIn].randomizeSymmetric();
		return;
		}
		if (in.find(L".reducedRowEchelonForm()") != std::wstring::npos) {
		ComplexMatrixContainer[tempIn].reducedRowEchelonForm();
		return;
		}
		if (in.find(L".removeColumn(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(16));
		ComplexMatrixContainer[tempIn].removeColumn(temp1);
		return;
		}
		if (in.find(L".removeRow(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(13));
		ComplexMatrixContainer[tempIn].removeRow(temp1);
		return;
		}
		if (in.find(L".row(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(7));
		std::wcout << in.substr(1) << std::endl;
		VectorContainer.push_back(Vector(ComplexMatrixContainer[tempIn].row(temp1), ComplexMatrixContainer[tempIn].columns));
		return;
		}
		if (in.find(L".rowDeviation(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(17));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".rowDeviation(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].rowDeviation(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".rowMean(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(10));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".rowMean(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].rowMean(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".rowNorm(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(10));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".rowNorm(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].rowNorm(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".rowSumSquares(") != std::wstring::npos) {
		double temp1 = ParseInputNumber(in.substr(16));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".rowSumSquares(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].rowSumSquares(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".rowVariance(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(16));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".rowVariance(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].rowVariance(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".sampleStandardDeviation()") != std::wstring::npos) {
		double std = ComplexMatrixContainer[tempIn].sampleStandardDeviation();
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".sampleStandardDeviation() = ");
		output.append(to_stringPrecision(std));
		output.append(L"\r\n\r\n");
		return;
		}
		*/	if (in.find(L".size()") != std::wstring::npos) {
			double sz = ComplexMatrixContainer[tempIn].size;
			output.append(L"CM");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".size() = ");
			output.append(to_stringPrecision(sz));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".submatrix(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(13));
			ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].submatrix(temp1[0], temp1[1], temp1[2], temp1[3]));
			return;
		}
		/*	if (in.find(L".sum()") != std::wstring::npos) {
		double sum = ComplexMatrixContainer[tempIn].sum();
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".sum() = ");
		output.append(to_stringPrecision(sum));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".sumColumn(") != std::wstring::npos) {
		double temp1 = ParseInputNumber(in.substr(13));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".sumColumn(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].sumColumn(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".sumRow(") != std::wstring::npos) {
		int temp1 = ParseInputNumber(in.substr(10));
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".sumRow(");
		output.append(to_stringPrecision(temp1));
		output.append(L") = ");
		output.append(to_stringPrecision(ComplexMatrixContainer[tempIn].sumRow(temp1)));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".sumSquares()") != std::wstring::npos) {
		double sumSquares = ComplexMatrixContainer[tempIn].sumSquares();
		output.append(L"M");
		output.append(to_stringPrecision(tempIn + 1));
		output.append(L".sumSquares() = ");
		output.append(to_stringPrecision(sumSquares));
		output.append(L"\r\n\r\n");
		return;
		}
		if (in.find(L".swapColumns(") != std::wstring::npos) {
		std::vector<double> temp1 = ParseNInputNumbers(in.substr(15));
		ComplexMatrixContainer[tempIn].swapColumns(temp1[0], temp1[1]);
		return;
		}
		if (in.find(L".swapRows(") != std::wstring::npos) {
		std::vector<double> temp1 = ParseNInputNumbers(in.substr(12));
		ComplexMatrixContainer[tempIn].swapRows(temp1[0], temp1[1]);
		return;
		}
		if (in.find(L".thinQR()") != std::wstring::npos) {
		Complexstd::vector<Matrix> qr = ComplexMatrixContainer[tempIn].thinQR();
		ComplexMatrixContainer.push_back(qr[0]);
		ComplexMatrixContainer.push_back(qr[1]);
		qr = NULL;
		return;
		}
		if (in.find(L".thinHouseholderQR()") != std::wstring::npos) {
		Complexstd::vector<Matrix> qr = ComplexMatrixContainer[tempIn].thinHouseholderQR();
		ComplexMatrixContainer.push_back(qr[0]);
		ComplexMatrixContainer.push_back(qr[1]);
		qr = NULL;
		return;
		}
		*/	if (in.find(L".transpose()") != std::wstring::npos) {
			ComplexMatrixContainer[tempIn] = ComplexMatrixContainer[tempIn].transpose();
			return;
		}
		/*	if (in.find(L".trim()") != std::wstring::npos) {
		ComplexMatrixContainer[tempIn].trim();
		return;
		}
		if (in.find(L".U()") != std::wstring::npos) {
		ComplexMatrixContainer.push_back(ComplexMatrixContainer[tempIn].U());
		return;
		}*/
	}

	//DUAL MATRICES
	//================
	//constructor
	if (in.find(L"DualNumberMatrix(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(14);
		if (tempstr.find(L",") != std::wstring::npos) {   //if there are both a row and column
			std::vector<double> temp = ParseNInputNumbers(tempstr);
			DualNumberMatrix t(temp[0], temp[1]);
			t.identity();
			DualNumberMatrixContainer.push_back(t);
		}
		else {
			double temp = ParseInputNumber(tempstr);
			DualNumberMatrix t(temp, temp);
			t.identity();
			DualNumberMatrixContainer.push_back(t);
		}
		return;
	}
	//Assorted class functions of dual matrices
	//============================================
	if (in[0] == 'D' && in[1] == 'M' && isNumber(in[2]) == true) {    //reference to saved DualNumberMatrix
		double tempIn = ParseInputNumber(in.substr(2, 2));
		tempIn--;


		if (in.find(L".add(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(7));
			DualNumberMatrixContainer.push_back(DualNumberMatrixContainer[tempIn].add(DualNumberMatrixContainer[temp1[0] - 1], DualNumberMatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".conjugateTranspose()") != std::wstring::npos) {
			DualNumberMatrixContainer[tempIn] = DualNumberMatrixContainer[tempIn].conjugateTranspose();
			return;
		}
		if (in.find(L".copy()") != std::wstring::npos) {
			DualNumberMatrixContainer.push_back(DualNumberMatrixContainer[tempIn]);
			return;
		}
		if (in.find(L".delete(") != std::wstring::npos) {
			int temp1 = ParseInputNumber(in.substr(10));
			DualNumberMatrixContainer[tempIn].element.clear();
			DualNumberMatrixContainer.erase(DualNumberMatrixContainer.begin() + temp1);
			return;
		}
		if (in.find(L".exponentiate(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(17));
			DualNumberMatrixContainer.push_back(DualNumberMatrixContainer[tempIn].exponent(DualNumberMatrixContainer[temp1[0] - 1], temp1[1]));
			return;
		}
		if (in.find(L".multiply(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(13));
			DualNumberMatrixContainer.push_back(DualNumberMatrixContainer[tempIn].multiply(DualNumberMatrixContainer[temp1[0] - 1], DualNumberMatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".size()") != std::wstring::npos) {
			double sz = DualNumberMatrixContainer[tempIn].size();
			output.append(L"DM");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".size() = ");
			output.append(to_stringPrecision(sz));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".submatrix(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(13));
			DualNumberMatrixContainer.push_back(DualNumberMatrixContainer[tempIn].submatrix(temp1[0], temp1[1], temp1[2], temp1[3]));
			return;
		}
		if (in.find(L".subtract(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(7));
			DualNumberMatrixContainer.push_back(DualNumberMatrixContainer[tempIn].subtract(DualNumberMatrixContainer[temp1[0] - 1], DualNumberMatrixContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".transpose()") != std::wstring::npos) {
			DualNumberMatrixContainer[tempIn] = DualNumberMatrixContainer[tempIn].transpose();
			return;
		}

	}//END DUAL MATRIX METHODS====================


	 //QUATERNIONS
	 //===========
	 //constructor
	if (in.find(L"Quaternion()") != std::wstring::npos) {
		QuaternionContainer.push_back(Quaternion());
		return;
	}
	if (in.find(L"Quaternion(") != std::wstring::npos) {
		std::vector<double> temp1 = ParseNInputNumbers(in.substr(11));
		Quaternion Q(temp1[0], temp1[1], temp1[2], temp1[3]);
		QuaternionContainer.push_back(Q);
		return;
	}
	//Assorted class functions of quaternions
	//=======================================
	if (in[0] == 'q' && isNumber(in[1]) == true) {    //reference to saved Quaternion
		double tempIn = ParseInputNumber(in.substr(1, 2));
		tempIn--;

		if (in.find(L".copy()") != std::wstring::npos) {
			QuaternionContainer.push_back(QuaternionContainer[tempIn]);
			return;
		}
		if (in.find(L".conjugate()") != std::wstring::npos) {
			QuaternionContainer.push_back(QuaternionContainer[tempIn].conjugate());
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			QuaternionContainer.erase(QuaternionContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".norm(") != std::wstring::npos) {
			double norm = QuaternionContainer[tempIn].norm();
			output.append(L"q");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".norm() = ");
			output.append(to_stringPrecision(norm));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".pNorm(") != std::wstring::npos) {
			double p = ParseInputNumber(in.substr(9));
			double pNorm = VectorContainer[tempIn].pNorm(p);
			output.append(L"q");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".pNorm() = ");
			output.append(to_stringPrecision(pNorm));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".randomize()") != std::wstring::npos) {
			QuaternionContainer[tempIn].randomize();
			return;
		}
		if (in.find(L".randomizeBoolean()") != std::wstring::npos) {
			QuaternionContainer[tempIn].randomizeBoolean();
			return;
		}
		if (in.find(L".randomizeInteger()") != std::wstring::npos) {
			QuaternionContainer[tempIn].randomizeInteger();
			return;
		}
		if (in.find(L".rotate(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(10));
			QuaternionContainer.push_back(QuaternionContainer[tempIn].rotate(QuaternionContainer[temp1[0] - 1], QuaternionContainer[temp1[1] - 1]));
			return;
		}
		if (in.find(L".set(") != std::wstring::npos) {
			std::vector<double> temp1 = ParseNInputNumbers(in.substr(7));
			if (temp1.size() >= 1) { QuaternionContainer[tempIn].element[temp1[0]] = temp1[1]; }
			return;
		}
		if (in.find(L".toMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(QuaternionContainer[tempIn].toMatrix());
			return;
		}
		if (in.find(L".toMatrix(") != std::wstring::npos) {
			std::vector<double> vec = ParseNInputNumbers(in.substr(12));
			MatrixContainer.push_back(QuaternionContainer[tempIn].toMatrix(vec[0], vec[1]));
			return;
		}
		if (in.find(L".toRotationMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(QuaternionContainer[tempIn].toRotationMatrix());
			return;
		}
		if (in.find(L".toScalingMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(QuaternionContainer[tempIn].toScalingMatrix());
			return;
		}
		if (in.find(L".toTranslationMatrix()") != std::wstring::npos) {
			MatrixContainer.push_back(QuaternionContainer[tempIn].toTranslationMatrix());
			return;
		}

	}

	//GROUPS
	//=======
	//constructor
	if (in.find(L"Group(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(6);
		std::vector<int> temp = STLDoubleToInt(ParseNInputNumbers(tempstr));
		Group G(temp);
		GroupContainer.push_back(G);
		return;
	}

	//class functions
	//===============
	if (in[0] == 'G' && isNumber(in[1]) == true) {    //reference to saved Vector
		double tempIn = ParseInputNumber(in.substr(1, 2));
		tempIn--;

		if (in.find(L".copy()") != std::wstring::npos) {
			GroupContainer.push_back(GroupContainer[tempIn]);
			return;
		}
		if (in.find(L".cycleLeft()") != std::wstring::npos) {
			GroupContainer[tempIn].cycleLeft();
			return;
		}
		if (in.find(L".cycleRight()") != std::wstring::npos) {
			GroupContainer[tempIn].cycleRight();
			return;
		}
		if (in.find(L".delete()") != std::wstring::npos) {
			GroupContainer.erase(GroupContainer.begin() + tempIn);
			return;
		}
		if (in.find(L".groupActionLeft(") != std::wstring::npos) {
			std::vector<int> v = STLDoubleToInt(ParseNInputNumbers(in.substr(19)));
			GroupContainer[tempIn].groupActionLeft(v);
			return;
		}
		if (in.find(L".groupActionRight(") != std::wstring::npos) {
			std::vector<int> v = STLDoubleToInt(ParseNInputNumbers(in.substr(20)));
			GroupContainer[tempIn].groupActionRight(v);
			return;
		}
		if (in.find(L".identity()") != std::wstring::npos) {
			GroupContainer[tempIn].identity();
			return;
		}
		if (in.find(L".order()") != std::wstring::npos) {
			double sz = GroupContainer[tempIn].order();
			output.append(L"G");
			output.append(to_stringPrecision(tempIn + 1));
			output.append(L".order() = ");
			output.append(to_stringPrecision(sz));
			output.append(L"\r\n\r\n");
			return;
		}
		if (in.find(L".swap(") != std::wstring::npos) {
			std::vector<int> v = STLDoubleToInt(ParseNInputNumbers(in.substr(15)));
			GroupContainer[tempIn].swap(v[0], v[1]);
			return;
		}
	}

	//==============================================================================================


	//ALL OTHER FUNCTIONS
	//===================
	//Ackermann's Function
	if (in.find(L"AckermannsFunction(") != std::wstring::npos) {
		std::wstring substr = in.substr(19);
		std::vector<double> tempIn = ParseNInputNumbers(substr);
		double x = AckermannsFunction(tempIn[0], tempIn[1]);
		output.append(L"AckermannsFunction(");
		output.append(std::to_wstring(tempIn[0]));
		output.append(L",");
		output.append(std::to_wstring(tempIn[1]));
		output.append(L") = ");
		output.append(std::to_wstring(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//Bell Polynomial
	if (in.find(L"BellPolynomial(") != std::wstring::npos) {
		std::wstring substr = in.substr(15);
		int tempIn = ParseInputNumber(substr);
		PolynomialContainer.push_back(BellPolynomial(tempIn));
		return;
	}
	//Bell Number
	if (in.find(L"BellNumber(") != std::wstring::npos) {
		std::wstring substr = in.substr(11);
		int tempIn = ParseInputNumber(substr);
		long long unsigned int x = BellNumber(tempIn);
		output.append(L"BellNumber(");
		output.append(std::to_wstring(tempIn));
		output.append(L") = ");
		output.append(std::to_wstring(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//Bernoulli Number
	if (in.find(L"BernoulliNumber(") != std::wstring::npos) {
		std::wstring substr = in.substr(16);
		int tempIn = ParseInputNumber(substr);
		double x = BernoulliNumber(tempIn);
		output.append(L"BernoulliNumber(");
		output.append(std::to_wstring(tempIn));
		output.append(L") = ");
		output.append(to_stringPrecision(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//Bernoulli Polynomial
	if (in.find(L"BernoulliPolynomial(") != std::wstring::npos) {
		std::wstring substr = in.substr(16);
		int tempIn = ParseInputNumber(substr);
		PolynomialContainer.push_back(BernoulliPolynomial(tempIn));
		return;
	}
	//Bessel Polynomial
	if (in.find(L"BesselPolynomial(") != std::wstring::npos && in.find(L"reverse") == std::wstring::npos) {
		std::wstring substr = in.substr(16);
		int tempIn = ParseInputNumber(substr);
		PolynomialContainer.push_back(BesselPolynomial(tempIn));
		return;
	}
	//CauchyCDF
	if (in.find(L"CauchyCDF(") != std::wstring::npos) {
		std::vector<double> tempIn = ParseNInputNumbers(in.substr(10));
		output.append(L"CauchyCDF(");
		output.append(to_stringPrecision(tempIn[0]));
		output.append(L", " + to_stringPrecision(tempIn[1]));
		output.append(L", " + to_stringPrecision(tempIn[2]));
		double x = CauchyCDF(tempIn[0], tempIn[1], tempIn[2]);
		output.append(L") = ");
		output.append(std::to_wstring(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//Chebyshev Polynomial of the 1st Kind
	if (in.find(L"ChebyshevPolynomial") != std::wstring::npos && in.find(L"2ndKind") == std::wstring::npos) {
		std::wstring tempstr = in.substr(27);
		while (tempstr.find(L")") != std::wstring::npos) { tempstr.pop_back(); }
		int tempIn = ParseInputNumber(tempstr);
		PolynomialContainer.push_back(ChebyshevPolynomial1stKind(tempIn));
		return;
	}
	//Chebyshev Polynomial of the 2nd Kind
	if (in.find(L"ChebyshevPolynomial2ndKind(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(27);
		while (tempstr.find(L")") != std::wstring::npos) { tempstr.pop_back(); }
		int tempIn = ParseInputNumber(tempstr);
		PolynomialContainer.push_back(ChebyshevPolynomial2ndKind(tempIn));
		return;
	}
	//chiSquareCDF
	if (in.find(L"chiSquareCDF(") != std::wstring::npos) {
		std::vector<double> tempIn = ParseNInputNumbers(in.substr(13));
		output.append(L"chiSquareCDF(");
		output.append(to_stringPrecision(tempIn[0]));
		output.append(L", " + to_stringPrecision(tempIn[1]));
		double x = chiSquareCDF(tempIn[0], tempIn[1]);
		output.append(L") = ");
		output.append(std::to_wstring(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//chiSquareInverse -- get test statistic for given probability
	if (in.find(L"chiSquareInverse(") != std::wstring::npos) {
		std::vector<double> tempIn = ParseNInputNumbers(in.substr(17));
		output.append(L"chiSquareInverse(");
		output.append(to_stringPrecision(tempIn[0]));
		output.append(L", " + to_stringPrecision(tempIn[1]));
		double x = chiSquareCDFInverseApproximation(tempIn[0], tempIn[1]);
		output.append(L") = ");
		output.append(std::to_wstring(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//convertRadiansToDegrees(double n);
	if (in.find(L"convertRadiansToDegrees(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(24);
		double temp = ParseInputNumber(tempstr);
		double result = convertRadiansToDegrees(temp);
		output.append(L"convertRadiansToDegrees(");
		output.append(tempstr);
		output.append(L" = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		return;
	}
	//convertDegreesToRadians(double n);
	if (in.find(L"convertDegreesToRadians(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(24);
		double temp = ParseInputNumber(tempstr);
		double result = convertDegreesToRadians(temp);
		output.append(L"convertDegreesToRadians(");
		output.append(tempstr);
		output.append(L" = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		return;
	}
	//convertToBinary(int n);
	if (in.find(L"convertToBinary(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(15);
		double temp = ParseInputNumber(tempstr);
		double result = convertToBinary(temp);
		output.append(L"convertToBinary(");
		output.append(tempstr);
		output.append(L" = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		return;
	}
	//Divisors of an integer
	if (in.find(L"divisors(") != std::wstring::npos) {
		std::wstring substr = in.substr(9);
		int tempIn = ParseInputNumber(substr);
		std::vector<int> x = divisors(tempIn);
		output.append(L"divisors of (");
		output.append(std::to_wstring(tempIn));
		output.append(L") = {");
		for (int i = 0; i < x.size(); ++i) {
			output.append(std::to_wstring(x[i]));
			if (i < x.size() - 1) { output.append(L","); }
		}
		output.append(L"}\r\n\r\n");
		return;
	}
	//Euler Polynomial
	if (in.find(L"EulerPolynomial(") != std::wstring::npos) {
		std::wstring substr = in.substr(16);
		int tempIn = ParseInputNumber(substr);
		PolynomialContainer.push_back(EulerPolynomial(tempIn));
		return;
	}
	//Euler Number
	if (in.find(L"EulerNumber(") != std::wstring::npos) {
		std::wstring substr = in.substr(12);
		int tempIn = ParseInputNumber(substr);
		double x = EulerNumber(tempIn);
		output.append(L"EulerNumber(");
		output.append(std::to_wstring(tempIn));
		output.append(L") = ");
		output.append(std::to_wstring(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//Euler's Totient Function
	if (in.find(L"EulersTotientFunction(") != std::wstring::npos) {
		std::wstring substr = in.substr(22);
		int tempIn = ParseInputNumber(substr);
		int x = EulersTotientFunction(tempIn);
		output.append(L"EulersTotientFunction(");
		output.append(to_stringPrecision(tempIn));
		output.append(L") = ");
		output.append(to_stringPrecision(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//Fibonacci Number
	if (in.find(L"FibonacciNumber(") != std::wstring::npos) {
		std::wstring substr = in.substr(16);
		int tempIn = ParseInputNumber(substr);
		double x = FibonacciNumber(tempIn);
		output.append(L"FibonacciNumber(");
		output.append(std::to_wstring(tempIn));
		output.append(L") = ");
		output.append(std::to_wstring(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//HermitePolynomial
	if (in.find(L"HermitePolynomial(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(18);
		while (tempstr.find(L")") != std::wstring::npos) { tempstr.pop_back(); }
		int tempIn = ParseInputNumber(tempstr);
		PolynomialContainer.push_back(HermitePolynomial(tempIn));
		return;
	}
	//Hermite Number
	if (in.find(L"HermiteNumber(") != std::wstring::npos) {
		std::wstring substr = in.substr(14);
		int tempIn = ParseInputNumber(substr);
		double x = HermiteNumber(tempIn);
		output.append(L"HermiteNumber(");
		output.append(to_stringPrecision(tempIn));
		output.append(L") = ");
		output.append(to_stringPrecision(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//inverseErf
	if (in.find(L"inverseErf(") != std::wstring::npos) {
		std::wstring substr = in.substr(10);
		double tempIn = ParseInputNumber(substr);
		output.append(L"inverseErf(");
		output.append(to_stringPrecision(tempIn));
		double x = inverseErfPadeApproximation(tempIn);
		output.append(L") = ");
		output.append(to_stringPrecision(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//inverseErfc
	if (in.find(L"inverseErfc(") != std::wstring::npos) {
		std::wstring substr = in.substr(10);
		double tempIn = ParseInputNumber(substr);
		output.append(L"inverseErfc(");
		output.append(to_stringPrecision(tempIn));
		double x = inverseErfc(tempIn);
		output.append(L") = ");
		output.append(to_stringPrecision(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//isPrime(int n);
	if (in.find(L"isPrime(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(7);
		unsigned long long int temp = ParseInputNumber(tempstr);
		bool result = isPrime(temp);
		output.append(L"isPrime(");
		output.append(tempstr);
		output.append(L" = ");
		if (result == 0) { output.append(L"false"); }
		if (result == 1) { output.append(L"true"); }
		output.append(L"\r\n\r\n");
		return;
	}
	//Lagrange Interpolating Polynomial
	if (in.find(L"LagrangeInterpolatingPolynomial(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(32);
		while (tempstr.find(L")") != std::wstring::npos) { tempstr.pop_back(); }
		std::vector<double> tempIn = ParseNInputNumbers(tempstr);
		std::vector<double> x;
		std::vector<double> y;
		for (int i = 0; i < tempIn.size(); i += 2) {
			x.push_back(tempIn[i]);
			y.push_back(tempIn[i + 1]);
		}
		PolynomialContainer.push_back(LagrangeInterpolatingPolynomial(x, y));
		return;
	}
	//Laguerre Polynomial
	if (in.find(L"LaguerrePolynomial(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(19);
		while (tempstr.find(L")") != std::wstring::npos) { tempstr.pop_back(); }
		int tempIn = ParseInputNumber(tempstr);
		PolynomialContainer.push_back(LaguerrePolynomial(tempIn));
		return;
	}
	//LegendrePolynomial
	if (in.find(L"LegendrePolynomial(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(19);
		while (tempstr.find(L")") != std::wstring::npos) { tempstr.pop_back(); }
		int tempIn = ParseInputNumber(tempstr);
		PolynomialContainer.push_back(LegendrePolynomial(tempIn));
		return;
	}
	//MoebiusFunction
	if (in.find(L"MoebiusFunction(") != std::wstring::npos) {
		std::wstring substr = in.substr(16);
		int tempIn = ParseInputNumber(substr);
		output.append(L"MoebiusFunction(");
		output.append(to_stringPrecision(tempIn));
		int x = MoebiusFunction(tempIn);
		output.append(L") = ");
		output.append(std::to_wstring(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//normalCDF
	if (in.find(L"normalCDF(") != std::wstring::npos) {
		std::wstring substr = in.substr(10);
		double tempIn = ParseInputNumber(substr);
		output.append(L"normalCDF(");
		output.append(to_stringPrecision(tempIn));
		double x = normalCDF(tempIn);
		output.append(L") = ");
		output.append(to_stringPrecision(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//Partition Number
	if (in.find(L"partitionNumber(") != std::wstring::npos) {
		std::wstring substr = in.substr(15);
		if (substr.find(',') != std::wstring::npos) {
			std::vector<double> temp = ParseNInputNumbers(substr);
			output.append(L"partitionNumber(");
			output.append(to_stringPrecision(temp[0]));
			output.append(L",");
			output.append(to_stringPrecision(temp[1]));
			int x = partitionNumber(temp[0], temp[1]);
			output.append(L") = ");
			output.append(to_stringPrecision(x));
			output.append(L"\r\n\r\n");
			return;
		}
		int tempIn = ParseInputNumber(substr);
		output.append(L"partitionNumber(");
		output.append(to_stringPrecision(tempIn));
		int x = partitionNumber(tempIn);
		output.append(L") = ");
		output.append(to_stringPrecision(x));
		output.append(L"\r\n\r\n");
		return;
	}
	//pentagonalNumber(double n);
	if (in.find(L"pentagonalNumber(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(17);
		double temp = ParseInputNumber(tempstr);
		double result = pentagonalNumber(temp);
		output.append(L"pentagonalNumber(");
		output.append(tempstr);
		output.append(L" = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		return;
	}
	//percentError(double n, double r);
	if (in.find(L"percentError(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(13);
		std::vector<double> temp = ParseNInputNumbers(tempstr);
		double result = percentError(temp[0], temp[1]);
		output.append(L"percentError(");
		output.append(tempstr);
		output.append(L" = ");
		output.append(to_stringPrecision(result));
		output.append(L"\r\n\r\n");
		return;
	}
	//Periodic Table(std::wstring Command, int n);
	if (in.find(L"PeriodicTable()") != std::wstring::npos) {
		output.append(PeriodicTable() + L"\r\n\r\n");
		return;
	}
	if (in.find(L"PeriodicTable(") != std::wstring::npos) {
		std::wstring tempstr = in.substr(14);
		std::vector<double> temp = ParseNInputNumbers(tempstr);
		if (temp.size() == 0) { temp.push_back(0); }
		if (temp.size() == 1) { temp.push_back(0); }
		std::wstring result = PeriodicTable(std::to_wstring((int)temp[0]), temp[1]);
		output.append(result);
		output.append(L"\r\n\r\n");
		return;
	}
	//Reverse Bessel Polynomial
	if (in.find(L"reverseBesselPolynomial(") != std::wstring::npos) {
		std::wstring substr = in.substr(23);
		int tempIn = ParseInputNumber(substr);
		PolynomialContainer.push_back(reverseBesselPolynomial(tempIn));
		return;
	}
	//Rogers-Szego Polynomial
	if (in.find(L"RogersSzegoPolynomial(") != std::wstring::npos) {
		std::wstring substr = in.substr(23);
		std::vector<double> tempIn = ParseNInputNumbers(substr);
		PolynomialContainer.push_back(RogersSzegoPolynomial(tempIn[0], tempIn[1]));
		return;
	}
	//studentsCDF
	if (in.find(L"studentsCDF(") != std::wstring::npos) {
		std::wstring substr = in.substr(12);
		std::vector<double> tempIn = ParseNInputNumbers(substr);
		output.append(L"studentsCDF(");
		output.append(to_stringPrecision(tempIn[0]));
		output.append(L", " + to_stringPrecision(tempIn[1]));
		double x = studentsCDF(tempIn[0], tempIn[1]);
		output.append(L") = ");
		output.append(to_stringPrecision(x));
		output.append(L"\r\n\r\n");
		return;
	}

	//END ORDINARY FUNCTIONS==============================

	//ATLAS OF FINITE GROUPS
	//(matrix representations)
	//========================

	//Sporadic Groups
	//===============
	if (in == L"M11") { M11(); return; }
	if (in == L"M12") { M12(); return; }
	if (in == L"M22") { M22(); return; }
	if (in == L"M23") { M23(); return; }
	if (in == L"M24") { M24(); return; }
	if (in == L"HS") { HS(); return; }
	if (in == L"J2") { J2(); return; }
	if (in == L"Co1") { Co1(); return; }
	if (in == L"Co2") { Co2(); return; }
	if (in == L"Co3") { Co3(); return; }
	if (in == L"McL") { McL(); return; }
	if (in == L"Suz") { Suz(); return; }
	if (in == L"He") { He(); return; }
	if (in == L"HN") { HN(); return; }
	if (in == L"Th") { Th(); return; }
	if (in == L"Fi22") { Fi22(); return; }
	if (in == L"Fi23") { Fi23();	 return; }
	if (in == L"Fi24") { Fi24(); return; }
	if (in == L"B") { B(); return; }
	if (in == L"M") { M(); return; }
	if (in == L"J1") { J1(); return; }
	if (in == L"ON") { ON(); return; }
	if (in == L"J3") { J3(); return; }
	if (in == L"Ru") { Ru();return; }
	if (in == L"J4") { J4();return; }
	if (in == L"Ly") { Ly();return; }
	if (in == L"T") { T();return; }

	//Alternating groups
	//==================
	if (in == L"A5") { A5();return; }
	if (in == L"A6") { A6();return; }
	if (in == L"A7") { A7();return; }
	if (in == L"A8") { A8();return; }
	if (in == L"A9") { A9();return; }
	if (in == L"A10") { A10();return; }
	if (in == L"A11") { A11();return; }
	if (in == L"A12") { A12();return; }
	if (in == L"A13") { A13();return; }
	if (in == L"A14") { A14();return; }
	if (in == L"A15") { A15();return; }
	if (in == L"A16") { A16();return; }
	if (in == L"A17") { A17();return; }
	if (in == L"A18") { A18();return; }
	if (in == L"A19") { A19();return; }
	if (in == L"A20") { A20();return; }
	if (in == L"A21") { A21();return; }
	if (in == L"A22") { A22();return; }
	if (in == L"A23") { A23();return; }

	//Linear Groups
	//=============
	if (in == L"L3(2)") { L32();return; }

	//Classical Groups
	//================

	//Exceptional Groups of Lie Type
	//==============================

	//Miscellaneous Groups
	//=====================
	if (in == L"M20") { M20(); return; }
	
	//Finally, if all else fails treat input like a normal calculator:
	//PARSE BY INPUT TYPE
	if (hasOperator(in) == true) {//handle numeric input like a normal calculator would, including functions like 'exp(2.34)'
		if (findComplexNumber(in) != std::wstring::npos) {//process expression w/complex numbers in it
			ComplexFunction fp(L"", in);
			ComplexNumber x(0);
			ComplexNumber result = fp.evaluateComplex(x);
			output.append(in);
			output.append(L" = ");
			output.append(result.toString());
			output.append(L"\r\n\r\n");
			return;
		}else {//process real-valued expression
			Function fp(L"", in);
			double x = 0;
			double result = fp.evaluate(x);
			output.append(in);
			output.append(L" = ");
			output.append(to_stringPrecision(result));
			output.append(L"\r\n\r\n");
			return;
		}
	}

	//====================================
	//if nothing is thrown, return nothing
	in = L"";
	output = L"";
}