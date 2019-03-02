#pragma once
#include "stdafx.h"

double FourierTransformAngularFrequency(double w, Function f) //w must be a real number
{   /* evaluates the FourierTF of some function
	in the time domain at a frequency, w.
	Assumes w is a real valued variable only.
	Integral bounds are [-infinity,infinity]

	FTF = integral [-inf,inf]: ( f(x)* [cos(2*PI*w*x) - i sin(2*PI*w*x)] ) dx     */

	const int n = 50 * 3;	//use n abscissa points, the roots of Hermite Polynomial P_n(x)
	std::vector<double> x = GaussHermiteTable();
	double answer = 0;
	//Hermite polynomial coefficients are symmetric -- we must run the loop once with negative values,
	//and once with negative values.  The following method is a neat way to perform both operations at once.
	for (int i = 0; i < n; i += 3) {
		double temp = x[i + 2] * f.evaluate(-x[i]) * cos(w*-x[i]);//sum w_i*f(x_i)*cos(-2Pi * wx)
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
		temp = x[i + 2] * f.evaluate(x[i]) * cos(w*x[i]);//sum w_i*f(x_i)*cos(-2Pi * wx)
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
	}
	x.clear();
	return answer;
}

double FourierTransform(double w, Function f) //w must be a real number
{   /* evaluates the FourierTF of some function
	in the time domain at a frequency, w.
	Assumes w is a real valued variable only.
	Integral bounds are [-infinity,infinity]

	FTF = integral [-inf,inf]: ( f(x)* [cos(2*PI*w*x) - i sin(2*PI*w*x)] ) dx     */

	const int n = 50 * 3;	//use n abscissa points, the roots of Hermite Polynomial P_n(x)
	std::vector<double> x = GaussHermiteTable();
	double answer = 0;
	//Hermite polynomial coefficients are symmetric -- we must run the loop once with negative values,
	//and once with negative values.  The following method is a neat way to perform both operations at once.
	for (int i = 0; i < n; i += 3) {
		double temp = x[i + 2] * f.evaluate(-x[i]) * cos(2 * PI*w*-x[i]);//sum w_i*f(x_i)*cos(-2Pi * wx)
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
		temp = x[i + 2] * f.evaluate(x[i]) * cos(2 * PI*w*x[i]);//sum w_i*f(x_i)*cos(-2Pi * wx)
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
	}
	x.clear();
	return answer;
}

double inverseFourierTransform(double w, Function f) {
	/*Note:	Since we have defined the Fourier over R as simply the cosine term of the
	Fourier TF integral, the same formula works for the inverse.  Thus, we do not need to
	calculate the '+ i sin(2*PI*w*x)' term since we don't care about the imaginary value.	*/
	return FourierTransform(w, f);
}

double LaplaceTransform(double s, Function f) {
	/* evaluates the LaplaceTF of some function
	in the time domain at a complex frequency, s.*/

	const int n = 100 * 3;	//use n abscissa points, the roots of Laguerre Polynomial P_n(x)
	std::vector<long double> x = GaussLaguerreTable();
	double answer = 0;
	for (int i = 0; i < n; i += 3) {
		double temp = x[i + 2] * f.evaluate(x[i]) * exp(-1 * s*x[i]);//sum w_i*f(x_i)*e^-sx
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
	}
	x.clear();
	return answer;
}

double twoSidedLaplaceTransform(double s, Function f) {
	/* evaluates the LaplaceTF of some function
	in the time domain at a complex frequency, s.*/

	const int n = 50 * 3;	//use n abscissa points, the roots of Hermite Polynomial P_n(x)
	std::vector<double> x = GaussHermiteTable();
	double answer = 0;
	for (int i = 0; i < n; i += 3) {
		double temp = x[i + 2] * f.evaluate(x[i]) * exp(-1 * s*x[i]);//sum w_i*f(x_i)*e^-sx
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
	}
	x.clear();
	return answer;
}

double MellinTransform(double s, Function f)
{   /* evaluates the MellinTF of some function
	in the time domain at a complex frequency, s.*/

	const int n = 100 * 3;	//use n abscissa points, the roots of Laguerre Polynomial P_n(x)
	std::vector<long double> x = GaussLaguerreTable();
	double answer = 0;
	for (int i = 0; i < n; i += 3) {
		double temp = x[i + 2] * f.evaluate(x[i]) * pow(x[i], s - 1);//sum w_i*f(x_i)*x^s-1
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
	}
	x.clear();
	return answer;
}

double twoSidedMellinTransform(double s, Function f)
{   /* evaluates the MellinTF of some function
	in the time domain at a complex frequency, s.*/

	const int n = 50 * 3;	//use n abscissa points, the roots of Laguerre Polynomial P_n(x)
	std::vector<double> x = GaussHermiteTable();
	double answer = 0;
	for (int i = 0; i < n; i += 3) {
		double temp = x[i + 2] * f.evaluate(x[i]) * pow(x[i], s - 1);//sum w_i*f(x_i)*x^s-1
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
	}
	x.clear();
	return answer;
}

double convolve1D(double t, Function f1, Function f2) {   //convolution is equivalent to multiplication in the frequency domain.
	const int n = 50 * 3;	//use n abscissa points, the roots of Hermite Polynomial P_n(x)
	std::vector<double> x = GaussHermiteTable();
	double answer = 0;
	//Hermite polynomial coefficients are symmetric -- we must run the loop once with negative values,
	//and once with negative values.  The following method is a neat way to perform both operations at once.
	for (int i = 0; i < n; i += 3) {
		double temp = x[i + 2] * f1.evaluate(t) * f2.evaluate(x[i] - t);//sum w_i*f(x_i)*cos(-2Pi * wx)
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
		temp = x[i + 2] * f1.evaluate(t) * f2.evaluate(x[i] - t);//sum w_i*f(x_i)*cos(-2Pi * wx)
		if (temp == temp) { answer += temp; }//reject crazy values (i.e. singularities)
	}
	x.clear();
	return answer;
}

std::vector<ComplexNumber> discreteFourierTransform(std::vector<double> data) {
	std::vector<ComplexNumber> c;
	for (int i = 0; i < data.size(); ++i) {
		c.push_back(ComplexNumber(data[i]));
	}
	return discreteFourierTransform(c);
}

std::vector<ComplexNumber> discreteFourierTransform(std::vector<ComplexNumber> data) {
	int N = data.size();
	std::vector<ComplexNumber> data2;
	for (int k = 0; k < data.size(); ++k) {
		std::complex<double> answer(0);
		for (int i = 0; i <= N-1; ++i) {
			std::complex<double> c = data[i].toComplex();
			double val = 2*PI*k*i / N;
			std::complex<double> z(cos(val),-sin(val));
			c *= z;
			answer += c;
		}
		data2.push_back(answer);
	}
	return data2;
}

std::vector<ComplexNumber> discreteFourierTransform(std::vector<double> x, Function f) {
	std::vector<double> vec;
	for (int i = 0; i < x.size(); ++i) {
		vec.push_back(f.evaluate(x[i]));
	}
	return discreteFourierTransform(vec);
}

std::vector<ComplexNumber> discreteFourierTransform(std::vector<ComplexNumber> x, ComplexFunction f) {
	std::vector<ComplexNumber> vec;
	for (int i = 0; i < x.size(); ++i) {
		vec.push_back(f.evaluateComplex(x[i]));
	}
	return discreteFourierTransform(vec);
}

ComplexMatrix discreteFourierTransform(Matrix A) {
	ComplexMatrix CM = A.toComplexMatrix();
	return discreteFourierTransform(CM);
}

ComplexMatrix discreteFourierTransform(ComplexMatrix A) {//computer DFT of input vector by matrix transform
	//pad matrix A if non-square
/*	if (A.rows != A.columns) {
		while (A.rows > A.columns) {
			A = A.extendColumns(1);
		}
		while (A.rows < A.columns) {
			A = A.extendRows(1);
		}
	}*/

	int N = A.rows;
	double val = 2 * PI / N;
	ComplexMatrix DFTMat(N, N);
	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j) {
				ComplexNumber temp = pow(ComplexNumber(cos(val),-sin(val)), i*j);
				DFTMat.element[i*N + j] = temp;
				if (i != j) { DFTMat.element[j*N + i] = temp; }
		}
	}
	return A.multiply(DFTMat,A);
}

std::vector<ComplexNumber> inverseDiscreteFourierTransform(std::vector<double> data) {
	std::vector<ComplexNumber> c;
	for (int i = 0; i < data.size(); ++i) {
		c.push_back(ComplexNumber(data[i]));
	}
	return inverseDiscreteFourierTransform(c);
}

std::vector<ComplexNumber> inverseDiscreteFourierTransform(std::vector<ComplexNumber> data) {
	double N = data.size();
	std::vector<ComplexNumber> data2;
	for (int k = 0; k < data.size(); ++k) {
		std::complex<double> answer(0);
		for (int i = 0; i <= N - 1; ++i) {
			std::complex<double> c = data[i].toComplex();
			double val = 2 * PI*k*i / N;
			std::complex<double> z(cos(val), sin(val));
			c *= z;
			answer += c;
		}
		data2.push_back(answer/N);
	}
	return data2;
}

ComplexMatrix inverseDiscreteFourierTransform(Matrix A) {
	ComplexMatrix CM = A.toComplexMatrix();
	return inverseDiscreteFourierTransform(CM);
}

ComplexMatrix inverseDiscreteFourierTransform(ComplexMatrix A) {//computer inverse DFT of input vector by matrix transform
	int N = A.rows;
	double val = 2 * PI / N;
	ComplexMatrix DFTMat(N,N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			DFTMat.element[i*N + j] = pow(ComplexNumber(cos(val), sin(val)), i*j);
		}
	}
	DFTMat = DFTMat.multiply(DFTMat, 1.0 / N);
	return A.multiply(DFTMat,A);
}

std::vector<ComplexNumber> fastFourierTransform(std::vector<ComplexNumber> data) {
	ComplexNumber* dta = new ComplexNumber[data.size()];
	ComplexNumber* temp = new ComplexNumber[data.size()];
	fft(dta, temp, data.size(), data.size());
	std::vector<ComplexNumber> vec = toSTLVector(dta, data.size());
	return vec;
}



void fft(ComplexNumber *data, ComplexNumber *temp, int n, int ndata) {
	ComplexNumber circle[((1 << 18) / 2) + 1];

	int seqlen;

	//knock off least significant digit of data
	seqlen = ndata;
	if (seqlen > 0) {
		while (!isPowerOfTwo(seqlen)) { seqlen &= (seqlen - 1); }	/* knock off least sig. 1 bit */
		while (!seqlen >= ndata) { seqlen <<= 1; }
		/* seqlen is now the least power of 2 .ge. ndata */
		while (ndata < seqlen) data[ndata++] = 0.0;
	}

	//set up circle
	for (int i = 0; i <= seqlen / 2; i++) { /* actually a semicircle; n.b. extra slot for window computation */
		double x = (2 * PI * i) / seqlen;
		circle[i].Re = cos(x);
		circle[(seqlen / 2) - i - 1].Re = cos(x);
		circle[i].Im = sin(x);
		circle[(seqlen / 2) - i - 1] = -sin(x);
	}

	if (n > 1) {
		int h = n / 2;
		for (int i = 0; i < h; i++) {
			int i2 = i * 2;
			temp[i] = data[i2];				/* even */
			temp[h + i] = data[i2 + 1];		/* odd	*/
		}
		fft(&temp[0], &data[0], h, seqlen);
		fft(&temp[h], &data[h], h, seqlen);
		int p = 0, t = seqlen / n;
		for (int i = 0; i < h; i++)
		{
			ComplexNumber wkt = circle[p] * temp[h + i];
			data[i] = temp[i] + wkt;
			data[h + i] = temp[i] - wkt;
			p += t;
		}
	}
}

Matrix discreteCosineTransform(std::vector<double> data) {
	//for examples, see: http://eeweb.poly.edu/~yao/EE3414/ImageCoding_DCT.pdf
	int N = data.size();
	if (N < 2) { return Matrix(); }
	std::vector<double> data2(N*N);	
	for (int k = 0; k <= N-1; ++k) {
		for (int i = 0; i <= N - 1; ++i) {
			double answer = 0;
			answer = cos(PI*k*(2 * i + 1) / (2 * N));
			if (k == 0) { answer *= sqrt(1.0 / N); }
			else { answer *= sqrt(2.0 / N); }
			data2[k*N + i] = answer;
		}
	}
	Matrix DCT(N,N,data2);
	Matrix vec(N, 1, data);
	return DCT.multiply(DCT, vec);
}

Matrix discreteCosineTransform(Matrix data) {
	//for examples, see: http://eeweb.poly.edu/~yao/EE3414/ImageCoding_DCT.pdf
	if (data.rows < 2 || data.columns < 1) { return Matrix(); }
	if (data.columns == 1) { 
		std::vector<double> vec;
		for (int i = 0; i < data.rows; ++i) {
			vec.push_back(data.element[i]);
		}
		return discreteCosineTransform(vec); 
	}

	//initialize all vectors u_n
	std::vector<std::vector<double>> vecs;
	int N = data.rows;
	for (int k = 0; k <= N - 1; ++k) {
		std::vector<double> data2(N);
		for (int i = 0; i <= N - 1; ++i) {
			double answer = 0;
			answer = cos(PI*k*(2 * i + 1) / (2 * N));
			if (k == 0) { answer *= sqrt(1.0 / N); }
			else { answer *= sqrt(2.0 / N); }
			data2[i] = answer;
		}
		vecs.push_back(data2);
	}

	//create matrices by calculating outerProduct(u_i, u_j), multiplying by data matrix, and summing the resulting values
	std::vector<double> elements;
	for (int k = 0; k <= N - 1; ++k) {
		for (int l = 0; l <= N - 1; ++l) {
			Vector v1(vecs[k]);
			Vector v2(vecs[l]);
			Matrix DCT = v1.outerProduct(v1, v2);
			Matrix Hadm = DCT.HadamardProduct(DCT, data);
			elements.push_back(Hadm.sumAll());
		}
	}
	return Matrix(data.rows, data.columns, elements);
}

Matrix inverseDiscreteCosineTransform(std::vector<double> X){
	//for examples, see: http://eeweb.poly.edu/~yao/EE3414/ImageCoding_DCT.pdf
	int N = X.size();
	if (N < 2) { return Matrix(); }
	std::vector<double> data2(N*N);
	for (int k = 0; k <= N - 1; ++k) {
		for (int i = 0; i <= N - 1; ++i) {
			double answer = 0;
			answer = cos(PI*k*(2 * i + 1) / (2 * N));
			if (k == 0) { answer *= sqrt(1.0 / N); }
			else { answer *= sqrt(2.0 / N); }
			data2[k*N + i] = answer;
		}
	}
	Matrix DCT(N, N, data2);
	Matrix DCT2 = DCT.transpose();
	Matrix vec(N, 1, X);
	return DCT.multiply(DCT2, vec);
}

Matrix inverseDiscreteCosineTransform(Matrix data) {
	//for examples, see: http://eeweb.poly.edu/~yao/EE3414/ImageCoding_DCT.pdf
/*	if (data.rows < 2 || data.columns < 1) { return Matrix(); }
	if (data.columns == 1) {
		std::vector<double> vec;
		for (int i = 0; i < data.rows; ++i) {
			vec.push_back(data.element[i]);
		}
		return inverseDiscreteCosineTransform(vec);
	}

	//initialize all vectors u_n
	std::vector<std::vector<double>> vecs;
	int N = data.rows;
	for (int k = 0; k <= N - 1; ++k) {
		std::vector<double> data2(N);
		for (int i = 0; i <= N - 1; ++i) {
			double answer = 0;
			answer = cos(PI*k*(2 * i + 1) / (2 * N));
			if (k == 0) { answer *= sqrt(1.0 / N); }
			else { answer *= sqrt(2.0 / N); }
			data2[i] = answer;
		}
		vecs.push_back(data2);
	}

	//create matrices by calculating outerProduct(u_i, u_j), multiplying by data matrix, and summing the resulting values
	std::vector<double> elements;
	for (int k = 0; k <= N - 1; ++k) {
		for (int l = 0; l <= N - 1; ++l) {
			Vector v1(vecs[k]);
			Vector v2(vecs[l]);
			Matrix DCT = v1.outerProduct(v2, v1);
			//Matrix DCT2 = DCT.inverse();	
			Matrix Hadm = DCT.HadamardProduct(DCT, data);
			elements.push_back(Hadm.sumAll());
		}
	}
	return Matrix(data.rows, data.columns, elements);	*/

	std::vector<std::vector<double>> elm;
	for (int m = 0; m < data.columns; ++m) {
		std::vector<double> vec = data.column(m);
		Matrix M = inverseDiscreteCosineTransform(vec);
		vec.clear();
		for (int i = 0; i < M.rows; ++i) {
			vec.push_back(M.element[i]);
		}
		elm.push_back(vec);
	}
	return Matrix(elm);
}