#pragma once
#include "stdafx.h"

ComplexNumber mean(std::vector<ComplexNumber> vec) {
	std::complex<double> c;
	for (int i = 0; i < vec.size(); ++i) {
		c += vec[i].toComplex();
	}
	return ComplexNumber(c / (double)vec.size());
}

double mean(int a, int b, Function f) {
	return ((1.0 / (b - a)) * integrateGaussLegendreQuadrature(a, b, f));
}

double logisticFunction(double n) {// (1+e^-f(x))^-1
	return 1.0 / (1 + exp(-n));
}

double logisticFunction(double x, Function f) {// (1+e^-f(x))^-1
	double n = f.evaluate(x);
	return 1.0 / (1 + exp(-n));
}

double logit(double x) {//inverse of logistic function
	return log(x / (1.0 - x));
}

double logit(double x, Function f) {//inverse of logistic function
	double n = f.evaluate(x);
	return log(n / (1.0 - n));
}

double Gaussian1D(double x, double o) {
	return exp(-0.5*pow(x/o,2))/(o*sqrt(2*PI));
}

double Gaussian1D(double x, double u, double o) {
	return exp(-0.5*pow((x-u) / o, 2))/(o*sqrt(2*PI));
}

double Gaussian2D(double x, double y, double o) {
	return (exp(-1.0*((pow(x, 2) + pow(y, 2)) / (2 * pow(o, 2)))) / ((2 * PI)*pow(o,2)));
}

double Gaussian2D(double x, double y, double u1, double u2, double o) {
	return (exp(-1.0*((pow(x-u1, 2) + pow(y-u2, 2)) / (2 * pow(o, 2)))) / ((2 * PI)*pow(o, 2)));
}

double GaussianND(std::vector<double> x, std::vector<double> o) {
	if (x.size() != o.size()) { return 0; }
	double val = 1;
	for (int i = 0; i < x.size(); ++i) {
		val *= Gaussian1D(x[i], o[i]);
	}
	return val;
}

Polynomial linearLeastSquaresRegression(Matrix A) {   //returns linear polynomial function
	if (A.columns == 2 && A.rows > 0) { A = A.transpose(); }
	if (A.rows != 2) { return Polynomial(); }//only accepts an input 2 x n matrix: one row for x values, the other for y values
	double sum_x = arrsum(A.row(0));
	double sum_y = arrsum(A.row(1));
	double x_mean = sum_x / A.columns;
	double y_mean = sum_y / A.columns;
	double sum_x_sqr = pow(arrNorm(A.row(0)), 2);
	double sum_xy = dotProduct(A.row(0), A.row(1));

	double m = (sum_xy - ((sum_x*sum_y) / A.columns)) / (sum_x_sqr - (sum_x*sum_x) / A.columns);
	double b = y_mean - (x_mean*m);

	std::vector<double> coef;
	coef.push_back(b);
	coef.push_back(m);
	std::vector<double> expo;
	expo.push_back(0);
	expo.push_back(1);

	Polynomial Regression(coef);
	return Regression;
}

Polynomial linearLeastSquaresRegression(Vector v1, Vector v2) {
	if (v1.size() != v2.size()) { return Polynomial(); }//to have a least squares regression, these vectors must have the same number of values in them

	std::vector<double> vnew = v1.vec;
	for (int i = 0; i < v2.size(); ++i) {
		vnew.push_back(v2.vec[i]);
	}
	Vector v3(vnew);
	Matrix M = v3.toMatrix();
	return linearLeastSquaresRegression(M);
}

double covariance(std::vector<double> v1, std::vector<double> v2) {//sample covariance
	int sz = v1.size();
	double mu1 = arrsum(v1) / sz;//mean v1
	double mu2 = arrsum(v2) / sz;//mean v2
	double answer = 0;
	for (int i = 0; i < sz; ++i) {
		answer += (v1[i] - mu1) * (v2[i] - mu2);
	}
	return answer / (sz - 1);
}

double covariance(double* v1, double* v2, int sz) {//sample covariance
	double mu1 = arrsum(v1, sz) / sz;//mean v1
	double mu2 = arrsum(v2, sz) / sz;//mean v2
	double answer = 0;
	for (int i = 0; i < sz; ++i) {
		answer += (v1[i] - mu1) * (v2[i] - mu2);
	}
	return answer / (sz - 1);
}

double zscore(double b, double mean, double stddev) { return ((b - mean) / stddev); }

Vector zscores(Vector v) {
	double mean = v.mean();
	double stddev = v.sampleStandardDeviation();
	for (int i = 0; i < v.size(); ++i) {
		v.vec[i] = (v.vec[i] - mean) / stddev;
	}
	return v;
}

double normalCDF(double z) {
	Function f(L"f(x) = exp(-1*0.5*x*x)");
	if (z == 0) { return 0.5; }
	if (z < 0) { z = std::abs(z); return 0.5 - ((1.0 / (sqrt(2 * PI))) * integrateGaussLegendreQuadrature(0, z, f)); }
	return ((1.0 / (sqrt(2 * PI))) * integrateGaussLegendreQuadrature(0, z, f)) + 0.5;
}

double normalCDFApproximation(double z) {
	if (z == 0) { return 0.5; }
	else { return (pow((exp((-358 * z / 23) + (111 * atan((37 * z) / 294))) + 1), -1)); }//Vasquez-Leal normal cdf approximation
}

double normalPDF(double z) { return ((1 / (2 * PI))*exp(-0.5*(z*z))); }

double normalCDFExact(double x)
{//this function will calculate probability using an integral of the normal distribution function
 // now, we integrate
	if (x = 0) { return 0.5; }
	bool isNeg = false;
	if (x < 0) { isNeg = true; }
	x = std::abs(x);
	double result = 0;
	double partitionwidth = .00001;
	double f = 0;
	for (double b = 0; b < x; b += partitionwidth) {
		f = normalPDF(b);
		result += f;
	}
	result *= partitionwidth;
	if (isNeg == true) { return 1 - result; }
	return result + 0.5;
}

Vector pvalues(Vector v) {
	v = zscores(v);
	/*for (int i = 0; i < v.size(); ++i) {
	v.vec[i] = normalCDFApproximation(v.vec[i]);
	}*/
	double ans = 0;
	double n = v.vec.size();
	for (int i = 0; i < n; ++i) {
		double part1 = ((2 * (i + 1)) - 1) / n;
		v.vec[i] += ((log(normalCDFApproximation(v.vec[i]))) + (log(1 - normalCDFApproximation(v.vec[n - 1 - i]))))*(part1);
		ans += v.vec[i];
	}
	ans *= (1.0 + (0.75 / n) + (2.25 / pow(n, 2)));//adjustment for unknown pop/variance

												   //calculate p-value from AD statistic
	double p = 0;
	if (ans >= 0.6) { p = exp(1.2937 - (5.709*ans) + (0.0186*pow(ans, 2))); }
	if (0.34 < ans < 0.6) { p = exp(0.9177 - (4.279*ans) - (1.38*ans*ans)); }
	if (0.2 < ans <= 0.34) { p = 1 - exp(-8.318 + (42.796*ans) - (59.938*ans*ans)); }
	if (ans <= 0.2) { p = 1 - exp(-13.436 + (101.14*ans) - (223.73*ans*ans)); }

	v.vec.push_back(ans);
	v.vec.push_back(p);
	return v;
}

double t_score(double b, double mean, double stddev) {
	double result = ((b - mean) / (stddev / (sqrt(b))));
	return result;
}

double errorFunction(double z) {
	Function f(L"f(x) = exp(-1*x*x)");
	return (2.0 / (sqrt(PI))) * integrateGaussLegendreQuadrature(0, z, f);
}

double errorFunctionComplement(double z) { return 1 - errorFunction(z); }

double inverseErfPadeApproximation(double z) {//reasonably close approximation
	std::vector<double> numCoef;
	std::vector<double> denomCoef;

	numCoef.push_back(1);
	numCoef.push_back(0);
	numCoef.push_back(-0.25340018441678192715537113877363*PI);
	numCoef.push_back(0);
	numCoef.push_back(0.00765295341280818459241695756219*PI*PI);

	denomCoef.push_back(1);
	denomCoef.push_back(0);
	denomCoef.push_back(-0.33673351775011526048870447210696*PI);
	denomCoef.push_back(0);
	denomCoef.push_back(0.0211307465586511229664756635711*PI*PI);

	Polynomial num(numCoef);
	Polynomial denom(denomCoef);

	return (z*0.88622692545275801364908374167057*(num.HornerEvaluate(z) / denom.HornerEvaluate(z)));
}

double inverseErfc(double z) {
	if (z < 0.1) {
		double temp = log(2.0 / (PI*z*z));
		return sqrt(0.5* (temp = log(temp)));
	}
	return 1 - inverseErfPadeApproximation(z);
}

double CauchyPDF(double x, double x0, double gamma) { return  (1 / (PI*gamma)) * (pow(gamma, 2) / (pow((x - x0), 2) + pow(gamma, 2))); }
double CauchyCDF(double x, double x0, double gamma) { return (0.5 + (atan(x - x0) / gamma*PI)); }

double studentsCDF(double t, double n) {
	if (t == 0) { return 0.5; }
	std::wstring n1 = to_stringPrecision(1.0 / n);
	std::wstring n2 = to_stringPrecision(-0.5*(n + 1));
	std::wstring funct = L"f(x) = (1+x*x*";
	funct.append(n1);
	funct += L")^";
	funct += n2;
	Function f(funct);

	if (t < 0) { t = std::abs(t);  return 0.5 - ((1.0 / (sqrt(n) * betaFunction(0.5, 0.5*n))) * integrateGaussLegendreQuadrature(0, t, f)); }
	return 	0.5 + ((1.0 / (sqrt(n) * betaFunction(0.5, 0.5*n))) * integrateGaussLegendreQuadrature(0, t, f));
}

double studentsCDFApproximation(double x, double n) {//works best from x = 0 to 3
	double a2[13] = {
		-4.633640051820906435864,
		1.237623470148502047294,
		-0.01726917489890150703147,
		0.6666490496381165531403,
		2.586906931031663692266,
		1.323688973798827772654,
		1.326852402001201136628,
		1.093394658115632411821,
		0.5884542805947355903484,
		0.1418372560478959321095,
		-0.5571843343020567163038,
		1.915990707020092109758,
		0.5075164021629131605451
	};
	return (1 / (exp((a2[0] * pow(x, a2[1])) + ((a2[2] - ((n - 2)*a2[12] / n)) *pow(x, a2[3])) + ((a2[4] - ((n - 2)*a2[12] / n)) *pow(x, a2[5])) + (a2[6] * pow(x, a2[7])) + (a2[8] * tanh(x*a2[9])) + ((a2[10] / (n*a2[12]))*atan(x*a2[11]))) + 1));
}

double studentsCDFBetaVersion(double t, double v) {
	if (t == 0) { return 0.5; }
	double x = (v / ((t*t) + v));
	return (1 - (0.5*incompleteBetaFunction(x, 0.5*v, 0.5)));
}

double studentsPDF(double x, double v) {
	if (v > 30) { return normalPDF(x); }
	double temp = (v + 1)*0.5;
	return ((tgamma(temp) / (tgamma(v*0.5)*sqrt(v*PI))) * (pow(v, temp) / pow((x*x) + v, temp)));
}

double studentsPDFApproximation(double x, double n) { //a not-quite-close-enough approximation of the students PDF
	double v = pow(0.94497, -n);
	double v2 = pow(1.022, -n);
	return (exp(pow(sech(1.52938*v*x), v2*0.47494) - pow(atan(0.9990151*v*x), v2*5.3973) + pow(tanh(0.848473*v*x), v2*7.104726) + pow(tanh(0.4581817*v*x), v2*12.750694) - 1));
}

double chiSquareCDF(double k, double x) { return regularizedGamma(0.5*k, 0.5*x); }

double chiSquareCDFAppoximation(int n, double X) {
	double a1 = -0.9911;
	double a2 = -0.6763;
	double b1 = 0.8055;
	double b2 = -1.2451;
	double z = sqrt(X) - sqrt(n);
	double answer = 0;
	if (z <= 0) { answer = 1 - (0.5 * (exp((b1*z) + (a1*z*z)))); }
	else { answer = 0.5 * exp((b2*z) + (a2*z*z)); }
	return answer;
}
double chiSquareCDFInverseApproximation(double p, int n) {
	if (p > 1) { return 0; }
	double a1 = -0.9911;
	double a2 = -0.6763;
	double b1 = 0.8055;
	double b2 = -1.2451;
	double c1 = -log(2 * (1 - p));
	double c2 = -log(2 * p);
	double z = (-b1 + sqrt((b1*b1) - (4 * a1*c1))) / (2 * a1);
	if (p < 0.5) { z = ((-b2 - sqrt((b2*b2) - (4 * a2*c2))) / (2 * a2)); }
	return  pow(z + sqrt(n), 2);
}

double ChiSquareTest(Matrix A) { return chiSquareCDF(A.ChiSquareDegreesFreedom(), A.ChiSquareTestStatistic()); }

double lognormalDistributionPDF(double x, double mean, double variance) {//for non-lognormally distributed data
																		 //generate log-normal mean and variance
	double mu = log(mean / sqrt(1 + (variance / (mean*mean))));
	double vari = log(1 + (variance / (mean*mean)));

	return (1.0 / (x*sqrt(vari * 2 * PI)))*exp((-0.5 * pow(log(x - mu), 2)) / vari);
}

double lognormalDistributionCDF(double x, double mean, double variance) {//for non-lognormally distributed data
																		 //generate log-normal mean and variance
	double mu = log(mean / sqrt(1 + (variance / (mean*mean))));
	double vari = log(1 + (variance / (mean*mean)));

	return 0.5*(1 + errorFunction(log(x - mu) / sqrt(2 * vari)));
}

double FDistributionApproximation(double x, double v1, double v2) {//x is test statistic, v1/v2 are degrees freedom
	double part1 = pow((((((2 * v2) + ((v1*x) / 3) + v1 - 2)*x) / ((2 * v2) + ((4 * v1*x) / 3)))), 0.333333333333333333333333333333333333333);
	double part2 = sqrt(2 / (9 * v1));
	double part3 = (1 - (2 / (9 * v1)));
	return normalCDF((part1 - part3) / part2);
}

double correlationCoefficient(Vector a, Vector b) {
	if (a.size() != b.size()) { return 0; }
	double n = a.size();
	double mu1 = a.mean();
	double mu2 = b.mean();
	double s1 = a.sampleStandardDeviation();
	double s2 = b.sampleStandardDeviation();
	double num = 0;
	double denom1 = 0;
	double denom2 = 0;
	for (int i = 0; i < n; ++i) {
		num += (a.vec[i] - mu1) * (b.vec[i] - mu2);
		denom1 += pow((a.vec[i] - mu1), 2);
		denom2 += pow((b.vec[i] - mu2), 2);
	}
	denom1 = sqrt(denom1);
	denom2 = sqrt(denom2);
	return num / (denom1*denom2);
}

double correlationCoefficient(double* n, double* p, int vals) {
	Vector v1(n, vals);
	Vector v2(p, vals);
	return correlationCoefficient(v1, v2);
}

double coefficientOfDetermination(Vector x, Vector y) {
	return pow(correlationCoefficient(x, y), 2);
}

Matrix correlationMatrix(Matrix A) {
	//this correlation matrix is the symmetric matrix of values for pearson correlations between variable i and variable j
	Matrix M(A.columns, A.columns);
	M.identity();
	M.matrixPrecision = 4;
	for (int i = 1; i < M.columns; ++i) {
		for (int j = 0; j < i; ++j) {
			if (i != j) {
				double temp = correlationCoefficient(A.column(i), A.column(j));
				M.set(i, j, temp);
				M.set(j, i, temp);
			}
		}
	}
	return M;
}

Matrix covarianceMatrix(Matrix M) {//returns a matrix of sample covariances between column vectors
	Matrix A(M.columns, M.columns);
	for (int i = 0; i < M.columns; ++i) {
		for (int j = 0; j < M.columns; ++j) {
			A.set(i, j, covariance(M.column(i), M.column(j)));
		}
	}
	return A;
}

std::vector<double> JarqueBeraTest(Matrix M, double criticalValue) {//tests normality
	std::vector<double> answer;
	int n = M.rows;
	//double Chi2 = chiSquareCDFInverseApproximation(criticalValue, 2);//Get X^2 value for critical region bound

	for (int i = 0; i < M.columns; ++i) {
		double skewness = 0;
		double kurtosis = 0;
		std::vector<double> x = M.column(i);
		double mean = arrmean(x);
		double s = arrSampleStandardDeviation(x);

		for (int j = 0; j < n; ++j) {
			skewness += pow(x[j] - mean, 3) / (n*pow(s, 3));
			kurtosis += pow(x[j] - mean, 4) / (n*pow(s, 4));
		}
		kurtosis -= 3;
		double JBtestStat = (n * ((pow(skewness, 2) / 6) + (pow(kurtosis, 2) / 24)));
		answer.push_back(chiSquareCDF(2, JBtestStat));
	}
	return answer;
}

std::vector<double> JarqueBeraSkewness(Matrix M) {//returns JB skewness vector
	std::vector<double> answer;
	int n = M.rows;

	for (int i = 0; i < M.columns; ++i) {
		double skewness = 0;
		std::vector<double> x = M.column(i);
		double mean = arrmean(x);
		double s = arrSampleStandardDeviation(x);

		for (int j = 0; j < n; ++j) { skewness += pow(x[j] - mean, 3) / (n*pow(s, 3)); }
		answer.push_back(skewness);
	}
	return answer;
}

std::vector<double> JarqueBeraKurtosis(Matrix M) {//returns JB kurtosis vector
	std::vector<double> answer;
	int n = M.rows;

	for (int i = 0; i < M.columns; ++i) {
		double kurtosis = 0;
		std::vector<double> x = M.column(i);
		double mean = arrmean(x);
		double s = arrSampleStandardDeviation(x);

		for (int j = 0; j < n; ++j) {
			kurtosis += pow(x[j] - mean, 4) / (n*pow(s, 4));
		}
		kurtosis -= 3;
		answer.push_back(kurtosis);
	}
	return answer;
}

std::vector<double> DAgostinoTest(Matrix M, double criticalValue) {//further refines JarqueBeraTest of normality
	std::vector<double> answer;
	double n = M.rows;

	std::vector<double> sk = JarqueBeraSkewness(M);
	std::vector<double> kt = JarqueBeraKurtosis(M);
	double skewDiv = ((6 * n*(n - 1)) / ((n - 2)*(n + 1)*(n + 3)));
	double kurtDiv = ((6 * n) / (((n - 2)*(n - 3)*(n + 3)*(n + 5))));
	skewDiv = sqrt(skewDiv);
	kurtDiv = sqrt(kurtDiv) * 2 * (n - 1);

	for (int i = 0; i < sk.size(); ++i) {
		double JBTestStat = pow(sk[i] / skewDiv, 2) + pow(kt[i] / kurtDiv, 2);
		answer.push_back(chiSquareCDF(2, JBTestStat));
	}
	return answer;
}

std::vector<double> DAgostinoSkewness(Matrix M) {//returns D'Agostino skewness vector
	std::vector<double> answer;
	double n = M.rows;

	std::vector<double> sk = JarqueBeraSkewness(M);
	double skewDiv = ((6 * n*(n - 1)) / ((n - 2)*(n + 1)*(n + 3)));
	skewDiv = sqrt(skewDiv);

	for (int i = 0; i < sk.size(); ++i) {
		double JBTestStat = sk[i] / skewDiv;
		answer.push_back(JBTestStat);
	}
	return answer;
}

std::vector<double> DAgostinoKurtosis(Matrix M) {//returns D'Agostino kurtosis vector
	std::vector<double> answer;
	double n = M.rows;

	std::vector<double> kt = JarqueBeraKurtosis(M);
	double kurtDiv = ((6 * n) / (((n - 2)*(n - 3)*(n + 3)*(n + 5))));
	kurtDiv = sqrt(kurtDiv) * 2 * (n - 1);

	for (int i = 0; i < kt.size(); ++i) {
		double JBTestStat = kt[i] / kurtDiv;
		answer.push_back(JBTestStat);
	}
	return answer;
}

std::vector<double> AndersonDarlingTest(Matrix M) {
	std::vector<double> answer;
	double n = M.rows;

	for (int i = 0; i < M.columns; ++i) {
		std::vector<double> x = M.column(i);//data vector
		sort(x);//sort data greatest-to-least
		x = reverse(x);
		double mean = arrmean(x);//sample mean
		double s = arraySampleStandardDeviation(x);//sample standard deviation

		double temp = 0;
		for (int j = 1; j <= M.rows; ++j) {//calculate S
			double z = (x[j - 1] - mean) / s;
			double Yi = normalCDFApproximation(z);
			double part1 = (2 * j) - 1;
			double part2 = (2 * (n - j)) + 1;
			temp += (part1*log(Yi)) + (part2*log(1 - Yi));
		}

		//calculate Anderson-Darling test statistic (AD)
		double AD = -n - (temp / n);
		AD *= (1.0 + (0.75 / n) + (2.25 / pow(n, 2)));//adjustment for unknown pop/variance

													  //calculate p-value from AD statistic
		double p = 0;
		if (AD >= 0.6) { p = exp(1.2937 - (5.709*AD) + (0.0186*pow(AD, 2))); }
		if (0.34 < AD < 0.6) { p = exp(0.9177 - (4.279*AD) - (1.38*AD*AD)); }
		if (0.2 < AD <= 0.34) { p = 1 - exp(-8.318 + (42.796*AD) - (59.938*AD*AD)); }
		if (AD <= 0.2) { p = 1 - exp(-13.436 + (101.14*AD) - (223.73*AD*AD)); }

		//answer.push_back(p);
		answer.push_back(AD);
		x.clear();
	}
	return answer;
}

std::vector<double> AndersonDarlingTest(Vector X) {
	std::vector<double> answer;

	std::vector<double> x = X.vec;//data vector
	sort(x);//sort data greatest-to-least
	x = reverse(x);
	double mean = arrmean(x);//sample mean
	double s = arraySampleStandardDeviation(x);//sample standard deviation
	double n = x.size();

	double temp = 0;
	for (int j = 1; j <= x.size(); ++j) {//calculate S
		double z = (x[j - 1] - mean) / s;
		double Yi = normalCDFApproximation(z);
		double part1 = (2 * j) - 1;
		double part2 = (2 * (n - j)) + 1;
		temp += (part1*log(Yi)) + (part2*log(1 - Yi));
	}

	//calculate Anderson-Darling test statistic (AD)
	double AD = -n - (temp / n);
	AD *= (1.0 + (0.75 / n) + (2.25 / pow(n, 2)));//adjustment for unknown pop/variance

												  //calculate p-value from AD statistic
	double p = 0;
	if (AD >= 0.6) { p = exp(1.2937 - (5.709*AD) + (0.0186*pow(AD, 2))); }
	if (0.34 < AD < 0.6) { p = exp(0.9177 - (4.279*AD) - (1.38*AD*AD)); }
	if (0.2 < AD <= 0.34) { p = 1 - exp(-8.318 + (42.796*AD) - (59.938*AD*AD)); }
	if (AD <= 0.2) { p = 1 - exp(-13.436 + (101.14*AD) - (223.73*AD*AD)); }

	//answer.push_back(p);
	answer.push_back(AD);
	x.clear();
	return answer;
}

double GompertzPDF(double x, double n, double b) {
	if (x < 0 || n <= 0 || b <= 0) { return 0; }
	return b*n*exp(n)*exp(b*x)*exp(-n*exp(b*x));
}

double GompertzCDF(double x, double n, double b) {
	if (x < 0 || n <= 0 || b <= 0) { return 0; }
	return 1 - exp(-n*(exp(b*x)-1));
}

double rootMeanSquare(std::vector<double> vec) {
	double answer = 0;
	for (int i = 0; i < vec.size(); ++i) {
		answer += pow(vec[i], 2);
	}
	answer *= 1.0 / (vec.size());
	return pow(answer,0.5);
}

ComplexNumber rootMeanSquare(std::vector<ComplexNumber> vec) {
	std::complex<double> answer = 0;
	for (int i = 0; i < vec.size(); ++i) {
		answer += pow(vec[i], 2).toComplex();
	}
	answer *= (1.0 / vec.size());
	ComplexNumber answer2(answer);
	return pow(answer2, 0.5);
}

double rootMeanSquare(int T1, int T2, Function f) {
	//format function of f(t) so that it is now abs(f(t)^2)
	std::wstring funct = L"(abs(";
	funct.append(f.function);
	funct.append(L")^2)");
	f.function = funct;
	return sqrt(integrateGaussLegendreQuadrature(T1,T2,f) / (T2 - T1));
}

ComplexNumber rootMeanSquare(int T1, int T2, ComplexFunction f) {
	
	// ADD SOMETHING LATER

	return 0;
}