#pragma once
#include "stdafx.h"

double EuclideanDistance(double p1, double p2) {
	return sqrt(pow(p2-p1,2));
}

double EuclideanDistance(std::vector<double> p1, std::vector<double> p2) {
	while (p1.size() < p2.size()) { p1.push_back(0); }
	while (p1.size() > p2.size()) { p2.push_back(0); }
	double ans = 0;
	for (int i = 0; i < p1.size(); ++i) {
		ans += pow(p2[i]-p1[i], 2);
	}
	return sqrt(ans);
}

double EuclideanDistance(ComplexNumber p1, ComplexNumber p2) {
	return sqrt( pow(p2.Re-p1.Re,2) + pow(p2.Im-p1.Im,2) );
}

double EuclideanDistance(std::vector<ComplexNumber> p1, std::vector<ComplexNumber> p2) {
	while (p1.size() < p2.size()) { p1.push_back(ComplexNumber(0)); }
	while (p1.size() > p2.size()) { p2.push_back(ComplexNumber(0)); }
	double ans = 0;
	for (int i = 0; i < p1.size(); ++i) {
		ans += pow(p2[i].Re-p1[i].Re,2);
		ans += pow(p2[i].Im-p1[i].Im,2);
	}
	return sqrt(ans);
}