#pragma once
#include "stdafx.h"

GeometricObject::GeometricObject(){}
GeometricObject::~GeometricObject(){}
//===================================
Curve::Curve() {}
Curve::~Curve() {}
Curve::Curve(std::vector<Function> f) {
	dimensions = f.size();
	functs = f;
}
//===================================
Line::Line(){}
Line::~Line() { }
Line::Line(Function f) {//case: line has two dimensions and is epitomized by a single equation
	p1.push_back(0);
	p1.push_back(f.evaluate(0));
	p2.push_back(1);
	p2.push_back(f.evaluate(1));
	dimensions = 2;
	objectType = L"Line";
}

Line::Line(double p1x, double p1y, double p2x, double p2y) {
	std::vector<double> pt1, pt2;
	pt1.push_back(p1x);
	pt1.push_back(p1y);
	pt2.push_back(p2x);
	pt2.push_back(p2y);
	p1 = pt1;
	p2 = pt2;
	dimensions = 2;
	objectType = L"Line";
}

Line::Line(std::vector<double> pt1, std::vector<double> pt2) {//case: n dimensional line
	p1 = pt1; p2 = pt2;
	while (p1.size() > p2.size()) { p2.push_back(0); }//make sure dimensions of points match
	while (p1.size() < p2.size()) { p1.push_back(0); }
	dimensions = pt1.size();//line equations are of the form a_0*t + a_1*t + ... + a_n = 0
	objectType = L"Line";
}

Matrix Line::toMatrix() {
	std::vector<double> coef = getSlope();
	Matrix A(p1.size(),2);
	for (int i = 0; i < p1.size(); ++i) {

	}
	return A;
}

VectorValuedFunction Line::toFunction() {
	std::vector<double> coef = getSlope();
	while (coef.size() != p1.size()) { coef.push_back(0); }
	std::vector<Function> functs;

	for (int i = 0; i < p1.size(); ++i) {
		std::wstring f = L"f(t) = ";
		f.append(to_stringPrecision(coef[i]));
		f.append(L"t ");
		if (p1[i] > 0) { f.append(L"+ "); }
		f.append(to_stringPrecision(p1[i]));
		functs.push_back(Function(f));
	}
	return VectorValuedFunction(functs);
}

std::wstring Line::toString() {
	VectorValuedFunction vvf = toFunction();
	std::wstring str = L"F(";
	str.append((vvf.variablesToString()));
	str.append(L") = ");
	str.append(vvf.toString());
	return str;
}

std::vector<double> Line::getSlope() {//gets vector of values P2-P1, such that P2 = P1 + (P2-P1)*t where t is a parameter
	std::vector<double> vec;
	for (int i = 0; i < p1.size(); ++i) {
		vec.push_back(p2[i]-p1[i]);
	}
	return vec;
}

std::vector<double> getIntersection(Line l1, Line l2) {
	std::vector<double> ans;
	

	return ans;
}

bool isPointOnLine(Line l, std::vector<double> point) {
	while (point.size() < l.dimensions) { point.push_back(0); }

	return false;
}

bool doLinesIntersect(Line l1, Line l2) {
	
	return false;
}
bool areParallel(Line l1, Line l2) {
	std::vector<double> v1, v2;
	v1 = l1.getSlope();
	v2 = l2.getSlope();
	if (v1 == v2 ||
		arrmultiply(v1,-1) == v2 ||
		v1 == arrmultiply(v2,-1)) 
	{ return true; }
	return false;
}
bool doCoincide(Line l1, Line l2) {
	
	return false;
}
//===================================
LineSegment::LineSegment() {}
LineSegment::~LineSegment() {}
LineSegment::LineSegment(double p1x, double p1y, double p2x, double p2y) {
	std::vector<double> pt1, pt2;
	pt1.push_back(p1x);
	pt1.push_back(p1y);
	pt2.push_back(p2x);
	pt2.push_back(p2y);
	p1 = pt1;
	p2 = pt2;
	end1 = p1;
	end2 = p2;
	length = EuclideanDistance(p1, p2);
	dimensions = 1;
	objectType = L"Line";
}
LineSegment::LineSegment(std::vector<double> endp1, std::vector<double> endp2) {
	dimensions = 1;
	p1 = endp1;
	p2 = endp2;
	end1 = endp1;
	end2 = endp2;
	length = EuclideanDistance(endp1,endp2);
	objectType = L"LineSegment";
}
//===================================
Triangle::Triangle() {}
Triangle::~Triangle(){}
Triangle::Triangle(LineSegment ls1, LineSegment ls2, LineSegment ls3) {
	line1 = ls1; line2 = ls2; line3 = ls3;
	area = getArea();
	circumference = getCircumference();
	GeometricObject::objectType = L"Triangle";
}
double Triangle::getArea() {return 0.5*line1.length*line2.length*sin(line3.length);}
double Triangle::getCircumference() { return line1.length + line2.length + line3.length; }
bool Triangle::isEquilateral() {
	if ((line1.length = line2.length) && (line2.length == line3.length)) { return true; }
	return false;
}
bool Triangle::isIsosceles() {
	
	return false;
}
//===================================
Sphere::Sphere() {}
Sphere::~Sphere() {};
Sphere::Sphere(int n) { dimensions = n; offcenter = false; radius = 1; }//unit sphere at center of basis declaration
Sphere::Sphere(int n, std::vector<double> x) { dimensions = n; offcenter = true; center = x; radius = 1; }//unit sphere centered off center declaration
Sphere::Sphere(int n, double r) { dimensions = n; offcenter = false; radius = r; }//centered non-unit sphere
Sphere::Sphere(int n, std::vector<double> x, double r) { dimensions = n; offcenter = true; center = x; radius = r; }//offcenter non-unit sphere
double Sphere::getVolume() { return (pow(radius, dimensions)*pow(PI, dimensions*0.5)) / tgamma(0.5*dimensions + 1); }
double Sphere::getCircumference(int p) { return getVolume() / pow(radius, dimensions); }
double Sphere::getSurfaceArea() {return dimensions * getCircumference(dimensions)*pow(radius, dimensions - 1);}

double Sphere::function(std::vector<double> x) {
	double xsumsq = arrsumexp(x, 2);
	if (radius*radius < xsumsq) { return INFINITY; }
	else { return sqrt(radius*radius - xsumsq); }
}
//==================================

Hyperplane::Hyperplane() {}
Hyperplane::~Hyperplane() {}

Hyperplane::Hyperplane(std::vector<double> p) {
	dimensions = p.size() - 1;
	vec = p;
}

bool Hyperplane::isCoplanar(std::vector<double> p1) {
	double sum = 0;
	for (int i = 0; i < p1.size(); ++i) {
		sum += p1[i] * vec[i];
	}
	if (abs(sum) < SYSTEM_EPSILON) { return true; }
	return false;
}