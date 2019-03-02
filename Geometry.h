#ifndef GEOMETRY_H
#define GEOMETRY_H
#pragma once
#include "stdafx.h"

class GeometricObject {
public:
	GeometricObject();
	~GeometricObject();
	int dimensions;
	std::wstring objectType;
};
class Curve : public GeometricObject {//an n-dimensional curve is a function defined by n parametric functions
public:
	std::vector<Function> functs;
	Curve();
	~Curve();
	Curve(std::vector<Function> f);
	std::wstring toString();
};
class Line : public GeometricObject {//a continuous subset of the reals between two points in N dimensions, parameterized by the variable 't'
public:
	std::vector<double> p1, p2;//points that define the line
	Line();
	~Line();
	Line(Function f);
	Line(double p1x, double p1y, double p2x, double p2y);
	Line(std::vector<double> pt1, std::vector<double> pt2);
	std::vector<double> getSlope();
	VectorValuedFunction toFunction();
	Matrix toMatrix();
	std::wstring toString();
};
std::vector<double> getIntersection(Line l1, Line l2);
bool doLinesIntersect(Line l1, Line l2);
bool areParallel(Line l1, Line l2);
bool doCoincide(Line l1, Line l2);

class LineSegment : public Line{//a line with two endpoints
public:	
	std::vector<double> end1;
	std::vector<double> end2;
	double length;
	LineSegment();
	LineSegment(double p1x, double p1y, double p2x, double p2y);
	LineSegment(std::vector<double> endp1, std::vector<double> endp2);
	~LineSegment();
	std::wstring toString();
};
class Triangle : public GeometricObject, public LineSegment {//defined as the intersection of 3 line segments
public:
	LineSegment line1, line2, line3;
	double area;
	double circumference;
	Triangle();
	Triangle(LineSegment ls1, LineSegment ls2, LineSegment ls3);
	~Triangle();
	double getArea();
	double getCircumference();
	bool isEquilateral();
	bool isIsosceles();
	std::wstring toString();
};
class Sphere : public GeometricObject {
public:
	std::vector<double> center;
	double radius;
	bool offcenter;
	Sphere();
	~Sphere();
	Sphere(int n);
	Sphere(int n, std::vector<double> x);
	Sphere(int n, double r);
	Sphere(int n, std::vector<double> x, double r);
	double getSurfaceArea();
	double getCircumference(int p);
	double getVolume();
	double function(std::vector<double> x);
	std::wstring toString();
};
class Circle : public Sphere {

	std::wstring toString();
};
class Hyperplane : public GeometricObject {//Hyperplanes are defined by a single linear equation, a1x1 + a2x2+...+ anxn - b = 0
public:
	std::vector<double> vec;
	Hyperplane();
	~Hyperplane();
	Hyperplane(std::vector<double> p1);
	bool isCoplanar(std::vector<double> p1);
};
#endif