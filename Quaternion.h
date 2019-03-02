#ifndef QUATERNION_H
#define QUATERNION_H
#pragma once
//#include "stdafx.h"

class Quaternion { 	// conventional form:  q = w+xi+yj+zk  | openGL quaternion:  q = xi+yj+zk+w
public: 
	std::vector<double> element;
	unsigned int QuaternionPrecision;
	Quaternion();
	Quaternion(std::vector<double> v);
	Quaternion(double x, double y, double z, double w);
	~Quaternion();
	void clear();
	Matrix CreateMatrix();
	void CreateFromAxisAngle(double x, double y, double z, double deg);
	void randomize();
	void randomizeBoolean();
	void randomizeInteger();
	Quaternion conjugate();
	double norm();
	double pNorm(double p);
	void normalize();	
	double* toArray();
	std::vector<double> toVector();
	double* toEulerAngles();
	std::vector<double> toEulerAnglesVector();
	Quaternion HamiltonProduct(Quaternion q, Quaternion p);
	Quaternion add(Quaternion q, double n);
	Quaternion add(Quaternion q, Quaternion p);
	Quaternion subtract(Quaternion q, double n);
	Quaternion subtract(Quaternion q, Quaternion p);
	Quaternion multiply(Quaternion q, double n);
	Quaternion multiply(Quaternion q, Quaternion p);
	Quaternion inverse();
	Quaternion operator*=(Quaternion rhs);
	Quaternion operator*=(double x);
	Quaternion operator*(Quaternion rhs);
	Quaternion operator*(double x);
	Quaternion operator+=(Quaternion rhs);
	Quaternion operator+=(double x);
	Quaternion operator+(Quaternion rhs);
	Quaternion operator+(double x);
	Quaternion operator-=(Quaternion rhs);
	Quaternion operator-=(double x);
	Quaternion operator-(Quaternion rhs);
	Quaternion operator-(double x);
	Quaternion convertEulerAnglesToQuaternion(double* v);
	Quaternion convertEulerAnglesToQuaternion(std::vector<double> v);
	Quaternion createRotationQuaternion(double* v);
	Quaternion createRotationQuaternion(std::vector<double> v);
	Quaternion rotate(Quaternion p, Quaternion q);
	Matrix toMatrix();
	Matrix toMatrix(int r, int c);
	Matrix toRotationMatrix();
	Matrix toScalingMatrix();
	Matrix toTranslationMatrix();
	void toOpenGLForm();
	std::wstring toString();
	void display();
};
//END QUATERNION CLASS=========================================

class DualQuaternion {
	/*	A quaternion is of the form q= w+xi+yj+zk,
	while a dual quaternion consists of two quaternions,
	Q = r + de, where r and d are the real and dual-valued
	quaternions.	*/
public:
	Quaternion R;
	Quaternion D;
	unsigned int DualQuaternionPrecision;
	DualQuaternion();
	~DualQuaternion();
	DualQuaternion(Quaternion r, Quaternion d);
	DualQuaternion(Quaternion r);
	DualQuaternion(std::vector<double> v, std::vector<double> w);
	DualQuaternion(double x, double y, double z, double w, double x2, double y2, double z2, double w2);
	void clear();
	DualQuaternion HamiltonProduct(DualQuaternion q, DualQuaternion p);
	DualQuaternion add(DualQuaternion q, double n);
	DualQuaternion add(DualQuaternion q, DualQuaternion p);
	DualQuaternion subtract(DualQuaternion q, double n);
	DualQuaternion subtract(DualQuaternion q, DualQuaternion p);
	DualQuaternion multiply(DualQuaternion q, double n);
	DualQuaternion multiply(DualQuaternion q, DualQuaternion p);
	DualQuaternion inverse(DualQuaternion q);
	DualQuaternion operator*=(DualQuaternion rhs);
	DualQuaternion operator*=(double x);
	DualQuaternion operator*(DualQuaternion rhs);
	DualQuaternion operator*(double x);
	DualQuaternion operator+=(DualQuaternion rhs);
	DualQuaternion operator+=(double x);
	DualQuaternion operator+(DualQuaternion rhs);
	DualQuaternion operator+(double x);
	DualQuaternion operator-=(DualQuaternion rhs);
	DualQuaternion operator-=(double x);
	DualQuaternion operator-(DualQuaternion rhs);
	DualQuaternion operator-(double x);
	std::wstring toString();
	void display();
};//========================================================

class Octonian {   // If using Octonians to represent rotations, O must be a unit Octonian.
public:
	std::vector<double> element;//element[0] is magnitude, element[n] is part of the vector
	int OctonianPrecision;
	Octonian();
	Octonian(std::vector<double> v);
	~Octonian();
	Matrix CreateMatrix();
	//void CreateFromAxisAngle(double x, double y, double z, double deg);
	Octonian conjugate();
	double norm();
	/*double* toArray();
	std::vector<double> toVector();
	double* toEulerAngles();
	Octonian add(Octonian q, double n);
	Octonian add(Octonian q, Octonian p);
	Octonian multiply(Octonian q, double n);
	Octonian inverse(Octonian q);
	Octonian HamiltonProduct(Octonian q, Octonian p);
	Octonian convertEulerAnglesToOctonian(double* v);
	Octonian rotateByOctonians(Octonian q, Octonian p);
	Octonian createRotationOctonian(double* v);*/
	std::wstring toString();
	void display(double* v);
};
//END Octonian CLASS=========================================


#endif