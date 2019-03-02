#pragma once
#include "stdafx.h"

Quaternion::Quaternion() { }
Quaternion::Quaternion(std::vector<double> v) {
		if (v.size() == 3) {
			element = v;
			element.push_back(sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2)));
			return;
		}
		if (v.size() == 4) { element = v; return; }
		for (int i = 0; i < 4; ++i) {
			if (i < v.size()) { element.push_back(v[i]); }
			else { element.push_back(0); }
		}
	}
Quaternion::Quaternion(double x, double y, double z, double w) {
		element.clear();
		element.push_back(x);
		element.push_back(y);
		element.push_back(z);
		element.push_back(w);
	}
Quaternion::~Quaternion() { element.clear(); };

	void Quaternion::clear() { element.clear(); }

	Matrix Quaternion::CreateMatrix() {
		Matrix A(4, 4);
		A.element[0] = element[0];
		A.element[4] = element[1];
		A.element[8] = element[2];
		A.element[12] = element[3];
		return A;
	}

	void Quaternion::CreateFromAxisAngle(double x, double y, double z, double deg) {
		double angle = double((deg / 180) * PI);
		double result = double(sin(angle / 2));
		element[0] = double(cos(angle / 2));
		// Calculate the x, y and z of the quaternion
		element[1] = double(x * result);
		element[2] = double(y * result);
		element[3] = double(z * result);
	}

	void Quaternion::randomize() {
		srand(time(0) + clock());
		std::vector<double> vec;
		for (int i = 0; i < 4; ++i) {
			vec.push_back(((double)rand() / 1000) * pow(-1, rand() % 2));
		}
		element = vec;
	}

	void Quaternion::randomizeBoolean() {
		srand(time(0) + clock());
		std::vector<double> vec;
		for (int i = 0; i < vec.size(); ++i) {
			vec.push_back(rand() % 2);
		}
		*this = Quaternion(vec);
	}

	void Quaternion::randomizeInteger() {
		srand(time(0) + clock());
		std::vector<double> vec;
		for (int i = 0; i < vec.size(); ++i) {
			vec.push_back((floor(rand() / 1000)) * pow(-1, rand() % 2));
		}
		*this = Quaternion(vec);
	}

	Quaternion Quaternion::conjugate() { return Quaternion(element[0], -element[1], -element[2], -element[3]); }

	double Quaternion::norm() { return sqrt((element[0] * element[0]) + (element[1] * element[1]) + (element[2] * element[2]) + (element[3] * element[3])); }

	double Quaternion::pNorm(double p) {
		double answer = 0;
		for (int i = 0; i < element.size(); ++i) {
			answer += pow(element[i], p);
		}
		return sqrt(answer);
	}


	void Quaternion::normalize() {
		double nrm = 1.0 / norm();
		element[0] *= nrm;
		element[1] *= nrm;
		element[2] *= nrm;
		element[3] *= nrm;
	}

	double* Quaternion::toArray() {
		double* v = new double[3];
		v[0] = element[1];
		v[1] = element[2];
		v[2] = element[3];
		return v;
	}

	std::vector<double> Quaternion::toVector() {
		std::vector<double> v;
		v.push_back(element[1]);
		v.push_back(element[2]);
		v.push_back(element[3]);
		return v;
	}

	double* Quaternion::toEulerAngles() {
		double v[3];
		double x = element[0];
		double y = element[1];
		double z = element[2];
		double w = element[3];
		v[0] = atan2(2 * (w*x + y*z), 1 - 2 * (x*x + y*y));
		v[1] = asin(2 * (w*y - z*x));
		v[2] = atan2(2 * (w*z + x*y), 1 - 2 * (y*y + z*z));
		return v;
	}

	std::vector<double> Quaternion::toEulerAnglesVector() {
		std::vector<double> v;
		double x = element[0];
		double y = element[1];
		double z = element[2];
		double w = element[3];
		v.push_back(atan2(2 * (w*x + y*z), 1 - 2 * (x*x + y*y)));
		v.push_back(asin(2 * (w*y - z*x)));
		v.push_back(atan2(2 * (w*z + x*y), 1 - 2 * (y*y + z*z)));
		return v;
	}

	//Quaternion arithmetic
	//=====================
	Quaternion Quaternion::HamiltonProduct(Quaternion q, Quaternion p) { //Hamilton Product 
		double w_val = q.element[0] * p.element[0] - q.element[1] * p.element[1] - q.element[2] * p.element[2] - q.element[3] * p.element[3];
		double x_val = q.element[0] * p.element[1] + q.element[1] * p.element[0] + q.element[2] * p.element[3] - q.element[3] * p.element[2];
		double y_val = q.element[0] * p.element[2] - q.element[1] * p.element[3] + q.element[2] * p.element[0] + q.element[3] * p.element[1];
		double z_val = q.element[0] * p.element[3] + q.element[1] * p.element[2] - q.element[2] * p.element[1] + q.element[3] * p.element[0];
		q.element[0] = w_val; q.element[1] = x_val; q.element[2] = y_val; q.element[3] = z_val;
		return q;
	}

	Quaternion Quaternion::add(Quaternion q, double n) { q.element[0] += n; q.element[1] += n; q.element[2] += n; q.element[3] += n; return q; }
	Quaternion Quaternion::add(Quaternion q, Quaternion p) { q.element[0] += p.element[0]; q.element[1] += p.element[1]; q.element[2] += p.element[2]; q.element[3] += p.element[3]; return q; }

	Quaternion Quaternion::subtract(Quaternion q, double n) { q.element[0] -= n; q.element[1] -= n; q.element[2] -= n; q.element[3] -= n; return q; }
	Quaternion Quaternion::subtract(Quaternion q, Quaternion p) { q.element[0] -= p.element[0]; q.element[1] -= p.element[1]; q.element[2] -= p.element[2]; q.element[3] -= p.element[3]; return q; }

	Quaternion Quaternion::multiply(Quaternion q, double n) { q.element[0] *= n; q.element[1] *= n; q.element[2] *= n; q.element[3] *= n; return q; }
	Quaternion Quaternion::multiply(Quaternion q, Quaternion p) { return HamiltonProduct(q, p); }

	Quaternion Quaternion::inverse() {
		Quaternion q_ = conjugate();
		q_ = multiply(q_, 1 / norm());
		return q_;
	}

	Quaternion Quaternion::operator*=(Quaternion rhs) { *this = multiply(*this, rhs); return*this; }
	Quaternion Quaternion::operator*=(double x) { *this = multiply(*this, x); return *this; }
	Quaternion Quaternion::operator*(Quaternion rhs) { return multiply(*this, rhs); }
	Quaternion Quaternion::operator*(double x) { return multiply(*this, x); }

	Quaternion Quaternion::operator+=(Quaternion rhs) { *this = add(*this, rhs); return*this; }
	Quaternion Quaternion::operator+=(double x) { *this = add(*this, x); return *this; }
	Quaternion Quaternion::operator+(Quaternion rhs) { return add(*this, rhs); }
	Quaternion Quaternion::operator+(double x) { return add(*this, x); }

	Quaternion Quaternion::operator-=(Quaternion rhs) { *this = subtract(*this, rhs); return*this; }
	Quaternion Quaternion::operator-=(double x) { *this = subtract(*this, x); return *this; }
	Quaternion Quaternion::operator-(Quaternion rhs) { return subtract(*this, rhs); }
	Quaternion Quaternion::operator-(double x) { return subtract(*this, x); }
	//==================================================================

	Quaternion Quaternion::convertEulerAnglesToQuaternion(double* v) {//    v = <pitch,yaw,roll> or <x-axis rotate, y rotate, z rotate>
		double q[4];

		double t0 = cos(v[2] * 0.5);
		double t1 = sin(v[2] * 0.5);
		double t2 = cos(v[1] * 0.5);
		double t3 = sin(v[1] * 0.5);
		double t4 = cos(v[0] * 0.5);
		double t5 = sin(v[0] * 0.5);

		q[0] = t0 * t2 * t4 + t1 * t3 * t5;
		q[1] = t0 * t3 * t4 - t1 * t2 * t5;
		q[2] = t0 * t2 * t5 + t1 * t3 * t4;
		q[3] = t1 * t2 * t4 - t0 * t3 * t5;

		return Quaternion(toSTLVector(q, 4));
	}

	Quaternion Quaternion::convertEulerAnglesToQuaternion(std::vector<double> v) {//   v = <pitch,yaw,roll> or <x-axis rotate, y rotate, z rotate>
		double q[4];

		double t0 = cos(v[2] * 0.5);
		double t1 = sin(v[2] * 0.5);
		double t2 = cos(v[1] * 0.5);
		double t3 = sin(v[1] * 0.5);
		double t4 = cos(v[0] * 0.5);
		double t5 = sin(v[0] * 0.5);

		q[0] = t0 * t2 * t4 + t1 * t3 * t5;
		q[1] = t0 * t3 * t4 - t1 * t2 * t5;
		q[2] = t0 * t2 * t5 + t1 * t3 * t4;
		q[3] = t1 * t2 * t4 - t0 * t3 * t5;

		return Quaternion(toSTLVector(q, 4));
	}

	Quaternion Quaternion::createRotationQuaternion(double* v) {
		/*  given a = angle to rotate
		[x, y, z] = axis to rotate around
		then rotation Quaternion R = [cos(a / 2), sin(a / 2)*x, sin(a / 2)*y, sin(a / 2)*z] */
		return Quaternion(cos(v[0] / 2), sin(v[0] / 2)*v[1], sin(v[0] / 2)*v[2], sin((v[0] / 2)*v[3]));
	}

	Quaternion Quaternion::createRotationQuaternion(std::vector<double> v) {
		/*  given a = angle to rotate
		[x, y, z] = axis to rotate around
		then rotation Quaternion R = [cos(a / 2), sin(a / 2)*x, sin(a / 2)*y, sin(a / 2)*z] */
		return Quaternion(cos(v[0] / 2), sin(v[0] / 2)*v[1], sin(v[0] / 2)*v[2], sin((v[0] / 2)*v[3]));
	}

	Quaternion Quaternion::rotate(Quaternion p, Quaternion q) { // p' = q p q'
		double nrm = norm();//save length of vector to apply later to final rotated vector
		q = createRotationQuaternion(q.element);
		Quaternion q_ = q;
		q_.element[1] *= -1;
		q_.element[2] *= -1;
		q_.element[3] *= -1;
		q = HamiltonProduct(HamiltonProduct(q, p), q_);
		if (abs(q.element[0]) < 0.00000000001) { q.element[0] = 0; }
		if (abs(q.element[1]) < 0.00000000001) { q.element[1] = 0; }
		if (abs(q.element[2]) < 0.00000000001) { q.element[2] = 0; }
		if (abs(q.element[3]) < 0.00000000001) { q.element[3] = 0; }
		return q;
	}

	Matrix Quaternion::toMatrix() {//uses homogeneous coordinates
		Matrix M(4, 1);//create column matrix
		M.set(0, 0, element[0]);
		M.set(1, 0, element[1]);
		M.set(2, 0, element[2]);
		M.set(3, 0, element[3]);
		return M;
	}

	Matrix Quaternion::toMatrix(int r, int c) {
		if (r < 4 && c < 4) { return Matrix(); }
		Matrix M(r, c);
		M.identity();
		for (int i = 0; i < element.size(); ++i) {
			M.set(i, i, element[i]);
		}
		return M;
	}

	Matrix Quaternion::toRotationMatrix() {
		if (element[1] == 1) {//x-axis rotation
			Matrix M(4, 4);
			M.identity();
			M.set(1, 1, cos(element[0]));
			M.set(1, 2, -sin(element[0]));
			M.set(2, 1, sin(element[0]));
			M.set(2, 2, cos(element[0]));
			return M;
		}
		if (element[2] == 1) {//y-axis rotation
			Matrix M(4, 4);
			M.identity();
			M.set(0, 0, cos(element[0]));
			M.set(0, 2, -sin(element[0]));
			M.set(2, 0, sin(element[0]));
			M.set(2, 2, cos(element[0]));
			return M;
		}
		if (element[3] == 1) {//z-axis rotation
			Matrix M(4, 4);
			M.identity();
			M.set(0, 0, cos(element[0]));
			M.set(0, 1, -sin(element[0]));
			M.set(1, 0, sin(element[0]));
			M.set(1, 1, cos(element[0]));
			return M;
		}
	}

	Matrix Quaternion::toScalingMatrix() {//uses homogeneous coordinates
		Matrix M(4, 4);
		M.identity();
		M.set(0, 0, element[0]);
		M.set(1, 1, element[1]);
		M.set(2, 2, element[2]);
		M.set(3, 3, element[3]);
		return M;
	}

	Matrix Quaternion::toTranslationMatrix() {//uses homogeneous coordinates
		Matrix M(4, 4);
		M.identity();
		M.set(0, 3, element[0]);
		M.set(1, 3, element[1]);
		M.set(2, 3, element[2]);
		M.set(3, 3, element[3]);
		return M;
	}

	void Quaternion::toOpenGLForm() {	//swap form from w,x,y,z to x,y,z,w
		double temp = element[0];
		element[0] = element[3];
		element[3] = element[0];
	}

	std::wstring Quaternion::toString() {
		std::wostringstream sstr;
		double w = element[0];
		double x = element[1];
		double y = element[2];
		double z = element[3];
		if (w != floor(w)) { QuaternionPrecision = precision; }
		else { QuaternionPrecision = 0; }
		sstr << L"<" << std::setprecision(QuaternionPrecision) << w << L",";
		if (x != floor(x)) { QuaternionPrecision = precision; }
		else { QuaternionPrecision = 0; }
		sstr << std::setprecision(QuaternionPrecision) << x << L",";
		if (y != floor(y)) { QuaternionPrecision = precision; }
		else { QuaternionPrecision = 0; }
		sstr << std::setprecision(QuaternionPrecision) << y << L",";
		if (z != floor(z)) { QuaternionPrecision = precision; }
		else { QuaternionPrecision = 0; }
		sstr << std::setprecision(QuaternionPrecision) << z << L">";
		std::wstring answer = sstr.str();
		sstr.clear();
		return answer;
	}

	void Quaternion::display() { std::wcout << toString() << std::endl; }
	//===================================================================

		DualQuaternion::DualQuaternion() { }
		DualQuaternion::~DualQuaternion() { R.clear(); D.clear(); };
		DualQuaternion::DualQuaternion(Quaternion r, Quaternion d) { R = r; D = d; }
		DualQuaternion::DualQuaternion(Quaternion r) { R = r; D = Quaternion(0, 0, 0, 0); }
		DualQuaternion::DualQuaternion(std::vector<double> v, std::vector<double> w) {
			std::vector<double> r;
			std::vector<double> d;
			if (v.size() < 3 && v.size() > 0) {
				for (int i = 0; i < 4; ++i) {
					if (i < v.size()) { r.push_back(v[i]); }
					else { r.push_back(0); }
				}
			}
			if (v.size() == 3) {
				r = v;
				r.push_back(sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2)));
			}
			if (v.size() == 4) { r = v; }
			//==============================
			if (w.size() < 3 && w.size() > 0) {
				for (int i = 0; i < 4; ++i) {
					if (i < w.size()) { d.push_back(w[i]); }
					else { d.push_back(0); }
				}
			}
			if (w.size() == 3) {
				d = w;
				d.push_back(sqrt(pow(w[0], 2) + pow(w[1], 2) + pow(w[2], 2)));
			}
			if (w.size() == 4) { d = w; }
			R = Quaternion(r);
			D = Quaternion(d);
		}

		DualQuaternion::DualQuaternion(double x, double y, double z, double w, double x2, double y2, double z2, double w2) {
			std::vector<double> r;
			std::vector<double> d;
			r.clear();
			r.push_back(x);
			r.push_back(y);
			r.push_back(z);
			r.push_back(w);
			//============
			d.clear();
			d.push_back(x2);
			d.push_back(y2);
			d.push_back(z2);
			d.push_back(w2);
			R = Quaternion(r);
			D = Quaternion(d);
		}

		void DualQuaternion::clear() { R.clear(); D.clear(); }

		//DualQuaternion arithmetic
		//=====================
		DualQuaternion DualQuaternion::HamiltonProduct(DualQuaternion q, DualQuaternion p) {
			//Hamilton Product of q and p, where q = a+be, p = c+de is ac+(ad+cb)e
			return DualQuaternion(q.R.HamiltonProduct(q.R, p.R), q.R.HamiltonProduct(q.R, p.D) + q.R.HamiltonProduct(p.R, q.D));
		}

		DualQuaternion DualQuaternion::add(DualQuaternion q, double n) { q.R.element[0] += n; q.R.element[1] += n; q.R.element[2] += n; q.R.element[3] += n; return q; }
		DualQuaternion DualQuaternion::add(DualQuaternion q, DualQuaternion p) {
			q.R.element[0] += p.R.element[0]; q.R.element[1] += p.R.element[1]; q.R.element[2] += p.R.element[2]; q.R.element[3] += p.R.element[3];
			q.D.element[0] += p.D.element[0]; q.D.element[1] += p.D.element[1]; q.D.element[2] += p.D.element[2]; q.D.element[3] += p.D.element[3];
			return q;
		}

		DualQuaternion DualQuaternion::subtract(DualQuaternion q, double n) {
			q.R.element[0] -= n; q.R.element[1] -= n; q.R.element[2] -= n; q.R.element[3] -= n;
			return q;
		}
		DualQuaternion DualQuaternion::subtract(DualQuaternion q, DualQuaternion p) {
			q.R.element[0] -= p.R.element[0]; q.R.element[1] -= p.R.element[1]; q.R.element[2] -= p.R.element[2]; q.R.element[3] -= p.R.element[3];
			q.D.element[0] -= p.D.element[0]; q.D.element[1] -= p.D.element[1]; q.D.element[2] -= p.D.element[2]; q.D.element[3] -= p.D.element[3];
			return q;
		}

		DualQuaternion DualQuaternion::multiply(DualQuaternion q, double n) {
			q.R.element[0] *= n; q.R.element[1] *= n; q.R.element[2] *= n; q.R.element[3] *= n;
			q.D.element[0] *= n; q.D.element[1] *= n; q.D.element[2] *= n; q.D.element[3] *= n;
			return q;
		}
		DualQuaternion DualQuaternion::multiply(DualQuaternion q, DualQuaternion p) { return HamiltonProduct(q, p); }

		DualQuaternion DualQuaternion::inverse(DualQuaternion q) {
			//For some DualQuaternion Q = p + qe, then Q^-1 = p^-1(1 - (q*p^-1)e)
			return DualQuaternion(q.R.inverse(), (q.R.inverse()*q.D*q.R.inverse())*-1);
		}

		DualQuaternion DualQuaternion::operator*=(DualQuaternion rhs) { *this = multiply(*this, rhs); return*this; }
		DualQuaternion DualQuaternion::operator*=(double x) { *this = multiply(*this, x); return *this; }
		DualQuaternion DualQuaternion::operator*(DualQuaternion rhs) { return multiply(*this, rhs); }
		DualQuaternion DualQuaternion::operator*(double x) { return multiply(*this, x); }

		DualQuaternion DualQuaternion::operator+=(DualQuaternion rhs) { *this = add(*this, rhs); return*this; }
		DualQuaternion DualQuaternion::operator+=(double x) { *this = add(*this, x); return *this; }
		DualQuaternion DualQuaternion::operator+(DualQuaternion rhs) { return add(*this, rhs); }
		DualQuaternion DualQuaternion::operator+(double x) { return add(*this, x); }

		DualQuaternion DualQuaternion::operator-=(DualQuaternion rhs) { *this = subtract(*this, rhs); return*this; }
		DualQuaternion DualQuaternion::operator-=(double x) { *this = subtract(*this, x); return *this; }
		DualQuaternion DualQuaternion::operator-(DualQuaternion rhs) { return subtract(*this, rhs); }
		DualQuaternion DualQuaternion::operator-(double x) { return subtract(*this, x); }

		//=======================
		std::wstring DualQuaternion::toString() { return R.toString() + L" + " + D.toString() + L"ɛ"; }
		void DualQuaternion::display() { std::wcout << toString() << std::endl; }
		//=========================================================================

		Octonian::Octonian() {}
		Octonian::Octonian(std::vector<double> v) { element = v; }
		Octonian::~Octonian() {};

			Matrix Octonian::CreateMatrix() {
				Matrix A(8, 8);
				for (int i = 0; i < 8; ++i) {
					A.set(i, i, element[i]);
				}
				return A;
			}

			/*void Octonian::CreateFromAxisAngle(double x, double y, double z, double deg) {
			double angle = double((deg / 180) * PI);
			double result = double(sin(angle / 2));
			w = double(cos(angle / 2));
			// Calculate the x, y and z of the Octonian
			x = double(x * result);
			y = double(y * result);
			z = double(z * result);
			}*/

			Octonian Octonian::conjugate() {
				std::vector<double> p = element;
				for (int i = 1; i < p.size(); ++i) { p[i] *= -1; }
				return Octonian(p);
			}

			double Octonian::norm() {
				double n = 0;
				for (int i = 0; i < 8; ++i) { n += pow(element[i], 2); }
				return sqrt(n);
			}

			/*
			double* Octonian::toArray() {
			double* v = new double[3];
			v[0] = x;
			v[1] = y;
			v[2] = z;
			return v;
			}

			std::vector<double> Octonian::toVector() {
			std::vector<double> v;
			v.push_back(x);
			v.push_back(y);
			v.push_back(z);
			return v;
			}

			double* Octonian::toEulerAngles() {
			double v[3];
			v[0] = atan2(2 * (w*x + y*z), 1 - 2 * (x*x + y*y));
			v[1] = asin(2 * (w*y - z*x));
			v[2] = atan2(2 * (w*z + x*y), 1 - 2 * (y*y + z*z));
			return v;
			}

			Octonian Octonian::add(Octonian q, double n) { q.w += n; q.x += n; q.y += n; q.z += n; return q; }
			Octonian Octonian::add(Octonian q, Octonian p) { q.w += p.w; q.x += p.x; q.y += p.y; q.z += p.z; return q; }
			Octonian Octonian::multiply(Octonian q, double n) { q.w *= n; q.x *= n; q.y *= n; q.z *= n; return q; }
			Octonian Octonian::inverse(Octonian q) {
			Octonian q_ = q.conjugate();
			q_ = multiply(q_, 1 / q.norm());
			return q_;
			}
			Octonian Octonian::HamiltonProduct(Octonian q, Octonian p) { //Hamilton Product
			double w_val = q.w*p.w - q.x*p.x - q.y*p.y - q.z*p.z;
			double x_val = q.w*p.x + q.x*p.w + q.y*p.z - q.z*p.y;
			double y_val = q.w*p.y + q.y*p.w + q.z*p.x - q.x*p.z;
			double z_val = q.w*p.z + q.z*p.w + q.x*p.y - q.y*p.x;
			q.w = w_val; q.x = x_val; q.y = y_val; q.z = z_val;
			return q;
			}

			Octonian Octonian::convertEulerAnglesToOctonian(double* v) {//    v = <pitch,yaw,roll> or <x-axis rotate, y rotate, z rotate>
			double q[4];

			double t0 = cos(v[2] * 0.5);
			double t1 = sin(v[2] * 0.5);
			double t2 = cos(v[1] * 0.5);
			double t3 = sin(v[1] * 0.5);
			double t4 = cos(v[0] * 0.5);
			double t5 = sin(v[0] * 0.5);

			q[0] = t0 * t2 * t4 + t1 * t3 * t5;
			q[1] = t0 * t3 * t4 - t1 * t2 * t5;
			q[2] = t0 * t2 * t5 + t1 * t3 * t4;
			q[3] = t1 * t2 * t4 - t0 * t3 * t5;

			return Octonian(q);
			}

			Octonian Octonian::rotateByOctonians(Octonian q, Octonian p) { // p' = q p q'
			Octonian q_ = inverse(q);
			q = HamiltonProduct(HamiltonProduct(q, p), q_);
			if (abs(q.w) < 0.00000000001) { q.w = 0; }
			if (abs(q.x) < 0.00000000001) { q.x = 0; }
			if (abs(q.y) < 0.00000000001) { q.y = 0; }
			if (abs(q.z) < 0.00000000001) { q.z = 0; }
			return q;
			}

			Octonian Octonian::createRotationOctonian(double* v) {
			/*  given a = angle to rotate
			[x, y, z] = axis to rotate around
			//then rotation Octonian R = [cos(a / 2), sin(a / 2)*x, sin(a / 2)*y, sin(a / 2)*z]
			return Octonian(cos(v[0] / 2), sin(v[0] / 2)*v[1], sin(v[0] / 2)*v[2], sin((v[0] / 2)*v[3]));

			}*/

			std::wstring Octonian::toString() {
				std::wostringstream sstr;
				sstr << L"<";
				for (int i = 0; i < 8; ++i) {
					if (element[i]<999999999999999999 && element[i] > -999999999999999999) {
						if (element[i] != floor(element[i])) { OctonianPrecision = precision; }
						else { OctonianPrecision = 0; }
						sstr << std::fixed << std::setprecision(OctonianPrecision) << element[i];
						if (i < element.size() - 1) { sstr << L","; }
					}
				}
				sstr << L">";
				std::wstring answer = sstr.str();
				sstr.clear();
				return answer;
			}

			void Octonian::display(double* v) { std::wcout << toString() << std::endl; }
		//END Octonian CLASS=========================================

