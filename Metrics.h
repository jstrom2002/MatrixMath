#ifndef METRICS_H
#define METRICS_H
#pragma once
#include "stdafx.h"
double EuclideanDistance(double p1, double p2);
double EuclideanDistance(std::vector<double> p1, std::vector<double> p2);
double EuclideanDistance(ComplexNumber p1, ComplexNumber p2);
double EuclideanDistance(std::vector<ComplexNumber> p1, std::vector<ComplexNumber> p2);
#endif