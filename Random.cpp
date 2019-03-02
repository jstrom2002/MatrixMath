#pragma once
#include "stdafx.h"

double randDouble() { return ((double)rand() / RAND_MAX)*pow(-1, rand() % 2); }//returns a random value between [-1,1]

double fRand(double fMin, double fMax, int y) {
	srand(y);                       // seeds the random number generator
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

double boundedRand(double a, double b){
	int rv = rand();
	return fRand(a, b, rv);   // the random number is within [x,k] and has a random integer as a seed
}

int pseudorandom(int x_prev, int a, int c, int m) {// generate a pseudorandom number in range [a,b] by linear congruence method
	/*
		For this generator to work, for X_n = (a * X_n-1 + c) mod m
		1) m and c must be relatively prime
		2) a-1 is divisible by all prime factors of m
		3) a-1 is divisible by 4 if m is divisible by 4.
	*/

	return (a * x_prev + c) % m;
}