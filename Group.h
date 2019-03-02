#ifndef GROUP_H
#define GROUP_H
#pragma once
#include "stdafx.h"

class Group {
	//Any group can be characterized as a set of integers paired with the "permute" operation.
	//This is due to Cayley's Theorem - "Every group G is isomorphic to a permutation group."
public:
	std::vector<int> element;
	int GroupPrecision;
	Group();
	~Group();
	Group(int n);
	Group(std::vector<int> v);
	double order();
	void identity();
	void swap(int a, int b);
	void cycleLeft();
	void cycleRight();
	void groupActionRight(std::vector<int> act);
	void groupActionLeft(std::vector<int> act);
	int get(int i);
	void set(int i, int x);
	Vector toVector();
	std::wstring toString();
	void display();
};
#endif