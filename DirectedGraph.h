#pragma once
#include "stdafx.h"
#include "Matrix.h"
#include "polynomial.h"

class DirectedGraph {
public:
	Matrix adjacencyMatrix;
	Matrix incidinceMatrix;
	int edges;
	int vertices;

public: DirectedGraph() {};
public: ~DirectedGraph() {};

public: DirectedGraph(Matrix Adj) {

}
		
	std::wstring toString() {
		std::wostringstream s;
		s << adjacencyMatrix.toString() << L"\n" << incidinceMatrix.toString();
		std::wstring answer = s.str();
		s.clear();
		return answer;
	}

	void display() { std::wcout << toString() << std::endl; }

};