#pragma once
#include "stdafx.h"
#include "Matrix.h"
#include "polynomial.h"

class UndirectedGraph {
public:
	Matrix adjacencyMatrix;
	Matrix incidenceMatrix;

public: UndirectedGraph() {};
public: ~UndirectedGraph() {};

public : UndirectedGraph(Matrix adj, Matrix inc) {;
		 adjacencyMatrix = adj;
		 incidenceMatrix = inc;
}

	int vertices() { return incidenceMatrix.rows; }
	int edges() { return adjacencyMatrix.rows; }

	std::wstring toString() {
		std::wostringstream s;
		s << adjacencyMatrix.toString() << L"\n" << incidenceMatrix.toString();
		std::wstring answer = s.str();
		s.clear();
		return answer;
	}

	void display() {std::wcout << toString() << std::endl;}

};