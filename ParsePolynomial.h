#pragma once
#include "stdafx.h"
#include "Polynomial.h"
#include "vector functions.h"

Polynomial ParsePolynomial(std::wstring input) { //turns a polynomial of form 'p1(x) = x^3 + 2x + -1 to a vector of values, then constructs a polynomial from that vector
	input = removeSpaces(input);
	//cut off anything proceeding the 'x^3 + ...' part like a 'p1(x) = ' statement
	if (input.find(L"=") != std::wstring::npos) {
		input = input.substr(input.find(L"=") + 1);
	}
	while (input[0] == ' ') { input = input.substr(1); }//remove blank space at the beginning of the string

	//insert a '+' in front of any '-' sign in case there isnt' one already there
	for (int i = 0; i < input.length(); ++i) {
		if (input[i] == '-') {
			if (isLetter(input[i + 1]) == true) { input.insert(input.begin() + i + 1, '1'); }
			for (int k = i; k >= 0; --k) {//go backwards from i and check if you hit a number.  if so, add a '+'
				if (isNumber(input[k]) == true) { input.insert(input.begin() + i, '+');  break; }
				if (input[k] == '+') { break; }
			}
		}
	}

	//SEPARATE VALUES
	//===============
	//cut up input string into smaller strings and place them in the "vals" string vector
	std::vector<std::wstring> vals;
	for (int i = 0; i < input.length(); ++i) {
		if (i > 0 && input[i] == '+') {
			std::wstring s = input.substr(0, i);
			if (s[0] == '+') { s[0] = ' '; }
			removeSpaces(s);
			while (s[0] == ' ') { s = s.substr(1); }//makes sure to remove blank space at the beginning of the string
			vals.push_back(s);
			input.erase(0, i);
			i = 0;
		}
	}
	//deal with the last of the input string
	if (input[0] == '+') { input[0] = ' '; }
	removeSpaces(input);
	while (input[0] == ' ') { input = input.substr(1); }//makes sure to remove blank space at the beginning of the string
	vals.push_back(input);


	//GET MATCHING EXPONENT-COEFFICIENT PAIRS
	//=======================================
	//now cut the strings in the vector into two parts:  the coefficient and the exponent
	std::vector<std::wstring> indices;
	std::vector<std::wstring> coefs;
	for (int i = 0; i < vals.size(); ++i) {
		if (vals[i].find(L"^") == std::wstring::npos) {
			if (vals[i].find(L"x") != std::wstring::npos) {//case: no "^" but an x.  i.e. '23.23x'
				indices.push_back(L"1");
				std::wstring temp = vals[i];
				bool onlyX = true;
				for (int j = 0; j < temp.length(); ++j) {//check to see if it's just 'x' or if it has a coefficient (i.e. '2.34x')
					if (isNumber(temp[j]) == true) { onlyX = false; }
				}
				if (onlyX == true) { temp.insert(temp.begin(), '1'); }//put a '1' at the beginning if there isn't one				
				std::replace(temp.begin(), temp.end(), 'x', ' ');
				removeSpaces(temp);
				coefs.push_back(temp);
			}
			if (vals[i].find(L"x") == std::wstring::npos) {//case: no o'^', no 'x'.  A single value alone (i.e. '12.32x^0')
				indices.push_back(L"0");
				coefs.push_back(vals[i]);
			}
		}
		else {
			for (int j = 0; j < vals[i].size(); ++j) {//case" a value with an 'x' and a '^' which needs separating in to coefficient and exponent
				if (j > 0 && vals[i][j] == '^') {
					std::wstring s = vals[i].substr(j + 1);
					std::wstring c = vals[i].substr(0, j-1);
					if (s[0] == '^') { s[0] = ' '; }
					removeSpaces(s);
					while (s[0] == ' ') { s = s.substr(1); }//makes sure to remove blank space at the beginning of the string
					indices.push_back(s);
					bool isCoefficient = false;//check now to see if there actually is a coefficient or if we've got a case like 'x^4'
					for (int q = 0; q < c.size(); ++q) { if (isNumber(c[q]) == true) { isCoefficient = true;	q = c.size(); } }
					if (isCoefficient == false) {coefs.push_back(L"1");	}//if there are no numbers, put a "1" string in there
					else { coefs.push_back(c); }//else, treat normally
				}
			}
		}
	}

/*	DEBUG CODE
	==========
	displayarr(coefs);
	std::wcout << L"\n";
	displayarr(indices);
	std::wcout << L"========\n";
	std::cin.ignore();	*/
	
	//convert all strings to double//int values
	//=========================================
	int largestExponent = ParseInputNumber(indices[0]);
	++largestExponent;
	double* v = new double[largestExponent+1];

	for (int i = 0; i < largestExponent+1; ++i) {
		v[i] = 0;
	}
	for (int i = 0; i < indices.size(); ++i) {
		v[(int)ParseInputNumber(indices[i])] = ParseInputNumber(coefs[i]);
	}
	std::vector<double> vec = toSTLVector(v, largestExponent);
	Polynomial p(vec);
	//p.display(); std::cin.ignore();//debug code
	return p;
}