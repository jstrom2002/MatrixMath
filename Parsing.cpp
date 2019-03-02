#pragma once
#include "Parsing.h"

std::vector<std::wstring> tokenizeByComma(std::wstring str) {
	std::vector<std::wstring> answer;
	while (str.find(L',') != std::wstring::npos) {
		answer.push_back(str.substr(0, str.find(L',')));
		str = str.substr(str.find(L',') + 1);
	}
	answer.push_back(str);
	return answer;
}

std::vector<std::wstring> tokenizeByChar(std::wstring str, wchar_t token) {
	std::vector<std::wstring> answer;
	while (str.find(token) != std::wstring::npos) {
		answer.push_back(str.substr(0, str.find(token)));
		str = str.substr(str.find(token) + 1);
	}
	answer.push_back(str);
	return answer;
}

double ParseInputNumber(std::wstring in) { return stringToDouble(in); }
double ParseInputNumber(std::string in) { return stringToDouble(s2ws(in)); }

std::vector<double> ParseNInputNumbers(std::string str) {
	return ParseNInputNumbers(s2ws(str));
}

std::vector<double> ParseNInputNumbers(std::wstring str) {
	str = removeSpaces(str);
	str = replaceConstants(str);	//replace a constant like "speedOfLight" in the string with the actual value for the speed of light

	std::vector<double> answer;
	std::vector<int> position;

	//get the number of values in the string by finding all ',' chars
	for (int i = 0; i < str.length(); ++i) {
		if (str[i] == ',') {
			position.push_back(i);
		}
		if (str[i] == ')') {//break if you've hit the end
			i = str.length();
		}
	}
	if (position.size() == 0) { answer.push_back(ParseInputNumber(str)); return answer; }//if there is only one single value in the string
	if (position.size() > 0) {//if there are multiple values in the string
		std::wstring tempStr = str.substr(0, position[0]);   //get first substring
		answer.push_back(stringToDouble(tempStr));
		for (int i = 1; i < position.size();++i) {
			tempStr = str.substr(position[i - 1] + 1, position[i] - position[i - 1] - 1);
			answer.push_back(stringToDouble(tempStr));
		}
		tempStr = str.substr(position[position.size() - 1] + 1);
		answer.push_back(stringToDouble(tempStr));
	}

	return answer;
}


ComplexNumber ParseComplexNumber(std::wstring in) {
	//if the input string is only real-valued and not complex, evaluate it as a real number
	if (in.find(L'i') == std::wstring::npos) {
		return ComplexNumber(ParseInputNumber(in));
	}

	//else if it is complex...
	while (isNumber(in[0]) == false && in[0] != '-') { in = in.substr(1); }
	if (in.find(L'i') != std::wstring::npos) {
		if (isOperator(in[in.find(L'i') - 1]) == true) {//case: '1-i' is the number
			in.insert(in.find(L'i'), L"1");
		}
		if (in.find(L"+") == false && in.find(L"-") == false) {//case: '3i' is the number
			in = replaceChar(in, L'i', L' ');
			in.append(L")");
			return ComplexNumber(0, ParseInputNumber(in));
		}
	}

	for (int i = 1; i < in.size(); ++i) {
		if (in[i] == L'-' && isNumber(in[i - 1]) == true) {
			in.insert(in.begin() + i, L',');
		}
	}
	in = replaceChar(in, L'+', L',');
	in = replaceChar(in, L'i', L' ');
	in.append(L")");
	std::vector<double> temp = ParseNInputNumbers(in);
	if (temp.size() == 1) { return ComplexNumber(temp[0]); }
	return ComplexNumber(temp[0], temp[1]);
}

std::vector<ComplexNumber> ParseNComplexNumbers(std::wstring str) {
	str = removeSpaces(str);
	std::vector<ComplexNumber> answer;
	std::vector<std::wstring> strs = tokenizeByComma(str);
	for (int i = 0; i < strs.size(); ++i) {
		if (findComplexNumber(strs[i]) != std::wstring::npos) {
			answer.push_back(ParseComplexNumber(strs[i]));
		}
		else if (findNearestDouble(strs[i]) != std::wstring::npos) {
			answer.push_back(ComplexNumber(ParseInputNumber(strs[i])));
		}
	}
	return answer;

	//OLD CODE -- DEPRECATED
	/*std::vector<int> position;

	//get the number of values in the string by finding all ',' chars
	for (int i = 0; i < str.length(); ++i) {
		if (str[i] == L',') {
			position.push_back(i);
		}
		if (str[i] == L')') {//break if you've hit the end
			i = str.length();
		}
	}
	if (position.size() == 0) { answer.push_back(ParseComplexNumber(str)); return answer; }//if there is only one single value in the string
	if (position.size() > 0) {//if there are multiple values in the string
		std::wstring tempStr = str.substr(0, position[0]);   //get first substring
		answer.push_back(ParseComplexNumber(tempStr));
		for (int i = 1; i < position.size();++i) {
			tempStr = str.substr(position[i - 1] + 1, position[i] - position[i - 1] - 1);
			answer.push_back(ParseComplexNumber(tempStr));
		}
		tempStr = str.substr(position[position.size() - 1] + 1);
		answer.push_back(ParseComplexNumber(tempStr));
	}
	return answer;*/
}

Matrix ParseMatrix(std::wstring inp) {//matrices should be declared like so, "[1,2,3;4,5,6;7,8,9]"
	std::vector<double> b;
	int rowCounter = 1;
	int colCounter = 1;

	inp = replaceString(inp, L"\t", L",    ");//replace all tabs with spaces

											  //Count # rows and columns
	bool firstCol = true;
	for (int i = 0; i < inp.length(); ++i) {
		if (inp[i] == ',' && firstCol == true) { ++colCounter; }
		if (inp[i] == ';' && i<(inp.length() - 2)) {
			++rowCounter;
			firstCol = false;
		}
	}

	inp = replaceChar(inp, ';', ',');
	inp = replaceChar(inp, '[', ' ');
	inp = replaceChar(inp, ']', ' ');
	if (inp.length() > 0) { inp = trim(inp); }
	while (isNumber(inp[inp.length() - 1]) == false) { inp.pop_back(); }//remove non-number final characters
	b = ParseNInputNumbers(inp);
	if (rowCounter*colCounter != b.size()) { std::wcout << L"ERROR!\n"; return Matrix(); }//check to make sure input is valid.
	return Matrix(rowCounter, colCounter, b);
}

Matrix ParseMatrixFromFile(std::wstring input) {
	std::vector<char> c(input.c_str(), input.c_str() + input.size() + 1u);
	std::vector<double> b;
	bool hasMore = true;
	double temp = 0;
	double tempDecimal = 0;
	int counter = 0;
	int rowCounter = 1; //assume at least a 1x1 matrix
	int colCounter = 1;
	int lengthCol = 0;
	bool pastDecimal = false;
	bool isNegative = false;

	input = replaceString(input, L"\t", L",   ");//replace all tabs with spaces
												 //===============================
	for (int i = 0; i<c.size(); ++i) {
		if (isNumber(c[i]) == true) {
			if (pastDecimal == false) {
				temp *= 10;
				temp += c[i] - '0';
			}
			if (pastDecimal == true) {
				tempDecimal *= 10;
				tempDecimal += c[i] - '0';
			}
		}
		if (c[i] == ',') {
			if (pastDecimal = true) {
				tempDecimal = convertIntToDecimalValue(tempDecimal);
				temp += tempDecimal;
				tempDecimal = 0;
			}
			if (isNegative == true) {
				temp *= -1;
				isNegative = false;
			}
			b.push_back(temp);
			temp = 0;
			pastDecimal = false;
			colCounter++;
		}
		if (c[i] == ';' && c[i + 1] != ']') {
			tempDecimal = convertIntToDecimalValue(tempDecimal);
			temp += tempDecimal;
			tempDecimal = 0;
			if (isNegative == true && temp>0) {
				temp *= -1;
				isNegative = false;
			}
			b.push_back(temp);
			temp = 0;
			lengthCol = colCounter;
			colCounter = 1;
			rowCounter++;
			pastDecimal = false;
		}
		if (c[i] == '.') {
			pastDecimal = true;
		}
		if (c[i] == '-') {
			isNegative = true;
		}
		counter++;
	}
	//====================end 'for' loop
	if (tempDecimal != 0) {
		tempDecimal = convertIntToDecimalValue(tempDecimal);
		temp += tempDecimal;
		tempDecimal = 0;
	}
	if (isNegative == true) { temp *= -1; }
	b.push_back(temp);
	temp = 0;

	if (rowCounter == 1) {
		lengthCol = colCounter;
	}

	return Matrix(rowCounter, lengthCol, b);
}

void ParseMatrixLoadFile(std::wstring filename) {
	std::wstring line;
	std::wifstream myfile(filename);
	bool isBoolean = false;
	if (myfile) {
		while (getline(myfile, line)) {
			//Matrices================
		topOfLoop:
			if (line[0] == 'M') {
				std::vector<std::wstring> finalStr;
				std::wstring nextLine;
				getline(myfile, nextLine);
				while (!myfile.eof() && nextLine.find(L"M") == std::wstring::npos) {   //while file isn't at end, and there's not a M in the next line
					finalStr.push_back(nextLine);
					if (line.find(L"Boolean") != std::wstring::npos) { isBoolean = true; }
					getline(myfile, nextLine);
				}
				if (nextLine.find(L"M") == std::wstring::npos && nextLine != L"") {
					finalStr.push_back(nextLine);
				}


				std::wstring f;
				//now format the string -- insert a ',' after the end of each number, then remove all the spaces in the string
				for (int k = 0; k < finalStr.size(); ++k) {
					bool firstSpace = true;
					for (int i = 0;i < finalStr[k].length();++i) {      //now format the entire string
						if (finalStr[k][i] == ' ' && firstSpace == true) { //the first space between matrix entries ought to be turned into a ','
							finalStr[k].erase(finalStr[k].begin() + i);
							finalStr[k].insert(i, L",");
							firstSpace = false;
						}
						if (isNumber(finalStr[k][i]) == true && finalStr[k][i] != L'.') {    //remove any spaces after the 1st space
							firstSpace = true;
						}
					}
					if (finalStr[k] != L"") {
						finalStr[k] = removeSpaces(finalStr[k]);
						while (finalStr[k][finalStr[k].length() - 1] == L' ') { finalStr[k].pop_back(); }
						if (finalStr[k][finalStr[k].length() - 1] == L',') { finalStr[k].pop_back(); }
						while (finalStr[k][finalStr[k].length() - 1] == L' ') { finalStr[k].pop_back(); }
						finalStr[k].append(L";");
						f.append(finalStr[k]);
						finalStr[k].clear();
					}
				}
				f.insert(0, L"[");
				f.append(L"]");
				Matrix M = ParseMatrix(f);
				if (isBoolean == true) { M.modulus = 2; }
				MatrixContainer.push_back(M);
				f.clear();
				finalStr.clear();
				line = nextLine;
				isBoolean = false;
				goto topOfLoop;
			}
		}
		myfile.close();
	}
	else {};
}

ComplexMatrix ParseComplexMatrix(std::wstring input) {//matrices should be declared like so, "[1,2,3;4,5,6;7,8,9;]"
	std::vector<char> c(input.c_str(), input.c_str() + input.size() + 1u);
	std::vector<double> b;
	bool hasMore = true;
	double temp = 0;
	double tempDecimal = 0;
	int counter = 0;
	int rowCounter = 1; //assume at least a 1x1 matrix
	int colCounter = 1;
	int lengthCol = 0;
	bool pastDecimal = false;
	bool isNegative = false;
	//===============================
	for (int i = 0; i<c.size(); ++i) {
		if (isNumber(c[i]) == true) {
			if (pastDecimal == false) {
				temp *= 10;
				temp += c[i] - '0';
			}
			if (pastDecimal == true) {
				tempDecimal *= 10;
				tempDecimal += c[i] - '0';
			}
		}
		if (c[i] == ',') {
			if (pastDecimal = true) {
				tempDecimal = convertIntToDecimalValue(tempDecimal);
				temp += tempDecimal;
				tempDecimal = 0;
			}
			if (isNegative == true) {
				temp *= -1;
				isNegative = false;
			}
			b.push_back(temp);
			temp = 0;
			pastDecimal = false;
			colCounter++;
		}
		if (c[i] == ';' && c[i + 1] != ']') {
			tempDecimal = convertIntToDecimalValue(tempDecimal);
			temp += tempDecimal;
			tempDecimal = 0;
			if (isNegative == true && temp>0) {
				temp *= -1;
				isNegative = false;
			}
			b.push_back(temp);
			temp = 0;
			lengthCol = colCounter;
			colCounter = 1;
			rowCounter++;
			pastDecimal = false;
		}
		if (c[i] == '.') {
			pastDecimal = true;
		}
		if (c[i] == '-') {
			isNegative = true;
		}
		counter++;
	}
	//====================end 'for' loop
	if (tempDecimal != 0) {
		tempDecimal = convertIntToDecimalValue(tempDecimal);
		temp += tempDecimal;
		tempDecimal = 0;
	}
	if (isNegative == true) { temp *= -1; }
	b.push_back(temp);
	temp = 0;

	if (rowCounter == 1) {
		lengthCol = colCounter;
	}

	std::vector<double> R;
	std::vector<double> C;
	for (int i = 0; i < b.size();++i) {
		if (i % 2 == 0) { R.push_back(b[i]); }
		if (i % 2 == 1) { C.push_back(b[i]); }
	}

	return ComplexMatrix(rowCounter, floor(lengthCol / 2), R, C);
}

DualNumberMatrix ParseDualNumberMatrix(std::wstring input) {
	std::vector<char> c(input.c_str(), input.c_str() + input.size() + 1u);
	std::vector<double> b;
	bool hasMore = true;
	double temp = 0;
	double tempDecimal = 0;
	int counter = 0;
	int rowCounter = 1; //assume at least a 1x1 matrix
	int colCounter = 1;
	int lengthCol = 0;
	bool pastDecimal = false;
	bool isNegative = false;
	//===============================
	for (int i = 0; i<c.size(); ++i) {
		if (isNumber(c[i]) == true) {
			if (pastDecimal == false) {
				temp *= 10;
				temp += c[i] - '0';
			}
			if (pastDecimal == true) {
				tempDecimal *= 10;
				tempDecimal += c[i] - '0';
			}
		}
		if (c[i] == ',') {
			if (pastDecimal = true) {
				tempDecimal = convertIntToDecimalValue(tempDecimal);
				temp += tempDecimal;
				tempDecimal = 0;
			}
			if (isNegative == true) {
				temp *= -1;
				isNegative = false;
			}
			b.push_back(temp);
			temp = 0;
			pastDecimal = false;
			colCounter++;
		}
		if (c[i] == ';' && c[i + 1] != ']') {
			tempDecimal = convertIntToDecimalValue(tempDecimal);
			temp += tempDecimal;
			tempDecimal = 0;
			if (isNegative == true && temp>0) {
				temp *= -1;
				isNegative = false;
			}
			b.push_back(temp);
			temp = 0;
			lengthCol = colCounter;
			colCounter = 1;
			rowCounter++;
			pastDecimal = false;
		}
		if (c[i] == '.') {
			pastDecimal = true;
		}
		if (c[i] == '-') {
			isNegative = true;
		}
		counter++;
	}
	//====================end 'for' loop
	if (tempDecimal != 0) {
		tempDecimal = convertIntToDecimalValue(tempDecimal);
		temp += tempDecimal;
		tempDecimal = 0;
	}
	if (isNegative == true) { temp *= -1; }
	b.push_back(temp);
	temp = 0;

	if (rowCounter == 1) {
		lengthCol = colCounter;
	}

	std::vector<double> R;
	std::vector<double> C;
	for (int i = 0; i < b.size();++i) {
		if (i % 2 == 0) { R.push_back(b[i]); }
		if (i % 2 == 1) { C.push_back(b[i]); }
	}

	std::vector<DualNumber> dvec;
	int sz = R.size();
	if (C.size() < R.size()) { sz = C.size(); }
	for (int i = 0; i < sz; ++i) { dvec.push_back(DualNumber(R[i], C[i])); }

	return DualNumberMatrix(rowCounter, floor(lengthCol / 2), dvec);
}

Tensor ParseTensor(std::wstring in) {
	// 'T[(1,2,3)(1e_1, 2e_2, 3e_3, 4e_4, 5e_5, 6e_6)]'
	while (in[0] != '(') { in = in.substr(1); }
	in = replaceChar(in, '[', ' ');
	in = replaceChar(in, ']', ' ');
	in = removeSpaces(in);

	std::wstring part1 = in.substr(0, in.find(L")") + 1);
	in.erase(0, in.find(L")") + 1);
	while (part1[0] != '(') { part1 = part1.substr(1); }
	part1 = part1.substr(1);

	std::wstring part2 = in.substr(0, in.find(L")") + 1);
	while (part2[0] != '(') { part2 = part2.substr(1); }
	part2 = part2.substr(1);

	bool isCoV = false;
	if (part2.find(L"^") != std::wstring::npos) { isCoV = true; }

	std::vector<double> vals;
	for (int i = 0; i < part2.size(); ++i) {
		if (part2.find(L",") != std::wstring::npos) {
			std::wstring temp = part2.substr(0, part2.find(L"e"));
			vals.push_back(ParseInputNumber(temp));
			i = part2.find(L",") + 2;
			part2 = part2.erase(0, part2.find(L",") + 1);
		}
	}
	std::wstring temp = part2.substr(0, part2.find(L"e"));//get very last value
	vals.push_back(ParseInputNumber(temp));

	std::vector<int> indices = STLDoubleToInt(ParseNInputNumbers(part1));
	return Tensor(vals, indices, isCoV);
}

Tensor NEWParseTensor(std::wstring in) {
	//initialization is of the form 'T[(1,2,3)(1,2,3,4,5,6)]' or T[(indices)(elements)]
	while (in[0] != L'(') { in = in.substr(1); }
	in = replaceChar(in, L'[', L' ');
	in = replaceChar(in, L']', L' ');
	in = removeSpaces(in);

	std::wstring part1 = in.substr(0, in.find(L")") + 1);
	in.erase(0, in.find(L")") + 1);
	part1 = replaceChar(part1, L'(', L' ');
	part1 = replaceChar(part1, L')', L' ');
	part1 = removeSpaces(part1);

	std::wstring part2 = in.substr(0, in.find(L")") + 1);
	part2 = replaceChar(part2, L'(', L' ');
	part2 = replaceChar(part2, L')', L' ');

	//determine if elements have '^' char determining that the tensor is covariant
	bool isCoV = false;
	if (part2.find(L"^") != std::wstring::npos) { isCoV = true; }
	part2 = replaceChar(part2, L'^', L' ');
	part2 = removeSpaces(part2);

	std::vector<double> vals = ParseNInputNumbers(part2);
	std::vector<int> indices = STLDoubleToInt(ParseNInputNumbers(part1));
	return Tensor(vals, indices, isCoV);
}

ComplexPolynomial ParseComplexPolynomial(std::wstring input) { //turns a polynomial of form 'p1(x) = x^3 + 2x + -1 to a vector of values, then constructs a polynomial from that vector
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
				if (j > 0 && vals[i][j] == L'^') {
					std::wstring s = vals[i].substr(j + 1);
					std::wstring c = vals[i].substr(0, j - 1);
					if (s[0] == L'^') { s[0] = L' '; }
					removeSpaces(s);
					while (s[0] == L' ') { s = s.substr(1); }//makes sure to remove blank space at the beginning of the string
					indices.push_back(s);
					bool isCoefficient = false;//check now to see if there actually is a coefficient or if we've got a case like 'x^4'
					for (int q = 0; q < c.size(); ++q) { if (isNumber(c[q]) == true) { isCoefficient = true;	q = c.size(); } }
					if (isCoefficient == false) { coefs.push_back(L"1"); }//if there are no numbers, put a "1" string in there
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

	//convert all strings to ComplexNumber//int values
	//=========================================
	int largestExponent = ParseInputNumber(indices[0]);
	++largestExponent;
	ComplexNumber* v = new ComplexNumber[largestExponent + 1];

	for (int i = 0; i < largestExponent + 1; ++i) {
		v[i] = 0;
	}
	for (int i = 0; i < indices.size(); ++i) {
		v[(int)ParseInputNumber(indices[i])] = ParseInputNumber(coefs[i]);
	}

	std::vector<ComplexNumber> vec;
	for (int i = 0; i < largestExponent; ++i) { vec.push_back(v[i]); }
	ComplexPolynomial p(vec);
	//p.display(); std::cin.ignore();//debug code
	return p;
}

DualNumber ParseDualNumber(std::wstring in) {
	while (isNumber(in[0]) == false && in[0] != L'-') { in = in.substr(1); }
	in = replaceChar(in, L'+', L',');
	in = replaceChar(in, L'e', L' ');
	in.append(L")");
	std::vector<double> temp = ParseNInputNumbers(in);
	if (temp.size() == 1) { return DualNumber(temp[0]); }
	return DualNumber(temp[0], temp[1]);
}

std::vector<int> ParseGroup(std::wstring input) {
	return STLDoubleToInt(ParseNInputNumbers(input));
}

ODE ParseODE(std::wstring in) {
	int colonpos = 0;
	int EQpos = 0;
	for (int i = 0; i < in.length(); i++) {
		if (in[i] == ':') { colonpos = i; }
		if (in[i] == '=') {
			EQpos = i;
			i = in.length();
		}
	}
	std::wstring lft = in.substr(colonpos + 1, EQpos - colonpos - 1);
	lft = trim(lft);
	std::wstring rt = in.substr(EQpos + 1);
	rt = trim(rt);
	ODE ode1(lft, rt);
	return ode1;
}

PDE ParsePDE(std::wstring in) {
	int colonpos = 0;
	int EQpos = 0;
	for (int i = 0; i < in.length(); i++) {
		if (in[i] == ':') { colonpos = i; }
		if (in[i] == '=') {
			EQpos = i;
			i = in.length();
		}
	}
	std::wstring lft = in.substr(colonpos + 1, EQpos - colonpos - 1);
	lft = trim(lft);
	std::wstring rt = in.substr(EQpos + 1);
	rt = trim(rt);
	PDE pde1(lft, rt);
	return pde1;
}

std::vector<double> ParseVector(std::wstring input) {
	return ParseNInputNumbers(input);
}

Quaternion ParseQuaternion(std::wstring input) {
	while (input[0] != '<') { input = input.substr(1); }
	input = replaceChar(input, '<', ' ');
	input = replaceChar(input, '>', ' ');
	input = removeSpaces(input);
	return Quaternion(ParseNInputNumbers(input));
}

DualQuaternion ParseDualQuaternion(std::wstring input) {
	while (input[0] != '<') { input = input.substr(1); }
	input = replaceChar(input, '<', ' ');
	input = replaceChar(input, '>', ' ');
	input = removeSpaces(input);
	std::vector<double> r;
	std::vector<double> d;
	std::vector<double> vals = ParseNInputNumbers(input);
	if (vals.size() > 7) {
		for (int i = 0; i < 4; ++i) {
			r.push_back(vals[i]);
			d.push_back(vals[i + 4]);
		}
	}
	return DualQuaternion(r, d);
}

SparseMatrix ParseSparseMatrix(std::wstring inp) {//sparse matrices should be declared like so, "S[1,2,3;4,5,6;7,8,9];rows,columns"
	std::vector<double> b;
	int rowCounter = 1;
	int colCounter = 1;

	inp = inp.substr(1);//knock off the 'S' at the beginning

						//separate the row and column values from the rest of the string
	int rws = 0;
	int cols = 0;
	std::vector<char> v;//rows
	std::vector<char> w;//columns
	bool pastComma = false;
	while (inp[inp.length() - 1] != ']') {
		if (isNumber(inp[inp.length() - 1]) == true && pastComma == false) { w.push_back(inp[inp.length() - 1]); }
		if (isNumber(inp[inp.length() - 1]) == true && pastComma == true) { w.push_back(inp[inp.length() - 1]); }
		if (inp[inp.length() - 1] == ',') { pastComma = true; }
		inp.pop_back();
	}
	std::wstring str(v.begin(), v.end());
	rws = ParseInputNumber(str);
	std::wstring str2(w.begin(), w.end());
	cols = ParseInputNumber(str2);

	//Count # rows and columns
	bool firstCol = true;
	int savedItems = 0;
	for (int i = 0; i < inp.length(); ++i) {
		if (inp[i] == ',' && firstCol == true) { ++colCounter; }
		if (inp[i] == ';' && i<(inp.length() - 2)) {
			++savedItems;
			++rowCounter;
			firstCol = false;
		}
	}

	inp = replaceChar(inp, ';', ',');
	inp = replaceChar(inp, '[', ' ');
	inp = replaceChar(inp, ']', ' ');
	inp = trim(inp);
	while (isNumber(inp[inp.length() - 1]) == false) { inp.pop_back(); }//remove non-number final characters
	b = ParseNInputNumbers(inp);
	if (rowCounter*colCounter != b.size()) { return SparseMatrix(); }//check to make sure input is valid.
	return SparseMatrix(rws, cols, b);
}

SparseMatrix ParseSparseMatrixFromFile(std::wstring input) {
	std::vector<char> c(input.c_str(), input.c_str() + input.size() + 1u);
	std::vector<double> b;
	bool hasMore = true;
	double temp = 0;
	double tempDecimal = 0;
	int counter = 0;
	int rowCounter = 1; //assume at least a 1x1 matrix
	int colCounter = 1;
	int lengthCol = 0;
	bool pastDecimal = false;
	bool isNegative = false;
	//===============================
	for (int i = 0; i<c.size(); ++i) {
		if (isNumber(c[i]) == true) {
			if (pastDecimal == false) {
				temp *= 10;
				temp += c[i] - '0';
			}
			if (pastDecimal == true) {
				tempDecimal *= 10;
				tempDecimal += c[i] - '0';
			}
		}
		if (c[i] == ',') {
			if (pastDecimal = true) {
				tempDecimal = convertIntToDecimalValue(tempDecimal);
				temp += tempDecimal;
				tempDecimal = 0;
			}
			if (isNegative == true) {
				temp *= -1;
				isNegative = false;
			}
			b.push_back(temp);
			temp = 0;
			pastDecimal = false;
			colCounter++;
		}
		if (c[i] == ';' && c[i + 1] != ']') {
			tempDecimal = convertIntToDecimalValue(tempDecimal);
			temp += tempDecimal;
			tempDecimal = 0;
			if (isNegative == true && temp>0) {
				temp *= -1;
				isNegative = false;
			}
			b.push_back(temp);
			temp = 0;
			lengthCol = colCounter;
			colCounter = 1;
			rowCounter++;
			pastDecimal = false;
		}
		if (c[i] == '.') {
			pastDecimal = true;
		}
		if (c[i] == '-') {
			isNegative = true;
		}
		counter++;
	}
	//====================end 'for' loop
	if (tempDecimal != 0) {
		tempDecimal = convertIntToDecimalValue(tempDecimal);
		temp += tempDecimal;
		tempDecimal = 0;
	}
	if (isNegative == true) { temp *= -1; }
	b.push_back(temp);
	temp = 0;

	if (rowCounter == 1) {
		lengthCol = colCounter;
	}

	return SparseMatrix(rowCounter, lengthCol, b);
}

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
					std::wstring c = vals[i].substr(0, j - 1);
					if (s[0] == '^') { s[0] = ' '; }
					removeSpaces(s);
					while (s[0] == ' ') { s = s.substr(1); }//makes sure to remove blank space at the beginning of the string
					indices.push_back(s);
					bool isCoefficient = false;//check now to see if there actually is a coefficient or if we've got a case like 'x^4'
					for (int q = 0; q < c.size(); ++q) { if (isNumber(c[q]) == true) { isCoefficient = true;	q = c.size(); } }
					if (isCoefficient == false) { coefs.push_back(L"1"); }//if there are no numbers, put a "1" string in there
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
	double* v = new double[largestExponent + 1];

	for (int i = 0; i < largestExponent + 1; ++i) {
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

FunctionMatrix ParseFunctionMatrix(std::wstring inp) {
	std::vector<std::wstring> vec;
	int rowCounter = 1;
	int colCounter = 1;

	//Count # rows and columns
	bool firstCol = true;
	for (int i = 0; i < inp.length(); ++i) {
		if (inp[i] == ',' && firstCol == true) { ++colCounter; }
		if (inp[i] == ';' && i<(inp.length() - 2)) {
			++rowCounter;
			firstCol = false;
		}
	}

	inp = replaceChar(inp, ';', ',');
	inp = replaceChar(inp, '[', ' ');
	inp = replaceChar(inp, ']', ' ');
	inp = trim(inp);

	while (inp.find(L",") != std::wstring::npos) {
		vec.push_back(trim(inp.substr(0, inp.find(L","))));
		inp = inp.substr(inp.find(L",") + 1);
	}
	vec.push_back(inp);

	return FunctionMatrix(rowCounter, colCounter, vec);
}

std::vector<std::wstring> ParseFunction(std::wstring in) {//parse input string of form "'f(x,...)='" into two strings " 'exp(x^2...)...etc', 'x,y,z,...etc' " 
	std::wstring str1;
	std::wstring str2;
	std::vector<std::wstring> answer;

	int leftbracket = 0;
	int rightbracket = 0;
	int EQsign = 0;
	bool needLB = true;
	bool needRB = true;
	bool needEQ = true;
	for (int i = 0; i < in.length(); ++i) {
		if (in[i] == ')' && needRB == true) {
			rightbracket = i;
			needRB = false;
		}
		if (in[i] == '(' && needLB == true) {
			leftbracket = i;
			needLB = false;
		}
		if (in[i] == '=' && needEQ == true) {
			EQsign = i;
			needEQ = false;
		}
	}
	if (needLB == true || needRB == true || needEQ == true) { return answer; }
	str1 = in.substr(EQsign + 1);//actual function
								 //removeSpaces(str1);//can't have any spaces in the function
	str2 = in.substr(leftbracket + 1, rightbracket - 2);//list of variables
	str1 = trim(str1);
	str2 = trim(str2);

	//Function is declared like so:  f[0] = function, f[1] = variables
	answer.push_back(str2);
	answer.push_back(str1);
	return answer;
}

std::vector<Function> ParseNFunctions(std::wstring str) {
	str = removeSpaces(str);
	std::vector<Function> answer;
	std::vector<int> position;

	//get the number of values in the string by finding all ',' chars
	for (int i = 0; i < str.length(); ++i) {
		if (str[i] == L',') {
			position.push_back(i);
		}
		if (str[i] == L')') {//break if you've hit the end
			i = str.length();
		}
	}
	if (position.size() == 0) { answer.push_back(Function(ParseFunction(str))); return answer; }//if there is only one single value in the string
	if (position.size() > 0) {//if there are multiple values in the string
		std::wstring tempStr = str.substr(0, position[0]);   //get first substring
		answer.push_back(Function(ParseFunction(tempStr)));
		for (int i = 1; i < position.size();++i) {
			tempStr = str.substr(position[i - 1] + 1, position[i] - position[i - 1] - 1);
			answer.push_back(Function(ParseFunction(tempStr)));
		}
		tempStr = str.substr(position[position.size() - 1] + 1);
		answer.push_back(Function(ParseFunction(tempStr)));
	}
	return answer;
}

std::vector<wchar_t> varsToChars(std::wstring vars) {
	std::vector<wchar_t> variables;
	for (int i = 0; i < vars.size(); ++i) {
		if (isLetter(vars[i]) == true && vars[i] != ',') { variables.push_back(vars[i]); }
	}
	return variables;
}