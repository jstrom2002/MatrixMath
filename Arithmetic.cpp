#pragma once
#include "stdafx.h"

bool isInteger(std::wstring n) {
	if (n.find(L'.') == std::wstring::npos) { return true; }
	return false;
}

std::wstring getMantissa(std::wstring n) {
	if (isInteger(n) == true) { return L""; }
	return n.substr(n.find(L'.')+1, n.length()-n.find(L'.')-1);
}

std::wstring getInteger(std::wstring n) {
	if (isInteger(n) == true) { return L""; }
	return n.substr(0, n.find(L'.'));
}

bool lessThan(std::wstring n, std::wstring m) {
	if (n == m) { return false; }
	if (n[0] == L'-' && m[0] != L'-') { return true; }
	if (n[0] != L'-' && m[0] == L'-') { return false; }
	if (n[0] == L'-' && m[0] == L'-') { 
		n = n.substr(1);
		m = m.substr(1);
		return toggleState(lessThan(n, m));
	}

	std::wstring mantN;
	std::wstring mantM;
	if (isInteger(n) == false) {
		mantN = getMantissa(n);
		n = getInteger(n);
	}
	if(isInteger(m)==false){
		mantM = getMantissa(m);
		m = getInteger(m);
	}

	//compare integer components
	if (n != m) {
		if (n.length() < m.length()) { return true; }
		if (n.length() > m.length()) { return false; }
		for (int i = 0; i < n.length(); ++i) {
			if (n[i] < m[i]) { return true; }
			if (n[i] > m[i]) { return false; }
		}
	}
	
	//is integer components are identical, compare mantissas
	while (mantN.length() < mantM.length()) { mantN.append(L"0"); }
	while (mantM.length() < mantN.length()) { mantM.append(L"0"); }
	for (int i = 0; i < mantN.length(); ++i) {
		if (mantN[i] < mantM[i]) { return true; }
		if (mantN[i] > mantM[i]) { return false; }
	}

	return false;
}

bool lessThanEqual(std::wstring n, std::wstring m) {
	if (n == m) { return true; }
	return lessThan(n, m);
}

std::wstring shiftLeft(std::wstring n, int p) {//shift the decimal point in a number to the left
	int decPlace = n.find(L'.');
	if (p <= 0 || decPlace-p < 0 || decPlace == std::wstring::npos) { return n; }
	n.insert(n.begin() + (decPlace - p), L'.');
	n.erase(decPlace + 1, 1);
	if (n[0] == L'.') { n.insert(n.begin(),L'0'); }
	return n;
}

std::wstring shiftRight(std::wstring n, int p) {//shift the decimal point in a number to the right
	int decPlace = n.find(L'.');
	if (p <= 0) { return NULL; }
	if (decPlace == std::wstring::npos) { 
		n.append(L"0");
		return n; 
	}
	while (n.length() < n.find(L'.') + p) { n.append(L"0"); }
	n.insert(n.begin() + (decPlace + p + 1), L'.');
	n.erase(decPlace, 1);
	if (n[n.length()-1] == L'.') { n.insert(n.end(), L'0'); }
	return n;
}

std::wstring abs(std::wstring n) {
	if (n[0] == L'-') { n = n.substr(1); }
	return n;
}

std::wstring add(std::wstring n, std::wstring m) {
	if (m == L"0") { return n; }
	//handle negative number input
	bool isNegative = false;
	if (n[0] == L'-' && m[0] != L'-') {	return subtract(m, n.substr(1)); }
	if (n[0] != L'-' && m[0] == L'-') { return subtract(n, m.substr(1)); }
	if (n[0] == L'-' && m[0] == L'-') { 
		isNegative = true; 
		n.erase(n.begin());
		m.erase(m.begin());
	}

	//handle numbers with floating-point decimal values
	std::wstring mantissa = L"";
	std::wstring mantRem = L"";
	if (n.find(L'.') != std::wstring::npos || m.find(L'.') != std::wstring::npos){//see if there's a mantissa on either number
		//condition: only one number has a mantissa. Then, just add that mantissa onto the result at the end.
		if (n.find(L'.') != std::wstring::npos && m.find(L'.') == std::wstring::npos) {
			mantissa = n.substr(n.find(L'.') + 1);
			n = n.substr(0, n.find(L'.'));
		}
		if(n.find(L'.') == std::wstring::npos && m.find(L'.') != std::wstring::npos){
			mantissa = m.substr(m.find(L'.') + 1);
			m = m.substr(0, m.find(L'.'));
		}
		if (n.find(L'.') != std::wstring::npos && m.find(L'.') != std::wstring::npos) {//if both numbers have mantissas, add them
			std::wstring a = n.substr(n.find(L'.') + 1);
			std::wstring b = m.substr(m.find(L'.') + 1);
			int len = a.length();
			if (b.length() > len) { 
				len = b.length(); 
				while (a.length() < len) {
					a.append(L"0");
				}
			}
			while (b.length() < len) { b.append(L"0"); }
			mantRem = add(a, b);
			for (int i = 0; i < len; ++i) {
				mantissa.insert(mantissa.begin(), mantRem[mantRem.length() - 1]);
				mantRem.pop_back();
			}
			n = n.substr(0, n.find(L'.'));
			m = m.substr(0, m.find(L'.'));
		}
	}

	//else, if only two integers, add normally
	std::wstring answer = L"";
	int carry = 0;
	bool secondLarger = false;
	int length = m.length();
	int endlength = n.length();
	if (m.length()>n.length()) { //determine which string is larger
		length = n.length();
		endlength = m.length();
		secondLarger = true;
	};
	if (secondLarger == true) {
		for (int i = 0; i<endlength - length; i++) {
			n.insert(n.begin(), L'0');
		}
	}
	if (secondLarger == false) {
		for (int i = 0; i<endlength - length; i++) {
			m.insert(m.begin(), L'0');
		}
	}
	for (int i = endlength; i>0; i--) {
		std::wstring temp = to_stringPrecision(charToInt(n[i - 1]) + charToInt(m[i - 1]) + carry); //add each value of strings
		if (temp.length()>1) {
			answer.insert(answer.begin(), (temp[1]));
			carry = charToInt(temp[0]);
			if (i == 1) { answer.insert(answer.begin(), (temp[0])); }
		}
		if (temp.length() == 1) {
			answer.insert(answer.begin(), (temp[0]));
			carry = 0;
		}
	}

	//add mantissa remainder to final value if necessary
	if (mantRem.length() > 0) {
		answer = add(answer, mantRem);
	}
	//tack on mantissa to the end if needed
	if (mantissa.size() > 0) {
		answer.append(L".");
		answer.append(mantissa);
	}
	//make negative if it must be
	if (isNegative == true) {
		std::wstring ans = L"-";
		ans.append(answer);
		answer = ans;
	}
	return answer;
}

std::wstring subtract(std::wstring n, std::wstring m) {
	if (m == L"0") { return n; }
	if (m[0] == L'-') { return add(n, m.substr(1)); }
	bool isNegative = false;
	if (n[0] == L'-') {
		std::wstring temp = subtract(m, n);
		temp.insert(temp.begin(), L'-');
		return temp;
	}

	if (n == m) { return L"0"; }

	//handle numbers with floating-point decimal values
	std::wstring mantN = L"";
	std::wstring intN = n;
	std::wstring mantM = L"";
	std::wstring intM = m;

	if (isInteger(n) == false) {
		mantN = getMantissa(n);
		intN = getInteger(n);
	}
	if(isInteger(m)==false){
		intM = getInteger(m);
		mantM = getMantissa(m);
	}

	while (intN.length() < intM.length()) { intN.insert(intN.begin(), L'0'); }
	while (intN.length() > intM.length()) { intM.insert(intM.begin(), L'0'); }
	while (mantN.length() < mantM.length()) { mantN.append(L"0"); }
	while (mantN.length() > mantM.length()) { mantM.append(L"0"); }

	n = intN;
	n.append(mantN);
	m = intM;
	m.append(mantM);

	int decPlace = mantN.length();

	std::wstring answer = L"";
	int borrow = 0;
	bool secondLarger = false;
	int length = m.length();
	int endlength = n.length();
	if (m.length()>n.length()) {  //determine which wstring is longer
		std::wstring tempSwap = n;
		n = m;
		m = tempSwap;
		length = m.length();
		endlength = n.length();
		secondLarger = false;
		isNegative = true;
	};
	if (secondLarger == true) {         //make both wstrings the same size by tacking on 0's in front of the smaller number
		for (int p = 0; p<endlength - length; p++) {
			n.insert(n.begin(), '0');
		}
	}
	if (secondLarger == false) {
		for (int p = 0; p<endlength - length; p++) {
			m.insert(m.begin(), '0');
		}
	}

	bool isLarger = false;
	for (int p = 0; p<n.length();p++) {     //determine if the second value is numerically larger
		if (n[p]>m[p]) { isLarger == true;  p = n.length(); }
		if (n[p]<m[p] && isLarger == false) {
			std::wstring tempSwap = n;
			n = m;
			m = tempSwap;
			isNegative = true;
			p = n.length();
		}
	}

	for (int p = endlength; p>0; p--) {
		if (charToInt(n[p - 1]) < charToInt(m[p - 1])) { // if necessary, borrow 10
			for (int j = p - 1; j>0; j--) {
				if (charToInt(n[j - 1]) > charToInt(m[j - 1])) {
					n[j - 1] -= 1;
					borrow = 10;
					for (int k = j; k<p - 1; k++) { n[k] += 9; }
					j = 0;
				}
			}
		}
		std::wstring temp = to_stringPrecision((charToInt(n[p - 1]) + borrow) - charToInt(m[p - 1])); //add each value of wstrings
		if (temp.length() == 1) { answer.insert(answer.begin(), temp[0]); }
		if (temp.length()>1) {
			answer.insert(answer.begin(), temp[1]);
			answer.insert(answer.begin(), temp[0]);
		}
		borrow = 0;
	}


	//now eliminate 0's at front of wstring
	while (answer[0] == '0') { 
		answer.erase(answer.begin());
	}
	//insert decimal point if necessary
	if (decPlace > 0) {
		answer.insert(answer.end()-decPlace ,L'.');
	}
	//make negative if it must be
	if (isNegative == true) {
		std::wstring ans = L"-";
		ans.append(answer);
		answer = ans;
	}
	return answer;
}

std::wstring multiply(std::wstring n, std::wstring m) {
	if (n == L"0" || m == L"0") { return L"0"; }
	if (m == L"1") { return n; }
	if (n == L"1") { return m; }
	
	std::wstring answer = abs(n);
	std::wstring mantissa = abs(n);
	bool isNegative = false;
	if (n[0] == L'-' && m[0] != L'-') { isNegative = true; }
	if (n[0] != L'-' && m[0] == L'-') { isNegative = true; }
	n = abs(n);
	m = abs(m);

	std::wstring mantN = L"";
	std::wstring intN = n;
	std::wstring mantM = L"";
	std::wstring intM = m;

	if (isInteger(n) == false) {
		mantN = getMantissa(n);
		intN = getInteger(n);
	}
	if (isInteger(m) == false) {
		intM = getInteger(m);
		mantM = getMantissa(m);
	}

	int decPlace = mantM.length();

	//multiply whole number component of m
	while (lessThan(L"1",intM)) {
		answer = add(answer, n);
		intM = subtract(intM, L"1");
	}
	//multiply mantissa of m
	while (lessThan(L"1", mantM)) {
		mantissa = add(mantissa, n);
		mantM = subtract(mantM, L"1");
	}
	//add the whole and mantissa values together
	if (decPlace > 0) {
		mantissa = shiftLeft(mantissa, decPlace);
		answer = add(answer, mantissa);
	}
	//make negative if necessary
	if (isNegative == true) {
		answer.insert(answer.begin(), L'-');
	}
	return answer;
}


std::wstring divide(std::wstring n, std::wstring m) {
	if (m == L"0") {return NULL;} 
	if (n == L"0") { return L"0"; }
	if (m == L"1") { return n; }

	std::wstring answer = L"0";
	std::wstring mantissa = L"";
	std::wstring remainder = L"";
	
	bool isNegative = false;
	if (n[0] == L'-' && m[0] != L'-') { isNegative = true; }
	if (n[0] != L'-' && m[0] == L'-') { isNegative = true; }
	n = abs(n);
	m = abs(m);

	std::wstring mantN = L"";
	std::wstring intN = n;
	std::wstring mantM = L"";
	std::wstring intM = m;

	if (isInteger(n) == false) {
		mantN = getMantissa(n);
		intN = getInteger(n);
	}
	if (isInteger(m) == false) {
		intM = getInteger(m);
		mantM = getMantissa(m);
	}

	//repeatedly subtract the divisor from the numerator until it cannot be done further
	while (lessThanEqual(m,n)==true) {
		n = subtract(n, m);
		answer = add(answer, L"1");
	}
	remainder = n;

	//divide the remainder to come up with a mantissa
	int offset = 0;
	while (lessThan(remainder,m) == true) {
		remainder = shiftRight(remainder, 1);
		++offset;
	}
	if (offset >= 1) { --offset; }

	//do repeated subtraction to calculate value
	std::wstring lastDecimal = L"";
	while (lastDecimal.length() < precision && remainder != L"0") {//set 3 as the max length of a decimal value
		while (lessThanEqual(m, remainder) == true) {
			remainder = subtract(remainder, m);
			mantissa = add(mantissa, L"1");
		}
		if (remainder != L"0") { 
			remainder = shiftRight(remainder, 1); 
			lastDecimal = mantissa;
			mantissa = shiftRight(mantissa, 1);
		}
		else { lastDecimal = mantissa; }
	}

	//put in offset 0's into mantissa
	mantissa = L"";
	while (offset>0) {
		mantissa.insert(mantissa.begin(), L'0');
		--offset;
	}
	mantissa.append(lastDecimal);

	//handle mantissa
	if (mantissa.length() > 0) {
		answer.append(L".");
		answer.append(mantissa);
	}
	return answer;
}
