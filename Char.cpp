#pragma once
#include "stdafx.h"

wchar_t getLetterByAlphabeticalOrder(int n) {
	wchar_t x[26] = {
		'a','b','c','d','e',
		'f','g','h','i','j',
		'k','l','m','n','o',
		'p','q','r','s','t',
		'u','v','w','x','y','z'
	};
	return x[n];
}

std::wstring getLetterByAlphabeticalOrderString(int n) {
	std::wstring x[26] = {
		L"a",L"b",L"c",L"d",L"e",
		L"f",L"g",L"h",L"i",L"j",
		L"k",L"l",L"m",L"n",L"o",
		L"p",L"q",L"r",L"s",L"t",
		L"u",L"v",L"w",L"x",L"y",L"z"
	};
	return x[n];
}

char ASCIItoChar(int a) {
	if (a >= 0 && a < 256) {
		char a1 = a;
		return a1;
	}
	else { return NULL; }
}

wchar_t ASCIItoWchar_t(int a) {
	if (a >= 0 && a < 256) {
		wchar_t a1 = a;
		return a1;
	}
	else { return NULL; }
}

bool isNumber(char a) {
	if (a > 47 && a < 58) { return true; }
	return false;
}

bool isNumber(wchar_t a) {
	if (a > 47 && a < 58) { return true; }
	return false;
}

bool isLetter(char a) {
	if ((a > 64 && a < 91) || (a > 96 && a <123)) { return true; }
	return false;
}

bool isLetter(wchar_t a) {
	if ((a > 64 && a < 91) || (a > 96 && a <123)) { return true; }
	return false;
}

bool isOperator(char c) {//operators are +,-,*,/,^,(,)
	if (
		c == '+' ||
		c == '-' ||
		c == '*' ||
		c == '/' ||
		c == '^' ||
		c == '%' ||
		c == '(' ||
		c == ')'
		) {
		return true;
	}
	else { return false; }
}

bool isOperator(wchar_t c) {//operators are +,-,*,/,^,(,)
	if (
		c == L'+' ||
		c == L'-' ||
		c == L'*' ||
		c == L'/' ||
		c == L'^' ||
		c == L'%' ||
		c == L'(' ||
		c == L')'
		) {
		return true;
	}
	else { return false; }
}

bool isOperatorNotParen(char c) {//operators are +,-,*,/,^
	if (
		c == '+' ||
		c == '-' ||
		c == '*' ||
		c == '/' ||
		c == '^' ||
		c == '%'
		) {
		return true;
	}
	else { return false; }
}

bool isOperatorNotParen(wchar_t c) {//operators are +,-,*,/,^
	if (
		c == L'+' ||
		c == L'-' ||
		c == L'*' ||
		c == L'/' ||
		c == L'^' ||
		c == L'%'
		) {
		return true;
	}
	else { return false; }
}

char randomChar() {
	srand(time(0) + clock());
	char c = (rand() % 255) + 33;
	while (c == 127) { c = (rand() % 255) + 33; }
	return c;
}

char randomWchar_t() {
	srand(time(0) + clock());
	wchar_t c = (rand() % 255) + 33;
	while (c == 127) { c = (rand() % 255) + 33; }
	return c;
}