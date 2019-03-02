#ifndef CHAR_H
#define CHAR_H
#pragma once
#include "stdafx.h"

wchar_t getLetterByAlphabeticalOrder(int n);
std::wstring getLetterByAlphabeticalOrderString(int n);
char ASCIItoChar(int a);
wchar_t ASCIItoWchar_t(int a);
bool isNumber(char a);
bool isNumber(wchar_t a);
bool isLetter(char a);
bool isLetter(wchar_t a);
bool isOperator(char c);
bool isOperator(wchar_t c);
bool isOperatorNotParen(char c);
bool isOperatorNotParen(wchar_t c);
char randomChar();
char randomWchar_t();
#endif