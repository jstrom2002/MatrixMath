#ifndef ARITHMETIC_H
#define ARITHMETIC_H
#pragma once
//this header is for arbitrary precision arithmetic using strings
bool isInteger(std::wstring n);
std::wstring getMantissa(std::wstring n);
std::wstring getInteger(std::wstring n);
bool lessThan(std::wstring n, std::wstring m);
bool lessThanEqual(std::wstring n, std::wstring m);
std::wstring shiftLeft(std::wstring n, int p);
std::wstring shiftRight(std::wstring n, int p);
std::wstring abs(std::wstring n);
std::wstring add(std::wstring n, std::wstring m);
std::wstring subtract(std::wstring n, std::wstring m);
std::wstring multiply(std::wstring n, std::wstring m);
std::wstring divide(std::wstring n, std::wstring m);
#endif