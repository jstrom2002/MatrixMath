#ifndef SECURITY_H
#define SECURITY_H
#pragma once
//#include "stdafx.h"
#include <string>

bool isEncryptedChar(char c);
std::wstring encrypt(std::wstring str);
std::wstring decrypt(std::wstring str);
std::wstring generateProductKey();
std::wstring generateSerialNumber();
bool testProductKey(std::wstring pk);
void writeToBinary(std::wstring file, std::wstring str);
std::wstring readFromBinary(std::wstring file, int stringLength, unsigned int size);
void clearFile(std::wstring file);
void writeToTxtFile(std::wstring file, std::wstring str);
std::wstring readFromTxtFile(std::wstring file);
std::vector<std::wstring> GetIPAddresses();
bool sendSecurityReport();
bool checkProductKey();
#endif