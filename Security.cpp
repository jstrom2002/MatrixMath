#pragma once
#include "stdafx.h"
#include "String.h"
#define _WIN32_WINNT 0x501
#include <ws2tcpip.h>
#include "Security.h"

//int encChars[10] = { 190,128,161,33,247,157,39,63,95,244 };

bool isEncryptedChar(char c) {
	if (c == 190 || c == 128 || c == 161 || c == 33 || c == 247 || c == 157 || c == 39 || c == 63 || c == 95 || c == 244) { return true; }
	return false;
}

std::wstring encrypt(std::wstring str) {
	std::wstring enc = L"";
	int counter = str.length() - 1;
	srand(time(0) + clock());
	while (counter >= 0) {//add char from str
		int r = rand();
		if ((r % 1271) == 0) {
			if (str[counter] == L'0') { enc.push_back((char)190); }
			if (str[counter] == L'1') { enc.push_back((char)128); }
			if (str[counter] == L'2') { enc.push_back((char)161); }
			if (str[counter] == L'3') { enc.push_back((char)33); }
			if (str[counter] == L'4') { enc.push_back((char)247); }
			if (str[counter] == L'5') { enc.push_back((char)157); }
			if (str[counter] == L'6') { enc.push_back((char)39); }
			if (str[counter] == L'7') { enc.push_back((char)63); }
			if (str[counter] == L'8') { enc.push_back((char)95); }
			if (str[counter] == L'9') { enc.push_back((char)244); }
			if (str[counter] == L'-') { enc.push_back('q'); }
			--counter;
		}
		else {//add "noise" chars, make sure they're not from the set of chars that will make up the SN
			char c = (r % 255) + 33;
			while (isEncryptedChar(c) == true || c == L'q' || c > 244 || c < 34 || c == 127) {
				c = (rand() % 255) + 33;
			}
			enc.append(charToString(c));
		}
	}
	enc = removeSpaces(enc);
	return enc;
}

std::wstring decrypt(std::wstring str) {
	std::wstring decrypted = L"";
	for (int i = 0; i < str.length(); ++i) {
		if (str[i] == (wchar_t)190) { decrypted.push_back('0'); }
		if (str[i] == (wchar_t)128) { decrypted.push_back('1'); }
		if (str[i] == (wchar_t)161) { decrypted.push_back('2'); }
		if (str[i] == (wchar_t)33) { decrypted.push_back('3'); }
		if (str[i] == (wchar_t)247) { decrypted.push_back('4'); }
		if (str[i] == (wchar_t)157) { decrypted.push_back('5'); }
		if (str[i] == (wchar_t)39) { decrypted.push_back('6'); }
		if (str[i] == (wchar_t)63) { decrypted.push_back('7'); }
		if (str[i] == (wchar_t)95) { decrypted.push_back('8'); }
		if (str[i] == (wchar_t)244) { decrypted.push_back('9'); }
		if (str[i] == L'q') { decrypted.push_back('-'); }
	}
	decrypted = reverseString(decrypted);
	return decrypted;
}

std::wstring generateProductKey() {
	srand(clock() + time(0));
	std::vector<std::wstring> valsStr;
	std::vector<double> coef;

	for (int i = 0; i < 6; ++i) {
		int rand1 = rand();
		while (rand1 > 10000) { rand1 = floor(rand1 / 10); }
		coef.push_back(rand1);
		std::wstring temp = std::to_wstring(rand1);
		while (temp.length() < 4) { temp.insert(temp.begin(), '0'); }
		valsStr.push_back(temp);
	}

	long double val = 0;
	val += stringToDouble(valsStr[3]) / 10000.0;
	val += stringToDouble(valsStr[4]) / 100000000.0;
	val += stringToDouble(valsStr[5]) / 1000000000000.0;

	coef.push_back(stringToDouble(valsStr[0]));
	coef.push_back(stringToDouble(valsStr[1]));
	coef.push_back(stringToDouble(valsStr[2]));
	long double testval = (coef[0] * pow(val, 16)) + (coef[1] * pow(val, 15)) + (coef[2] * pow(val, 14));

	int temp = precision;
	precision = 20;
	std::wstring chkstr = to_stringPrecision(testval);
	precision = temp;
	while (chkstr[0] != L'.') { chkstr = chkstr.substr(1); }//remove everything up to the mantissa for the test string
	chkstr = chkstr.substr(1);
	for (int i = 0; i < 3; ++i) {
		std::wstring temp = L"";
		for (int j = 0; j < 4; ++j) {
			temp.append(charToString(chkstr[j]));
			chkstr = chkstr.substr(1);
		}
		valsStr.push_back(temp);
	}

	std::wstring answer = L"";
	for (int i = 0; i < valsStr.size() - 1; ++i) {
		answer.append(valsStr[i]);
		answer.append(L"-");
	}
	answer.append(valsStr[valsStr.size() - 1]);
	return answer;
}

std::wstring generateSerialNumber() {
	const std::wstring nextSerialNumber = L"0000-0000-0002";
	return nextSerialNumber;
}

bool testProductKey(std::wstring pk) {
	pk = trim(pk);
	std::vector<std::wstring> strs;
	while (pk.length() > 0) {
		strs.push_back(pk.substr(0, pk.find(L"-")));
		pk = pk.substr(pk.find(L"-") + 1);
		if (pk.find(L"-") == std::wstring::npos) {
			strs.push_back(pk);
			pk = L"";
		}
	}

	//setup test polynomial
	std::vector<double> coef;
	coef.push_back(stringToDouble(strs[0]));
	coef.push_back(stringToDouble(strs[1]));
	coef.push_back(stringToDouble(strs[2]));


	//get value for evaluation
	long double val = 0;
	val += stringToDouble(strs[3]) / 10000.0;
	val += stringToDouble(strs[4]) / 100000000.0;
	val += stringToDouble(strs[5]) / 1000000000000.0;

	//evaluate polynomial
	long double testval = (coef[0] * pow(val, 16)) + (coef[1] * pow(val, 15)) + (coef[2] * pow(val, 14));

	//convert result to string for checking
	int temp = precision;
	precision = 20;
	std::wstring chkstr = to_stringPrecision(testval);
	precision = temp;
	while (chkstr[0] != L'.') { chkstr = chkstr.substr(1); }//remove everything up to the mantissa for the test string
	chkstr = chkstr.substr(1);

	//make the last 3 strings into a single string to check
	std::wstring last3 = L"";
	last3.append(strs[6]);
	last3.append(strs[7]);
	last3.append(strs[8]);

	//do comparison of evaluated polynomial vs last third of pk
	for (int i = 0; i < last3.length(); ++i) {
		if (chkstr[i] != last3[i]) { return false; }
	}
	return true;
}

void writeToBinary(std::wstring file, std::wstring str) {
	std::wstring ref = str;
	unsigned int size = str.size();

	std::wofstream supp_info_output(file, std::ios::out | std::ios::binary); // saving file
	unsigned int stringLength = ref.length();
	supp_info_output.write((wchar_t*)(&stringLength), sizeof(stringLength));
	supp_info_output.write(ref.c_str(), ref.length());
	supp_info_output.write((wchar_t*)(&size), sizeof(size));
	supp_info_output.close();
}

std::wstring readFromBinary(std::wstring file, int stringLength, unsigned int size) {
	std::wstring ref2;

	std::ifstream supp_info_input(file, std::ios::in | std::ios::binary); // loading file
	supp_info_input.read((char*)(&stringLength), sizeof(stringLength));
	ref2.resize(stringLength);
	supp_info_input.read((char*)ref2.c_str(), stringLength);
	supp_info_input.read((char*)(&size), sizeof(size));
	supp_info_input.close();

	return ref2;
}

void clearFile(std::wstring file) {//deletes all data saved in the .txt file
	std::wofstream outf(file);
	outf.open(file, std::wofstream::out | std::wofstream::trunc);
	outf.close();
}

void writeToTxtFile(std::wstring file, std::wstring str) {
	std::wofstream outf2(file);
	outf2 << str;
	outf2.close();
}

std::wstring readFromTxtFile(std::wstring file) {
	std::wstring answer = L"";
	std::wstring line;
	std::wifstream myfile(file);
	if (myfile) {
		while (myfile.eof() == false) {
			getline(myfile, line);
			answer.append(line);
		}
		myfile.close();
	}
	else {
		std::wstring s = L"";
		return s;
	};
	return line;
}

std::vector<std::wstring> GetIPAddresses() {
	std::vector<std::wstring> result;

	std::wostringstream os;
	WSAData wsaData;
	if (WSAStartup(MAKEWORD(1, 1), &wsaData) != 0) { return result; }
	int retval;

	//get hostname
	//============
	char ac[80];
	if (gethostname(ac, sizeof(ac)) == SOCKET_ERROR) { return result; }
	os << ac;
	hostname = os.str();
	os.clear();
	os.flush();

	// We use getaddrinfo (gethostbyname has been deprecated)
	struct addrinfo hints = { 0 };
	hints.ai_family = AF_UNSPEC;    // Want both IPv4 and IPv6
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_protocol = IPPROTO_TCP;

	struct addrinfo *paddrinfo = NULL;

	if (getaddrinfo(ws2s(hostname).c_str(), NULL, &hints, &paddrinfo) != 0) { return result; }

	// The machine can have multiple IP addresses (IPv4, IPv6, etc.)
	for (struct addrinfo *ptr = paddrinfo; ptr != NULL;ptr = ptr->ai_next)
	{
		// inet_ntop is not available for all versions of Windows, we implement our own
		char ipaddress[NI_MAXHOST] = { 0 };
		if (ptr->ai_family == AF_INET)
		{
			if (getnameinfo(ptr->ai_addr, sizeof(struct sockaddr_in), (PCHAR)ipaddress, _countof(ipaddress) - 1, NULL, 0, NI_NUMERICHOST) == 0)
				result.push_back(s2ws(ipaddress));
		}
		else if (ptr->ai_family == AF_INET6)
		{
			if (getnameinfo(ptr->ai_addr, sizeof(struct sockaddr_in6), (PCHAR)ipaddress, _countof(ipaddress) - 1, NULL, 0, NI_NUMERICHOST) == 0)
				result.push_back(s2ws(ipaddress));
		}
	}

	freeaddrinfo(paddrinfo);

	WSACleanup();
	return result;
}

bool sendSecurityReport() {
	//add function later to email IP address and user details to my gmail
	return false;
}

bool checkProductKey() {
	std::wstring dec = decrypt(readFromTxtFile(L"sec.txt"));
	IPs = GetIPAddresses();
	serialnumber = dec.substr(dec.length() - 14);
	productKey = dec.substr(0, dec.length() - 15);
	if (dec.length() > 0) {
		return testProductKey(dec);
	}
	return false;
}
