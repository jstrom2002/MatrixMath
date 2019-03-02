/* this header is for operating upon integers of the 'unsigned int'
data type, treating them as binary values. */
#pragma once
#include "stdafx.h"

bool toggleState(bool n) {
	if (n == true) { return false; }
	return true;
}

bool XOR(bool A, bool B) {
	A = !A;
	B = !B;
	return (A + B) % 2;
}

int getValue(int x, int n)
{//returns the value at position n of int x
	int temp = 0;
	int temp2 = 0;
	int counter = 0;
	int length = 0;

	//get length of n
	temp = n;
	while (temp >= 1) {
		temp /= 10;
		length++;
	}
	//handle case where 'x' is a single digit
	if (length == 1 && n == 1) {
		return x;
	}

	//remove each value of 'x' on the LHS of 'n'
	temp = n;
	for (int i = length - 1; i >= n; --i) {
		temp2 = temp;
		while (temp2 > 1) {
			temp2 /= 10;
		}
		temp *= pow(10, i);
		temp -= temp;
		printf("%d\n", temp);
	}

	//shift the result to the right until only the desired
	//digit remains
	while (temp >= 10) {
		temp /= 10;
		printf("%d\n", temp);

	}

	return temp;
}

unsigned int getBit(unsigned int x, int n)
{
	unsigned int temp = 1;
	temp = temp << (n - 1);
	x = x & temp;
	while (x > 1) {
		x = x >> 1;
	}
	return x;
}

void setBit(unsigned int* x, int n)
{
	//set bit n of binary number x
	unsigned int temp = 1;
	temp = temp << (n - 1);
	*x = *x | temp;
}

void clearBit(unsigned int* x, int n)
{
	//set bit n of binary number x
	unsigned int temp = 1;
	temp = temp << (n - 1);
	temp = ~temp;
	*x = *x & temp;
}

unsigned int getBitByte(unsigned int x, int byte, int bit)
{
	//returns x bit at n byte
	unsigned int temp = 1;
	temp = temp << (((byte - 1) * 8) + (bit - 1));
	x = x & temp;
	while (x > 1) {
		x = x >> 1;
	}
	return x;
}

void setBitByte(unsigned int* num, int byte, int bit)
{
	//sets x bit at n byte of a binary number
	unsigned int temp = 1;
	temp = temp << ((byte - 1) * 8);
	temp = temp << (bit - 1);
	*num = *num | temp;
}

unsigned int bitLength(unsigned int x)
{
	//returns the length of a binary number
	unsigned int len = 0;
	while (x != 0x0) {
		x = x >> 1;
		len++;
	}
	return len;
}

void printBinary(unsigned int x)
{
	//prints a binary number in proper format
	unsigned int len = bitLength(x);
	int counter = 0;

	//pad byte with extra 0's
	if (len % 4 != 0 || x == 0) {
		for (int i = 0; i < 4 - (len % 4); ++i) {
			printf("0");
			counter++;
		}
	}
	//print out bits in number
	for (int i = len; i>0; --i) {
		if (i % 4 == 0 && len>4 && i != len) {
			printf(" ");
		}
		printf("%x", getBit(x, i));
	}
	printf("\n");
}

unsigned int rotate_right(unsigned int x, int n)
{
	unsigned int temp, temp2, len;
	temp = temp2 = len = 0;

	len = bitLength(x); //store length of input binary number

						// save the values that will be knocked off the number
						// after it is shifted to the right n times.
	for (int i = n - 1; i >= 0; --i) {
		temp2 = 1;
		temp2 = temp2 << i;
		temp2 = temp2 & x;
		temp = temp | temp2;
	}
	temp = temp & x;

	// shift input number to the right n times
	x = x >> n;

	// shift the knocked off numbers to the left fully,
	// and then 'OR' the results of this shift with the
	// now right-shifted input number
	temp = temp << (len - 1);
	x = x | temp;
	return x;
}

unsigned int rotate_left(unsigned int x, int n)
{
	unsigned int temp, len;
	temp = len = 0;

	len = bitLength(x); // store length of input binary number

						// get numbers knocked off by the left shift,
						// so that they can be appended to the number
	temp = x >> (len - n);

	for (int i = len; i > len - n; --i) {
		clearBit(&x, i);
	}

	// shift input number to the left n times
	x = x << n;

	// add on the knocked off numbers back to the right side of the binary number
	x = x | temp;
	return x;
}

int charToInt(char c) {return c - '0';}

int charToInt(wchar_t c) {return c - L'0';}

int charToInt(char* str, int sz)
{
	int i;
	int val = 0;
	int temp = sz;

	for (i = 0; i<temp; ++i) {
		if (str[i] == '\0') { i = sz + 1; }
		else if (str[i] != '0') {
			val += ((int)str[i] - '0') * (int)pow(10, temp - i - 1);
		}
	}
	return val;
}

int charToInt(wchar_t* str, int sz)
{
	int i;
	int val = 0;
	int temp = sz;

	for (i = 0; i<temp; ++i) {
		if (str[i] == L'\0') { i = sz + 1; }
		else if (str[i] != L'0') {
			val += ((int)str[i] - L'0') * (int)pow(10, temp - i - 1);
		}
	}
	return val;
}

char* intToChar(int n)
{
	char* str;
	int temp = 0;
	int counter = 0;
	int length = 0;

	// get length of n
	temp = n;
	while (temp >= 1) {
		temp /= 10;
		length++;
	}
	str = (char*)malloc(sizeof(char)*(length + 1));

	//handle case where 'n' is a single digit
	if (length == 1) {
		str[0] = ((char)getValue(n, 1) + '0');
		str[1] = '\0';
		return str;
	}


	// get each individual value at each postiion of n,
	// then save them each as a char in str
	temp = n;
	for (int i = length - 1; i >= 0; --i) {
		temp = getValue(n, i + 1);
		str[i] = temp + '0';
	}
	str[length] = '\0';
	return str;
}

wchar_t* intToWchar(int n)
{
	wchar_t* str;
	int temp = 0;
	int counter = 0;
	int length = 0;

	// get length of n
	temp = n;
	while (temp >= 1) {
		temp /= 10;
		length++;
	}
	str = (wchar_t*)malloc(sizeof(wchar_t)*(length + 1));

	//handle case where 'n' is a single digit
	if (length == 1) {
		str[0] = ((wchar_t)getValue(n, 1) + L'0');
		str[1] = '\0';
		return str;
	}

	// get each individual value at each postiion of n,
	// then save them each as a char in str
	temp = n;
	for (int i = length - 1; i >= 0; --i) {
		temp = getValue(n, i + 1);
		str[i] = temp + L'0';
	}
	str[length] = L'\0';
	return str;
}

unsigned int parityTest(unsigned int x)
{
	unsigned int counter = 0;
	for (int i = 0; i < bitLength(x); ++i) {
		counter += getBit(x, i);
	}
	return counter % 2;
}

unsigned int flipBinary(unsigned int x)
{
	unsigned int temp = 0;
	unsigned int len = bitLength(x);
	for (int i = len; i >= 1; --i) {
		temp = temp | (getBit(x, i) << (len - i));
	}
	return temp;
}

unsigned int flipByteLength(unsigned int x, unsigned int sz)
{
	unsigned int temp = 0;
	unsigned int len = bitLength(x);
	for (int i = len; i >= 1; --i) {
		temp = temp | (getBit(x, i) << (len - i));
	}
	temp = (temp << (sz - len));
	return temp;
}

unsigned int getByte(unsigned int x, int n)
{
	if (n == 0) { return 0; }
	unsigned int temp = 0;
	unsigned int temp2 = 0;
	for (int i = 8; i >= 1; --i) {
		temp2 = getBit(x, ((n - 1) * 8) + i);
		if (temp2 == 1) {
			temp = temp | (temp2 << (i - 1));
		}
	}
	return temp;
}

unsigned int swapBytes(unsigned int x, unsigned int bytes)
{
	unsigned int temp = 0;
	for (int i = bytes; i >= 1; --i) {
		temp = temp | (getByte(x, i) << ((bytes - i) * 8));
	}

	return temp;
}

unsigned int swapReverseBytes(unsigned int x)
{
	unsigned int bytes = bitLength(x) / 8;
	unsigned int temp = 0;
	for (int i = bytes; i >= 1; --i) {
		temp = temp | (flipByteLength(getByte(x, i), 8) << ((i - 1) * 8));
	}

	return temp;
}

unsigned char circshift(unsigned char x, int n) {
	return (x << n) | (x >> (8 - n));
}

/*unsigned short ntohs(unsigned short x) {return swapBytes(x, 2);}
unsigned short ntohl(unsigned short x) {return swapBytes(x, 4);}
unsigned short htons(unsigned short x) {return swapBytes(x, 2);}
unsigned short htonl(unsigned short x) {return swapBytes(x, 4);}*/