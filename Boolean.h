#ifndef BOOLEAN_H
#define BOOLEAN_H
/* this header is for operating upon integers of the 'unsigned int' 
data type, treating them as binary values. */
#pragma once

bool toggleState(bool n);
bool XOR(bool A, bool B);
int getValue(int x, int n);
unsigned int getBit(unsigned int x, int n);
void setBit(unsigned int* x, int n);
void clearBit(unsigned int* x, int n);
unsigned int getBitByte(unsigned int x, int byte, int bit);
void setBitByte(unsigned int* num, int byte, int bit);
unsigned int bitLength(unsigned int x);
void printBinary(unsigned int x);
unsigned int rotate_right(unsigned int x, int n);
unsigned int rotate_left(unsigned int x, int n);
int charToInt(char c);
int charToInt(wchar_t c);
int charToInt(char* str, int sz);
int charToInt(wchar_t* str, int sz);
char* intToChar(int n);
wchar_t* intToWchar(int n);
unsigned int parityTest(unsigned int x);
unsigned int flipBinary(unsigned int x);
unsigned int flipByteLength(unsigned int x, unsigned int sz);
unsigned int getByte(unsigned int x, int n);
unsigned int swapBytes(unsigned int x, unsigned int bytes);
unsigned int swapReverseBytes(unsigned int x);
unsigned char circshift(unsigned char x, int n);
/*unsigned short ntohs(unsigned short x);
unsigned short ntohl(unsigned short x);
unsigned short htons(unsigned short x);
unsigned short htonl(unsigned short x);*/
#endif