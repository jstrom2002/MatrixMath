#ifndef IMAGEPROCESSING_H
#define IMAGEPROCESSING_H
#pragma once
#include "stdafx.h"

/*
	NOTE:
	This header is for image processing, and will display image files as n x m x 3
	tensors	representing each entry as a 3D vector of RGB values, where the number of
	rows = the image height, and the number of columns = the image width,
	and the depth = 3 (one entry for each R,G,B value).
*/

class IMAGE {
public:
	int colorDepth; //read color depth to determine color array
	int compressionMethod;//tells what method of compression is used (0 = uncompressed) 
	double gamma;//gamma value for the immage, Vout = A*Vin^(gamma)
	int headerSize;//size of header before pixel array
	int hgt;//gives image width
	int horizontalResolution;
	int imageSize;//gives image size
	int numberOfBytes;//read number of bytes in file
	int numberOfColorsInPalatte;//colors in palatte, 0 is default (2^n colors)
	int numberOfPixels;//number of total pixel vectors, ie for 24-bit arrays it would be 3 * wdt * hgt
	int pixelsPerRow;//number of pixels in a row, width*3
	int padBytes;//number of bytes padding each row
	int reservedBytes;//field of the header
	int rowSize;//number of bytes necessary to store one row
	int verticalResolution;
	int wdt;//gives image height 
	std::wstring filename;
	std::wstring filetype;//the extension of the image file (ie BMP, JPEG, GIF, etc)
	std::wstring headertype;
	std::wstring colorOrdering;//displays whether the array values are RGB, BGR, grayscale, etc
	std::vector<unsigned char> header;//byte vector of the data for the header to BMP file
	std::vector<std::vector<unsigned char>> pixelArray;//an array of all the pixels in the BMP file 

	IMAGE();
	IMAGE(std::wstring fn);
	IMAGE(Matrix A);
	IMAGE(Matrix A, std::wstring fileType);
	~IMAGE();
	void convertToGrayscale();
	Matrix getBlackAndWhiteMatrix();
	void saveImage();
	void saveImage(std::wstring name);
	void printToTextFile();
	void printToTextFile(std::wstring outputFile);
	Vector toVector();
	std::vector<unsigned char> get(int i, int j);
	unsigned char get(int i, int j, int n);
	void set(int i, int j, std::vector<unsigned char> val);
	void set(int i, int j, int n, unsigned char val);
	Matrix getRedMatrix();
	Matrix getGreeMatrix();
	Matrix getBlueMatrix();
	void maskOutRed();
	void maskOutGreen();
	void maskOutBlue();
	void invertColor();
	void transpose();
	void convertBGRtoRBG();
	void adjustGamma(double gam);
	std::vector<unsigned char> makeHeader(std::wstring fileType);
	void updateHeader();
	int calculateRowSize(std::wstring type);
	int calculateNumberOfBytes(std::wstring type);
	std::wstring getHeadertype();
	Tensor toTensor();
	std::wstring toString();
	void display();
};
bool isBMP(std::wstring filename);
bool isGIF(std::wstring filename);
bool isJPEG(std::wstring filename);
bool isPNG(std::wstring filename);

void fileToText(std::wstring filename);
void BMPtoText(std::wstring filename, std::wstring outputName);
void BMPtoText(std::wstring filename);
unsigned long* make_crc_table();

Matrix identityFilter(unsigned int sz);
Matrix sharpeningFilter(unsigned int sz);
Matrix edgeDetectorFilter(unsigned int sz);
Matrix boxFilter(unsigned int sz);
Matrix GaussianKernel1D(int sz);
Matrix GaussianKernel1D(int sz,double sigma);
Matrix GaussianKernel2D(int sz);
Matrix GaussianKernel2D(int sz, double sigma);
Matrix gradientFilter(std::wstring type);
Matrix gradientFilter(std::wstring type, unsigned int wrt);
Matrix testMatrix();
#endif