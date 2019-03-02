#pragma once
#include "stdafx.h"
/*
	NOTE:
	This header is for image processing, and will display image files as 2x2x3
	tensors	representing each entry as a 3D vector of RGB values, where the number of 
	rows = the image height, and the number of columns = the image width,
	and the depth = 3 (one entry for each R,G,B value).
*/

IMAGE::IMAGE() {}
IMAGE::~IMAGE() {}

IMAGE::IMAGE(std::wstring fn) {

	if (isBMP(fn) == true) {
		filename = fn;
		filetype = L"BMP";
		colorOrdering = L"BGR";

		//open BMP file
		FILE *f = fopen(ws2s(filename).c_str(), "rb"); 

		//read preliminary file data -- 14 bytes
		unsigned char prelimData[14];
		fread(prelimData, sizeof(unsigned char), 14, f);
		for (int i = 0; i < 14; ++i) { header.push_back(prelimData[i]); }
		numberOfBytes = (header[5] << 24) ^ (header[4] << 16) ^ (header[3] << 8) ^ header[2];//read number of bytes in file
		reservedBytes = (header[9] << 24) ^ (header[8] << 16) ^ (header[7] << 8) ^ header[6];//read reserved data
		headerSize = (header[13] << 24) ^ (header[12] << 16) ^ (header[11] << 8) ^ header[10];//read starting address
		headertype = getHeadertype();

		//read and interpret file header data
		unsigned char* headerData = new unsigned char[headerSize-14];
		fread(headerData, sizeof(unsigned char), headerSize-14, f);	 //read the 54-byte header
		for (int i = 0; i < headerSize-14; ++i) { header.push_back(headerData[i]); }

		//initialize class variables;
		wdt = (header[21] << 24) ^ (header[20] << 16) ^ (header[19] << 8) ^ header[18];//gives image height 
		hgt = (header[25] << 24) ^ (header[24] << 16) ^ (header[23] << 8) ^ header[22];//gives image width
		horizontalResolution = (header[41] << 24) ^ (header[40] << 16) ^ (header[39] << 8) ^ header[38];
		verticalResolution = (header[45] << 24) ^ (header[44] << 16) ^ (header[43] << 8) ^ header[42];
		colorDepth = (header[27] << 8) ^ header[28]; //read color depth to determine color array
		compressionMethod = (header[33] << 24) ^ (header[32] << 16) ^ (header[31] << 8) ^ header[30];//tells what method of compression is used (0 = uncompressed) 
		numberOfColorsInPalatte = (header[49] << 24) ^ (header[48] << 16) ^ (header[47] << 8) ^ header[46];//gives number of colors in color palatte, 0 = 2^n (default)
		gamma = 1;
		updateHeader();

		std::vector<unsigned char> data;
		unsigned char data2;
		while (fread(&data2, sizeof(unsigned char), 1, f)) {
			data.push_back(data2);
		}
		fclose(f);		

		int pixelCount = 0;//counts the number of pixel vectors read
		std::vector<unsigned char> tempVec;
			if (colorDepth == 24) {
				for (int i = 0; i < data.size(); ++i) {
					if (i%rowSize < rowSize-padBytes) {
						if (i > 0 && i % 3 == 0) {
							pixelArray.push_back(tempVec);
							tempVec.clear();
						}
						tempVec.push_back(data[i]);
					}
				}
		}
		fclose(f);
	}

	if (isPNG(fn) == true) {
		filename = fn;
		filetype = L"PNG";

		FILE *f = fopen(ws2s(filename).c_str(), "rb");
		int colorType = 0;
		int filterMethod = 0;
		int interlaceMethod = 0;

		//read preliminary file header -- 8 bytes
		unsigned char prelimData[8];
		fread(prelimData, sizeof(unsigned char), 8, f);
		for (int i = 0; i < 8; ++i) { header.push_back(prelimData[i]); }		

		//read chunks
		bool moreChunks = true;
		while (moreChunks==true) {			
			//read chunk
			unsigned char length[4];
			fread(length, sizeof(unsigned char), 4, f);
			for (int i = 0; i < 4; ++i) { header.push_back(length[i]); }
			unsigned int len = (length[0] << 24) ^ (length[1] << 16) ^ (length[2] << 8) ^ length[3];

			unsigned char type[4];
			fread(type, sizeof(unsigned char), 4, f);
			for (int i = 0; i < 4; ++i) { header.push_back(type[i]); }
			std::string chunkType = { (char)type[0], (char)type[1], (char)type[2], (char)type[3] };

			unsigned char* data = new unsigned char[len];
			fread(data, sizeof(unsigned char), len, f);
			for (int i = 0; i < len; ++i) { header.push_back(data[i]); }

			unsigned char CRC[4];
			fread(CRC, sizeof(unsigned char), 4, f);
			for (int i = 0; i < 4; ++i) { header.push_back(CRC[i]); }

			//handle chunk data
			if (chunkType == "IHDR" && len == 13) {
			/*  IHDR Byte/Data Fields:
				---------------------
				Width:				4 bytes
				Height:             4 bytes
				Bit depth:          1 byte
				Color type:         1 byte
				Compression method: 1 byte
				Filter method:      1 byte
				Interlace method:   1 byte   
				
				Color   Allowed     Interpretation
				Type    Bit Depths
				0       1,2,4,8,16  Each pixel is a grayscale sample.
				2       8,16        Each pixel is an R,G,B triple.
				3       1,2,4,8     Each pixel is a palette index; a PLTE chunk must appear.
				4       8,16        Each pixel is a grayscale sample, followed by an alpha sample.
				6       8,16        Each pixel is an R,G,B triple, followed by an alpha sample.			*/

				wdt = (data[0] << 24) ^ (data[1] << 16) ^ (data[2] << 8) ^ data[3];
				hgt = (data[4] << 24) ^ (data[5] << 16) ^ (data[6] << 8) ^ data[7];
				colorDepth = data[8];
				colorType = data[9];
				compressionMethod = data[10];
				filterMethod = data[11];
				interlaceMethod = data[12];
			}
			if (chunkType == "IDAT") {
				if (colorDepth == 8) {
					//skip the first two bytes and the last four
					for (int i = 2; i < len-4; ++i) {

					}
				}
			}
			if (chunkType == "IEND") {
				moreChunks = false;
			}

			//Auxillary chunks
			if (chunkType == "gAMA") {
				gamma = ((data[0] << 24) ^ (data[1] << 16) ^ (data[2] << 8) ^ data[3]) / 100000.0;
			}
			if (chunkType == "bKGD") {

			}
			if (chunkType == "cHRM") {
			/*	The cHRM chunk contains:
				White Point x: 4 bytes
				White Point y: 4 bytes
				Red x:         4 bytes
				Red y:         4 bytes
				Green x:       4 bytes
				Green y:       4 bytes
				Blue x:        4 bytes
				Blue y:        4 bytes	*/

			}
			if (chunkType == "hIST") {

			}
			if (chunkType == "pHYs") {

			}
			if (chunkType == "sBIT") {

			}
			if (chunkType == "tEXt") {

			}
			if (chunkType == "tIME") {

			}
			if (chunkType == "tRNS") {

			}
			if (chunkType == "zTXt") {

			}
		}		

	}

	if (isJPEG(fn) == true) {
		filename = fn;
		filetype = L"JPEG";

		FILE *f = fopen(ws2s(filename).c_str(), "r");

	}
}

IMAGE::IMAGE(Matrix A) {//default file type is BMP
	*this = IMAGE(A, L"BMP");
}

IMAGE::IMAGE(Matrix A, std::wstring type) {//save matrix as grayscale image
	
	if (type == L"BMP") {
		for (int i = 0; i < A.size(); ++i) {
			std::vector<unsigned char> vec(3, 0);
			pixelArray.push_back(vec);
		}

		for (int i = 0; i < A.rows; ++i) {
			for (int j = 0; j < A.columns; ++j) {
				pixelArray[i*A.columns + j][0] = pixelArray[i*A.columns + j][1] = pixelArray[i*A.columns + j][2] = std::floor(A.get(i, j));
			}
		}

		colorDepth = 24;
		compressionMethod = 0;
		gamma = 1.0;
		horizontalResolution = verticalResolution = 2835;
		headerSize = 54;//size of header before pixel array
		hgt = A.rows;
		wdt = A.columns;
		numberOfColorsInPalatte = 0;//colors in palatte, 0 is default (2^n colors)
		numberOfPixels = 3 * wdt * hgt;//number of total pixel vectors, ie for 24-bit arrays it would be 3 * wdt * hgt
		reservedBytes = 0;
		rowSize = calculateRowSize(L"BMP");//number of bytes needed to store one row
		pixelsPerRow = wdt * 3;
		padBytes = rowSize - pixelsPerRow;//if the width is not divisible by 4, then it is filled in by padding (ie 10 mod 4 = 2 padded entries)		
		numberOfBytes = calculateNumberOfBytes(L"BMP");
		imageSize = numberOfBytes - headerSize;
		filename = L"matrix.bmp";
		filetype = type;
		header = makeHeader(L"BMP");
	}

}

void IMAGE::convertToGrayscale() {
	for (int i = 0; i < pixelArray.size(); ++i) {
		/* this method uses the same process as the
		cv2.cvtColor(img, COLOR_RGB2Gray) method */
		//get values for RGB in range [0,1]
		double R = pixelArray[i][0] / 255.0;
		double G = pixelArray[i][1] / 255.0;
		double B = pixelArray[i][2] / 255.0;
		double Y = (R*0.299 + G * 0.587 + B * 0.114);//calculate luminance from formula
		pixelArray[i][0] = pixelArray[i][1] = pixelArray[i][2] = (Y*255);
	}
}

Matrix IMAGE::getBlackAndWhiteMatrix() {
	std::vector<double> vals;
	for (int i = 0; i < pixelArray.size(); ++i) {
		/* this method uses the same process as the
		cv2.cvtColor(img, COLOR_RGB2Gray) method */

		//get values for RGB in range [0,1]
		double R = pixelArray[i][0]/255.0;
		double G = pixelArray[i][1]/255.0;
		double B = pixelArray[i][2]/255.0;
		double Y = (R*0.299+G*0.587+B*0.114);//calculate luminance from formula
		vals.insert(vals.begin(), std::round(Y*255.0));//denormalize and save grayscale value
	}
	return Matrix(hgt, wdt, vals);
}

std::vector<unsigned char> IMAGE::get(int i, int j) {
	return pixelArray[i*wdt + j];
}

unsigned char IMAGE::get(int i, int j, int n) {
	return pixelArray[i*wdt + j][n];
}

void IMAGE::set(int i, int j, std::vector<unsigned char> val) {
	pixelArray[i*wdt + j] = val;
}

void IMAGE::set(int i, int j, int n, unsigned char val) {
	pixelArray[i*wdt + j][n] = val;
}

void IMAGE::saveImage() {
	saveImage(filename);
}

void IMAGE::saveImage(std::wstring name) {
	FILE *f = fopen(ws2s(name).c_str(), "wb"); //write file

	if (filetype == L"BMP") {
		//write header
		unsigned char* headerData = new unsigned char[headerSize];
		for (int i = 0; i < header.size(); ++i) {
			headerData[i] = header[i];
		}
		fwrite(headerData, sizeof(unsigned char), headerSize, f);

		//write pixel data
		int pixelCount = 0;
		while (pixelCount < pixelArray.size()) {
			unsigned char* row = new unsigned char[rowSize];
			for (int i = 0; i < rowSize; i += 3) {
				if (padBytes == 0 || (i < pixelsPerRow && pixelCount<pixelArray.size()) ){
					row[i + 0] = pixelArray[pixelCount][0];
					row[i + 1] = pixelArray[pixelCount][1];
					row[i + 2] = pixelArray[pixelCount][2];
					++pixelCount;
				}else {
					while (i<rowSize) {
						row[i+1] = 0x00;
						++i;
					}
				}
			}
			fwrite(row, sizeof(unsigned char), rowSize, f);
		}
	}

	if (filetype == L"PNG") {

	}

	if (filetype == L"JPEG") {

	}

	fclose(f);
}

void IMAGE::printToTextFile() {
	std::wstring str = filename;
	str = str.substr(0,str.find(L"."));
	str.append(L".txt");
	printToTextFile(str);
}

void IMAGE::printToTextFile(std::wstring outputFile) {	
	std::wofstream txt(ws2s(outputFile));//create test.txt and write all lines to the file
	std::wstring line = L"";

	if (filetype == L"JPEG"){}

	if (filetype == L"PNG") {
	}

	if (filetype == L"BMP") {
		//write header
		txt << L"Header\n------------\n";
		for (int i = 0; i < 54; i += 4) {
			//write 4 bytes per line in .txt
			if (i > 0 && i % 4 == 0) {
				txt << line << L"\n";
				line.clear();
			}
			int temp = header[i];
			std::wstring tmpStr = toHex(temp);
			while (tmpStr.size() < 6) { tmpStr.append(L" "); }
			line.append(tmpStr);
			line.append(L" ");
		}
		txt << L"\n";
		line.clear();

		//write pixel Array
		txt << toString() << L"\n";
	}
	txt.close();
}

Vector IMAGE::toVector() {
	std::vector<double> vec;
	for (int i = 0; i < pixelArray.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			vec.push_back(pixelArray[i][j]);
		}
	}
	return vec;
}

Matrix IMAGE::getRedMatrix() {
	std::vector<double> vals;
	for (int i = 0; i < pixelArray.size(); ++i) {
		vals.push_back(pixelArray[i][0]);
	}
	return Matrix(hgt, wdt, vals);
}

Matrix IMAGE::getGreeMatrix() {
	std::vector<double> vals;
	for (int i = 0; i < pixelArray.size(); ++i) {
		vals.push_back(pixelArray[i][1]);
	}
	return Matrix(hgt, wdt, vals);
}

Matrix IMAGE::getBlueMatrix() {
	std::vector<double> vals;
	for (int i = 0; i < pixelArray.size(); ++i) {
		vals.push_back(pixelArray[i][2]);
	}
	return Matrix(hgt, wdt, vals);
}

void IMAGE::maskOutRed() {
	for (int i = 0; i < pixelArray.size(); ++i) {
		pixelArray[i][1] = pixelArray[i][2] = 0;
	}
}

void IMAGE::maskOutGreen() {
	for (int i = 0; i < pixelArray.size(); ++i) {
		pixelArray[i][0] = pixelArray[i][2] = 0;
	}
}

void IMAGE::maskOutBlue() {
	for (int i = 0; i < pixelArray.size(); ++i) {
		pixelArray[i][1] = pixelArray[i][0] = 0;
	}
}

void IMAGE::invertColor() {
	for (int i = 0; i < pixelArray.size(); ++i) {
		pixelArray[i][0] = 255-pixelArray[i][0];
		pixelArray[i][1] = 255-pixelArray[i][1];
		pixelArray[i][2] = 255-pixelArray[i][2];
	}
}

void IMAGE::convertBGRtoRBG() {
	for (int i = 0; i < pixelArray.size(); ++i) {
		unsigned char temp = pixelArray[i][0];
		pixelArray[i][0] = pixelArray[i][2];
		pixelArray[i][2] = temp;
	}
}

void IMAGE::adjustGamma(double gam) {
	for (int i = 0; i < pixelArray.size(); ++i) {
		pixelArray[i][0] = pow(pixelArray[i][0],gam);
		pixelArray[i][1] = pow(pixelArray[i][1],gam);
		pixelArray[i][2] = pow(pixelArray[i][2],gam);
	}
}

void IMAGE::transpose() {
	std::vector<std::vector<unsigned char>> arr;	
	for (int i = 0; i < wdt; ++i) {
		for (int j = 0; j < hgt; ++j) {
			arr.push_back(pixelArray[j*wdt+i]);
		}
	}	
	pixelArray = arr;

	//swap wdt/hgt
	int temp = wdt;
	wdt = hgt;
	hgt = temp;
	updateHeader();
}

void IMAGE::updateHeader() {
	rowSize = calculateRowSize(filetype);
	numberOfBytes = calculateNumberOfBytes(filetype);
	int temp = verticalResolution;
	verticalResolution = horizontalResolution;
	horizontalResolution = temp;
	imageSize = hgt * rowSize;
	numberOfPixels = 3 * hgt*wdt;
	pixelsPerRow = 3 * wdt;
	padBytes = rowSize - pixelsPerRow;
}

std::vector<unsigned char> IMAGE::makeHeader(std::wstring fileType) {
	if (header.size() > 0) {
		header.clear();
	}
	std::vector<unsigned char> hdr;

	if (fileType == L"BMP") {
		//make signature 8-bit header - note: this header is little-endian
		unsigned char dataHeader[54];
		dataHeader[0] = 'B';
		dataHeader[1] = 'M';

		dataHeader[2] = numberOfBytes & 0xFF;
		dataHeader[3] = (numberOfBytes & 0xFF00) >> 8;
		dataHeader[4] = (numberOfBytes & 0xFF0000) >> 16;
		dataHeader[5] = (numberOfBytes & 0xFF000000) >> 24;

		dataHeader[6] = reservedBytes & 0xFF;
		dataHeader[7] = (reservedBytes & 0xFF00) >> 8;
		dataHeader[8] = (reservedBytes & 0xFF0000) >> 16;
		dataHeader[9] = (reservedBytes & 0xFF000000) >> 24;

		dataHeader[10] = headerSize & 0xFF;
		dataHeader[11] = (headerSize & 0xFF00) >> 8;
		dataHeader[12] = (headerSize & 0xFF0000) >> 16;
		dataHeader[13] = (headerSize & 0xFF000000) >> 24;

		//make BITMAPINFOHEADER	
		dataHeader[14] = 40;//the size of this header (40 bytes)
		dataHeader[15] = 0;
		dataHeader[16] = 0;
		dataHeader[17] = 0;

		dataHeader[18] = wdt & 0xFF;//the bitmap width in pixels (signed integer) 
		dataHeader[19] = (wdt & 0xFF00) >> 8;
		dataHeader[20] = (wdt & 0xFF0000) >> 16;
		dataHeader[21] = (wdt & 0xFF000000) >> 24;

		dataHeader[22] = hgt & 0xFF;//the bitmap height in pixels (signed integer) 
		dataHeader[23] = (hgt & 0xFF00) >> 8;
		dataHeader[24] = (hgt & 0xFF0000) >> 16;
		dataHeader[25] = (hgt & 0xFF000000) >> 24;

		dataHeader[26] = 0x01;//the number of color planes (must be 1) 
		dataHeader[27] = 0;

		dataHeader[28] = colorDepth;//the number of bits per pixel, which is the color depth of the image. Typical values are 1, 4, 8, 16, 24 and 32. 
		dataHeader[29] = 0;

		dataHeader[30] = compressionMethod;//the compression method being used. See the next table for a list of possible values 
		dataHeader[31] = 0;
		dataHeader[32] = 0;
		dataHeader[33] = 0;

		dataHeader[34] = imageSize & 0xFF;//the image size. This is the size of the raw bitmap data; a dummy 0 can be given for BI_RGB bitmaps. 
		dataHeader[35] = (imageSize & 0xFF00) >> 8;
		dataHeader[36] = (imageSize & 0xFF0000) >> 16;
		dataHeader[37] = (imageSize & 0xFF000000) >> 24;

		dataHeader[38] = horizontalResolution & 0xFF;//the horizontal resolution of the image. (pixel per metre, signed integer)	
		dataHeader[39] = (horizontalResolution & 0xFF00) >> 8;
		dataHeader[40] = (horizontalResolution & 0xFF0000) >> 16;
		dataHeader[41] = (horizontalResolution & 0xFF000000) >> 24;

		dataHeader[42] = verticalResolution & 0xFF;//the vertical resolution of the image. (pixel per metre, signed integer) 
		dataHeader[43] = (verticalResolution & 0xFF00) >> 8;
		dataHeader[44] = (verticalResolution & 0xFF0000) >> 16;
		dataHeader[45] = (verticalResolution & 0xFF000000) >> 24;

		dataHeader[46] = numberOfColorsInPalatte & 0xFF;//the number of colors in the color palette, or 0 to default to 2^n
		dataHeader[47] = (numberOfColorsInPalatte & 0xFF00) >> 8;
		dataHeader[48] = (numberOfColorsInPalatte & 0xFF0000) >> 16;
		dataHeader[49] = (numberOfColorsInPalatte & 0xFF000000) >> 24;

		dataHeader[50] = 0;//the number of important colors used, or 0 when every color is important; generally ignored 
		dataHeader[51] = 0;
		dataHeader[52] = 0;
		dataHeader[53] = 0;

		for (int i = 0; i < 54; ++i) {
			hdr.push_back(dataHeader[i]);
		}
	}
	return hdr;
}

int IMAGE::calculateRowSize(std::wstring type) {
	if (type == L"BMP") {return (std::floor(((colorDepth*wdt) + 31.0) / 32.0) * 4);}
	return 0;
}

int IMAGE::calculateNumberOfBytes(std::wstring type) {
	if (type == L"BMP") {return headerSize + (rowSize*hgt);}
}

std::wstring IMAGE::getHeadertype() {
	if (filetype == L"BMP") {
		std::wstring str = L"";
		switch (headerSize) {
			case 12+14:  str = L"BITMAPCOREHEADER"; break;
			case 64+14:  str = L"OS22XBITMAPHEADER";break;
			case 16+14:  str = L"OS22XBITMAPHEADER";break;
			case 40+14:  str = L"BITMAPINFOHEADER";break;
			case 52+14:  str = L"BITMAPV2INFOHEADER";break;
			case 56+14:  str = L"BITMAPV3INFOHEADER";break;
			case 10+14:  str = L"BITMAPV4HEADER ";break;
			case 124+14: str = L"BITMAPV5HEADER";break;
		}
		return str;
	}
}

Tensor IMAGE::toTensor() {
	std::vector<double> val;
	for (int i = 0; i < hgt; ++i) {
		for (int j = 0; j < wdt; ++j) {
			for (int n = 0; n < 3; ++n) {
				val.push_back(pixelArray[i*wdt + j][n]);
			}
		}
	}
	std::vector<int> dim;
	dim.push_back(hgt);
	dim.push_back(wdt);
	dim.push_back(3);
	Tensor T(val, dim, false);
	return T;
}

std::wstring IMAGE::toString() {
	std::wstringstream txt;
	std::wstring line = L"";

	if (colorDepth == 24) {
		//write pixel Array
		txt << L"Pixels\n------------\n";
		int counter = 0;
		//read pixel array
		unsigned char dat;
		while (counter<pixelArray.size()) {
			for (int n = 0; n < 3; ++n) {
				unsigned char dat = pixelArray[counter][n];
				int temp = dat;
				std::wstring tmpStr = s2ws(std::to_string(temp));
				while (tmpStr.size() < 4) { tmpStr.append(L" "); }//pad strings to fit them nicely in the file
				line.append(tmpStr);								
			}
			++counter;
			line.append(L"|");
			if (counter%wdt == 0 && counter > 0) { txt << line << L";\n\n";	line.clear(); }
		}
	}
	return txt.str();
}

void IMAGE::display() { std::wcout << toString(); }

//==================================================

bool isBMP(std::wstring filename) {
	FILE *f = fopen(ws2s(filename).c_str(), "rb");

	//check to make sure file has more than 26 bytes in it
	fseek(f, 0, SEEK_END);
	long len = ftell(f);
	fseek(f, 0, SEEK_SET);
	if (len < 26) {
		fclose(f);
		return false;
	}
	unsigned char buf[26];
	fread(buf, 1, 24, f);
	fclose(f);
	if (buf[0] == 'B' && buf[1] == 'M') {
		return true;
	}
	return false;
}

bool isGIF(std::wstring filename) {
	FILE *f = fopen(ws2s(filename).c_str(), "rb");

	//check to make sure file has more than 26 bytes in it
	fseek(f, 0, SEEK_END);
	long len = ftell(f);
	fseek(f, 0, SEEK_SET);
	if (len < 26) {
		fclose(f);
		return false;
	}
	unsigned char buf[24];
	fread(buf, 1, 24, f);
	fclose(f);
	if (buf[0] == 'G' && buf[1] == 'I' && buf[2] == 'F') {
		return true;
	}
	return false;
}

bool isJPEG(std::wstring filename) {
	FILE *f = fopen(ws2s(filename).c_str(), "rb");

	//check to make sure file has more than 26 bytes in it
	fseek(f, 0, SEEK_END);
	long len = ftell(f);
	fseek(f, 0, SEEK_SET);
	if (len < 26) {
		fclose(f);
		return false;
	}
	unsigned char buf[24];
	fread(buf, 1, 24, f);
	fclose(f);
	if (buf[0] == 0xFF && buf[1] == 0xD8 && buf[2] == 0xFF && buf[3] == 0xE0 && buf[6] == 'J'
		&& buf[7] == 'F' && buf[8] == 'I' && buf[9] == 'F') {
		return true;
	}
	return false;
}

bool isPNG(std::wstring filename) {
	FILE *f = fopen(ws2s(filename).c_str(), "rb");

	//check to make sure file has more than 26 bytes in it
	fseek(f, 0, SEEK_END);
	long len = ftell(f);
	fseek(f, 0, SEEK_SET);
	if (len < 26) {
		fclose(f);
		return false;
	}

	unsigned char buf[8];
	fread(buf, 1, 8, f);
	fclose(f);
	if (buf[0] == 0x89 && buf[1] == 'P' && buf[2] == 'N' && buf[3] == 'G'){
		return true;
	}
	return false;
}

void fileToText(std::wstring filename) {
	FILE *f = fopen(ws2s(filename).c_str(), "rb");
	std::wofstream txt(L"HexDump.txt");//create HexDump.txt and write all lines to the file
	std::wstring line = L"";

	int counter = 0;
	while (!feof(f)) {		
		unsigned char line[4];
		while (counter < 4) {
			unsigned char buf = 0;
			fread(&buf, sizeof(unsigned char), 1, f);
			line[counter] = buf;

			int temp = buf;
			std::wstring tmpStr = toHex(temp);
			txt << tmpStr << L" ";
			
			++counter;
		}
		txt << L"\n";
		counter = 0;
	}
	fclose(f);
}

void BMPtoText(std::wstring filename) { BMPtoText(filename, L"HexDump.txt"); }

void BMPtoText(std::wstring filename, std::wstring outputName) {
	FILE *f = fopen(ws2s(filename).c_str(), "rb");
	std::wofstream txt(outputName);
	std::wstring line = L"";

	//read the 54-byte header
	unsigned char header[54];
	fread(header, sizeof(unsigned char), 54, f);	 

	int numberOfBytes = (header[5] << 24) ^ (header[4] << 16) ^ (header[3] << 8) ^ header[2];//read number of bytes in file
	int wdt = (header[21] << 24) ^ (header[20] << 16) ^ (header[19] << 8) ^ header[18];//gives image height 
	int hgt = (header[25] << 24) ^ (header[24] << 16) ^ (header[23] << 8) ^ header[22];//gives image width
	int horizontalResolution = (header[41] << 24) ^ (header[40] << 16) ^ (header[39] << 8) ^ header[38];
	int verticalResolution = (header[45] << 24) ^ (header[44] << 16) ^ (header[43] << 8) ^ header[42];
	int imageSize = (header[37] << 24) ^ (header[36] << 16) ^ (header[35] << 8) ^ header[34];
	int colorDepth = (header[27] << 8) ^ header[28]; //read color depth to determine color array
	int compressionMethod = (header[33] << 24) ^ (header[32] << 16) ^ (header[31] << 8) ^ header[30];//tells what method of compression is used (0 = uncompressed) 
	int numberOfColorsInPalatte = (header[49] << 24) ^ (header[48] << 16) ^ (header[47] << 8) ^ header[46];//gives number of colors in color palatte, 0 = 2^n (default)
	int numberOfPixels = 3 * wdt * hgt; //number of bytes in pixel array
	int rowSize = std::floor(((colorDepth*wdt) + 31.0) / 32.0) * 4;//number of bytes needed to store one row
	int pixelsPerRow = wdt * 3;
	int padBytes = rowSize - pixelsPerRow;//if the width is not divisible by 4, then it is filled in by padding (ie 10 mod 4 = 2 padded entries)


	//write header
	txt << L"Header\n------------\n";
	for (int i = 0; i < 54; i += 4) {
		//write 4 bytes per line in .txt
		if (i > 0 && i % 4 == 0) {
			txt << line << L"\n";
			line.clear();
		}
		int temp = header[i];
		std::wstring tmpStr = toHex(temp);
		while (tmpStr.size() < 6) { tmpStr.append(L" "); }
		line.append(tmpStr);
		line.append(L" ");
	}
	txt << L"\n";
	line.clear();

	if (colorDepth == 24) {
		//write pixel Array
		txt << L"Pixels\n------------\n";
		int counter = 0;
		//read pixel array
		unsigned char dat;
		while (!feof(f)) {
			fread(&dat, sizeof(unsigned char), 1, f);
			int temp = dat;
			std::wstring tmpStr = s2ws(std::to_string(temp));
			while (tmpStr.size() < 4) { tmpStr.append(L" "); }//pad strings to fit them nicely in the file
			line.append(tmpStr);
			++counter;
			if (counter > 0 && counter % 3 == 0) { line.append(L"|"); }
			if (counter%rowSize == 0 && counter > 0) { txt << line << L";\n\n";	line.clear(); }
		}
	}
	fclose(f);
	txt.close();
}

unsigned long* make_crc_table(){//make table for calculating CRC code used in PNG file format
	unsigned long crc_table[256];
	unsigned long c;
	int n, k;

	for (n = 0; n < 256; n++) {
		c = (unsigned long)n;
		for (k = 0; k < 8; k++) {
			if (c & 1)
				c = 0xedb88320L ^ (c >> 1);
			else
				c = c >> 1;
		}
		crc_table[n] = c;
	}
	return crc_table;
}

Matrix identityFilter(unsigned int sz) {
	if (sz % 2 == 0) { return Matrix(); }//filters must be of odd size
	std::vector<double> vals(sz*sz, 0);
	vals[floor((sz*sz) / 2)] = 1;
	Matrix Mat(sz, sz, vals);
	return Mat;
}

Matrix sharpeningFilter(unsigned int sz) { 
	if (sz % 2 == 0) { return Matrix(); }//filters must be of odd size
	std::vector<double> vals;
	if (sz == 3) {
		double vals2[9] = {
			 0, -1,  0,
			-1,  5, -1,
			 0, -1,  0
		};
		vals = toSTLVector(vals2, 9);
	}
	Matrix Mat(sz, sz, vals);
	return Mat;
}

Matrix edgeDetectorFilter(unsigned int sz) {
	if (sz % 2 == 0) { return Matrix(); }//filters must be of odd size
	std::vector<double> vals;
	if (sz == 3) {
		double vals2[9] = {
			-1, -1, -1,
			-1,  8, -1,
			-1, -1, -1
		};
		vals = toSTLVector(vals2, 9);
	}
	Matrix Mat(sz, sz, vals);
	return Mat;
}

Matrix boxFilter(unsigned int sz) {
	if (sz % 2 == 0) { return Matrix(); }//filters must be of odd size
	std::vector<double> vals(sz*sz,1);
	vals[floor((sz*sz) / 2) + 1] = 1;
	Matrix Mat(sz, sz, vals);
	return Mat;
}

Matrix GaussianKernel1D(int sz) { return GaussianKernel1D(sz,1); }

Matrix GaussianKernel1D(int sz,double sigma) {
	double sigmaX = sigma > 0 ? sigma : ((sz - 1)*0.5 - 1)*0.3 + 0.8;
	double scale2X = -0.5 / (sigmaX*sigmaX);
	double sum = 0;

	std::vector<double> cd(sz);

	const int SMALL_GAUSSIAN_SIZE = 7;
	static const float small_gaussian_tab[][SMALL_GAUSSIAN_SIZE] = {
		{1.f},
		{0.25f, 0.5f, 0.25f},
		{0.0625f, 0.25f, 0.375f, 0.25f, 0.0625f},
		{0.03125f, 0.109375f, 0.21875f, 0.28125f, 0.21875f, 0.109375f, 0.03125f}
	};

	const float* fixed_kernel = sz % 2 == 1 && sz <= SMALL_GAUSSIAN_SIZE && sigma <= 0 ? small_gaussian_tab[sz >> 1] : 0;
	int i;
	for (i = 0; i < sz; i++) {
		double x = i - (sz - 1)*0.5;
		double t = fixed_kernel ? (double)fixed_kernel[i] : std::exp(scale2X*x*x);
		cd[i] = t;
		sum += cd[i];
	}
	sum = 1.0 / sum;
	for (i = 0; i < sz; i++) { cd[i] *= sum; }
	return Matrix(sz, 1, cd);
}

Matrix GaussianKernel2D(int sz) { return GaussianKernel2D(sz, 0.3*((sz - 1)*0.5 - 1) + 0.8); }

Matrix GaussianKernel2D(int sz, double sigma) {
	Matrix Mat1 = GaussianKernel1D(sz, sigma);
	Matrix MatT = Mat1.transpose();
	return Mat1.multiply(Mat1, MatT);
}

Matrix gradientFilter(std::wstring type) { return gradientFilter(type,1); }

Matrix gradientFilter(std::wstring type, unsigned int wrt) {
	int sz1, sz2;
	double* vals = NULL;

	if (type == L"finite difference" || type == L"") {
		double vals2[3*1] = {0,-1,1};
		sz1 = 3;
		sz2 = 1;
		vals = new double[sz1*sz2];
		vals = vals2;			
	}
	if (type == L"Sobel") {
		double vals2[3 * 3] = {
			 -1, 0, 1,
			 -2, 0, 2,
			 -1, 0, 1
		};
		sz1 = 3;
		sz2 = 3;
		vals = new double[sz1*sz2];
		vals = vals2;
	}
	if (type == L"Schurr") {
		double vals2[3 * 3] = {
			 -3, 0, 3,
			-10, 0,10,
			 -3, 0, 3
		};
		sz1 = 3;
		sz2 = 3;
		vals = new double[sz1*sz2];
		vals = vals2;
	}
	if (type == L"Laplacian") {
		double vals2[3 * 3] = {
			0, 1, 0,
			1,-4, 1,
			0, 1, 0
		};
		sz1 = 3;
		sz2 = 3;
		vals = new double[sz1*sz2];
		vals = vals2;
	}
	std::vector<double> vals3 = toSTLVector(vals, sz1*sz2);
	Matrix Mat(sz1, sz2, vals3);
	if (wrt > 0) { Mat = Mat.transpose(); }
	return Mat;
}

Matrix testMatrix() {
	double vals[15 * 12] = {
	255,255,255,   255,255,255,  255,255,255,  255,255,255,   255,255,255,
	  0,  0,  0,     0,  0,  0,    0,  0,  0,    0,  0,  0,     0,  0,  0,
	255,255,255,   255,255,255,  255,255,255,  255,255,255,   255,255,255,

	255,255,255,     0,255,255,    0,255,255,    0,255,255,    255,255,255,
	255,255,255,     0,255,255,    0,255,255,    0,255,255,    255,255,255,
	255,255,255,     0,255,255,    0,255,255,    0,255,255,    255,255,255,

	255,255,255,     0,255,255,    0,255,255,    0,255,255,    255,255,255,
	255,255,255,     0,255,255,    0,255,255,    0,255,255,    255,255,255,
	255,255,255,     0,255,255,    0,255,255,    0,255,255,    255,255,255,

	255,255,255,   255,255,255,  255,255,255,   255,255,255,   255,255,255,
	  0,  0,  0,     0,  0,  0,    0,  0,  0,    0,  0,  0,     0,  0,  0,
	255,255,255,   255,255,255,  255,255,255,   255,255,255,   255,255,255
	};
	Matrix Mat(12, 15, toSTLVector(vals, 15 * 12));
	return Mat;
}