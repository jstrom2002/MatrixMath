#pragma once
#include "FileIO.h"

std::vector<double> readTxtFile(std::string filename) {
	return readTxtFile(s2ws(filename));
}

std::vector<double> readTxtFile(std::wstring filename) {
	std::ifstream file(filename);
	std::vector<double> vec;

	std::string str;
	while (std::getline(file, str)) {
		std::vector<double> tempVec = ParseNInputNumbers(str);
		for (int i = 0; i < tempVec.size(); ++i) {
			vec.push_back(tempVec[i]);
		}
	}
	return vec;
}