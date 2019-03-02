#ifndef FILEIO_H
#define FILEIO_H
#pragma once
#include "stdafx.h"

std::vector<double> readTxtFile(std::string filename);
std::vector<double> readTxtFile(std::wstring filename);

#endif
