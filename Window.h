#ifndef WINDOW_H
#define WINDOW_H
#include "stdafx.h"

std::vector<ComplexNumber> BartlettWindow(std::vector<ComplexNumber> data);
std::vector<ComplexNumber> GaussianWindow(std::vector<ComplexNumber> data);
std::vector<ComplexNumber> HammingWindow(std::vector<ComplexNumber> data);
std::vector<ComplexNumber> WelchWindow(std::vector<ComplexNumber> data);

#pragma once
#endif