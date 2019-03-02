#pragma once
#include "stdafx.h"

ODE::ODE() {}
ODE::~ODE() {}
ODE::ODE(std::wstring a, std::wstring b) {
	lhs = a;
	rhs = b;
}
//========================================
PDE::PDE() {}
PDE::~PDE() {}
PDE::PDE(std::wstring a, std::wstring b) {
	lhs = a;
	rhs = b;
}