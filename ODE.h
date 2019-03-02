#pragma once

class ODE {
public:
	std::wstring lhs;
	std::wstring rhs;
	ODE();
	~ODE();
	ODE(std::wstring a, std::wstring b);
};

class PDE{
public:
	std::wstring lhs;
	std::wstring rhs;
	PDE();
	~PDE();
	PDE(std::wstring a, std::wstring b);
};