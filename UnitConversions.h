#pragma once
#include "stdafx.h"

double ATMToPascal(double x){}
double CelsiusToFahrenheit(double x) { return ((9/5.0)*x) + 32; }
double CelsiusToKelvin(double x) {return x + 273.15;}
double FahrenheitToCelsius(double x) { return (5.0/9)*(x - 32); }
double KelvinToCelsius(double x) { return x - 273.15; }
double mmHgToPascal(double x) { return (x*PascalsPermmHg); }
double PascalTommHg(double x) { return (x*mmHgPerPascal); }