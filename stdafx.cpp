#pragma once
#include "stdafx.h"

//ENVIRONMENT VARIABLES
//=====================
unsigned int precision = 4;
double SYSTEM_EPSILON = 0.00000001;											//value for which any smaller value is considered to be 0

//CONSTANTS
//=========
const double PI = 3.14159265358979323846264338328;							//Pi to as many digits as is useful
const double EMconstant = 0.57721566490153286060651209008240243;			//Euler-Mascheroni constant
const double EulersNumber = 2.71828182845904523536028747135266249;			//Euler's number, e
const double SoldnersConstant = 1.45136923488;
const double AperysConstant = 1.20205690315959428540;
const double R_constant = 1.097373156853955e+7;								//Rydberg's constant, R
const int speedOfLight = 299792458;                                     	//Speed of light in m/s
const double speedOfSound = 340.29;											//speed of sound in m/s
const double gravitationalConstant = 6.67408e-11;							//Newton's gravitational constant, G in units of Nm^2/kg^2
const double g_constant = 9.80665;											//Earth's gravitational constant (g) in units of m/s^2
const double molarVolumeIdealGas = 2.271094713e-2;							//at T = 273.15 Kelvin and p = 100 kPascals
const double idealGasConstantATMs = 0.082057338;                            //ideal gas constant in atm*Liters/K*mol
const double idealGasConstantKiloPascals = 8.314459848;                     //ideal gas constant in kPascals*Liters/K*mol also Joules/K*mol
double EPSILON = 1.0 / DBL_MAX;												//smallest values for double data type, set as the "epsilon" of our system by default
const double PHI = (1.0 + sqrt(5)) / 2;										//PHI, half of Fibonacci number calculation
const double PSI = (1.0 - sqrt(5)) / 2;										//PSI, other half of Fibonnaci number calculation

//		unit conversions
//=================================
const double PascalsPerATM = 101325;
const double ATMsPerPascal = 9.86923e-6;
const double mmHgPerPascal = 0.0075006156130264;
const double PascalsPermmHg = 133.32239;

//constants from particle physics
//===============================
const double PlancksConstantJoules = 6.626070040e-34;
const double PlancksConstanteVs = 4.135667662e-15;
const double vacuumPermitivitty = 8.854187817e-12;
const double vacuumPermeability = 4 * PI;
const double elementaryCharge = 1.6021766208e-19;							//in Coulombs
const double electronMasskg = 9.10938356e-31;								//electron mass in kg
const double electronMassMeV = 0.5109989461;
const double electronMassAMU = 0.000548579909070;
const double protonMasskg = 1.672621898e-27;								//proton mass in kg
const double protonMassMeV = 938.2720813;
const double protonMassAMU = 1.007276466879;
const double neutronMasskg = 1.674927471e-27;								//neutron mass in kg
const double neutronMassMeV = 939.5654133;
const double neutronMassAMU = 1.00866491588;
const double AvogadrosNumber = 6.022140857e23;
const double BoltzmannsConstant = 1.38064852e-23;
const double StefanBoltzmannConstant = 5.670367e-8;							//units = W/m^2 K^4
const double WeinDisplacementConstantmm = 2.8977729;						//units = mm/K
const double WeinDisplacementConstantGHz = 58.789238;						//units = GHz/K
const double atomicMassUnit = 1.66053892173e-27;							//i.e. 1 AMU = C-12 /12

//misc constants
//==============
const double HubbleConstantkm = 69.3;										//units = km/s/Mpc
const double HubbleConstantsec = 2.25e-18;									//units = sec^-1

//audio/electrical engineering constants
//=======================================
const double defaultAudioSamplingRate = 44100.0;							//default audio sample rate = 44.1k
const double HammingAlpha = 0.54;											// for Hamming window