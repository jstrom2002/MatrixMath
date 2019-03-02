#ifndef EXTERN_H
#define EXTERN_H
#pragma once
//Classes
//=======
//local declaration of extern objects
extern class Circle;
extern class ComplexFunction;
extern class ComplexInfiniteSeries;
extern class ComplexMatrix;
extern class ComplexNumber;
extern class ComplexPolynomial;
extern class ComplexPolynomialFraction;
extern class ComplexVector;
extern class Curve;
extern class DirectedGraph;
extern class DualNumber;
extern class DualNumberMatrix;
extern class DualNumberPolynomial;
extern class DualQuaternion;
extern class Fraction;
extern class Function;
extern class FunctionMatrix;
extern class GeometricObject;
extern class Group;
extern class Hyperplane;
extern class InfiniteSeries;
extern class Line;
extern class LineSegment;
extern class Matrix;
extern class Octonian;
extern class ODE;
extern class ParametricFunction;
extern class PDE;
extern class Polynomial;
extern class PolynomialFraction;
extern class Quaternion;
extern class SparseMatrix;
extern class Sphere;
extern class Tensor;
extern class Triangle;
extern class UndirectedGraph;
extern class Vector;
extern class VectorValuedFunction;

//Vectors of objects to save/display
//==================================
extern std::vector<ComplexFunction> ComplexFunctionContainer;
extern std::vector<ComplexInfiniteSeries> ComplexInfiniteSeriesContainer;
extern std::vector<ComplexMatrix> ComplexMatrixContainer;
extern std::vector<ComplexNumber> ComplexNumberContainer;
extern std::vector<ComplexPolynomial> ComplexPolynomialContainer;
extern std::vector<ComplexPolynomialFraction> ComplexPolynomialFractionContainer;
extern std::vector<ComplexVector> ComplexVectorContainer;
extern std::vector<DirectedGraph> DirectedGraphContainer;
extern std::vector<DualNumber> DualNumberContainer;
extern std::vector<DualNumberMatrix> DualNumberMatrixContainer;
extern std::vector<DualNumberPolynomial> DualNumberPolynomialContainer;
extern std::vector<DualQuaternion> DualQuaternionContainer;
extern std::vector<Function> FunctionContainer;
extern std::vector<FunctionMatrix> FunctionMatrixContainer;
extern std::vector<Group> GroupContainer;
extern std::vector<InfiniteSeries> InfiniteSeriesContainer;
extern std::vector<Matrix> MatrixContainer;
extern std::vector<Octonian> OctonianContainer;
extern std::vector<ODE> ODEContainer;
extern std::vector<ParametricFunction> ParametricFunctionContainer;
extern std::vector<PDE> PDEContainer;
extern std::vector<Polynomial> PolynomialContainer;
extern std::vector<PolynomialFraction> PolynomialFractionContainer;
extern std::vector<Quaternion> QuaternionContainer;
extern std::vector<SparseMatrix> SparseMatrixContainer;
extern std::vector<Tensor> TensorContainer;
extern std::vector<UndirectedGraph> UndirectedGraphContainer;
extern std::vector<Vector> VectorContainer;
extern std::vector<VectorValuedFunction> VectorValuedFunctionContainer;

//program environment variables
//=============================
extern bool changeScreen;
extern std::wstring filename;	//holds the name of the output file where the data will be saved to or loaded from
extern bool first;				//variable indicating if this is the first time the program boots up
extern bool fullscreen;			//variable determining whether or not the main console window is fullscreen or windowed
extern std::wstring hostname;    //name of internet provider host
extern std::vector<std::wstring> IPs;	//IP addresses in all formats(IPv6, IPv4, etc.)
extern std::wstring productKey;			//the key for the program, used to unlock it when it runs
extern std::wstring password;	//user-defined password
extern std::wstring path;		//directory path of where the .exe is running from
extern unsigned int precision;	//length of the mantissa to be dispalyed
extern bool running;			//variable determining whether or not the program ought to still be running
extern std::wstring serialnumber;//contains the program's serial number


extern unsigned int precision;
extern double SYSTEM_EPSILON;
extern const double PI;//pi to as many digits as is useful
extern const double EMconstant;//Euler-Mascheroni constant
extern const double EulersNumber;//Euler's number, e
extern const double SoldnersConstant;
extern const double AperysConstant;
extern const double R_constant;//Rydberg's constant, R
extern const int speedOfLight;//Speed of light in m/s
extern const double speedOfSound;//speed of sound in m/s
extern const double gravitationalConstant;//Newton's gravitational constant, G in units of Nm^2/kg^2
extern const double g_constant;//Earth's gravitational constant (g) in units of m/s^2
extern const double molarVolumeIdealGas;//at T = 273.15 Kelvin and p = 100 kPascals
extern const double idealGasConstantATMs;//ideal gas constant in atm*Liters/K*mol
extern const double idealGasConstantKiloPascals;//ideal gas constant in kPascals*Liters/K*mol also Joules/K*mol
extern double EPSILON;
extern const double PHI;
extern const double PSI;
extern const double PascalsPerATM;
extern const double ATMsPerPascal;
extern const double mmHgPerPascal;
extern const double PascalsPermmHg;
extern const double PlancksConstantJoules;
extern const double PlancksConstanteVs;
extern const double vacuumPermitivitty;
extern const double vacuumPermeability;
extern const double elementaryCharge;
extern const double electronMasskg;
extern const double electronMassMeV;
extern const double electronMassAMU;
extern const double protonMasskg;
extern const double protonMassMeV;
extern const double protonMassAMU;
extern const double neutronMasskg;
extern const double neutronMassMeV;
extern const double neutronMassAMU;
extern const double AvogadrosNumber;
extern const double BoltzmannsConstant;
extern const double StefanBoltzmannConstant;
extern const double WeinDisplacementConstantmm;
extern const double WeinDisplacementConstantGHz;
extern const double atomicMassUnit;
extern const double HubbleConstantkm;
extern const double HubbleConstantsec;
extern const double defaultAudioSamplingRate;
extern const double HammingAlpha;

//windows enviroment variables
//============================
extern HDC hdc;				//parent device context
extern MSG messages;        //parent message container
extern HWND hwnd;           //main window
extern HWND outputPanel;    //output text panel
extern HWND hEdit;			//user input bar
//screen size/variables:
extern int fullscreenWidth;
extern int fullscreenHeight;
extern int colourBits;
extern int refreshRate;

//IO strings
//==========
extern std::wstring input;   //holds string from input text box
extern std::wstring output;	//holds display to be shown to user
#endif