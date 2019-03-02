#pragma once
#include "stdafx.h"

std::wstring PeriodicTable(std::wstring Command, int n) {
	/*A matrix of rows containing elements and their
	Name, Number, Grp, Period, Formula, Charge, Mass, Special information*/
	std::wstring Name, Formula, Charge, Special;
	int Number, Grp, Period;
	double Mass;
	if ((Command == L"1") || (Command == L"hydrogen") || (Command == L"H"))
	{
		Name = L"hydrogen";
		Number = 1;
		Grp = 1;
		Period = 1;
		Formula = L"H";
		Charge = L"1+";
		Mass = 1.01;
		Special = L"Diatomic";
	}
	else if ((Command == L"2") || (Command == L"helium") || (Command == L"He"))
	{
		Name = L"helium";
		Number = 2;
		Grp = 18;
		Period = 1;
		Formula = L"He";
		Charge = L"0";
		Mass = 4.00;
		Special = L"None";
	}
	else if ((Command == L"3") || (Command == L"lithium") || (Command == L"Li"))
	{
		Name = L"lithium";
		Number = 3;
		Grp = 1;
		Period = 2;
		Formula = L"Li";
		Charge = L"1+";
		Mass = 6.94;
		Special = L"None";
	}
	else if ((Command == L"4") || (Command == L"beryllium") || (Command == L"Be"))
	{
		Name = L"beryllium";
		Number = 4;
		Grp = 2;
		Period = 2;
		Formula = L"Be";
		Charge = L"2+";
		Mass = 9.01;
		Special = L"None";
	}
	else if ((Command == L"5") || (Command == L"boron") || (Command == L"B"))
	{
		Name = L"boron";
		Number = 5;
		Grp = 13;
		Period = 2;
		Formula = L"B";
		Charge = L"3+";
		Mass = 10.81;
		Special = L"None";
	}
	else if ((Command == L"6") || (Command == L"carbon") || (Command == L"C"))
	{
		Name = L"carbon";
		Number = 6;
		Grp = 14;
		Period = 2;
		Formula = L"C";
		Charge = L"4+";
		Mass = 12.01;
		Special = L"None";
	}
	else if ((Command == L"7") || (Command == L"nitrogen") || (Command == L"N"))
	{
		Name = L"nitrogen";
		Number = 7;
		Grp = 15;
		Period = 2;
		Formula = L"N";
		Charge = L"3-";
		Mass = 14.01;
		Special = L"Diatomic";

	}
	else if ((Command == L"8") || (Command == L"oxygen") || (Command == L"O"))
	{
		Name = L"oxygen";
		Number = 8;
		Grp = 16;
		Period = 2;
		Formula = L"O";
		Charge = L"2-";
		Mass = 16.00;
		Special = L"Diatomic";

	}
	else if ((Command == L"9") || (Command == L"flourine") || (Command == L"F"))
	{
		Name = L"fluorine";
		Number = 9;
		Grp = 17;
		Period = 2;
		Formula = L"F";
		Charge = L"1-";
		Mass = 19.00;
		Special = L"Diatomic";

	}
	else if ((Command == L"10") || (Command == L"neon") || (Command == L"Ne"))
	{
		Name = L"neon";
		Number = 10;
		Grp = 18;
		Period = 2;
		Formula = L"Ne";
		Charge = L"0";
		Mass = 20.18;
		Special = L"Noble Gas";

	}
	else if ((Command == L"11") || (Command == L"sodium") || (Command == L"Na"))
	{
		Name = L"sodium";
		Number = 11;
		Grp = 1;
		Period = 3;
		Formula = L"Na";
		Charge = L"1+";
		Mass = 22.99;
		Special = L"None";

	}
	else if ((Command == L"12") || (Command == L"magnesium") || (Command == L"Mg"))
	{
		Name = L"magnesium";
		Number = 12;
		Grp = 2;
		Period = 3;
		Formula = L"Mg";
		Charge = L"2+";
		Mass = 24.31;
		Special = L"None";

	}
	else if ((Command == L"13") || (Command == L"aluminum") || (Command == L"Al"))
	{
		Name = L"aluminum";
		Number = 13;
		Grp = 13;
		Period = 3;
		Formula = L"Al";
		Charge = L"3+";
		Mass = 26.98;
		Special = L"None";

	}
	else if ((Command == L"14") || (Command == L"silicon") || (Command == L"Si"))
	{
		Name = L"silicon";
		Number = 14;
		Grp = 14;
		Period = 3;
		Formula = L"Si";
		Charge = L"4+";
		Mass = 28.09;
		Special = L"None";

	}
	else if ((Command == L"15") || (Command == L"phosphorous") || (Command == L"P"))
	{
		Name = L"phosphorus";
		Number = 15;
		Grp = 15;
		Period = 3;
		Formula = L"P";
		Charge = L"3-";
		Mass = 30.97;
		Special = L"Diatomic";

	}
	else if ((Command == L"16") || (Command == L"sulfur") || (Command == L"S"))
	{
		Name = L"sulfur";
		Number = 16;
		Grp = 16;
		Period = 3;
		Formula = L"S";
		Charge = L"2-";
		Mass = 32.06;
		Special = L"Diatomic";

	}
	else if ((Command == L"17") || (Command == L"chlorine") || (Command == L"Cl"))
	{
		Name = L"chlorine";
		Number = 17;
		Grp = 17;
		Period = 3;
		Formula = L"Cl";
		Charge = L"1-";
		Mass = 35.45;
		Special = L"Diatomic";

	}
	else if ((Command == L"18") || (Command == L"argon") || (Command == L"Ar"))
	{
		Name = L"argon";
		Number = 18;
		Grp = 18;
		Period = 3;
		Formula = L"Ar";
		Charge = L"0";
		Mass = 39.95;
		Special = L"Noble Gas";

	}
	else if ((Command == L"19") || (Command == L"potassium") || (Command == L"K"))
	{
		Name = L"potassium";
		Number = 19;
		Grp = 1;
		Period = 4;
		Formula = L"K";
		Charge = L"1+";
		Mass = 39.10;
		Special = L"None";

	}
	else if ((Command == L"20") || (Command == L"calcium") || (Command == L"Ca"))
	{
		Name = L"calcium";
		Number = 20;
		Grp = 2;
		Period = 4;
		Formula = L"Ca";
		Charge = L"2+";
		Mass = 40.08;
		Special = L"None";

	}
	else if ((Command == L"21") || (Command == L"scandium") || (Command == L"Sc"))
	{
		Name = L"scandium";
		Number = 21;
		Grp = 3;
		Period = 4;
		Formula = L"Sc";
		Charge = L"3+";
		Mass = 44.96;
		Special = L"None";

	}
	else if ((Command == L"22") || (Command == L"titanium") || (Command == L"Ti"))
	{
		Name = L"titanium";
		Number = 22;
		Grp = 4;
		Period = 4;
		Formula = L"Ti";
		Charge = L"4+ \\ 3+";
		Mass = 47.88;
		Special = L"None";

	}
	else if ((Command == L"23") || (Command == L"vanadium") || (Command == L"V"))
	{
		Name = L"vanadium";
		Number = 23;
		Grp = 5;
		Period = 4;
		Formula = L"V";
		Charge = L"5+ \\ 4+";
		Mass = 50.94;
		Special = L"None";

	}
	else if ((Command == L"24") || (Command == L"chromium") || (Command == L"Cr"))
	{
		Name = L"chromium";
		Number = 24;
		Grp = 6;
		Period = 4;
		Formula = L"Cr";
		Charge = L"3+ \\ 2+";
		Mass = 52.00;
		Special = L"None";

	}
	else if ((Command == L"25") || (Command == L"manganese") || (Command == L"Mn"))
	{
		Name = L"manganese";
		Number = 25;
		Grp = 7;
		Period = 4;
		Formula = L"Mn";
		Charge = L"2+ \\ 4+";
		Mass = 54.94;
		Special = L"None";

	}
	else if ((Command == L"26") || (Command == L"iron") || (Command == L"Fe"))
	{
		Name = L"iron";
		Number = 26;
		Grp = 8;
		Period = 4;
		Formula = L"Fe";
		Charge = L"3+ \\ 2+";
		Mass = 55.85;
		Special = L"None";

	}
	else if ((Command == L"27") || (Command == L"cobalt") || (Command == L"Co"))
	{
		Name = L"cobalt";
		Number = 27;
		Grp = 9;
		Period = 4;
		Formula = L"Co";
		Charge = L"2+ \\ 3+";
		Mass = 58.93;
		Special = L"None";

	}
	else if ((Command == L"28") || (Command == L"nickel") || (Command == L"Ni"))
	{
		Name = L"nickel";
		Number = 28;
		Grp = 10;
		Period = 4;
		Formula = L"Ni";
		Charge = L"2+ \\ 3+";
		Mass = 58.69;
		Special = L"None";

	}
	else if ((Command == L"29") || (Command == L"copper") || (Command == L"Cu"))
	{
		Name = L"copper";
		Number = 29;
		Grp = 11;
		Period = 4;
		Formula = L"Cu";
		Charge = L"2+ \\ 1+";
		Mass = 63.55;
		Special = L"None";

	}
	else if ((Command == L"30") || (Command == L"zinc") || (Command == L"Zn"))
	{
		Name = L"zinc";
		Number = 30;
		Grp = 12;
		Period = 4;
		Formula = L"Zn";
		Charge = L"2+";
		Mass = 65.38;
		Special = L"None";

	}
	else if ((Command == L"31") || (Command == L"gallium") || (Command == L"Ga"))
	{
		Name = L"gallium";
		Number = 31;
		Grp = 13;
		Period = 4;
		Formula = L"Ga";
		Charge = L"3+";
		Mass = 69.72;
		Special = L"None";

	}
	else if ((Command == L"32") || (Command == L"germanium") || (Command == L"Ge"))
	{
		Name = L"germanium";
		Number = 32;
		Grp = 14;
		Period = 4;
		Formula = L"Ge";
		Charge = L"4+";
		Mass = 72.61;
		Special = L"None";

	}
	else if ((Command == L"33") || (Command == L"arsonic") || (Command == L"As"))
	{
		Name = L"arsonic";
		Number = 33;
		Grp = 15;
		Period = 4;
		Formula = L"As";
		Charge = L"3-";
		Mass = 74.92;
		Special = L"None";

	}
	else if ((Command == L"34") || (Command == L"selenium") || (Command == L"Se"))
	{
		Name = L"selenium";
		Number = 34;
		Grp = 16;
		Period = 4;
		Formula = L"Se";
		Charge = L"2-";
		Mass = 78.96;
		Special = L"None";

	}
	else if ((Command == L"35") || (Command == L"bromine") || (Command == L"Br"))
	{
		Name = L"bromine";
		Number = 35;
		Grp = 17;
		Period = 4;
		Formula = L"Br";
		Charge = L"1-";
		Mass = 79.90;
		Special = L"Diatomic";

	}
	else if ((Command == L"36") || (Command == L"krypton") || (Command == L"Kr"))
	{
		Name = L"krypton";
		Number = 36;
		Grp = 18;
		Period = 4;
		Formula = L"Kr";
		Charge = L"0";
		Mass = 83.80;
		Special = L"Noble Gas";

	}
	else if ((Command == L"37") || (Command == L"rubidium") || (Command == L"Rb"))
	{
		Name = L"rubidium";
		Number = 37;
		Grp = 1;
		Period = 5;
		Formula = L"Rb";
		Charge = L"1+";
		Mass = 85.47;
		Special = L"None";

	}
	else if ((Command == L"38") || (Command == L"stronthum") || (Command == L"Sr"))
	{
		Name = L"stronthum";
		Number = 38;
		Grp = 2;
		Period = 5;
		Formula = L"Sr";
		Charge = L"2+";
		Mass = 87.62;
		Special = L"None";

	}
	else if ((Command == L"39") || (Command == L"ythrium") || (Command == L"Y"))
	{
		Name = L"ythrium";
		Number = 39;
		Grp = 3;
		Period = 5;
		Formula = L"Y";
		Charge = L"3+";
		Mass = 88.91;
		Special = L"None";

	}
	else if ((Command == L"40") || (Command == L"zirconium") || (Command == L"Zr"))
	{
		Name = L"zirconium";
		Number = 40;
		Grp = 4;
		Period = 5;
		Formula = L"Zr";
		Charge = L"4+";
		Mass = 91.22;
		Special = L"None";

	}
	else if ((Command == L"41") || (Command == L"niobium") || (Command == L"Nb"))
	{
		Name = L"niobium";
		Number = 41;
		Grp = 5;
		Period = 5;
		Formula = L"Nb";
		Charge = L"5+ \\ 3+";
		Mass = 92.91;
		Special = L"None";

	}
	else if ((Command == L"42") || (Command == L"molybdenum") || (Command == L"Mo"))
	{
		Name = L"molybdenum";
		Number = 42;
		Grp = 6;
		Period = 5;
		Formula = L"Mo";
		Charge = L"6+";
		Mass = 95.94;
		Special = L"None";

	}
	else if ((Command == L"43") || (Command == L"techenium") || (Command == L"Tc"))
	{
		Name = L"techenium";
		Number = 43;
		Grp = 7;
		Period = 5;
		Formula = L"Tc";
		Charge = L"7+";
		Mass = 98.91;
		Special = L"None";

	}
	else if ((Command == L"44") || (Command == L"ruthenium") || (Command == L"Ru"))
	{
		Name = L"ruthenium";
		Number = 44;
		Grp = 8;
		Period = 5;
		Formula = L"Ru";
		Charge = L"3+ \\ 4+";
		Mass = 101.07;
		Special = L"None";

	}
	else if ((Command == L"45") || (Command == L"rhodium") || (Command == L"Rh"))
	{
		Name = L"rhodium";
		Number = 45;
		Grp = 9;
		Period = 5;
		Formula = L"Rh";
		Charge = L"3+";
		Mass = 102.91;
		Special = L"None";

	}
	else if ((Command == L"46") || (Command == L"palladium") || (Command == L"Pd"))
	{
		Name = L"palladium";
		Number = 46;
		Grp = 10;
		Period = 5;
		Formula = L"Pd";
		Charge = L"2+ \\ 4+";
		Mass = 106.42;
		Special = L"None";

	}
	else if ((Command == L"47") || (Command == L"silver") || (Command == L"Ag"))
	{
		Name = L"silver";
		Number = 47;
		Grp = 11;
		Period = 5;
		Formula = L"Ag";
		Charge = L"1+";
		Mass = 107.87;
		Special = L"None";

	}
	else if ((Command == L"48") || (Command == L"cadmium") || (Command == L"Cd"))
	{
		Name = L"cadmium";
		Number = 48;
		Grp = 12;
		Period = 5;
		Formula = L"Cd";
		Charge = L"2+";
		Mass = 112.41;
		Special = L"None";

	}
	else if ((Command == L"49") || (Command == L"indium") || (Command == L"In"))
	{
		Name = L"indium";
		Number = 49;
		Grp = 13;
		Period = 5;
		Formula = L"In";
		Charge = L"3+";
		Mass = 114.82;
		Special = L"None";

	}
	else if ((Command == L"50") || (Command == L"tin") || (Command == L"Sn"))
	{
		Name = L"tin";
		Number = 50;
		Grp = 14;
		Period = 5;
		Formula = L"Sn";
		Charge = L"4+ \\ 2+";
		Mass = 118.69;
		Special = L"None";

	}
	else if ((Command == L"51") || (Command == L"antimony") || (Command == L"Sb"))
	{
		Name = L"antimony";
		Number = 51;
		Grp = 15;
		Period = 5;
		Formula = L"Sb";
		Charge = L"3+ \\ 5+";
		Mass = 121.75;
		Special = L"None";

	}
	else if ((Command == L"52") || (Command == L"tellurium") || (Command == L"Te"))
	{
		Name = L"tellurium";
		Number = 52;
		Grp = 16;
		Period = 5;
		Formula = L"Te";
		Charge = L"2-";
		Mass = 127.60;
		Special = L"None";

	}
	else if ((Command == L"53") || (Command == L"iodine") || (Command == L"I"))
	{
		Name = L"iodine";
		Number = 53;
		Grp = 17;
		Period = 5;
		Formula = L"I";
		Charge = L"1-";
		Mass = 126.90;
		Special = L"Diatomic";

	}
	else if ((Command == L"54") || (Command == L"xenon") || (Command == L"Xe"))
	{
		Name = L"xenon";
		Number = 54;
		Grp = 18;
		Period = 5;
		Formula = L"Xe";
		Charge = L"0";
		Mass = 131.29;
		Special = L"Noble Gas";

	}
	else if ((Command == L"55") || (Command == L"cosium") || (Command == L"Cs"))
	{
		Name = L"cosiumn";
		Number = 55;
		Grp = 1;
		Period = 6;
		Formula = L"Cs";
		Charge = L"1+";
		Mass = 132.91;
		Special = L"None";

	}
	else if ((Command == L"56") || (Command == L"barium") || (Command == L"Ba"))
	{
		Name = L"barium";
		Number = 56;
		Grp = 2;
		Period = 6;
		Formula = L"Ba";
		Charge = L"2+";
		Mass = 137.33;
		Special = L"None";

	}
	else if ((Command == L"57") || (Command == L"lanthanum") || (Command == L"La"))
	{
		Name = L"lanthanum";
		Number = 57;
		Grp = 3;
		Period = 6;
		Formula = L"La";
		Charge = L"3+";
		Mass = 138.91;
		Special = L"None";

	}
	else if ((Command == L"58") || (Command == L"cerium") || (Command == L"Ce"))
	{
		Name = L"cerium";
		Number = 58;
		Grp = 5;
		Period = 6;
		Formula = L"Ce";
		Charge = L"3+";
		Mass = 140.12;
		Special = L"None";

	}
	else if ((Command == L"59") || (Command == L"pruseodymium") || (Command == L"Pr"))
	{
		Name = L"pruseodymium";
		Number = 59;
		Grp = 6;
		Period = 6;
		Formula = L"Pr";
		Charge = L"3+";
		Mass = 140.91;
		Special = L"None";

	}
	else if ((Command == L"60") || (Command == L"neodymium") || (Command == L"Nd"))
	{
		Name = L"neodymium";
		Number = 60;
		Grp = 7;
		Period = 6;
		Formula = L"Nd";
		Charge = L"3+";
		Mass = 144.24;
		Special = L"None";

	}
	else if ((Command == L"61") || (Command == L"promethium") || (Command == L"Pm"))
	{
		Name = L"promethium";
		Number = 61;
		Grp = 8;
		Period = 6;
		Formula = L"Pm";
		Charge = L"3+";
		Mass = 145;
		Special = L"None";

	}
	else if ((Command == L"62") || (Command == L"samarium") || (Command == L"Sm"))
	{
		Name = L"samarium";
		Number = 62;
		Grp = 9;
		Period = 6;
		Formula = L"Sm";
		Charge = L"3+ \\ 2+";
		Mass = 150.40;
		Special = L"None";

	}
	else if ((Command == L"63") || (Command == L"europium") || (Command == L"Eu"))
	{
		Name = L"europium";
		Number = 63;
		Grp = 10;
		Period = 6;
		Formula = L"Eu";
		Charge = L"3+ \\ 2+";
		Mass = 151.97;
		Special = L"None";

	}
	else if ((Command == L"64") || (Command == L"gadolinium") || (Command == L"Gd"))
	{
		Name = L"gadolinium";
		Number = 64;
		Grp = 11;
		Period = 6;
		Formula = L"Gd";
		Charge = L"3+";
		Mass = 157.25;
		Special = L"None";

	}
	else if ((Command == L"65") || (Command == L"terbium") || (Command == L"Tb"))
	{
		Name = L"terbium";
		Number = 65;
		Grp = 12;
		Period = 6;
		Formula = L"Tb";
		Charge = L"3+";
		Mass = 158.93;
		Special = L"None";

	}
	else if ((Command == L"66") || (Command == L"dysprosium") || (Command == L"Dy"))
	{
		Name = L"dysprosium";
		Number = 66;
		Grp = 13;
		Period = 6;
		Formula = L"Dy";
		Charge = L"3+";
		Mass = 162.50;
		Special = L"None";

	}
	else if ((Command == L"67") || (Command == L"helmium") || (Command == L"Ho"))
	{
		Name = L"helmium";
		Number = 67;
		Grp = 14;
		Period = 6;
		Formula = L"Ho";
		Charge = L"3+";
		Mass = 164.93;
		Special = L"None";

	}
	else if ((Command == L"68") || (Command == L"erbium") || (Command == L"Er"))
	{
		Name = L"erbium";
		Number = 68;
		Grp = 15;
		Period = 6;
		Formula = L"Er";
		Charge = L"3+";
		Mass = 167.26;
		Special = L"None";

	}
	else if ((Command == L"69") || (Command == L"thulium") || (Command == L"Tm"))
	{
		Name = L"thulium";
		Number = 69;
		Grp = 16;
		Period = 6;
		Formula = L"Tm";
		Charge = L"3+";
		Mass = 168.94;
		Special = L"None";

	}
	else if ((Command == L"70") || (Command == L"ytlerhium") || (Command == L"Yb"))
	{
		Name = L"ytlerhium";
		Number = 70;
		Grp = 17;
		Period = 6;
		Formula = L"Yb";
		Charge = L"3+ \\ 2+";
		Mass = 173.04;
		Special = L"None";

	}
	else if ((Command == L"71") || (Command == L"lutelium") || (Command == L"Lu"))
	{
		Name = L"lutelium";
		Number = 71;
		Grp = 18;
		Period = 6;
		Formula = L"Lu";
		Charge = L"3+";
		Mass = 174.97;
		Special = L"None";

	}
	else if ((Command == L"72") || (Command == L"hefnium") || (Command == L"Hf"))
	{
		Name = L"hefnium";
		Number = 72;
		Grp = 4;
		Period = 6;
		Formula = L"Hf";
		Charge = L"4+";
		Mass = 178.49;
		Special = L"None";

	}
	else if ((Command == L"73") || (Command == L"tantalum") || (Command == L"Ta"))
	{
		Name = L"tantalum";
		Number = 73;
		Grp = 5;
		Period = 6;
		Formula = L"Ta";
		Charge = L"5+";
		Mass = 180.95;
		Special = L"None";

	}
	else if ((Command == L"74") || (Command == L"wolfrum") || (Command == L"tungsten") || (Command == L"W"))
	{
		Name = L"wolfrum (tungsten)";
		Number = 74;
		Grp = 6;
		Period = 6;
		Formula = L"W";
		Charge = L"6+";
		Mass = 183.85;
		Special = L"None";

	}
	else if ((Command == L"75") || (Command == L"rhenium") || (Command == L"Re"))
	{
		Name = L"rhenium";
		Number = 75;
		Grp = 7;
		Period = 6;
		Formula = L"Re";
		Charge = L"7+";
		Mass = 186.21;
		Special = L"None";

	}
	else if ((Command == L"76") || (Command == L"osmium") || (Command == L"Os"))
	{
		Name = L"osmium";
		Number = 76;
		Grp = 8;
		Period = 6;
		Formula = L"Os";
		Charge = L"4+";
		Mass = 190.2;
		Special = L"None";

	}
	else if ((Command == L"77") || (Command == L"iridium") || (Command == L"Ir"))
	{
		Name = L"iridium";
		Number = 77;
		Grp = 9;
		Period = 6;
		Formula = L"Ir";
		Charge = L"4+";
		Mass = 192.22;
		Special = L"None";

	}
	else if ((Command == L"78") || (Command == L"platinum") || (Command == L"Pt"))
	{
		Name = L"platinum";
		Number = 78;
		Grp = 10;
		Period = 6;
		Formula = L"Pt";
		Charge = L"4+ \\ 2+";
		Mass = 195.08;
		Special = L"None";

	}
	else if ((Command == L"79") || (Command == L"gold") || (Command == L"Au"))
	{
		Name = L"gold";
		Number = 79;
		Grp = 11;
		Period = 6;
		Formula = L"Au";
		Charge = L"3+ \\ 1+";
		Mass = 196.97;
		Special = L"None";

	}
	else if ((Command == L"80") || (Command == L"mercury") || (Command == L"Mg"))
	{
		Name = L"mercury";
		Number = 80;
		Grp = 12;
		Period = 6;
		Formula = L"Hg";
		Charge = L"2+ \\ 1+";
		Mass = 200.59;
		Special = L"None";

	}
	else if ((Command == L"81") || (Command == L"thallium") || (Command == L"Tl"))
	{
		Name = L"thallium";
		Number = 81;
		Grp = 13;
		Period = 6;
		Formula = L"Tl";
		Charge = L"1+ \\ 3+";
		Mass = 204.38;
		Special = L"None";

	}
	else if ((Command == L"82") || (Command == L"lead") || (Command == L"Pb"))
	{
		Name = L"lead";
		Number = 82;
		Grp = 14;
		Period = 6;
		Formula = L"Pb";
		Charge = L"2+ \\ 4+";
		Mass = 207.20;
		Special = L"None";

	}
	else if ((Command == L"83") || (Command == L"bismuth") || (Command == L"Bi"))
	{
		Name = L"bismuth";
		Number = 83;
		Grp = 15;
		Period = 6;
		Formula = L"Bi";
		Charge = L"3+ \\ 5+";
		Mass = 208.98;
		Special = L"None";

	}
	else if ((Command == L"84") || (Command == L"polonium") || (Command == L"Po"))
	{
		Name = L"polonium";
		Number = 84;
		Grp = 16;
		Period = 6;
		Formula = L"Po";
		Charge = L"2+ \\ 4+";
		Mass = 209;
		Special = L"None";

	}
	else if ((Command == L"85") || (Command == L"asiatine") || (Command == L"At"))
	{
		Name = L"asiatine";
		Number = 85;
		Grp = 17;
		Period = 6;
		Formula = L"At";
		Charge = L"1-";
		Mass = 210;
		Special = L"None";

	}
	else if ((Command == L"86") || (Command == L"radon") || (Command == L"Rn"))
	{
		Name = L"radon";
		Number = 86;
		Grp = 18;
		Period = 6;
		Formula = L"Rn";
		Charge = L"0";
		Mass = 222;
		Special = L"Noble Gas";

	}
	else if ((Command == L"87") || (Command == L"fruncium") || (Command == L"Fr"))
	{
		Name = L"fruncium";
		Number = 87;
		Grp = 1;
		Period = 7;
		Formula = L"Fr";
		Charge = L"1+";
		Mass = 223;
		Special = L"None";

	}
	else if ((Command == L"88") || (Command == L"radium") || (Command == L"Ra"))
	{
		Name = L"radium";
		Number = 88;
		Grp = 2;
		Period = 7;
		Formula = L"Ra";
		Charge = L"2+";
		Mass = 226.03;
		Special = L"None";

	}
	else if ((Command == L"89") || (Command == L"actinium") || (Command == L"Ac"))
	{
		Name = L"actinium";
		Number = 89;
		Grp = 3;
		Period = 7;
		Formula = L"Ac";
		Charge = L"3+";
		Mass = 227.03;
		Special = L"None";

	}
	else if ((Command == L"90") || (Command == L"thorlum") || (Command == L"Th"))
	{
		Name = L"thorlum";
		Number = 90;
		Grp = 5;
		Period = 7;
		Formula = L"Th";
		Charge = L"4+";
		Mass = 232.04;
		Special = L"None";

	}
	else if ((Command == L"91") || (Command == L"protactinium") || (Command == L"Pa"))
	{
		Name = L"protactinium";
		Number = 91;
		Grp = 6;
		Period = 7;
		Formula = L"Pa";
		Charge = L"5+ \\ 4+";
		Mass = 231.04;
		Special = L"None";

	}
	else if ((Command == L"92") || (Command == L"uranium") || (Command == L"U"))
	{
		Name = L"uranium";
		Number = 92;
		Grp = 7;
		Period = 7;
		Formula = L"U";
		Charge = L"6+ \\ 4+";
		Mass = 238.03;
		Special = L"None";

	}
	else if ((Command == L"93") || (Command == L"neplunium") || (Command == L"Np"))
	{
		Name = L"neplunium";
		Number = 93;
		Grp = 8;
		Period = 7;
		Formula = L"Np";
		Charge = L"5+";
		Mass = 237.05;
		Special = L"None";

	}
	else if ((Command == L"94") || (Command == L"plutonium") || (Command == L"Pu"))
	{
		Name = L"plutonium";
		Number = 94;
		Grp = 9;
		Period = 7;
		Formula = L"Pu";
		Charge = L"4+ \\ 6+";
		Mass = 244;
		Special = L"None";

	}
	else if ((Command == L"95") || (Command == L"americium") || (Command == L"Am"))
	{
		Name = L"americium";
		Number = 95;
		Grp = 10;
		Period = 7;
		Formula = L"Am";
		Charge = L"3+ \\ 4+";
		Mass = 244;
		Special = L"None";

	}
	else if ((Command == L"96") || (Command == L"curium") || (Command == L"Cm"))
	{
		Name = L"curium";
		Number = 96;
		Grp = 11;
		Period = 7;
		Formula = L"Cm";
		Charge = L"3+";
		Mass = 247;
		Special = L"None";

	}
	else if ((Command == L"97") || (Command == L"borkelium") || (Command == L"Bk"))
	{
		Name = L"borkelium";
		Number = 97;
		Grp = 12;
		Period = 7;
		Formula = L"Bk";
		Charge = L"3+ \\ 4+";
		Mass = 247;
		Special = L"None";

	}
	else if ((Command == L"98") || (Command == L"californium") || (Command == L"Cf"))
	{
		Name = L"californium";
		Number = 98;
		Grp = 13;
		Period = 7;
		Formula = L"Cf";
		Charge = L"3+";
		Mass = 251;
		Special = L"None";

	}
	else if ((Command == L"99") || (Command == L"einsteinium") || (Command == L"Es"))
	{
		Name = L"einsteinium";
		Number = 99;
		Grp = 14;
		Period = 7;
		Formula = L"Es";
		Charge = L"3+";
		Mass = 252;
		Special = L"None";

	}
	else if ((Command == L"100") || (Command == L"formium") || (Command == L"Fm"))
	{
		Name = L"formium";
		Number = 100;
		Grp = 15;
		Period = 7;
		Formula = L"Fm";
		Charge = L"3+";
		Mass = 257;
		Special = L"None";

	}
	else if ((Command == L"101") || (Command == L"mendelevium") || (Command == L"Md"))
	{
		Name = L"mendelevium";
		Number = 101;
		Grp = 16;
		Period = 7;
		Formula = L"Md";
		Charge = L"2+ \\ 3+";
		Mass = 258;
		Special = L"None";

	}
	else if ((Command == L"102") || (Command == L"nebelium") || (Command == L"No"))
	{
		Name = L"nebelium";
		Number = 102;
		Grp = 17;
		Period = 7;
		Formula = L"No";
		Charge = L"2+ \\ 3+";
		Mass = 259;
		Special = L"None";

	}
	else if ((Command == L"103") || (Command == L"lawrencium") || (Command == L"Lr"))
	{
		Name = L"lawrencium";
		Number = 103;
		Grp = 18;
		Period = 7;
		Formula = L"Lr";
		Charge = L"3+";
		Mass = 260;
		Special = L"None";

	}
	else if ((Command == L"104") || (Command == L"rutherfordium") || (Command == L"Rf"))
	{
		Name = L"rutherfordium";
		Number = 104;
		Grp = 4;
		Period = 7;
		Formula = L"Rf";
		Charge = L"4+";
		Mass = 267;
		Special = L"None";

	}
	else if ((Command == L"105") || (Command == L"dubnium") || (Command == L"Db"))
	{
		Name = L"dubnium";
		Number = 105;
		Grp = 5;
		Period = 7;
		Formula = L"Db";
		Charge = L"5+";
		Mass = 268;
		Special = L"None";

	}
	else if ((Command == L"106") || (Command == L"seaborgium") || (Command == L"Sg"))
	{
		Name = L"seaborgium";
		Number = 106;
		Grp = 6;
		Period = 7;
		Formula = L"Sg";
		Charge = L"6+";
		Mass = 269;
		Special = L"None";

	}
	else if ((Command == L"107") || (Command == L"bohrium") || (Command == L"Bh"))
	{
		Name = L"bohrium";
		Number = 107;
		Grp = 7;
		Period = 7;
		Formula = L"Bh";
		Charge = L"7+";
		Mass = 278;
		Special = L"None";

	}
	else if ((Command == L"108") || (Command == L"hassium") || (Command == L"Hs"))
	{
		Name = L"hassium";
		Number = 108;
		Grp = 8;
		Period = 7;
		Formula = L"Hs";
		Charge = L"?";
		Mass = 277;
		Special = L"None";

	}
	else if ((Command == L"109") || (Command == L"meitnerium") || (Command == L"Mt"))
	{
		Name = L"meitnerium";
		Number = 109;
		Grp = 9;
		Period = 7;
		Formula = L"Mt";
		Charge = L"?";
		Mass = 278;
		Special = L"None";

	}
	else if ((Command == L"110") || (Command == L"darmstatium") || (Command == L"Ds"))
	{
		Name = L"darmstatium";
		Number = 110;
		Grp = 10;
		Period = 7;
		Formula = L"Ds";
		Charge = L"?";
		Mass = 281;
		Special = L"None";

	}
	else if ((Command == L"111") || (Command == L"roentgenium") || (Command == L"Rg"))
	{
		Name = L"roentgenium";
		Number = 111;
		Grp = 11;
		Period = 7;
		Formula = L"Rg";
		Charge = L"?";
		Mass = 282;
		Special = L"None";

	}

	else if ((Command == L"112") || (Command == L"copernicium") || (Command == L"Cn"))
	{
		Name = L"copernicium";
		Number = 112;
		Grp = 12;
		Period = 7;
		Formula = L"Cn";
		Charge = L"?";
		Mass = 285;
		Special = L"None";

	}
	else if ((Command == L"113") || (Command == L"nihonium") || (Command == L"Nh"))
	{
		Name = L"nihonium";
		Number = 113;
		Grp = 13;
		Period = 7;
		Formula = L"Nh";
		Charge = L"?";
		Mass = 286;
		Special = L"None";

	}
	else if ((Command == L"114") || (Command == L"flerovium") || (Command == L"Fl"))
	{
		Name = L"flerovium";
		Number = 114;
		Grp = 14;
		Period = 7;
		Formula = L"Fl";
		Charge = L"?";
		Mass = 289;
		Special = L"None";

	}
	else if ((Command == L"115") || (Command == L"moscovium") || (Command == L"Mc"))
	{
		Name = L"moscovium";
		Number = 115;
		Grp = 15;
		Period = 7;
		Formula = L"Mc";
		Charge = L"?";
		Mass = 290;
		Special = L"None";

	}
	else if ((Command == L"116") || (Command == L"livermorium") || (Command == L"Lv"))
	{
		Name = L"livermorium";
		Number = 116;
		Grp = 16;
		Period = 7;
		Formula = L"Lv";
		Charge = L"?";
		Mass = 293;
		Special = L"None";

	}
	else if ((Command == L"117") || (Command == L"tennessine") || (Command == L"Ts"))
	{
		Name = L"tennesine";
		Number = 117;
		Grp = 17;
		Period = 7;
		Formula = L"Ts";
		Charge = L"?";
		Mass = 294;
		Special = L"None";

	}
	else if ((Command == L"118") || (Command == L"oganesson") || (Command == L"Og"))
	{
		Name = L"oganesson";
		Number = 118;
		Grp = 18;
		Period = 7;
		Formula = L"Og";
		Charge = L"?";
		Mass = 294;
		Special = L"None";

	}


	//====================================
	else if (stringToDouble(Command)>103) {
		std::wstring nl = L"";
		return nl;
	}
	if (n <= 0 || n>8) {
		return std::wstring(
			L"Element name: " + Name + L",\n" +
			L"Element atomic number: " + std::to_wstring(Number) + L",\n" +
			L"Element group: " + std::to_wstring(Grp) + L",\n" +
			L"Element period: " + std::to_wstring(Period) + L",\n" +
			L"Element formula: " + Formula + L",\n" +
			L"Element common ionic charge: " + Charge + L",\n" +
			L"Element atomic mass: " + to_stringPrecision(Mass) + L",\n" +
			L"Element special info: " + Special);
	}
	if (n == 2) { return std::wstring(L"Element atomic number: " + std::to_wstring(Number)); }
	if (n == 3) { return std::wstring(L"Element group: " + std::to_wstring(Grp)); }
	if (n == 4) { return std::wstring(L"Element period: " + std::to_wstring(Period)); }
	if (n == 5) { return std::wstring(L"Element formula: " + Formula); }
	if (n == 6) { return std::wstring(L"Element common ionic charge: " + Charge); }
	if (n == 7) { return std::wstring(L"Element atomic mass: " + to_stringPrecision(Mass)); }
	if (n == 8) { return std::wstring(L"Element special info: " + Special); }
	std::wstring nl = L"";
	return nl;
}


std::wstring PeriodicTable() {//prints whole periodic table as a string
	std::wstring str = L"";
	for (int i = 1; i < 118; ++i) {
		str.append(PeriodicTable(std::to_wstring(i), 0) + L"\n\n");
	}
	str.append(PeriodicTable(std::to_wstring(118), 0));
	return str;
}