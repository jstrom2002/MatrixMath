#pragma once
#include "stdafx.h"
#include <direct.h>

std::wstring namesOfSporadicGroups = L"M11,M12,M22,M23,M24,HS,J2,Co1,Co2,Co3,McL,Suz,He,HN,Th,Fi22,Fi23,Fi24,B,M,J1,ON,J3,Ru,J4,Ly,T";

//Note: when possible, the representations over Z are used here.  If that's not available, then
//a representation over the boolean field GF(2) is used.  

//Sporadic Simple Groups
//======================
void M11() { ParseMatrixLoadFile(L"Atlas/M11.txt"); }
void M12() { ParseMatrixLoadFile(L"Atlas/M12.txt"); }
void M22() { ParseMatrixLoadFile(L"Atlas/M22.txt"); }
void M23() { ParseMatrixLoadFile(L"Atlas/M23.txt"); }
void M24() { ParseMatrixLoadFile(L"Atlas/M24.txt"); }
void HS() { ParseMatrixLoadFile(L"Atlas/HS.txt"); }
void J2() { ParseMatrixLoadFile(L"Atlas/J2.txt"); }
void Co1() { ParseMatrixLoadFile(L"Atlas/Co1.txt"); }
void Co2() { ParseMatrixLoadFile(L"Atlas/Co2.txt"); }
void Co3() { ParseMatrixLoadFile(L"Atlas/Co3.txt"); }
void McL() { ParseMatrixLoadFile(L"Atlas/McL.txt"); }
void Suz() { ParseMatrixLoadFile(L"Atlas/Suz.txt"); }
void He() { ParseMatrixLoadFile(L"Atlas/He.txt"); }
void HN() { ParseMatrixLoadFile(L"Atlas/HN.txt"); }
void Th() { ParseMatrixLoadFile(L"Atlas/Th.txt"); }
void Fi22() { ParseMatrixLoadFile(L"Atlas/Fi22.txt"); }
void Fi23() { ParseMatrixLoadFile(L"Atlas/Fi23.txt"); }
void Fi24() { ParseMatrixLoadFile(L"Atlas/Fi24.txt"); }
void B() { ParseMatrixLoadFile(L"Atlas/B.txt"); }
void M() { ParseMatrixLoadFile(L"Atlas/M.txt"); }
void J1() { ParseMatrixLoadFile(L"Atlas/J1.txt"); }
void ON() { ParseMatrixLoadFile(L"Atlas/ON.txt"); }
void J3() { ParseMatrixLoadFile(L"Atlas/J3.txt"); }
void Ru() { ParseMatrixLoadFile(L"Atlas/Ru.txt"); }
void J4() { ParseMatrixLoadFile(L"Atlas/J4.txt"); }
void Ly() { ParseMatrixLoadFile(L"Atlas/Ly.txt"); }
void T() { ParseMatrixLoadFile(L"Atlas/T.txt"); }

//Alternating Groups
//==================
void A5() { ParseMatrixLoadFile(L"Atlas/A5.txt"); }
void A6() { ParseMatrixLoadFile(L"Atlas/A6.txt"); }
void A7() { ParseMatrixLoadFile(L"Atlas/A7.txt"); }
void A8() { ParseMatrixLoadFile(L"Atlas/A8.txt"); }
void A9() { ParseMatrixLoadFile(L"Atlas/A9.txt"); }
void A10() { ParseMatrixLoadFile(L"Atlas/A10.txt"); }
void A11() { ParseMatrixLoadFile(L"Atlas/A11.txt"); }
void A12() { ParseMatrixLoadFile(L"Atlas/A12.txt"); }
void A13() { ParseMatrixLoadFile(L"Atlas/A13.txt"); }
void A14() { ParseMatrixLoadFile(L"Atlas/A14.txt"); }
void A15() { ParseMatrixLoadFile(L"Atlas/A15.txt"); }
void A16() { ParseMatrixLoadFile(L"Atlas/A16.txt"); }
void A17() { ParseMatrixLoadFile(L"Atlas/A17.txt"); }
void A18() { ParseMatrixLoadFile(L"Atlas/A18.txt"); }
void A19() { ParseMatrixLoadFile(L"Atlas/A19.txt"); }
void A20() { ParseMatrixLoadFile(L"Atlas/A20.txt"); }
void A21() { ParseMatrixLoadFile(L"Atlas/A21.txt"); }
void A22() { ParseMatrixLoadFile(L"Atlas/A22.txt"); }
void A23() { ParseMatrixLoadFile(L"Atlas/A23.txt"); }

//Linear Groups
//=============
void L32() { ParseMatrixLoadFile(L"Atlas/L3(2).txt"); }


//Classical Groups
//================

//Exceptional Groups of Lie Type
//==============================

//Miscellaneous Groups
//=====================
void M20() { ParseMatrixLoadFile(L"Atlas/M20.txt"); }
