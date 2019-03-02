//MATRIX MATH v1.2 by JH Strom
//copyright 2017, Canada

#include "stdafx.h"
#include "vector.h"
#include "MainParser.h"
#include "Output.h"
#include "Security.h"
#include "Resource.h"
#include "String.h"
#include <direct.h>
#include <io.h>
#include <fcntl.h>

//global variables
bool changeScreen;
bool first;
bool fullscreen;
bool needOutputWindow;
bool running;
bool showTime;
std::wstring filename;
std::wstring hostname;
std::vector<std::wstring> IPs;
std::wstring productKey;
std::wstring path;
std::wstring serialnumber;

//WinAPI specific variables/objects
LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam);
bool CALLBACK SetFont(HWND child, LPARAM font);
LPCWSTR inputLPCWSTR; //converted LPCSTR from input box for use in message boxes
LPCWSTR outputLPCWSTR;
std::vector<LPCWSTR> outputLPCWSTRarr;
LPCWSTR errorMessage = L"Error!";
TEXTMETRIC tm;

int fullscreenWidth;
int fullscreenHeight;
int colourBits;
int refreshRate;

HDC hdc;					 //parent device context
HINSTANCE hInstance;		 //current instance
MSG messages;			     //parent message container
HWND hwnd;				     //main window
HWND outputPanel;		     //output text panel
HWND hEdit;					 //input bar for user-entered strings
HWND errorBox;				 //error dialog box

const TCHAR szClassName[] = TEXT("MatrixMath"); //main window class name = szClassName
std::wstring szClassName2 = L"MatrixMath";
#define ID_Save 1
#define ID_Load 2
#define ID_Exit 3
#define ID_Undo 4
#define ID_Redo 5
#define ID_HELP 6
#define IDC_MAIN_BUTTON  104         
#define IDC_MAIN_EDIT    105    // Edit box identifier -- this panel displays output text
#define IDC_COMMAND_LINE 106		
#define IDC_LOGO_PANE    107	
#define ICD_ERROR_BOX	 108	// Error box identifier -- shows a pop box message when an error is thrown

bool enterFullscreen(HWND hwnd, int fullscreenWidth, int fullscreenHeight, int colourBits, int refreshRate);
bool exitFullscreen(HWND hwnd, int windowX, int windowY, int windowedWidth, int windowedHeight, int windowedPaddingX, int windowedPaddingY);
void SetTransparency(HWND hwnd, std::uint8_t Transperancy);

/*===============================================
//					MAIN  METHOD				//
================================================*/
int WINAPI WinMain(HINSTANCE hThisInstance, HINSTANCE hPrevInstance, LPSTR lpszArgument, int nFunsterStil) {

	//on load, run security checks, make fullscreen, initialize variables, etc.
	bool SN = checkProductKey();
	if (SN == false) { 
		LPCTSTR Caption = L"ERROR!!!";
		MessageBox(NULL,
			L"Pirated copy of MatrixMath detected.\n"
			L"Shutting down and notifying the developer.",
			Caption,
			MB_ICONERROR);
		return 0;
		exit(0); 			
	}
	if (SN == true) {
		filename = L"MatrixMath.txt";
		first = true;
		fullscreen = true;
		running = true;
		showTime = false;
		changeScreen = false;
		precision = 4;
		loadData();
		needOutputWindow = true;
	}


	hdc = GetDC(NULL);			//get parent device context
    hInstance = hThisInstance;	//get instance

	//Set up main window class
	//=========================
	WNDCLASSEX wincl;			//main window class
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = (szClassName2.c_str());
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof(WNDCLASSEX);
    wincl.hIcon = LoadIcon(hThisInstance, MAKEINTRESOURCE(IDI_ICON1));
    wincl.hIconSm = LoadIcon(hThisInstance, MAKEINTRESOURCE(IDI_ICON1));
    wincl.hCursor = LoadCursor(NULL, IDC_ARROW);				//use standard mouse pointer
    wincl.lpszMenuName = NULL;									//add menu later
    wincl.cbClsExtra = 0;										/* No extra bytes after the window class */
    wincl.cbWndExtra = 0;										/* structure or the window instance */                                             
    wincl.hbrBackground = CreatePatternBrush((HBITMAP)LoadImage(0, L"bg2.bmp",IMAGE_BITMAP, 0, 0,LR_CREATEDIBSECTION | LR_LOADFROMFILE));
 
    /* Register the window class, and if it fails quit the program */
	if (!RegisterClassEx(&wincl)) { return 0; }

	//get window size, refresh rate, etc
	fullscreenWidth = GetDeviceCaps(hdc, HORZRES);
	fullscreenHeight = GetDeviceCaps(hdc, VERTRES);
	colourBits = GetDeviceCaps(hdc, BITSPIXEL);
	refreshRate = GetDeviceCaps(hdc, VREFRESH);

	//Create main window	
	hwnd = CreateWindowEx(
		0,                   /* Extended possibilites for variation */
		szClassName,         /* Classname */
		L"MatrixMath v1.1",       /* Title Text */
		WS_OVERLAPPEDWINDOW, /* default window */
		CW_USEDEFAULT,       /* Windows decides the position */
		CW_USEDEFAULT,       /* where the window ends up on the screen */
		GetSystemMetrics(SM_CXSCREEN) + 10,   /* The programs width */
		GetSystemMetrics(SM_CYSCREEN),  /* and height in pixels */
		0,        /* The window is a child-window to desktop */
		0,                /* No menu yet, this is added later*/
		hThisInstance,       /* Program Instance handler */
		0                 /* No Window Creation data */
	);
	ShowWindow(hwnd, SW_MAXIMIZE); //Make the window visible on the screen 

	// SET UP GUI
	//===========
	//set up menu
	HMENU hMenubar = CreateMenu();
	HMENU hFileMenu = CreateMenu();
	HMENU hEditMenu = CreateMenu();
	HMENU hHelpMenu = CreateMenu();

	AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hFileMenu, TEXT("File"));
	AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hEditMenu, TEXT("Edit"));
	AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hHelpMenu, TEXT("Help"));

	AppendMenu(hFileMenu, MF_STRING, ID_Load, TEXT("Load"));
	AppendMenu(hFileMenu, MF_STRING, ID_Save, TEXT("Save"));
	AppendMenu(hFileMenu, MF_STRING, ID_Exit, TEXT("Exit"));

	AppendMenu(hEditMenu, MF_STRING, ID_Undo, TEXT("Undo"));
	AppendMenu(hEditMenu, MF_STRING, ID_Redo, TEXT("Redo"));

	AppendMenu(hHelpMenu, MF_STRING, ID_HELP, TEXT("Help file"));
	SetMenu(hwnd, hMenubar);

	//=====================
	
	//Create "enter command: " text box
	HWND enterCommand = CreateWindowEx(
		0,
		L"STATIC", 
		L"Enter command:",
		WS_VISIBLE | WS_CHILD,
		10, 13, 110, 16, 
		hwnd, 
		0, 
		hThisInstance, 
		0);      


	//====================
	// Create the command text field
	hEdit = CreateWindowEx(
		WS_EX_CLIENTEDGE,
		TEXT("EDIT"),
		TEXT(""),
		WS_CHILD | WS_VISIBLE,
		150, 9, fullscreenWidth - 350, 24,     //x,y,w,h
		hwnd,
		(HMENU)IDC_COMMAND_LINE, 
		GetModuleHandle(NULL), 
		0);
	HGDIOBJ hfDefault = GetStockObject(DS_FIXEDSYS);
	SendMessage(hEdit, WM_SETFONT, (WPARAM)hfDefault, MAKELPARAM(FALSE, 0));
	SendMessage(hEdit, EM_SETLIMITTEXT, 2147483646, 0);
	//=========================

	// Create the static control for the logo bar at the top
	HWND logo = CreateWindowEx(
		0,
		TEXT("STATIC"),
		TEXT(""),
		WS_CHILD | WS_VISIBLE | SS_BITMAP | SS_CENTERIMAGE,
		(fullscreenWidth/2)-300, 35, 525, 100,     //x,y,w,h
		hwnd,
		(HMENU)IDC_LOGO_PANE,
		GetModuleHandle(NULL),
		0);
	//load image for button
	HBITMAP logoBMP = (HBITMAP)LoadImageW(NULL, L"logobar.bmp", IMAGE_BITMAP, 0, 0, LR_DEFAULTCOLOR | LR_DEFAULTSIZE | LR_LOADFROMFILE);
	SendMessageW(logo, STM_SETIMAGE, IMAGE_BITMAP, (LPARAM)logoBMP);

	//=========================
	// Create the push button
	HWND hWndButton = CreateWindowEx(
		0,
		TEXT("BUTTON"),
		TEXT(""),
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_PUSHBUTTON | BS_BITMAP,
		fullscreenWidth - 145, 7, 95, 29,	//x,y,w,h
		hwnd,
		(HMENU)IDC_MAIN_BUTTON,
		GetModuleHandle(NULL),
		NULL);
	SendMessage(hWndButton,
		WM_SETFONT,
		(WPARAM)hfDefault,
		MAKELPARAM(FALSE, 0));
	
	//load image for button
	HBITMAP hbutton = (HBITMAP)LoadImageW(NULL, L"button2Crop.bmp", IMAGE_BITMAP, 0, 0, LR_DEFAULTCOLOR | LR_DEFAULTSIZE | LR_LOADFROMFILE);
	SendMessageW(hWndButton, BM_SETIMAGE, IMAGE_BITMAP, (LPARAM)hbutton);

	//==========================
	//create an output panel for displaying text
	outputPanel = CreateWindowEx(
		WS_EX_CLIENTEDGE, 
		L"EDIT", 
		outputLPCWSTR, 
		WS_CHILD | WS_VISIBLE | ES_LEFT | ES_MULTILINE | WS_VSCROLL | WS_HSCROLL | ES_AUTOHSCROLL | ES_AUTOVSCROLL,
		12, 136, fullscreenWidth - 30, fullscreenHeight - 160,    //offset height/width so that the program doesn't overlap the taskbar
		hwnd, 
		(HMENU)IDC_MAIN_EDIT,
		hThisInstance, 
		0);
	ShowWindow(outputPanel, nFunsterStil);

	//set font for output window
	HFONT hFont = CreateFont(11, 5, 0, 0, FW_DONTCARE, FALSE, FALSE, FALSE, ANSI_CHARSET,
		OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY,
		DEFAULT_PITCH | FF_DONTCARE, TEXT("Fixedsys"));             //fixedsys is a monospaced font (so it will display matrices properly) and that most older computers will recognize
	SendMessage(outputPanel, WM_SETFONT, (WPARAM)hFont, TRUE);
	SendMessage(outputPanel, EM_SETLIMITTEXT, 2147483646, 0);
	
	//===============================

	//get path of where .exe is running
	char cCurrentPath[FILENAME_MAX];
	if (!_getcwd(cCurrentPath, sizeof(cCurrentPath))) {	return errno;}
	path = s2ws(cCurrentPath);

	//===========//
	// MAIN LOOP //
	//===========//
	while (GetMessage(&messages, NULL, 0, 0) && running==true) {   //Run the message loop. It will run until GetMessage() returns 0
		TranslateMessage(&messages);            /* Translate virtual-key messages into character messages */
		DispatchMessage(&messages);     /* Send message to WindowProcedure */
		if (GetAsyncKeyState(VK_RETURN) < 0) { SendMessage(hWndButton, BM_CLICK, 0, 0); Sleep(150); }    //if user hits "enter," process input
		if (GetAsyncKeyState(VK_LCONTROL) < 0 && GetAsyncKeyState('A') < 0) { SendMessage(hEdit, EM_SETSEL, 0, 99999999); Sleep(200); }    //if user hits "enter," process input


		//toggle fullscreen if necessary
		if (changeScreen == true) {
			if (fullscreen == false) {exitFullscreen(hwnd, 0, 0, fullscreenWidth, fullscreenHeight, 0, 0);}
			if (fullscreen == true) {enterFullscreen(hwnd, fullscreenWidth, fullscreenHeight, colourBits, refreshRate);}
			changeScreen = false;
		}

		//UPDATE TEXT FIELD LOOP
		//======================
		if (needOutputWindow == true) {
			output.clear();
			clock_t timer1, timer2;
			std::vector <std::wstring> ins;
			
			if (first != true) {
				if (input.find(L";") != std::wstring::npos && input.find(L"[") == std::wstring::npos) {//if input is multi-command demarcated with ';' char...
					while (input.find(L";") != std::wstring::npos) {
						std::wstring tempStr = input.substr(0, input.find(L";"));
						tempStr = removeSpaces(tempStr);
						ins.push_back(tempStr);
						input.erase(0, input.find(L";") + 1);
					}
					if (input.length() > 0) { ins.push_back(input); }
				}
				
				if (showTime == true) { timer1 = clock(); }//Parse input and handle clocking of calculation

				if (ins.size() > 0) {
					for (int i = 0; i < ins.size(); ++i) {//parse each line
						MainParser(ins[i]);
					}
				}else { MainParser(input); }//parse input string and do calculations
				
				if (showTime == true && input != L"showTime") {
					timer2 = clock();
					float diff((float)timer2 - (float)timer1);
					float sc = (diff / CLOCKS_PER_SEC);
					output.insert(0, L" seconds\n\n");
					output.insert(0, to_stringPrecision(sc));
					output.insert(0, L"Calculation Time: ");
				}
				if (input == L"showTime") {
					if (showTime == true) {
						showTime = false;
						output.append(L"show calculation time: OFF.\n\n");
					}
					if (showTime == false) {
						showTime = true;
						output.append(L"show calculation time: ON.\n\n");
					}
				}
			}
	
			output.append(printOutput());//display saved variables and objects
			//output = formatStringAsWString(output);
			outputLPCWSTR = output.c_str();//convert output string to an LPCSTR for display on a window
			SetWindowText(outputPanel, outputLPCWSTR);
			SetFocus(hEdit);
			if (first == true) { first = false; }
			needOutputWindow = false;	//stop this loop from running after it is finished
		}//==============================
	}
	return messages.wParam;
}
/*===============================================
//				END  MAIN  METHOD				//
================================================*/
	
LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam){
	switch (message) {
	case WM_CREATE: {break;}
	case WM_COMMAND: {
		if(LOWORD(wParam) == IDC_MAIN_BUTTON)       //catch if user presses button
		{
			output.clear();
			outputLPCWSTR = L"";
			wchar_t buffer[256];
			SendMessage(hEdit,
				WM_GETTEXT,
				sizeof(buffer) / sizeof(buffer[0]),
				reinterpret_cast<LPARAM>(buffer));
			std::wstring ws1 = buffer;
			input = ws1;
			inputLPCWSTR = input.c_str();
			needOutputWindow = true;
		}//===========================

		if (LOWORD(wParam) == ID_Load) {//if user clicks the 'load file' menu option
			//set load menu condition
			OPENFILENAME ofn;
			char szFileName[MAX_PATH] = "";
			std::wstring tmp = s2ws(szFileName);
			LPCWSTR tmp2 = tmp.c_str();

			ZeroMemory(&ofn, sizeof(ofn));

			//set file picker attributes
			ofn.lStructSize = sizeof(ofn);
			ofn.hwndOwner = hwnd;
			ofn.lpstrFilter = TEXT("Text Files (*.txt)\0*.txt\0");
			ofn.lpstrFile = (LPWSTR)szFileName;
			ofn.nMaxFile = MAX_PATH;
			ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
			ofn.lpstrDefExt = TEXT("txt");

			if (GetOpenFileName(&ofn))
			{
				ClearContainers();
				std::wstring tps = filename;
				filename = s2ws(szFileName);
				loadData();
				filename = tps;
				input = L"load";
				output = formatStringAsWString(output);
				std::wstring tmp = output;    //convert output string to an LPCSTR for display on a window
				outputLPCWSTR = tmp.c_str();
				SetWindowText(outputPanel, outputLPCWSTR);
				SetFocus(hEdit);
				needOutputWindow = true;
			}
			break;
		}//===========================

		if (LOWORD(wParam) == ID_Save) {//if user clicks the 'save file' menu option
			//set save condition
			OPENFILENAME ofn;
			wchar_t szFileName[MAX_PATH] = L"";
			ZeroMemory(&ofn, sizeof(ofn));

			//options for file picker
			ofn.lStructSize = sizeof(ofn);
			ofn.hwndOwner = hwnd;
			ofn.lpstrFilter = TEXT("Text Files (*.txt)\0*.txt\0");
			ofn.lpstrFile = szFileName;
			ofn.nMaxFile = MAX_PATH;
			ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
			ofn.lpstrDefExt = TEXT("txt");

			if (GetSaveFileName(&ofn))
			{
				// Now save the new file name in the 'filename' variable and run saveData()
				int positn = 0;
				for (int i = 0; i < sizeof(szFileName); ++i) { if (szFileName[i] == '\\') { positn = i; } }//find last '\' char and remove everything before that, leaving only the new filename in the string (not the path)
				std::wstring temp = szFileName;
				if (positn < temp.length()) { temp = temp.substr(positn + 1); }
				filename = temp; //set the default file name to the new file name
				input = L"save";  //let the MainParser method do the work of saving the file
				needOutputWindow = true;
			}
			break;
		}//===========================

		if (LOWORD(wParam) == ID_Exit) {
			exit(0);
		}
		if (LOWORD(wParam) == ID_Undo) {
			//set undo option
		}
		if (LOWORD(wParam) == ID_Redo) {
			//set redo option
		}
		if (LOWORD(wParam) == ID_HELP) {
			ShellExecute(NULL, TEXT("Open"), TEXT("https://matrixmathofficial.wordpress.com/help"), NULL, NULL, SW_SHOWNORMAL);
		}
		break;
	}//===========================

	case WM_DESTROY:
	{
		PostQuitMessage(0);       /* send a WM_QUIT to the message queue */
		break;
	}
	default:                      /* for messages that we don't deal with */
		return DefWindowProc(hwnd, message, wParam, lParam);
	}//===========================
	return 0;
}//END WindowProcedure================================================


bool CALLBACK SetFont(HWND child, LPARAM font) {
	SendMessage(child, WM_SETFONT, font, true);
	return true;
}

bool enterFullscreen(HWND hwnd, int fullscreenWidth, int fullscreenHeight, int colourBits, int refreshRate) {
	DEVMODE fullscreenSettings;
	bool isChangeSuccessful;
	RECT windowBoundary;

	EnumDisplaySettings(NULL, 0, &fullscreenSettings);
	fullscreenSettings.dmPelsWidth = fullscreenWidth;
	fullscreenSettings.dmPelsHeight = fullscreenHeight;
	fullscreenSettings.dmBitsPerPel = colourBits;
	fullscreenSettings.dmDisplayFrequency = refreshRate;
	fullscreenSettings.dmFields = DM_PELSWIDTH |
		DM_PELSHEIGHT |
		DM_BITSPERPEL |
		DM_DISPLAYFREQUENCY;

	SetWindowLongPtr(hwnd, GWL_EXSTYLE, WS_EX_APPWINDOW | WS_EX_TOPMOST);
	SetWindowLongPtr(hwnd, GWL_STYLE, WS_POPUP | WS_VISIBLE);
	SetWindowPos(hwnd, HWND_NOTOPMOST, 0, 0, fullscreenWidth, fullscreenHeight, SWP_SHOWWINDOW);
	isChangeSuccessful = ChangeDisplaySettings(&fullscreenSettings, CDS_FULLSCREEN) == DISP_CHANGE_SUCCESSFUL;
	ShowWindow(hwnd, SW_MAXIMIZE);

	return isChangeSuccessful;
}

bool exitFullscreen(HWND hwnd, int windowX, int windowY, int windowedWidth, int windowedHeight, int windowedPaddingX, int windowedPaddingY) {
	bool isChangeSuccessful;

	SetWindowLongPtr(hwnd, GWL_EXSTYLE, WS_EX_LEFT);
	SetWindowLongPtr(hwnd, GWL_STYLE, WS_OVERLAPPEDWINDOW | WS_VISIBLE);
	isChangeSuccessful = ChangeDisplaySettings(NULL, CDS_RESET) == DISP_CHANGE_SUCCESSFUL;
	SetWindowPos(hwnd, HWND_NOTOPMOST, windowX, windowY, windowedWidth + windowedPaddingX, windowedHeight + windowedPaddingY, SWP_SHOWWINDOW);
	ShowWindow(hwnd, SW_RESTORE);

	return isChangeSuccessful;
}

void SetTransparency(HWND hwnd, std::uint8_t Transperancy){
	long wAttr = GetWindowLong(hwnd, GWL_EXSTYLE);
	SetWindowLong(hwnd, GWL_EXSTYLE, wAttr | WS_EX_LAYERED);
	SetLayeredWindowAttributes(hwnd, 0, Transperancy, 0x02);
}