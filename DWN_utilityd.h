
#ifndef FIXIT_JOSH
#define FIXIT_JOSH

#include <iostream>

using namespace std;

#include "DWN_utilityconstants.h"
#include "DWN_utilitys.h"
#include "DWN_utilitym.h"
#include "DWN_utilityc.h"

//Pseudo-global variables---

void Inc_TempA();
void Inc_TempB();
void Inc_TempC();
 //Increase the value of a global variable
 //by one.
void Show_TempA();
void Show_TempB();
void Show_TempC();
 //Prints the global variable to screen.

//Debug functions---

void Say_Hi();
 //Outputs "Hello!" to screen.
void Show_String(const char*);
 //Outputs a cstring.
void Show_Num(int);
void Show_Num(double);
 //Outputs a number to screen.
void Show_Nums(const int*);
void Show_Nums(const double*);
 //Outputs three numbers to screen.
 
void Show_File(const char*, bool);
 //Shows a file to the screen, effectively showing EVERY character, using
 //ASCII integers to represent the characters of otherwise invisible
 //characters.

#endif