
#ifndef TALKY_MISA
#define TALKY_MISA

#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>

using namespace std;

#include "DWN_utilityconstants.h"

const int MAX_STATEMENT_LENGTH = 100;
 //Maximum to-display string statement to be expected.

//String operations---

//Simple character analysis---

bool Is_Number(char);
 //Returns true if the character passed is a numerical character
 //i.e. between 0 and 9.
bool Has_Number(const char*);
 //Returns true if the string has a number within it (at any location).
bool Has_Number_At_End(const char*);
 //Returns true if the string is terminated by a number (e.g. "Me100 "),
 //ignoring insigificant characters such as white space.

bool Is_Text(char);
 //Returns true if the character passed is a lowercase or uppercase character.
bool Is_Capital_Letter(char);
 //Returns true if the character passed is an uppercase character.
bool Is_Uppercase_Letter(char);
 //Returns true if the character passed is an uppercase character.
bool Is_Undercase_Letter(char);
 //Returns true if the character passed is an undercase character.
char Make_Undercase(char);
 //Returns the passed character with no changes, unless the passed character
 //is a capital letter, in which case it is made to an undercase letter.
bool Contains_At_Least_Two_Capitals(const char*);
 //Returns true if the string contains at least two capital characters.
bool Contains_Two_Capitals_At_Start(const char*);
 //Returns true if the first two characters of the string are capitals.
bool Contains_Capital_Then_Number_At_Start(const char*);
 //Returns true if the first two characters of the string are a 
 //capital and then number.
bool Is_Phrase_End(char);
 //Returns true if the character passed is a phrase termination
 //character, i.e. space, end-of-line, or NULL character.

int Find_First_Meaningful_Char_From_End(const char* ze_string);
 //Returns the index of the first character from the string end that is not
 //a space, end-of-line, or NULL. Returns 0 if no meaningful character is found.
int Find_First_Integer_From_End(const char* ze_string);
 //Returns the index of the first number character from the string end.
 //Returns -1 if no number character is found.
int Find_First_FPChar_From_End(const char* ze_string);
 //Returns the index of the first number character or period from the string end.
 //Returns -1 if no number character or period character is found.

void String_Copy_StopAt_NotInitial_Captial(char*, const char*);
 //Copies string over, but stops string copying at the first captial letter
 //that is not the first character in the string.

//Number retrieval/insertion---

int Get_Positive_Number(const char*);
 //Returns the positive number that terminates a string in integer form.
 //Returns 0 if no number is found. 
int Get_Number(const char*);
 //Returns the number that terminates a string in integer form.
 //Returns 0 if no number is found.
double Get_FP_Number(const char*);
 //Returns the floating-point number that terminates a string.
 //Returns 0.00 if no number is found.
int Get_Numbers(const char*, int*, int);
 //Does as above, but obtains multiple numbers in a row (e.g. 2 2 2).
 //Returns the number of arguments obtained.
int Get_FP_Numbers(const char*, double*, int);
 //Does as above, but obtains multiple numbers in a row. Returns
 //the number of arguments obtained.
void Add_Number(char*, int);
 //Adds the passed number to the end of the character string, if
 //there is no number already there.
void Remove_Number(char*);
 //Removes the number that terminates a string.
int String_Size(int); 
 //Returns the size of an integer in terms of c-string characters needed
 //to represent it. e.g. -1004 has a string size of 5.

//Number alteration---

void Change_Integer(char*, int);
 //Finds the first integer in the passed string and modifies it by the passed value.
 //Note: WILL NOT add digits to the line.
void Change_Integers(char*, int*, int);
 //Alters integers in the passed string by the passed values.
 //Note: WILL NOT add digits to the line.
void Change_FP(char*, double);
 //Finds the first number (or decimal point) in the passed string and modifies 
 //the FP value by the passed value. Note: WILL NOT add digits to the line.
void Change_FPS(char*, double*, int);
 //Alters FP values in the passed string by the passed values.
 //Note: WILL NOT add digits to the line.

//Extension logic---
 
void Get_Extension(char*, const char*);
 //Stores the extension of a file name. Result = NULL if there is no extension.
void Remove_Extension(char*);
 //Removes the file extension, if present, from the passed string.
void Remove_Root(char*);
 //Removes the folder root at the start of the string.
void Copy_With_No_Extension(char*, const char*);
 //Copies a string (latter) into another string space (former), not including
 //file extensions. 
void Add_File_Extension(char*, const char*);
 //Adds an extension string to the end of a name string. Pre-extension
 //period should not be present.
void Make_File_Name(char*, const char*);
 //Add an extension string to the end of a name string.
 //Also adds "Results/" to the start of the file name.
void Make_File_Name(char*, const char*, const char*);
 //Adds an extension string to the end of a name string and places the result
 //in another string. Also adds "Results/" to the start of the file name.
void Change_File_Type(char*, const char*);
 //This function replaces the file extension in a string with a different extension.
 //e.g. "File.txt" ----> "File.goc"

//Command checking---

bool Strict_Command_Check(const char*, const char*);
 //Compares two strings for equality, up to the length of the command string.
 //Requires the next character in the test string to be insignificant 
 //(e.g. a space). For example, COMMAND = "HIP" works with "HIP " but not "HIPS".
bool Command_Check(const char*, const char*);
 //Compares two strings for equality, up to the length of the command string.
 //Returns FALSEV if the last call of this function processed a '/' character.
bool Command_Check_Set_BooleanT(const char*, const char*, bool&);
 //Performs a command check. If the result is TRUEV, the boolean is set
 //to TRUEV. No change otherwise.
bool Command_Check_Set_BooleanF(const char*, const char*, bool&);
 //Performs a command check. If the result is TRUEV, the boolean is set
 //to FALSEV. No change otherwise.
bool Name_Check(const char*, const char*);
 //Returns TRUEV if the entire first string passed matches the
 //entire second string passed. "Entire content" here means all characters up
 //to (and not including) a space or end-of-line or null character, for both strings.
bool Name_Check_To_Index(const char*, const char*);
 //Performs a name check with versions of the passed strings that have terminating
 //numbers removed from them. e.g. "Au3" = "Au5" because "Au" = "Au."
int Find_String_Match(const char*, char**, int);
 //Looks through an array of character strings for a name check match to
 //the passed string. Returns that string's index or -1 if no string is found.

//Read pointer movement---

void Move_Pointer_NextLine(const char*&);
 //Moves a string pointer to the beginning of the next line.
void Move_Pointer_NextLine(const char*&, char);
 /*Moves a string pointer to the next line of text in the character array
   and then moves the pointer to the first instance of the passed character.
   If that character is not found, the pointer is placed at the end-of-line
   on the current line of text (i.e. not the next line discussed).*/
void Move_Pointer_NextLine(const char*&, char, int);
 //Moves the pointer to the next line and then moves it to the first instance
 //of a passed character and then an integer number of characters forward.

void Move_Pointers_NextLine(const char*&, char*&);
 //Moves both pointers to the beginning of the next line.
void Move_Pointers_NextLine(const char*&, char*&, char, int);
 //Moves both points to the beginning of the next line, and then assigns
 //the second pointer to an integer number of characters ahead of a sought character.
 //If character is not found, both pointers stay at the current line.

void MPNL_Get_Positive_Number(const char*&, char, int&);
 //Performs an MPNL operation and then attempts to read a positive integer
 //from the next line of text (if search character is found).
void MPNL_Get_Number(const char*&, char, int&);
 //Performs an MPNL operation and then attempts to read an integer
 //from the next line of text (if search character is found).
void MPNL_Get_FP_Number(const char*&, char, double&);
 //Performs an MPNL operation and then attempts to read a floating-point value
 //from the next line of text (if search character is found).

void MPNL_Get_Numbers(const char*&, char, int*, int);
 //Performs an MPNL operation and then attempts to read integers
 //from the next line of text (if search character is found).
void MPNL_Get_FP_Numbers(const char*&, char, double*, int);
 //Performs an MPNL operation and then attempts to read floating-point values
 //from the next line of text (if search character is found).
void MPNL_Get_Names(const char*&, char, int, char**, int);
 //Performs an MPNL operation and then attempts to read a series
 //of phrases (if search character is found).

void MPSNL_Get_Numbers(const char*&, char*&, char, int, int*, int);
 //Performs an MPSNL (Move PointerS Next Line) operation and then attempts to 
 //read integers from the line (if search character is found).
void MPSNL_Get_FP_Numbers(const char*&, char*&, char, int, double*, int);
 //Performs an MPSNL (Move PointerS Next Line) operation and then attempts to 
 //read floating-point values from the line (if search character is found).
 //Example command --> "Add Si3N4 0.5 4.1 3.2"

void Find_Next_Line_Start(const char*&);
 //Moves a string pointer forward to the start of the next line, that being
 //the first character after an end-of-line character. 
void Find_Char_Next_Line(const char*&, char); 
 /*Moves a string pointer forward to the next occurence of a character occuring
   after a line feed character, ignoring any NULL terminating characters.
   If the character is not found, the pointer is moved to the end of the 
   line it had been on before the function call.
   If argument is '\n', the end of the current line is returned.*/
void Skip_Phrase(const char*&);
 //Moves the pointer to one character in front of the
 //next space character. Search starts at one character
 //in front of the passed reference point.
 
//Combined reader movement with command checking---

bool MPNL_CommandCheck(const char*&, const char*);
 //Performs an MPNL operation and then looks for the passed command. If not found,
 //the pointer is left at the end of the initial line. If found, it remains at
 //the first character of the next line.
bool MPNL_CommandCheck(const char*&, const char*, int);
 //Performs an MPNL operation and then looks for the passed command. If not found,
 //the pointer is left at the end of the initial line. If found, the pointer
 //is moved an integer value of characters ahead of the first character of the command.
bool MPNL_CommandCheck_SC(const char*&, const char*);
 //Performs an MPNL operation and then looks for the passed command. If not found,
 //the pointer is left at the end of the initial line. If found, the pointer
 //is moved to two spaces in front of the end of the command.

bool MPSNL_CommandCheck(const char*&, char*&, const char*);
 //Performs an MPSNL operation and then looks for the passed command. If not found,
 //the pointers are left at the end of the initial line.
bool MPSNL_CommandCheck(const char*&, char*&, const char*, int);
 //Performs an MPSNL operation and then looks for the passed command. If not found,
 //the pointers are left at the end of the initial line. If found, the second pointer
 //is moved an integer value of characters ahead of the command.
bool MPSNL_CommandCheck_SC(const char*&, char*&, const char*);
 //Performs an MPSNL operation and then looks for the passed command. If not found,
 //the pointer is left at the end of the initial line. If found, the second pointer
 //is moved to two spaces in front of the end of the command.

//The following operations combine an MPNL or MPSNL operation wtih retrieval
//of numbers or names to be found after the command---

bool MPNL_CommandCheck_GetPosNum(const char*&, const char*, int&);
bool MPNL_CommandCheck_GetNum(const char*&, const char*, int&);
bool MPNL_CommandCheck_GetFPNum(const char*&, const char*, double&);
bool MPNL_CommandCheck_GetFPNum(const char*&, const char*, double&, double);
bool MPNL_CommandCheck_GetNums(const char*&, const char*, int*, int);
bool MPNL_CommandCheck_GetFPNums(const char*&, const char*, double*, int);
bool MPNL_CommandCheck_GetNames(const char*&, const char*, int, char**, int);
bool MPSNL_CommandCheck_GetNums(const char*&, char*&, const char*, int, int*, int);
bool MPSNL_CommandCheck_GetFPNums(const char*&, char*&, const char*, int, double*, int);

//File-to-character logic---

void Get_Characterized_File(ifstream&, char*, const int, int&);
 //Dumps a file into a character array, changing line feeds to line
 //feeds preceded by NULL characters.
void String_Replace(char*, int, int, const char*, int);
  //Replaces n-characters of a string with m-characters of another string.
  //Assumes there is enough space in the character string to deal with an 
  //increase in overall string size.
 
//File logic---

void Skip_Phrases(ifstream&, int);
 //Skips past an integer number of phrases during a file read operation.
 //For example, 5 skips would move past "HEY MY NAME IS DAVID!".
void Load_Values(ifstream&, int*, int);
 //Reads an integer number of integers during a file read operation.
void Load_Values(ifstream&, double*, int);
 //Reads an integer number of floating-point values during a file read operation.
void Write_Values(ofstream&, const int*, int);
 //Stores an integer number of integers during a file write operation.
void Write_Values(ofstream&, const double*, int);
 //Stores an integer number of floating-point values during a file write operation.
void Write_Intro(ofstream&, const char*);
 //Stores a file introduction string during a file write operation.

//Character output logic---

//General statement output---

void Show_Statement(const char*);
 //Show a string to the user.
void Show_Statement(const char*, int);
 //Show a string with appended number to the user.
void Show_Statement(const char*, double);
 //Show a string with appended number to the user.
void Show_Statement(const char*, int, const char*);
 //Show a string with appended number and then string to the user.
void Show_Statement(const char*, const char*);
 //Show a set of connected strings to the user.
void Show_Statement(const char*, const char*, const char*);
 //Show a set of connected strings to the user.

void Show_Warning(const char*);
 //Shows a string to the user, as a warning statement.
void Show_Warning(const char*, const char*, const char*);
 //Shows a set of connected strings to the user, as a warning statement.
void Show_Congratulations(const char*);
 //Special congratulations statement.

//Number output---

void Show_XYZ(const double*, double);
 //Shows a set of coordinates on screen, multiplying by the scaling factor.
void Write_XYZ(ofstream&, const double*, double);
 //Writes a set of coordinates to file, multiplying by the scaling factor.
void Read_XYZ(ifstream&, double*, double);
 //Reads in a set of coordinates from file, multiplying by the scaling factor.

//Number output with column justification---

void Output_Spaces(ofstream&, int);
 //Outputs a passed number of spaces to file.
void Left_Justify(ofstream&, double, int, int);
 //Outputs a value with specificed precision and left justifies it.
void Left_Justify(ofstream&, int, int);
 //Outputs an integer and left justifies it.
void Right_Justify(ofstream&, double, int, int);
 //Outputs a value with specificed precision and right justifies it.
void Right_Justify(ofstream&, int, int);
 //Outputs an integer and right justifies it.

#endif
