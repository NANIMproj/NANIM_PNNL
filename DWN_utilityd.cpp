
#include "stdafx.h"

#include "DWN_utilityd.h"


//Psuedo-Global variables

int TEMPA = 0;
int TEMPB = 0;
int TEMPC = 0;
 //Global variables just for the purpose of debugging.

void Inc_TempA()
  { 
  ++TEMPA;
  }

void Show_TempA()
  {
  cout << endl << "TEMPORARY A'S VALUE IS: " << TEMPA;
  }

void Inc_TempB()
  { 
  ++TEMPB;
  }

void Show_TempB()
  {
  cout << endl << "TEMPORARY B'S VALUE IS: " << TEMPB;
  }

void Inc_TempC()
  { 
  ++TEMPC;
  }

void Show_TempC()
  {
  cout << endl << "TEMPORARY C'S VALUE IS: " << TEMPC;
  }

//Debug function---

void Say_Hi()
     {
     cout << endl << "HEY THERE! HOPE YOU SEE ME ELSE THE SHIT HAS HIT THE FAN.";        
     }

void Show_String(const char* cstring)
     {
     cout << endl << cstring;                  
     }

void Show_Num(int num)
     {
     cout << endl << num;
     }

void Show_Num(double num)
     {
     cout << endl << num;                
     }

void Show_Nums(const int* nums)
     {
     cout << endl << nums[0] << " " << nums[1] << " " << nums[2];               
     }

void Show_Nums(const double* nums)
     {
     cout << endl << nums[0] << " " << nums[1] << " " << nums[2]; 
     }


void Show_File(const char* file_name, bool pause_mode)
	 {
	 ifstream show_file(file_name, ios::in);
     if (show_file.is_open())
		{
		char test;
		bool is_letter;
		bool is_number;
		bool is_special;

		cout << endl << "File Read Starts on Next Line: " << endl;
		while (!show_file.eof())
			{
			show_file.get(test);
			is_letter = Is_Text(test);
			is_number = Is_Number(test);
			is_special =   ( (test >= 32) && (test <= 47) ) 
				        || ( (test >= 58) && (test <= 64) )
						|| ( (test >= 91) && (test <= 96) )
						|| ( (test >= 123) && (test <= 126) );
			 //Additional ASCII characters that would be
			 //output as characters, not ASCII integers.
			if (is_letter || is_number || is_special)
			   {
			   cout << test;
			   }
			else
			   {
			   cout << int(test) << " ";
			   }

			if (test == '\n')
			 //Process the end of line.
			   {
			   cout << endl;
			   if (pause_mode)
			     {
				 system("PAUSE");
			     }
			   }

			}
		cout << endl << "End File Read!";
	    show_file.close();
		}
	 }