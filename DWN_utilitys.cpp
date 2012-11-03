
#include "stdafx.h"

#include "DWN_utilitys.h"

//String operations---

bool Is_Number(char num_character)
     {
     bool is_number = ( (num_character >= '0') && (num_character <= '9') );              
     return is_number;
     }

bool Has_Number(const char* ze_string)
       {
       bool has_number = FALSEV;
       int char_index = Find_First_Integer_From_End(ze_string);
       if (char_index != -1)
          {
          has_number = TRUEV;            
          }
       return has_number;             
       }

bool Has_Number_At_End(const char* ze_string)
       {
       bool has_number = FALSEV;
       int char_index = Find_First_Meaningful_Char_From_End(ze_string);
       if (Is_Number(ze_string[char_index]))
          {
          has_number = TRUEV;            
          }
       return has_number;    
       }

bool Is_Text(char text_character)
     {
     bool is_lowercase =  Is_Undercase_Letter(text_character);
     bool is_uppercase =  Is_Uppercase_Letter(text_character);
     bool is_text = is_lowercase || is_uppercase;
     return is_text;           
     }

bool Is_Capital_Letter(char ze_char)
	{
    bool is_cap = ( (ze_char >= 'A') && (ze_char <= 'Z') );              
    return is_cap;
	}

bool Is_Uppercase_Letter(char ze_char)
	{
    bool is_cap = ( (ze_char >= 'A') && (ze_char <= 'Z') );              
    return is_cap;
	}

bool Is_Undercase_Letter(char ze_char)
	{
    bool is_und = ( (ze_char >= 'a') && (ze_char <= 'z') );              
    return is_und;
	}

char Make_Undercase(char ze_char)
	{
	char undercase_char = ze_char;
	if (Is_Uppercase_Letter(ze_char))
	   {
	   undercase_char -= ('A' - 'a' );
	   }
	return undercase_char;
	}

bool Contains_At_Least_Two_Capitals(const char* ze_string)
	{
	bool contains_two = FALSEV;
	int string_length = strlen(ze_string);
	int cap_count = 0;
	for (int a = 0; a < string_length; ++a)
	    {
		if (Is_Capital_Letter(ze_string[a]))
			{ 
			++cap_count;
			if (cap_count > 1)
			   {
			   contains_two = TRUEV;
			   a = string_length;
			   }
			}
	    }
	return contains_two;
	}

bool Contains_Two_Capitals_At_Start(const char* ze_string)
    {
    bool contains_two = FALSEV;
	if ( (strlen(ze_string) > 1) && Is_Capital_Letter(ze_string[0])
		  && Is_Capital_Letter(ze_string[1]) )
	   {
	   contains_two = TRUEV;
	   }
	return contains_two;
    }   

bool Contains_Capital_Then_Number_At_Start(const char* ze_string)
	{
	bool contains_capnum = FALSEV;
	if ((strlen(ze_string) > 1))
	   {
	   if (Is_Capital_Letter(ze_string[0]) && Is_Number(ze_string[1]))
	      {
		  contains_capnum = TRUEV;
	      }
	   }
	return contains_capnum;
	}

bool Is_Phrase_End(char ze_char)
     {
     bool is_end = ( (ze_char == ' ') || (ze_char == '\0') || 
                     (ze_char == '\n') ); 
     return is_end;                      
     }

int Find_First_Meaningful_Char_From_End(const char* ze_string)
    {
    int string_size = strlen(ze_string);
	int char_index = 0;
	if (string_size > 0)
	   {
       char_index = string_size - 1;
       char test_char = ze_string[char_index];
       while ( Is_Phrase_End(test_char) && (char_index > 0) )
        //Move past white space.
          {  
          --char_index;
          test_char = ze_string[char_index];
          }
	   }
    return char_index;                                         
    }

int Find_First_Integer_From_End(const char* ze_string)
    {
    int string_size = strlen(ze_string);
    int char_index = -1;
	if (string_size > 0)
	   {
	   char_index = string_size - 1;
	   char test_char = ze_string[char_index];
	   while (!Is_Number(test_char) && (char_index > 0))
        //Move past white space.
          {
          --char_index;
          test_char = ze_string[char_index];
          }   
	    if (!Is_Number(test_char))
	      {
		  char_index = -1;
	      }
	    }
    return char_index;                                         
    }

int Find_First_FPChar_From_End(const char* ze_string)
    {
    int string_size = strlen(ze_string);
    int char_index = -1;
	if (string_size > 0)
	   {
	   char_index = string_size - 1;
	   char test_char = ze_string[char_index];
	   while (!( Is_Number(test_char) || (test_char == '.') ) && (char_index > 0))
        //Move past white space.
          {
          --char_index;
          test_char = ze_string[char_index];
          }   
	    if (!( Is_Number(test_char) || (test_char == '.') ))
	      {
		  char_index = -1;
	      }
	    }
    return char_index;                                    
    }

void String_Copy_StopAt_NotInitial_Captial(char* string1, const char* string2)
    {
	int char_index = 0;
	if (string2[0] != '\0')
	 //Copy string over until the string end or a capital is encountered.
	  {
	  do
	    {
	    string1[char_index] = string2[char_index];
		++char_index;
	    }
	  while ( (string2[char_index] != '\0') 
		             && !Is_Capital_Letter(string2[char_index]) );
      }  

	string1[char_index] = '\0';
	 //Place the final string-terminating character.
    }

//Number logic---

int Get_Positive_Number(const char* ze_string)
       {
       if (!Has_Number(ze_string))
          {
          return 0;
          }         
       
       int start_index = Find_First_Integer_From_End(ze_string);
       
       int number = 0;
       int multi_factor = 1;
       int index;
       for (int a = 0; a < MAX_POI_DIGITS; ++a)
           {
           index = start_index - a;
           
           if (index < 0)
              {
              break; 
              }
           
           if (Is_Number(ze_string[index]))
              {
              number += multi_factor * int(ze_string[index] - '0');
              multi_factor *= 10;     
              }
           else
              {
              break;      
              }
           }     
       
       return number;
       }    

int Get_Number(const char* ze_string)
       {
       if (!Has_Number(ze_string))
          {
          return 0;
          }         
       
       int start_index = Find_First_Integer_From_End(ze_string);
       
       int number = 0;
       int multi_factor = 1;
       int index;
       for (int a = 0; a < MAX_I_DIGITS; ++a)
           {
           index = start_index - a;
           
           if (index < 0)
              {
              break; 
              }
           
           if (Is_Number(ze_string[index]))
              {
              number += multi_factor * int(ze_string[index] - '0');
              multi_factor *= 10;     
              }
           else
              {
              break;      
              }
           }     
       
       if ( (index >= 0) && (ze_string[index] == '-') )
        //Number is meant to be negative.
          {
          number *= -1;                  
          }    
       
       return number;
       }    
       
double Get_FP_Number(const char* ze_string)
       {
       if (!Has_Number(ze_string))
          {
          return 0.00;               
          }         
       
       int start_index = Find_First_FPChar_From_End(ze_string);
       
	   //Find decimal point, if present---

       int index, decimal_index;
       bool has_decimal = FALSEV;
       for (int a = 0; a < MAX_FP_DIGITS; ++a)
        //Find the decimal point, if present.
           {
           index = start_index - a; 
           
           if (index < 0)
            //Don't want to explore negative space.
              {
              break;  
              }
           else if (ze_string[index] == '.')
              {
              decimal_index = index;
              has_decimal = TRUEV;
			  break;
              }    
		   else if (!Is_Number(ze_string[index]))
		    //No period found.
			  {
			  break;
			  }
           }	
	          
	   //Evaluate the base 10 exponent for the
	   //first number in the string---

       double number = 0.00;
       double multi_factor = 1.0;
       if (has_decimal)
          {
          int num_dec_digits = start_index - decimal_index;
		   //Get number of digits after the decimal place.
          multi_factor = pow(10.0, -1*num_dec_digits);
           //Determine base 10 power for the last digit in the FP number string.
          }
       
	   //Determine the floating-point value of the string number---

       for (int a = 0; a < MAX_FP_DIGITS; ++a)
           {
           index = start_index - a;
           
		   if (index < 0)
			  {
			  break;
			  }

           if (has_decimal && (index == decimal_index))
            //Do not mathematically evalute the decimal place.
              {
              continue;              
              }
			
           if (Is_Number(ze_string[index]))
              {
              number += multi_factor * double(ze_string[index] - '0');
              multi_factor *= 10.0;     
              }
           else
              {
              break;  
              }
           }     
       
       if ( (index >= 0) && (ze_string[index] == '-'))
        //Number is negative.
          {
          number *= -1.0;                  
          }    
       
       return number;
       }    
  
int Get_Numbers(const char* ze_string, int* numbers, int variable_max)
       {
       //Prepare the test string for integer retrieval---

       char temp_string[MAX_PHRASE_LENGTH];
       char* temp_pointer;
       strcpy(temp_string, ze_string);
       int string_length = strlen(ze_string);
       
       bool is_num_back;
       bool is_not_num_here;
       for (int index = 1; index < string_length; ++index)
        //Separate the string numbers with NULL characters.
           {
           is_num_back = Is_Number(temp_string[index - 1]);
           is_not_num_here = !(Is_Number(temp_string[index]));
           if (is_num_back && is_not_num_here)
            //Found the end of a number string. Separate it with
			//a null-terminating character---
              {
              temp_string[index] = '\0';             
              }
           }

	   //Grab the numbers, one by one, from the prepared string---
           
       bool is_NULL;
       numbers[0] = Get_Number(temp_string);
	    //Grab the first number from the string.
       int count_temp = 1;
       for (int index = 0; index < (string_length - 1); ++index)
        //Process each separated number.
           {
		   if (count_temp == variable_max)
              {  
			  break;
              }

           is_NULL = (temp_string[index] == '\0');
		    //Found the start of a new number in the string.
           if (is_NULL)
              {
              temp_pointer = &(temp_string[index + 1]);
              numbers[count_temp] = Get_Number(temp_pointer);
              ++count_temp;         
              }     
           }
      
       return count_temp;
       }  
       
int Get_FP_Numbers(const char* ze_string, double* numbers, int variable_max)
       {
	   //Prepare the test string for floating-point retrieval---

       char temp_string[MAX_PHRASE_LENGTH];
       char* temp_pointer;
       strcpy(temp_string, ze_string);
       int string_length = strlen(ze_string);
       
       bool is_num_back;
       bool is_not_num_here;
       for (int index = 1; index < string_length; ++index)
        //Separate the numbers with NULL characters.
           {
           is_num_back = Is_Number(temp_string[index - 1]);
           is_not_num_here = !(Is_Number(temp_string[index]) || (temp_string[index] == '.') );
           if (is_num_back && is_not_num_here)
            //Found the end of a number string. 
              {
              temp_string[index] = '\0';             
              }
           }
           
	   //Grab the numbers, one by one, from the prepared string---

       bool is_NULL;
       numbers[0] = Get_FP_Number(temp_string);
	    //Grab the first number.
       int count_temp = 1;
       for (int index = 0; index < (string_length - 1); ++index)
        //Process each separated number.
           {
		   if (count_temp == variable_max)
              {
              break;        
              }

           is_NULL = (temp_string[index] == '\0');
           if (is_NULL)
              {
              temp_pointer = &(temp_string[index + 1]);
              numbers[count_temp] = Get_FP_Number(temp_pointer);
              ++count_temp;         
              }     
           }

       return count_temp;
       }

void Add_Number(char* ze_string, int number)
       {
       if (Has_Number_At_End(ze_string))
		//Don't add a number to a number.
          {
          return;               
          }         

	   //Transform number into a character string---
       
       char num_string[MAX_I_DIGITS + 1];
	   Initialize_Array(num_string, MAX_I_DIGITS + 1, ' ');
        //Extra space is for an eventual NULL terminating character.
       int highest_order_digit = 0;
       int digit, index;
	   int test_number = number;
       for (int a = 0; a < MAX_I_DIGITS; ++a)
        //Create a string-version of the number.
           {
           index = MAX_I_DIGITS - 1 - a;
           if (index < 0)
            //Safety check to prevent breaking memory bounds.
              {
              break;       
              }
           digit = test_number % 10;
           num_string[index] = '0' + digit;
           ++highest_order_digit;
		   test_number /= 10;
           if (test_number == 0)
            //Done converting the number, minus a minus sign.
              {              
              if ( (number < 0) && (index > 0) )
                 {
                 num_string[index - 1] = '-';
                 ++highest_order_digit;         
                 }
              a = MAX_I_DIGITS;  
              }
           }
      
	   //Position number string at the front of the character array---

       int pushback_factor = MAX_I_DIGITS - highest_order_digit;    
       for (int a = 0; a < highest_order_digit; ++a)
        //Move the number into place (at start of C-string).
           {
           num_string[a] = num_string[a + pushback_factor];     
           }	
       num_string[highest_order_digit] = '\0';
        //Add NULL terminating character.
       strcat(ze_string, num_string);
        //Copy number onto the target string.
       }      

void Remove_Number(char* ze_string)
     {
     if (!Has_Number_At_End(ze_string))
          {
          return;               
          }         
       
     int string_size = strlen(ze_string);
     int offset = 1;
     int index;
     char current_digit;     
     do
      //Find the first character, moving back from the end of the string, 
      //that is not a number or minus sign.
      {
      ++offset;
      index = string_size - offset;
      current_digit = ze_string[index];               
      }
     while ( (Is_Number(current_digit)) || (current_digit == '-') );
       
     ze_string[index + 1] = '\0';
      //Remove number by placing a NULL terminating character before it.                            
     } 

int String_Size(int value)
	 {
	 int size = 1;
	 while ( (value / 10) != 0)
	      {
		  ++size;
		  value /= 10;
	      }
	 if (value < 0) ++size;
	 return size;
	 }
    
//Number alteration--

void Change_Integer(char* line, int change_val)
     {
	 //Determine the new value of the integer---

	 int orig_val = Get_Number(line);
	 int new_val = orig_val + change_val;

	 //Characterize the integer string---

	 int val_size = 0;
	 int number_location;
	 bool num_found = FALSEV;
	 for (int a = 0; a < int(strlen(line)); ++a)
	     {
		 if (Is_Number(line[a]) && !num_found)
		  //Found starting index for the integer.
		    {
			number_location = a;
			num_found = TRUEV;
		    }
		 if (Is_Number(line[a]) && num_found)
		  //Get string length of the integer.
		    {
			++val_size;
		    }
		 if (!Is_Number(line[a]) && num_found)
		    {
			break;
		    }
	     }

	 //Change the integer string to the new value---

	 int division_factor = 1;
	 int abs_val = abs(new_val);
	 for (int b = 0; b < val_size; ++b)
	     {
		 line[number_location + val_size - 1 - b] = '0' + 
			  char(abs_val/division_factor % 10);
		 division_factor *= 10;
	     }

	 //Check for a minus sign that needs to be removed---

	 if ( (number_location != 0) && (line[number_location - 1] == '-') && (new_val > 0) )
	  //If number has been changed from negative value to positive value,
	  //the minus sign must be removed. Reverse case is not allowed here.
	     {
         line[number_location - 1] = '0';
		  //Replace with a zero. Keeps string length of the number the same.
	     }

	 //Check for a new value that is larger in string space than the original value---

	 int num_temp = Get_Number(line);
	 if (new_val != num_temp)
	     {
		 Show_Warning("AUTOMATION OF FILE TEXT REQUIRES MORE SPACE IN INTEGER NUMBER STRINGS!");
		 Change_Integer(line, -1*num_temp);
		  //Zero the integer string in the case of overflow.
	     }
	 }

void Change_Integers(char* line, int* change_vals, int num_vals)
     {
	 int num_count = 0;
	 bool found_num, set_space_to_null, is_minus;
	 char* line_pointer;
	 int end_number_location;

	 //Go through string, finding integers to apply changes to, and applying
	 //those changes---

	 for (int a = 0; a < int(strlen(line)); ++a)
		{
		if (a < (int(strlen(line)) - 1))
		   {
		   is_minus = (line[a] == '-') && (Is_Number(line[a + 1]));
		   }
		else
		   {
		   is_minus = FALSEV;
		   }
		if (Is_Number(line[a]) || is_minus)
		   {
		   found_num = TRUEV;
		   if ( (a != 0) )
			//Only accept numbers that are set apart by white space.
		      {
			  if (line[a - 1] != ' ')
			     {
				 found_num = FALSEV;
			     }
		      }
		   if (found_num)
			//Found a number to apply the change to.
		      {
			  line_pointer = &(line[a]);
			  set_space_to_null = FALSEV;
			  for (int b = a; b < int(strlen(line)); ++b)
			     {
			     if (line[b] == ' ')
				    {
					end_number_location = b;
					set_space_to_null = TRUEV;
					break;
				    }
			     }
			  if (set_space_to_null)
			     {
				 line[end_number_location] = '\0';
			     }
			  Change_Integer(line_pointer, change_vals[num_count]);   
			  if (set_space_to_null)
				//Undo the change.
			     {
				 line[end_number_location] = ' ';
			     }
			  ++num_count;
		      }
		   }

		if (num_count == num_vals)
		   {
		   break;
		   }
	    }
     }

void Change_FP(char* line, double change_val)
     {
	 //Determine the new value of the variable---

	 double orig_val = Get_FP_Number(line);
	 double new_val = orig_val + change_val;

	 //Characterize the floating-point string---

	 int val_size = 0;
	 int after_dec = 0;
	 int before_dec = 0;
	 int number_location;
	 bool num_found = FALSEV;
	 bool period_found = FALSEV;
	 bool is_num_dec;
	 for (int a = 0; a < int(strlen(line)); ++a)
	     {	
		 is_num_dec = (line[a] == '.') && (Is_Number(line[a + 1]));
		 if (is_num_dec)
		    {
			period_found = TRUEV;
		    }
		 if ( (Is_Number(line[a]) || period_found) && !num_found)
		  //Start of the number string has been found.
		    {
			number_location = a;
			num_found = TRUEV;
		    }
		 if (Is_Number(line[a]) && num_found)
	      //Determine how many digits are before and after the decimal place.
		    {
			++val_size;
			if (period_found)
			  {
			  ++after_dec;
			  }
			else
			  {
			  ++before_dec;
			  }
		    }
		 if (!(Is_Number(line[a]) || is_num_dec) && num_found)
		  //End of the number string has been encountered.
		    {
			break;
		    }
	     }

	 //Change the floating-point string to the new value---

	 int dif_from_ones_place = 0;
	 int char_index;
	 double temp_shift;
	 double abs_val = abs(new_val);
	 for (int b = 0; b < val_size; ++b)
	     {
		 char_index = b;
		 dif_from_ones_place = (before_dec - 1) - b;
		 temp_shift = pow(10.0, double(dif_from_ones_place));
		 if (b >= before_dec)
		  //Skip the decimal point.
		    {
			++char_index;
		    } 
		 line[number_location + char_index] = '0' + 
			  char( int(abs_val/temp_shift + FP_ERROR_FIX) % 10 );
	     }

	 //Check for the need to remove a minus sign---

	 if ( (number_location != 0) && (line[number_location - 1] == '-') && (new_val > 0.0) )
	  //If number has been changed from negative value to positive value,
	  //the minus sign must be removed. Reverse case is not allowed here.
	     {
         line[number_location - 1] = '0';
		  //Replace with a zero. Keeps the string length the same.
	     }

	 //Check for a failed change to the FP value---

	 double num_temp = Get_FP_Number(line);
	 if (!Check_FP_Equality(num_temp, new_val))
	     {
		 Show_Warning("AUTOMATION OF FILE TEXT REQUIRES MORE SPACE FOR FP VALUES!");
		 Change_FP(line, -1.0*num_temp);
		  //Zero the FP value.
	     }
     }

void Change_FPS(char* line, double* change_vals, int num_vals)
     {
	 int num_count = 0;
	 char* line_pointer;
	 int end_number_location;
	 bool found_num, set_space_to_null, is_num_dec, is_num_minus_sign;

	 //Go through string, finding integers to apply changes to, and applying
	 //those changes---

	 for (int a = 0; a < int(strlen(line)); ++a)
		{
		if (a < (int(strlen(line)) - 1))
		   {
           is_num_dec = (line[a] == '.') && (Is_Number(line[a + 1]));
		   is_num_minus_sign = (line[a] == '-') && 
			                   ( Is_Number(line[a + 1]) || (line[a + 1] == '.') );
		   }
		else
		   {
		   is_num_dec = is_num_minus_sign = FALSEV;
		   }
		if (Is_Number(line[a]) || is_num_dec || is_num_minus_sign)
		   {
		   found_num = TRUEV;
		   if ( (a != 0) )
			//Only consider a number that is preceded by white space (and thus set apart).
		      {
			  if (line[a - 1] != ' ')
			     {
				 found_num = FALSEV;
			     }
		      }
		   if (found_num)
		      {
			  line_pointer = &(line[a]);
			  set_space_to_null = FALSEV;
			  for (int b = a; b < int(strlen(line)); ++b)
			     {
			     if (line[b] == ' ')
				    {
					end_number_location = b;
					set_space_to_null = TRUEV;
					break;
				    }
			     }
			  if (set_space_to_null)
			     {
				 line[end_number_location] = '\0';
			     }
			  Change_FP(line_pointer, change_vals[num_count]);   
			  if (set_space_to_null)
				//Undo the change.
			     {
				 line[end_number_location] = ' ';
			     }
			  ++num_count;
		      }
		   }

		if (num_count == num_vals)
		   {
		   break;
		   }
	    }
     }

//Extension logic--

void Get_Extension(char* extension, const char* file_name)
     {
     int string_length = strlen(file_name);
     if (file_name[string_length - 4] == '.')
        {
        strcpy(extension, &(file_name[string_length - 3]));                     
        }   
     else
        {
        extension[0] = '\0'; 
        }
     }

void Remove_Extension(char* ze_string)
     {
     char done1 = '\0';
     char done2 = '.';
      //Characters that indicate the string extension is found.
      
     int index = 0; 
     while ( (ze_string[index] != done1) && (ze_string[index] != done2) )
      {
      ++index;     
      } 
     
     ze_string[index] = '\0';
     }

void Remove_Root(char* ze_string)
     {
	 char search_char = '\\';
	 int string_length = strlen(ze_string);
	 int search_index = -1;
	 for (int a = 0; a < string_length; ++a)
	  //Find last instance of the search character (\).
	     {
		 if (ze_string[a] == search_char)
		    {
		    search_index = a;
		    }
	     }

	 int pushback = search_index + 1;
	 for (int a = 0; a < (string_length - pushback); ++a)
	     {
		 ze_string[a] = ze_string[a + pushback];
	     }

	 ze_string[string_length - pushback] = '\0';
     }

void Copy_With_No_Extension(char* dest, const char* source)
     {
     strcpy(dest, source);
     Remove_Extension(dest);                       
     }

void Add_File_Extension(char* ze_string, const char* extension)
     {
	 int string_length = strlen(ze_string);
     if (ze_string[string_length - 4] == '.')
      //Don't want to add an extension to an extension.
        {
        return;                         
        }
     ze_string[string_length] = '.';
     ze_string[string_length + 1] = '\0';
     strcat(ze_string, extension);    
     }
     
void Make_File_Name(char* ze_string, const char* extension)
     {
	 Add_File_Extension(ze_string, extension);

	 char temp_string[100];
	 temp_string[0] = '\0';
	 if (!Command_Check(ze_string, "Results\\"))
	  //Get proper root added.
	    {
	    strcpy(temp_string, "Results\\\0");
		}
	 strcat(temp_string, ze_string);
	 strcpy(ze_string, temp_string);
     }

void Make_File_Name(char* ze_string, const char* file, const char* extension)
     {
	 strcpy(ze_string, file);
	 Add_File_Extension(ze_string, extension);

	 char temp_string[100];
	 temp_string[0] = '\0';
	 if (!Command_Check(ze_string, "Results\\"))
	  //Get proper root added.
	    {
	    strcpy(temp_string, "Results\\\0");
		}
	 strcat(temp_string, ze_string);
	 strcpy(ze_string, temp_string);
     }


void Change_File_Type(char* ze_string, const char* new_extension)
     {
     Remove_Extension(ze_string);
     Make_File_Name(ze_string, new_extension);
     }

//Command checking---

bool Strict_Command_Check(const char* file_pointer, const char* command_pointer)
     {
     bool check = FALSEV;                    
     int command_length = strlen(command_pointer);
     if(strncmp(file_pointer, command_pointer, command_length) == 0 )
         {
         check = TRUEV;                      
         }  
	 bool phrase_end = Is_Phrase_End(file_pointer[command_length]);
	 if (!phrase_end)
	    {
		check = FALSEV;
	    }
     return check;      
     }

bool Command_Check(const char* file_pointer, const char* command_pointer)
     {
     static char last_character = '+';
     bool check = FALSEV;                    
     int command_length = strlen(command_pointer);
     if(strncmp(file_pointer, command_pointer, command_length) == 0 )
         {
         check = TRUEV;                      
         }         
     if (last_character == '/')
      //Don't want to process commented lines.
         {
         check = FALSEV;               
         }    
     last_character = file_pointer[0];
     return check;      
     }

bool Command_Check_Set_BooleanT(const char* file_pointer, const char* command_pointer,
		                       bool& boolean_pam)
     {
	 bool check = Command_Check(file_pointer, command_pointer);
	 if (check)
	    {
		boolean_pam = TRUEV;
	    }
     return check;      
     }

bool Command_Check_Set_BooleanF(const char* file_pointer, const char* command_pointer,
		                       bool& boolean_pam)
     {
	 bool check = Command_Check(file_pointer, command_pointer);
	 if (check)
	    {
		boolean_pam = FALSEV;
	    }
     return check;      
     }

bool Name_Check(const char* name_pointerA, const char* name_pointerB)
     {
     bool check = FALSEV;
     int index = -1;
     char test_charA, test_charB;
     bool same_char, phrase_end, phrase_endA, phrase_endB;
     do
	       {
           ++index;
           test_charA = name_pointerA[index];
           test_charB = name_pointerB[index];
            
           same_char = (test_charA == test_charB);
           phrase_endA = Is_Phrase_End(test_charA);
           phrase_endB = Is_Phrase_End(test_charB);
           phrase_end = phrase_endA || phrase_endB;
           }
     while (same_char && !phrase_end);
      //Loop through the strings until unlike characters or a phrase end
      //(space, end-of-line, or NULL character) is encountered.
     
     if (phrase_endA && phrase_endB)
      //Must find a phrase-ender for both strings at the same location for this
      //function to accept the names as truly equivalent!
        {
        check = TRUEV;            
        }
             
     return check;        
     }

bool Name_Check_To_Index(const char* name_pointerA, const char* name_pointerB)
	 {
	 char test_stringA[100];
	 strcpy(test_stringA, name_pointerA);
	 Remove_Number(test_stringA);
	 char test_stringB[100];
	 strcpy(test_stringB, name_pointerB);
	 Remove_Number(test_stringB);
	 
	 bool check = Name_Check(test_stringA, test_stringB);
	 return check;
	 }

int Find_String_Match(const char* string, char** array, int array_size)
    {
    int index = -1;
    for (int a = 0; a < array_size; ++a)
        {
        if (Name_Check(string, array[a]))
           {
           index = a;  
           a = array_size;
           }
        }   
    return index;            
    }

//Read pointer movement---

void Move_Pointer_NextLine(const char*& string_pointer)
     {
     Find_Next_Line_Start(string_pointer);                            
     }

void Move_Pointer_NextLine(const char*& string_pointer, char search)
     {
     Find_Char_Next_Line(string_pointer, search);
     }

void Move_Pointer_NextLine(const char*& string_pointer, char search, int index)
     {
     Move_Pointer_NextLine(string_pointer, search);
	 if ( (string_pointer[0] == '\n') && (search != '\n') )
      //If search character is not found, keep the pointer at the end of line.
        {
        index = 0;                   
        }
     string_pointer = &(string_pointer[index]);                                
     }

void Move_Pointers_NextLine(const char*& string_pointer, char*& second_pointer)
     {
     Move_Pointer_NextLine(string_pointer);    
     second_pointer = const_cast<char*>(string_pointer);                        
     }

void Move_Pointers_NextLine(const char*& string_pointer, char*& assign_pointer, 
                            char search, int index)
     {
     Move_Pointer_NextLine(string_pointer, search);
     if ( (string_pointer[0] == '\n') && (search != '\n') )
      //If  search character is not found, place the assignment pointer
      //at the end of line.
        {
        index = 0;                   
        }
     assign_pointer = const_cast<char*>(&(string_pointer[index]));                        
     }

void MPNL_Get_Positive_Number(const char*& string_pointer, char search, int& num)
     {
     Move_Pointer_NextLine(string_pointer, search);
     num = Get_Positive_Number(string_pointer);                               
     }

void MPNL_Get_Number(const char*& string_pointer, char search, int& num)
     {
     Move_Pointer_NextLine(string_pointer, search);
     num = Get_Number(string_pointer);                         
     }
     
void MPNL_Get_FP_Number(const char*& string_pointer, char search, double& num)
     {
     Move_Pointer_NextLine(string_pointer, search);
     num = Get_FP_Number(string_pointer);  
     }

void MPNL_Get_Numbers(const char*& string_pointer, char search, int* nums, int count)
     {
	 Zero_Array(nums, count);
     Move_Pointer_NextLine(string_pointer, search);
     Get_Numbers(string_pointer, nums, count);
     }
     
void MPNL_Get_FP_Numbers(const char*& string_pointer, char search, double* nums, int count)
     {
	 Zero_Array(nums, count);
     Move_Pointer_NextLine(string_pointer, search);
     Get_FP_Numbers(string_pointer, nums, count);
     }
     
void MPNL_Get_Names(const char*& string_pointer, char search, int index, char** name_array, int num_names)
     {
     Move_Pointer_NextLine(string_pointer, search, index);
     int name_index = 0;
     int character_index = 0;

	 for (int a = 0; a < num_names; ++a)
	  //Initialize name array.
		 {
		 name_array[a] = const_cast<char*>(&(string_pointer[0]));
		 }
     while ( (string_pointer[character_index] != '\n') && (name_index < num_names) )
           {
           ++character_index;
           if ( Is_Text(string_pointer[character_index]) && (string_pointer[character_index - 1] == ' ') )
              {
              name_array[name_index] = const_cast<char*>(&(string_pointer[character_index]));
              ++name_index;
              }
           }
     }

void MPSNL_Get_Numbers(const char*& string_pointer, char*& assign_pointer, 
                       char search, int index, int* nums, int count)
     {
	 Zero_Array(nums, count);
     Move_Pointers_NextLine(string_pointer, assign_pointer, search, index);
     const char* temp_pointer = assign_pointer;
     if (temp_pointer[0] != '\n')
      //If the next-line move failed (due to lack of input), don't bother
      //analyzing any further.
        {
        Skip_Phrase(temp_pointer);
        Get_Numbers(temp_pointer, nums, count);           
        }
     }

void MPSNL_Get_FP_Numbers(const char*& string_pointer, char*& assign_pointer, 
                          char search, int index, double* nums, int count)
     {
	 Zero_Array(nums, count);
     Move_Pointers_NextLine(string_pointer, assign_pointer, search, index);
     const char* temp_pointer = assign_pointer;
     if (temp_pointer[0] != '\n')
        {
        Skip_Phrase(temp_pointer);
        Get_FP_Numbers(temp_pointer, nums, count);           
        }                    
     }

void Find_Next_Line_Start(const char*& string_pointer)
     {
     int index = 0;
     while (string_pointer[index] != '\n')
      //Find end of current line.
           {
           ++index;                       
           }      
     string_pointer = &(string_pointer[index + 1]);                                  
     }

void Find_Char_Next_Line(const char*& string_pointer, char search)
     {
     int index = 0;
     while (string_pointer[index] != '\n')
      //Find end of current line.
           {
           ++index;                       
           }       
     int temp_index = index;
     ++index;   
     while ( (string_pointer[index] != search) && (string_pointer[index] != '\n') )
      //Look for the search character, but do not go to the line after the next line.  
           {
           ++index;                         
           }
     if (string_pointer[index] == '\n')
      //If end of line is found again, just go back to the former end of line.
           {
           index = temp_index;                    
           }
     string_pointer = &(string_pointer[index]);  
      //Assign the new position (end of line 1 or the search character).      
     }

void Skip_Phrase(const char*& string_pointer)
     {
     int index = 1;
     while (!Is_Phrase_End(string_pointer[index]))
      //Find end of current line.
           {
           ++index;                       
           }      
     string_pointer = &(string_pointer[index + 1]);                  
     }

//Read pointer movement + command checking---

bool MPNL_CommandCheck(const char*& string_pointer, const char* command_pointer)
	 {
	 bool c_check = FALSEV;
     Find_Next_Line_Start(string_pointer); 
	 if (Command_Check(string_pointer, command_pointer))
		{
		c_check = TRUEV;
	    }
	 else
	    {
		string_pointer = &(string_pointer[-1]);  
		 //Return pointer to the previous line.
	    }
	 return c_check;
     }

bool MPNL_CommandCheck(const char*& string_pointer, const char* command_pointer, 
		               int index)
	 {
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 if (c_check)
      //Move pointer ahead in string by the passed amount.
		{
	    string_pointer = &(string_pointer[index]);  
	    }
	 return c_check;
     }

bool MPNL_CommandCheck_SC(const char*& string_pointer, const char* command_pointer)
	 {
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 if (c_check)
		{
	    Skip_Phrase(string_pointer); 
	    }
	 return c_check;
     }

bool MPSNL_CommandCheck(const char*& string_pointer, char*& second_pointer, 
		                const char* command_pointer)
	 {
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 second_pointer = const_cast<char*>(string_pointer);   
	 return c_check;
     }

bool MPSNL_CommandCheck(const char*& string_pointer, char*& second_pointer, 
		                const char* command_pointer, int index)
	 {
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 if (c_check)
		{
		second_pointer = const_cast<char*>(&(string_pointer[index]));   
	     //Second pointer gets the boost, while the first pointer remains at the
	     //beginning of the line.
		}
	 else
		{
		second_pointer = const_cast<char*>(string_pointer);
	    }
	 return c_check;
     }

bool MPSNL_CommandCheck_SC(const char*& string_pointer, char*& second_pointer, 
						   const char* command_pointer)
	 {
	 bool c_check = MPNL_CommandCheck_SC(string_pointer, command_pointer);
	 second_pointer = const_cast<char*>(string_pointer);
	 return c_check;
     }

bool MPNL_CommandCheck_GetPosNum(const char*& string_pointer, const char* command_pointer,
								 int& num)
	 {
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 num = Get_Positive_Number(string_pointer);
	 return c_check;
     }

bool MPNL_CommandCheck_GetNum(const char*& string_pointer, const char* command_pointer,
							  int& num)
	 {
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 num = Get_Number(string_pointer);
	 return c_check;
     }

bool MPNL_CommandCheck_GetFPNum(const char*& string_pointer, const char* command_pointer,
							    double& num)
	 {
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 num = Get_FP_Number(string_pointer);
	 return c_check;
     }

bool MPNL_CommandCheck_GetFPNum(const char*& string_pointer, const char* command_pointer,
							    double& num, double default_value)
	 {
	 num = default_value;
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 if (c_check)
		{
		num = Get_FP_Number(string_pointer);
	    }
	 return c_check;
     }

bool MPNL_CommandCheck_GetNums(const char*& string_pointer, const char* command_pointer,
							   int* nums, int count)
	 {
	 Zero_Array(nums, count);
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 if (c_check)
		{
		Get_Numbers(string_pointer, nums, count);
	    }
	 return c_check;
     }

bool MPNL_CommandCheck_GetFPNums(const char*& string_pointer, const char* command_pointer,
							     double* nums, int count)
	 {
	 Zero_Array(nums, count);
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 if (c_check)
		{
		Get_FP_Numbers(string_pointer, nums, count);
	    }
	 return c_check;
     }
     
bool MPNL_CommandCheck_GetNames(const char*& string_pointer, const char* command_pointer,
							    int index, char** name_array, int name_count)
	 {
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 for (int a = 0; a < name_count; ++a)
	  //Initialize.
		 {
		 name_array[a] = const_cast<char*>(&(string_pointer[0]));
		 }
	 if (c_check)
		{
	    string_pointer = &(string_pointer[index]);  
		int name_index = 0;
		int character_index = 0;
		while ( (string_pointer[character_index] != '\n') && (name_index < name_count) )
           {
           ++character_index;
           if ( Is_Text(string_pointer[character_index]) && 
			    (string_pointer[character_index - 1] == ' ') )
              {
              name_array[name_index] = const_cast<char*>(&(string_pointer[character_index]));
              ++name_index;
              }
           }
	    }
	 return c_check;
	 }

bool MPSNL_CommandCheck_GetNums(const char*& string_pointer, char*& assign_pointer,
								const char* command_pointer, int index, int* nums, int count)
	 {
     Zero_Array(nums, count);
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 if (c_check)
		{
		assign_pointer = const_cast<char*>(&(string_pointer[index]));  
		const char* temp_pointer = assign_pointer;
        Skip_Phrase(temp_pointer);
        Get_Numbers(temp_pointer, nums, count);           
	    }
	 else
	    {
	    assign_pointer = const_cast<char*>(string_pointer);  
	    }
	 return c_check;
     }

bool MPSNL_CommandCheck_GetFPNums(const char*& string_pointer, char*& assign_pointer,
							      const char* command_pointer, int index, double* nums, int count)
	 {
	 Zero_Array(nums, count);
	 bool c_check = MPNL_CommandCheck(string_pointer, command_pointer);
	 if (c_check)
		{
		assign_pointer = const_cast<char*>(&(string_pointer[index]));  
		const char* temp_pointer = assign_pointer;
        Skip_Phrase(temp_pointer);
        Get_FP_Numbers(temp_pointer, nums, count);           
	    }
	 else
	    {
	    assign_pointer = const_cast<char*>(string_pointer);  
	    }
	 return c_check;
     }


//File-to-character logic---

void Get_Characterized_File(ifstream& file_reader, char* file_text, 
                            int array_size, int& file_size)
    {
    const char TERM_CHAR = '#';
	 //Indicates the end-of-input-file character.
    file_size = 0;
    while ( !file_reader.eof() && (file_size < ( array_size - MAX_COMMAND_LENGTH ) ) 
            && (file_text[file_size - 1] != TERM_CHAR) )
          {
          file_reader.get(file_text[file_size]);
          if (file_text[file_size] == '\n')
           //Line feed character. Format by adding NULL character.
                  {
                  file_text[file_size] = '\0';
                  ++file_size;
                  file_text[file_size] = '\n';     
                  }
          ++file_size;                    
          }            
    for (int a = file_size; a < ( array_size - 3 ); ++a)
     //Add final terminating characters that indicate the end of a file.
        {               
        file_text[a] = TERM_CHAR; 
        }            
    file_text[array_size - 3] = '\0';
    file_text[array_size - 2] = '\n';
	 //Add normal c-string terminating characters to the string.
    }

void String_Replace(char* change_string, int change_text_length, 
	                int total_character_size, const char* rep_text, int rep_text_length)
    {
	int length_change = rep_text_length - change_text_length;
	int start_index = 0;

	//Shift the entire character array as in vector insertion logic---

	if (length_change < 0)
	 //Negative shift of characters.
	   {
	   for (int a = 0; a < (total_character_size + length_change); ++a)
	     {
		 change_string[a] = change_string[a - length_change];
	     }
	   change_string[total_character_size + length_change] = '\0';
	   }

    if (length_change > 0)
	 //Positive shift of characters.
	   {
	   for (int a = total_character_size + length_change - 1; a >= length_change; --a)
	     {
		 change_string[a] = change_string[a - length_change];
	     }
	   change_string[total_character_size + length_change] = '\0';
	   }

	for (int a = 0; a < rep_text_length; ++a)
	 //Final change/exchange of strings.
	    {
		change_string[a] = rep_text[a];
	    }

    }


//File logic---

void Skip_Phrases(ifstream& file_read, int num_skips)
     {
     char temp_string[MAX_PHRASE_LENGTH];
     for (int a = 0; a < num_skips; ++a)
         {
         file_read >> temp_string;    
         }                        
     }

void Load_Values(ifstream& file_read, int* values, int num_reads)    
     {
     for (int a = 0; a < num_reads; ++a)
         {
         file_read >> values[a];
         }                      
     }

void Load_Values(ifstream& file_read, double* values, int num_reads)    
     {
     for (int a = 0; a < num_reads; ++a)
         {
         file_read >> values[a];
         }                      
     }

void Write_Values(ofstream& file_write, const int* values, int num_writes)    
     {
     for (int a = 0; a < num_writes; ++a)
         {
         file_write << " " << values[a];
         }                      
     }
     
void Write_Values(ofstream& file_write, const double* values, int num_writes)
     {
     for (int a = 0; a < num_writes; ++a)
         {
         file_write << " " << values[a];
         }                             
     }
     
void Write_Intro(ofstream& file_write, const char* intro_string)
     {
     file_write << endl << intro_string << endl << endl;                       
     }


//General output statements---

void Show_Statement(const char* statement)
     {
	 cout << endl << statement;
     }

void Show_Statement(const char* first_part, int num)
	 {
	 char statement[MAX_STATEMENT_LENGTH];
	 strcpy(statement, first_part);
	 Add_Number(statement, num);
	 Show_Statement(statement);
	 }

void Show_Statement(const char* first_part, double num)
	 {
	 cout << endl << first_part << num; 
	  //Work in progress: No function for adding FP values to
	  //strings exist in my library.
	 }

void Show_Statement(const char* first_part, int num, const char* final_part)
	 {
	 char statement[MAX_STATEMENT_LENGTH];
	 strcpy(statement, first_part);
	 Add_Number(statement, num);
	 strcat(statement, final_part);
	 Show_Statement(statement);
	 }

void Show_Statement(const char* first_part, const char* second_part)
	 {
	 char statement[MAX_STATEMENT_LENGTH];
	 strcpy(statement, first_part);
	 strcat(statement, second_part);
	 Show_Statement(statement);
	 }

void Show_Statement(const char* first_part, const char* second_part, const char* third_part)
	 {
	 char statement[MAX_STATEMENT_LENGTH];
	 strcpy(statement, first_part);
	 strcat(statement, second_part);
	 strcat(statement, third_part);
	 Show_Statement(statement);
	 }

void Show_Warning(const char* warning) 
     {
     cout << endl << "WARNING: " << warning;
     }

void Show_Warning(const char* first_part, const char* second_part, const char* third_part)
	 {
	 char statement[MAX_STATEMENT_LENGTH];
	 strcpy(statement, first_part);
	 strcat(statement, second_part);
	 strcat(statement, third_part);
	 Show_Warning(statement);
	 }

void Show_Congratulations(const char* congrats)
	 {
	 cout << endl << "CONGRATULATIONS: " << congrats;
	 }

//Number output---

void Show_XYZ(const double* coors, double scale)
     {
	 cout << coors[0]*scale << " " << coors[1]*scale << " " << coors[2]*scale;
     }

void Write_XYZ(ofstream& file_writer, const double* coors, double scale)
     {
     file_writer << coors[0]*scale << " " << coors[1]*scale << " " << coors[2]*scale;
     }

void Read_XYZ(ifstream& file_reader, double* coors, double scale)
     {
	 file_reader >> coors[0] >> coors[1] >> coors[2];
	 Multiply_XYZ(coors, scale);
     }

//Number output with column justification---

void Output_Spaces(ofstream& writer, int num_spaces)
      {
	  for (int a = 0; a < num_spaces; ++a)
	      {
		  writer << ' ';
	      }
      }

void Left_Justify(ofstream& writer, double num, int prec, int total_space)
      {
	  writer.precision(prec);
	  int value_string_length = 1 + prec + String_Size(int(num));
	   //Decimal point + decimal places + pre-decimal number size.
	  writer << fixed << num;
	  Output_Spaces(writer, total_space - value_string_length);
      }

void Left_Justify(ofstream& writer, int num, int total_space)
      {
	  int value_string_length = String_Size(int(num));
	  writer << fixed << num;
	  Output_Spaces(writer, total_space - value_string_length);
      }

void Right_Justify(ofstream& writer, double num, int prec, int total_space)
      {
	  writer.precision(prec);
	  int value_string_length = 1 + prec + String_Size(int(num));
	   //Decimal point + decimal places + pre-decimal number size.
	  Output_Spaces(writer, total_space - value_string_length);
      writer << fixed << num;
      }

void Right_Justify(ofstream& writer, int num, int total_space)
      {
	  int value_string_length = String_Size(int(num));
	  Output_Spaces(writer, total_space - value_string_length);
      writer << fixed << num;
      }
