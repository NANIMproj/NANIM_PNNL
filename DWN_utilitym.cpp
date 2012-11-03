
#include "stdafx.h"

#include "DWN_utilitym.h"
 

//Array initialization---

void Zero_Array(double* arr, int num_elements)
	{
	Initialize_Array(arr, num_elements, 0.0);
	}

void Zero_Array(int* arr, int num_elements)
	{
	Initialize_Array(arr, num_elements, 0);
	}

void False_Array(bool* arr, int num_elements)
	{
	Initialize_Array(arr, num_elements, FALSEV);
	}

void Initialize_Array(double* arr, int num_elements, double init_value)
	{
	for (int a = 0; a < num_elements; ++a)
		{
		arr[a] = init_value;
		}
	}

void Initialize_Array(int* arr, int num_elements, int init_value)
	{
	for (int a = 0; a < num_elements; ++a)
		{
		arr[a] = init_value;
		}
	}

void Initialize_Array(bool* arr, int num_elements, bool init_value)
	{
	for (int a = 0; a < num_elements; ++a)
		{
		arr[a] = init_value;
		}
	}

void Initialize_Array(char* arr, int num_elements, char init_value)
	{
	for (int a = 0; a < num_elements; ++a)
		{
		arr[a] = init_value;
		}
	}

//Relative math---
     
int Get_Largest(int a, int b, int c)
     {
     int largest = a;
     if (b > a)
        {
        largest = b;     
        }
     if (c > largest)
        {
        largest = c;        
        }
     return largest;
       }

double Get_Largest(double a, double b, double c)
     {
     double largest = a;
     if (b > a)
        {
        largest = b;     
        }
     if (c > largest)
        {
        largest = c;        
        }
     return largest;
     }

int Get_Smallest(int a, int b, int c)
     {
     int smallest = a;
     if (b < a)
        {
        smallest = b;     
        }
     if (c < smallest)
        {
        smallest = c;        
        }
     return smallest;
     }  
       
double Get_Smallest(double a, double b, double c)
     {
     double smallest = a;
     if (b < a)
        {
        smallest = b;     
        }
     if (c < smallest)
        {
        smallest = c;        
        }
     return smallest;
     }

int Get_Smallest_NonZero(int a, int b, int c)
     { 
     int smallest = a;
     if ( (b < a) || (a == 0) )
        {
        smallest = b;     
        }
     if ( ( c < smallest) || (smallest == 0) )
        {
        smallest = c;        
        }
     return smallest;                    
     }
     
double Get_Smallest_NonZero(double a, double b, double c)
     {
     double smallest = a;
     if ( (b < a) || (Check_FP_Equality(a, 0.0)) )
        {
        smallest = b;     
        }
     if ( (c < smallest) || (Check_FP_Equality(smallest, 0.0)) )
        {
        smallest = c;        
        }
     return smallest;                               
     }

int Get_Largest_3DVec_Component(const double* vec)
     {
	 int comp = 0;
	 double avec[3];
	 Set_XYZ(avec, vec);
	 Abs_XYZ(avec);
	 double largest = avec[0];
	 if (avec[1] > largest)
		{
		comp = 1;
		largest = avec[1];
		}
	 if (avec[2] > largest)
	    {
		comp = 2;
		}
	 return comp;
     }

int Get_Largest_RelVec_Component(const double* vec, const double* rel_vec)
	 {
	 int comp = 0;
	 double avec[3];
	 avec[0] = abs(vec[0]/rel_vec[0]);
	 avec[1] = abs(vec[1]/rel_vec[1]);
	 avec[2] = abs(vec[2]/rel_vec[2]);
	 double largest = avec[0];
	 if (avec[1] > largest)
		{
		comp = 1;
		largest = avec[1];
		}
	 if (avec[2] > largest)
	    {
		comp = 2;
		}
	 return comp;
	 }

int Abs_Difference(int val1, int val2)
    {
    int dif = val1 - val2;
    int abs_dif = abs(dif);
    return abs_dif;                    
    }
    
double Abs_Difference(double val1, double val2)
    {
    double dif = val1 - val2;
    double abs_dif = abs(dif);
    return abs_dif;                         
    }

//Numerical range enforcement---

void Set_Max(int& vari, int max)
    {
    if (vari > max)
       {
       vari = max;      
       }                 
    }
    
void Set_Max(double& vari, double max)
    {
    if (vari > max)
       {
       vari = max;      
       }                 
    }

void Set_Min(int& vari, int min)
    {
    if (vari < min)
       {
       vari = min;      
       }                 
    }
     
void Set_Min(double& vari, double min)
    {
    if (vari < min)
       {
       vari = min;      
       }                 
    }
     
void Set_Maxs(double* varis, const double* max, int vari_count)
    {
    for (int a = 0; a < vari_count; ++a)
        {
        Set_Max(varis[a], max[a]);     
        }                  
    }
     
void Set_Mins(double* varis, const double* min, int vari_count)
    {
    for (int a = 0; a < vari_count; ++a)
        {
        Set_Min(varis[a], min[a]);
        }              
    }

//Floating-point math logic---

void Check_FP_Zero(double& variable)
    {
    if ( (variable < ZERO_UPP_BOUND) && (variable > ZERO_LOW_BOUND) )
     //Set values that are very close to zero to zero.
      {
      variable = 0.0;                          
      }                        
    }

void Check_FP_Zeros(double* variables)
    {
    for (int a = 0; a < 3; ++a)
        {
		Check_FP_Zero(variables[a]);
        }
    }   
     
bool Check_FP_Equality(double variable1, double variable2)
    {
    bool equal_check = FALSEV;
    double dif = variable1 - variable2;
    Check_FP_Zero(dif);
    if (dif == 0.0)
       {
       equal_check = TRUEV;           
       }                          
    return equal_check;
    }

bool Check_FP_LessOrEqual(double variable1, double variable2)
    {
    bool equal_check = FALSEV;
    double dif = variable1 - variable2;
    Check_FP_Zero(dif);
    if (dif <= 0.0)
       {
       equal_check = TRUEV;           
       }                          
    return equal_check;
    }
     
bool Check_FP_GreaterOrEqual(double variable1, double variable2)
    {
    bool equal_check = FALSEV;
    double dif = variable1 - variable2;
    Check_FP_Zero(dif);
    if (dif >= 0.0)
       {
       equal_check = TRUEV;           
       }                          
    return equal_check;
    }
     
bool Check_FP_LessThan(double variable1, double variable2)
    {
    bool equal_check = FALSEV;
    double dif = variable1 - variable2;
    Check_FP_Zero(dif);
	if (dif == 0.0)
	   {
	   equal_check = FALSEV;
	   }
    else if (dif < 0.0)
       {
       equal_check = TRUEV;           
       }                          
    return equal_check;
    }
     
bool Check_FP_GreaterThan(double variable1, double variable2)
    {
    bool equal_check = FALSEV;
    double dif = variable1 - variable2;
    Check_FP_Zero(dif);
	if (dif == 0.0)
	   {
	   equal_check = FALSEV;
	   }
    else if (dif > 0.0)
       {
       equal_check = TRUEV;           
       }                          
    return equal_check;
    }

void Protect_Truncation(double& variable)
	{
	if (variable > 0.00) 
	    {
        variable += FP_ERROR_FIX;
	    }
	 else	  
 	    {
        variable -= FP_ERROR_FIX;
	    }
	}

void Safe_Truncation(int& division_result, double numer, double denom)
	{
	Protect_Truncation(numer);
	division_result = int(numer/denom);
	}

//PEAK ANALYSIS---

void Analyze_Peaks(double* data, int num_points, double min_x,
	               double precision, double** analysis_res, 
				   int& num_peaks, const int MAX_PEAKS)
     {
	 num_peaks = 0;
	 for (int a = 0; a < MAX_PEAKS; ++a)
	  //Initialize peak array.
	     {
		 Zero_Array(analysis_res[a], 4);
	     }

	 //Find the first occurrence of a positive slope to begin analysis---

	 double new_min_x;
	 int new_start = 0;
	 double x_val = min_x;
	 for (int a = 0; a < (num_points - 1); ++a)
	     {
		 if (data[a + 1] > data[a])
		    {
			new_min_x = x_val;
			new_start = a;
			break;
		    }
		 x_val += precision;
	     }
	  double x_range = double(num_points - new_start) * precision; 

	 //Characterize the lowest, highest, and average values
	 //in the data set---

	 double lowest_value = data[new_start];
	 double highest_value = data[new_start];
	 double average = data[new_start];
	 for (int a = new_start + 1; a < num_points; ++a)
	     {
		 if (data[a] > highest_value)
		    {
			highest_value = data[a];
		    }
		 if (data[a] < lowest_value)
		    {
			lowest_value = data[a];
		    }
		 average += data[a]/double(num_points - new_start);
	     }

	 double value_range = highest_value - lowest_value;
	 double amp_criterion = 0.05*value_range;
	  //Peaks are determined as locations of maximum values
	  //where the amplitude decays to at least 5% of the y-range
	  //of values on both the left and right sides of the peak.
	 bool left_test, right_test;
	  //Booleans used in that test.
	 double required_value;
	  //The amplitude that a peak must decay to on both sides in
	  //order to truly be considered a peak of interest.

	 //Find and evalute the peaks---

	 double target_intensity, left_half_maximum, right_half_maximum;
	 bool found_left_half_maximum, found_right_half_maximum;
	 double right_contribution, left_contribution;
	 int steps_to_half_maximum;
	  //Variables used in the analysis of the full-width at half
	  //maximum (FWHM) of the peaks.
	 for (int a = new_start; a < num_points; ++a)
	  //Test each data point as a potential maximum.
	     {
		 right_test = FALSEV;
		 left_test = FALSEV;
		 required_value = data[a] - amp_criterion;
		 for (int b = (a + 1); b < num_points; ++b)
		  //Right side.
		     {
			 if (data[b] < required_value)
			    {
			    b = num_points;
				right_test = TRUEV;
			    }
			 else if (data[b] > data[a])
			    {
				b = num_points;
			    }
		     }
		 if (right_test)
		     {
		     for (int c = (a - 1); c >= new_start; --c)
		      //Left side.
		        {
			    if (data[c] < required_value)
			       {
			       c = new_start - 1;
				   left_test = TRUEV;
			       }
			    else if (data[c] > data[a])
			       {
			 	   c = new_start - 1;
			       }
			    }
		     }
		 if (left_test && right_test)
		  //Found a peak!
		     {
			 analysis_res[num_peaks][0] = new_min_x + double(a - new_start) * precision;
			 analysis_res[num_peaks][1] = data[a];
			  //Record peak location and intensity at peak maximum.

			 if (analysis_res[num_peaks][1] < 0.0)
			  //Peaks located at negative intensities can not be fully analyzed,
			  //so stop analysis here.
			    {
			    ++num_peaks;
			    continue;
			    }

			 target_intensity = data[a]/2.0;
			  //Assumes the peak arises out of a zero-intensity baseline.
			  //Thus, the FWHM occurs at half the maximum intensity.
			 right_contribution = data[a]*precision;
			 found_right_half_maximum = FALSEV;
			 for (int b = a + 1; b < num_points; ++b)
			     {
				 right_contribution += data[b]*precision;
			     if (!found_right_half_maximum && (data[b] < target_intensity))
				    {
				    found_right_half_maximum = TRUEV;
					steps_to_half_maximum = b - a;
					right_half_maximum = analysis_res[num_peaks][0] 
					                   + double(steps_to_half_maximum - 1)*precision;
					 //Set right half maximum to the point before crossing the half-maximum
					 //intensity.
					right_half_maximum += precision * ( target_intensity - data[b - 1] ) 
						                  / ( data[b] - data[b - 1] );
				    }
				 else if (found_right_half_maximum && ( (b - a) > 2*steps_to_half_maximum))
				  //Limit for integrated intensity calculation: 2 * steps to half maximum
				  //of peak from peak center.
				    {
					b = num_points;
				    }
				 else if (data[b] > data[a])
				   //FWHM can not be found, due to a neighboring, intense peak. Abort mission.
				    {
					b = num_points;
				    }
			     }

			 left_contribution = 0.0;
			  //Don't want to double-count peak maximum's contribution to integrated intensity.
			 found_left_half_maximum = FALSEV;
			 for (int c = a - 1; c >= new_start; --c)
			     {
				 left_contribution += data[c]*precision;
			     if (!found_left_half_maximum && (data[c] < target_intensity))
				    {
				    found_left_half_maximum = TRUEV;
					steps_to_half_maximum = a - c;
					left_half_maximum = analysis_res[num_peaks][0] 
					                   - (steps_to_half_maximum - 1)*precision;
					left_half_maximum -= precision * ( target_intensity - data[c + 1] ) 
						                  / ( data[c] - data[c + 1] );
				    }
				 else if (found_left_half_maximum && ( (a - c) > 2*steps_to_half_maximum))
				  //Limit for integrated intensity calculation: 2 * steps to half maximum
				  //of peak from peak center.
				    {
					c = new_start - 1;
				    }
				 else if (data[c] > data[a])
			      //FWHM not discovered, due to neighboring, intense peak.
				    {
					c = new_start - 1;
				    }
			     }
				 
			 if (found_left_half_maximum && found_right_half_maximum)
			     {
				 analysis_res[num_peaks][2] = left_contribution + right_contribution; 
				 analysis_res[num_peaks][3] = right_half_maximum - left_half_maximum;
			     }

			 ++num_peaks;
			 if (num_peaks == MAX_PEAKS)
			    {
				Show_Warning("TOO MANY PEAKS TO PROCESS IN DATA ANALYSIS!");
			    a = num_points;
				}
		     }
	     }
     }
    
//MEMORY---

void Get_Memory(double**& ptr, int a, int b)
    {
    ptr = new double*[a];
    for (int index = 0 ; index < a; ++index)
       {
       ptr[index] = new double[b];     
       }
    }

void Free_Memory(double**& ptr, int a)
    {
    for (int index = 0; index < a; ++index)
       {
       delete[] ptr[index];      
       }
    delete[] ptr;
    }

void Get_Memory(int**& ptr, int a, int b)
    {
    ptr = new int*[a];
    for (int index = 0 ; index < a; ++index)
       {
       ptr[index] = new int[b];     
       }
    }

void Free_Memory(int**& ptr, int a)
    {
    for (int index = 0; index < a; ++index)
       {
       delete[] ptr[index];      
       }
    delete[] ptr;
    }

void Get_Memory(char**& ptr, int a, int b)
    {
    ptr = new char*[a];
    for (int index = 0 ; index < a; ++index)
       {
       ptr[index] = new char[b];     
       }
    }

void Free_Memory(char**& ptr, int a)
    {
    for (int index = 0; index < a; ++index)
       {
       delete[] ptr[index];      
       }
    delete[] ptr;
    }
