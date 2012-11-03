
#ifndef MATH_CHUCK
#define MATH_CHUCK

#include <cstdlib>
#include <cmath>

using namespace std;

#include "DWN_utilityconstants.h"

//Array initialization---

void Zero_Array(double*, int);
 //Sets all values in an array to zero.
void Zero_Array(int*, int);
 //Sets all values in an array to zero.
void False_Array(bool*, int);
 //Sets all values in an array to false (0 = FALSEV).

void Initialize_Array(double*, int, double);
 //Initializes an array to the passed value.
void Initialize_Array(int*, int, int);
 //Initializes an array to the passed value.
void Initialize_Array(bool*, int, bool);
 //Initializes an array to the passed value.
void Initialize_Array(char*, int, char);
 //Initializes an array to the passed value.

//Relative math---

int Get_Largest(int, int, int);
 //Returns the largest of the three passed integers.
double Get_Largest(double, double, double);
 //Returns the largest of the three passed values.
int Get_Smallest(int, int, int);
 //Returns the smallest of the three passed integers.
double Get_Smallest(double, double, double);
 //Returns the smallest of the three passed values.
int Get_Smallest_NonZero(int, int, int);
 //Returns the smallest non-zero of the three passed integers.
double Get_Smallest_NonZero(double, double, double);
 //Returns the smallest non-zero of the three passed values.

int Get_Largest_3DVec_Component(const double*);  
 //Returns the index of the component in the passed 3-D vector
 //that has the largest magnitude.
int Get_Largest_RelVec_Component(const double*, const double*);
 //Returns the index of the components in the passed relative
 //3-D vector that has the largest magnitude.

int Abs_Difference(int, int);
 //Returns the absolute difference between the integers.
double Abs_Difference(double, double);
 //Returns the absolute difference between the values.

//Numerical range enforcement---

void Set_Max(int&, int);
 //If passed variable is greater than the maximum value,
 //the variable is set to that maximum value.
void Set_Max(double&, double);
 //If passed variable is greater than the maximum value,
 //the variable is set to that maximum value.
void Set_Min(int&, int);
 //If passed variable is greater than the minimum value,
 //the variable is set to that maximum value.
void Set_Min(double&, double);
 //If passed variable is greater than the minimum value,
 //the variable is set to that maximum value.
void Set_Maxs(double*, const double*, int);
 //Enforces maximum values on a set of variables.
void Set_Mins(double*, const double*, int);
 //Enforces minimum values on a set of variables.
 
//Safe floating-point math--- 
 
void Check_FP_Zero(double&);
 //Makes sure a floating-point value of zero is actually stored as zero
 //in computational floating-point terms. See the "utilityconstants.h"
 //file for the FP-range that is considered effectively zero.
void Check_FP_Zeros(double*);
 //Makes sure the three passed FP values are stored as zero if they
 //effectively are zero.
bool Check_FP_Equality(double, double);
 //Determines if the passed FP variables are effectively equal to
 //each other, i.e. the difference between them is effectively zero.
bool Check_FP_LessOrEqual(double, double);
 //<= operation between FP variables. Ignores insignificant differences.
 //Here, that means effectively equal FP variables will return TRUEV.
bool Check_FP_GreaterOrEqual(double, double);
 //>= operation between FP variables. Ignores insignificant differences.
 //Here, that means effectively equal FP variables will return TRUEV.
bool Check_FP_LessThan(double, double);
 //< operation between FP variables. Ignores insignificant differences.
 //Here, that means effectively equal FP variables will return FALSEV.
bool Check_FP_GreaterThan(double, double);
 //> operation between FP variables. Ignores insignificant differences.
 //Here, that means effectively equal FP variables will return FALSEV.

void Protect_Truncation(double&);
 //Increases the magnitude of the passed variable by a very small value, such
 //that truncation will return the desired integer value.
void Safe_Truncation(int&, double, double);
 //Obtains the integer truncation of the result of an FP division.

//Peak analysis---

void Analyze_Peaks(double*, int, double, double, double**, int&, const int);
 /*Analyzes a data set and determines peak locations, maximums,
   FWHMs, and integrated intensity are determined. Assumes
   a zero-baseline for determination of FWHM and integrated intensity.*/

//Memory management---

//Dynamic two-dimensional memory handling---

void Get_Memory(double**&, int, int);
 //Dynamically allocates memory of two-dimensional FP array.
void Free_Memory(double**&, int);
 //Dynamic deallocation of memory. 
void Get_Memory(int**&, int, int);
 //Dynamically allocates memory of two-dimensional integer array.
void Free_Memory(int**&, int);
 //Dynamic deallocation of memory. 
void Get_Memory(char**&, int, int);
 //Dynamically allocates memory of two-dimensional character array.
void Free_Memory(char**&, int);
 //Dynamic deallocation of memory. 


#endif