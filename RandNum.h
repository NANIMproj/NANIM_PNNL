//Random number generator

#ifndef RANDVAL_XXX
#define RANDVAL_XXX

#include <cstdlib>
#include <time.h>

using namespace std;

//Random number seeding---

void Seed_RandNums();
 //Seeds the random number generator.

//Random number retrieval---

short int GetRandPosNum(short int);
 //Obtains a positive number between zero and the passed integer.
short int GetRandInteger(short int, short int);
 //Obtains an integer within the passed range.
long int GetRandPosNum(long int);
 //Obtains an integer between zero and the passed integer. Max = 1,073,741,823.
long int GetRandInteger(long int, long int);
 //Obtains an integer within the passed range.
double GetRandPosDec(double);
 //Gets a random floating point number between 0.00 and the passed floating point value.
double GetRandDec(double, double);
 //Gets a random floating point number within the passed range.

//Random number sets---

void Get_Rand_Direction(double*, int);
 //Gets a random vector (values of each coordinate vary from -1.0 to 1.0)
 //of n-dimensionality. Note that the vector magnitude is not set to a constant.
void Get_Rand_Box_Position(double*, const double*, const double*);
 //Gets a random position within the box identified by its corner and size.

#endif