//Random number generator

#include "stdafx.h"

#include "RandNum.h"

//Random number seeding---

void Seed_RandNums()
 //Seeds the random number generator
{
//timeval ze_time;
//gettimeofday(&ze_time, NULL);
//srand(ze_time.tv_usec * ze_time.tv_sec);
 //Use microseconds in seeding random numbers.
srand((unsigned)time(NULL));     
}

//Random number retrieval---

short int GetRandPosNum(short int max)
 //Obtains a positive number between zero and the passed integer.
{
return (rand() % (max + 1));          
}

short int GetRandInteger(short int min, short int max)
 //Obtains an integer within a certain range.
{
return ((rand() % (max + 1 - min)) + min);           
}

long int GetRandPosNum(long int max)
{
long int randA = long int(rand()) % 32768;
long int randB = long int(rand()) % 32768;
long int rand_val = randA * 32768 + randB;
rand_val = rand_val % (max + 1);
return rand_val;
}

long int GetRandInteger(long int min, long int max)
{
long int randA = long int(rand()) % 32768;
long int randB = long int(rand()) % 32768;
long int rand_val = randA * 32768 + randB;
rand_val = ( rand_val % (max + 1 - min) ) + min;
return rand_val;
}

double GetRandPosDec(double maximum)
 //Gets a random floating point number between 0.00 and passed float
{
return ( (double(rand()) / double(RAND_MAX) ) * maximum);      
}

double GetRandDec(double minimum, double maximum)
 //Same as above, but for obtaining a number between min and max passed.
{
return ( (double(rand()) / double(RAND_MAX)) * (maximum-minimum) + minimum);         
}

//Random number sets---

void Get_Rand_Direction(double* direc, int dimensions)
	{
	for (int a = 0; a < dimensions; ++a)
		{
		direc[a] = GetRandDec(-1.0, 1.0);
		}

	}

void Get_Rand_Box_Position(double* pos, const double* box_loc, const double* box_size)
	{
	for (int a = 0; a < 3; ++a)
	      {
		  pos[a] = box_loc[a] + GetRandDec(0.0, 1.0) * box_size[a];
		  }
	}



