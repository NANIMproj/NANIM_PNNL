
#ifndef FINAL_JOE
#define FINAL_JOE

#include <iostream>
#include <cstring>
#include <cstdio>

using namespace std;

#include "DWN_atom.h"
#include "DWN_molecule.h"
#include "DWN_nanos.h"
#include "DWN_nanop.h"
#include "DWN_MD.h"
#include "DWN_scope.h"
#include "DWN_File.h"
#include "DWN_utility.h"

const int MAX_SCRIPT_FILE_SIZE = 1000000;
 //Maximum number of characters in user-defined text input file.
const int MAX_AUTO_SIZE = 100000;
 //Maximum number of characters in user-defined automation input file.

void Start();
 //Starts the program!
 
void Script_Reader(const char*, int, const char*, int, const char*, int, int, int);
 //Script engine. Reads and processes program request files.
 
//Serial reading functions---

void MoleculeAdd_Reader(const char*, double*, double*, char**); 
 //Aids in reading molecule construction information.
void Convert_ShapeFactor(char, int&);
 //Aids in reading molecular shape-factors.
void Set_Faceting_Planes(crystal_system&, const char*, double);
 //Adds a set of faceting planes to a crystal system
 //e.g. cubo-octahedral facets.

//Second order automation programming---

void Process_NAutomation_Commands(char*, int, const char*, int, int);
 //Alters the main script file text via an automation command file.
void Process_SAutomation_Commands(char*, int, char*, int, const char*, int, int);
 //Third-order level of automation. Alters the script file AND the automation command
 //file text via another automation command file.


//Analysis functions---

void Analysis_Out(const char*, const char*, int, bool, char*, bool, double, bool, double, double,
				  bool, double, double, bool, double, double, bool, double, bool, double*,
				  bool, bool, bool, double, bool, double);
 //Adds system information to the analysis file.
 //Information may include the following: Atom count, Composition,
 //Energy, AB distance average, AA distance average, General distance average,
 //Sphere fit, Box fits, Strangeness parameter, and Electrostatic potential.

#endif