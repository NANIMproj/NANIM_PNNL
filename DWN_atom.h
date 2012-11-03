
#ifndef ATOMIC_PING
#define ATOMIC_PING

#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "DWN_utility.h"

///Atom structure-------

const int MAX_ANAME_SIZE = 10;
 //Maximum number of characters in the atom's name.
 
const char DUMMY_NAME[5] = "XXX"; 
 //Invisible atom label.
const char RGROUP_NAME[3] = "R";
 //Label for an organic group represented by a single atom, e.g. amino acid side chains/substituents.
 //Used in organic molecule construction.

const double RADII_SCALING_FACTOR = 1.201;
 /*Indicates the distance cutoff between covalent and non-covalent interaction.
   E.g. If a radii sum (A + B) is 2 A and the scaling factor is 1.2, any distance
   less than 2.4 A will be considered covalent. This cut-off is the same 
   used in spatial overlap analysis to determine if atoms are "too close for 
   comfort" in simulation (excluding intended molecular bonds).*/

const int MAX_ATOMIC_NUMBER = 118;
 //Maximum atomic number that can be assigned to an atom.

class atom
       {
         
       friend class atom_collection; 
	    //Atom collection class contains a variety of functions for how
	    //atoms interact with each other. Due to this, the "friend"
	    //exception is made here.
           
       public:
       
       //Constructors and destructors---

       atom();
	    //Basic constructor. Sets atom to DUMMY ATOM, which has the name "NOATOM."
       atom(const char*, int, double, double, double);
        //Set atom name, -ic number, -ic relative mass, -ic charge, and -ic radius,
	    //in that order.
       atom(const char*, const double*, int, double, double, double);
        //Set atom name, location, -ic number, -ic relative mass, -ic charge, and -ic radius,
	    //in that order.
	   atom (const char*, const double*, int, double, double, double, int);
	    //Set atom name, location, -ic number, -ic relative mass, -ic charge, -ic radius,
	    //and group tag, in that order.
	   atom(const atom&);
	    //Sets all properties of atom to that of the passed atom.
       
	   ~atom();      
	    //Destructor.

       //Properties settings (non-spatial)---

       void Set_Atom_Specifics(const char*, const double*);
	    //Set atom name and location.
	   void Set_Atom_Properties(int, double, double, double);
        //Set atomic number, relative mass, charge, and radius.

       void Set_Atom_Name(const char*);
	    //Sets the atom name.
	   void Set_Atom_Name(const char*, bool);
	    //Sets atom name to the passed string, but stops string copying
	    //at the first capital letter beyond the first string character
	    //if boolean is set to TRUEV. (i.e. "MG" --> "M").
	   void Index_Atom_Name(int);
	    //Adds a number to the end of the atom name.
	    //If there is a number there before this function is called,
	    //that number will be removed to make way.

	   void Set_Atom_Number(int);
	    //Sets the atomic number.
	   void Set_Atom_Rel_Mass(double);
	    //Sets the relative atomic mass of the element.
	   void Set_Atom_Charge(double);
	    //Sets the atomic charge.
	   void Set_Atom_Radius(double);
	    //Sets the effective (covalent) radius.
	   void Set_Atom_Group_Tag(int);
	    //Sets the group tag (molecule bonding index).
	   
	   //Properties setting (spatial)---

       void Set_Atom_Location(double, double, double);
       void Set_Atom_Location(const double*);
        //Set atom location (3-D array).
	   void AddTo_Atom_Location(double, double, double);
       void AddTo_Atom_Location(const double*);
        //Adds the passed vector to the atom location.
       void SubFrom_Atom_Location(double, double, double);
       void SubFrom_Atom_Location(const double*);
        //Subtracts the passed vector from the atom location.
        
       //Information retriveal-- 
        
       void Get_Atom_Name(char*) const;
	   bool Same_Name(const atom&) const;
        //Returns TRUEV if this atom has the same name as the passed atom.
       int Get_Atom_Number() const;
       double Get_Atom_Rel_Mass() const;
       double Get_Atom_Charge() const;
	   double Get_Atom_Radius() const;
       void Get_Atom_Location(double&, double&, double&) const;
       void Get_Atom_Location(double*) const;
	   int Get_Atom_Group_Tag() const;
	    //Returns various properties of the atom in read-only form.

	   //Distance logic---

       double Get_Distance_To_Origin() const;
	    //Returns the distance of the atom location from the (0, 0, 0) position.
       double Get_Distance(const atom&) const;
        //Returns the distance between this atom and the passed atom.
       double Get_Distance(const atom&, int, const double*) const;
        //Returns the periodic boundary-based distance.
        //(Second parameter = periodicity, third parameter = orthorhombic box sizes).
	   double Get_Distance(const double*) const;
	    //Returns the distance between this atom and the passed coordinates.
	   double Get_Distance(const double*, int, const double*) const;
	    //Returns the PB-based distance.
       bool Covalent_Overlap(const atom&) const;
	    //Returns TRUEV if the two atoms are covalently bonded.
	    //Test = 1.201 * covalent radii sum.
	   bool Covalent_Overlap(const atom&, int, const double*) const;
	    //Returns the PB-based bond boolean.
	   bool Same_Molecule(const atom&) const;
	    //Returns TRUEV if the two atoms have the same group tag.

       //Specific identity functions---
       
       bool Is_Hydrogen() const;
        //Returns TRUEV if the atom's atomic number is 1.
       bool Is_Dummy_Atom() const;
        //Returns TRUEV if the atom is assigned the dummy atom identity shown above.
       bool Is_Invisible(bool) const;
        //Returns TRUEV if the atom should not be included in output.
        //Boolean parameter indicates if hydrogen inclusion is turned on.
	    //Tests thus for dummy atoms and hydrogen inclusion.
       
       //Operator overloads---
       
       void operator = (const atom&);
	    //Atom takes on all properties of another atom.
       bool operator == (const atom&);
        //Note: Equality comparison is done in terms of spatial location only.
       
       //File I/O (units: Angstroms, amu, elementary charge)---
       
       void Print_Atom_Location() const; 
        //Prints atom name and location on a DOS prompt (in Angstroms).
       void Store_Atom_Location(ofstream&) const;
        //Stores atomic coordinates in a file (in Angstroms).
       void Store_Atom_Location_With_Charge(ofstream&) const; 
        //Stores atomic coordinates with atomic charge as well 
        //(in Angstroms/elementary charge units).
	   void Store_Atom_Location_With_CoreShell(ofstream&) const;
	    //Writes a line for atom core and atom shells (2 lines per atom).
	    //Does not include charge.
       void Store_Atom_Info(ofstream&) const;
	    //Saves all properties of the atom to an output file.
       void Load_Atom_Info(ifstream&);
        //Reads all properties of the atom from an input file.
       
       private:
       
       //Specific atomic information-- 
             
       char atom_name[MAX_ANAME_SIZE];
        //Atomic name, i.e. "Au." This can be indexed, e.g. "Au23," if needed.
       double atom_location[3];
        //Spatial location of the atom.   
        
       //General atomic information that generally belongs to all atom of a certain type--
       
       int atomic_number;
        //Number of protons.
       double rel_atomic_mass; 
        //Relative atomic mass. 
	    //IUPAC defines this as ratio to carbon-12 mass multipled by 12.
       double atomic_charge;
        //Charge of the atom.
	   double atomic_radius;
	    //Effective radius of the atom.


	   //Molecule/bonding matrix index---

	   int group_tag;
	    //Index which represents the molecule that the atom belongs to.

       };

#endif