
#ifndef ORGANIC_RACHEL
#define ORGANIC_RACHEL

#include <fstream>

using namespace std;

#include "DWN_atom.h"
#include "DWN_atomcoll.h"

const int INIT_MOLECULE_SIZE = 100;
 //Initial amount of memory reserved for atoms in the molecule.

const bool MOLECULE_CENTER = 1;
const bool FIRST_ATOM = 0;
 //Methods of molecule placement at spatial locations, 
 //either by its coordinate average or the first atom defined.

//Standard bond lengths for automatic molecule construction---

const double OH_BOND = 0.100;
const double CH_BOND = 0.109;
const double CC_BOND = 0.1544;
const double CN_BOND = 0.147;
const double NH_BOND = 0.101;
const double CO_BOND = 0.143;
const double SC_BOND = 0.181;
const double SH_BOND = 0.135;
 //Single bond.
const double CO_DBOND = 0.123;
 //Double bond.
const double R_BOND = CC_BOND;
 //For initial placement of the R-groups 
 //during peptide creation.

class molecule
      {
      public:
      
      //Constructors, destructors, and initialization functions---

      molecule();
      molecule(const molecule&);
      ~molecule();
      
      void Initialize();
	   //Initializes molecule properties.
      
	  //Information retrieval---

	  int Size() const;
	   //Returns number of atoms in the molecule.
	  double Get_Charge() const;
	   //Returns the charge of the molecule.

	  //Properties setting---

      void operator = (const molecule&);
       //Atomic collection becomes that
	   //of the passed molecule.

      //Building functions---
      
      void Add_Atom(const atom&, const double*);
	   //Adds an atom to the molecule at the passed coordinates.
      void Add_Atom(const atom&, double, double, double);
       //Adds an atom to the molecule at the passed coordinates.
      void Add_Atom(const atom&, const double*, double, int, int, bool);
	   /*Adds an atom based on geometrical variables and the location of a
	     different atom that the atom is to be bonded to.
	     e.g. this can be used for following the tetrahedral atom placement in 
	     a carbon backbone.*/
            
      void Place_Water(const double*, const atom&, const atom&);
       //Adds a water molecule at the passed oxygen atom coordinates.
       //Last two atoms are the water and hydrogen atomic information. 
       
      //Adding specific chemical group types--- 
       
      void Add_CH3(const double*, const double*, const atom&, const atom&);
      void Add_CH2(const double*, const double*, const atom&, const atom&);
      void Add_SH(const double*, const double*, const atom&, const atom&);
      void Add_OH(const double*, const double*, const atom&, const atom&);
	   //Adds a chemical group to another atom, e.g. X-CH3, by knowledge 
	   //of the location of the X atom and that of the central atom (e.g. C).
      
      //Adding general chemical group types---
      
      void Add_AX(const double*, const double*, const atom&, const atom&, 
                  double, int);
       //Adds an AX chemical unit which is bonded to a previous atom, based
	   //on bond vector and geometric shape arguments.
      void Add_AXY(const double*, const double*, const atom&, const atom&, 
                   const atom&, double, double, int);
	   //Adds an AXY chemical unit.
      void Add_AXYZ(const double*, const double*, const atom&, const atom&, 
                    const atom&, const atom&, double, double, double, int);
	   //Adds an AXYZ chemical unit.

      //Adding pieces of molecules---
      
	  void Add_Amino_Acid(const double*, const double*, double*, const atom&, 
		                  const atom&, const atom&, const atom&, const atom&, 
						  const atom&, bool);
	   //Adds an amino acid residue with the R group represented with
	   //a dummy "R" atom. Starts with carbon atom of the carbonyl group
	   //in amino acid addition.
      void Add_Carbon_Chain(const double*, const double*, double*, double*, 
                            int, const atom&, const atom&);
	   //Adds a carbon chain in terms of number of -CH2 units.
	  void Place_R_Group(const double*, const double*, char, const atom&, const atom&, 
		                 const atom&, const atom&, const atom&, const atom&);
	   //Places an amino acid substituent group at the passed position.

       //Building specific molecules---
       
      void Build_Alkanethiol(int, const atom&, const atom&, const atom&, const atom&);
	    //Constructs an alkanethiol molecule of desired length.
      void Build_Peptide(int, const atom&, const atom&, const atom&, const atom&, 
		                 const atom&, const atom&, const atom&, const atom&);
	    //Constructs a peptide chain of desired length with all substituent groups
	    //replaced with "R" atoms.
	  void Peptide_Substituent_Addition(const char*, int, const atom&, const atom&, 
		                                const atom&, const atom&, const atom&, const atom&);
	    //Replaces "R" groups in a basic peptide chain with amino acid substituent groups.

      //Utility functions---
      
      void Return_Molecule_Center(double*) const;
       //Returns the atomic coordinate average.
      void Translate_Molecule(const double*);
       //Shifts the molecule location.
      void ReAssign_Center(const double*);
       //Moves the center of the molecule to a new location.
      
	  void Rotate_Molecule(const double*, const double*);
       //Rotates the molecule by vec A --> vec B logic.
	  void Rotate_Molecule(const double*, const double*, const double*);
       //Rotates the molecule by vec A --> vec B logic with respect
	   //to an origin.
	  void Rotate_Molecule(const double*, double);
	   //Rotates the molecule by a rotation vector and angle.
      void Rotate_Molecule(const double*, double, const double*);
       //Rotates the molecule by a rotation vector and angle with
	   //respect to an origin.
	  void Rotate_MoleculeXYZ(const double*, bool);
	   //Rotate the molecule by an X-Y-Z rotation, clockwise (TRUEV)
	   //or counter-clockwise (FALSEV).
	  void Rotate_MoleculeXYZ(const double*, bool, const double*);
	   //Rotates the molecule by an X-Y-Z rotation with respect
	   //to an origin, clockwise (TRUEV) or counter-clockwise (FALSEV).

	  void Orient_Molecule(const double*);
	   //Re-orients the molecule such that its backbone vector
	   //(last atom location - first atom location) matches the passed vector.
	  void Orient_Molecule(const double*, const double*);
	   //Re-orients the backbone vector of the molecule and then reorients
	   //its first bond (second atom location - first atom location) via
	   //rotation about the backbone vector.
	  void Orient_Molecule(double);
	   //Rotates a molecule about the backbone vector by the passed angle.
	  void Standard_Orientation();
	   //Orients the molecule by standard orientation, that being:
	   //(1) backbone along the -z axis, and (2) first bond along the +x direction
	   //as much as is possible without molecular distortion.
      
      //void Get_Formula(char*);
       //Returns the full chemical formula of the molecule.
     // double Calculate_IVolume_Occupied();
       //Returns effective initial volume occupied by the molecule.
      
	  //Molecule copying---
      
	  int Copy_Molecule(atom_collection&) const;
       //Adds a copy of the molecule into an atomic collection.
      int Copy_Molecule(atom_collection&, const double*, bool) const;
       //Makes a copy of the molecule into an atomic collection, placing
	   //the molecule by its center or its first atom.
      int Copy_Molecule(atom_collection&, const double*,
                        const double*, bool) const;
       //Makes a copy of the molecule into an atomic collection, placing 
	   //the molecule by its center or its first atom. Also performs X-Y-Z rotation.
	  int Copy_Molecule(atom_collection&, const double*,
                         const double*, const double*, bool) const;
       //Makes a copy of the molecule into an atomic collection, placing 
	   //the molecule by its center or its first atom. Also performs 
	   //start vec-->end vec rotation.

      //I/O (units: Angstroms, amu)---
      
      void Store_Atom_Locations(ofstream&) const;
       //Stores atomic information in a file.
      
      void Save_Molecule(const char*) const;
	   //Saves a molecule to file.
      void Load_Molecule(const char*);
       //Loads a molecule from file.
      
      private:
      
	  void Evaluate_Addition(double*, int, double*, int);
	   //Assists chain-building functions.

      atom_collection molecule_atoms;
       //Atoms of the molecule.
      };


#endif
