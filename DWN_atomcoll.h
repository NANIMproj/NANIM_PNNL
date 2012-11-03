
#ifndef GLUE_PATRICIA
#define GLUE_PATRICIA

#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

#include "DWN_atom.h"
#include "DWN_utilitys.h"
#include "DWN_utilityconstants.h"

const int DEFAULT_INITIAL_MEMORY = 100;
 //Default initial memory allocation for the array of 
 //atoms in the collection. 

const int MAX_ATOM_TYPES = 1000;
 //Maximum number of unique types of atoms (e.g. Au, Au2, Au3, C12 isotope, C13 isotope, O, etc.)
 //expected for one atom collection. Used in creation of unique-atom arrays.

class atom_collection
     {
     public:

	 //Constructors and destructors---
     
     atom_collection();
     atom_collection(int);
     ~atom_collection();
     
	 //Memory management---

     void Set_Initial_Allocation(int);
     
     //Size and atom member retrieval---
     
     int Size() const;
      //Returns atom count.
     int Real_Size(bool) const;
	  //Returns the number of atoms that are not dummy atoms
	  //or hydrogen atoms (if boolean = FALSEV, for hydrogen inclusion).

     atom& operator [] (int);
      //Returns the atom of interest. Allows overwriting.
     atom operator [] (int) const;
      //Constant return operator. Does not perform memory check
	  //and thus allows quick read-only access.
     const atom& Get_Atom(int) const;
      //Alternative access method for read-only purposes.

     //Copying and addition of atoms---
     
     void operator = (const atom_collection&);
      //Copy operator. Takes the capacity of the passed atom collectionn
	  //as its new size and copies the properties of all atoms.

	 void Add_Atom(const atom&);
	  //Adds an atom to the end of the atom collection array.
	 void Add_Atom(const atom&, const double*, bool);
	  //Adds an atom to the end of the array. Also uses the passed
	  //spatial coordinates to change the location of the atom.
	  //TRUEV = shift original coordinates, FALSEV = place at coordinates.

	 //Liberal use of this class may not include automatic change 
	 //to the atom count value. Use these functions if you must---

     atom_collection& operator ++ ();
      //Increases the atom count of the array by one.
     void Inc_Coll_Size(int);
      //Increases the atom count of the array by the passed parameter.

	 //Group tag handling---

	 int Get_Group_Tag(int) const;
	  //Returns the group tag of the atom at the passed index.

	 void Set_Group_Tag(int, int, int);
	  //Sets the group tag of atoms in the passed index range to
	  //the passed group tag.
	 int Next_Tag();
	  //Determines the highest integer value of group tags in the 
	  //collection. Returns that integer plus one.
	 int Number_Of_Groups();
	  //Determines the number of groups in the collection.
	  //Note: Not necessarily related to the maximum.
	 int Molecule_Size(int);
	  //Returns the size of a given group (i.e. the number of atoms
	  //with a certain group tag).

	 void Assign_Group_Tags();
	  //Uses covalent overlap analysis (with atomic radii) to determine
	  //the collection's bonding network. Then proceeds to assign all
	  //atom group tag values according to that bonding network.
	 void Assign_Group_Tags(int*);
	  //Copies an array of group tag values into the group tag values
	  //of the atoms in the collection.
	 void Defrag_Group_Tags();
	  //Removes "holes" in group tag values caused by atom deletion.
	  //E.g. Molecules 0 1 3 5 8 ---> Molecules 0 1 2 3 4.

	 //Atom removal---

     void Delete_Atom(int);
	  //Removes the atom at the passed index.
	 void Delete_Atoms(int, int);
	  //Deletes atoms in the passed index range.
	 void Delete_Atoms(const double*, const double*);
	  //Deletes atoms in a given spatial region (i.e. orthorhombic box).
	  //Currently does not wrap around PBC boundaries.
	 void Make_Orthorhombic_Void(const double*, const double*);
	  //Same function as above Delete_Atoms(...) function. Nicer name.
	 void Make_Spherical_Void(const double*, double, bool);
	  //Makes a spherical cut in the atom collection as defined by
	  //the passed origin and radius. Boolean parameters is for
	  //the observation of PBC boundaries in void-creation logic.

	 //Molecule removal---

	 void Delete_Molecule(int);
	  //Deletes all atoms with the passed group tag.
	 void Delete_Molecule_Containing_Atom(int);
	  //Deletes all atoms that have the group tag belonging to
	  //atom at the passed index.

	 //Periodic box handling---

	 void Set_Periodicity(int);
	  //Sets the periodicity of the box which holds the atoms.
	 void No_Periodicity();
	  //Sets periodicity to zero. Alternative to previous function.
	 void Set_Box_Vectors(const double*);
	  //Sets the orthorhombic box parameters. Must be a 3-D array passed.
	 void Set_Box(int, const double*);
	  //Simultaneous periodicity/orthorhombic box parameter setting.

	 int Get_Periodicity() const;
	  //Returns the periodicity of the collection.
	 void Get_Box_Vectors(double*) const;
	  //Returns the orthorhombic periodic box dimensions.

	 //Spatial information---

	 void Get_Bond_Vec(int, int, double*) const;
	   //Returns a bond vector of: (b - a). Follows PBC boundaries.

	 void Get_Coordinate_Mins(double*) const;
	  //Returns the lowest points in coordinate space found (may belong to
	  //different atoms).
	 void Get_Coordinate_Mins(int, int, double*) const;
	  //Returns the lowest points in coordinate space found for atoms
	  //in a given index range.

	 void Get_Coordinate_Maxs(double*) const;
	  //Returns the highest points in coordinate space found (may belong to
	  //different atoms).
     void Coordinate_Average(double*) const;
       //Returns the average over all atomic coordinates. No PBC boundaries involved
	   //as that would not make sense here.
	 void Coordinate_Average(int, double*) const;
	   //Does the above for a certain number of atoms at the end of the collection array.
	 void Coordinate_Average(int, int, double*) const;
	   //For an (indexed) range of atoms.
	 void Fit_Sphere(double, bool, int, double&, int, bool) const;
	  /*Fits a sphere to the surface of the atomic collection. Takes the square-root
	    of the average square-distance of surface atoms from the coordinate average
	    of the system. Analysis in terms of a specific atom type (identified by the passed
	    atomic number) can be performed, but coordinate averaging is not influenced by
	    this parameter. Use a value of -1 to not use a specific atom. Also note that
		this is meant for a contained group of atoms, not crossing PBC boundaries. */
	 void Fit_Box(double, bool, int, double*, int, bool) const;
	  //Fits an orthorhombic box to the atomic collection. Although this takes many
	  //parameters, this function does not use them at the moment. It simply uses the
	  //formula: Box fit = Max coordinates (per dimension) - Min coordinatnes (per dimension).


	 //Collection --> collection copying---

	 void Copy_Atomic_Group(const atom_collection&);
       //Copies passed collection of atoms into this collection of atoms.
     void Copy_Atomic_Group(const atom_collection&, int, int, int, const double*); 
       //Copies passed collection of atoms into this collection of atoms.
       //Position shift for copied atoms is used. 
	 void Copy_Atomic_Group(const atom_collection&, int, int, int, const double*, 
		                    const double*);
	  //Also performs an X-Y-Z rotation about the coordinate average of the atoms copied.
	 void Copy_Atomic_Group(const atom_collection&, int, int, int, const double*, 
		                    const double*, bool);
	  //Also performs an X-Y-Z rotation about the coordinate average of the atoms copied.
	  //Allows for the prevention of hydrogen atoms and dummy atoms during copying.
            
     //Atomic name formatting---
       
     void Index_Similar_Atoms();
       //Adds a number to the end of all similar atom names.
       //E.g. Au, Au, Au becomes Au1, Au2, Au3, etc.
     void Remove_Indices();
       //Removes any numbers from the ends of all atom names.     
	 void Rename_Atoms(const char*, const char*);
	   //Finds atoms with the first passed name and replaces them with the second passed name.
          
     //Charge based operations---
       
	 double Return_Total_Charge() const;
	   //Returns the total charge of all atoms in the collection.
	 double Return_Total_Charge(const double*, const double*) const; 
	   //Returns the total charge of atoms in a given spatial region.
	 double Calculate_ES_Potential(bool) const;
	  //Sums up the electrostatic interactions in the atomic collection.
	  //Currently limited to a direct sum.
	 
	 //Translational operations---

     void Shift_Atomic_Coordinates(const double*);
       //Shift a collection of atoms by the passed spatial coordinates.
	 void Shift_Atomic_Coordinates(int, const double*);
	   //Shift only a certain number of atoms at the end of the collection array.
	 void Shift_Atomic_Coordinates(int, int, const double*);
	   //Shifts a certain (index) range of atoms.

	 void Inverse_Shift_Atomic_Coordinates(const double*);
       //Shift a collection of atoms by the passed spatial coordinates.
	   //Moves in the opposite direction (negative direction).
	 void Inverse_Shift_Atomic_Coordinates(int, const double*);
	   //Shift only a certain number of atoms at the end of the collection array.
	   //Moves in the opposite direction (negative direction).
	 void Inverse_Shift_Atomic_Coordinates(int, int, const double*);
	   //Shifts a certain (index) range of atoms.
	   //Moves in the opposite direction (negative direction).

	 void Translate_To_New_Coordinate_Mins(const double*);
	  //Shifts collection of atoms by a certain traslation vector such that
	  //the minimum coordinates change to the passed values.
	 void Translate_To_New_Coordinate_Mins(int, const double*);
	  //End-set limited version.
	 void Translate_To_New_Coordinate_Mins(int, int, const double*);
	  //For an (indexed) range of atoms.

     void Translate_To_New_Coordinate_Average(const double*);
       //Shifts collection of atoms by a certain translation vector such that
       //the coordinate average changes to the passed value.
	 void Translate_To_New_Coordinate_Average(int, const double*);
	   //End-set limited version.
	 void Translate_To_New_Coordinate_Average(int, int, const double*);
	   //For an (indexed) range of atoms.

	 void Move_All_Positions_Into_Box();
	   //Moves all atoms to be located within an orthorhombic periodic box
	   //by translations along box dimensions. This can tear bonds apart
	   //in terms of absolute coordinates, so be cautious in use.

	 //Rotational operations (vec 1 --> vec 2)---

     void Rotate_Atoms(const double*, const double*);
       //Rotate the atoms by vector1 --> vector 2 logic.
	 void Rotate_Atoms(int, const double*, const double*);
	   //Rotates the atoms from the end of the collection.
     void Rotate_Atoms(int, int, const double*, const double*);
	   //Rotates an indexed range of atoms.

	 void Rotate_Atoms(const double*, const double*, const double*);
      //Same as above, but does a rotation about a non-zero origin point.
     void Rotate_Atoms(int, const double*, const double*, const double*);
	  //Rotates the atoms from the end of the collection about an origin.
	 void Rotate_Atoms(int, int, const double*, const double*, const double*);
	  //Rotates the atoms in an indexed range about an origin.

	 //Rotational operations (vector and angle)---

	 void Rotate_Atoms(const double*, double);
	  //Rotates atoms by a rotation vector and rotation angle.
	 void Rotate_Atoms(int, const double*, double);
	  //From the end of the atom collection.
	 void Rotate_Atoms(int, int, const double*, double);
      //For an indexed range of atoms.

	 void Rotate_Atoms(const double*, double, const double*);
      //Rotates atoms by a rotation vector and rotation angle,
	  //with an origin of rotation.
     void Rotate_Atoms(int, const double*, double, const double*);
	  //From the end of the atom collection, with a center of
	  //rotation.
	 void Rotate_Atoms(int, int, const double*, double, const double*);
	  //For an indexed range of atoms, with a center of rotation.

	 //Rotational operations (set of axis rotation-based angles)---

     void Rotate_AtomsXYZ(const double*, bool);
	  //Rotates by a set of three angles (X-Y-Z) by 
	  //clockwise rotation (boolean parameter). 
	 void Rotate_AtomsXYZ(int, const double*, bool);
	  //From the end of the collection.
	 void Rotate_AtomsXYZ(int, int, const double*, bool);
	  //For an indexed range of atoms.
	 
	 void Rotate_AtomsXYZ(const double*, bool, const double*);
	  //Rotates by a set of three angles (X-Y-Z) by clockwise
	  //rotation (boolean parameter).
	 void Rotate_AtomsXYZ(int, const double*, bool, const double*);
	  //From the end of the collection, with an origin of rotation.
	 void Rotate_AtomsXYZ(int, int, const double*, bool, const double*);
	  //For an indexed range of atoms, with an origin of rotation.

	 //Atom searching---

     int Find_Atom_Name(const char*, bool) const;
      /*Looks through a collection of atoms for the first atom with a name
        that matches the passed string. Returns the corresponding atom index.
        Returns -1 if no index is found. Will first look for an exact name
	    match and then precede to consider atom names that only differ by 
		indices contained in the atom name. Finally looks for differences
		to be expected from PDB input (e.g. "MG" will equal "Mg"). Boolean
		parameter indicates if the passed name is an atomic symbol.*/
	 int Find_Atom_Name(int, int, const char*, bool) const;
      //Looks through a collection of atoms for an atom name. Has min/max
	  //arguments.
     int Find_Closest_Atom(const double*) const;
      //Looks through a collection of atoms and finds the atom that is closest
      //to the passed spatial coordinates.
	 int Find_Closest_Atom(int, int, const double*) const;
      //Looks through a collection of atoms and finds the atom that is closest
      //to the passed spatial coordinates. Has min/max index arguments.

	 //Category lists---

     void Get_Name_List(char**, int&, const int) const;
      //Determines an array of unique atom names (e.g. Au, Ag, Cu, etc.). 
      //Names should not be completely indexed before use of this function
	  //as this could cause array overflow.
	 void Get_Index_List(char**, int*, int) const;
      //Determines a reference index for a set of atom names.
      //Returns -1 if no index is found AND outputs a warning string.
	 void Get_Index_List(int*, int&) const;
	  //Determines a reference index for all unique atom types in the collection.
	  //Returns the number of atom types found.
	 void Get_Surface_Atom_List(int*, const int, int&, double, bool, int, bool) const;
	  //Returns a set of indices that correspond to all the surface atoms in the system.
	  //Does so by analyzing which atoms do not have a coordination sphere like that in
	  //the bulk. Removes invisible atoms from consideration.

	 void Get_Charged_Atom_List(int*, const int, int&, bool) const;
	  //Returns an array of atoms that carry non-zero charge. Also removes
	  //invisible atoms from consideration. Meant for use in atomistic simulation.
	 void Get_Molecule_List(int*, const int, int, int&, bool) const;
	  //Returns an array of atoms that have a certain group tag value.

	 //Atom counting---

	 void Get_Composition(char*, const int) const;
	  //Returns the elemental composition (e.g. Si_303,N_404) of the atomic collection.
	  //Second parameter is the length of the passed c-string.
	 int Count_Atom_Type(const char*) const;
	  //Takes an atom name and counts the number of atoms with that name.
	 int Count_Atom_Type(const char*, int*, int) const;
	  //Takes an atom name and counts the number of atoms within a molecular
	  //matrix contained in a bonding matrix.
	 int Count_Unique_Atom_Types() const;
	  //Returns the number of unique atom names in the system.

	 //Coordination analysis---

	 int Get_Coordination_Number(int, double, bool) const;
	  //Determines the coordination number for a given atom, using a
	  //coordination cut-off and an unlike-atom only boolean.
	 int Get_Coordination_Number(int, double, double, bool) const;
	  //Allows a minimum coordination number to be included in the analysis.
	 bool Is_Surface_Atom(int, double, bool, int) const;
	  //Returns TRUEV if the indicated atom is a surface atom. To determine this, the
	  //following assumption is used: Surface atom = Atom with CN less than bulk CN.
	 void Coordex_Atoms(int, int, double, double);
	  //Assigns the coordination number to a range of atoms in the collection based on the 
	  //min/max coordination sphere parameters passed. E.g. "Au" + 3-CN --> "Au3".

	 double Calculate_Surface_Fraction(double, bool, int);
	  //Determines the fraction of atoms that are surface atoms (i.e. that are
	  //experiencing undercoordination).

	 //Large-scale spatial analysis.

     void Calculate_RDF(const char*, const char*, double*, double) const;
      //Calculates the radial distribution function for the atomic collection
	  //for a given atom-atom pair (e.g. O-O RDF) with a given data point
	  //spacing.
	 int Pair_Count(double, double, bool, bool, bool) const;
	  //Returns the number of atomic pairs within a distance range. Allows for 
	  //consideration of like-only and unlike-only interactions.
	 double Distance_Average(double, double, bool, bool, bool) const;
	  //Gives the average distance of atom pairs within a certain distance range.
	  //Also allows for like-only and unlike-only interactions to be considered.
	 double Coordination_Average(double, double, bool, bool, bool) const;
	  //Determines the average coordination number over atoms in the collection.
	 void Coordination_Count(double, double, bool, bool, bool, int*, const int) const;
	  //Counts the number of atoms with each possible coordination number.
	  //E.g. 2 CN-5 atoms, 3 CN-15 atoms, 4 CN-40 atoms, etc. A maximum coordination 
	  //number (1 greater than that used in calculation) is specified.

	 //X-ray diffraction calculation---

	 void Get_Debye_Intensity(double*, double, double, int, double, 
							  const double**, double, bool, double, 
							  double, int) const;
	  /*Get the Debye scattering intensity over a range of angles. Uses
	    the traditional fit of 4 Gaussians to describe atomic scattering 
		factor-scattering angle relationship. Also includes the Debye-Waller 
		temperature factor, which describes thermal attentuation of intensity as:
		Relative attenuation = exp(-2*B*s*s), where s is the scattering vector.
		Also allows for the use of different Debye-Waller factors in describing
		surface atoms*/
      
	 //Spatial violation checks---

     bool No_Covalent_Overlap(int);
      //Given an atom specified by an index, this returns TRUEV if the atom is
	  //not covalently bonded to another atom.
	 bool No_Covalent_Overlap(int, int);
	  //Returns TRUEV if no atom in the passed index range is covalently bonded
	  //to atoms outside that range.
	 bool No_Covalent_Overlap(int, int, int, int);
	  //Checks that atoms in one passed index range are not covalently bonded
	  //to any atoms in another range of atoms (= TRUEV). This function 
	  //will ignore index range overlap of the 2nd passed range by the 1st range.
	 bool No_Covalent_Overlap(int, int, int, int, int&);
	  /*Index range-index range covalent overlap check in which the atom
	    that is in a covalent bond (if found) is stored in the passed parameter.
	    That parameter index comes from the 1st passed index range. Gives -1
	    if no such overlap is found. Also ignores 1st range overlap with 2nd range.*/
	 bool No_Molecule_Overlap();
	   //Returns TRUEV if no atoms from different molecules are covalently bonded.
	   //Ignores atoms with the tag of -1.
	 bool No_Molecule_Overlap(int, int);
	   //Returns TRUEV if no atoms from two different molecules are covalently bonded.
	   //This is an O(n) function so don't use it unless you really need it.

	 void Monolayer_Overlap_Prevention(int, int, int, int, double, int);
	  /*Meant for addition of monolayers, where molecular chains are both
	  added sequentially and must not overlap.
	  Checks for covalent overlap. If found, this function will attempt
	  to reorient the atoms in the second passed range by 
	  rotation along the vector that connects the atoms at the 
	  start and end of the index range (origin in rotation = start of that range). 
	  A reorientation rotation angle and the number of reorientations
	  to attempt to perform must be given. If spatial overlap in the system can
	  not be prevented, both the indicated group and the first group to
	  trigger overlap is removed*/
	 void Remove_Atoms_In_The_Way(int, int, int, int);
	  //Removes atoms in the first passed range that are
	  //exhibiting any covalent overlap with atoms of the 
	  //second passed range.

	 //Neighbor-list generation---

	 //Bond network functions---

	 void Generate_Bond_Network(int*, int, int, int&) const;
	  //Fills the passed integer array with numerical tags that
	  //indicate the molecular groups each atom is a part of.
	  //Number of molecules found is returned.
	 void Get_Molecule_Composition_String(char*, int*, int, const int) const;
	  //Takes a bonding network array and determines a specific bond
	  //network index's molecule name.
	 void Clean_Up_Them_Molecules(int);
	  /*Removes molecular fragments by determining the molecular bonding
	    network and then looking for molecules in small population that
	    appear to be fragments of more commonly found molecules. Parameter
	    passed is the required difference between fragment count and
	    like molecule counts for elimination (recommended = 0).*/

     //Special functions---

	 bool Check_Strange_Location(double, double, bool) const;
	  /*Looks for peculiarities in the atomic spacings in the collections. If 
	    1) any atom has any other atom at a distance less than the passed minimum
		distance (i.e. atoms are too close) or 2) any atom has not even one atom 
		within the second distance parameter (i.e. an atom is off on its own 
		in coordinate space), this function returns TRUEV.*/

     //I/O---
     
	 void Prepare_Collection_For_Output();
	   //Prepares the data for output. Checks for effective floating-point 
	   //zero values and sets them exactly equal to zero. Makes data look much nicer.

     void Print_Atoms_Location() const;
	  //Outputs the spatial locations to DOS prompt.
     void Store_Atoms_Location(ofstream&, bool) const;
	  //Outputs the spatial locations of atoms in the atomic collection.
     void Store_Atoms_Location_With_Charge(ofstream&, bool) const;      
	  //Same as above, but also stores atomic charges (with each atom).
	 void Store_Atoms_Location_With_CoreShell(ofstream&, bool, bool, bool) const;
	  //Stores atom location with cation and/or anion core-shell terms included.

	 void Store_Atom_Types_With_Masses(ofstream&, bool) const;
	   //Stores a list of the atom types with their atomic masses.
         
     void Atoms_Storage(ofstream&) const;
	   //Save a collection of atoms.
     void Atoms_Retrieval(ifstream&);
       //Load a collection of atoms. Appends the atoms to the collection.
       
     private: 

     void Monolayer_Cutter(int, int, int);
	  //Cuts the group of atoms in the range of 2nd to 3rd parameters (one molecule
	  //in the monolayer) and then figures out the molecule containing the first
	  //atom index and cuts that one as well.

     void Memory_Check();
      //Ensures that the initial allocation of memory has taken place.
     void Memory_Check(int);
      //Same as above, but also guarantees that the atom vector is larger 
	  //in size than the passed index.
     
     vector<atom> atom_array;
      //Vector of atoms in the atom collection.
     int num_atoms; 
      //Number of atoms in the collection.
     int initial_allocation;
      //Initial amount of memory to reserve for atom vector.
     bool memory_initialized;
      //Inidicates if initial memory allocation has taken place.

	 int periodicity;
	 double box_dimensions[3];
	  //Periodicity and size parameters of the orthorhombic/rectangle/line
	  //translational periodicity that describes the atomic collection
	  //in space.

     };



#endif