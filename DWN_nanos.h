
#ifndef FLAT_MASA
#define FLAT_MASA

#include <iostream>
#include <fstream>

using namespace std;

#include "DWN_atom.h"
#include "DWN_atomcoll.h"
#include "DWN_cryst.h"
#include "DWN_molecule.h"
#include "DWN_utility.h"

//Nanosurface class---

const long int INIT_SATOM_COUNT = 500;
 //Initial amount of memory space reserved for atoms in the surface.

const double SIZE_VARI_CONST = 0.01;
 //Statistical parameter in the minimization of surface charge by
 //3D surface-thickness variation (not ion removal).

const bool BIND = 1;
const bool FREE = 0;
 //Monolayer coating modes. Bind ensures that molecules in the monolayer
 //are attached/bonded to a specific atom on the surface.

class nanosurface
      {
      public:

      //Constructors, destructors, and initialization functions---

      nanosurface();
      nanosurface(const crystal_system&, double, double, double, 
		          int, int, int, double, bool, bool, bool);
       //Detailed constructor. Determines the atoms of the nanosurface.
      ~nanosurface();  
      
      void Initialize();
      
      //Properties setting---
      
      void Set_Properties(const crystal_system&, double, double, double, int, int, int, 
		                  double, bool, bool, bool);
       //Sets all properties of the nanosurface. This function should be used
	   //prior to calling Construct_Surface().

	  //Nanosurface construction---

	  void Construct_Surface();
       //After the nanosurface properties are set,
       //the nanosurface may be created by calling this function.
      
      //Utility functions---
      
      void Get_Corner(double*) const;
       //Returns the upper-left (minimum) coordinates of the nanosurface.
      void Get_Corner_Atom_Loc(double*) const;
       //Returns the position of the atom closest to the surface corner, with a guarantee
	   //that this atom is in the top surface plane.
      void Get_OpposingCorner_Atom_Loc(double*) const;
       //Gets the position of the atom closest to the opposing surface corner, with a guarantee
	   //that this atom is in the top surface plane.
      void Translate_Surface(const double*);
       //Shifts the location of the nanosurface.
      void ReAssign_Corner(const double*);
       //Reassigns the location of the nanosurface in terms of the new location
       //of the minimum-coordinates corner.
	  void Staircase_Shift(double, double);
	   //Progessively shifts all atomic coordinates as in a staircase by 
	   //Val = [ (z position)/height parameter) * relative shift parameter ]
	   //along the x and y axes.
	  void Random_Shift(const double*);
	   //Shifts all atomic coordinates using ranges for random values.
	  void Slice(double, double);
	   /*Forces all atomic coordinates to have a z-coordinate exactly between
	     2-D slice positions. Useful for image simulation. Assumes the origin of 
	     the surface box to be the location of the first slice. Offset parameter
	     should ideally be less than half the slice thickness and never ever
	     more than that.*/

	  void Rotate_Along_Norm(double);
	   //Rotates the surface slab along its surface normal by the passed angle,
	   //in radians.

      void Get_XY_Repeat(double&, double&) const;
       //Returns the crystal repeat direction along the x and y axes of the surface
       //plane, that is the minimum distance traveled independently along the x or y 
	   //axes from the surface corner atom that finds itself at another atom's position.
      
      //Coating functions---
      
      void Coat_Surface(const molecule&, double, double, double, bool);
       //Coat the nanosurface's top layer with a monolayer.
       
	  //Nanosurface copying---

	  int Copy_NanoSurface(atom_collection&) const;
       //Copies the nanosurface atomic information into an array of atoms.
      int Copy_NanoSurface(atom_collection&, const double*) const;
       //Copies the nanosurface atomic information into an array of atoms.
       //Reassigns the spatial location of that copy.
	  int Copy_NanoSurface(atom_collection&, const double*, const double*) const;
	   //Copies the nanosurface atomic information into an array of atoms.
	   //Reassigns the spatial location of that copy and performs an
	   //X-Y-Z rotation of the surface.
	  int Copy_NanoSurface(atom_collection&, const double*, const double*, bool) const;
	   //Copies the nanosurface atomic information into an array of atoms.
	   //Reassigns the spatial location of that copy and performs an
	   //X-Y-Z rotation of the surface. Allows for the prevention of hydrogen atom
	   //copying into the new nansourface.


      //I/O (units: Angstroms, amu)---
       
      void Print_Atom_Locations() const;
        //Output the names and locations of all atoms to a DOS prompt.
      void Store_Atom_Locations(ofstream&) const;  
        //Output the names and locations of all atoms to file.      
      void Save_NanoSurface(const char*) const;
	   //Saves the nanosurface properties to file.
      void Load_NanoSurface(const char*);
       //Loads the nanosurface properties from file.
      
      private:
            
	  void Neutralize_Surface();
	   /*Attempts to bring the total surface charge as close as possible
	     to zero. There are two modes user can turn on for this:
	     Surface thickness reduction (in all three dimensions) or
	     random ion removal (from top layer).*/
              
      atom_collection surface_atoms;
       //Atomic collection of surface atoms.
      int num_core_atoms;
       //Number of atoms before a monolayer coating is applied to the surface.

      crystal_system cryst_struct;
       //Crystal system that is used to build the nanosurface structure.
      int term_surface_indices[3];
       //Miller indices that define the crystallographic plane that lies along
	   //the Cartesian (0, 0, -1) plane in space.
	  double surface_angle;
	   //The angle of rotation performed along the surface normal during
	   //nanosurface creation.
      double dimensions[3];
       //Desired dimensions of the orthorhombic nanosurface slab.
      bool stuff_space;
       //Indicates if surface creation should try to fit all available space
       //in the orthorhombic box. If set to FALSEV, nanosurface creation will 
	   //not allow the partial addition of basis sets (e.g. 2 out of 4 atoms).
	  bool variable_size;
	   //Indicates if surface creation should vary the dimensions
	   //of the surface until  net charge minimization results.
	  bool neutralize_by_surface;
	   //Indicates if the surface is to be neutralized by random removal of ions
	   //from the top (0, 0, -1) surface layer.
      
      molecule layer_molecule;
       //Holds the prototype of the molecule used for adding monolayers to
	   //the surface.
      };

#endif