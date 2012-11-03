
#ifndef NANOPART_MISA
#define NANOPART_MISA

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#include "DWN_atom.h"
#include "DWN_atomcoll.h"
#include "DWN_cryst.h"
#include "DWN_molecule.h"
#include "DWN_utility.h"

//Nanoparticle class------

const long int INITIAL_PATOM_COUNT = 500;
 //Initial amount of memory reserved for atoms in the nanoparticle.

const int CRY_MODEL = 0;
const int CON_MODEL = 1;
const int SPH_MODEL = 2;
 /*Values for the continuous, crystalline, and spherical models of nanoparticles:
   Continuous model --> Spherical area filled with random close-packed atoms.
   Spherical model --> Spherical area filled with crystalline material
   Crystalline model --> Faceted area filled with crystalline material */

const bool NO_FIRST_ATOM_METHOD = TRUEV;
 //Turns off the "first atom removed is from the (1, 1, 1) corner" method
 //of charge neutralization of crystalline nanoparticles.

class nanoparticle
      {
      public:
       
      //Constructors, destructrors, and initialization functions---

      nanoparticle();
      nanoparticle(double, double, int, const crystal_system&, bool, bool);
	   //Constructor with all properties defined.
      ~nanoparticle();   

	  void Initialize();
       //Initializes nanoparticle properties.
      
      //Properties setting---
      
      void Set_Properties(double, double, int, const crystal_system&, bool, bool);
       //Sets all properties of the nanoparticle. Carries out construction
       //of the list of atoms in the nanoparticle as well.

	  //Particle construction---

	  void Construct_Particle();
	   //Constructs the nanoparticle (should be called after the properties are
	   //set).       

	  //Coating functions---
      
      void Coat_Spherical_Particle(const molecule&, double, double, bool); 
       /*Coat the particle of a nanosurface with a monolayer.
       Requires the particle-head spacing and head-head spacing.
       Here, head-head spacing is less clearly enforced as with surface
       coating. The algorithm attempts to find an atom after moving an arc
       length equal to the head-head spacing when placing new molecules if the 
       last parameter is set to TRUE.*/
      void Coat_Faceted_Particle(const molecule&, double, double, bool);
       //Similiar as above, but coats the faces of a faceted nanoparticle.
      
	  //Utility functions---
      
      void Get_Center(double*) const;
       //Obtains the coordinate average of nanoparticle atoms.
      void Translate_Particle(const double*);
       //Shifts the location of the nanoparticle.
      void Zero_Center();
	   //Zeros the coordinate average of the nanoparticle.
	  void ReAssign_Center(const double*);
       //Assigns a new coordinate average of the nanoparticle.
	  void Rotate(const double*);
	   //Performs an X-Y-Z rotation.

	  //Nanoparticle copying---

	  int Copy_NanoParticle(atom_collection&) const;
       //Adds the nanoparticle to an array of atoms.
      int Copy_NanoParticle(atom_collection&, const double*) const;
       //Adds the nanoparticle to an array of atoms, reassigning the nanoparticle
	   //coordinate average.
	  int Copy_NanoParticle(atom_collection&, const double*, const double*) const;
	   //Adds the nanoparticle to an array of atoms, reassigning the nanoparticle
	   //coordinate average and performing an X-Y-Z rotation.

      //I/O (units: Angstroms, amu)---
       
      void Print_Atom_Locations() const;
       //Output the names and locations of all atoms to a DOS prompt.
      void Store_Atom_Locations(ofstream&) const;  
       //Output atomic information to file.      
      void Save_NanoParticle(const char*) const;
	   //Saves a nanoparticle to file.
      void Load_NanoParticle(const char*);
       //Loads a nanoparticle from file.  
       
      private:
	  
      void Construct_Conti_Particle();      
      void Construct_Cryst_Particle();
       //Particle construction functions for different model types.

	  void Neutralize_Particle();
	    /*Neutralize the particle charge. Tries to remove the first atom 
	     from the (1, 1, 1) particle corner if model = CRYST, 
		 and then at random from the particle's surface.
	     Credit for idea: Josh Deetz*/	
      
      atom_collection particle_atoms;
       //Atoms that make up the nanoparticle.
      
      double effective_diameter;
       //Effective diameter of the nanoparticle.
	  double hole_diameter;
	   //Effective diameter of the central hole within the nanoparticle 
	   //(no hole = 0.0 size).
      double effective_radius;
       //Effective radius of the nanoparticle.
      bool stuff_space;
       //Indicates if particle creation can partially add basis sets
	   //during particle construction (e.g. adding 2 out of 4 atoms
	   //in the basis).
	  bool neutralize_particle;
	   //Charge neutralization boolean.
      crystal_system cryst_struct;
       //Crystal structure used to build crystalline nanoparticles.
      int model_type;
	   //Type of model used for the atomic arrangement (amorphous, crystalline,
	   //or spherical).
       
      molecule layer_molecule;
       //Molecule used in coating the nanoparticle surface.               
      };
      
#endif
