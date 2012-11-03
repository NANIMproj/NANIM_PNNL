
#ifndef CRYSTAL_JOSH
#define CRYSTAL_JOSH

//Crystal system class used for handling ordering of atoms when
//constructing systems such as nanoparticles and nanosurfaces---

#include <cmath>
#include <vector>

using namespace std;

#include "DWN_atom.h"
#include "DWN_atomcoll.h"
#include "DWN_utility.h"
#include "RandNum.h"

const int INIT_BASIS_SIZE = 100;
 //Initial number of atoms to reserve in memory space.
const int MAX_FACETING_PLANES = 500;
 //Maximum number of planes to be used in faceting crystalline nanoparticles.
const int MAX_AUTO_PLANE_INDEX = 2;
 //Largest value of a Miller indice used in automatic crystal plane analysis
 //for the crystal system. Not relevant to the user-handled setting of faceting planes.

const double DEFAULT_MESH_VARI = 0.10;
 //Default variation parameter in the mesh grid generation used in construction
 //of amorphous nanoparticles.


class crystal_system
      {
      public:            
      
	  //Constructors and initialization functions---

      crystal_system();

	  void Initialize();
	   //Reserves initial memory for crystal system and 
	   //sets initial values to data members.

	  //Destructors---

      ~crystal_system();
      
      //Lattice vector settings---
      
      void Set_Cubic_Vectors(double);
       // a = a = a
      void Set_Tetragonal_Vectors(double, double);
       // a = a != c
      void Set_Orthorhombic_Vectors(double, double, double);
       // a != b != c
      void Set_Hexagonal_Vectors(double, double);
       // a = a != c, a-> c is 90 degrees and a1 -> a2 is 120 degrees.
      void Set_Vectors(const double*);
	   //Sets the lattice vectors with a 1-dimensional 9 member matrix.
	  void Set_Vectors(const double*, const double*, const double*);
       //Sets the lattice vectors with 3 1-dimensional 3 member matrices.
	  void Set_Vectors(const double*, const double*);
	   //Sets the lattice vectors according to six parameters: three
	   //vector magnitudes and three vector angles.
      void Set_Amor_Vectors(double, double);
       //Sets the mesh grid and variation parameters for making a random close-packing
	   //amorphous structure.
       
      //Basis creation---
      
      void Add_Atom_To_Basis(const atom&, const double*);
	   //Adds the passed atom to the basis at the passed location
	   //which is relative to the basis origin.
	  void Add_Atom_To_Basis(const atom&, double, double, double);
	   //Adds the passed atom to the basis at the passed location
	   //which is relative to the basis origin.
      void Add_Atom_To_Basis_Rel(const atom&, const double*);
       //Adds the passed atom to the basis at the passed relative
	   //location, relative being to the crystal lattice vectors.

	  //Set specific crystal system--- 
       
      void Set_SC_System(const atom&, double); 
	   //Defines the crystal system for a basic simple cubic system
	   //with a one-atom basis.
      void Set_FCC_System(const atom&, double);
	   //Defines the crystal system for a basic face-centered cubic
	   //system with a four-atom basis.
      void Set_BCC_System(const atom&, double); 
       //Defines the crystal system for a basis body-centered cubic
	   //system with a two-atom basis.
	        
      //Information retriveal---
      
	  int Get_Basis_Size() const;
       //Returns the number of atoms in the basis of the crystal system.
	  bool Basis_Contains_Charged_Atom() const;
       //Returns TRUEV if at least one atom in the basis has non-zero charge.

      void Get_Vectors(double*) const;
       //Returns a one-dimensional matrix with nine parameters (see above).
      void Get_Vectors(double*, double*, double*) const;
       //Returns the lattice vectors.
      void Get_Vectors_Mag(double*) const;
       //Returns a one-dimensional matrix with the three lattice vector 
       //magnitudes.
	  void Get_Vector_Mag(double&, int) const;
	   //Returns the magnitude of the indexed vector.
      double Size_A_Step() const;
      double Size_B_Step() const;
      double Size_C_Step() const;
      double Smallest_NonZero_Step() const;
      double Biggest_Step() const;
       //Returns the magnitude of one of the three (A, B, or C) lattice 
       //vectors.

	  //Lattice vector tests/operations---

	  bool Is_OrthoRhombic() const;
	   //Returns true if the crystal system is orthorhombic
	   //(i.e. all lattice vector angles are ninety degrees).
	  bool Is_Full_Step_Length(double, int) const;
	   //Tests if the passed parameter is equal to an integer multiple
	   //of the magnitude of the indicated lattice vector.
	   //For TRUEV: n * lattice_vec_mag = passed parameter, within 1%.
	  bool Is_Half_Step_Length(double, int) const;
	   //Tests if the passed parameter is equal to an integer multiple
	   //plus 1/2 of the magnitude of the indicated lattice vector.
	   //For TRUEV: (n + 1/2) * lattice_vec_mag = passed parameter, within 1%.

	  void Get_Abs_Coors(double*, const double*) const;
       //Transforms relative coordinates into absolute coordinates via
	   //lattice vectors.

	  //Basis addition---

	  void Get_Basis_Shift(double*, const double*) const;
	   //Returns the coordinate shift to move the first atom in the
	   //basis set to the passed coordinates. Returns a zero-vector
	   //if there are no atoms in the basis.
	  void Partial_Basis_Check(atom_collection&, int, int, bool) const;
	   //Checks for a partial addition of a basis set (i.e. only some
	   //atoms were added). If boolean parameter is FALSEV, that
	   //partial addition will be undone.

	  void Add_Basis(atom_collection&, const double*) const;
       //Adds the atom basis to an atomic collection such that the first atom 
	   //in the basis is at the passed coordinates.
	  void Add_Basis_ToPointSet(atom_collection&, const double**, int&) const;
	   //Adds the atoms of the basis to an atomic collection, using one set of
	   //spatial coordinates per atom in the basis.
	   //Also updates the passed index variable by adding the basis size to it.
	  void Add_Basis_ToPointSet_ImposeSphere(atom_collection&, 
		                         const double**, int&, double, double, bool) const;
	   /*Adds the atoms of the basis to an atomic collection, using one set of
	     spatial coordinates per atom in the basis.
	     Updates the index variable and imposes a spherical distance minimum/maximum
		 range upon basis atom addition.*/

      void Add_Basis_ToSParticle(atom_collection&, const double*, double, double, bool) const;
       //Adds the atoms of the basis to an atomic collection, imposing a spherical
	   //distance min/max range requirement.
      void Add_Basis_ToFParticle(atom_collection&, const double*, double, double, bool) const;
	   //Adds the atoms of the basis to an atomic collection, imposing a multi-faceted
	   //surface within which all atoms must reside (or be removed).
      void Add_Basis_ToCParticle(atom_collection&, double, double, bool) const;
	   //Adds the atoms of the basis to an atomic collection, using semi-random placement
	   //of atoms to generate the random-close packed configuration of a spherical particle.
	   //The particle will potentially grow off of any atoms in the atomic collection passed.
      void Add_Basis_ToSurface(atom_collection&, const double*, const double*, bool) const;
	   //Adds the atoms of the basis to an atomic collection, imposing an orthorhombic
	   //box in which all atoms added must be located.
	   
      //Lattice vector stepping---
      
      void A_Step(double*) const;
      void B_Step(double*) const;
      void C_Step(double*) const;
       //Apply one periodic step along a lattice vector to the passed 
       //Cartesian coordinates.        
      void A_Step(double*, int) const;
      void B_Step(double*, int) const;
      void C_Step(double*, int) const;
       //Apply integer number of steps along a lattice vector to the passed 
       //Cartesian coordinates.
      void Multi_Step(double*, int, int, int) const;
       //Apply integer number of steps along the three lattice vectors to
	   //the passed Cartesian coordinates.
	  void Multi_Step(double*, int*) const;
	   //Apply integer number of steps along the three lattice vectors to
	   //the passed Cartesian coordinates.
       
     //Cell rotation---  
       
      void Rotate_Planes(const int*, const int*);
       //Rotates the cell such that the plane indicated by the second set of Miller
	   //indices spatially occupies (after rotation) the planar space of the plane 
	   //indicated by the first set of Miller indices.
      void Rotate_Cell(const double*, double);
       //Rotates the cell along a rotation axis by a rotation angle. 
       
      //Faceting plane determination---
       
      int Get_Facet_Count() const; 
	   //Returns the number of faceting planes declared for use
	   //in faceting atomic surfaces.
      void Set_Faceting_Plane(int, int, int, double);
       //Sets the n-th faceting plane.
	  void Get_Faceting_Plane(int, int*, double&) const;
       //Returns the n-th declared faceting plane.
      void Determine_Normal_Vector(const int*, double*) const;
       //Determines the normal vector to the plane indicated 
       //by Miller indices. 
       
      //Supercell logic---
      
      int Add_SuperCell(atom_collection&, int, int, int) const;
       //Adds a supercell out of this crystal system to the atomic collection. 
	   //Returns the number of atoms in the supercell.
       
      //Operator overloads---
      
      void operator = (const crystal_system&); 
	   //Crystal system takes on all properties of another crystal system.
      
      //File I/O---
      
	  void Save_Crystal_System(const char*) const;
	   //Saves the crystal system to file, as specified
	   //by file name.
      void Load_Crystal_System(const char*);
	   //Loads the crystal system from file, as
	   //specified by file name.
	  void Save_Crystal_System(ofstream&) const;
       //Saves the crystal system to file.
      void Load_Crystal_System(ifstream&);
	   //Loads the crystal system from file.

      void Save_SuperCell(const char*, int, int, int) const;
       //Makes and saves an a x b x c supercell of atoms to file.
       
      private:
                    
      void Calculate_ClosePacked_Planes();
       /*Determines the close-packed planes for the crystal system
	     and assigns them to the faceting planes. This is performed
	     only if faceting is requested with no faceting planes yet declared.
	     Note: Only uses one basis of atoms in analysis. This is meant
	     to be something more interesting than practical. It is not
		 currently a part of main program usage.*/
               
      atom_collection basis_atoms;
       //The basis set of atoms for the crystal system.
      
      double lattice_vectors[3][3];
       //The vectors that define the crystallographic unit cell.  
      
      int faceting_planes[MAX_FACETING_PLANES][3];
       //Miller indices for faceting planes.
      int num_faceting_planes;
       //Number of faceting planes declared.
      double scaling_factors[MAX_FACETING_PLANES];
	   //Relative placement of faceting planes when creating atomic surfaces.
	   //A higher scaling factor corresponds to a faceting plane farther 
	   //away from the origin of the system being faceted.

      };      
      
#endif