
#ifndef ELECTRON_SHEEP
#define ELECTRON_SHEEP

#include <vector>

using namespace std;

#include "DWN_atom.h"
#include "DWN_atomcoll.h"
#include "DWN_cryst.h"
#include "DWN_nanop.h"
#include "DWN_nanos.h"
#include "DWN_molecule.h"

const int INIT_CNUM_ATOMS = 2000;
 //Initial allocation of memory space for atoms in the
 //composite = 2000 atoms.

class general_composite
      {
      public:
             
	  //Constructors and initialization functions---

      general_composite();
      general_composite(const double*);
      general_composite(const double*, int);
       //Constructors.
      void Initialize();
       //Data initialization.
      
	  //Destructors---

      ~general_composite();
       //Destructor.        

      //Information retriveal---

	  int Get_Atom_Count() const;
	   //Returns the number of atoms in the composite.
	  double Get_Charge() const;
	   //Returns the total charge of all atoms in the composite.
       
	  void Get_Box_Size(double*) const;
       //Gets the composite's orthorhombic box size parameters. 
	  int Get_Periodicity() const;
       //Returns the periodicity of the box.

      //Addition of systems (e.g. nanoparticles, nanosurfaces, etc.)
	  //to the composite---
      
      void Add_Atom(const atom&, const double*);
       //Adds a single atom at the passed coordinates.
      void Add_Atom_Collection(const atom_collection&, const double*, bool);
       //Adds a collection of atoms and shifts all atomic locations by the
	   //passed coordinate shift.
      void Add_NanoParticle(const nanoparticle&, const double*, const double*, bool);
       //Nanoparticle addition with particle center positioned at passed coordinates.
      void Add_NanoSurface(const nanosurface&, const double*, const double*, bool, bool);
       //Nanosurface addition with effective (0, 0, 0) corner placed at passed coordinates.
	   //Prevention of hydrogen addition is enabled.
      void Add_Molecule(const molecule&, const double*, const double*, bool, bool);
       //Adds a molecule with center placed at passed coordinates.
       //Also rotates the molecule by the three angles passed (X-Y-Z rotation).
	  void Add_Molecule(const molecule&, const double*, double, double, bool, bool, double);
	   //Same as above, but includes a minimum distance parameter to be used in basic
	   //spatial violation analysis. If any atoms outside the molecule are within that
	   //distance from any atoms within the molecule, no addition takes place.
      void Add_SuperCell(const crystal_system&, const double*, int, int, int, bool); 
       //Adds a supercell of a crystal system with its effective (0, 0, 0) corner 
	   //being positioned at the passed coordinates. 
	  void FillSpace_With_Solvent(const molecule&, const double*,  
                                  const double*, const double*, bool); 
	   /*Fills a specified spatial region with solvent molecules placed
	     at a certain distance from each other within an orthorhombic box
		 positioned with passed corner/box size values. Minimum distance
		 between solvent molecule center and non-solvent molecules must also
		 be specified*/
       
	  //Complementary atom addition/removal logic---

	  void Charge_Balance(const molecule&, const double*, const double*, bool);
	   //Fills a box-like area with the charged molecule passed until charge balance 
	   //(i.e. total charge of composite = 0) is best achieved. 
	   //Uses random placement with covalent overlap checking.
	  void Charge_Balance(const double*, const double*);
	   //Balances the composite's total charge by removing atoms found within a 
	   //given box at random until the total absolute charge is less than 0.3.
	  void Clean_Composite(int);
	   //Removes any fragmented molecules in the composite.
	   //E.g. If a composite has 2 OH molecules and 100 H2O molecules, the 2 OH molecules
	   //are removed.
	  void Make_Void(const double*, const double*, bool);
	   //Deletes all atoms in the space specified by the passed origin and sizes.
	   //If size parameters 2 and 3 are set to 0, a spherical void is created
	   //with that radius. Else, an orthorhombic void is created.
      
	  //Non-normal component addition modes--- 
             
	  void Temp_Set_Sub_Mode();
	   /*Turns on substitution mode: the next (and only the next) addition
	     command will remove atoms in the composite that would covalently 
	     overlap with the group of atoms being added. In other words,
	     atoms "in the way" of the next addition are simply removed.*/
	  void Temp_Set_Bulldoze_Mode();
	   //Turns on bulldoze mode: This is just like substitution mode except
	   //that the next addition does not even add the new group of atoms.
	  void Temp_Set_InvSub_Mode();
	   //Turns on inverse substitution mode: The next addition command will
	   //remove atoms in the new group of atoms that would covalently overlap
	   //with atoms in the composite, i.e. the opposite of substitution.

	  //Box control---
              
	  void Resize(double, double, double);
	   //Resizes the box with a 3 parameter vector.
      void Resize(const double*);
       //Resizes the box with a 3 parameter vector.
      void Set_Periodicity(int);
       //Sets the periodicity (0 to 3) of box.

	  //Special Effects---

	  void Coordex_Me(int, int, double, double);
	   //A range of atoms have their atomic name indexed with their
	   //coordination numbers. E.g. "Au" + 3 atoms nearby = "Au3."

	  //Molecular mechanics/dynamics initiation---

	  void Rigid_Minimization(int, int, bool, double);
	   //Starts rigid body energy minimization (which uses translation and rotation
	   //of the body to find a minimum via conjugate gradient).
      
      //File I/O---
      
      void Save_Composite(const char*) const;
	   //Saves the composite to a file.
      void Load_Composite(const char*, const double*, const double*, bool);
	   //Loads a composite from file. Positions composite by passed coordinates,
	   //which can serve as the new minimum-coordinates (per dimension) or new
	   //coordinate center of the composite. Also rotates the composite.

	  void Get_File_Runners(atom_collection&) const;
       //Get composite box size and atomic collection for use in data analysis.
      
             
      private:

	  //Private functions---
          
      void Add_Group(int, int, bool);
	   //Assigns a group tag value (determined by composite) to the atoms
	   //in the passed index range, or gives a group tag of -1 to the 
	   //atoms if they are not to be part of a group.  
      
      //Data members---
              
      atom_collection atoms;
	   //The collection of all atoms in the composite.
      int group_counter;
       //Number of atom groups/molecular groups added to the composite.

	  bool sub_mode_set;
	   //Substitution mode active boolean.
	  bool bulldoze_mode_set;
	   //Bulldoze mode active boolean.
	  bool invsub_mode_set;
       //Inverse substitution mode active boolean.
      };

#endif
