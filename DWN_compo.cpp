
#include "stdafx.h"

#include "DWN_compo.h"

//Constructors and initialization---

general_composite::general_composite()
      {
      Initialize();                            
      }
      
general_composite::general_composite(const double* size)
      {
      Initialize();
      Resize(size);                              
      }
      
general_composite::general_composite(const double* size, int period)
      {
      Initialize();
      Resize(size);   
      Set_Periodicity(period);                             
      }   

void general_composite::Initialize()
      {
      atoms.Set_Initial_Allocation(INIT_CNUM_ATOMS);
      Resize(10.0, 10.0, 10.0);
      Set_Periodicity(0);
	   //Assume a OD system for initialization.
      group_counter = 0;    
	   //First group/molecule added should have an index of 0.

	  sub_mode_set = FALSEV;
	  bulldoze_mode_set = FALSEV;
	  invsub_mode_set = FALSEV;
	   //All alternative addition styles should be set to FALSEV.
      }

//Destructor---
      
general_composite::~general_composite()
      {                                
      }   

//Information retriveal--

int general_composite::Get_Atom_Count() const
     {
	 return (atoms.Size());
     }

double general_composite::Get_Charge() const
     {
	 return (atoms.Return_Total_Charge());
     }

void general_composite::Get_Box_Size(double* sizes) const
     {
	 atoms.Get_Box_Vectors(sizes);
     }

int general_composite::Get_Periodicity() const
     {
     return atoms.Get_Periodicity();                                           
     }

//Addition of systems (e.g. nanoparticles, nanosurfaces, etc.)
//to the composite---

void general_composite::Add_Atom(const atom& ze_atom, const double* pos)
     {
	 atoms.Add_Atom(ze_atom, pos, PLACE_MODE);
     Add_Group(atoms.Size() - 1, atoms.Size(), FALSEV); 
	  //Atom gets a group tag value of -1 (not a group).
     }

void general_composite::Add_Atom_Collection(const atom_collection& ze_atoms, 
                                            const double* pos, bool is_group)
     {	 
     int old_composite_size = atoms.Size();
	 double min_vals[3];
	 ze_atoms.Get_Coordinate_Mins(min_vals);
	 double shift_vals[3];
	 SetSub_XYZ(shift_vals, pos, min_vals);
	  //Now the shift vector is: new corner - old corner.
     for (int a = 0; a < ze_atoms.Size(); ++a)
         {
		 atoms.Add_Atom(ze_atoms[a], shift_vals, SHIFT_MODE);
         }
     Add_Group(old_composite_size, atoms.Size(), is_group);                                             
     }

void general_composite::Add_NanoParticle(const nanoparticle& parti, 
                                         const double* pos, const double* rot_angles, 
										 bool is_group)
      {
      int old_composite_size = atoms.Size();
      parti.Copy_NanoParticle(atoms, pos, rot_angles);  
       //Use nanoparticle's copier.
      Add_Group(old_composite_size, atoms.Size(), is_group);
       //Define the group identity for the atoms added.
      }
      
void general_composite::Add_NanoSurface(const nanosurface& surf, const double* pos, 
                                        const double* rot_angles, bool is_group, bool stop_hydro)
      {
      int old_composite_size = atoms.Size();
      surf.Copy_NanoSurface(atoms, pos, rot_angles, stop_hydro); 
	   //Use nanosurface's copier.
      Add_Group(old_composite_size, atoms.Size(), is_group);    
      }
      
void general_composite::Add_Molecule(const molecule& mol, const double* pos, 
                                     const double* rot_angles, bool is_group,
                                     bool use_first_atom)
      {
      int old_composite_size = atoms.Size();
      mol.Copy_Molecule(atoms, pos, rot_angles, use_first_atom);  
	   //Use molecule's copier.
      Add_Group(old_composite_size, atoms.Size(), is_group);                        
      }
      
void general_composite::Add_SuperCell(const crystal_system& cryst, const double* pos, 
                                      int x, int y , int z, bool is_group)
      {
      int old_composite_size = atoms.Size();
      int size_change = cryst.Add_SuperCell(atoms, x, y, z);
	   //Supercell maker returns the number of atoms in the supercell.
      atoms.Shift_Atomic_Coordinates(old_composite_size, old_composite_size + size_change, pos);
	   //Shift supercell atoms into requested position.
      Add_Group(old_composite_size, atoms.Size(), is_group);
      } 
       
void general_composite::FillSpace_With_Solvent(const molecule& mol, 
                                               const double* box_location,
                                               const double* box_size,
											   const double* solvent_spacing,
											   bool is_group)
     {
     double tp[3];
      //Test position for solvent molecule placement.
	 double angles[3] = { 0.0, 0.0, 0.0 };
      //Angles for rotating molecule as it is placed down, which is
	  //always zero degrees here.

	 double end_corner[3];
	 SetAdd_XYZ(end_corner, box_location, box_size);
	  //Get the opposing corner of the solvation box.
     
     double start_vals[3], end_vals[3];
     Get_Spatial_Sampling_Pams(box_location, end_corner, solvent_spacing, 3, start_vals, end_vals);
	  //Determine the correct start/end positions for having a solvent placement that is symmetrical
	  //with respect to the solvation box and the solvent molecule centers.
     
     int original_size = atoms.Size();
      //Size of the composite before adding any solvent molecules.
	 int molecule_size = mol.Size();
     for (tp[0] = start_vals[0]; tp[0] < end_vals[0]; tp[0] += solvent_spacing[0])
         {
         for (tp[1] = start_vals[1]; tp[1] < end_vals[1]; tp[1] += solvent_spacing[1])        
             {
             for (tp[2] = start_vals[2]; tp[2] < end_vals[2]; tp[2] += solvent_spacing[2])
			  //Check each point in the orthorhombic solvent placement box...
                 {
				 Add_Molecule(mol, tp, angles, is_group, MOLECULE_CENTER);       
                 if (!atoms.No_Covalent_Overlap(0, original_size, atoms.Size() - molecule_size, atoms.Size()))
				  //Check for covalent overlap between the original composite atoms (no solvent considered)
				  //and the molecule just added. If there is such overlap, delete the molecule.
					{
					atoms.Delete_Atoms(atoms.Size() - molecule_size, atoms.Size());
					}
                 }        
             }
         }                                    
     }

//Complementary atom addition/removal logic---

void general_composite::Charge_Balance(const molecule& balancer, const double* box_pos, 
	                                   const double* box_size, bool is_group)
     {
     double charge_to_add = -1.0*Get_Charge();
	  //Balance by adding the opposite charge of the composite to the composite.
	 int number_to_add;
	 Safe_Truncation(number_to_add, charge_to_add, balancer.Get_Charge());
	  //Get number of charge balancing molecules to add in order to obtain neutral charge:
	  //Number to add = Charge to add/Molecule charge.
	 if (number_to_add < 0)
	  //Trying to balance with a molecule that has the wrong charge,
	  //e.g. using Cl- to balance a negative system.
	    {
		return;
	    }

	 double random_point[3];
	  //Used to store a random Cartesian position to attempt
	  //to place a charge-balancing molecule at.
	 int try_count;
	  //Number of tries made at placing a specific molecule.
	 const int MAX_TRY_COUNT = 10000;
	  //Maximum number of tires at adding one molecule
	  //before giving up (i.e. there is no open space for 
	  //the charge-balancing molecule if this happens).

	 double rot_angles[3] = { 0.0, 0.0, 0.0 };
	  //Rotation for adding molecules to the composite,
	  //which is here set to zero.

	 bool overlap_problem;
	 int molecule_start_index;
	 for (int a = 0; a < number_to_add; ++a)
	    {
	    overlap_problem = TRUEV;
		try_count = 0;
        do
		  {
		  Get_Rand_Box_Position(random_point, box_pos, box_size);
		  molecule_start_index = atoms.Size();
		  Add_Molecule(balancer, random_point, rot_angles, is_group, MOLECULE_CENTER);     
          if (!atoms.No_Covalent_Overlap(0, molecule_start_index, molecule_start_index, atoms.Size()))
			//Check for covalent overlap problems with the balancing molecule being added.
			    {
			    atoms.Delete_Atoms(molecule_start_index, atoms.Size());
				++try_count;
				}
		  else
			//No covalent overlap problem. Keep the charge-balancing molecule
			//and move on to the next one to be added.
				{
				overlap_problem = FALSEV;
				}
		  if (try_count == MAX_TRY_COUNT)
		   //Abandon this loop if open space for charge-balancing molecule addition
		   //can not be found.
			 {
			 Show_Warning("CHARGE BALANCER COULD NOT FIND SPACE FOR ADDING COUNTER-IONS");
			 overlap_problem = FALSEV;
			 a = number_to_add;
			 }
		  }
		  while (overlap_problem);
        }
     }

void general_composite::Charge_Balance(const double* box_pos, const double* box_size)
	 {
	 double box_max[3];
	 Set_XYZ(box_max, box_pos);
	 Add_XYZ(box_max, box_size);
	  //Obtain opposite corner of the box.

	 int random_index;
	  //Index for random selection of an atom within the box
	  //to eliminate.
	 int try_count = 0;
	  //Number of tries at finding an atom within the box to
	  //be removed.
	 const int MAX_TRY_COUNT = 10000;
	  /*Maximum number of tries at atom removal before giving up.
	    Reasons for failure on a test include: 1) atom not
	    within the box and 2) atom does not have the right sign
	    of chare*/

	 const double CHARGE_BALANCE_DONE = 0.3; 
	  //Absolute charge of the composite needed for the charge balancing
	  //to be considered finished.
	 double atom_loc[3];
	 double current_charge = Get_Charge();
	 while (abs(current_charge) > CHARGE_BALANCE_DONE)
	    {
		random_index = GetRandPosNum(long int(atoms.Size() - 1));
		atoms[random_index].Get_Atom_Location(atom_loc);
		++try_count;
	    if (Check_Boundaries(atom_loc, box_pos, box_max)
			&& Check_FP_GreaterThan(atoms[random_index].Get_Atom_Charge()/current_charge, 0.0 ))
		 //If atom is in the requested removal zone and has the right sign of charge, cut it.
           {
		   current_charge -= atoms[random_index].Get_Atom_Charge();
           atoms.Delete_Atom(random_index);
		   try_count = 0;
           }
		if (try_count == MAX_TRY_COUNT)
		 //Abandon this loop if an atom for removal can not be found.
		   {
		   Show_Warning("CHARGE BALANCER COULD NOT FIND ATOMS FOR REMOVAL WITHIN THE SPECIFIED BOX!");
		   current_charge = 0.0;
		   }
        }
     }

void general_composite::Clean_Composite(int clean_pam)
     {
	 atoms.Clean_Up_Them_Molecules(clean_pam);
     }


void general_composite::Make_Void(const double* origin, const double* sizes, bool pbc)
     {
	 bool is_sphere = FALSEV;
	 if (Check_FP_Equality(sizes[1], 0.0) && Check_FP_Equality(sizes[2], 0.0) )
	  //Spherical void is requested by specifying only one parameter in sizes
	  //array.
	    {
	    is_sphere = TRUEV;
	    atoms.Make_Spherical_Void(origin, sizes[0], pbc);
	    }
	 else
	    {
		atoms.Make_Orthorhombic_Void(origin, sizes);
	    }
     }

//Non-normal component addition modes--- 

void general_composite::Temp_Set_Sub_Mode()
	 {
	 sub_mode_set = TRUEV;
	 }

void general_composite::Temp_Set_Bulldoze_Mode()
     {
	 bulldoze_mode_set = TRUEV;
     }

void general_composite::Temp_Set_InvSub_Mode()
     {
	 invsub_mode_set = TRUEV;
     }

//Box control---

void general_composite::Resize(double x, double y, double z)
     {
	 double vecs[3] = { x, y, z };
	 Resize(vecs);
     }

void general_composite::Resize(const double* size)
     {
	 atoms.Set_Box_Vectors(size);
     }

void general_composite::Set_Periodicity(int period)
     {
	 if (period > 3)
		{
		Show_Congratulations("YOU ARE WANDERING INTO THE FOURTH-DIMENSION!");
		period = 3;
		}
	 if (period < 0)
		{
		Show_Warning("I AM SORRY! I don't know how to do negative dimensions.");
		period = 0;
		}
	 atoms.Set_Periodicity(period);
     }

//Special Effects---

void general_composite::Coordex_Me(int start, int end, double cn_min, double cn_max)
     {
	 atoms.Coordex_Atoms(start, end, cn_min, cn_max);
     }

//Molecular Mechanics/Dynamics Initiation---

void general_composite::Rigid_Minimization(int start, int end, bool include_hydro, 
	                                       double precision)
     {
	 //atoms.Rigid_Minimization(start, end, include_hydro, precision);
     }

//File I/O---

void general_composite::Save_Composite(const char* save_file) const
   {
   ofstream ze_file(save_file, ios::out);
   if (ze_file.is_open())
      {
      Write_Intro(ze_file, "File Type = Composite Description");
      ze_file << atoms.Get_Periodicity();
	  double box_vec[3];
	  atoms.Get_Box_Vectors(box_vec);
      Write_Values(ze_file, box_vec, 3);
	   //Save box parameters.
      ze_file << endl;
      atoms.Atoms_Storage(ze_file); 
	   //Store atoms.
	  ze_file.close();
      }          
   }
   
void general_composite::Load_Composite(const char* load_file, const double* pos, 
	                                   const double* rot_angles, bool at_center)
   {
   ifstream ze_file(load_file, ios::in);
   int start_index = Get_Atom_Count();
   if (ze_file.is_open())
      {
      Skip_Phrases(ze_file, 5);
       //Get rid of introductory text.
	  int periodicity;
      ze_file >> periodicity; 
	  Set_Periodicity(periodicity);
	   //Get periodicity.
	  double box_vec[3];
      Load_Values(ze_file, box_vec, 3);
	  Resize(box_vec);
	   //Get composite box size.
      atoms.Atoms_Retrieval(ze_file);  
      Skip_Phrases(ze_file, 2);
	  ze_file.close(); 
      }

   int end_index = Get_Atom_Count();

   atoms.Rotate_AtomsXYZ(start_index, end_index, rot_angles, CLOCKWISE);

   if (at_center)
    //Center composite request.
      {
	  atoms.Translate_To_New_Coordinate_Average(start_index, end_index, pos);
      }
   else
	 //Shift composite corner request.
      {
	  atoms.Translate_To_New_Coordinate_Mins(start_index, end_index, pos);
      }

   Show_Statement("Done loading composite of size: ", atoms.Size());
    //Inform user of the successful load!
   }

void general_composite::Get_File_Runners(atom_collection& ze_atoms) const
     {
     ze_atoms = atoms;        
	  //Get a copy of the atomic collection.

	 ze_atoms.Prepare_Collection_For_Output();
	  //Check floating-point values to ensure good file output.
	  //Also defragment the group tags.
     }
              
//Private functions---

void general_composite::Add_Group(int start, int end, bool is_group)
     {
	 if (is_group)
		{
		atoms.Set_Group_Tag(start, end, group_counter);
		++group_counter;
		}
	 else
	  //Make sure group tag is set to -1 (no group) if
	  //grouping is not desired.
		{
		atoms.Set_Group_Tag(start, end, -1);
		}

	 if ( (sub_mode_set == TRUEV) || (bulldoze_mode_set == TRUEV) )
		 {
		 atoms.Remove_Atoms_In_The_Way(0, start, start, end);
		 sub_mode_set = FALSEV; 
		 if (bulldoze_mode_set == TRUEV)
		  //Remove the atoms just added as the bulldozer does 
		  //not stick around.
		    {
		    int number_to_cut = end - start;
			atoms.Delete_Atoms(atoms.Size() - number_to_cut, atoms.Size());
			 //Atoms of the bulldozer are at the end of the atomic collection
			 //before removal.
		    bulldoze_mode_set = FALSEV;
		    }
		 }

	 if (invsub_mode_set == TRUEV)
	     {
		 atoms.Remove_Atoms_In_The_Way(start, end, 0, start);
		 invsub_mode_set = FALSEV;
	     }
     }
