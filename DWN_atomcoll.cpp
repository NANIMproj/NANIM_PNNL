
#include "stdafx.h"

#include "DWN_atomcoll.h"

//Constructors/destructors---

atom_collection::atom_collection()
     {
     num_atoms = 0;
     memory_initialized = FALSEV;
	 No_Periodicity();
	 Set_Initial_Allocation(DEFAULT_INITIAL_MEMORY);           
     }
     
atom_collection::atom_collection(int memory_allocation)
     {
     num_atoms = 0;
     memory_initialized = FALSEV;
	 No_Periodicity();
     Set_Initial_Allocation(memory_allocation);
     }
     
atom_collection::~atom_collection()
     {                         
     }

//Memory allocation---

void atom_collection::Set_Initial_Allocation(int val)
     {
     initial_allocation = val; 
     Memory_Check();                      
     }
     
void atom_collection::Memory_Check()
     { 
     if(!memory_initialized)
       {
       atom_array.resize(initial_allocation);
       memory_initialized = TRUEV;
       }                              
     }

void atom_collection::Memory_Check(int access_value)
     {    
     Memory_Check();
	  //Initial allocation check.
     if (access_value >= int(atom_array.size()))
      //Memory needs to be expanded. 
       {
       int new_space = atom_array.size() * 2 + 1;
       while (access_value >= new_space)
		//Double the available memory until the size
		//of the vector is greater than the requested
		//access value.
             {
             new_space *= 2;           
             }
       atom_array.resize(new_space);
       }
     }

//Size and atom retrieval---
     
int atom_collection::Size() const
    {              
    return num_atoms;                                      
    }

int atom_collection::Real_Size(bool include_hydro) const
    {
    int num_real_atoms = 0;
	for (int a = 0; a < num_atoms; ++a)
	    {
        if (!atom_array[a].Is_Invisible(include_hydro))
		    {
			++num_real_atoms;
		    }
		}
	return num_real_atoms;
    }

atom& atom_collection::operator [] (int index)
    {
    Memory_Check(index);
    return (atom_array[index]);              
    }

atom atom_collection::operator [] (int index) const
    {
    return (atom_array[index]);              
    }

const atom& atom_collection::Get_Atom(int index) const
    {
    return (atom_array[index]);                                  
    }

//Copying and addition of atoms---

void atom_collection::operator = (const atom_collection& ze_coll)
     {
     initial_allocation = ze_coll.atom_array.capacity();
     atom_array.resize(initial_allocation);
     memory_initialized = TRUEV;
     for (int a = 0; a < ze_coll.num_atoms; ++a)
         {
         atom_array[a] = ze_coll.atom_array[a];     
         }
	 if (num_atoms > ze_coll.num_atoms)
	     {
		 Delete_Atoms(ze_coll.num_atoms, num_atoms);
		 }
	 num_atoms = ze_coll.num_atoms;
	 Set_Box(ze_coll.periodicity, ze_coll.box_dimensions);
     }

void atom_collection::Add_Atom(const atom& ze_atom)
     {
	 Memory_Check(2*num_atoms);
     atom_array[num_atoms] = ze_atom;
	  //Place new atom at end of atom array.
	 ++num_atoms;
     }

void atom_collection::Add_Atom(const atom& ze_atom, const double* coors, bool mode)
     {
	 Memory_Check(2*num_atoms);
	 atom_array[num_atoms] = ze_atom;
	 if (mode == SHIFT_MODE)
	  //Shift coordinates.
	    {
		atom_array[num_atoms].AddTo_Atom_Location(coors);
	    }
	 else
	  //Atom placement.
	    {
	    atom_array[num_atoms].Set_Atom_Location(coors);
	    }
	 ++num_atoms;
     }

atom_collection& atom_collection::operator ++ ()
    {
    ++num_atoms;         
	return *this;
    }

void atom_collection::Inc_Coll_Size(int inc)
    {
    num_atoms += inc;
    }

//Group tag handling---

int atom_collection::Get_Group_Tag(int index) const
	{
	return atom_array[index].group_tag;
    }

void atom_collection::Set_Group_Tag(int start, int end, int tag_val)
	{
	for (int a = start; a < end; ++a)
		{
		atom_array[a].Set_Atom_Group_Tag(tag_val);
		}
	}	
	 
int atom_collection::Next_Tag()
	{
	int highest_tag = -1;
	 //-1 should be the lowest possible value for a group tag.
	for (int a = 0; a < num_atoms; ++a)
		{
		if (atom_array[a].group_tag > highest_tag)
		   {
		   highest_tag = atom_array[a].group_tag;
		   }
		}
	return (highest_tag + 1);
	}

int atom_collection::Number_Of_Groups()
	{
	int next_tag = Next_Tag();

	//Get array of booleans, one for each group tag
	//that can exist (i.e. is less than next_tag)---

	bool* tag_finder = new bool[next_tag]; 
	False_Array(tag_finder, next_tag);

	for (int a = 0; a < num_atoms; ++a)
	 //Find which tag values are found on
	 //at least one atom.
		{
		if (atom_array[a].group_tag > -1)
		   {
		   tag_finder[atom_array[a].group_tag] = TRUEV;
		   }
		}

	//Count the number of group tag booleans turned on----

	int total_tags = 0;
	for (int a = 0; a < next_tag; ++a)
		{
		if (tag_finder[a])
		   {
		   ++total_tags;
		   }
		}

	delete[] tag_finder;
	return total_tags;
	}

int atom_collection::Molecule_Size(int mol_index)
	{
	int size = 0;
	for (int a = 0; a < num_atoms; ++a)
		{
		if (atom_array[a].group_tag == mol_index)
		   {
		   ++size;
		   }
		}
	return size;
	}

void atom_collection::Assign_Group_Tags()
	{
	int* bond_network = new int[num_atoms];
	int molecule_count;
	Generate_Bond_Network(bond_network, 0, num_atoms, molecule_count);
	Assign_Group_Tags(bond_network);
	delete[] bond_network;
	}

void atom_collection::Assign_Group_Tags(int* bond_network)
	{
	for (int a = 0; a < num_atoms; ++a)
		{
		atom_array[a].group_tag = bond_network[a];
		}
	}

void atom_collection::Defrag_Group_Tags()
	{
	int next_tag = Next_Tag();

	//Get array of booleans, one for each group tag
	//that can exist (i.e. is less than next_tag)---

	bool* tag_finder = new bool[next_tag]; 
	False_Array(tag_finder, next_tag);
	 //Initialize array to FALSEV.

	for (int a = 0; a < num_atoms; ++a)
	 //Find which tag values are found on
	 //at least one atom.
		{
		if (atom_array[a].group_tag > -1)
		   {
		   tag_finder[atom_array[a].group_tag] = TRUEV;
		   }
		}

	//Determine shifting values---

	int* shift_values = new int[next_tag];
	int cur_shift_value = 0;
	for (int a = 0; a < next_tag; ++a)
		{
		shift_values[a] = cur_shift_value;
		if (tag_finder[a] == FALSEV)
		 //"Hole" found.
			{
			++cur_shift_value;
			}
		}
	
	//Apply shifting values/defrag---

	int tag;
	for (int a = 0; a < num_atoms; ++a)
		{
		tag = atom_array[a].group_tag;
		if (tag > -1)
		   {
		   atom_array[a].group_tag -= shift_values[tag];
		   }
		}

	
	delete[] shift_values;
	delete[] tag_finder;
	}

//Atom removal---

void atom_collection::Delete_Atom(int index)
    {
    for (int a = index; a < (num_atoms - 1); ++a)
		{
		atom_array[a] = atom_array[a + 1];
		}
	--num_atoms;
    }

void atom_collection::Delete_Atoms(int index1, int index2)
    {
    int deletion_range = index2 - index1;
	for (int a = index1; a < (num_atoms - deletion_range); ++a)
		{
		atom_array[a] = atom_array[a + deletion_range];
		}
	num_atoms -= deletion_range;
	}

void atom_collection::Delete_Atoms(const double* corner, const double* box_size)
    {
    int org_size = num_atoms;
	double loc[3];
	for (int a = 0; a < num_atoms; ++a)
 	    {
        atom_array[a].Get_Atom_Location(loc);
		Sub_XYZ(loc, corner);
		if (Check_BoundariesZeroToMax(loc, box_size))
		   {
		   Delete_Atom(a);
		   --a;
		   }
	    }
    atom_array.resize(org_size);
	 //Don't lose memory reservation.
	}

void atom_collection::Make_Orthorhombic_Void(const double* corner, const double* box_size)
    {
    Delete_Atoms(corner, box_size);
    }  

void atom_collection::Make_Spherical_Void(const double* origin, double radius, bool pbc)
    {
    int org_size = num_atoms;
	double dist;
	for (int a = 0; a < num_atoms; ++a)
	    {
		if (pbc)
		   {
		   dist = atom_array[a].Get_Distance(origin, periodicity, box_dimensions);
		   }
		else
		   {
		   dist = atom_array[a].Get_Distance(origin);
		   }
		if (radius > (dist - FP_ERROR_FIX) )
		   {
		   Delete_Atom(a);
		   --a;
		   }
	    }
	atom_array.resize(org_size);
    }

//Molecule removal---

void atom_collection::Delete_Molecule(int group_tag)
	{
	for (int a = 0; a < num_atoms; ++a)
		{
		if (atom_array[a].group_tag == group_tag)
		   {
		   Delete_Atom(a);
		   --a;
		   }
		}
	}

void atom_collection::Delete_Molecule_Containing_Atom(int index)
	{
	int tag = atom_array[index].group_tag;
	Delete_Molecule(tag);
	}

//Periodic box handling---

void atom_collection::Set_Periodicity(int period)
	 {
	 periodicity = period;
     }

void atom_collection::No_Periodicity()
	 {
	 periodicity = 0;
     }
	 
void atom_collection::Set_Box_Vectors(const double* box_vec)
	 {
	 Set_XYZ(box_dimensions, box_vec);
	 }
	 
void atom_collection::Set_Box(int period, const double* box_vec)
	 {
	 periodicity = period;
	 Set_XYZ(box_dimensions, box_vec);
	 }

int atom_collection::Get_Periodicity() const
     {
	 return periodicity;
     }

void atom_collection::Get_Box_Vectors(double* vec) const
     {
	 Set_XYZ(vec, box_dimensions);
     }

//Spatial information---

void atom_collection::Get_Bond_Vec(int ind1, int ind2, double* bond_vec) const
    {
	Get_ClosestVec_OrthoPBC(bond_vec, atom_array[ind1].atom_location, 
		                    atom_array[ind2].atom_location, periodicity, box_dimensions);
    }

void atom_collection::Get_Coordinate_Mins(double* mins) const
     {
	 Get_Coordinate_Mins(0, num_atoms, mins);
     }

void atom_collection::Get_Coordinate_Mins(int start, int end, double* mins) const
	 {
	 double locs[3];
	 Zero_XYZ(mins);
	 if (num_atoms > start)
      //Initial minimum.
	    {
		atom_array[start].Get_Atom_Location(mins);
	    } 
	 for (int a = start + 1; a < end; ++a)
	  //Look through all atoms for the minimum
	  //coordinates.
         {
         atom_array[a].Get_Atom_Location(locs);
		 for (int b = 0; b < 3; ++b)
		  //Compare coordinates, per dimension.
		    {
		    if (locs[b] < mins[b])
		      {
			  mins[b] = locs[b];
		      }
		    }
         }      
	 }

void atom_collection::Get_Coordinate_Maxs(double* maxs) const
     {
     double locs[3];
	 Zero_XYZ(maxs);
	 if (num_atoms > 0)
	  //Initial maximum.
	    {
	    atom_array[0].Get_Atom_Location(maxs);
	    } 
	 for (int a = 1; a < num_atoms; ++a)
	  //Look through all atoms for the minimum
	  //coordinates.
         {
         atom_array[a].Get_Atom_Location(locs);
		 for (int b = 0; b < 3; ++b)
		  //Compare coordinates, per dimension.
		    {
		    if (locs[b] > maxs[b])
		      {
			  maxs[b] = locs[b];
		      }
		    }
         }      
     }

void atom_collection::Coordinate_Average(double* avg_coors) const
     {
	 Coordinate_Average(num_atoms, avg_coors);                             
     } 

void atom_collection::Coordinate_Average(int end_atoms, double* avg_coors) const
     {
	 Coordinate_Average(num_atoms - end_atoms, num_atoms, avg_coors);                             
     }

     
void atom_collection::Coordinate_Average(int start, int end, double* avg_coors) const
     {
	 double temp_coors[3];
     Zero_XYZ(temp_coors);   
     for (int a = start; a < end; ++a)
         {
         Add_XYZ(temp_coors, atom_array[a].atom_location);
         }
     Divide_XYZ(avg_coors, temp_coors, end - start);                            
     } 


void atom_collection::Fit_Sphere(double dist_max, bool AB_only, 
								 int internal_CN, double& sphere_radius,
								 int atomic_number_specific, bool anal_hydrogen) const
    {
	int num_surf_atoms;
	int* surf_atom_indices = new int[num_atoms];
	 //Storage space for the indices of the atoms that are at the surface.
	Get_Surface_Atom_List(surf_atom_indices, num_atoms, num_surf_atoms, 
		                  dist_max, AB_only, internal_CN, anal_hydrogen);
	double effective_sphere_center[3];
	Coordinate_Average(effective_sphere_center);
	 //Center of the sphere is to be placed at the coordinate average.
	 //For any intended uses of this function, the dimensionality of
	 //the system should not be significant.
	double square_distance_sum = 0.0;
	 //Holds the sum of square distances from the coordinate center
	 //taken over all surface atoms.
	int weight = 0;
	 //Number of surface atoms considered in the following analysis.
    double temp_distance;
	for (int a = 0; a < num_surf_atoms; ++a)
	    {
		if (atomic_number_specific != -1)
		 //Allow for an analysis based on atoms of a specific atomic number.
		   {
		   if (atom_array[surf_atom_indices[a]].Get_Atom_Number() 
			   != atomic_number_specific)
		      {
			  continue;
		      }
		   }

		if (atom_array[surf_atom_indices[a]].Is_Invisible(anal_hydrogen))
		 //Ignore dummy atoms and hydrogens if requested.
		   {
		   continue;
		   }

		temp_distance = atom_array[surf_atom_indices[a]].
			            Get_Distance(effective_sphere_center, periodicity, box_dimensions);
		square_distance_sum += temp_distance*temp_distance;
		++weight;
	    }
	sphere_radius = sqrt(square_distance_sum/double(weight));
	delete[] surf_atom_indices;
    }

void atom_collection::Fit_Box(double dist_max, bool AB_only, 
							  int internal_CN, double* box_fit,
							  int atomic_number_specific,
							  bool anal_hydrogen) const
    {
	double min_coors[3];
	Get_Coordinate_Mins(min_coors);
	Get_Coordinate_Maxs(box_fit);
	Sub_XYZ(box_fit, min_coors);
    }

//Collection --> collection copying---

void atom_collection::Copy_Atomic_Group(const atom_collection& atoms_to_copy)
     {
	 Memory_Check(num_atoms + atoms_to_copy.num_atoms);
     for (int a = 0; a < atoms_to_copy.num_atoms; ++a)
       {
       atom_array[num_atoms + a] = atoms_to_copy.atom_array[a];   
       }
	 num_atoms += atoms_to_copy.num_atoms;
     } 

void atom_collection::Copy_Atomic_Group(const atom_collection& atoms_to_copy, int index1, int index2, 
                                        int num_to_copy, const double* posi_shift)
     {
     Memory_Check(num_atoms + num_to_copy);
     for (int a = 0; a < num_to_copy; ++a)
       {
       atom_array[index1 + a] = atoms_to_copy.atom_array[index2 + a];   
       Add_XYZ(atom_array[index1 + a].atom_location, posi_shift);
       }
	 int number_added = num_to_copy - (num_atoms - index1);
	 if (number_added < 0)
      //No change to number of atoms in the atomic collection.
	    {
		number_added = 0;
	    }
	 num_atoms += number_added;
     }

void atom_collection::Copy_Atomic_Group(const atom_collection& atoms_to_copy, int index1, int index2, 
                                        int num_to_copy, const double* posi_shift, const double* rot_angles)
     {
     Memory_Check(num_atoms + num_to_copy);
	 int starting_size = num_atoms;
     for (int a = 0; a < num_to_copy; ++a)
	  //Copy the atoms over and apply translation.
       {
       atom_array[index1 + a] = atoms_to_copy.atom_array[index2 + a];   
       Add_XYZ(atom_array[index1 + a].atom_location, posi_shift);
       }
	 int number_added = num_to_copy - (num_atoms - index1);
	 if (number_added < 0)
      //No change to number of atoms in the atomic collection.
	    {
		number_added = 0;
	    }
	 num_atoms += number_added;
	 double group_center[3];
	 Coordinate_Average(num_atoms - starting_size, group_center);
	  //Record the coordinate average for the rotation operation.
	 Rotate_AtomsXYZ(num_atoms - starting_size, rot_angles, 
		             TRUEV, group_center);
	  //Clockwise rotation of the added atoms.
     }

void atom_collection::Copy_Atomic_Group(const atom_collection& atoms_to_copy, int index1, int index2, 
                                        int num_to_copy, const double* posi_shift, const double* rot_angles,
										bool stop_hydro)
     {
     Memory_Check(num_atoms + num_to_copy);
	 int starting_size = num_atoms;
	 int num_copied = 0;
	 bool include_hydro = !stop_hydro;
     for (int a = 0; a < num_to_copy; ++a)
	  //Copy the atoms over and apply translation.
       {
	   if (!atoms_to_copy.atom_array[index2 + a].Is_Invisible(include_hydro))
		  {
          atom_array[index1 + num_copied] = atoms_to_copy.atom_array[index2 + a];  
          Add_XYZ(atom_array[index1 + num_copied].atom_location, posi_shift);
		  ++num_copied;
		  }
       }
	 int number_added = num_copied - (num_atoms - index1);
	 if (number_added < 0)
      //No change to number of atoms in the atomic collection.
	    {
		number_added = 0;
	    }
	 num_atoms += number_added;
	 double group_center[3];
	 Coordinate_Average(num_atoms - starting_size, group_center);
	  //Record the coordinate average for the rotation operation.
	 Rotate_AtomsXYZ(num_atoms - starting_size, rot_angles, 
		             TRUEV, group_center);
	  //Clockwise rotation of the added atoms.
     }

//Name formatting---

void atom_collection::Index_Similar_Atoms()
       {
	   Remove_Indices();
	    //Atom names must not be terminated with a number for the 
	    //following indexation to work.

	   //Find list of unique atom names--- 

       char** atom_names;
	   Get_Memory(atom_names, MAX_ATOM_TYPES, MAX_ANAME_SIZE);
	    //List of atom names found after removing indices.
       int num_atom_types;
        //Number of atom names in the system.
       Get_Name_List(atom_names, num_atom_types, MAX_ATOM_TYPES);
		
	   //Prepare array of current indices---

	   int* current_indices = new int[num_atom_types];
	   Initialize_Array(current_indices, num_atom_types, 1);
		 //First index = 1.

	   //Index the atom names---

       char* name_pointer;
       for (int a = 0; a < num_atoms; ++a)
           {
           name_pointer = atom_array[a].atom_name;
           for (int b = 0; b < num_atom_types; ++b)
               {
               if (strcmp(name_pointer, atom_names[b]) == 0)
                  //If true, the strings match each other.
                  {
                  Add_Number(name_pointer, current_indices[b]);
                  ++current_indices[b];     
                  b = num_atom_types;
                  }       
               }
           }       

	   //Memory release---

	   delete[] current_indices;
	   Free_Memory(atom_names, MAX_ATOM_TYPES);
       }
       
void atom_collection::Remove_Indices()
       {
       for (int a = 0; a < num_atoms; ++a)
           {
           Remove_Number(atom_array[a].atom_name);
           }                         
       }

void atom_collection::Rename_Atoms(const char* old_name, const char* new_name)
	   {
	   for (int a = 0; a < num_atoms; ++a)
		//Find atoms with the old name and replace them with the new.
		   {
		   if (Name_Check(atom_array[a].atom_name, old_name))
			  {
			  atom_array[a].atom_name[0] = '\0';
			  strcpy(atom_array[a].atom_name, new_name);
			  }
		   }
       }
       
//Charge-based operations---

double atom_collection::Return_Total_Charge() const
     {
	 double total_charge = 0.0;
     for (int a = 0; a < num_atoms; ++a)
	  //Sum up the charges.
	     {
         total_charge += atom_array[a].Get_Atom_Charge();
	     }
	 return total_charge;
     }

double atom_collection::Return_Total_Charge(const double* corner, 
	                                        const double* size) const
     {
	 double total_charge = 0.0;
	 double maxs[3];
	 Set_XYZ(maxs, corner);
	 Add_XYZ(maxs, size);
     for (int a = 0; a < num_atoms; ++a)
	     {
	    if (Check_Boundaries(atom_array[a].atom_location, corner, maxs))
		   {
		   total_charge += atom_array[a].Get_Atom_Charge();
		   }
	     }
	 return total_charge;
     }

double atom_collection::Calculate_ES_Potential(bool include_hydro) const
     {

     //Get array of charged atoms for electrostatic summation---

	 int num_charged_atoms;
	 int* charged_atom_indices = new int[num_atoms];
	 Get_Charged_Atom_List(charged_atom_indices, num_atoms, num_charged_atoms,
		                   include_hydro);

	 //Electrostatic summation---

	 double total_potential = 0.0;
	 double distance_temp, charge1, charge2;
	 int index1_temp, index2_temp;
	 for (int a = 0; a < num_charged_atoms; ++a)
	     {
		 index1_temp = charged_atom_indices[a];
		 charge1 = atom_array[index1_temp].Get_Atom_Charge();
		 for (int b = a + 1; b < num_charged_atoms; ++b)
		     {
			 index2_temp = charged_atom_indices[b];
			 charge2 = atom_array[index2_temp].Get_Atom_Charge();
			 distance_temp = atom_array[index1_temp].Get_Distance(atom_array[index2_temp], 
																  periodicity, box_dimensions);
			 total_potential += charge1*charge2/distance_temp;
		     }
	     }

	 delete[] charged_atom_indices;

	 return (ELEC_POT_CONSTANT_EV*total_potential);
     }

//Translational operations---

void atom_collection::Shift_Atomic_Coordinates(const double* shift)
     {
	 Shift_Atomic_Coordinates(0, num_atoms, shift);                   
     }

void atom_collection::Shift_Atomic_Coordinates(int end_atoms, const double* shift)
     {
	 Shift_Atomic_Coordinates(num_atoms - end_atoms, num_atoms, shift);                            
     }
     
void atom_collection::Shift_Atomic_Coordinates(int start, int end, const double* shift)
     {
	 for (int a = start; a < end; ++a)
         {
         atom_array[a].AddTo_Atom_Location(shift);     
         }     
     }

void atom_collection::Inverse_Shift_Atomic_Coordinates(const double* shift)
     {
	 Inverse_Shift_Atomic_Coordinates(0, num_atoms, shift);                   
     }

void atom_collection::Inverse_Shift_Atomic_Coordinates(int end_atoms, const double* shift)
     {
	 Inverse_Shift_Atomic_Coordinates(num_atoms - end_atoms, num_atoms, shift);                            
     }
     
void atom_collection::Inverse_Shift_Atomic_Coordinates(int start, int end, const double* shift)
     {
	 double neg_shift[3];
	 Zero_XYZ(neg_shift);
	 Sub_XYZ(neg_shift, shift);
	 for (int a = start; a < end; ++a)
         {
         atom_array[a].AddTo_Atom_Location(neg_shift);     
         }     
     }

void atom_collection::Translate_To_New_Coordinate_Mins(const double* new_min)
     {
	 Translate_To_New_Coordinate_Mins(0, num_atoms, new_min);
     }

void atom_collection::Translate_To_New_Coordinate_Mins(int end_atoms, const double* new_min)
	 {
	 Translate_To_New_Coordinate_Average(num_atoms - end_atoms, num_atoms, new_min);
     }

void atom_collection::Translate_To_New_Coordinate_Mins(int start, int end, const double* new_min)
	 {
	  double min_vals[3];
	  Get_Coordinate_Mins(start, end, min_vals);
	   //Get current minimum coordinates.
	  double shift_vals[3];
	  SetSub_XYZ(shift_vals, new_min, min_vals);
	   //Determine shift to get the new set of coordinate minimums.
	  Shift_Atomic_Coordinates(start, end, shift_vals);
	   //Perform the shift.
     }

void atom_collection::Translate_To_New_Coordinate_Average(const double* new_avg_coors)
     {
	 Translate_To_New_Coordinate_Average(0, num_atoms, new_avg_coors);
     }

void atom_collection::Translate_To_New_Coordinate_Average(int end_atoms,
	                                                      const double* new_avg_coors)
     {
	 Translate_To_New_Coordinate_Average(num_atoms - end_atoms, num_atoms, new_avg_coors);
     }

void atom_collection::Translate_To_New_Coordinate_Average(int start, int end,
														  const double* new_avg_coors)
     {
     double old_avg_coors[3];
     Coordinate_Average(start, end, old_avg_coors); 
	  //Get coordinate average to determine how far it is from the desired result.
         
     double shift_coors[3];
	 Set_XYZ(shift_coors, new_avg_coors);
     Sub_XYZ(shift_coors, old_avg_coors);
	  //Get the translational shift.

     Shift_Atomic_Coordinates(start, end, shift_coors);  
	  //Peform the shift.
     }

void atom_collection::Move_All_Positions_Into_Box()
	 {
	 for (int a = 0; a < num_atoms; ++a)
	     {
		 Get_PositiveVec_OrthoPBC(atom_array[a].atom_location, periodicity, box_dimensions);
	     }
	 }

//Rotational operations (vec 1 --> vec 2)---

void atom_collection::Rotate_Atoms(const double* start_vec, 
								   const double* final_vec)
    {
	Rotate_Atoms(0, num_atoms, start_vec, final_vec);
    }

void atom_collection::Rotate_Atoms(int atoms_to_rotate, const double* start_vec, 
                                   const double* final_vec)
    {
	Rotate_Atoms(num_atoms - atoms_to_rotate, num_atoms, start_vec, final_vec);                    
    }

void atom_collection::Rotate_Atoms(int start, int end, const double* start_vec, 
                                   const double* final_vec)
    {
	double rot_matrix[9];
	Calc_Rotation_Matrix(start_vec, final_vec, rot_matrix);
    for (int a = start; a < end; ++a)
       {
       Rotate_Vector(rot_matrix, atom_array[a].atom_location);   
       }                           
    }

void atom_collection::Rotate_Atoms(const double* start_vec, 
                        const double* final_vec, const double* origin)
    {
	Rotate_Atoms(0, num_atoms, start_vec, final_vec, origin);                     
    }

void atom_collection::Rotate_Atoms(int atoms_to_rotate, const double* start_vec, 
                                   const double* final_vec, const double* origin)
    {
    Rotate_Atoms(num_atoms - atoms_to_rotate, num_atoms, start_vec, final_vec, origin);                   
    }

void atom_collection::Rotate_Atoms(int start, int end, const double* start_vec,
	                               const double* final_vec, const double* origin)
	{
	Inverse_Shift_Atomic_Coordinates(origin);
    Rotate_Atoms(start, end, start_vec, final_vec); 
    Shift_Atomic_Coordinates(origin);
	}

//Rotational operations (vector and angle)---

void atom_collection::Rotate_Atoms(const double* rotation_vec, double rotation_angle)
    {
	Rotate_Atoms(0, num_atoms, rotation_vec, rotation_angle);                        
    }

void atom_collection::Rotate_Atoms(int num_to_rotate, const double* rotation_vec, 
	                               double rotation_angle)
    {
	Rotate_Atoms(num_atoms - num_to_rotate, num_atoms, rotation_vec, rotation_angle);                    
    }

void atom_collection::Rotate_Atoms(int start, int end, 
	                               const double* rotation_vec, double rotation_angle)
	{
    double rotation_matrix[9];
    Calc_Rotation_Matrix(rotation_vec, rotation_angle, rotation_matrix);
    for (int a = start; a < end; ++a)
       {
       Rotate_Vector(rotation_matrix, atom_array[a].atom_location);   
       }  
	}

void atom_collection::Rotate_Atoms(const double* rotation_vec, double rotation_angle, const double* origin)
    {
	Rotate_Atoms(0, num_atoms, rotation_vec, rotation_angle, origin);                    
    }

void atom_collection::Rotate_Atoms(int num_to_rotate, const double* rotation_vec, 
	                               double rotation_angle, const double* origin)
    {
	Rotate_Atoms(num_atoms - num_to_rotate, num_atoms, rotation_vec, 
				 rotation_angle, origin);                        
    }

void atom_collection::Rotate_Atoms(int start, int end, const double* rotation_vec,
	                               double rotation_angle, const double* origin)
	{
	Inverse_Shift_Atomic_Coordinates(origin);
    Rotate_Atoms(start, end, rotation_vec, rotation_angle); 
    Shift_Atomic_Coordinates(origin);      
	}

//Rotational operations (set of axis rotation-based angles)---

void atom_collection::Rotate_AtomsXYZ(const double* angles, bool clockwise)
    {
	Rotate_AtomsXYZ(0, num_atoms, angles, clockwise); 
    }

void atom_collection::Rotate_AtomsXYZ(int num_to_rotate, const double* angles, 
	                                  bool clockwise)
    {
	Rotate_AtomsXYZ(num_atoms - num_to_rotate, num_atoms, angles, clockwise);
    }

void atom_collection::Rotate_AtomsXYZ(int start, int end, const double* angles, 
	                                  bool clockwise)
    {
	double rangles[3];
	Set_XYZ(rangles, angles);
	if (!clockwise)
	   {
	   Neg_XYZ(rangles);
	   }
	double rotation_matrix[9];
	Calc_Rotation_Matrix(rangles[0], rangles[1], rangles[2], rotation_matrix);
	for (int a = start; a < end; ++a)
       {
       Rotate_Vector(rotation_matrix, atom_array[a].atom_location);   
       }       
    }

void atom_collection::Rotate_AtomsXYZ(const double* angles,
	                                  bool clockwise, const double* origin)
    {
	Rotate_AtomsXYZ(0, num_atoms, angles, clockwise, origin); 
    }

void atom_collection::Rotate_AtomsXYZ(int num_to_rotate, const double* angles, 
	                                   bool clockwise, const double* origin)
    {
	Rotate_AtomsXYZ(num_atoms - num_to_rotate, num_atoms, angles, 
		            clockwise, origin); 
    }

void atom_collection::Rotate_AtomsXYZ(int start, int end, const double* angles, 
	                                  bool clockwise, const double* origin)
    {
	Inverse_Shift_Atomic_Coordinates(origin);
    Rotate_AtomsXYZ(start, end, angles, clockwise); 
    Shift_Atomic_Coordinates(origin);    
    }

//Atom searching---
     
int atom_collection::Find_Atom_Name(const char* name_string, bool is_symbol) const
    {
    int name_index = Find_Atom_Name(0, num_atoms, name_string, is_symbol);
	return name_index;
    }

int atom_collection::Find_Atom_Name(int min_index, int max_index, 
	                                const char* name_string, bool is_symbol) const
    {
    int index = -1; 

	//Search for an exact match between passed string and atom names---

    for (int a = min_index; a < max_index; ++a)
        {
        if (Name_Check(name_string, atom_array[a].atom_name))
           {
           index = a;  
           a = max_index;
           }
        }      

	//Check (if previous search failed) for atom names that differ
	//from the passed string by just an index at the string end
	//E.g. atom name "Au3" = search string "Au5"---

	if (index == -1)
	 //Check for atom name without regard for numbers at string ends.
	   {
	   for (int a = min_index; a < max_index; ++a)
          {
          if (Name_Check_To_Index(name_string, atom_array[a].atom_name))
            {
            index = a;  
            a = max_index;
            }
          }    
	   }

	//If atomic symbol is set to TRUEV, look for a capitalization difference.
	//E.g. "ZN" vs. "Zn." This is useful particularly for PDB file analysis---

	if ( (index == -1) && is_symbol)
	   { //cout << endl << "HERE WITH: " << name_string;
	   bool two_capitals = Contains_Two_Capitals_At_Start(name_string);
	    //E.g. "ZN", "MG" , etc. 
	   if (two_capitals)
	     { 
		 char temp_string[MAX_ANAME_SIZE];
		 temp_string[0] = name_string[0];
		 temp_string[1] = Make_Undercase(name_string[1]);
		 temp_string[2] = '\0';
		 for (int a = min_index; a < max_index; ++a)
		   {
           if (Name_Check_To_Index(temp_string, atom_array[a].atom_name))
             {
             index = a;  
             a = max_index;
             }
		   }
		  }
        } 

    return index;  
    }

int atom_collection::Find_Closest_Atom(const double* coors) const
       {
       int atom_index;
       atom_index = Find_Closest_Atom(0, num_atoms, coors);
       return atom_index;
       }
 
int atom_collection::Find_Closest_Atom(int min_index, int max_index, 
	                                   const double* coors) const
       {
	   int closest_atom_index = min_index;
	   double least_distance = atom_array[min_index].Get_Distance(
		                       coors, periodicity, box_dimensions);
	   	 //Initialize test to closest atom = first atom investigated.
       double dist_val;
       for (int a = min_index + 1; a < max_index; ++a)
           {
		   dist_val = atom_array[a].Get_Distance(coors, periodicity, box_dimensions);
           if (dist_val < least_distance)
              {
              least_distance = dist_val;
              closest_atom_index = a;    
              }
           }     
       return closest_atom_index;
       }

//Category lists---

void atom_collection::Get_Name_List(char** name_array, 
                                    int& name_index, 
									const int MAX_NAME_COUNT) const
    {
    name_index = 0;
    bool found_new_name;
    const char* temp_string;
    for (int a = 0; a < num_atoms; ++a)
        {
        temp_string = atom_array[a].atom_name;
        found_new_name = TRUEV;
        for (int b = 0; b < name_index; ++b)
		 //Look through name list to see if the atom name has
		 //been found previously.
            {
            if (strcmp(temp_string, name_array[b]) == 0)  
               {
               found_new_name = FALSEV;
               b = name_index;                        
               }   
            }
        if (found_new_name)
         //Found a new atom name!
           {
           strcpy(name_array[name_index], temp_string);
           ++name_index;
           if (name_index == MAX_NAME_COUNT)
              {
			  Show_Warning("TOO MANY ATOM TYPES IN SYSTEM FOR ATOM NAME LIST CREATION!");
              a = num_atoms;
              }
           }     
        }
    }

void atom_collection::Get_Index_List(char** atom_names, int* reference_array,
                                     int num_names) const
    {
    for (int a = 0; a < num_names; ++a)
     //Find an example atom for every atom name of interest.
       {  
       reference_array[a] = Find_Atom_Name(atom_names[a], FALSEV);    
       if (reference_array[a] == -1)
          {
		  Show_Warning("DEFINE ALL ATOMS! CAN NOT FIND ", atom_names[a], " ATOM!");                    
          }                           
       }                     
    }

void atom_collection::Get_Index_List(int* reference_array, int& num_types) const
    {
	char** atom_names;
	Get_Memory(atom_names, MAX_ATOM_TYPES, MAX_ANAME_SIZE);
	 //List of atom names, one per atom type.
	
	Get_Name_List(atom_names, num_types, MAX_ATOM_TYPES);
	Get_Index_List(atom_names, reference_array, num_types);
	 //Find the name list and then find atoms to be 
	 //associated with those names.

	Free_Memory(atom_names, MAX_ATOM_TYPES);              
    }

void atom_collection::Get_Surface_Atom_List(int* index_list, 
	                  const int MAX_INDEX_NUM, int& num_surf_atoms, 
					  double dist_max, bool AB_only, int internal_CN,
					  bool include_hydro) const
	{
	num_surf_atoms = 0;
	for (int a = 0; a < num_atoms; ++a)
		{
		if (atom_array[a].Is_Invisible(include_hydro))
		 //Ignore invisible atoms.
		   {
		   continue;
		   }

		if (Is_Surface_Atom(a, dist_max, AB_only, internal_CN))
		 //Test atom for a surface atom coordination environment.
			{
			if (num_surf_atoms == MAX_INDEX_NUM)
			   {
			   Show_Warning("NOT ENOUGH MEMORY IN SURFACE ATOM INDEX ARRAY!");
			   break;
			   }

			index_list[num_surf_atoms] = a;
			++num_surf_atoms;
			}
		}
	}

void atom_collection::Get_Charged_Atom_List(int* index_list, 
	                  const int MAX_INDEX_NUM, int& num_charged_atoms, 
					  bool include_hydro) const
	{
	num_charged_atoms = 0;
	for (int a = 0; a < num_atoms; ++a)
		{
		if (atom_array[a].Is_Invisible(include_hydro))
		 //Ignore invisible atoms.
		   {
		   continue;
		   }
		if (!Check_FP_Equality(atom_array[a].atomic_charge, 0.0))
			{
			if (num_charged_atoms == MAX_INDEX_NUM)
		     //Check for maximum array capacity violation.
			   {
			   Show_Warning("NOT ENOUGH MEMORY IN CHARGED ATOM INDEX ARRAY!");
			   break;
			   }

			index_list[num_charged_atoms] = a;
			++num_charged_atoms;
			}
		}
	}

void atom_collection::Get_Molecule_List(int* index_list, const int MAX_INDEX_NUM, int molecule_index, 
	                                    int& number_atoms, bool include_hydro) const
	{
	number_atoms = 0;
	for (int a = 0; a < num_atoms; ++a)
		{
		if (atom_array[a].Is_Invisible(include_hydro))
		 //Ignore invisible atoms.
		   {
		   continue;
		   }
		if (atom_array[a].group_tag == molecule_index)
			{
			if (number_atoms == MAX_INDEX_NUM)
		     //Check for maximum array capacity violation.
			   {
			   Show_Warning("NOT ENOUGH MEMORY IN MOLECULE ATOM INDEX ARRAY!");
			   break;
			   }

			index_list[number_atoms] = a;
			++number_atoms;
			}
		}	
	}

//Atom counting---

void atom_collection::Get_Composition(char* ele_comp_string, const int MAX_STRING_LENGTH) const
    {
	char** names;
	int name_count;
	Get_Memory(names, MAX_ATOM_TYPES, MAX_ANAME_SIZE);
	Get_Name_List(names, name_count, MAX_ATOM_TYPES);
	 //Obtain a list of the unique atom types in the collection.

	int char_index = 0;
	 //Position in the composition string.
	int atom_count;
	for (int a = 0; a < name_count; ++a)
	    {
		//Add name of atom type---

		strcpy(&(ele_comp_string[char_index]), names[a]);
		char_index += strlen(names[a]);
		
		//Name-number spacing character, e.g. "Si_300"---

		ele_comp_string[char_index] = '_';
		++char_index;

		//Addition of atom count to string---

		ele_comp_string[char_index] = '\0';
		atom_count = Count_Atom_Type(names[a]);
		Add_Number(ele_comp_string, atom_count);
		char_index += String_Size(atom_count);
		
		//Number-name spacing character, e.g. "Si_300,N_400"---

		ele_comp_string[char_index] = ',';
		++char_index;

		//Memory-violation check---

		if ( (char_index > (MAX_STRING_LENGTH - MAX_ANAME_SIZE - 10) ) )
			{
		    Show_Warning("NOT ENOUGH STRING SPACE FOR ELEMENTAL COMPOSITION!");
			break;
			}

	    }

	ele_comp_string[char_index - 1] = '\0';
	 //Final string-terminating character.
	Free_Memory(names, MAX_ATOM_TYPES);
    }

int atom_collection::Count_Atom_Type(const char* atom_name) const
	{
	int count = 0;
	for (int a = 0; a < num_atoms; ++a)
	 //Look for atoms with the passed atom name.
		{
		if (Name_Check(atom_name, atom_array[a].atom_name))
		   {
		   ++count;
		   }
		} 
	return count;
	}

int atom_collection::Count_Atom_Type(const char* atom_name, int* b_net, 
	                                 int index) const
	{
	int count = 0;
	for (int a = 0; a < num_atoms; ++a)
	 //Look for atoms with the passed atom name.
		{
	    if (b_net[a] != index)
		 //Skip atoms that are not in the passed bonding matrix.
		   {
		   continue;
		   }
		if (Name_Check(atom_name, atom_array[a].atom_name))
		   {
		   ++count;
		   }
		} 
	return count;
	}


int atom_collection::Count_Unique_Atom_Types() const
    {
    int num_types = 0;

	bool* test_array = new bool[num_atoms];
     //Function uses a series of booleans that mirror the atomic array
	 //in order to assess which atoms have been accounted for
	 //in the atom type counter.
	False_Array(test_array, num_atoms);

	for (int a = 0; a < num_atoms; ++a)
	 //Move through the array, looking for unencountered
	 //atom types, as indicated by an on-switch (= FALSEV).
	    {
	    if (test_array[a] == FALSEV)
		 //If FALSEV, the atom type is a new one.
		   {
		   for (int b = a; b < num_atoms; ++b)
			//Switch off all booleans corresponding to atoms
			//with the same atom name as the one just discovered.
		       {
			   if (strcmp(atom_array[a].atom_name, atom_array[b].atom_name) == 0)
			      {
				  test_array[b] = TRUEV;
			      }
		       }
           ++num_types;
		   }
	    }

	delete[] test_array;
	return num_types;
    }

//Coordination analysis---

int atom_collection::Get_Coordination_Number(int index, 
	                 double dist_max, bool AB_only) const
	//Returns number of neighbor atoms within a cut-off distance 
	//from an atom.
    {
	int coordination_count = 
	    Get_Coordination_Number(index, 0.0, dist_max, AB_only);
	return coordination_count;
    }

int atom_collection::Get_Coordination_Number(int index, 
	                 double dist_min, double dist_max, bool AB_only) const
	{
	int coordination_count = 0;
	double temp_distance;
	int atomic_num1, atomic_num2;
    atomic_num1 = atom_array[index].Get_Atom_Number();
    for (int a = 0; a < num_atoms; ++a)
      { 
	  if ( (a == index) || atom_array[a].Is_Invisible(TRUEV) )
		//Check for dummy atoms.
	     {
		 continue;
	     }

	  atomic_num2 = atom_array[a].Get_Atom_Number();
	  if (AB_only && (atomic_num1 == atomic_num2) )
		//If AB_only is set, like atom pairs can
		//not be considered further.
		 {
		 continue;
	     }

      temp_distance = atom_array[index].Get_Distance(atom_array[a], periodicity, box_dimensions);
	  if ( (temp_distance < dist_max) && (temp_distance > dist_min) )
		 {
		 ++coordination_count;
		 }
	  } 
	
	return coordination_count;
    }

bool atom_collection::Is_Surface_Atom(int index, double dist_max, 
	                  bool AB_only, int internal_CN) const
    {
	bool is_surf_atom = FALSEV;
	if (!atom_array[index].Is_Invisible(TRUEV))
	 //Don't analyze dummy atoms.
	   {
	   int CN = Get_Coordination_Number(index, dist_max, AB_only);
	   is_surf_atom = (CN < internal_CN);
	   }	 
	return is_surf_atom;
    }

void atom_collection::Coordex_Atoms(int start, int end, double dist_min, double dist_max)
    {
    int cn_temp;
    for (int a = start; a < end; ++a)
	 //For each atom, grab the coordination number and append it
	 //to the atom name.
	    {
		cn_temp = Get_Coordination_Number(a, dist_min, dist_max, FALSEV);
		atom_array[a].Index_Atom_Name(cn_temp);
	    }
    }

double atom_collection::Calculate_Surface_Fraction(double dist_max, 
	                    bool AB_only, int internal_CN)
	{
	int total_count = Real_Size(TRUEV);
	 //Get total number of atoms (surface and bulk).
	int surface_count = 0;
	for (int a = 0; a < num_atoms; ++a)
	  {
	  if (!atom_array[a].Is_Invisible(TRUEV))
	   //Don't analyze dummy atoms.
	     {
	     int CN = Get_Coordination_Number(a, dist_max, AB_only);
	     if (CN < internal_CN)
		  //Found a surface atom!
			{
			++surface_count;
			}
	     }	 
	  }
	return (double(surface_count)/double(total_count));
	}
		
//Large-scale spatial analysis---

void atom_collection::Calculate_RDF(const char* atom_name1, const char* atom_name2, 
                                    double* RDF_array, double step_size) const
     {
     const double START_DIST = 0.0;
     const double END_DIST = 1.0;    
	  //RDF spans from 0 to 1 nm along the distance axis.
     int RDF_max_index = int(END_DIST/step_size) + 1; 
	  //Number of data points in the RDF to be calculated.
	 Zero_Array(RDF_array, RDF_max_index);
      //Make sure RDF array is initialized to zero.                
     
     const char* temp_stringA;
	 const char* temp_stringB;
	   //Atom name holders.
     bool name_test1, name_test2;
	   //For checking if the correct atoms are being considered.
     const double FOUR_PI = 4.0 * PI_CONST;
     double temp_distance, RDF_denominator, effective_data_val;
     int step_point;
	  //Data point reference index.
     int total_pair_count = 0;
	  //Number of pairs found in a certain distance range.
     for (int a = 0; a < num_atoms; ++a)
         {
         temp_stringA = atom_array[a].atom_name;
		 name_test1 = ( Name_Check(temp_stringA, atom_name1) || 
						Name_Check(temp_stringA, atom_name2) );
		 if (!name_test1)
		   //If atom is not one of the atoms in the pair being considered,
		   //don't bother analyzing it.
			{
			continue;
			}
         for (int b = a + 1; b < num_atoms; ++b)
             {

			 //Determine if the a-b pair is relevant to the RDF being calculated---

             temp_stringB = atom_array[b].atom_name;
             name_test1 = ( Name_Check(temp_stringA, atom_name1) && 
							Name_Check(temp_stringB, atom_name2) );
             name_test2 = ( Name_Check(temp_stringA, atom_name2) && 
							Name_Check(temp_stringB, atom_name1) );
             if (!(name_test1 || name_test2) )
              //If we are not looking at the pair of interest, ignore it.
                {
                continue;
                }
                
             temp_distance = atom_array[a].Get_Distance(atom_array[b], 
                                           periodicity, box_dimensions);
             if ( (temp_distance < START_DIST) || (temp_distance > END_DIST) )
              //Must consider the limited range.
                {
                continue;                 
                } 
         
			 //At this point, a pair interaction of interest has been found---

             ++total_pair_count;
             step_point = int(temp_distance/step_size);
              //Determine the data point (distance value) that this pair
			  //corresponds to.

			 //Calculate the contribution to the RDF---

			 effective_data_val = (double(step_point) + 0.5) * step_size;
             RDF_denominator = FOUR_PI * effective_data_val * effective_data_val * step_size;
             RDF_array[step_point] += 1.0/RDF_denominator;
             }     
         }
     
	 //Normalize the radial distribution function---

     double volume_considered = ( (FOUR_PI/3.0) * END_DIST * END_DIST * END_DIST)
                              - ( (FOUR_PI/3.0) * START_DIST * START_DIST * START_DIST);
     double normalization_factor = double(total_pair_count)/(volume_considered);
	  //Normalization factor = Average density over "whole" system.
     for (int a = 0; a < RDF_max_index; ++a)
         {
         RDF_array[a] /= normalization_factor;    
         }
     }

int atom_collection::Pair_Count(double min, double max, 
							    bool unlike_only, bool like_only,
								bool include_hydro) const
     {
	 int pair_count = 0;
	 double temp_distance;
	 int atomic_num1, atomic_num2;
	 for (int a = 0; a < num_atoms; ++a)
         {
		 if (atom_array[a].Is_Invisible(include_hydro))
		  //Check for dummy atoms and unwanted hydrogens.
		    {
			continue;
		    }
		 atomic_num1 = atom_array[a].Get_Atom_Number();
         for (int b = (a + 1); b < num_atoms; ++b)
		     {	
			 if (atom_array[b].Is_Invisible(include_hydro))
			  //Check other atom.
				 {
				 continue;
				 }

			 //Check for unlike only/like only pair requests.

			 atomic_num2 = atom_array[b].Get_Atom_Number();
			 if ( unlike_only && (atomic_num1 == atomic_num2) )
			    {
				continue;
			    }
			 if ( like_only && (atomic_num1 != atomic_num2) )
			    {
				continue;
			    }

			 //Add up the pair interaction if it is within
			 //the requested distance range---

             temp_distance = atom_array[a].Get_Distance(atom_array[b], 
                                           periodicity, box_dimensions);
			 if ( (temp_distance > min) && (temp_distance < max) )
			    {
				++pair_count;
			    }
		     }
	      }

	 return pair_count;
     }

double atom_collection::Distance_Average(double min, double max, 
										 bool unlike_only, bool like_only, 
										 bool include_hydro) const
     {
	 int pair_count = 0;
	 double distance_sum = 0.0;
	 double temp_distance;
	 int atomic_num1, atomic_num2;
	 for (int a = 0; a < num_atoms; ++a)
         {
		 if (atom_array[a].Is_Invisible(include_hydro))
		  //Check for dummy atoms and unwanted hydrogens.
		    {
			continue;
		    }

		 atomic_num1 = atom_array[a].Get_Atom_Number();
         for (int b = (a + 1); b < num_atoms; ++b)
		     {
			 if (atom_array[b].Is_Invisible(include_hydro))
			  //Check other atom.
				{
				continue;
				}

			 //Check for unlike-only/like-only interaction requests.

			 atomic_num2 = atom_array[b].Get_Atom_Number();
			 if (unlike_only && (atomic_num1 == atomic_num2) )
			    {
				continue;
			    }
			 if (like_only && (atomic_num1 != atomic_num2)	)
			    {
				continue;
			    }
			 if (unlike_only && like_only)
			    {
				Show_Congratulations("YOU HAVE UNLOCKED SECRET MESSAGE #1");
			    } 

			 //Add up the pair interaction if it is within
			 //the requested distance range---

             temp_distance = atom_array[a].Get_Distance(atom_array[b], 
                                           periodicity, box_dimensions);
			 if ( (temp_distance > min) && (temp_distance < max) )
			    {
				++pair_count;
				distance_sum += temp_distance;
			    }
		     }
	     }
	 return (distance_sum/double(pair_count));
     }

double atom_collection::Coordination_Average(double min, double max,
										     bool unlike_only, bool like_only,
										     bool include_hydro) const
	 {
	 int pair_count = 2*Pair_Count(min, max, unlike_only, like_only, include_hydro);
	  //Each pair counts as one coordination for EACH atom. Hence,
	  //the multiplying by 2.
	 int atoms_in_pair_count = Real_Size(include_hydro);
	  //Determine the number of atoms involved in the just performed
	  //pair count. Real_Size(...) prevents invisible atoms from
	  //being considered in the average.
	 return (double(pair_count)/double(atoms_in_pair_count));
	 }

void atom_collection::Coordination_Count(double min, double max, 
										 bool unlike_only, bool like_only, 
										 bool include_hydro, int* coordination_counts, 
										 const int MAX_COOR) const
     {
	 int pair_count_sum;
	 double temp_distance;
	 int atomic_num1, atomic_num2;

	 Zero_Array(coordination_counts, MAX_COOR);
	  //Initialize the coordination count array to zero.

	 for (int a = 0; a < num_atoms; ++a)
         {
		 if (atom_array[a].Is_Invisible(include_hydro))
		  //Check for dummy atoms and unwanted hydrogens.
		    {
			continue;
		    }

		 atomic_num1 = atom_array[a].Get_Atom_Number();
		 pair_count_sum = 0;
         for (int b = 0; b < num_atoms; ++b)
		  //Get the coordination number for each atom and then add 1
		  //to the corresponding coordination number count.
		     {
			 if (b == a) continue;
			  //Ignore self-interaction.

			 if (atom_array[b].Is_Invisible(include_hydro))
			  //Check other atom.
				{
				continue;
				}

			 //Like-only/unlike-only request consideration---

			 atomic_num2 = atom_array[b].Get_Atom_Number();
			 if (unlike_only && (atomic_num1 == atomic_num2) )
			    {
				continue;
			    }
			 if (like_only && (atomic_num1 != atomic_num2)	)
			    {
				continue;
			    }

			 //If pair interaction is in distance range, increase
			 //the atom's coordination number---

             temp_distance = atom_array[a].Get_Distance(atom_array[b], 
                                           periodicity, box_dimensions);
			 if ( (temp_distance > min) && (temp_distance < max) )
			    {
				++pair_count_sum;
			    }
		     }

		 if (pair_count_sum < MAX_COOR)
			{
		    ++coordination_counts[pair_count_sum];
			}
		 else
			{
			Show_Warning("COORDINATION NUMBERS GREATER THAN YOUR MAXIMUM VALUE FOUND!");
			}
	     }
     }

//X-ray diffraction calculation---

void atom_collection::Get_Debye_Intensity(double* xrd_data, double min_angle, 
	                                      double step_size, int num_steps, 
										  double wavelength, const double** scat_fit,
										  double dw_factor, bool surf_effect,
										  double surface_dw_enhance, 
										  double coor_cutoff, int coor_num) const
     {
	 double dist, as_factor, scat_factor, angle;
	 double s_value, neg_s_square, sin_angle, dw_atten, dw_effect, surf_atten;
	  //Various variables involved in the calculation of peak intensity.
	 int an_1, an_2;
	  //Atomic number temps.
	 double pre_factor = 4.0*PI_CONST/wavelength;
	 double pre_factor_sin;
	  //Factors in the Debye formula.

	 Zero_Array(xrd_data, num_steps);
	  //Initialize XRD pattern to zero intensity at all angles.

	 //Get a boolean array of surface atoms in the system if needed.
	 //This array essentially runs alongside the atom collection array,
	 //indicating which atoms are surface atoms---

	 int* surface_atoms;
	 bool* surface_list;
	 int surface_atom_count;
	 if (surf_effect)
		{
		surface_atoms = new int[num_atoms];
		Get_Surface_Atom_List(surface_atoms, num_atoms, surface_atom_count, 
		                      coor_cutoff, FALSEV, coor_num, TRUEV);
		 //Obtain a list of atom indices that correspond to surface atoms.

		surface_list = new bool[num_atoms];
		False_Array(surface_list, num_atoms);
		for (int a = 0; a < surface_atom_count; ++a)
	     //Set all surface atoms in the list array to TRUEV.
		 //Remaining atoms are FALSEV, and this array can be used
		 //to quickly assess which atoms are surface atoms.
		    {
		    surface_list[surface_atoms[a]] = TRUEV;
		    }
		
		delete[] surface_atoms;
		 //Don't need this anymore.
	    }

	 //Carry out calculation of pattern, performing it in the 
	 //loop order angle-atom-atom, not atom-atom-angle, as to speed
	 //up the calculation--

	 double* scat_factors;
	 scat_factors = new double[MAX_ATOMIC_NUMBER];
	  //Atomic scattering coefficient array.
	 
	 angle = min_angle;
	 for (int c = 0; c < num_steps; ++c)
	  //Do sampling, angle by angle.
	   {
	   //Determine values for a specific angle and wavelength to 
	   //be used in the Debye formula---

	   angle += step_size;
	   sin_angle = sin(angle); 
	   pre_factor_sin = pre_factor * sin_angle;

	   s_value = sin_angle/wavelength;
	    //s = sin(theta)/lambda.
	   neg_s_square = -1.0*s_value*s_value;
	   dw_atten = exp(dw_factor*neg_s_square);
	     //Debye-Waller attenuation of signal = exp(-2*DWF*s^2).
	     //Per atom value is given here.
	   surf_atten = exp(dw_factor*surface_dw_enhance*neg_s_square);
	     //Surface-enhanced attenuation of signal, if needed.
	   
	   for (int d = 0; d < MAX_ATOMIC_NUMBER; ++d)
	    //Determine the atomic scattering factors at this angle.
		//Fit = four Gaussians of the form: A*exp(-B*s^2)
	       {
	       scat_factors[d] = scat_fit[d][0]*exp(scat_fit[d][1]*neg_s_square)
					       + scat_fit[d][2]*exp(scat_fit[d][3]*neg_s_square)
						   + scat_fit[d][4]*exp(scat_fit[d][5]*neg_s_square)
						   + scat_fit[d][6]*exp(scat_fit[d][7]*neg_s_square);
	       }
	   for (int a = 0; a < num_atoms; ++a)
	     {	
		 an_1 = atom_array[a].Get_Atom_Number();
		 if (Check_FP_Equality(scat_fit[an_1][0], 0.0))
	      //Check to make sure scattering factor is defined for the atom.
			{
			Show_Warning("SCATTERING FACTOR IS UNDEFINED OR ZERO!");
			}
		 for (int b = a; b < num_atoms; ++b)
		  //For each atom pair, determine the contribution to the XRD pattern.
		     {
			 dist = atom_array[a].Get_Distance(atom_array[b], periodicity, box_dimensions);
			 scat_factor = pre_factor_sin*dist;
			  //Term in the Debye formula: (4*PI*sin(angle)*r/wavelength)

			 dw_effect = dw_atten*dw_atten;
			  //Temperature attenutation of the signal.
			 if (surf_effect)
				  //Surface enhancement of DW factors.
				 {
				 if (surface_list[a] && surface_list[b])
				  //2 surface atoms in interaction.
				    {
				    dw_effect = surf_atten*surf_atten;
				    }
				 else if (surface_list[a] || surface_list[b])
				  //1 surface atom in interaction.
				    {
				    dw_effect = surf_atten*dw_atten;
				    }
				 }

			  an_2 = atom_array[b].Get_Atom_Number();
			  as_factor = scat_factors[an_1] * scat_factors[an_2] * dw_effect;
			  if (a != b)
				    {
			        xrd_data[c] += as_factor * sin(scat_factor)/scat_factor;
				    }
			  else
				  //Self-scattering, which is weighted half as much as a-b
				  //interactions.
				    {
				    xrd_data[c] += as_factor/2.0;
				    }
		     }
	     }
	   }

	 //Memory release---

	 delete[] scat_factors;
	 if (surf_effect)
		{
		delete[] surface_list;
		}
     }
    
//Spatial violations check---

bool atom_collection::No_Covalent_Overlap(int index)
     {
	 int temp;
     bool overlap = No_Covalent_Overlap(0, num_atoms, index, index + 1, temp);
	 return overlap;     
     }

bool atom_collection::No_Covalent_Overlap(int index_min, int index_max)
     {
	 int temp;
     bool overlap = No_Covalent_Overlap(0, num_atoms, index_min, index_max, temp);
	 return overlap;
     }

bool atom_collection::No_Covalent_Overlap(int comp_index_min, int comp_index_max,
	                                      int index_min, int index_max)
     {
	 int temp;
     bool overlap = No_Covalent_Overlap(comp_index_min, comp_index_max, index_min, 
		                                index_max, temp);
	 return overlap;
     }

bool atom_collection::No_Covalent_Overlap(int comp_index_min, int comp_index_max,
	                                      int index_min, int index_max, int& overlap_index)
     {
	 overlap_index = -1;
	  //Default result if no covalent overlap is found.
	 
	 bool no_overlap = TRUEV; 
	 for (int a = index_min; a < index_max; ++a)
	  //2nd passed range of indices (e.g. atoms being added to a composite).
	    {
		for (int b = comp_index_min; b < comp_index_max; ++b)
		 //1st passed range of indices.
		   {
           if ( (b >= index_min) && (b < index_max) )
		    //Prevent index range overlap in the analysis.
		      {
		      continue;
		      }

			 if (atom_array[a].Covalent_Overlap(atom_array[b], periodicity, box_dimensions))
			  {
			  no_overlap = FALSEV;
              a = index_max;
			  b = comp_index_max;
			  overlap_index = b;
              }
		   }
	    }

	 return no_overlap;
     }

bool atom_collection::No_Molecule_Overlap()
	 {
	 bool no_overlap = TRUEV;

	 for (int a = 0; a < num_atoms; ++a)
		 {
		 for (int b = a + 1; b < num_atoms; ++b)
			 {
			 if (!atom_array[a].Same_Molecule(atom_array[b]))
			  //Only compare atoms not in the same molecule.
				{
		        if (atom_array[a].Covalent_Overlap(atom_array[b], periodicity, box_dimensions))
					{
					no_overlap = FALSEV;
					a = num_atoms;
					b = num_atoms;
					}
				}
			 }	 
		 }
	 
	 return no_overlap;
	 }

bool atom_collection::No_Molecule_Overlap(int ind1, int ind2)
	 {
	 bool no_overlap = TRUEV;

	 //Get list of atoms belonging to those two molecules---

	 int num_atoms_in_mol1, num_atoms_in_mol2;
	 int mol_size1 = Molecule_Size(ind1);
	 int mol_size2 = Molecule_Size(ind2);
	 int* molecule1 = new int[mol_size1];
	 int* molecule2 = new int[mol_size2];
	 Get_Molecule_List(molecule1, mol_size1, ind1, num_atoms_in_mol1, TRUEV);
	 Get_Molecule_List(molecule2, mol_size2, ind2, num_atoms_in_mol2, TRUEV);

	 //Search for covalent overlap---
	 
	 int index_temp1, index_temp2;
	 for (int a = 0; a < num_atoms_in_mol1; ++a)
	    {
		index_temp1 = molecule1[a];
		for (int b = 0; b < num_atoms_in_mol2; ++b)
		   {
		   index_temp2 = molecule2[b];
		   if (atom_array[index_temp1].Covalent_Overlap
			   (atom_array[index_temp2], periodicity, box_dimensions))
			  {
			  no_overlap = FALSEV;
              a = num_atoms_in_mol1;
			  b = num_atoms_in_mol2;
              }
		   }
	    }

	 delete[] molecule2;
	 delete[] molecule1;

	 return no_overlap;
	 }

void atom_collection::Monolayer_Overlap_Prevention(int comp_index_min, int comp_index_max, int index_min, 
	                                               int index_max, double angle, int num_tries)
     {
     int overlap_index;
	 No_Covalent_Overlap(comp_index_min, comp_index_max, index_min, index_max, overlap_index);
	 if (overlap_index != -1)
	  //Overlap problem exists...
	     {
	     bool covalent_overlap = TRUEV;
		 double group_axis[3];
         Get_Bond_Vec(index_min, index_max - 1, group_axis);
		  //Rotation is done along the connecting vector over the range of atoms
		  //being considered.
		 double origin[3];
		 atom_array[index_min].Get_Atom_Location(origin);
		  //Origin in rotation is the first atom in the group being considered.
		 int molecule_size = index_max - index_min;
	     for (int a = 0; a < num_tries; ++a)
		  //Try reorientation along the molecular axis by (angle), (angle*2), etc. to prevent
		  //a spatial overlap problem.
			     {
				 Rotate_Atoms(molecule_size, group_axis, angle, origin);
			     if (No_Covalent_Overlap(comp_index_min, comp_index_max, index_min, index_max))
				    {
                    covalent_overlap = FALSEV;
					break;
				    }
			     }

		 if (covalent_overlap)
		  //Covalent overlap was not able to be prevented, so the group must be removed.
		    {
            Monolayer_Cutter(overlap_index, index_min, index_max); 
		    }
	     }
	 }

void atom_collection::Remove_Atoms_In_The_Way(int remove_range_min, int remove_range_max,
											  int range_min, int range_max)
	 {
	 bool remove_range_caution = FALSEV;
	 if (remove_range_min < range_min)
	  //Determine if the atoms to be removed are behind (in memory) the
	  //other atoms being considered. Important to how vector deletion is done.
		{
		remove_range_caution = TRUEV;
		}
	 for (int a = remove_range_min; a < remove_range_max; ++a)
	  //For each atom that might be removed...
		 {
		 if (!No_Covalent_Overlap(a, a + 1, range_min, range_max))
		  //Check the atom's covalent overlap to see if it is
		  //"in the way" of the considered range of atoms.
			{
		    Delete_Atom(a);
			--remove_range_max;
		     //Due to vector deletion, it is necessary to decrement
			 //the index range maximum.
			--a;
			 //Accordingly, the search index is also decreased.
			if (remove_range_caution)
			 //Check if that deletion requires the range of comparisons to
			 //also be updated.
				{
				--range_min;
				--range_max;
				}
			}
		 }
	 }

void atom_collection::Monolayer_Cutter(int group_index, int index_min, int index_max)
     {

	 if (group_index != -1)
	    {
		int molecule_size = index_max - index_min;
		int addition_parameter;
		if (group_index > index_min)
		   {
           addition_parameter = ( ( (group_index - index_max)/molecule_size) + 1) * molecule_size;
		   }
		else
		   {
           addition_parameter = ( ( (group_index - index_min + 1)/molecule_size) - 1) * molecule_size;
		   }
	    int group_start_index = index_min + addition_parameter;
		atom_array.erase(atom_array.begin() + group_start_index, 
			             atom_array.begin() + group_start_index + molecule_size);
		atom_array.erase(atom_array.begin() + index_min, 
			             atom_array.begin() + index_max);
		num_atoms -= 2*molecule_size;
	    }
     }

//Neighbor list generation---

//Bond network functions---

void atom_collection::Generate_Bond_Network(int* network, int start_index, 
	                  int max_index, int& molecule_count) const
     {
	 for (int a = start_index; a < max_index; ++a)
      /*Initialize molecule matrix (no molecule = -1).
	    This matrix is a list of indices in which all atoms with the same index 
		are a part of the same molecule.
	    E.g. Array = [ 2 2 1 1 5 1 1 2 ] --> 1st/2nd/last atom are a part of the 
		same molecule and the fifth atom is not a part of any molecule.*/
	     {
		 network[a] = -1;
	     }
	 
	 int bonding_index = -1;
	  //Count of the number of molecules found minus one.
	 int index_temp;
	 for (int a = start_index; a < max_index; ++a)
	     {
		 for (int b = a + 1; b < max_index; ++b)
		  //Search through all atom interactions, creating
		  //the overall bonding matrix.
		     {
			 if (atom_array[a].Covalent_Overlap(atom_array[b]))
			  //Found a bond from A to B.
			    {
				if ( (network[a] == -1) && (network[b] == -1) )
				 //Both atoms are not yet part of any bonding matrix.
				 //Hence, a new molecule has been found.
				      {
					  ++bonding_index;
				      network[a] = bonding_index;
					  network[b] = bonding_index;
				      }
				else if ( (network[a] == -1) && (network[b] != -1) )
				 //Atom A is not a part of a molecule matrix, but B is.
				 //Atom A is thus to be absorbed into atom B's bonding matrix.
				      {
					  network[a] = network[b];
				      }
				else if ( (network[a] != -1) && (network[b] == -1) )
				 //Atom B is to be absorbed into atom A's bonding matrix.
				      {
					  network[b] = network[a];
				      }
				else if ( (network[a] != network[b]) )
				 /*Special case: Both atoms are part of a bonding matrix, but
				   not the same matrix. This means that the program has found
				   two molecules that are actually only one molecule, and
				   a merger of the two molecule matrices must be performend.*/
				      {
				      index_temp = network[b];
					  for (int c = start_index; c < max_index; ++c)
					      {
						  if (network[c] == index_temp)
						   //Molecule B becomes a part of molecule A.
						     {
							 network[c] = network[a];
						     }
						  if (network[c] > index_temp)
						   //Due to the loss of a molecule in indexing, all
						   //molecule indices greater than that of the
						   //molecule index being removed must be decremented.
						     {
                             --network[c];
						     }
					      } 
					  --bonding_index;
				      }
			    }
		     }
	     }

	 for (int a = start_index; a < max_index; ++a)
	   //Take care of lone atoms. Assign them
	   //a dedicated molecule index.
	   {
	   if (network[a] == -1)
	      {
		  ++bonding_index;
          network[a] = bonding_index;
	      }
	   }

	 molecule_count = bonding_index + 1;
	  //First molecule has index = 0, so the molecule count is one greater.
     }

void atom_collection::Get_Molecule_Composition_String(char* molecule_comp, int* b_net,
	                                                  int index, const int MAX_STRING_LENGTH
													  ) const
     {
	 char** names;
	 int name_count;
	 Get_Memory(names, MAX_ATOM_TYPES, MAX_ANAME_SIZE);
	 Get_Name_List(names, name_count, MAX_ATOM_TYPES);
	  //Obtain a list of atom types in the atomic collection.

	 int char_index = 0;
	  //Position in molecule composition string.
	 int atom_count;
	 for (int a = 0; a < name_count; ++a)
	    {
		atom_count = Count_Atom_Type(names[a], b_net, index);
	     //Detemine the number of atoms in the molecule.
		if (atom_count != 0)
		   {
		   
		   //Addition of atom name---

		   strcpy(&(molecule_comp[char_index]), names[a]);
		   char_index += strlen(names[a]);

		   //Addition of name-number spacer---

		   molecule_comp[char_index] = '_';
		   ++char_index;

		   //Addition of atom count---

		   molecule_comp[char_index] = '\0';
		   Add_Number(molecule_comp, atom_count);
		   char_index += String_Size(atom_count);

		   //Addition of number-name spacer---

		   molecule_comp[char_index] = ',';
		   ++char_index;

		   if ( (char_index > (MAX_STRING_LENGTH - MAX_ANAME_SIZE - 10) ) )
		    //String memory space violation check.
			 {
			 Show_Warning("NOT ENOUGH STRING SPACE FOR MOLECULE COMPOSITIONS!");
			 break;
			 }
		   }
	    }

	 molecule_comp[char_index - 1] = '\0';
	  //Final string-terminating character.
	 Free_Memory(names, MAX_ATOM_TYPES);
	  //Memory release.
     }

void atom_collection::Clean_Up_Them_Molecules(int dif_pam)
     {

	 //Determine the bonding network of the atomic collection---

	 int* b_net = new int[num_atoms];
	 int original_molecule_count;
	 Generate_Bond_Network(b_net, 0, num_atoms, original_molecule_count);

	 //Get a list of the unique atom types---

	 char** names;
	 int name_count;
	 Get_Memory(names, MAX_ATOM_TYPES, MAX_ANAME_SIZE);
	 Get_Name_List(names, name_count, MAX_ATOM_TYPES);

	 //Obtain a molecule atom count array, which contains
	 //the number of each atom type in each molecule----

	 int** molecule_atom_count;
	 Get_Memory(molecule_atom_count, original_molecule_count, name_count + 2);
	  //Extra space (+ 2) is used for another purpose: comparisons between molecules.

	 bool equal_molecule_found, bigger_molecule_found, smaller_molecule_found;
	  //Molecule fragment variables, e.g. OH is smaller molecule than H2O
	  //and H3O is larger molecule than H2O.
	 for (int a = 0; a < original_molecule_count; ++a)
	  //Fills a list of the number of atoms of each type in each molecule.
	     {
		 for (int b = 0; b < name_count; ++b)
		     {
             molecule_atom_count[a][b] = Count_Atom_Type(names[b], b_net, a);
		     }
		 molecule_atom_count[a][name_count] = 0;
		 molecule_atom_count[a][name_count + 1] = 0;
		  //Initialize last two variables in the array, which are used
		  //for a different purposes: molecular comparisons.
	     }

	 //Determine, for each molecule, how many other molecules
	 //that exist that are larger or smaller versions
	 //of that molecule---

	 //To clarify on the molecule counters with examples:
	 //Fragment counter: OH molecule has 5 molecules larger than it (e.g. H2O).
	 //Like molecule counter: OH molecule has a population of 3.

	 for (int a = 0; a < original_molecule_count; ++a)
	    {
	 	for (int b = a + 1; b < original_molecule_count; ++b)
		     {
			 bigger_molecule_found = TRUEV;
			 equal_molecule_found = TRUEV;
			 smaller_molecule_found = TRUEV;
			 for (int c = 0; c < name_count; ++c)
			  /*In order for a smaller or larger molecule to be truly found,
			    all atom counts must compare accordingly.
			    E.g. CH4 is nothing compared to CH2O as one atom count
			    is smaller and one is larger.*/
			     {
				 if (molecule_atom_count[b][c] != molecule_atom_count[a][c])
				    {
					equal_molecule_found = FALSEV;
				    }
				 if (molecule_atom_count[b][c] < molecule_atom_count[a][c])
				    {
				    bigger_molecule_found = FALSEV;
				    }
				 if (molecule_atom_count[b][c] > molecule_atom_count[a][c])
				    {
				    smaller_molecule_found = FALSEV;
				    }
			     }
			 if (equal_molecule_found)
			  //Add up like molecule count (e.g. 100 H2O molecules) 
			  //for both molecules.
				 {
				 ++molecule_atom_count[b][name_count];
				 ++molecule_atom_count[a][name_count];
			     }
			 else if (bigger_molecule_found)
			  //Add up fragment-counter (e.g. 5 molecules smaller than H2O)
			  //for smaller molecule.
			     {
				 ++molecule_atom_count[a][name_count + 1];
			     }
			 else if (smaller_molecule_found)
			  //Add up fragment-counter for smaller molecule.
			     {
				 ++molecule_atom_count[b][name_count + 1];
			     }
		     }
	    }

	 //Now that like molecule counts and fragment counts have been totaled,
	 //the next step is to look for molecules with a high fragment
	 //count and lower like molecule count and eliminate them---

	 int molecule_count = original_molecule_count;
	 int orig_size = num_atoms;
	 int num_atoms_removed = 0;
	 for (int a = 0; a < original_molecule_count; ++a)
	     {
	     if (molecule_atom_count[a][name_count + 1] > 
			 ( molecule_atom_count[a][name_count] + dif_pam ) )
			//Remove criterion: Fragment count > Like molecule count + dif_pam
		    {
			--molecule_count;
	        for (int b = 0; b < num_atoms; ++b)
	           {
		       if (b_net[b] == a)
				//Found an atom that is a part of the fragment to be
				//removed, so delete it.
		          {
				  for (int c = b; c < (num_atoms - 1); ++c)
				   //Update change to bonding network by the deletion
				   //of the atom.
				      {
					  b_net[c] = b_net[c + 1];
				      }
			      Delete_Atom(b);
				   //Remove an atomic piece of the fragment.
			      ++num_atoms_removed;
				  --b;
				   //Vector deletion requires the b index to remain the same
				   //for the next run through of the loop.
		          }
	           }
		    }
	     }

	 //Memory release---

	 Free_Memory(molecule_atom_count, original_molecule_count);
	 Free_Memory(names, MAX_ATOM_TYPES);
	 delete[] b_net;
     }

//Special functions---

bool atom_collection::Check_Strange_Location(double dc_min, double dc_required, bool include_hydro) const
 //Note: Does not change atom collection cast's data members.
	 {
	 bool is_strange = FALSEV;	
	 double distance_temp;
	 bool atom_is_weird_by_isolation, atom_is_weird_by_closeness;
	 for (int a = 0; a < num_atoms; ++a)
	    {
		atom_is_weird_by_isolation = TRUEV;
		atom_is_weird_by_closeness = FALSEV;
		 //Reset weirdness booleans. The following
		 //atom distance search will see if either
		 //boolean need be switched.
		if (atom_array[a].Is_Invisible(include_hydro))
		  //Check for invisible atoms.
		     {
			 continue;
		     }

		for (int b = 0; b < num_atoms; ++b)
		   {
		   if (a == b)
		     {
			 continue;
		     }

		   if (atom_array[b].Is_Invisible(include_hydro))
		  //Check for invisible atoms.
		     {
			 continue;
		     }

		   distance_temp = atom_array[a].Get_Distance(atom_array[b], 
			               periodicity, box_dimensions);
           if (distance_temp < dc_required)
		    //Isolation strangeness not a factor for this atom.
              {
			  atom_is_weird_by_isolation = FALSEV;
              }
		   if (distance_temp < dc_min)
			 //Atoms are too close!
		      {
			  atom_is_weird_by_closeness = TRUEV;
			  b = num_atoms;
		      }
		   }
		if (atom_is_weird_by_isolation || atom_is_weird_by_closeness)
		   {
		   is_strange = TRUEV;
		   a = num_atoms;
		   }
	    }

	 return is_strange;
	 }

//I/O----

void atom_collection::Prepare_Collection_For_Output()
    {
    for (int a = 0; a < num_atoms; ++a)
	    {
	    Check_FP_Zeros(atom_array[a].atom_location);
		 //Prevent atoms located at the zero coordinate along
		 //a Cartesian axis from creating ugly output
		 //of FP values in the data files.
		Check_FP_Zero(atom_array[a].atomic_charge);
		Check_FP_Zero(atom_array[a].atomic_radius);
	    }
	Defrag_Group_Tags();
	 //Organize group tags so that no numbers in the group
	 //tag range are left unused (e.g. 1 2 5 --> 1 2 3).
    }

void atom_collection::Print_Atoms_Location() const
       {
       for (int a = 0; a < num_atoms; ++a)
           {
           atom_array[a].Print_Atom_Location();      
           }                    
       }

void atom_collection::Store_Atoms_Location(ofstream& out_file, bool include_hydro) const
       {
       for (int a = 0; a < num_atoms; ++a)
           {
           if (!atom_array[a].Is_Invisible(include_hydro))
		    //Don't include invisible atoms.
              {
              atom_array[a].Store_Atom_Location(out_file);      
              }
           }                                          
       }
       
void atom_collection::Store_Atoms_Location_With_Charge(ofstream& out_file, 
                                                       bool include_hydro) const
       {
       for (int a = 0; a < num_atoms; ++a)
           {
           if (!atom_array[a].Is_Invisible(include_hydro))
		    //Don't include invisible atoms.
              {
              atom_array[a].Store_Atom_Location_With_Charge(out_file);      
              }
           }                                       
       }

void atom_collection::Store_Atoms_Location_With_CoreShell(ofstream& out_file, 
                                                          bool include_hydro,
														  bool cat_coreshell,
														  bool ani_coreshell) const
       {
	   double atom_charge;
	   bool needs_shell;
       for (int a = 0; a < num_atoms; ++a)
           {
           if (!atom_array[a].Is_Invisible(include_hydro))
			//Don't include charges.
              {
			  needs_shell = FALSEV;
			  atom_charge = atom_array[a].Get_Atom_Charge();
			  if (cat_coreshell && Check_FP_GreaterOrEqual(atom_charge, 0.00))
			   //If charge is greater than 0 and cation shells are on... 
				 {
				 needs_shell = TRUEV;
				 }
			  else if (ani_coreshell && Check_FP_LessThan(atom_charge, 0.00))
			   //Anion case.
				 {
				 needs_shell = TRUEV;
				 }
			  if (needs_shell)
			   //Output a line for both the core and shell of the atom.
				 {
				 atom_array[a].Store_Atom_Location_With_CoreShell(out_file);      
				 }
			  else
			   //Output a line for the atom only (one line).
				 {
			     atom_array[a].Store_Atom_Location(out_file);
				 }
			  }
           }                                       
       }

void atom_collection::Store_Atom_Types_With_Masses(ofstream& out_file, 
	                                               bool include_hydro) const
     {
	 int indices[MAX_ATOM_TYPES];
	 int num_types;
     Get_Index_List(indices, num_types);
		//Get a list of all unique atom types in terms of representative
	    //atoms.
	 for (int a = 0; a < num_types; ++a)
	    {
		if (!atom_array[indices[a]].Is_Invisible(include_hydro))
		 //Ignore invisible atoms.
		  {
	      out_file << "mass " << atom_array[indices[a]].atom_name << " "
	               << atom_array[indices[a]].rel_atomic_mass << endl;
		  }
		}
     }

void atom_collection::Atoms_Storage(ofstream& out_file) const
     {
     out_file << endl << num_atoms;
	  //Store number of atoms on a distinct line in the file.
     for (int a = 0; a < num_atoms; ++a)
	  //Store each atom, line by line.
         {
         atom_array[a].Store_Atom_Info(out_file);  
         }
     }
     
void atom_collection::Atoms_Retrieval(ifstream& in_file)
     {
	 int temp_atom_count;
     in_file >> temp_atom_count;
	 Memory_Check(temp_atom_count + num_atoms);
	  //Reserve memory space for the atom loading.
     for (int a = num_atoms; a < temp_atom_count + num_atoms; ++a)
         {
         atom_array[a].Load_Atom_Info(in_file);   
         }        
	 num_atoms += temp_atom_count;
     }

