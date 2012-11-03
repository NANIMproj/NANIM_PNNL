
#include "stdafx.h"

#include "DWN_molecule.h"

molecule::molecule()
   {
   Initialize();
   }

molecule::molecule(const molecule& molec)
   {
   Initialize();
   *this = molec;               
   }

molecule::~molecule()
   {
   }
   
void molecule::Initialize()
   {
   molecule_atoms.Set_Initial_Allocation(INIT_MOLECULE_SIZE);
   }

//Information retrieval---

int molecule::Size() const
   {
   return molecule_atoms.Size();
   }

double molecule::Get_Charge() const
   {
   return molecule_atoms.Return_Total_Charge();
   }

//Properties retrieval---

void molecule::operator = (const molecule& molec)
   {
   molecule_atoms = molec.molecule_atoms;               
   }

//Building functions---

void molecule::Add_Atom(const atom& ze_atom, const double* pos)
   {
   if (!ze_atom.Is_Dummy_Atom())
	//Dummy atoms likely used by user for small changes
	//to the molecule.
      {
	  molecule_atoms.Add_Atom(ze_atom, pos, PLACE_MODE);                         
      }
   }

void molecule::Add_Atom(const atom& ze_atom, double x, double y, double z)
   {
   double coors[3] = { x, y, z };
   Add_Atom(ze_atom, coors);
   }
   
void molecule::Add_Atom(const atom& ze_atom, const double* pos, double bond_length, 
                        int geometry_index, int index, bool invert)
   {
   double pos_coors[3];
   double change_coors[3];
   for (int a = 0; a < 3; ++a)
	//Get the coordinates for the requested geometry.
       {
       if (geometry_index == LINEAR)
          {
          change_coors[a] = (LINEAR_POS[index][a] * bond_length);               
          }
       else if (geometry_index == TRIG)
          {
          change_coors[a] = (TRIG_POS[index][a] * bond_length); 
          }
       else if (geometry_index == TETRA)
          {     
          change_coors[a] = (TETRA_POS[index][a] * bond_length);
          }
       else
          {
		  Show_Warning("NON-PROGRAMMED GEOMETRY REQUESTED IN MOLECULE CONSTRUCTION!");
          }
       }

   if (invert)
	//Change to the opposite direction with the geometrical logic,
    //e.g. going left vs. right along a carbon backbone.
      {
	  Neg_XYZ(change_coors);
      }

   SetAdd_XYZ(pos_coors, pos, change_coors);
   Add_Atom(ze_atom, pos_coors);   
   }  

void molecule::Place_Water(const double* pos, const atom& water_atomO, 
                           const atom& water_atomH)
   {      
   Add_Atom(water_atomO, pos);
   Add_Atom(water_atomH, pos, OH_BOND, TETRA, 0, FALSEV);
   Add_Atom(water_atomH, pos, OH_BOND, TETRA, 1, FALSEV);
    //H2O is the nectar of life, baby.              
   }

//AXYZ addition---
   
void molecule::Add_CH3(const double* prev, const double* pos, const atom& hc_atomC, 
                       const atom& hc_atomH)
   {      
   Add_AXYZ(prev, pos, hc_atomC, hc_atomH, hc_atomH, hc_atomH, 
            CH_BOND, CH_BOND, CH_BOND, TETRA);                    
   }  

void molecule::Add_AXYZ(const double* prev, const double* pos, const atom& atomA, 
                        const atom& atomX, const atom& atomY, const atom& atomZ, 
                        double AX_bond_length, double AY_bond_length, 
                        double AZ_bond_length, int shape_factor)
   {
   //Get R-A bond description (pos vector - prev vector)---

   double prev_bond[3];
   double first_bond_length;
   Get_Vec_And_Mag(prev_bond, first_bond_length, prev, pos);
    
   //Set down an initial AXYZ group; a final rotation 
   //will get things into place---

   double atomA_pos[3];
   double temp_RA_bond[3];
   Add_Atom(atomA, prev, first_bond_length, shape_factor, VERTEX_A, FALSEV);
    //R-A atom placement goes along the opposite direction as does the
    //A-X/Y/Z atom placements within the specified geometric shape.
   molecule_atoms[molecule_atoms.Size() - 1].Get_Atom_Location(atomA_pos);
   SetSub_XYZ(temp_RA_bond, prev, atomA_pos);
   Add_Atom(atomX, atomA_pos, AX_bond_length, shape_factor, VERTEX_B, TRUEV); 
   Add_Atom(atomY, atomA_pos, AY_bond_length, shape_factor, VERTEX_C, TRUEV);
   Add_Atom(atomZ, atomA_pos, AZ_bond_length, shape_factor, VERTEX_D, TRUEV);   
   if (Check_FP_GreaterThan(first_bond_length, 0.0))
    //Passing the same values for "prev" and "pos" turns addition to a molecule
    //into simple atom placement. If that is the case, no rotation is needed.
      { 
      molecule_atoms.Rotate_Atoms(4, temp_RA_bond, prev_bond, prev);                                 
      }                

   }

//A-X-Y addition (A is the central atom for molecule construction)---

void molecule::Add_CH2(const double* prev, const double* pos, 
                       const atom& me_atomC, const atom& me_atomH)
   {
   Add_AXY(prev, pos, me_atomC, me_atomH, me_atomH, CH_BOND, CH_BOND, TETRA);                    
   }

void molecule::Add_AXY(const double* prev, const double* pos, const atom& atomA, 
                       const atom& atomX, const atom& atomY, double AX_bond_length, 
                       double AY_bond_length, int shape_factor)
   {
   //Get R-A bond description (pos vector - prev vector)---

   double prev_bond[3];
   double first_bond_length;
   Get_Vec_And_Mag(prev_bond, first_bond_length, prev, pos);
   
   //Set down an initial AXY group; a final rotation 
   //will get things into place---

   double atomA_pos[3];
   double temp_RA_bond[3];
   Add_Atom(atomA, prev, first_bond_length, shape_factor, VERTEX_A, FALSEV);
   molecule_atoms[molecule_atoms.Size() - 1].Get_Atom_Location(atomA_pos);
   SetSub_XYZ(temp_RA_bond, prev, atomA_pos);
   Add_Atom(atomX, atomA_pos, AX_bond_length, shape_factor, VERTEX_B, TRUEV); 
   Add_Atom(atomY, atomA_pos, AY_bond_length, shape_factor, VERTEX_C, TRUEV);  
   if (Check_FP_GreaterThan(first_bond_length, 0.0))
    //Passing the same values for "prev" and "pos" turns addition to a molecule
    //into simple atom placement. If that is the case, no rotation is needed.
      { 
      molecule_atoms.Rotate_Atoms(3, temp_RA_bond, prev_bond, prev);                            
      }                   
   }

//A-X addition (A is the central atom for molecule construction)---

void molecule::Add_SH(const double* prev, const double* pos, const atom& sh_atomS, 
                      const atom& sh_atomH)
   {
   Add_AX(prev, pos, sh_atomS, sh_atomH, SH_BOND, TETRA);                           
   }

void molecule::Add_OH(const double* prev, const double* pos, const atom& oh_atomO, const atom& oh_atomH)
   {   
   Add_AX(prev, pos, oh_atomO, oh_atomH, OH_BOND, TETRA);                      
   }

void molecule::Add_AX(const double* prev, const double* pos, const atom& atomA, const atom& atomX, 
                      double AX_bond_length, int shape_factor)
   {
   //Get R-A bond description (pos vector - prev vector)---

   double prev_bond[3];
   double first_bond_length;
   Get_Vec_And_Mag(prev_bond, first_bond_length, prev, pos);

   //Set down an initial AX group; a final rotation 
   //will get things into place---
   
   double atomA_pos[3];
   double temp_RA_bond[3];
   Add_Atom(atomA, prev, first_bond_length, shape_factor, VERTEX_A, FALSEV);
   molecule_atoms[molecule_atoms.Size() - 1].Get_Atom_Location(atomA_pos);
   SetSub_XYZ(temp_RA_bond, prev, atomA_pos);
   Add_Atom(atomX, atomA_pos, AX_bond_length, shape_factor, VERTEX_B, TRUEV);  
   if (Check_FP_GreaterThan(first_bond_length, 0.0))
    //Passing the same values for "prev" and "pos" turns addition to a molecule
    //into simple atom placement. If that is the case, no rotation is needed.
      { 
      molecule_atoms.Rotate_Atoms(2, temp_RA_bond, prev_bond, prev);                            
      }
   }

//Large portion additions---

void molecule::Add_Carbon_Chain(const double* prev, const double* pos, 
                                double* final_pos, double* next_pos, 
                                int numb_carbons, const atom& hc_atomC, 
                                const atom& hc_atomH)
   {

   //Characterize the first bond to be added in the carbon chain
   //(pos - prev) vector---

   double initial_pos_diff[3];
   double first_bond_length;
   Get_Vec_And_Mag(initial_pos_diff, first_bond_length, prev, pos);
   double first_bond_ratio = first_bond_length/CC_BOND;
    //X-C bond / C-C bond ratio, which will make later logic
    //a little easier.

   //Get the stepping parameters that describe how to step along the
   //carbon backbone (C-C-C-C-etc.) of the chain as in
   //tetrahedral-tetrahedral-tetrahedral-etc. connectivity---
   
   double steps[4][3];
    //Steps that take the carbon placement algorithm along a carbon backbone.
   Set_XYZ(steps[0], TETRA_POS[VERTEX_A]);
   Set_XYZ(steps[1], TETRA_POS[VERTEX_D]);
   Neg_XYZ(steps[1]);
   Set_XYZ(steps[2], TETRA_POS[VERTEX_B]);
   Set_XYZ(steps[3], TETRA_POS[VERTEX_C]);
   Neg_XYZ(steps[3]);
   for (int a = 0; a < 4; ++a)
       {
       for (int b = 0; b < 3; ++b)
           {
           steps[a][b] *= CC_BOND;  
            //Set the steps relative to the carbon-carbon bond length.   
           }     
       } 

   //Get the stepping parameters for the very first X-C bond---

   double true_first_steps[3];
   Set_XYZ(true_first_steps, steps[0]);
   Multiply_XYZ(true_first_steps, first_bond_ratio);

   //Get parameters for hydrogen placement (two hydrogens for each
   //methyl group) as we move along the carbon backbone. The following
   //parameters describe how to use a tetrahedral description to do so.---

   int hydrogen_pos[4][2];
   bool inversion[4][2];
   hydrogen_pos[0][0] = 1;
   hydrogen_pos[0][1] = 2;
   inversion[0][0] = inversion[0][1] = TRUEV;
   hydrogen_pos[1][0] = 0;
   hydrogen_pos[1][1] = 2;
   inversion[1][0] = inversion[1][1] = FALSEV;
   hydrogen_pos[2][0] = 3;
   hydrogen_pos[2][1] = 0;
   inversion[2][0] = inversion[2][1] = TRUEV;
   hydrogen_pos[3][0] = 1;
   hydrogen_pos[3][1] = 3;
   inversion[3][0] = inversion[3][1] = FALSEV;

   //Create the backbone---

   double backbone_pos[3];
   Set_XYZ(backbone_pos, prev);
    //Start the backbone at the X group.
   int step_index;
    //Index that indicates which stepping/placement parameters to use.
   int first_atom_index = molecule_atoms.Size();
    //Record the initial size of the molecule before 
    //carbon backbone addition.
   for (int a = 0; a < ( numb_carbons + 1); ++a)
    //Add one more carbon than that requested to make the assignment of the "next
    //atom position" easier.
       {
       step_index = a % 4;
       if (a != 0)
          {
          Add_XYZ(backbone_pos, steps[step_index]); 
          }
       else
		//First addition --> X--(CH2) group.
          {
          Add_XYZ(backbone_pos, true_first_steps);
          }
       Add_Atom(hc_atomC, backbone_pos);
       Add_Atom(hc_atomH, backbone_pos, CH_BOND, TETRA, 
                hydrogen_pos[step_index][0], inversion[step_index][0]);     
       Add_Atom(hc_atomH, backbone_pos, CH_BOND, TETRA, 
                hydrogen_pos[step_index][1], inversion[step_index][1]); 
        //Add the methyl group.
       }
       
   //Finally, rotate the initial X--C bond such that it is in the
   //requested location, (pos - prev) vector and get rid of that
   //extra methyl group added just for spatial analysis---
   
   int num_atoms_added = molecule_atoms.Size() - first_atom_index;
   molecule_atoms.Rotate_Atoms(num_atoms_added, true_first_steps, 
	                           initial_pos_diff, prev);
    //Rotate the carbon chain into place.
   molecule_atoms[molecule_atoms.Size() - 6].Get_Atom_Location(final_pos);
    //Set the final position to the last carbon atom added.     
   molecule_atoms[molecule_atoms.Size() - 3].Get_Atom_Location(next_pos);
    //Set the position where a next carbon atom would be added.
   molecule_atoms.Delete_Atoms(molecule_atoms.Size() - 3, molecule_atoms.Size());
    //Undo the addition of that extra, last methyl group. It was only added
    //in order to figure out the value of "next_pos."
   }

//Specific molecule building---

void molecule::Build_Alkanethiol(int num_methyl_groups, const atom& s_atomS, 
                                 const atom& CH2_atomC, const atom& CH3_atomC, 
								 const atom& hc_atomH)
   {
   double sulfur_pos[3] = { 0.0, 0.0, 0.0 };
   Add_Atom(s_atomS, sulfur_pos);
    //Add the initial sulfur head group.
   double first_carbon_pos[3]; 
   Set_XYZ(first_carbon_pos, TETRA_POS[VERTEX_C]);
   Multiply_XYZ(first_carbon_pos, SC_BOND);
   double last_carbon_pos[3];
   double next_atom_pos[3];
   Add_Carbon_Chain(sulfur_pos, first_carbon_pos, last_carbon_pos, 
                    next_atom_pos, num_methyl_groups, CH2_atomC, hc_atomH);
    //Add the carbon backbone.
   Add_CH3(last_carbon_pos, next_atom_pos, CH3_atomC, hc_atomH); 
    //Add the final methyl group.
   Standard_Orientation();        
   }

void molecule::Add_Amino_Acid(const double* prev, const double* pos, 
							  double* next_pos, const atom& co_atomC, 
							  const atom& co_atomO, const atom& cr_atomC, 
							  const atom& cr_atomH, const atom& n_atomN, 
							  const atom& n_atomH, bool is_terminal)
   {
   double new_pos[3];
   double old_pos[3];

   //Set up dummy atom R, which is to be replaced with the amino acid signature 
   //side-chain at a later point---

   atom R_atom;
   R_atom.Set_Atom_Name(RGROUP_NAME);

   //CO group----

   Add_AXY(prev, pos, co_atomC, co_atomO, cr_atomC, CO_DBOND, CC_BOND, TRIG);
   Evaluate_Addition(new_pos, 1, old_pos, 3);

   //CRHN group----

   Add_AXYZ(old_pos, new_pos, cr_atomC, R_atom, cr_atomH, n_atomN,
	        R_BOND, CH_BOND, CN_BOND, TETRA);
   Evaluate_Addition(new_pos, 1, old_pos, 4);

   //NH2 group----
   
   if (!is_terminal)
	  {
	  Add_AXY(old_pos, new_pos, n_atomN, n_atomH, n_atomH, NH_BOND, CN_BOND, TETRA);
	  molecule_atoms[molecule_atoms.Size() - 1].Get_Atom_Location(next_pos);
	  molecule_atoms.Delete_Atom(molecule_atoms.Size() - 1);
	  }
   else
	//Terminal group has one extra hydrogen atom.
      {
      Add_AXY(old_pos, new_pos, n_atomN, n_atomH, n_atomH, NH_BOND, NH_BOND, TETRA);
	  Zero_XYZ(next_pos);
	   //No next position if builder is at a terminal group.
      }
   }

void molecule::Evaluate_Addition(double* new_pos, int new_coor, 
	                             double* old_pos, int old_coor)
 //Gets spatial coordinates of desired atoms and deletes the last atom
 //placed in the molecule. Common need in my molecule construction algorithms.
   {
   molecule_atoms[molecule_atoms.Size() - new_coor].Get_Atom_Location(new_pos);
   molecule_atoms[molecule_atoms.Size() - old_coor].Get_Atom_Location(old_pos);
   molecule_atoms.Delete_Atom(molecule_atoms.Size() - 1);
   }

void molecule::Build_Peptide(int num_amino_acids, const atom& co_atomC, 
                             const atom& co_atomO, const atom& co_atomOT, 
							 const atom& co_atomHT, const atom& cr_atomC, 
							 const atom& cr_atomH, const atom& n_atomN, 
							 const atom& n_atomH)
   {
   double start[3] = { 0.0, 0.0, 0.0 }; 
   Add_AX(start, start, co_atomOT, co_atomHT, OH_BOND, TETRA);
    //Add terminal -OH acid group.
   double next_carbon_pos[3];
   Set_XYZ(next_carbon_pos, TETRA_POS[VERTEX_C]);
   Multiply_XYZ(next_carbon_pos, CO_BOND);
   double prev_atom_pos[3];
   Set_XYZ(prev_atom_pos, start);
   bool final_add = FALSEV;
   for (int a = 0; a < num_amino_acids; ++a)
       {
	   if (a == (num_amino_acids - 1))
	      {
          final_add = TRUEV;
	      }
       Add_Amino_Acid(prev_atom_pos, next_carbon_pos, next_carbon_pos, co_atomC,
		              co_atomO, cr_atomC, cr_atomH, n_atomN, n_atomH, final_add);
	   molecule_atoms[molecule_atoms.Size() - 2].Get_Atom_Location(prev_atom_pos);
	    //Nitrogen atom added is built off of in the next amino acid
	    //addition step.
	   }
   Standard_Orientation();
   }

void molecule::Peptide_Substituent_Addition(const char* substituent_list, 
	                                        int num_amino_acids, const atom& c_atomC,
											const atom& c_atomH, const atom& s_atomS,
											const atom& s_atomH, const atom& o_atomO,
											const atom& o_atomH)
   {
   double C_prev[3];
   double R_pos[3];
   char name_temp[MAX_ANAME_SIZE];
   int num_added = 0;
   for (int a = 0; a < molecule_atoms.Size(); ++a)
      {
	  molecule_atoms[a].Get_Atom_Name(name_temp);
	  if (Name_Check(name_temp, const_cast<char*>(RGROUP_NAME)))
	     {
         molecule_atoms[a].Get_Atom_Location(R_pos);
		 molecule_atoms[a - 1].Get_Atom_Location(C_prev);
		   //Get spatial details of C-R bond, which indicates how to 
		   //place the substituent.
		 molecule_atoms.Delete_Atom(a);
		   //Erase dummy R atom.
		 if (substituent_list[num_added] == 'P')
		   //Special case: Proline. Hydrogen atom is removed from nitrogen.   
		    {
            molecule_atoms.Delete_Atom(a + 2);   
		    }

		 Place_R_Group(C_prev, R_pos, substituent_list[num_added], c_atomC, c_atomH, s_atomS, 
			           s_atomH, o_atomO, o_atomH);
         ++num_added;
	     if (num_added == num_amino_acids)
	        {
		    break;
	        }
	     }
      }
   }

void molecule::Place_R_Group(const double* prev, const double* pos, 
	                         char substituent_label, const atom& c_atomC,
							 const atom& c_atomH, const atom& s_atomS,
							 const atom& s_atomH, const atom& o_atomO,
							 const atom& o_atomH)
   {
   double new_pos[3];
   double old_pos[3];
   double new2_pos[3];
   double old2_pos[3];
   double temp_pos[3];

   switch (substituent_label)
      {
   
	  //Amino acids with electrically charged side chains---	   

      case 'R':
	    //Arginine-(CH2)3-NH-C-(NH2)2, positively charged end

	 	

	  case 'H': 
		//Histidine

	  case 'K':
       //Lysine

	  case 'D':
	   //Aspartic acid

	  case 'E': 
	   //Glutamic acid

	  case 'S':
	   //Serine

	  Add_AXYZ(prev, pos, c_atomC, c_atomH, c_atomH, o_atomO, 
		       CH_BOND, CH_BOND, CO_BOND, TETRA);
	  Evaluate_Addition(new_pos, 1, old_pos, 4);
	  Add_AX(new_pos, old_pos, o_atomO, o_atomH, OH_BOND, TETRA);

	  case 'T':
	   //Threonine

	  case 'N': 
	   //Asparagine

	  case 'Q':
	   //Glutamine

      case 'C':
	   //Cysteine

	  Add_AXYZ(prev, pos, c_atomC, c_atomH, c_atomH, s_atomS, 
		       CH_BOND, CH_BOND, SC_BOND, TETRA);
	  Evaluate_Addition(new_pos, 1, old_pos, 4);
	  Add_AX(new_pos, old_pos, s_atomS, s_atomH, SH_BOND, TETRA);

	  case 'U': 
       //Selenocysteine

	  case 'G':
       //Glycine

	  SetSub_XYZ(temp_pos, pos, prev);
	  Multiply_XYZ(temp_pos, CH_BOND/CC_BOND);
	  Add_XYZ(temp_pos, prev);
	   //Change bond length to be that of a CH bond.
	  Add_Atom(c_atomH, temp_pos);

	  case 'P':
	   //Proline

	  case 'A': 
	   //Alanine

	  Add_AXYZ(prev, pos, c_atomC, c_atomH, c_atomH, c_atomH, 
		       CH_BOND, CH_BOND, CH_BOND, TETRA);

	  case 'V':
	   //Valine

      Add_AXYZ(prev, pos, c_atomC, c_atomH, c_atomC, c_atomC, 
		       CH_BOND, CC_BOND, CC_BOND, TETRA);
	  Evaluate_Addition(new_pos, 1, old_pos, 4);
	  Evaluate_Addition(new2_pos, 1, old2_pos, 3);

	  Add_AXYZ(old_pos, new_pos, c_atomC, c_atomH, c_atomH, c_atomH, 
		       CH_BOND, CH_BOND, CH_BOND, TETRA);
	  Add_AXYZ(old2_pos, new2_pos, c_atomC, c_atomH, c_atomH, c_atomH, 
		       CH_BOND, CH_BOND, CH_BOND, TETRA);

	  case 'I':
	   //Isoleucine	    
 
	  case 'L': 
	   //Leucine

	  case 'M':
	   //Methionine

	  case 'F':
       //Phenylalanine

	  case 'Y': 
       //Tyrosine

	  case 'W':
       //Tryptophan

	  default: Show_Statement("PLEASE DO NOT MAKE UP AMINO ACIDS!");

      }
   }

//Utility functions---

void molecule::Return_Molecule_Center(double* cen_coors) const
   {
   molecule_atoms.Coordinate_Average(cen_coors);                            
   } 

void molecule::Translate_Molecule(const double* shift_coors)
   {
   molecule_atoms.Shift_Atomic_Coordinates(shift_coors);                   
   }

void molecule::ReAssign_Center(const double* new_cen_coors)
   {
   molecule_atoms.Translate_To_New_Coordinate_Average(new_cen_coors);                                                                       
   }     

void molecule::Rotate_Molecule(const double* start_vec, const double* final_vec)
   {  
   molecule_atoms.Rotate_Atoms(start_vec, final_vec);            
   }

void molecule::Rotate_Molecule(const double* start_vec, const double* final_vec, const double* origin)
   {
   molecule_atoms.Rotate_Atoms(start_vec, final_vec, origin);                                   
   }

void molecule::Rotate_Molecule(const double* rot_vec, double angle)
   {  
   molecule_atoms.Rotate_Atoms(rot_vec, angle);            
   }

void molecule::Rotate_Molecule(const double* rot_vec, double angle, const double* origin)
   {  
   molecule_atoms.Rotate_Atoms(rot_vec, angle, origin);            
   }

void molecule::Rotate_MoleculeXYZ(const double* rot_angles, bool clockwise)
   {  
   molecule_atoms.Rotate_AtomsXYZ(rot_angles, clockwise);            
   }

void molecule::Rotate_MoleculeXYZ(const double* rot_angles, bool clockwise, const double* origin)
   {
   molecule_atoms.Rotate_AtomsXYZ(rot_angles, clockwise, origin);
   }
      

void molecule::Orient_Molecule(double rot_angle)
   {
   double first_atom_pos[3];
   molecule_atoms[0].Get_Atom_Location(first_atom_pos);
   double mol_vec[3];
   molecule_atoms.Get_Bond_Vec(0, molecule_atoms.Size() - 1, mol_vec);
   Rotate_Molecule(mol_vec, rot_angle, first_atom_pos);
   }


void molecule::Orient_Molecule(const double* orientation)
   {
   double mol_vec[3];
   molecule_atoms.Get_Bond_Vec(0, molecule_atoms.Size() - 1, mol_vec);
   double first_atom_pos[3];
   molecule_atoms[0].Get_Atom_Location(first_atom_pos);
   Rotate_Molecule(mol_vec, orientation, first_atom_pos);
   }

void molecule::Orient_Molecule(const double* mol_orientation, const double* firstbond_orientation)
   {
   Orient_Molecule(mol_orientation);

   double first_atom_pos[3];
   molecule_atoms[0].Get_Atom_Location(first_atom_pos);
   double bond_vec[3];
   molecule_atoms.Get_Bond_Vec(0, 1, bond_vec);
   double rot_angle;
   Get_Angle(bond_vec, firstbond_orientation, rot_angle);
   Rotate_Molecule(mol_orientation, rot_angle, first_atom_pos);
   }

void molecule::Standard_Orientation()
   {
   double NEG_Z_AXIS[3] = { 0.0, 0.0, -1.0 };
    //Negative z-axis, e.g. as in coating of a nanosurface.
   Orient_Molecule(NEG_Z_AXIS);
   double bond_vec[3];
   molecule_atoms.Get_Bond_Vec(0, 1, bond_vec);
    //Get bond vector connecting the first two atoms in the molecule.
   bond_vec[0] = sqrt(bond_vec[0]*bond_vec[0] + bond_vec[1]*bond_vec[1]);
   bond_vec[1] = 0.0; 
    //Transform x/y components of bond vector into pure x component.
   Orient_Molecule(NEG_Z_AXIS, bond_vec);
   }
 
//Molecule copying---

int molecule::Copy_Molecule(atom_collection& atomic_array) const
   {
   int ind = atomic_array.Size();
   atomic_array.Copy_Atomic_Group(molecule_atoms);
   return molecule_atoms.Size();                            
   }
   
int molecule::Copy_Molecule(atom_collection& atomic_array, const double* pos, 
                            bool use_center) const
   {
   int ind = atomic_array.Size();
   molecule copy(*this);
   if (use_center)
      {
      copy.ReAssign_Center(pos);
      }
   if (!use_center)
      {
      double temp_coors[3];
      double start_coors[3];
      copy.molecule_atoms[0].Get_Atom_Location(start_coors);
      SetSub_XYZ(temp_coors, pos, start_coors);
      copy.molecule_atoms.Shift_Atomic_Coordinates(temp_coors);              
      }
   atomic_array.Copy_Atomic_Group(copy.molecule_atoms);
   return molecule_atoms.Size();                            
   }
   
int molecule::Copy_Molecule(atom_collection& atomic_array, const double* pos, 
                            const double* rot_angles, bool use_center) const
   {
   int ind = atomic_array.Size();
   molecule copy(*this);
   copy.Rotate_MoleculeXYZ(rot_angles, CLOCKWISE);
   if (use_center)
      {
      copy.ReAssign_Center(pos);
      }
   if (!use_center)
      {
      double temp_coors[3];
      double start_coors[3];
      copy.molecule_atoms[0].Get_Atom_Location(start_coors);
      SetSub_XYZ(temp_coors, pos, start_coors);
      copy.molecule_atoms.Shift_Atomic_Coordinates(temp_coors);              
      } 
   atomic_array.Copy_Atomic_Group(copy.molecule_atoms);
   return molecule_atoms.Size();                            
   }

int molecule::Copy_Molecule(atom_collection& atomic_array, const double* pos, 
                            const double* start_vec, const double* end_vec, 
							bool use_center) const
   {
   int ind = atomic_array.Size();
   molecule copy(*this);
   copy.Rotate_Molecule(start_vec, end_vec);
   if (use_center)
      {
      copy.ReAssign_Center(pos);
      }
   if (!use_center)
      {
      double temp_coors[3];
      double start_coors[3];
      copy.molecule_atoms[0].Get_Atom_Location(start_coors);
      SetSub_XYZ(temp_coors, pos, start_coors);
      copy.molecule_atoms.Shift_Atomic_Coordinates(temp_coors);              
      } 
   atomic_array.Copy_Atomic_Group(copy.molecule_atoms);
   return molecule_atoms.Size();                            
   }

//I/O---

void molecule::Store_Atom_Locations(ofstream& out_file) const
   {
   molecule_atoms.Store_Atoms_Location(out_file, TRUEV);                        
   }
   
void molecule::Save_Molecule(const char* save_file) const
   {
   ofstream ze_file(save_file, ios::out);
   if (ze_file.is_open())
      {
      Write_Intro(ze_file, "File Type = Molecule Description");
      molecule_atoms.Atoms_Storage(ze_file); 
	  ze_file.close();    
      }              
   }
   
void molecule::Load_Molecule(const char* load_file)
   {
   ifstream ze_file(load_file, ios::in);
   if (ze_file.is_open())
      {
      Skip_Phrases(ze_file, 5);
      molecule_atoms.Atoms_Retrieval(ze_file);  
	  ze_file.close();  
      }
   }
