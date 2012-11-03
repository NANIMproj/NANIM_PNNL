
#include "stdafx.h"

#include "DWN_nanos.h"

//Constructors, destructors, and initialization functions---

nanosurface::nanosurface()
      {
      Initialize();      
      }

nanosurface::nanosurface(const crystal_system& ze_cryst, double x_dim, double y_dim, 
                         double z_dim, int x_sur, int y_sur, int z_sur, double angle,
						 bool stuffed, bool vary_size, bool neut_surf)
      {
      Initialize();
      Set_Properties(ze_cryst, x_dim, y_dim, z_dim, x_sur, y_sur, z_sur, 
		             angle, stuffed, vary_size, neut_surf);                         
      }
      
nanosurface::~nanosurface()
      {
      }

void nanosurface::Initialize()
	  {
      surface_atoms.Set_Initial_Allocation(INIT_SATOM_COUNT);
	  num_core_atoms = 0;
      term_surface_indices[0] = 1;
      term_surface_indices[1] = 0;
      term_surface_indices[2] = 0;   
	  surface_angle = 0.0;
	  dimensions[0] = 0.0;
      dimensions[1] = 0.0;
      dimensions[2] = 0.0;
      stuff_space = TRUEV;
	  variable_size = FALSEV;
	  neutralize_by_surface = FALSEV;
	  }

//Properties setting---

void nanosurface::Set_Properties(const crystal_system& ze_cryst, double x_dim, double y_dim, 
                                 double z_dim, int x_sur, int y_sur, int z_sur, double surf_angle,
								 bool stuffed, bool size_vary, bool neut_surf)
      {
      cryst_struct = ze_cryst;
      dimensions[0] = x_dim;
      dimensions[1] = y_dim;
      dimensions[2] = z_dim;
      term_surface_indices[0] = x_sur;
      term_surface_indices[1] = y_sur;
      term_surface_indices[2] = z_sur;   
	  surface_angle = surf_angle;
      stuff_space = stuffed;
	  variable_size = size_vary;
	  neutralize_by_surface = neut_surf;
      }  
 
//Nanosurface construction--- 
     
void nanosurface::Construct_Surface()  
      {
      /*Surface construction follows this algorithm:
       1) Rotate the unit cell as to have the desired plane in proper Cartesian orientation
       2) Search through coordinate space for all lattice points that will be
          in the nanosurface volume
       3) Define the nanosurface through addition of crystal basis atoms
       4) Undo the unit cell rotation
	   5) Charge balance the surface, if requested
	   6) Position the surface corner at the (0, 0, 0) position
      */
      
      //Do the unit cell rotation---
      
      const int NEW_PLANE_INDICES[3] = { 0, 0, -1 };
	   //Terminating surface plane of the nanosurface is, by default,
	   //the upper plane of the orthorhombic surface box. The terminating
	   //crystal plane desired must reside in this plane position.
       
      bool h_check = (NEW_PLANE_INDICES[0] != term_surface_indices[0]);
      bool k_check = (NEW_PLANE_INDICES[1] != term_surface_indices[1]);
      bool l_check = (NEW_PLANE_INDICES[2] != term_surface_indices[2]); 
      bool do_rotation = (h_check || k_check || l_check);
       //If the user requests a (00-1) terminating surface, no rotation is necessary.
      if (do_rotation)  
         {
         cryst_struct.Rotate_Planes(NEW_PLANE_INDICES, term_surface_indices);
          //Rotate the unit cell so that the nanosurface terminating plane is along
          //the (00-1) plane of the original Cartesian system.
         }
	  double surface_normal[3] = { 0.0, 0.0, -1.0 };
	  if (!Check_FP_Equality(surface_angle, 0.0))
		//Final rotation about the z-axis of the crystal system.
	     {
		 cryst_struct.Rotate_Cell(surface_normal, surface_angle);
	     }
         
	  //Determine the points in three-dimensional space to place crystal
	  //basis atoms at---

	  double** points;
	  int point_index = 0;
	  double volume = dimensions[0]*dimensions[1]*dimensions[2];
	  double mem_scaling_factor = 10000.0/8.0;
	   //Assume that a 8 nm^3 surface would not require more than 10000 lattice points.
	  int max_points = int(mem_scaling_factor * (volume + 1.0));
	   //Scale that factor to the actual volume of the surface to be created.
	  Get_Memory(points, max_points, 3);
	   //Reserve memory for lattice point sampling.

	  double lattice_vecs[3][3];
	  cryst_struct.Get_Vectors(lattice_vecs[0], lattice_vecs[1], lattice_vecs[2]);
	  double largest_vector = cryst_struct.Biggest_Step();
	  double extra_space[3] = { 4.0*largest_vector, 4.0*largest_vector, 4.0*largest_vector };
	   //Extra space in spatial sampling parameters ensures that edges of the orthorhombic box are 
	   //not missed out upon.

	  double sampling_box_size[3];
	  SetAdd_XYZ(sampling_box_size, dimensions, extra_space);
	   //Get the region within which to obtain a point sampling = nanosurface dimensions + extra space.
	  Get_OrthoRhombic_Sampling(lattice_vecs[0], lattice_vecs[1], lattice_vecs[2], 
		                        sampling_box_size, points, max_points, point_index); 

	  //Shift those points such that: (0, 0, 0) origin ---> box center. This will work much
	  //better with the algorithm to come that ensures an orthorhombic box is created---

	  for (int a = 0; a < point_index; ++a)
	      {
          points[a][0] += dimensions[0]/2.0;
		  points[a][1] += dimensions[1]/2.0;
		  points[a][2] += dimensions[2]/2.0;
	      }

	  //Add the basis atoms to those points---

	  for (int a = 0; a < point_index; ++a)
	      {
	      cryst_struct.Add_Basis_ToSurface(surface_atoms, points[a],
                                           dimensions, stuff_space);
	      }

	  Free_Memory(points, max_points);

      //Undo the cell rotation---

	  if (surface_angle != 0.0)
	     {
		 cryst_struct.Rotate_Cell(surface_normal, -1.0*surface_angle);
		  //Undo the z-axis rotation.
	     }
      if (do_rotation)
         {
         cryst_struct.Rotate_Planes(term_surface_indices, NEW_PLANE_INDICES);
          //Undo the unit cell plane A --> plane B rotation. 
         }

	  //Charge neutralization, if requested---

	  double total_charge = surface_atoms.Return_Total_Charge();
	  if ( (variable_size || neutralize_by_surface) && !Check_FP_Equality(total_charge, 0.0))
		  {
		  Neutralize_Surface();
		  }   
      
      //Set the corner of the nanosurface to be exactly at the origin---
      
      double INITIAL_CORNER[3] = { 0.0, 0.0, 0.0 };
      ReAssign_Corner(INITIAL_CORNER);

	  Show_Statement("Done creating nanosurface of atomic size: ", surface_atoms.Size(), "!");
      }

void nanosurface::Neutralize_Surface()
	  {	  
	  //Make sure corner is at the origin---	

      double INITIAL_CORNER[3] = { 0.0, 0.0, 0.0 };
      ReAssign_Corner(INITIAL_CORNER);

	  //Parameters for the region of the nanosurface
	  //that is trimmed during one neutralization step---

      double region_loc[3];
	   //Region for trimming.
	  double sizes[3];
	   //Size of the region for trimming.
	  double boundaries[3];
	  const double SHIFT[3] = { FP_ERROR_FIX, FP_ERROR_FIX, FP_ERROR_FIX };
	  SetAdd_XYZ(boundaries, dimensions, SHIFT);
	   //Opposite corner (trimming takes place from all sides of the surface).

	  double sys_charge;
	   //System charge, which is ideally to be made zero.
	  if (variable_size)
	   {	  
	   const int MAX_CYCLES = 5;
	    //Maximum number of thickness reductions performed
	    //before stopping.
	   const int NUM_DIMENSIONS = 3;
	    //Controls if x (1-3), y(2-3), and z(3) directions
	    //are used in size shrinking.

	   double region_backup[3];
	   double sizes_backup[3];
	    //Temporary backups.
	   double removal_slice_size;
	    //Thickness of slice that is trimmed off the surface.
	   double slice_charge;
	    //Charge of the slice removed.
	   double charge_dif_old, charge_dif_new;
	    //Comparison of charges with and without cutting.
	   int step_count;

	   //Initiate the trimming routine---

	   for (int a = 0; a < MAX_CYCLES; ++a)
		 {
		 for (int b = 0; b < NUM_DIMENSIONS; ++b)
		  //Axes of movement (x, y, z)
			{
			for (int c = 0; c < 2; ++c)
			 //Approach from the negative or positive side of the surface.
				{
				sys_charge = surface_atoms.Return_Total_Charge();
				 //Get current charge.
				Set_XYZ(region_loc, INITIAL_CORNER);
				Set_XYZ(sizes, boundaries);
				 //Initialize region and region size variables. 
				removal_slice_size = 0.0;
				 //Initialize slice thickness for cutting.
				charge_dif_old = 0.0;
				step_count = -1;
				do
			     //Increase value of thickness variation (cutting thickness)
				 //until the charge post-cutting would be closer to zero.
					{
					++step_count;
					if (step_count != 0)
					   {
					   charge_dif_old = charge_dif_new;
					   }

					Set_XYZ(region_backup, region_loc);
					Set_XYZ(sizes_backup, sizes);

					removal_slice_size += SIZE_VARI_CONST;
					sizes[b] = removal_slice_size;
					if (c == 1)
					 //Approaching from the positive side...
					   {
					   region_loc[b] = boundaries[b] - removal_slice_size;
					   }
					slice_charge = surface_atoms.Return_Total_Charge(region_loc, sizes);
					charge_dif_new = abs(sys_charge) - abs(sys_charge - slice_charge);
					}
				while (Check_FP_GreaterOrEqual(charge_dif_new, charge_dif_old));
				 //Above conditional indicates that the new trimming conditions
				 //would remove even more net charge.

				Set_XYZ(sizes, sizes_backup);
				Set_XYZ(region_loc, region_backup);
				 //Restores sizes and region location to their value before the
				 //last step of the previous loop.

				if (Check_FP_GreaterThan(charge_dif_old, 0.0))
				 //Make sure trimming would bring total charge closer to zero.
				   {
				   surface_atoms.Delete_Atoms(region_loc, sizes);
				   dimensions[b] -= (removal_slice_size - SIZE_VARI_CONST);
				   }
				}
			}	
		 }
	   }

	  //Second mode: Charge neutralization by ion removal---

	  if (neutralize_by_surface)
	   {
	   Set_XYZ(region_loc, INITIAL_CORNER);
	   Set_XYZ(sizes, boundaries);
	   sizes[2] = SIZE_VARI_CONST;
	    //Set initial box parameters. Ion removal takes
	    //place off of the top layer of the surface.

	   int try_count = 0;
	   sys_charge = surface_atoms.Return_Total_Charge();
	   double test_pos[3];
	   int test_index;
	   const int MAX_TRIES = int(200.0 * (0.01/SIZE_VARI_CONST));
	    //In the worst case, this algorithm will being increasing
	    //the depth into the surface within which it looks for
	    //atoms to remove. This value gives a max depth of 10 Angstroms.

	   while (abs(sys_charge) > 0.3)
	     {
		 Get_Rand_Box_Position(test_pos, region_loc, sizes);
		 test_index = surface_atoms.Find_Closest_Atom(test_pos);
		 ++try_count;
		 if (surface_atoms[test_index].Get_Atom_Charge()/sys_charge > 0.0)
		  //If atom is of the right charge, remove it.
			{
			sys_charge -= surface_atoms[test_index].Get_Atom_Charge();
			surface_atoms.Delete_Atom(test_index);
			try_count = 0; 
			}
		 else if (try_count > MAX_TRIES/2)
		  //After repeated failures looking for an atom to removal,
		  //start expanding the search box into the depth of the surface.
			{
			sizes[2] += SIZE_VARI_CONST;
			}
		 
		 if (try_count > MAX_TRIES)
		    {
			Show_Warning("WARNING: CHARGE REMOVAL FROM SURFACE WAS A FAILURE!");
			sys_charge = 0.0;
		    }
		 }	 
	   }

	  }


//Utility functions---

void nanosurface::Get_Corner(double* corner_coors) const
 //Obtains the "(0, 0, 0) corner" by finding the lowest x/lowest y/lowest z 
 //coordinates within the slab.
      {
	  surface_atoms.Get_Coordinate_Mins(corner_coors);           
      }
 
void nanosurface::Get_Corner_Atom_Loc(double* atom_coors) const
     {
     double corner_coors[3];
     Get_Corner(corner_coors);

     Zero_XYZ(atom_coors);
     double temp_coors[3];
     double corner_distance;
	 double least_distance;
	 if (surface_atoms.Size() != 0)
		{
        least_distance = surface_atoms[0].Get_Distance(corner_coors);
		surface_atoms[0].Get_Atom_Location(atom_coors);
	    }
     const double ALLOWED_SURF_PLANE_THICKNESS = 0.04;
	  //Allowed variation in the z-coordinate, relative to the
	  //minimum-coordinates z-value.
     for (int a = 1; a < surface_atoms.Size(); ++a)
          {
          surface_atoms[a].Get_Atom_Location(temp_coors);
          Sub_XYZ(temp_coors, corner_coors);
          if (temp_coors[2] < ALLOWED_SURF_PLANE_THICKNESS)
           //If atom is along the surface plane.
             {
		     corner_distance = Get_VecMag(temp_coors);
             if (corner_distance < least_distance)
                {
                least_distance = corner_distance;
                surface_atoms[a].Get_Atom_Location(atom_coors);
                }               
             }
          }                                  
     }

void nanosurface::Get_OpposingCorner_Atom_Loc(double* atom_coors) const
     {
     double corner_coors[3];
     Get_Corner(corner_coors);

     Zero_XYZ(atom_coors);
     double temp_coors[3];
     double corner_distance;
	 double worst_distance;
     if (surface_atoms.Size() != 0)
		{
        worst_distance = surface_atoms[0].Get_Distance(corner_coors);
		surface_atoms[0].Get_Atom_Location(atom_coors);
	    }
     const double ALLOWED_SURF_PLANE_THICKNESS = 0.04;
     for (int a = 0; a < surface_atoms.Size(); ++a)
          {
          surface_atoms[a].Get_Atom_Location(temp_coors);
          Sub_XYZ(temp_coors, corner_coors);
          if (temp_coors[2] < ALLOWED_SURF_PLANE_THICKNESS)
           //If atom is along the surface plane
             {
		     corner_distance = Get_VecMag(temp_coors);
             if (corner_distance > worst_distance)
                {
                worst_distance = corner_distance;
                surface_atoms[a].Get_Atom_Location(atom_coors);  
                }               
             }
          }                                  
     }

void nanosurface::Translate_Surface(const double* shift_coors)
      {
      surface_atoms.Shift_Atomic_Coordinates(shift_coors);                      
      }
     
void nanosurface::ReAssign_Corner(const double* new_corner_coors)
      {
	  surface_atoms.Translate_To_New_Coordinate_Mins(new_corner_coors);          
      }

void nanosurface::Staircase_Shift(double shift_space, double shift_val)
      {
	  double loc[3];
	  double corner[3];
	  Get_Corner(corner);

	  double far_corner[3];
	  Set_XYZ(far_corner, corner);
	  Add_XYZ(far_corner, dimensions);
	  const double TOLERANCE = 0.0;
	   //How far out of the orthorhombic surface box an atom can be pushed
	   //before it is moved back into the box.

	  double column_ratio, total_shift;
	  for (int a = 0; a < surface_atoms.Size(); ++a)
	      {
		  surface_atoms[a].Get_Atom_Location(loc);
		  column_ratio = (loc[2]-corner[2])/shift_space;
		  total_shift = column_ratio*shift_val;
		  for (int b = 0; b < 2; ++b)
		      {
			  loc[b] += total_shift;
			  while (loc[b] > (far_corner[b] + TOLERANCE))
		         {
			     loc[b] -= (far_corner[b] - corner[b]);
			     if (loc[b] < far_corner[b])
			      //Try to have flow to the other side not produce
			      //a local density increase.
			       {
				   loc[b] -= TOLERANCE;
			       }
		         }
		      }
		  surface_atoms[a].Set_Atom_Location(loc);
	      }
      }

void nanosurface::Random_Shift(const double* rand_parameters)
      {
	  double corner[3];
	  Get_Corner(corner);

	  double far_corner[3];
	  Set_XYZ(far_corner, corner);
	  Add_XYZ(far_corner, dimensions); 

	  double loc[3];
	  double rand_shifts[3];
	  Seed_RandNums();
	  for (int a = 0; a < surface_atoms.Size(); ++a)
	      {
		  surface_atoms[a].Get_Atom_Location(loc);
		  for (int b = 0; b < 3; ++b)
		      {
			  rand_shifts[b] = GetRandDec(-1.0*rand_parameters[b], rand_parameters[b]);
			  loc[b] += rand_shifts[b];
			  if (loc[b] > far_corner[b])
			     {
				 loc[b] -= (far_corner[b] - corner[b]);
			     }
			  else if (loc[b] < corner[b])
			     {
				 loc[b] += (far_corner[b] - corner[b]);
			     }
		      }
		  surface_atoms[a].Set_Atom_Location(loc);
	      }
      }

void nanosurface::Slice(double slice_thickness, double pot_offset)
	  {
	  double loc[3];
	  double z_temp;
	  int slice_index;
	  double original_origin[3];
	  Get_Corner(original_origin);

	  for (int a = 0; a < surface_atoms.Size(); ++a)
	      {
		  surface_atoms[a].Get_Atom_Location(loc);
		  z_temp = (loc[2] + pot_offset)/slice_thickness;
		  slice_index = int(z_temp);
		  if (slice_index < 0) slice_index = 0;
		  loc[2] = original_origin[2] + 
			       (double(slice_index) + 0.5)*slice_thickness - pot_offset;
		  if (loc[2] < original_origin[2])
		   //Prevent atoms from being placed outside the box, in the event
		   //this function is used with a large offset parameter.
			{
			loc[2] = original_origin[2] + 1.5*slice_thickness - pot_offset;
			}
		  surface_atoms[a].Set_Atom_Location(loc);
	      }
	  atom fake_atom;
	  fake_atom.Set_Atom_Location(original_origin);
	  surface_atoms.Add_Atom(fake_atom);
	   //Dummy atom put at original origin.
	   //Prevents placement logic in composite from
	   //being affected by 2D slicing.
	  }

void nanosurface::Rotate_Along_Norm(double angle)
	  {
	  double rot_vec[3] = { 0.0, 0.0, -1.0 };
	  surface_atoms.Rotate_Atoms(rot_vec, angle);
	  }

void nanosurface::Get_XY_Repeat(double& x_repeat, double& y_repeat) const
      {
      double corner_coors[3];
      Get_Corner_Atom_Loc(corner_coors);
      
	  double oppcorner_coors[3];
	  Get_OpposingCorner_Atom_Loc(oppcorner_coors);
      x_repeat = oppcorner_coors[0] - corner_coors[0];
      y_repeat = oppcorner_coors[1] - corner_coors[1];
      
	  //Look for atoms that are connected to the corner atom
	  //by an x or y Cartesian axis. The closest atoms thus
	  //define the x/y repeat distances---

      double temp_difs[3];
      for (int a = 0; a < surface_atoms.Size(); ++a)
          {
          surface_atoms[a].Get_Atom_Location(temp_difs);
          Sub_XYZ(temp_difs, corner_coors);
          if (Check_FP_Equality(temp_difs[1], 0.0) && Check_FP_Equality(temp_difs[2], 0.0))
		  //Check to see if the atoms are connected via the x-axis.
            {
            if (( temp_difs[0] < x_repeat ) && !Check_FP_Equality(temp_difs[0], 0.0))
		     //If TRUEV, function has found a closer atom.
              {
              x_repeat = temp_difs[0];       
              }                                 
            }
          if (Check_FP_Equality(temp_difs[0], 0.0) && Check_FP_Equality(temp_difs[2], 0.0))
		   //Check to see if the atoms are connected via the y-axis.
            {
            if (( temp_difs[1] < y_repeat ) && !Check_FP_Equality(temp_difs[1], 0.0))
              {
              y_repeat = temp_difs[1];       
              }                                 
            }        
          }                            
      }

//Coating---

void nanosurface::Coat_Surface(const molecule& coat_mol, double sh_dist, 
                               double hh_dist_x, double hh_dist_y, 
                               bool place_on_atom)
     {
     layer_molecule = coat_mol;
      //Record the coating molecule. 
     layer_molecule.Standard_Orientation();
	  //Orient the molecular backbone along the surface normal.

     num_core_atoms = surface_atoms.Size();
      //Record the number of atoms in the surface before addition of
	  //monolayer.

	 //Initial spatial analysis for monolayer placement---
     
     double surf_corner[3];
     double opp_corner[3];
     Get_Corner_Atom_Loc(surf_corner); 
     Get_OpposingCorner_Atom_Loc(opp_corner);
      //Get corners of the "sheet" upon which the monolayer is
	  //to be placed.

     double spread[3];
     Set_XYZ(spread, opp_corner);
     Sub_XYZ(spread, surf_corner);
      //Get x/y length of that sheet.

     bool no_x_define = Check_FP_Equality(hh_dist_x, 0.0); 
     bool no_y_define = Check_FP_Equality(hh_dist_y, 0.0);
	  //Values of 0 from user = Default = determine monolayer
	  //molecule-molecule spacing.
     if (no_x_define && no_y_define)
        {
        Get_XY_Repeat(hh_dist_x, hh_dist_y);                              
        }
     else if (no_x_define || no_y_define)
      //If one dimension only is undefined, set the monolayer to be just
      //a row of molecules along one dimension.
        {                        
        if (no_x_define)     
           {
           hh_dist_x = spread[0]; 
		    //Molecule-molecule spacing = spatial length of surface sheet
		    //will result in only one row of molecules along the other direction.
           }
        if (no_y_define)
           {
           hh_dist_y = spread[1];             
           }
        }    

     double step_pams[2];
     step_pams[0] = hh_dist_x;
     step_pams[1] = hh_dist_y;   
     double start_vals[3];
     double end_vals[3];
     Get_Spatial_Sampling_Pams(surf_corner, opp_corner, step_pams, 2, start_vals, end_vals);
	  //Obtain the proper starting location for molecule placement such that the molecules
	  //look evenly placed across the whole sheet.
     
     double atom_pos[3];
     int molecule_size;
     int nearest_atom_index;
     for (double a = start_vals[0]; a < end_vals[0]; a += step_pams[0])
         {
         for (double b = start_vals[1]; b < end_vals[1]; b += step_pams[1])
             {
             atom_pos[0] = a;
             atom_pos[1] = b;
             atom_pos[2] = surf_corner[2] - sh_dist;
             if (place_on_atom)
                {
                nearest_atom_index = surface_atoms.Find_Closest_Atom(0, num_core_atoms, atom_pos);
                surface_atoms[nearest_atom_index].Get_Atom_Location(atom_pos);
                atom_pos[2] -= sh_dist;
                }
             molecule_size = layer_molecule.Copy_Molecule(surface_atoms, atom_pos, FIRST_ATOM);
			 int mol_start = surface_atoms.Size() - molecule_size;
			 int mol_end = surface_atoms.Size();
			 surface_atoms.Monolayer_Overlap_Prevention(num_core_atoms, mol_start, mol_start, mol_end, PI_CONST/18.0, 36);
              //Spatial overlap check. If molecule-molecule overlap can not be stopped, both molecules will be removed.
             }     
         } 
     }

//Nanosurface copying---

int nanosurface::Copy_NanoSurface(atom_collection& atomic_array) const
   {
   atomic_array.Copy_Atomic_Group(surface_atoms);       
   return surface_atoms.Size();                       
   }
   
int nanosurface::Copy_NanoSurface(atom_collection& atomic_array,
                                  const double* position) const
   {
   double shift_coors[3];
   double corner[3];
   Get_Corner(corner);
   SetSub_XYZ(shift_coors, position, corner);
   atomic_array.Copy_Atomic_Group(surface_atoms, atomic_array.Size(), 0, 
	                              surface_atoms.Size(), shift_coors);    
   return surface_atoms.Size();                     
   }

int nanosurface::Copy_NanoSurface(atom_collection& atomic_array,
                                  const double* position, const double* rot_angles) const
   {
   double shift_coors[3];
   double corner[3];
   Get_Corner(corner);
   SetSub_XYZ(shift_coors, position, corner);
   atomic_array.Copy_Atomic_Group(surface_atoms, atomic_array.Size(), 0, 
	                              surface_atoms.Size(), shift_coors, rot_angles); 
   return surface_atoms.Size();                     
   }

int nanosurface::Copy_NanoSurface(atom_collection& atomic_array,
                                  const double* position, const double* rot_angles, bool stop_hydro) const
   {
   double shift_coors[3];
   double corner[3];
   Get_Corner(corner);
   SetSub_XYZ(shift_coors, position, corner);
   atomic_array.Copy_Atomic_Group(surface_atoms, atomic_array.Size(), 0, 
	                              surface_atoms.Size(), shift_coors, rot_angles, stop_hydro); 
   return surface_atoms.Size();                     
   }

//I/O---

void nanosurface::Print_Atom_Locations() const
      {
	  cout << endl << "HERE COMES THE MIGHTY FLOW OF INFORMATION!!!" << endl;
      surface_atoms.Print_Atoms_Location();
      }

void nanosurface::Store_Atom_Locations(ofstream& out_file) const
   {                                              
   surface_atoms.Store_Atoms_Location(out_file, TRUEV);
   }
   
void nanosurface::Save_NanoSurface(const char* save_file) const
   {
   ofstream ze_file(save_file, ios::out);
   if (ze_file.is_open())
      {
      Write_Intro(ze_file, "File Type = NanoSurface Description");
      Write_Values(ze_file, term_surface_indices, 3);
      ze_file << endl;
      Write_Values(ze_file, dimensions, 3);
      ze_file << endl << stuff_space << ' ' << variable_size << ' ' 
		      << neutralize_by_surface << ' ' << surface_angle << endl;
      cryst_struct.Save_Crystal_System(ze_file);
      ze_file << endl << endl << "Atomic Information List";
      surface_atoms.Atoms_Storage(ze_file);  
	  ze_file.close();
      }    
   }
   
void nanosurface::Load_NanoSurface(const char* load_file)
   {
   ifstream ze_file(load_file, ios::in);
   if (ze_file.is_open())
      {
      Skip_Phrases(ze_file, 5);
      Load_Values(ze_file, term_surface_indices, 3);
      Load_Values(ze_file, dimensions, 3);
      ze_file >> stuff_space >> variable_size >> neutralize_by_surface >> surface_angle;
      cryst_struct.Load_Crystal_System(ze_file);
      Skip_Phrases(ze_file, 3);
      surface_atoms.Atoms_Retrieval(ze_file);
	  ze_file.close();
      }
   }