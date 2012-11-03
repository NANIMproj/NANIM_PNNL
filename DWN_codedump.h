//Point sampling bins---


	 //Free_Memory(bins, total_bin_count);
	 //Free_Memory(points_placed, total_bin_count);

			 /*  bin_index = try_pos[0]/length_1D 
				         + (try_pos[1]/length_1D)*bin_count_1D 
				         + (try_pos[2]/length_1D)*bin_count_2D;
			   for (int a = 0; a < points_per_volume; ++a)
			      {
                  if (points_placed[bin_index][a] == FALSEV)
				     {
					 bins[bin_index][a] = point_index;
                     points_placed[bin_index][a] = TRUEV;
				     }
			      }*/

	/* double volume = (4.0/3.0)*PI_CONST*(radius*radius*radius);
	 int points_per_volume = int(double(max_points)/volume).
      //Number of points allowed per unit volume.
	 double length_1D = pow(volume, 1.0/3.0);
     int bin_count_1D = int(length_1D) + 2;
	 int bin_count_2D = bin_count_1D*bin_count_1D;
	 int total_bin_count = bin_count_1D*bin_count_1D*bin_count_1D;
	 int** bins;
	 bool** points_placed;
     Get_Memory(bins, total_bin_count, points_per_volume);
	 Get_Memory(bins, total_bin_count, points_per_volume); 
	  //Bin location is calculated as follows:
	  //Bin num = int(x_cor/length_1D) + int(y_cor/length_1D)*bin_count_1D 
	  //        + int(z_cor/length_1D)*bin_count_2D
	  //For sphere, ratio values must be added to half the bin_count_1D
	  //Neighboring bins = bin num +/- 1, bin_count_1D, bin_count_2D
	  //Safety: If bun count < 0 or > max, set equal to 0. */

//Rotations for plane sampling---
//BCC_GRID ~0.953 rotation angle
	   // rotation_factor = PI_CONST/4.0;
		 //Rotate to either 110 or 1-10 direction, depending on the plane.
		//if ( (ind[1] + ind[0]) == -2)
	      //Correction factor.
		 // {
         // step_arrowA[0] = step_arrowB[0];
		 // step_arrowA[1] = -1.0*step_arrowB[1];
		 // }

/*bool Check_Boundaries_PBC(const double* coor, const double* min, const double* max, 
					      int periodicity, const double* box_dimensions)
	 {
	 bool in_bounds = FALSEV;
	 for (int a = 0; a < 3; ++a)
		 {
		 
		 }
	 return in_bounds;
	 }
	 */
//Bin logic---

/*void Get_Bin_Memory(int**& bins, bool**& points_placed, double volume, int total_point_count, 
	                double& 1D_length, int& bin_count_1D, int& bin_count_2D, 
					int& total_bin_count, int& points_per_volume)
   {
   int points_per_volume = int(double(total_point_count)/volume).
      //Number of points allowed per unit volume.
   1D_length = pow(volume, 1.0/3.0);
   bin_count_1D = int(length_1D) + 2;
   bin_count_2D = bin_count_1D*bin_count_1D;
   total_bin_count = bin_count_1D*bin_count_1D*bin_count_1D;
   Get_Memory(bins, total_bin_count, points_per_volume);
   Get_Memory(bins, total_bin_count, points_per_volume); 
   }*/


//Bin location is calculated as follows:
//Bin num = int(x_cor/length_1D) + int(y_cor/length_1D)*bin_count_1D 
//        + int(z_cor/length_1D)*bin_count_2D
//For sphere, ratio values must be added to half the bin_count_1D
//Neighboring bins = bin num +/- 1, bin_count_1D, bin_count_2D
//Safety: If bun count < 0 or > max, set equal to 0. 

/*void Release_Bin_Memory(int**& bins, bool**& points_placed, int total_bin_count) 
   {
   Free_Memory(bins, total_bin_count);
   Free_Memory(points_placed, total_bin_count);
   }
*/		


/*Problem solvers---


//Problem solvers developed along the way---

void Coordinate_Solver();
 //A highly specific function that is used only by me to help figure out coordinates
 //of atoms in big unit cells.
void Solver_Function(const double*, const char*); 
 //Solves a coordinate function in terms of the 3 passed parameters, storing the
 //six results in the second array passed.
 
void Rotation_Problem(const double*, const double*, const double*, double*);
 //Rotates set B into set C for point D with origin A.
void Rotation_Problem2(); 

//Specific functions---

void Coordinate_Solver()
     {
     cout << endl << "AtomR N 0.0 0.0 0.0" 
          << endl << "AtomR N 0.0 0.0 0.5000"
          << endl << "AtomR N 0.3333 0.6666 0.7500"
          << endl << "AtomR N 0.6666 0.3333 1.2500";
     
     double pam[3];
     
     pam[0] = 0.33333;
     pam[1] = 0.000;
     pam[2] = 0.000;
     char nitro[3] = "N";
     Solver_Function(pam, nitro);
     
     pam[0] = 0.33333;
     pam[1] = 0.33333;
     pam[2] = 0.25;
     Solver_Function(pam, nitro);
     
     pam[0] = 0.5;
     pam[1] = 1.0/12.0;
     pam[2] = 0.25;
     char silicon[4] = "Si";
     Solver_Function(pam, silicon);
     
     pam[0] = 1.0/6.0;
     pam[1] = 0.25;
     pam[2] = 0.00;
     Solver_Function(pam, silicon);
     }

void Solver_Function(const double* pams, const char* atom_name)
     {
     double results[6][3];
     results[0][0] = pams[0];
     results[0][1] = pams[1];
     results[0][2] = pams[2];
     
     results[1][0] = -1.0 * pams[1];
     results[1][1] = pams[0] - pams[1];
     results[1][2] = pams[2];
     
     results[2][0] = pams[1] - pams[0];
     results[2][1] = -1.0 * pams[0];
     results[2][2] = pams[2];
     
     results[3][0] = pams[1];
     results[3][1] = pams[0];
     results[3][2] = 0.5 + pams[2];
     
     results[4][0] = -1.0 * pams[0];
     results[4][1] = pams[1] - pams[0];
     results[4][2] = 0.5 + pams[2];
     
     results[5][0] = pams[0]-pams[1];
     results[5][1] = -1.0 * pams[1];
     results[5][2] = 0.5 + pams[2];
     
     for (int a = 0; a < 6; ++a)
         {
         cout << endl << "AtomR " << atom_name << " ";
         for (int b = 0; b < 3; ++b)
             {
             Check_FP_Zero(results[a][b]); 
             cout << results[a][b] << " ";                       
             }
         }
     }
     
void Rotation_Problem(const double* A, const double* B, const double* C, double* D)
     {
     double rotation_vec[3];
     double angle;
     Calculate_Rotation_Parameters(B, C, rotation_vec, angle);
     double rot_matrix[9];
     Calc_Rotation_Matrix(rotation_vec, angle, rot_matrix);
     rot_matrix[0] *= 1.0;
     rot_matrix[1] *= 1.0;
     rot_matrix[2] *= 1.0;
     rot_matrix[3] *= 1.0;
     rot_matrix[4] *= 1.0;
     rot_matrix[5] *= 1.0;
     rot_matrix[6] *= 1.0;
     rot_matrix[7] *= 1.0;
     rot_matrix[8] *= 1.0;
     Set_XYZ(D, B);
     cout << endl << endl << "START " << D[0] << " " << D[1] << " " << D[2];
     Rotate_Vector(rot_matrix, D);        
     cout << endl << "END " << D[0] << " " << D[1] << " " << D[2];  
     cout << endl << "DESIRED " << C[0] << " " << C[1] << " " << C[2] << endl;                        
     }
     
void Rotation_Problem2()
     {
     double rot_vec[3] = { 0.0, 2.0, 2.0 };
     double vec[3] = { 1.0, 1.0, -1.0 }; 
     double angle = 3.141/2.0;
     double rot_matrix[9];
     Calc_Rotation_Matrix(rot_vec, angle, rot_matrix);
     Rotate_Vector(rot_matrix, vec);
     
     cout << endl << endl << endl << vec[0] << " " << vec[1] << " " << vec[2] << endl << endl;  
     
     }*/

	 	   // char sec_file_name[50] = "BIG_TEST.cfg";
	  /*double* test_vec[3] = { 1.0, 0.0, 0.0 };
	  double* rot_vec[3] = { -1.0, 2.0, 0.0 };
	  double* rot_vec_b[3] = { 1.0, -2.0, 0.0 };
	  rot_angle = PI_CONST/2.0;
	  Rotate_Vector(test_vec, rot_vec, rot_angle);*/

	/* char test_script[20] = "Hehe 500  ";
	  char test_FPscript[20] = "Hm 16.09 ";
	  char test_script2[40] = "Integers 500 71 9004 ";
	  char test_FPscript2[40] = "FPS  500.78 20.0078 19.05";
	  int change_vals[3] = { 178, 5, 703 };
	  double changeFP_vals[3] = {100.00, 4.11, 0.007};
	  Change_Integer(test_script, 215);
	  Change_FP(test_FPscript, 5.0);
	  Change_Integers(test_script2, change_vals, 3);
	  Change_FPS(test_FPscript2, changeFP_vals, 3);*/


//{111} faceting planes for cubic-rounded---

                         //const double f3 = 1.73185;
						 //Slightly less than sqrt(3).
						//const double f3 = 2.0;
						 //Currently have this turned off via a large value.

						/*new_crystal_systems[struct_index].Set_Faceting_Plane(1, 1, 1, f3);
						new_crystal_systems[struct_index].Set_Faceting_Plane(-1, 1, 1, f3);
						new_crystal_systems[struct_index].Set_Faceting_Plane(-1, -1, 1, f3);
						new_crystal_systems[struct_index].Set_Faceting_Plane(1, -1, 1, f3);
						new_crystal_systems[struct_index].Set_Faceting_Plane(1, 1, -1, f3);
						new_crystal_systems[struct_index].Set_Faceting_Plane(-1, 1, -1, f3);
						new_crystal_systems[struct_index].Set_Faceting_Plane(-1, -1, -1, f3);
						new_crystal_systems[struct_index].Set_Faceting_Plane(1, -1, -1, f3);*/
					      //{111} family.

//Nanosurface code---

    //  vector<atom> slab_atoms(INITIAL_SATOM_COUNT);
       /*Initial list of positions to place the crystal basis set at. 
       For example, with bcc, two atoms would be placed for each coordinate
       set in this list*/
      /*int atom_index = 0;
        //Index for the atom list.     

      double lattice_vec_mags[3];
      cryst_struct.Get_Vectors_Mag(lattice_vec_mags);
       //Get crystal system information.

      double largest_dim = Get_Largest(dimensions[0], dimensions[1], dimensions[2]);

      const double search_factor = 2.0;
      int max_a_steps = int(big_enough_factor * largest_dim/lattice_vec_mags[0]);
      int max_b_steps = int(big_enough_factor * largest_dim/lattice_vec_mags[1]);
      int max_c_steps = int(big_enough_factor * largest_dim/lattice_vec_mags[2]);*/
       /*Maximum steps along the crystal system vectors needed to make a big 
       enough slab to cut the surface out of. Based on largest slab dimension
       and the BEF factor must be fairly large to ensure that the desired slab
       can be found within the greater slab initially made.*/
      

  /*  
    for (int a = 0; a < max_a_steps; ++a)
          {
          for (int b = 0; b < max_b_steps; ++b)
              {
              for (int c = 0; c < max_c_steps; ++c)
                  {  
                  const int new_memory_size = 1000;
                  if ( slab_atoms.size() < (atom_index + new_memory_size) )
                   //Make sure memory is available for new atoms. 
                     {
                     slab_atoms.resize(atom_index + new_memory_size);                  
                     }       
                  
                  slab_atoms[atom_index].Set_Atom_Location(slab_x, slab_y, slab_z);
                  cryst_struct.Add_Basis(slab_atoms, atom_index);
                  cryst_struct.C_Step(slab_x, slab_y, slab_z);
                   //Incrementation of coordinates along the c vector by one step.
                  }
              cryst_struct.B_Step(slab_x, slab_y, slab_z);
               //Incrementation of coordinates along the b vector.
              cryst_struct.C_Step(slab_x, slab_y, slab_z, -1*max_c_steps);
               //Reset any incrementation from the c vector from the previous vector.   
              } 
          cryst_struct.A_Step(slab_x, slab_y, slab_z);
          cryst_struct.B_Step(slab_x, slab_y, slab_z, -1*max_b_steps);      
          }   

      
      //Adjust relative coordinates---
      
      double x_min = 1000.0;
      double y_min = 1000.0;
      double xy_min = 1000.0; 
      double z_min = 1000.0;
       //Initialize to arbitrary and large values.
      double x_temp, y_temp, z_temp, xy_temp, dif_temp;
      
      for (int a = 0; a < atom_index; ++a)
       //Find lowest z-coordinate, which should be the approximate z-coordinate 
       //location of the top surface plane on the nanosurface.
          {  
          slab_atoms[a].Get_Atom_Location(x_temp, y_temp, z_temp); 
          if (z_temp < z_min)
             {
             z_min = z_temp;                           
             }
          }
      
      for (int a = 0; a < atom_index; ++a)
       //Find low xy-coordinate corner of the top surface plane.
          {
          slab_atoms[a].Get_Atom_Location(x_temp, y_temp, z_temp);    
          dif_temp = z_temp - z_min;
          const double rotation_error_bar = 0.05;
          if (abs(dif_temp) < rotation_error_bar)
           //Looking for atoms in the top surface plane while leaving some
           //room for error in the ability of the previous rotation to make
           //that surface truly flat.
             {
             xy_temp = x_temp + y_temp;
             if (xy_temp < xy_min)
                {
                xy_min = xy_temp;
                x_min = x_temp;
                y_min = y_temp;         
                }       
             }
          }
          
       double neg_x_min = -1.0*x_min;
       double neg_y_min = -1.0*y_min;
       double neg_z_min = -1.0*z_min;
       for (int a = 0; a < atom_index; ++a)
        //Finally, subtract the x, y, and z minimum coordinates from all atoms.
           {
           slab_atoms[a].AddTo_Atom_Location(neg_x_min, neg_y_min, neg_z_min);      
           }
      
      
      //Trim surface slab down to requested sizes---
      //Also: Define the final array of atoms to do further calculation with.
  
      for (int a = 0; a < atom_index; ++a)
          {    
          slab_atoms[a].Get_Atom_Location(x_temp, y_temp, z_temp);
          if ( (x_temp < dimensions[0]) && (y_temp < dimensions[1])
                && (z_temp < dimensions[2]) )
           //If within the requested dimensions of the nanosurface, include the
           //atom in the final array.     
             {
             surface_atoms[num_atoms] = slab_atoms[a];
             ++num_atoms;
             }
          }*/



    /*      need_to_skip = FALSEV;
                  for (int a = 0; a < index_temp; ++a)
                      {
                      x_check = Check_FP_Equality(slab_coors[0], memory_visits[a][0]);  
                      y_check = Check_FP_Equality(slab_coors[1], memory_visits[a][1]); 
                      z_check = Check_FP_Equality(slab_coors[2], memory_visits[a][2]);
                      if (x_check && y_check && z_check)
                       //Found a match.
                        {
                        need_to_skip = TRUEV;
                        a = index_temp;        
                        }     
                      }
                  if (need_to_skip)
                     { 
                     continue;
                     } */


//Composite based storage calls---

/*void

general_composite::Store_Atom_Locations(ofstream& out_file, 
                                             bool include_hydro) const
     {
     atoms.Store_Atoms_Location(out_file, include_hydro);
     }

void general_composite::Store_Atom_Locations_With_Charge(ofstream& out_file, 
                                                         bool include_hydro) const
     {
     atoms.Store_Atoms_Location_With_Charge(out_file, include_hydro);                                                  
     }

void general_composite::Store_Atom_Locations_With_CoreShell(ofstream& out_file, 
                                                            bool include_hydro,
															bool cat_coreshell,
															bool ani_coreshell) const
     {
     atoms.Store_Atoms_Location_With_CoreShell(out_file, include_hydro, 
		                                       cat_coreshell, ani_coreshell);                                                  
     }

void general_composite::Store_Atom_Types_With_Masses(ofstream& out_file, 
	                                                 bool include_hydro) const
     {
     atoms.Store_Atom_Types_With_Masses(out_file, include_hydro);
     }*/
/* Molecular dynamics code, to be released when the time is ready---

//Potentials---

void Calc_Buck_Pot(DW_data& energy, const DW_data dist, const DW_data A, 
	               const DW_data B, const DW_data C)
	{
	DW_data dist_const1 = -1.0*dist/B;
	 //Negative ratio of the distance to the constant B.
	DW_data dist_const2 = (1.0/dist);
	dist_const2 = (dist_const2 * dist_const2 * dist_const2);
	dist_const2 = (dist_const2 * dist_const2);
	 //= 1/r^6
	energy = A*exp(dist_const1) - C*dist_const2;
	}

void Calc_126LJ_Pot(DW_data& energy, const DW_data dist, const DW_data A, const DW_data B)
	{
	DW_data dist_const1 = (1.0/dist);
	dist_const1 = (dist_const1 * dist_const1 * dist_const1);
	dist_const1 = (dist_const1 * dist_const1);
	 // = 1/r^6
	DW_data dist_const2 = (dist_const1 * dist_const1);
	 // = 1/r^12
	energy = A*dist_const2 - B*dist_const1;
	}

void Calc_96LJ_Pot(DW_data& energy, const DW_data dist, const DW_data A, const DW_data B)
	{
	DW_data dist_const = (1.0/dist);
	dist_const = (dist_const * dist_const * dist_const);
	 // = 1/r^3
	DW_data dist_const1 = (dist_const * dist_const);
	 // = 1/r^6
	DW_data dist_const2 = (dist_const * dist_const * dist_const);
	 // = 1/r^9
	energy = A*dist_const2 - B*dist_const1;
	}

void Calc_Morse_Potential(DW_data& energy, const DW_data dist, const DW_data De,
		                  const DW_data alpha, const DW_data equi_dist)
	{
	DW_data exp_term = -1.0*alpha*(dist - equi_dist);
	DW_data rel_term = 1.0 - exp(exp_term);
	energy = De*rel_term*rel_term;
	}

void Calc_Stretch_Potential(DW_data& energy, const DW_data dist, 
	                        const DW_data force_const, const DW_data equi_dist)
	{
	DW_data bond_distort = dist - equi_dist;
	energy = 0.5*force_const*bond_distort*bond_distort;
	}

void Calc_Bend_Potential(DW_data& energy, const DW_data angle, 
	                     const DW_data force_const, const DW_data equi_angle)
	{
	DW_data ang_distort = angle - equi_angle;
	energy = 0.5*force_const*ang_distort*ang_distort;
	}

void Calc_Torsion_Potential(DW_data& energy, const DW_data angle, const DW_data angle0,
							const DW_data a1, const DW_data a2, const DW_data a3)
	{
	DW_data term1 = a1*(1 + cos(angle - angle0));
	DW_data term2 = a2*(1 + cos(2.0*angle - angle0));
	DW_data term3 = a3*(1 + cos(3.0*angle - angle0));
	energy = term1 + term2 + term3;
	}				  

void Calc_CoreShell_Potential(DW_data& energy, const DW_data cs_dist,
							  const DW_data force_const1, const DW_data force_const2) 
	{
	DW_data dist_square = cs_dist*cs_dist;
	DW_data term1 = 0.5*force_const1*dist_square;
	DW_data term2 = (1.0/24.0)*force_const2*dist_square*dist_square;
	energy = term1 + term2;
	}

//Forces---

//Note: Connection vectors should be passed such that the force returned
//is expected for the 2nd atom in the definition of the connection vector.
//(I.e. for vector 12, the resulting force is on atom 2).

void Calc_Buck_Force(DW_data* force, const DW_data* con_vec, 
					 const DW_data dist, const DW_data A, const DW_data B, const DW_data C)
	{
	DW_data dist_const = 1.0/dist;
	DW_data dist_const1 = -1.0*dist/B;
	 //Negative ratio of the distance to the constant B.
	DW_data dist_const2 = (dist_const * dist_const * dist_const);
	dist_const2 = (dist_const2 * dist_const2) * dist_const;
	 //= 1/r^7
	DW_data force_scale = A/B*exp(dist_const1) - 6.0*C*dist_const2;
	Get_NewMag_Vec(con_vec, force, force_scale);
	}

void Calc_126LJ_Force(DW_data* force, const DW_data* con_vec, 
					 const DW_data dist, const DW_data A, const DW_data B)
	{
	DW_data dist_const = (1.0/dist);
	DW_data dist_const1 = (dist_const * dist_const * dist_const);
	dist_const1 = (dist_const1 * dist_const1) * dist_const;
	 // = 1/r^7
	DW_data dist_const2 = (dist_const1 * dist_const1) * dist;
	 // = 1/r^13
	DW_data force_scale = 12.0*A*dist_const2 - 6.0*B*dist_const1;
	Get_NewMag_Vec(con_vec, force, force_scale);
	}

void Calc_96LJ_Force(DW_data* force, const DW_data* con_vec, 
					 const DW_data dist, const DW_data A, const DW_data B)
	{
	DW_data dist_const = (1.0/dist);
	DW_data dist_const_r3 = (dist_const * dist_const * dist_const);
	DW_data dist_const1 = (dist_const_r3 * dist_const_r3) * dist_const;
	 // = 1/r^7
	DW_data dist_const2 = dist_const1 * dist_const_r3;
	 // = 1/r^10
	DW_data force_scale = 9.0*A*dist_const2 - 6.0*B*dist_const1;
	Get_NewMag_Vec(con_vec, force, force_scale);
	}

void Calc_Morse_Force(DW_data* force, const DW_data* con_vec, 
					 const DW_data dist, const DW_data De,
		             const DW_data alpha, const DW_data equi_dist)
	{
	DW_data exp_term1 = -1.0*alpha*(dist - equi_dist);
	DW_data exp_term2 = 2.0*exp_term1;
	DW_data al2 = 2.0*alpha;

	DW_data force_scale = De*al2*(exp(exp_term2) - exp(exp_term1));
	Get_NewMag_Vec(con_vec, force, force_scale);
	}

void Calc_Stretch_Force(DW_data* force, const DW_data* con_vec, 
					    const DW_data dist, const DW_data force_const, 
					    const DW_data equi_dist)
	{
	DW_data bond_distort = dist - equi_dist;
	DW_data force_scale = -1.0*force_const*bond_distort;
	Get_NewMag_Vec(con_vec, force, force_scale);
	}

void Calc_Bend_Potential(DW_data& energy, const DW_data angle, 
	                     const DW_data force_const, const DW_data equi_angle)
	{
	DW_data ang_distort = angle - equi_angle;
	DW_data force_scale = -1.0*force_const*ang_distort;
	}

void Calc_Torsion_Potential(DW_data& energy, const DW_data angle,
							const DW_data a1, const DW_data a2, const DW_data a3)
	{
	DW_data term1 = a1*(1 + cos(7angle));
	DW_data term2 = a2*(1 + cos(2.0*angle));
	DW_data term3 = a3*(1 + cos(3.0*angle));
	energy = term1 + term2 + term3;
	}				  

void Calc_CoreShell_Force(DW_data* force, const DW_data* con_vec, const DW_data cs_dist,
					      const DW_data force_const1, const DW_data force_const2) 
	//Determines force on the core or shell, whichever one is point 2 in the connecting
	//vector.
	{
	DW_data term1 = -1.0*force_const1*cs_dist;
	DW_data term2 = -(1.0/6.0)*force_const2*cs_dist*cs_dist*cs_dist;
	
	DW_data force_scale = term1 + term2;
	Get_NewMag_Vec(con_vec, force, force_scale);
	}


//Integrators---

void Take_Verlet_Step(DW_data* new_pos, const DW_data* cur_pos, 
		              const DW_data* old_pos, const DW_data* acceleration, 
					  const MD_param& parameters)
	{
	for (int a = 0; a < parameters.dimensionality; ++a)
		{
		new_pos[a] = (2*cur_pos[a] - old_pos[a]) + 
			       acceleration[a]*parameters.time_step*parameters.time_step;
		}
	}

void Calc_Verlet_Speed(DW_data& speed, const DW_data* new_pos, 
			           const DW_data* old_pos, const MD_param& parameters)
	{
	double cur_velocity[3];
	Zero_XYZ(cur_velocity);
	for (int a = 0; a < parameters.dimensionality; ++a)
		{
		cur_velocity[a] = (new_pos[a] - old_pos[a]) / (2*parameters.time_step);
		}
	speed = Get_VecMag(cur_velocity);
	}

void Calc_Verlet_SquareSpeed(DW_data& square_speed, const DW_data* new_pos, 
			                 const DW_data* old_pos, const MD_param& parameters)
	{
	double cur_velocity[3];
	Zero_XYZ(cur_velocity);
	for (int a = 0; a < parameters.dimensionality; ++a)
		{
		cur_velocity[a] = (new_pos[a] - old_pos[a]) / (2*parameters.time_step);
		}
	square_speed = Get_CoorSquare_Sum(cur_velocity);
	}

void Take_Euler_Step(DW_data* new_pos, const DW_data* cur_pos, 
		             const DW_data* velocity, const DW_data* acceleration, 
					 const MD_param& parameters)
	{
	for (int a = 0; a < parameters.dimensionality; ++a)
		{
		new_pos[a] = cur_pos[a] + velocity[a]*parameters.time_step + 
			         0.5*acceleration[a]*parameters.time_step*parameters.time_step;
		}
	}


//Physics---

void Calc_KinEnergy(DW_data& KE, const DW_data mass, const DW_data square_speed)
    {
	KE = 0.5*mass*square_speed;
	}

//Molecular mechanics---

void molecule_sim::Rigid_Minimization(const atom_collection&, int start, int end, bool include_hydro,
	                                  double precision)
	 {
	 double reference_axis[3];
	  //Set a standard reference by indicating the directionality of one 
	  //of the bonds in the material.

	 double step_size_rotate = 1.0 * (PI_CONST/180.0) * (PI_CONST/180.0);
	  //Change in angle per eV/radians gradient, in radians.
	 double step_size_displace = 0.001;
	  //Change in position (nm) per eV/nm gradient.
	 double max_rotation_step = 3.0 * PI_CONST/180.0;
	 double max_translation_step = 0.03;
	  //Prevent steps from getting absurdly large, just in case
	  //a huge gradient does end up being encountered.
	  //Limits = 3.0 degrees rotation (1-D), 0.3 A translation (1-D).
	 double trial_translation, trial_rotation;
	 trial_translation = 0.0001;
	 trial_rotation = 0.001 * PI_CONST/180.0;
	 //double trial_translation_max = 0.001;
	 //double trial_rotation_max = 0.01 * PI_CONST/180.0;
	 //double trial_translation_min = 0.00001;
	 //double trial_rotation_min = 0.0001 * PI_CONST/180.0;
	  //Distance group of atoms are moved during the potential gradient calculation.
	  //NOTE: TRY THIS WITH SMALL VALUES IN THE FIRST PLACE.

	 if (Check_FP_Equality(precision, 0.0))
	    {
	    precision = 0.01;
	     //Energy difference from a steepest descent step that is viewed as 
	     //being "close enough" to an energy minimum.
		 //Default = 1 kJ/mol.
	    }
	 const int MAX_STEPS = 10000;
	  //Maximum number of steps taken before the algorithm gives up.
	 int num_steps_taken = 0;
	 double start_energy, end_energy, change_energy;
	 start_energy = Calculate_Fast_Total_Potential(include_hydro);
	 double translation_step[3];
	 double rotation_step[3];
	 double group_center[3];

	 bool done_minimization = FALSEV;
	 bool translation_done = FALSEV;
	 bool rotation_done = FALSEV;
	 while (!done_minimization)
	       {
		   ++num_steps_taken;

		   //Calculate translation move---


		   Trial_Move_Translate_Get_Gradient(start, end, translation_step, 
			                                 trial_translation);
		   Multiply_XYZ(translation_step, -1.0*step_size_displace);
		   Max_Bound_XYZ(translation_step, max_translation_step);

		   //Execute translation move---
			
		   Shift_Atomic_Coordinates(start, end, translation_step);

		   //Check new energies---

		   end_energy = Calculate_Fast_Total_Potential(include_hydro);
		   change_energy = end_energy - start_energy;
		   start_energy = end_energy;
		    //Get the move along the potential energy surface.
		   if ( (-1.0*change_energy) < precision)
		      {
			  translation_done = TRUEV;
		      }
		   else
		      {
			  translation_done = FALSEV;
		      }

		   //Calculate rotation move---

		   Trial_Move_Rotate_Get_Gradient(start, end, rotation_step, trial_rotation);
		   Multiply_XYZ(rotation_step, -1.0*step_size_rotate);
		   Max_Bound_XYZ(rotation_step, max_rotation_step);

		   //Execute rotation move---

		   Coordinate_Average(start, end, group_center);
		    //Check where the center of the group of atoms is (for rotation).
		   Rotate_AtomsXYZ(start, end, rotation_step, group_center, TRUEV);

		   //Check new energies---

		   end_energy = Calculate_Fast_Total_Potential(include_hydro);
		   change_energy = end_energy - start_energy;
		   start_energy = end_energy;
		   if ( (-1.0*change_energy) < precision)
		      {
			  rotation_done = TRUEV;
		      }
		   else
		      {
			  rotation_done = FALSEV;
		      }

		   //Check if local minimum has been found---

		   if (rotation_done && translation_done)
			//Done all translation and rotation minimizations.
		       {
			   cout << endl << "Minimization completed with: " << num_steps_taken
				    << " steps taken!";
			   Get_Bond_Vec(3, 4, reference_axis);
			   cout << endl << "Atom bond 1 reference direction is now along: "
					<< reference_axis[0] << " " << reference_axis[1] 
					<< " " << reference_axis[2];
			   Get_Bond_Vec(5, 1, reference_axis);
			    //Reverse order should be removed if you are standardizing this.
			   cout << endl << "Atom bond 2 reference direction is now along: "
					<< reference_axis[0] << " " << reference_axis[1] 
					<< " " << reference_axis[2];
			   Get_Bond_Vec(2, 5, reference_axis);
			   cout << endl << "Atom bond 3 reference direction is now along: "
					<< reference_axis[0] << " " << reference_axis[1] 
					<< " " << reference_axis[2];
			   done_minimization = TRUEV;
		       }

		   //Check if we are giving up on finding the minimum---

		   //Add change in precision according to changes in energy, look for zig-zags. 

		   if (num_steps_taken == MAX_STEPS)
		      {
			  cout << endl << "WARNING: RIGID BODY MINIMIZATION DID NOT FIND MINIMUM!";
			  break;
		      }

	       }

     }

//Trial moves---

void molecule_sim::Trial_Move_Translate_Get_Gradient(int start, int end, 
	                                                    double* translate_matrix,
	                                                    double move_size)
 
     {

     }

void molecule_sim::Trial_Move_Rotate_Get_Gradient(int start, int end, 
	                                                 double* rotate_matrix,
	                                                 double move_size)
     {


     }

//Potential calculations---

double molecule_sim::Calculate_Fast_Total_Potential(bool include_hydro) const
     {
	 int num_charged_atoms;
	 int* charged_atom_indices = new int[num_atoms];
	 Get_Charged_Atom_List(charged_atom_indices, num_atoms, num_charged_atoms,
		                   include_hydro);
	  //Only include charged atoms in following summation.

	 double total_potential = 0.0;
	 double pair_potential, distance_temp, charge1, charge2;
	 int index1_temp, index2_temp;
	 for (int a = 0; a < num_charged_atoms; ++a)
	     {
		 index1_temp = charged_atom_indices[a];
		 charge1 = atom_array[index1_temp].Get_Atom_Charge();
		 for (int b = a + 1; b < num_charged_atoms; ++b)
		     {
			 index2_temp = charged_atom_indices[b];
			 charge2 = atom_array[index2_temp].Get_Atom_Charge();
			 distance_temp = atom_array[index1_temp].Get_Distance(atom_array[index2_temp]);
			 pair_potential = charge1*charge2/distance_temp;
			 total_potential += pair_potential;
		     }
	     }
	 return (ELEC_POT_CONSTANT_EV*total_potential);
     }

double molecule_sim::Calculate_Two_Body_Potential(bool include_hydro) const
     {
	 double es_potential = Calculate_ES_Potential(include_hydro);
	 double sr_potential = Calculate_SR_Two_Body_Potential(include_hydro);
	 double total_potential = es_potential + sr_potential;
	 return total_potential;
     }

double molecule_sim::Calculate_SR_Two_Body_Potential(bool include_hydro) const
     {
	 double total_potential = 0.0;
     return total_potential;

	 for (int a = 0; a < num_atoms; ++a)
	     {
		 for (int b = 0; b < num_atoms; ++b)
          //Temporary until neighbor lists are added.
		     {
            // if (Buckingham_Used(a, b))
			    {

			    }
		     }
	     }


     }
*/

//PDB file writting---

	//Final TER line----

		/*out_file << line_feed << "TER";
        index_value_length = String_Size(atom_count + 1);
	    Output_Spaces(out_file, 8 - index_value_length);
		  //Right-justify the atomic index.
	    out_file << atom_count << "     AAA A";
		value_string_length = String_Size(group_tag_index);
		Output_Spaces(out_file, 4 - value_string_length);
	    out_file << group_tag_index; */

/*

void Write_GRO_File(const general_composite& comp, const char* Wfile_name,
					bool include_hydro)
     {
	 ofstream out_file(Wfile_name, ios::out);
	 if (out_file.is_open())
		{
		char system_name[50];
        Copy_With_No_Extension(system_name, Wfile_name);
         //Use the file name without the extension as a system name in the JEMS file.
		out_file << system_name << " , t= 0.0 ";

		double cell_temps[9];
		atom_collection atoms;
		comp.Get_File_Runners(cell_temps, atoms);

		out_file << endl << atoms.Real_Size(include_hydro); 

		double coors[3];
		char atom_name[MAX_ANAME_SIZE];
		for (int a = 0; a < atoms.Size(); ++a)
			{
			if (!atoms[a].Is_Invisible(include_hydro))
			   {
			   atoms[a].Get_Atom_Name(atom_name);
			   out_file << endl << "\t" << comp.Get_Group_Tag(a) << " GROUP" 
				        << atom_name << " " << a; 
			   atoms[a].Get_Atom_Location(coors);
			   out_file << fixed << coors[0] << " " << coors[1] << " " << coors[2];
			   }
			}

		if (comp.Get_Periodicity() > 0)
		   {
		   out_file << endl << cell_temps[0] << " ";
		   }
	    if (comp.Get_Periodicity() > 1)
		   {
		   out_file << cell_temps[4] << " ";
		   }
	    if (comp.Get_Periodicity() > 2)
		   {
		   out_file << cell_temps[8];
		   }
		}
	 out_file.close();	 
     }

	 */

/*

	 void Write_GULP_File(ofstream& out_file, const nanosurface& surface)
     {
     out_file << "conp opti" << endl << endl << "Cartesian region 1";
      //GULP file requires these lines before the addition of the atomic coordinates.
     surface.Store_Atom_Locations(out_file);
      //Write the atomic coordinates.
     out_file << endl << endl << "dump resultsS.gin" ;
      //GULP file requirement.                          
     }

void Write_GULP_File(ofstream& out_file, const nanoparticle& particle)
     {
     out_file << "conp opti" << endl << endl << "Cartesian region 1";
      //GULP file requires these lines before the addition of the atomic coordinates.
     particle.Store_Atom_Locations(out_file);
      //Write the atomic coordinates.
     out_file << endl << endl << "dump resultsP.gin" ;
      //GULP file requirement.                          
     }
     
void Write_GULP_File(const nanosurface& surface, const char* Wfile_name, const char* Dfile_name)
     {
     ofstream out_file(Wfile_name, ios::out);
       //Open up file to write to.
     if (out_file.is_open())
        {
        out_file << "conp opti" << endl << endl << "Cartesian region 1";
         //GULP file requires these lines before the addition of the atomic coordinates.
        surface.Store_Atom_Locations(out_file);
         //Write the atomic coordinates.
        out_file << endl << endl << "dump " << Dfile_name;
         //GULP file requirement.   
        }   
     out_file.close();
       //Done writing.                    
     }

void Write_GULP_File(const nanoparticle& particle, const char* Wfile_name, 
                     const char* Dfile_name)
     {
     ofstream out_file(Wfile_name, ios::out);
       //Open up file to write to.
     if (out_file.is_open())
        {
        out_file << "conp opti" << endl << endl << "Cartesian region 1";
         //GULP file requires these lines before the addition of the atomic coordinates.
        particle.Store_Atom_Locations(out_file);
         //Write the atomic coordinates.
        out_file << endl << endl << "dump " << Dfile_name;
         //GULP file requirement.      
        }
     out_file.close();
       //Done writing.                        
     }

	*/

    /* FOR CLOSEPACKED PLANE CALCULATIONS----
    
     atom basis_copy[INIT_BASIS_SIZE];
     for (int a = 0; a < num_basis_atoms; ++a)
      //Make a copy of the basis.
         {
         basis_copy[a] = basis_atoms[a];     
         }*/
    /* for (int a = 0; a < CPP_LIMIT; ++a)
         {
         for (int b = 0; b < CPP_LIMIT; ++b)
             {
             for (int c = 0; c < CPP_LIMIT; ++c)
                 {
                 packing_test_results[a][b][c] = 0;     
                 }     
             }     
         }*/

/*void crystal_system::Add_Basis(vector<atom>& ze_atom_array, int& index) const
     {
     double temps[3];
     ze_atom_array[index].Get_Atom_Location(temps[0], temps[1], temps[2]);
      //Location of basis origin relative to the system being constructed.
     for (int a = 0; a < basis_atoms.Size(); ++a)
         {
         ze_atom_array[index] = basis_atoms[a];
         ze_atom_array[index].AddTo_Atom_Location(temps[0], temps[1], temps[2]);    
          //Get the basis placed relative to the location passed in the variable
          //arguments (at the initial index).
         ++index;
         }          
     ze_atom_array.Set_Collection_Size(index);        
     }*/
     
/*

         for (int b = 1; b < num_CPP_planes; ++b)
          //Find the closest plane to be cutting the surface with.
             {
             dist_temp = Get_Dist(atom_loc, plane_intersect[b]);
             if (dist_temp < closest_dist)
                {
                closest_normal = b; 
                closest_dist = dist_temp;          
                }
             }
             
             */


//Nanoparticle logic---


/*void nanoparticle::Construct_Particle()
      {
      vector<atom> box_list(INTIAL_ATOM_COUNT);
       Initial array of atoms used in creating the nanoparticle. This array
       is used in the following algorithm to hold the atoms located at lattice 
       point positions (i.e. the following search is purely for lattice points 
       and the basis of the crystal system is considered later)
      int box_index = 0;
       //Index for the atom array.   
           
      const double big_enough_factor = 2.0;
      int max_a_steps = int(big_enough_factor * effective_diameter/Size_A_Step());
      int max_b_steps = int(big_enough_factor * effective_diameter/Size_B_Step());
      int max_c_steps = int(big_enough_factor * effective_diameter/Size_C_Step());
       Maximum steps along the crystal system vectors needed to make a box to 
       cut the nanoparticle from. Based on effective diameter of the particle. 
       The "big enough factor" is included to ensure that there is a large 
       enough box for dealing with any (i.e. non-orthogonal) crystal systems
       
      double box_x = 0.0;
      double box_y = 0.0;
      double box_z = 0.0;
       //Coordinates for moving through the box.
      for (int a = 0; a < max_a_steps; ++a)
          {
          for (int b = 0; b < max_b_steps; ++b)
              {
              for (int c = 0; c < max_c_steps; ++c)
                  {                  
                  const int new_memory_size = 1000;
                  if ( box_list.size() < (box_index + new_memory_size) )
                   //Make sure memory is available for new atoms. 
                     {
                     box_list.resize(box_index + new_memory_size);                  
                     }  
                                    
                  box_list[box_index].Set_Atom_Location(box_x, box_y, box_z);
                  ++box_index;
                  cryst_struct.C_Step(box_x, box_y, box_z);
                   //Incrementation of coordinates along the c vector.
                  }
              cryst_struct.B_Step(box_x, box_y, box_z);
               //Incrementation of coordinates along the b vector.
              cryst_struct.C_Step(box_x, box_y, box_z, -1*max_c_steps);
               //Reset any incrementation from the c vector from the previous vector.   
              } 
          cryst_struct.A_Step(box_x, box_y, box_z);
          cryst_struct.B_Step(box_x, box_y, box_z, -1*max_b_steps);     
          }
        //Box is now defined.  
        
       double center_coordinates[3]; 
       cryst_struct.Multi_Step(center_coordinates, max_a_steps/2, max_b_steps/2, 
                              max_c_steps/2);
        //Look out for the center-of-box basis point.
       double temp_x, temp_y, temp_z;
       double new_x, new_y, new_z;
       for (int a = 0; a < box_index; ++a)
           {
           box_list[a].Get_Atom_Location(temp_x, temp_y, temp_z);
           temp_x -= center_coordinates[0];
           Check_FP_Zero(temp_x);
               //Check floating-point math, specifically for zeros that are not
               //represented as exactly zero (which makes for ugly looking
               //output).    
           temp_y -= center_coordinates[1];
           Check_FP_Zero(temp_y);
           temp_z -= center_coordinates[2];
           Check_FP_Zero(temp_z);
           box_list[a].Set_Atom_Location(new_x, new_y, new_z);
           }
        //Set the coordinate list values relative to the center of the box.
              
      double eff_radius = effective_diameter/2.0; 
      double eff_radius_squared = eff_radius*eff_radius;
      double breathing_room = eff_radius + cryst_struct.Largest_Vector_Mag();     
      double eff_radius_squared_plusBR = breathing_room*breathing_room;
      double temp_dist_square;
      int number_atoms_added;
      
      if (model_type == con_model)
         {
         for (int a = 0; a < box_index; ++a)
             {
             box_list[a].Get_Atom_Location(temp_x, temp_y, temp_z);
             temp_dist_square = Get_CoorSquare_Sum(temp_x, temp_y, temp_z);
             if (temp_dist_square <= eff_radius_squared_plusBR)
              //Make sure basis origin is close to the nanoparticle space.
              //Note: An additional check is done below for the individual
              //atoms of the basis, which is particularly needed for the atoms
              //pushing the limits of the space.
                {
                particle_atoms[num_atoms].Set_Atom_Location(temp_x, temp_y, temp_z);  
                cryst_struct.Add_Basis(particle_atoms, num_atoms, eff_radius);
                } 
             }                      
         }
      else if (model_type == cry_model)
         {
         cout << "NOT READY!";                 
         }    
      }*/

       /*double center_coordinates[3];     
       //Center of the box being created.
       double temp_x = 0.0;
       double temp_y = 0.0;
       double temp_z = 0.0;   
       for (int a = 0; a < box_index; ++a)
           {
           temp_x += box_list[a][0];
           temp_y += box_list[a][1];
           temp_z += box_list[a][2];     
           }
       center_coordinates[0] = temp_x/box_index;
       center_coordinates[1] = temp_y/box_index;
       center_coordinates[2] = temp_z/box_index;
        //Estimate the center of the box by averaging all coordinates.*/


/*

Old size violation logic, pre-covalent radii analysis---

	  bool Size_Violation(const double*, double, int) const;
       /*Takes a position and a size parameter. Returns TRUEV if any atoms
         present in the box are within a distance of one size parameter from
         the position being tested. Only looks at atoms up to the index
         passed to this function (e.g. so that solvent addition need not
         worry about solvent-solvent overlap, which is handled elsewhere).  


bool general_composite::Size_Violation(const double* pos, double min_dist, 
                                       int atom_range) const
     {
     double atom_pos[3];
     double dist_val;
     bool violation = FALSEV;
     for (int a = 0; a < atom_range; ++a)
         {
         atoms[a].Get_Atom_Location(atom_pos);
         dist_val = Get_Dist_OrthoPBC(atom_pos, pos, periodicity, comp_box_size);
         if (dist_val < (min_dist - 0.002))
            {
            violation = TRUEV;
            a = atom_range;          
            }        
         }            
     return violation;                      
     }
	 
*/
     
     /*
     
     void general_composite::Add_SAMParticle(SAM_nanoparticle& part, double* pos, bool is_group)
      {
      int old_composite_size = num_atoms;
      int size_change = part.Copy_SAMParticle(atoms, old_composite_size, pos);  
      num_atoms = old_composite_size + size_change;
      Add_Group(old_composite_size, num_atoms, is_group);                              
      }
      
      */


/*

GROUP TAGS----

int general_composite::Get_Group_Tag(int atom_index) const
     {
	 int tag = -1;
	 if (atom_index < atoms.Size())
	    {
	    tag = group_tag[atom_index];
		}
	 return tag;
	 }*/

/* Molecule cleaning uping---

		    Get_Molecule_Composition_String(mol_fit_string, b_net, a, 
				                            MOL_STRING_LENGTH);
			if (num_cut == 0)
			   {
			   cout << endl << "Molecules have been cut from your system during cleaning"
				    << endl << "Cut list: " << mol_fit_string;
			   }
			else
			   {
			   cout << " | " << mol_fit_string;
			   }
			++num_cut;
*/


//Molecule deletion---

/*void atom_collection::Delete_Molecule(int* bonding_network, int index)
     {
	 int orig_size = num_atoms;
	 int num_removed = 0;
	 for (int a = 0; a < orig_size; ++a)
	     {
		 if (bonding_network[a] == index)
		    {
			Delete_Atom(a - num_removed);
			++num_removed;
		    }
	     }
     }*/

/* Approximate Debye XRD calculation---

double Get_Approx_Debye_Intensity(double, double) const;
	  /*Returns the x-ray diffraction intensity for a given scattering angle and
	    X-ray wavelength. Calculated by the Debye formula, assuming the atomic scattering 
		factors are proportional to Z and (1 - sin(scattering angle)).
	    This is VERY APPROXIMATE and is not meant for obtaining professional results.

double atom_collection::Get_Approx_Debye_Intensity(double angle, double wavelength) const
     {
	 double intensity = 0.0;
	 double dist, as_factor, scat_factor;
	 for (int a = 0; a < num_atoms; ++a)
	     {
		 for (int b = a; b < num_atoms; ++b)
		     {
			 dist = atom_array[a].Get_Distance(atom_array[b]);
			 as_factor = double(atom_array[a].Get_Atom_Number()*
				                atom_array[b].Get_Atom_Number());
			  //Approximate scattering factor vs. atomic number relationship.
			 scat_factor = 4*PI_CONST*dist*sin(angle)/wavelength;
			 if (a != b)
			    {
			    intensity += as_factor * sin(scat_factor)/(scat_factor);
			    }
			 else
			  //Self-scattering logic. a-a interaction has half the weight as
			  //that of a-b (also b-a) interaction.
			    {
			    intensity += as_factor/2.0;
			    }
		     }
	     }

	 intensity *= ( 1 - sin(angle) ) * ( 1 - sin(angle) );
	  //Scale intensities by the decay in scattering intensity with 
	  //increasing scattering angle. This is very approximate.
	 return intensity;
     }
*/

//Utility string get numbers stuff---

/*

    char TERM_CHAR = '#';
    file_reader.get(file_text, array_size, TERM_CHAR);
     //Dump the script file into a character array.
          
    file_size = 0;
    hile (file_text[file_size] != TERM_CHAR)
      //Get size of the character string taken from the whole text file.
      //Also convert line feed characters to null characters
      //for later string logic.
         {
         if (file_text[file_size] == '\n')
             {
             file_text[file_size] = '\0';                      
             }
         ++file_size;
         }  
         
         
int Get_FP_Numbers(char* ze_string, double* numbers, int variable_max)
       {
       char temp_string[100];
       char* temp_pointer;
       strcpy(temp_string, ze_string);
       int string_length = strlen(ze_string);
       int count_temp = 0;
       int first_digit_index;
       bool found_first_digit = FALSEV;
       bool is_num;
       
       
       for (int index = 0; index < string_length; ++index)
           {
           if (temp_string[index] == '.')
            //Special treatment for decimal points.
              {
              if (!found_first_digit)
               //Dealing with numbers shown like ".409"
                 {
                 found_first_digit = TRUEV;
                 first_digit_index = index;
                 if (index > 0)
                  //Look for minus signs to be included.
                    {
                    if (temp_string[index - 1] == '-')
                       {
                       first_digit_index = index - 1;                   
                       }
                    }
                  }
              continue;
              }    
                
           is_num = Is_Number(temp_string[index]);
           if (!found_first_digit && is_num)
              {
              found_first_digit = TRUEV;    
              first_digit_index = index;  
              if (index > 0)
               //Look for minus signs to be included.
                {
                if (temp_string[index - 1] == '-')
                    {
                    first_digit_index = index - 1;                   
                    }
                }                      
              }     
           else if (!is_num)
            //If true, the loop is going through one of the numbers and just found the end of the number
              {
              temp_string[index] = '\0'; 
              temp_pointer = &(temp_string[first_digit_index]);
              numbers[count_temp] = Get_FP_Number(temp_pointer);
              found_first_digit = FALSEV;
              ++count_temp;     
              if (count_temp == variable_max)
               //If true, hit the maximum allowed number of variable reads.
                 {
                 index = string_length;            
                 }
              }
           }                
       return count_temp;
       }
         
               
*/


//RDF calculation---

			       /*double tempA[3];
                   atom_array[a].Get_Atom_Location(tempA);
                   double tempB[3];
                   atom_array[b].Get_Atom_Location(tempB);
			 if (temp_distance < 0.0543)
				{
				cout << endl << a << " " << b << " " << atom_name1 << " " << atom_name2 << endl 
					 << tempA[0] << " " << tempA[1] << " " << tempA[2] << " " << endl << 
					   tempB[0] << " " << tempB[1] << " " << tempB[2];
				}  */

             /*if ( (atom_name1[1] == 'u') && (atom_name2[0] == 'S') )
                {
                if (temp_distance < 0.3)
                   {
                   double tempA[3];
                   atom_array[a].Get_Atom_Location(tempA);
                   double tempB[3];
                   atom_array[b].Get_Atom_Location(tempB);
                   
                   if (b == 2187)
                      {
                      cout << endl << tempA[0] - 3.25 << " " << tempA[1]- 3.25 << " " << tempA[2]- 3.25;
                      cout << endl << tempB[0]- 3.25 << " " << tempB[1]- 3.25 << " " << tempB[2]- 3.25;   
                      cout << endl << atom_array[a].atom_name << " " << atom_array[b].atom_name << " " << temp_distance;         
                      }
                   }                 
                } */

//PBC coding---

     //distance_coors[0] = Abs_Difference(coor1[0], coor2[0]);
     //distance_coors[1] = Abs_Difference(coor1[1], coor2[1]);
     //distance_coors[2] = Abs_Difference(coor1[2], coor2[2]);

       /*  if (coor1[a] < coor2[a])
            {
            temp_difference = Abs_Difference(coor1[a], coor2[a] - box_dimensions[a]);
            Set_Max(distance_coors[a], temp_difference);
            }
         else
            {
            temp_difference = Abs_Difference(coor1[a] - box_dimensions[a], coor2[a]);
            Set_Max(distance_coors[a], temp_difference);       
            }  */

//PDB READING---

/*	   if (read_length > 3)
			//Due to PDB file formatting, its possible for the components of the atom
			//description to blend together with no spaces between them. A real pain
			//in the butt for typical C++ reading logic.
		      {
			  temp_atom_name[3] = '\0';
		      }*/

//

/*const double DEFAULT_BOND_COMPRESSION = 0.95;
 //Indicates the range that is considered an acceptable bond length.
 //Useful for spatial overlap analysis in system construction.
 //E.g. For bond length = 2 A, allowed compression = 0.95, any distance greater than
 //1.90 A will be considered OK in a system.
const double DEFAULT_BOND_ELONGATION = 1.0/DEFAULT_BOND_COMPRESSION;
 //Similiar meaning to above.
const double DEFAULT_VARIATION_PARAMETER = 1.0 - DEFAULT_BOND_COMPRESSION;*/

//Atom rotation---


/*void atom::Rotate_Atom(double* start_vec, double* final_vec)
       {
       double rotation_vec[3];
       double angle;   
       Calculate_Rotation_Parameters(start_vec, final_vec, rotation_vec, angle);
       double rot_matrix[9];
       Calc_Rotation_Matrix(rotation_vec, angle, rot_matrix);
       Rotate_Vector(rot_matrix, rel_atom_location);                       
       }*/


//Findings----

//2.0 nm nanoparticle: Amorphous has ~1250 atoms (variation on the order of +/- 10 atoms) and spherical particle
//has 1985 atoms...63% density

//Nanosurface for-loop surface slab creation---

      //Set up search parameters--- 

     // double largest_dim = Get_Largest(dimensions[0], dimensions[1], dimensions[2]);
     // double smallest_vec = cryst_struct.Smallest_NonZero_Step();
      
     // double search_parameter = largest_dim/smallest_vec;
       /*Number of smallest lattice vector steps required to move across the 
       largest slab dimension. This is a first guess at how many steps will be 
       needed to explore coordinate space. This is very crude as the unit cell
       vector orientation may look nothing like the surface slab vectors in. 
       Thus, the following factor increases the amount of space explored as to 
       ensure that the full nanosurface slab is created */
     /* const double DEEP_SPACE = 4.0;
      int fsp = int(search_parameter * DEEP_SPACE);
       //Final Search Parameter. This is used in the space exploration algorithm
       //below.
      int range = fsp/2;
      Set_Min(range, 3);

      //Carry out the search and define the surface--- 
      double slab_coors[3];
       //Coordinates for finding points of interest in the nanosurface.  
      
      int current_size = 0; 
      for (int a = 1 - range; a < range; ++a)
          {
          for (int b = 1 - range; b < range; ++b)
              {
              for (int c = 1 - range; c < range; ++c)
                  { 
                  slab_coors[0] = dimensions[0]/2.0;
                  slab_coors[1] = dimensions[1]/2.0;
                  slab_coors[2] = dimensions[2]/2.0;
                   //Search around the slab center and assume a lattice point is
                   //at that center.
                  cryst_struct.Multi_Step(slab_coors, a, b, c);
                  current_size = surface_atoms.Size();
                  surface_atoms[current_size].Set_Atom_Location(slab_coors);
                  cryst_struct.Add_Basis_ToSurface(surface_atoms, current_size,
                                                   dimensions, stuff_space); 
                   //Adds the basis atoms at the coordinate point being
                   //explored, not adding outside the requested dimensions of 
                   //the slab.                 
                  }
              } 
          } 
          */

//Nanoparticle version---

      /*int max_a_steps = 1;
      int max_b_steps = 1;
      int max_c_steps = 1;
      
      const double EXTRA_ROOM = 2.01;
      if (!Check_FP_Equality(cryst_struct.Size_A_Step(), 0.0))
         {
         max_a_steps = int(EXTRA_ROOM * effective_diameter/cryst_struct.Size_A_Step());
         }
      if (!Check_FP_Equality(cryst_struct.Size_B_Step(), 0.0))
         {
         max_b_steps = int(EXTRA_ROOM * effective_diameter/cryst_struct.Size_B_Step());
         }
      if (!Check_FP_Equality(cryst_struct.Size_C_Step(), 0.0))
         {
         max_c_steps = int(EXTRA_ROOM * effective_diameter/cryst_struct.Size_C_Step());
         }
       /*Maximum steps along the crystal lattice vectors needed to explore a large
       enough coordinate space to place the nanoparticle within. 
       The "EXTRA ROOM" is included to ensure that there is a large 
       enough space for dealing with any (i.e. non-orthogonal) crystal systems
       */
       
      /*double part_coor[3];
       //Coordinates for moving through coordinate space.
      double center_coordinates[3] = {0.0, 0.0, 0.0}; 
      cryst_struct.Multi_Step(center_coordinates, max_a_steps/2, max_b_steps/2, 
                              max_c_steps/2);
       //Look out for the center-of-box lattice point.
       //Choose an atom in the middle of the lattice space being explored to
       //serve as the effective center.
      Neg_XYZ(center_coordinates);
      int current_size = 0;
      for (int a = 0; a < max_a_steps; ++a)
          {
          for (int b = 0; b < max_b_steps; ++b)
              {
              for (int c = 0; c < max_c_steps; ++c)
                  {
                  Set_XYZ(part_coor, center_coordinates);
                  cryst_struct.Multi_Step(part_coor, a, b, c);
                   //Get the lattice point position being explored.
                  current_size = particle_atoms.Size();
                  particle_atoms[current_size].Set_Atom_Location(part_coor);
                  if (model_type == SPH_MODEL)
                     { 
                     cryst_struct.Add_Basis_ToSParticle(particle_atoms, current_size, 
                                                        effective_radius, stuff_space);
                      //Add a basis of atoms at that lattice point, not allowing
                      //any atoms outside the spherical volume of the nanoparticle
                      //to be added.
                     }
                  else
                     {
                     cryst_struct.Add_Basis_ToFParticle(particle_atoms, current_size,
                                                        effective_radius, stuff_space);
                      //Add a basis of atoms at that lattice point, not allowing
                      //any atoms outside the internal volume of the faceted nanoparticle.     
                     }
                  } 
              }    
          }    */

//Recursion for orthorhombic sampling---


/*void Ortho_Recursion(const double* A_vec, const double* B_vec, const double* C_vec,  
	                 const double* sizes, double** point_array, int max_points, int& point_index,
					 const double* test_pos)
	//Recursive function for the following function to use.
     {
	 static int num_times = 0;
	 ++num_times;
	 cout << endl << num_times;
     bool memory_okay = (point_index < max_points);
	 double abs_coors[3];
	 Set_XYZ(abs_coors, test_pos);
	 Abs_XYZ(abs_coors);
	 bool in_box = Check_BoundariesMax(abs_coors, sizes);
     if (memory_okay && in_box)
	    {
        Set_XYZ(point_array[point_index], test_pos);
		++point_index;

		double try_pos[3];
		SetAdd_XYZ(try_pos, test_pos, A_vec);
	    Ortho_Recursion(A_vec, B_vec, C_vec, sizes, point_array, max_points, point_index, try_pos);
		SetSub_XYZ(try_pos, test_pos, A_vec);
	    Ortho_Recursion(A_vec, B_vec, C_vec, sizes, point_array, max_points, point_index, try_pos);
		 //A.

		SetAdd_XYZ(try_pos, test_pos, B_vec);
	    Ortho_Recursion(A_vec, B_vec, C_vec, sizes, point_array, max_points, point_index, try_pos);
		SetSub_XYZ(try_pos, test_pos, B_vec);
	    Ortho_Recursion(A_vec, B_vec, C_vec, sizes, point_array, max_points, point_index, try_pos);
		 //B.

		SetAdd_XYZ(try_pos, test_pos, C_vec);
	    Ortho_Recursion(A_vec, B_vec, C_vec, sizes, point_array, max_points, point_index, try_pos);
		SetSub_XYZ(try_pos, test_pos, C_vec);
	    Ortho_Recursion(A_vec, B_vec, C_vec, sizes, point_array, max_points, point_index, try_pos);
		 //C.

	    }
     }
	 */


/////

/*

void Get_Random_Mesh_Sampling(double mesh_pam, double vari_pam, long int max_tries,
	                          double** mesh_grid, long int& num_points, double box_size)
     {
     double start_pos[3] = { 0.0, 0.0, 0.0 };  
	 mesh_grid[0][0] = start_pos[0];
	 mesh_grid[0][1] = start_pos[1];
	 mesh_grid[0][2] = start_pos[2];
     num_points = 1;
     //Expand_Mesh_Grid(mesh_pam, vari_pam, max_tries, mesh_grid, box_size, 0, num_points);
	 
	 
	 Order_Coordinate_List(mesh_grid, num_points);
     } */


//Charge neutral nanosurface---


	/*  for (int a = 0; a < point_index; ++a)
         //Move the origin in the point sampling by the (box center - box corner) vector. 
	     //This will be better for the later basis addition downs.
	      {
          points[a][0] += dimensions[0]/2.0;
		  points[a][1] += dimensions[1]/2.0;
		  points[a][2] += dimensions[2]/2.0;
	      }

	  bool size_variation_criterion_met = FALSEV;
	  int vari_index = 0;
	  double prev_charge, total_charge;
	  bool direction_change = FALSEV;
	  int vari_count = 0;
	  int cycle_count = 0;
	  const int MAX_CYCLES = 5;
	  const int MAX_VARI_INDEX = 3;
	   //Only cut from the x and y directions.
	  int current_size;
	  while (!size_variation_criterion_met)
		{
	    if (vari_count != 0)
			{
			if (!direction_change)
				{
				prev_charge = total_charge;
				}
			else
				{
				direction_change = FALSEV;
				}
			}
		else
			{
			prev_charge = HUGE_VALUE;
			 //Give dummy value for later logic.
			}

	    current_size = 0;
	    for (int a = 0; a < point_index; ++a)
	      { 
	      surface_atoms[current_size].Set_Atom_Location(points[a]);
	      cryst_struct.Add_Basis_ToSurface(surface_atoms, current_size,
                                          dimensions, stuff_space);
	      }
	    total_charge = surface_atoms.Return_Total_Charge();
		if (!variable_size || Check_FP_Equality(total_charge, 0.0))
		  {
		  size_variation_criterion_met = TRUEV;
		  }
		else
		  {
		  cout << endl << "Charge reduction is: " << abs(total_charge) - abs(prev_charge)
		       << " for: " << vari_count << " direction " << vari_index
			   << endl << " with prev of " << prev_charge << " and new of " << total_charge;
		  if (Check_FP_LessThan(abs(prev_charge), abs(total_charge)))
			{
			dimensions[vari_index] += SIZE_VARI_CONST;
			++vari_index;
			direction_change = TRUEV;
			if (vari_index == MAX_VARI_INDEX)
		     //Done one cycle once we finish the z coordinate.
			   {
			   ++cycle_count;
			   if (cycle_count == MAX_CYCLES)
				  {
				  break;
				  }
			   vari_index = 0;
			   }
			}
		  dimensions[vari_index] -= SIZE_VARI_CONST;
		  surface_atoms.Delete_Atoms(0, current_size);
		  ++vari_count;
		  }
	    }*/

/*

	int num_surf_atoms;
	int* surf_atom_indices = new int[num_atoms];
	 //Can't have more surface atoms than the actual atom count.
	Get_Surface_Atom_List(surf_atom_indices, num_atoms, num_surf_atoms, 
		                  dist_max, AB_only, periodicity, box_sizes, internal_CN);

	double low_pos_sums[3] = { 0.0, 0.0, 0.0 };
	double high_pos_sums[3] = { 0.0, 0.0, 0.0 };
	double high_count[3] = { 0.0, 0.0, 0.0};
	double low_count[3] = { 0.0, 0.0, 0.0};
	double box_center[3];
	Coordinate_Average(box_center);
	double location_temp[3];
	int largest_component;
	int coor_number, under_coordination;
	double SA_contribution_factor;
	for (int a = 0; a < num_surf_atoms; ++a)
	    {
		if (atomic_number_specific != -1)
		   {
		   if (atom_array[surf_atom_indices[a]].Get_Atom_Number() 
			   != atomic_number_specific)
		      {
			  continue;
		      }
		   }
		atom_array[surf_atom_indices[a]].Get_Atom_Location(location_temp);
		Sub_XYZ(location_temp, box_center);
		largest_component = Get_Largest_RelVec_Component(location_temp, box_center);
		coor_number = Get_Coordination_Number(surf_atom_indices[a], dist_max, AB_only, 
			          periodicity, box_sizes);
		under_coordination = internal_CN - coor_number;
		SA_contribution_factor = sqrt(double(under_coordination));
		if (location_temp[largest_component] < 0.0)
			   {
			   low_pos_sums[largest_component] += location_temp[largest_component]
			                                      *SA_contribution_factor;
			   low_count[largest_component] += SA_contribution_factor;
			   }
		else
			   {
			   high_pos_sums[largest_component] += location_temp[largest_component]
			                                       *SA_contribution_factor;
			   high_count[largest_component] += SA_contribution_factor;
			   }
		}

	for (int a = 0; a < 3; ++a)
	    {
		if ( (high_count[a] > 0.0) && (low_count[a] > 0.0) )
			{
			box_fit[a] = high_pos_sums[a]/double(high_count[a]) 
			           - low_pos_sums[a]/double(low_count[a]);
			}
		else
			{
			box_fit[a] = 0.0;
			}
	    }

	delete[] surf_atom_indices; 

	*/


     
     /* 
     
     CRYSTAL SYSTEM INITIALIZATIONS
     
         //  my_cryst.Set_Orthorhombic_Vectors(0.2, 0.1, 0.2);
    //  my_cryst.Add_Basis_Atom(test_atom, 0.0, 0.0 0.0);
    //  my_cryst.Add_Basis_Atom(test_atomB, 0.1, 0.0, 0.1);
     // my_cryst.Set_BCC_System(test_atom, 0.324);
     
     CLOSE-PACKED PLANE TESTING
     
     for (int a = 0; a < 20; ++a)
          {
          int ind1, ind2, ind3;
          my_cryst.Get_ClosePacked_Plane(a, ind1, ind2, ind3);
          cout << a << " " << ind1 << " " << ind2 << " " << ind3 << endl;
          }
          system("PAUSE");*/
      
     // nanoparticle piece_of_nano(2.0, con_model, my_cryst); 
     // piece_of_nano.Print_Atom_Locations();




/*

      //INTRODUCTION---      
            
      cout << "WELCOME TO THE NEW PROGRAM!";
          
      //DECLARE TEST ATOMS---    
          
      char atom_nameAu[30] = "Au";
      char atom_nameAg[30] = "Ag";
      atom test_atomAu(atom_nameAu, 79, 197.001);
      atom test_atomAg(atom_nameAg, 47, 108.001);
      
     //DECLARE CRYSTAL SYSTEMS--- 
      
      crystal_system my_crystAu;
      my_crystAu.Set_FCC_System(test_atomAu, 0.408);   
      
      crystal_system my_crystAg;
      my_crystAg.Set_FCC_System(test_atomAg, 0.409);  

      //CREATE NANOSURFACE SYSTEM---
       
      nanosurface surface(my_crystAg, 4.0, 4.0, 2.0, 1, 1, 1);
      
      //CREATE NANOPARTICLES---
      
      nanoparticle part1(2.0, con_model, my_crystAu); 
      nanoparticle part2(1.0, con_model, my_crystAu); 
      nanoparticle part3(1.5, con_model, my_crystAu); 
      nanoparticle part4(0.5, con_model, my_crystAu); 
      nanoparticle part5(1.0, con_model, my_crystAu); 
      nanoparticle part6(2.0, con_model, my_crystAu); 
      nanoparticle part7(1.0, con_model, my_crystAu); 
      
      //CREATE WATER MOLECULE---
      
      molecule water;
      water.Add_Water(0.1, 0.1, 0.1);
       //H2O is the nectar of life, baby.
       
      //PUT IT ALL TOGETHER---
      
      double box_dimensions[3] = {6.0, 6.0, 6.0};
      general_composite comp(box_dimensions);
      
      double surface_coordinates[3] = {1.0, 1.0, 4.0};
      comp.Add_NanoSurface(surface, surface_coordinates, TRUEV);
   
      double nanop_coordinates1[3] = {2.0, 2.0, 2.0};
      comp.Add_NanoParticle(part1, nanop_coordinates1, TRUEV); 
      double nanop_coordinates2[3] = {5.0, 5.0, 1.23};
      comp.Add_NanoParticle(part2, nanop_coordinates2, TRUEV);    
      double nanop_coordinates3[3] = {4.0, 4.0, 1.5};
      comp.Add_NanoParticle(part3, nanop_coordinates3, TRUEV); 
      double nanop_coordinates4[3] = {6.0, 6.0, 2.0};
      comp.Add_NanoParticle(part4, nanop_coordinates4, TRUEV); 
      double nanop_coordinates5[3] = {0.0, 5.0, 1.5};
      comp.Add_NanoParticle(part5, nanop_coordinates5, TRUEV); 
      double nanop_coordinates6[3] = {0.0, 0.0, 1.0};
      comp.Add_NanoParticle(part6, nanop_coordinates6, TRUEV); 
      double nanop_coordinates7[3] = {4.0, 3.0, 2.2};
      comp.Add_NanoParticle(part7, nanop_coordinates7, TRUEV); 
      
      comp.FillSpace_With_Solvent(water, 0.400, TRUEV);
      
      //OUTPUT ATOMIC COORDINATES---
      
      char file_name[15] = "Coors.txt\0";
      char dump_name[15] = "Dumper.txt\0";
      Write_GULP_File(comp, file_name, dump_name);
      
      
      char QSTEM_Name[20] = "firsttest.cfg\0";
      comp.Get_QSTEM_File(QSTEM_Name);
      */


					  /*cout << endl << a << " " << b << endl << network[a] << " " << network[b];
					  char temp_name[3]; char tempp_name[3]; 
					  atom_array[a].Get_Atom_Name(temp_name);
					  atom_array[b].Get_Atom_Name(tempp_name);
					  cout << endl << temp_name << " " << tempp_name; system("PAUSE");
					  cout << endl << "WARNING: SOMETHING FUNKY WITH BOND NETWORK CALC!";*/


/*void atom_collection::Get_Debye_Intensity(double* xrd_data, double min_angle, 
	                                      double step_size, int num_steps, 
										  double wavelength, double** scat_fit,
										  double dw_factor, bool surf_effect,
										  double surface_dw_enhance, 
										  double coor_cutoff, int coor_num) const
     {
	 double dist, as_factor, scat_factor, angle, pre_factor;
	 double s_value, neg_s_square, sin_angle, scat_1, scat_2, dw_atten, surf_atten;
	 int an_1, an_2;

	 for (int a = 0; a < num_steps; ++a)
	  //Initialize to zero intensity.
	     {
	     xrd_data[a] = 0.0;
	     }

	 bool* surface_atom = new bool[num_atoms];
	 int surface_atom_count = 0;
	 double surface_dw_factor = surface_dw_enhance*dw_factor;
	 if (surf_effect)
	   {
	   for (int a = 0; a < num_atoms; ++a)
	     {
		 if (Is_Surface_Atom(a, coor_cutoff, TRUEV, coor_num))
		    {
			surface_atom[a] = TRUEV;
			++surface_atom_count;
		    }
		 else
		    {
		    surface_atom[a] = FALSEV;
		    }
	     }
	   cout << endl << "Surface atom count: " << surface_atom_count
		    << " out of " << num_atoms;
	   }

	   for (int a = 0; a < num_atoms; ++a)
	     {
		 for (int b = a; b < num_atoms; ++b)
		  //For each pair...
		     {
			 dist = atom_array[a].Get_Distance(atom_array[b]);
			  //Approximate scattering factor vs. atomic number relationship.
			 angle = min_angle;
			 pre_factor = 4*PI_CONST*dist/wavelength;

			 an_1 = atom_array[a].Get_Atom_Number();
			 an_2 = atom_array[b].Get_Atom_Number();

			 if (Check_FP_Equality(scat_fit[an_1][0], 0.0)
				 || Check_FP_Equality(scat_fit[an_2][0], 0.0) ) 
			     {
				 cout << endl << "WARNING: SCATTERING FACTOR IS UNDEFINED OR ZERO!";
			     }

			 for (int c = 0; c < num_steps; ++c)
			  //Calculate the contribution of the given pair to intensities at 
			  //all angles analyzed.
			     {
				 angle += step_size;
				 sin_angle = sin(angle);
				 s_value = sin_angle/wavelength;
				 neg_s_square = s_value*s_value*-1.0;
				 scat_1 = scat_fit[an_1][0]*exp(scat_fit[an_1][1]*neg_s_square)
					    + scat_fit[an_1][2]*exp(scat_fit[an_1][3]*neg_s_square)
						+ scat_fit[an_1][4]*exp(scat_fit[an_1][5]*neg_s_square)
						+ scat_fit[an_1][6]*exp(scat_fit[an_1][7]*neg_s_square);
				 scat_2 = scat_1;
				 if (an_2 != an_1)
				    {
				    scat_2 = scat_fit[an_2][0]*exp(scat_fit[an_2][1]*neg_s_square)
					       + scat_fit[an_2][2]*exp(scat_fit[an_2][3]*neg_s_square)
					       + scat_fit[an_2][4]*exp(scat_fit[an_2][5]*neg_s_square)
						   + scat_fit[an_2][6]*exp(scat_fit[an_2][7]*neg_s_square);
				    }
				 dw_atten = exp(2.0*dw_factor*neg_s_square);
				 if (surf_effect)
				  //Surface enhancement of DW factors.
				    {
					surf_atten = exp(surface_dw_factor*neg_s_square);
					if (surface_atom[a])
					   {
					   dw_atten *= surf_atten;
					   }
					if (surface_atom[b])
					   {
					   dw_atten *= surf_atten;
					   }
				    }
				 as_factor = scat_1 * scat_2 * dw_atten;
				 if (a != b)
				    {
					scat_factor = pre_factor*sin_angle;
			        xrd_data[c] += as_factor * sin(scat_factor)/scat_factor;
				    }
				 else
				  //Self-scattering, which would not work with computer math
				  //in the above equation (limit of sin(x)/x -> 0 is 1).
				  //Also, the a-a interaction is weighted half as much
				  //as other pair interactions (counted once, not twice).
				    {
				    xrd_data[c] += as_factor/2.0;
				    }
				 }
		     }
	     }
	 delete[] surface_atom;
     }
    


	 void atom_collection::Get_Debye_Intensity(double* xrd_data, double min_angle, 
	                                      double step_size, int num_steps, 
										  double wavelength, double** scat_fit,
										  double dw_factor, bool surf_effect,
										  double surface_dw_enhance, 
										  double coor_cutoff, int coor_num) const
     {
	 double dist, as_factor, scat_factor, angle;
	 double s_value, neg_s_square, sin_angle, scat_1, scat_2, dw_atten, surf_atten;
	 int an_1, an_2;

	 for (int a = 0; a < num_steps; ++a)
	  //Initialize to zero intensity.
	     {
	     xrd_data[a] = 0.0;
	     }

	 bool* surface_atom = new bool[num_atoms];
	 int surface_atom_count = 0;
	 double surface_dw_factor = surface_dw_enhance*dw_factor;
	 if (surf_effect)
	   {
	   for (int a = 0; a < num_atoms; ++a)
	     {
		 if (Is_Surface_Atom(a, coor_cutoff, TRUEV, coor_num))
		    {
			surface_atom[a] = TRUEV;
			++surface_atom_count;
		    }
		 else
		    {
		    surface_atom[a] = FALSEV;
		    }
	     }
	   cout << endl << "Surface atom count: " << surface_atom_count
		    << " out of " << num_atoms;
	   }

	 double* scat_factors;
	 scat_factors = new double[MAX_ATOMIC_NUMBER];
	 
	 angle = min_angle;
	 double pre_factor = 4.0*PI_CONST/wavelength;
	 for (int c = 0; c < num_steps; ++c)
	  //Do sampling, angle by angle.
	   {
	   angle += step_size;
	   sin_angle = sin(angle); 
	   s_value = sin_angle/wavelength;
	   neg_s_square = s_value*s_value*-1.0;
	   for (int d = 0; d < MAX_ATOMIC_NUMBER; ++d)
	    //Determine the atomic scattering factors at this angle.
	       {
	       scat_factors[d] = scat_fit[d][0]*exp(scat_fit[d][1]*neg_s_square)
					       + scat_fit[d][2]*exp(scat_fit[d][3]*neg_s_square)
						   + scat_fit[d][4]*exp(scat_fit[d][5]*neg_s_square)
						   + scat_fit[d][6]*exp(scat_fit[d][7]*neg_s_square);
	       }
	   for (int a = 0; a < num_atoms; ++a)
	     {
		 for (int b = a; b < num_atoms; ++b)
		  //For each pair...
		     {
			 dist = atom_array[a].Get_Distance(atom_array[b]);
			  //Approximate scattering factor vs. atomic number relationship.
			 scat_factor = pre_factor*sin_angle*dist;

			 an_1 = atom_array[a].Get_Atom_Number();
			 an_2 = atom_array[b].Get_Atom_Number();

			 if (Check_FP_Equality(scat_fit[an_1][0], 0.0)
				 || Check_FP_Equality(scat_fit[an_2][0], 0.0) ) 
			     {
				 cout << endl << "WARNING: SCATTERING FACTOR IS UNDEFINED OR ZERO!";
			     }

		     scat_1 = scat_factors[an_1];
		     scat_2 = scat_factors[an_2];
			 dw_atten = exp(2.0*dw_factor*neg_s_square);
			 if (surf_effect)
				  //Surface enhancement of DW factors.
				 {
			     surf_atten = exp(surface_dw_factor*neg_s_square);
				 if (surface_atom[a])
				    {
				    dw_atten *= surf_atten;
				    }
				 if (surface_atom[b])
				    {
				    dw_atten *= surf_atten;
				    }
				 }
			  as_factor = scat_1 * scat_2 * dw_atten;
			  if (a != b)
				    {
			        xrd_data[c] += as_factor * sin(scat_factor)/scat_factor;
				    }
			  else
				  //Self-scattering, which would not work with computer math
				  //in the above equation (limit of sin(x)/x -> 0 is 1).
				  //Also, the a-a interaction is weighted half as much
				  //as other pair interactions (counted once, not twice).
				    {
				    xrd_data[c] += as_factor/2.0;
				    }
		     }
	     }
	   }
	 delete[] scat_factors;
	 delete[] surface_atom;
     }*/