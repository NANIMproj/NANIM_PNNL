
#include "stdafx.h"

#include "DWN_cryst.h"

//Constructors and initialization functions---

crystal_system::crystal_system()
      {
      Initialize();
      }

void crystal_system::Initialize()
	  {
      basis_atoms.Set_Initial_Allocation(INIT_BASIS_SIZE);                         
      
      for (int a = 0; a < 3; ++a)
       {
	   Zero_XYZ(lattice_vectors[a]);
       }
      num_faceting_planes = 0;
	  }

//Destructors---

crystal_system::~crystal_system()
      {                     
      }

     
//Lattice vector settings--- 
        
void crystal_system::Set_Cubic_Vectors(double lattice_pam)
      {
      for (int a = 0; a < 3; ++a)
         {
	     Zero_XYZ(lattice_vectors[a]);
	     lattice_vectors[a][a] = lattice_pam;
		  //Diagonal components are set to the lattice parameter.
         }        
      }

void crystal_system::Set_Tetragonal_Vectors(double lattice_pamA, double lattice_pamC)
      {
      for (int a = 0; a < 3; ++a)
		 {
		 Zero_XYZ(lattice_vectors[a]);
		 }
	  lattice_vectors[0][0] = lattice_pamA;
	  lattice_vectors[1][1] = lattice_pamA;
	  lattice_vectors[2][2] = lattice_pamC;
	   //Diagonal components = A A C.
      }
     
void crystal_system::Set_Orthorhombic_Vectors(double lattice_pamA, double lattice_pamB, 
                                              double lattice_pamC)
      {
      for (int a = 0; a < 3; ++a)
		 {
		 Zero_XYZ(lattice_vectors[a]);
		 }
	  lattice_vectors[0][0] = lattice_pamA;
	  lattice_vectors[1][1] = lattice_pamB;
	  lattice_vectors[2][2] = lattice_pamC;
	   //Diagonal components = A B C.
      }

void crystal_system::Set_Hexagonal_Vectors(double lattice_pamA, double lattice_pamC)
     {
     const double GAMMA_ANGLE = (2.0/3.0) * PI_CONST;
      //Vector angles = 90/90/120 degrees.          
     const double cos_angle = cos(GAMMA_ANGLE);
     const double sin_angle = sin(GAMMA_ANGLE);
      //Trigonmetry must be used for the second lattice vector's determination.
	 for (int a = 0; a < 3; ++a)
		 {
		 Zero_XYZ(lattice_vectors[a]);
		 }
     lattice_vectors[0][0] = lattice_pamA;
      //Define a-vector 1 along the x direction.
     lattice_vectors[1][0] = lattice_pamA * cos_angle;
     lattice_vectors[1][1] = lattice_pamA * sin_angle;
      //Define a-vector 2 along the x/y directions.       
     lattice_vectors[2][2] = lattice_pamC;
      //Define c-vector along the z direction.     
     }

void crystal_system::Set_Vectors(const double* parameters)
      {
      for (int a = 0; a < 9; ++a)
          {
          lattice_vectors[a / 3][a % 3] = parameters[a];   
		   //Relate [9] matrix to [3][3] matrix.
          }                   
      }  

void crystal_system::Set_Vectors(const double* vec1, const double* vec2, 
                                 const double* vec3)
     {
     for (int a = 0; a < 3; ++a)
          {
          lattice_vectors[0][a] = vec1[a];   
          lattice_vectors[1][a] = vec2[a]; 
          lattice_vectors[2][a] = vec3[a]; 
          }                        
     }

void crystal_system::Set_Vectors(const double* vecm, const double* veca)
	//Uses following conventions:
	//a vector is along the x-axis
	//b vector is in the xy plane
	//c vector makes a right-handed set with a and b.
      {
	  lattice_vectors[0][0] = vecm[0];
	  lattice_vectors[0][1] = 0.0;
	  lattice_vectors[0][2] = 0.0;
	   //a vector is purely along the x-axis.
	  lattice_vectors[1][0] = vecm[1]*cos(veca[2]);
	  lattice_vectors[1][1] = vecm[1]*sin(veca[2]);
	  lattice_vectors[1][2] = 0.0;
	   //b vector is determined by angle with the "a" vector
	   //such that it lies in the xy plane.

	  //c vector, the complicated one---

	  lattice_vectors[2][0] = vecm[2]*cos(veca[1]);
	  lattice_vectors[2][1] = vecm[2]*(cos(veca[0])-
		                      cos(veca[2])*cos(veca[1]))/sin(veca[2]);
	  double temp = 1.0 - cos(veca[1])*cos(veca[1]) - 
		            pow((cos(veca[0]) - cos(veca[2])*cos(veca[1]))/sin(veca[2]), 2.0);
	  lattice_vectors[2][2] = vecm[2]*sqrt(temp);
      }

void crystal_system::Set_Amor_Vectors(double mesh_pam, double vari_pam)
     {
	 for (int a = 0; a < 3; ++a)
		 {
		 Zero_XYZ(lattice_vectors[a]);
		 }
     lattice_vectors[0][0] = mesh_pam;
     lattice_vectors[0][1] = vari_pam;                
     }

//Basis creation---
                         
void crystal_system::Add_Atom_To_Basis(const atom& ze_atom, const double* coors)
      {
      basis_atoms[basis_atoms.Size()] = ze_atom;
      basis_atoms[basis_atoms.Size()].Set_Atom_Location(coors);
      ++basis_atoms;                              
      }

void crystal_system::Add_Atom_To_Basis(const atom& ze_atom, double x, double y, double z)
	  {
	  basis_atoms[basis_atoms.Size()] = ze_atom;
      basis_atoms[basis_atoms.Size()].Set_Atom_Location(x, y, z);
      ++basis_atoms;   
	  }

void crystal_system::Add_Atom_To_Basis_Rel(const atom& ze_atom, const double* rel_coors)
      {
      double abs_coors[3];   
      Get_Abs_Coors(abs_coors, rel_coors);
      Add_Atom_To_Basis(ze_atom, abs_coors);
      }

//Set specific crystal system---

void crystal_system::Set_SC_System(const atom& ze_atom, double lattice_pam)
 //SC = one atom at the basis origin
     {
     Add_Atom_To_Basis(ze_atom, 0.0, 0.0, 0.0);
     Set_Cubic_Vectors(lattice_pam);
     } 

void crystal_system::Set_BCC_System(const atom& ze_atom, double lattice_pam)
 //BCC = two atoms, one at the basis origin and one at (0.5, 0.5, 0.5)*lattice_pam.
     {
     Add_Atom_To_Basis(ze_atom, 0.0, 0.0, 0.0);
     double atom2_loc = lattice_pam/2.0;
     Add_Atom_To_Basis(ze_atom, atom2_loc, atom2_loc, atom2_loc);
     Set_Cubic_Vectors(lattice_pam);                    
     }

void crystal_system::Set_FCC_System(const atom& ze_atom, double lattice_pam)
 //FCC = four atoms at (0, 0, 0), (0.5, 0.5, 0.0), 
 //(0.5, 0.0, 0.5), and (0.0, 0.5, 0.5)*lattice_pam.
     {
     Add_Atom_To_Basis(ze_atom, 0.0, 0.0, 0.0);
     double atom234_loc = lattice_pam/2.0;
     Add_Atom_To_Basis(ze_atom, atom234_loc, atom234_loc, 0.0);
     Add_Atom_To_Basis(ze_atom, atom234_loc, 0.0, atom234_loc);
     Add_Atom_To_Basis(ze_atom, 0.0, atom234_loc, atom234_loc);
     Set_Cubic_Vectors(lattice_pam);                       
     }

//Information retriveal---

int crystal_system::Get_Basis_Size() const
     {
     return basis_atoms.Size();
     }

bool crystal_system::Basis_Contains_Charged_Atom() const
     {
     bool charged_atom = FALSEV;
     for (int a = 0; a < basis_atoms.Size(); ++a)
         {
         if (!Check_FP_Equality(basis_atoms[a].Get_Atom_Charge(), 0.0))
		  //Look for at least one atom with non-zero charge.
            { 
            charged_atom = TRUEV;
            a = basis_atoms.Size();                        
            }     
         }
     return charged_atom;                                       
     }

void crystal_system::Get_Vectors(double* parameters) const
     {
      for (int a = 0; a < 9; ++a)
          {
          parameters[a] = lattice_vectors[a / 3][a % 3];    
		   //Relate the 3 1-D matrices back to the single 1-D
		   //matrix representation.
          }                           
     } 

void crystal_system::Get_Vectors(double* vec1, double* vec2, double* vec3) const
     {
     for (int a = 0; a < 3; ++a)
          {
          vec1[a] = lattice_vectors[0][a];   
          vec2[a] = lattice_vectors[1][a];  
          vec3[a] = lattice_vectors[2][a];  
          }                           
     }   
     
void crystal_system::Get_Vectors_Mag(double* magnitudes) const
     {
     for (int a = 0; a < 3; ++a)
         {
		 Get_Vector_Mag(magnitudes[a], a);
         }                                        
     }  

void crystal_system::Get_Vector_Mag(double& magnitude, int index) const
	 {
	 magnitude = Get_VecMag(lattice_vectors[index]);
	 }

double crystal_system::Size_A_Step() const
     {
     double A_size = Get_VecMag(lattice_vectors[0]);  
     return A_size;                             
     }
 
double crystal_system::Size_B_Step() const
     {
     double B_size = Get_VecMag(lattice_vectors[1]);  
     return B_size;                             
     } 
 
double crystal_system::Size_C_Step() const
     {
     double C_size = Get_VecMag(lattice_vectors[2]);  
     return C_size;                             
     }
     
double crystal_system::Biggest_Step() const
     {
     double A_size = Size_A_Step();
     double B_size = Size_B_Step();
     double C_size = Size_C_Step();
     double biggest_size = Get_Largest(A_size, B_size, C_size);
     return biggest_size;                             
     }
     
double crystal_system::Smallest_NonZero_Step() const
     {
     double A_size = Size_A_Step();
     double B_size = Size_B_Step();
     double C_size = Size_C_Step();
     double smallest_size = Get_Smallest_NonZero(A_size, B_size, C_size);
     return smallest_size;  
     }     

//Lattice vector tests/operations---

bool crystal_system::Is_OrthoRhombic() const
	{
	bool is_ortho = TRUEV;
	for (int a = 0; a < 0; ++a)
		{
		for (int b = 0; b < 3; ++b)
			{
			if (!(a == b) && !Check_FP_Equality(lattice_vectors[a][b], 0.0))
		     //Check off-diagonal components of the lattice vector matrix, which
			 //should all be equal to zero, if the system is orthogonal.
				{
				is_ortho = FALSEV;
				a = 3;
				b = 3;
				}
			}
		}
	return is_ortho;
	}	

bool crystal_system::Is_Full_Step_Length(double length, int index) const
	{
	bool is_full = FALSEV;
	double step_length;
	Get_Vector_Mag(step_length, index);
	 //Get magnitude of vector being considered.
	double ratio = length/step_length;
	 //Determine ratio of the passed parameter to that magnitude.
	int n_int = int(ratio);
	double multiple_base = double(n_int) * step_length;
	 //Determine the integer multiple of the magnitude that is
	 //closest in value to the passed parameter.
	double dif = length - multiple_base;
	 //Get the difference from that integer multiple.
	 //Criterion: Result = TRUEV if within 1% of magnitude.
	double step_frac = dif/step_length;
	if ( (step_frac < 0.01) || (step_frac > 0.99) )
	   {
	   is_full = TRUEV;
	   }
	return is_full;
	}

bool crystal_system::Is_Half_Step_Length(double length, int index) const
	{
	bool is_half = FALSEV;
	double step_length;
	Get_Vector_Mag(step_length, index);
	 //Get magnitude of vector being considered.
	double ratio = length/step_length;
	 //Determine ratio of the passed parameter to that magnitude.
	int n_int = int(ratio);
	double multiple_base = double(n_int) * step_length;
	 //Determine the integer multiple of the magnitude that is
	 //closest in value to the passed parameter.
	double dif = length - multiple_base;
	 //Get the difference from that integer multiple.
	 //Criterion: Result = TRUEV if within 1% of half-magnitude.
	double step_frac = dif/step_length;
	if ( (step_frac > 0.49) && (step_frac < 0.51) )
	   {
	   is_half = TRUEV;
	   }
	return is_half;
	}

void crystal_system::Get_Abs_Coors(double* abs_coors, const double* rel_coors) const
     {
	 Convert_RelCoor_To_AbsCoor(abs_coors, rel_coors, lattice_vectors[0], 
		                        lattice_vectors[1], lattice_vectors[2]);
     }

//Basis addition---

void crystal_system::Get_Basis_Shift(double* shift, const double* new_pos) const
	 {
	 if (basis_atoms.Size() > 0)
	    {
		double temps[3];
		basis_atoms[0].Get_Atom_Location(temps);
		SetSub_XYZ(shift, new_pos, temps);
		 //Coordinate shift values determined here can move the first atom
		 //in the basis to be at the passed position.
	    }
	 else
	  //Set shift vector to zero if there no basis atoms.
	    {
	    Zero_XYZ(shift);
	    }
	 }

void crystal_system::Partial_Basis_Check(atom_collection& ze_atom_array, int original_size, 
	                                     int number_added, bool partial_add) const
	 {
     if ( !partial_add && (number_added != 0) && (number_added != basis_atoms.Size()) )
      //If adding the whole basis is a requirement, check to make sure
      //the whole basis was added. If not, undo any atom additions.
          {
          ze_atom_array.Delete_Atoms(original_size, ze_atom_array.Size());          
          }  
	 }
     
void crystal_system::Add_Basis(atom_collection& ze_atom_array, const double* pos) const
     {
	 double shift_vals[3];
	 Get_Basis_Shift(shift_vals, pos);
	 for (int a = 0; a < basis_atoms.Size(); ++a)
           {
		   ze_atom_array.Add_Atom(basis_atoms[a], shift_vals, SHIFT_MODE); 
           }
     }

void crystal_system::Add_Basis_ToPointSet(atom_collection& ze_atom_array, 
	                                      const double** point_set, int& index) const
     {
     for (int a = 0; a < basis_atoms.Size(); ++a)
         {
		 ze_atom_array.Add_Atom(basis_atoms[a], point_set[index], PLACE_MODE);
         ++index;
         }           
     }

void crystal_system::Add_Basis_ToPointSet_ImposeSphere(atom_collection& ze_atom_array,
	                                       const double** point_set, int& index, double dist_min, 
										   double dist_max, bool partial_add) const
     {
     double dist_temp;
	 int original_size = ze_atom_array.Size();
	  //Size of atomic collection before any basis addition.
	 int number_added = 0;
	  //Number of basis atoms successfully added.
     for (int a = 0; a < basis_atoms.Size(); ++a)
         {
         dist_temp = Get_Distance_FromOrigin(point_set[index]);  
         if (Check_FP_LessOrEqual(dist_temp, dist_max) &&
			 Check_FP_GreaterOrEqual(dist_temp, dist_min))
		  //Spherical minimum/maximum distance enforcement.
			{
			ze_atom_array.Add_Atom(basis_atoms[a], point_set[index], PLACE_MODE);
			++number_added;
		    }
		 ++index;
         }         

	 Partial_Basis_Check(ze_atom_array, original_size, number_added, partial_add);
	  //Remove partial basis set addition if appropriate.
     }



void crystal_system::Add_Basis_ToSParticle(atom_collection& ze_atom_array, const double* pos, 
                                           double dist_min, double dist_max,
										   bool partial_add) const
     {
     double dist_temp;
	 int original_size = ze_atom_array.Size();
	  //Size of atomic collection before any basis addition.
	 int number_added = 0;
	  //Number of basis atoms successfully added.

	 Add_Basis(ze_atom_array, pos);
	  //Add the basis and THEN check on which atoms added are outside
	  //of the spherical boundaries.
     for (int a = original_size; a < ze_atom_array.Size(); ++a)
         {
		 dist_temp = ze_atom_array[a].Get_Distance_To_Origin();
		  //Determine the distance of the atom from the origin of the system,
		  //which here is taken to be the center of the desired sphere.
         if (Check_FP_LessOrEqual(dist_temp, dist_max) &&
			 Check_FP_GreaterOrEqual(dist_temp, dist_min))
		  //Spherical minimum/maximum distance enforcement.
			{
			++number_added;
			 //Keep the atom, which is within spherical bounds.
		    }
		 else
		    {
			ze_atom_array.Delete_Atom(a);
			--a;
			 //Delete the atom just added.
			}
         }         

	 Partial_Basis_Check(ze_atom_array, original_size, number_added, partial_add);
	  //Remove partial basis set addition if appropriate.           
     }

void crystal_system::Add_Basis_ToFParticle(atom_collection& ze_atom_array, const double* pos, 
                                           double dist_min, double target_radius, 
										   bool partial_add) const
     {

	 //Get a geometric description of the planes that make up the faceted surface---

     int plane_indices[MAX_FACETING_PLANES][3];
     double plane_intersect[MAX_FACETING_PLANES][3];
     double normal_vecs[MAX_FACETING_PLANES][3];
     double scaling_vals[MAX_FACETING_PLANES];
     double normal_mag;
     for (int a = 0; a < num_faceting_planes; ++a)
      //Determine the normal vector of each plane that connects the plane 
	  //to the origin of the particle.
		 {
         Get_Faceting_Plane(a, plane_indices[a], scaling_vals[a]);
		  //Obtain the indices of the plane as well as its scaling factor.
         Determine_Normal_Vector(plane_indices[a], normal_vecs[a]);
		  //Translates the Miller indices of the plane to a Cartesian vector which
		  //is normal to the plane and intersects the origin. Next step is to scale
		  //this vector such that the vector terminates at the actual facet.
         normal_mag = target_radius*scaling_vals[a];
		  //Distance of point-of-intersection of plane with origin-connecting
		  //normal vector (i.e. how far away the plane is from the particle origin).
         Get_NewMag_Vec(normal_vecs[a], plane_intersect[a], normal_mag);
		  //Scale the plane's normal vector's magnitude such that the faceting
		  //plane at which the normal vector terminates is effectively at the 
		  //right distance from the origin.
         }
         
	 //Add the basis of atoms, rejecting any atoms that are outside the faceted surface---

     double atom_loc[3];
     double term_sum, dist_temp;
     bool on_origin_side;

	 int original_size = ze_atom_array.Size();
	 int number_added = 0;
	 Add_Basis(ze_atom_array, pos);
	  //Add the basis and THEN remove any atoms
	  //not residing in the multi-faceted surface.
     for (int a = original_size; a < ze_atom_array.Size(); ++a)
         { 
         ze_atom_array[a].Get_Atom_Location(atom_loc);  
         on_origin_side = TRUEV;
         for (int b = 0; b < num_faceting_planes; ++b)
             {
             term_sum = Solve_Plane_Equation(plane_indices[b], atom_loc, plane_intersect[b]);
               //For points in the plane/facet: h(x - x0) + k(y - y0) + l(z - z0) = 0. 
             if (term_sum > 0.0)
               //Values over 0 when solving the plane equation indicate that the atom
			   //is not on the particle origin side of the plane.
                {
                on_origin_side = FALSEV;   
                b = num_faceting_planes;       
                }  
             }

		 dist_temp = ze_atom_array[a].Get_Distance_To_Origin();  
         if (on_origin_side && Check_FP_GreaterOrEqual(dist_temp, dist_min))
            {
            ++number_added;
			 //Atom passes the faceting and distance minimum check.
            }
		 else
		    {
		    ze_atom_array.Delete_Atom(a);
		    --a;
			 //Delete the atom just added.
		    } 
         }                  
	 
	 Partial_Basis_Check(ze_atom_array, original_size, number_added, partial_add);
	  //Remove partial basis set addition if appropriate.        
     }

void crystal_system::Add_Basis_ToCParticle(atom_collection& ze_atom_array, 
                                           double dist_min, double dist_max,
										   bool partial_add) const
     {                                                          
     double mesh_pam = lattice_vectors[0][0];
	 double vari_pam = lattice_vectors[0][1];
	  //First two values of lattice vectors represent the mean distance between
	  //atoms and the variation in this distance within the amorphous structure.
	 double min_mesh_pam = lattice_vectors[0][0] * (1.0 - vari_pam);
	 double max_mesh_pam = lattice_vectors[0][0] * (1.0 + vari_pam);
	  /*When adding an atom nearby another atom in the "growth-like" creation of
	    the amorphous nanoparticle, these represent the minimum and maximum distances
	    of that atom from the other atom. The minimum distance is also the lowest
	    atom-atom bond distance that will be allowed in the final structure*/
     double start_pos[3] = { 0.0, 0.0, 0.0 };  
	 int original_size = ze_atom_array.Size();
     if (original_size > 0) 
      //Look for an atom to place the new basis of atoms nearby.
        {   
        int rand_atom_index = GetRandPosNum(long int(original_size - 1));
        ze_atom_array[rand_atom_index].Get_Atom_Location(start_pos);
        }
	 int number_added = 0;
     double test_pos[3], move_coors[3], final_move[3];
     double bond_length, dist_temp;
	 bool overlap_check;

	 Add_Basis(ze_atom_array, start_pos);
	  //Add the basis and THEN reposition the atoms to create the
	  //random-close packing configuration.
     for (int a = original_size; a < ze_atom_array.Size(); ++a)
		 {
		 //Determine the semi-random position to place each atom of the
		 //basis at. To do this, a relative displacement vector from
		 //the test starting position must be determined.

         move_coors[0] = GetRandDec(-1.0, 1.0);
         move_coors[1] = GetRandDec(-1.0, 1.0);
         move_coors[2] = GetRandDec(-1.0, 1.0);
		  //Determine a random direction.

		 if (!Check_FP_Equality(vari_pam, 0.0))
		  //Determine a semi-random magnitude.
			{
			bond_length = GetRandDec(min_mesh_pam, max_mesh_pam);
			}
		 else
		    {
			bond_length = mesh_pam;
		    }

		 Get_NewMag_Vec(move_coors, final_move, bond_length);
		  //Combine magnitude and direction to get the displacement
		  //vector (i.e. where to put the basis atom at).
		 SetAdd_XYZ(test_pos, start_pos, final_move);
		  //Get the new location.
         ze_atom_array[a].Set_Atom_Location(test_pos);   
          //Set the basis atom to that location.
		 dist_temp = ze_atom_array[a].Get_Distance_To_Origin();  
		 overlap_check = TRUEV;
		 for (int b = 0; b < a; ++b)
			 {
			 if (ze_atom_array[b].Get_Distance(ze_atom_array[a]) < min_mesh_pam)
				{
				overlap_check = FALSEV;
				b = ze_atom_array.Size();
				}
			 }
         if (Check_FP_LessOrEqual(dist_temp, dist_max) 
			 && Check_FP_GreaterOrEqual(dist_temp, dist_min)
             && overlap_check)
            { 
            ++number_added;
            }
		 else
		    {
            ze_atom_array.Delete_Atom(a);
		    --a;
			 //Delete the atom just added.
		    }
         }    

     Partial_Basis_Check(ze_atom_array, original_size, number_added, partial_add);
	  //Remove partial basis set addition if appropriate.      
     }

void crystal_system::Add_Basis_ToSurface(atom_collection& ze_atom_array, const double* pos,       
                                         const double* dimensions, 
                                         bool partial_add) const
     { 
	 bool bounds_check;
	 double new_loc[3];
	 int original_size = ze_atom_array.Size();
	  //Size of atomic collection before any basis addition.
	 int number_added = 0;
	  //Number of basis atoms successfully added.

	 Add_Basis(ze_atom_array, pos);
	  //Add the basis and THEN check on which atoms added are outside
	  //of the orthorhombic boundaries.
     for (int a = original_size; a < ze_atom_array.Size(); ++a)
	  //For each atom just added to the atomic collection...
        {
        ze_atom_array[a].Get_Atom_Location(new_loc);
        bounds_check = Check_BoundariesZeroToMax(new_loc, dimensions);
		 //Check to see if atom is within the orthorhombic box.
        if (bounds_check)
		 //Atom is inside the orthorhombic bounds.
           {
		   ++number_added;
           }
		else
		   {
           ze_atom_array.Delete_Atom(a);
		   --a;
			 //Delete the atom just added.
		   }
        }   
	 Partial_Basis_Check(ze_atom_array, original_size, number_added, partial_add);
	  //Remove partial basis set addition if appropriate.                
     } 

//Lattice vector stepping---

void crystal_system::A_Step(double* coor) const
     {
	 Add_XYZ(coor, lattice_vectors[0]);
     } 

void crystal_system::B_Step(double* coor) const
     {   
     Add_XYZ(coor, lattice_vectors[1]);                 
     }

void crystal_system::C_Step(double* coor) const
     {
     Add_XYZ(coor, lattice_vectors[2]);              
     }
     
void crystal_system::A_Step(double* coor, int num_steps) const
     {
	 double total_step[3];
	 Set_XYZ(total_step, lattice_vectors[0]);
	 Multiply_XYZ(total_step, double(num_steps));
	 Add_XYZ(coor, total_step);   
     } 

void crystal_system::B_Step(double* coor, int num_steps) const
     {
     double total_step[3];
	 Set_XYZ(total_step, lattice_vectors[1]);
	 Multiply_XYZ(total_step, double(num_steps));
	 Add_XYZ(coor, total_step);        
     } 

void crystal_system::C_Step(double* coor, int num_steps) const
     {
	 double total_step[3];
	 Set_XYZ(total_step, lattice_vectors[2]);
	 Multiply_XYZ(total_step, double(num_steps));
	 Add_XYZ(coor, total_step);         
     } 

void crystal_system::Multi_Step(double* coor, int a, int b, int c) const
     {
     A_Step(coor, a);
     B_Step(coor, b);
     C_Step(coor, c);                                    
     }

void crystal_system::Multi_Step(double* coor, int* step_counts) const
     {
     A_Step(coor, step_counts[0]);
     B_Step(coor, step_counts[1]);
     C_Step(coor, step_counts[2]);                                    
     }

//Cell rotation---

void crystal_system::Rotate_Planes(const int* dest_plane, const int* sou_plane)
     {     
	 double old_normal[3];
     Determine_Normal_Vector(sou_plane, old_normal);
     double new_normal[3];
     Determine_Normal_Vector(dest_plane, new_normal); 
	  //Determine normal vectors for the two planes involved in rotation.
     double rot_vec[3];
     double rot_angle;
     Calculate_Rotation_Parameters(old_normal, new_normal, rot_vec, rot_angle);
      //Use plane normal vectors to determine the rotation vector and angle
      //needed.                        
     Rotate_Cell(rot_vec, rot_angle);
     } 


void crystal_system::Rotate_Cell(const double* rotation_vector, double angle)
     {
     double rot_matrix[9];
     Calc_Rotation_Matrix(rotation_vector, angle, rot_matrix);                                   
      //Calculates the rotation matrix for the following rotation operation.
     
     //Rotate the lattice vectors---
      
     Rotate_Vector(rot_matrix, lattice_vectors[0]); 
     Rotate_Vector(rot_matrix, lattice_vectors[1]); 
     Rotate_Vector(rot_matrix, lattice_vectors[2]);  
 
	 //Rotate the basis atom location vectors
	 //in the same way---

     double vec_temp[3];
     for (int a = 0; a < basis_atoms.Size(); ++a)
         {
         basis_atoms[a].Get_Atom_Location(vec_temp);
         Rotate_Vector(rot_matrix, vec_temp);
         basis_atoms[a].Set_Atom_Location(vec_temp);
         }
     }

//Faceting plane determination---     
     
int crystal_system::Get_Facet_Count() const
    {
    return num_faceting_planes;                       
    }    

void crystal_system::Set_Faceting_Plane(int index1, int index2, int index3, 
                                        double scaling_factor)
     {
     faceting_planes[num_faceting_planes][0] = index1;
     faceting_planes[num_faceting_planes][1] = index2;
     faceting_planes[num_faceting_planes][2] = index3;   
     scaling_factors[num_faceting_planes] = scaling_factor;
     ++num_faceting_planes;           
	 if (num_faceting_planes == MAX_FACETING_PLANES)
		{
		Show_Warning("TOO MANY FACETING PLANES DEFINED!");
		}
     }  

void crystal_system::Get_Faceting_Plane(int rank, int* indices, double& scale) const
     {
     if (num_faceting_planes == 0)
        {    
		Show_Warning("FACETS FOR CRYSTAL SYSTEM NOT DECLARED!");
        }
	 indices[0] = faceting_planes[rank][0];
	 indices[1] = faceting_planes[rank][1];
	 indices[2] = faceting_planes[rank][2];
     scale = scaling_factors[rank];                        
     }

void crystal_system::Determine_Normal_Vector(const int* ind, double* vec) const
     {
	 double rel_comp[3];
	 Zero_XYZ(rel_comp);
	 for (int a = 0; a < 3; ++a)
		 {
		 if (ind[a] != 0)
	      //Miller index of 0 = no component in normal vector. Otherwise,
	      //the inverse of the Miller index must be used.
			{
			rel_comp[a] = 1.0/double(ind[a]);
			}
		 }
	 Get_Abs_Coors(vec, rel_comp);                                                                
     }   
 
     
//Supercell logic---

int crystal_system::Add_SuperCell(atom_collection& ze_collection, int x, int y, int z) const
     {
     int atom_count = 0;
     double pos[3];

     for (int a = 0; a < x; ++a)
         {
         for (int b = 0; b < y; ++b)
             {
             for (int c = 0; c < z; ++c)
                 {
                 Zero_XYZ(pos);
                 Multi_Step(pos, a, b, c);
                 Add_Basis(ze_collection, pos);     
				 atom_count += basis_atoms.Size();
                 }     
             }     
         } 

     return atom_count;                     
     } 

//Operator overloads---

void crystal_system::operator = (const crystal_system& ze_syst)
     {        
	 //Assign lattice vectors---

     for (int a = 0; a < 3; ++a)
        {
	    Set_XYZ(lattice_vectors[a], ze_syst.lattice_vectors[a]);
        }

	 //Assign basis atoms---

	 basis_atoms.Copy_Atomic_Group(ze_syst.basis_atoms);

	 //Assign faceting planes---
     
     num_faceting_planes = ze_syst.num_faceting_planes;
     for (int a = 0; a < num_faceting_planes; ++a)
        {
		faceting_planes[a][0] = ze_syst.faceting_planes[a][0];
		faceting_planes[a][1] = ze_syst.faceting_planes[a][1];
		faceting_planes[a][2] = ze_syst.faceting_planes[a][2];
        scaling_factors[a] = ze_syst.scaling_factors[a]; 
        }
     }

//File I/O---

void crystal_system::Save_Crystal_System(const char* file_name) const
     {
     ofstream save_file(file_name, ios::out);
     if (save_file.is_open())
        {
        Save_Crystal_System(save_file);   
		save_file.close();  
        }                          
     }

void crystal_system::Load_Crystal_System(const char* file_name)
     {       
     ifstream load_file(file_name, ios::in);
     if (load_file.is_open())
        {
        Load_Crystal_System(load_file); 
		load_file.close();
        }
     }
      
void crystal_system::Save_Crystal_System(ofstream& save_file) const
     {
     Write_Intro(save_file, "Crystal Structure Information");
     
	 //Store lattice vectors and basis atoms---
	 
	 for (int a = 0; a < 3; ++a)
         {
         Write_Values(save_file, lattice_vectors[a], 3);
         save_file << endl;
         }
     basis_atoms.Atoms_Storage(save_file);   
     
	 //Store faceting planes---

     save_file << endl << endl << "Planes for Faceting:" << endl;
     save_file << num_faceting_planes;
     for (int a = 0; a < num_faceting_planes; ++a)
         {
         save_file << endl;
		 Write_Values(save_file, faceting_planes[a], 3);
         save_file << " " << scaling_factors[a];
         }               
     }   


void crystal_system::Load_Crystal_System(ifstream& load_file)
     {   
     Skip_Phrases(load_file, 3);

	 //Load lattice vectors and basis atoms---

     for (int a = 0; a < 3; ++a)
         {
         Load_Values(load_file, lattice_vectors[a], 3);
         }
     basis_atoms.Atoms_Retrieval(load_file);
     
	 //Load faceting planes---

     Skip_Phrases(load_file, 3);
     load_file >> num_faceting_planes;
     for (int a = 0; a < num_faceting_planes; ++a)
          {
		  Load_Values(load_file, faceting_planes[a], 3);
          load_file >> scaling_factors[a];   
          } 
     }

void crystal_system::Save_SuperCell(const char* file_name, int x, int y, 
                                    int z) const
     {
     atom_collection temp_atoms(x*y*z*basis_atoms.Size());
     Add_SuperCell(temp_atoms, x, y, z);
	  //Create the desired supercell, reserving an initial
	  //memory that is equal to the expected number of
	  //atoms in the supercell.

	 //Store the supercell as a list of atoms with Cartesian locations---

     ofstream ze_file(file_name, ios::out);
     if (ze_file.is_open())
        {
        temp_atoms.Store_Atoms_Location(ze_file, TRUEV);
		ze_file.close();   
        }           
     } 

//Automatic close-packed plane calculation (not currently in usage)---

void crystal_system::Calculate_ClosePacked_Planes()
     { 
     const int CPP_LIMIT = MAX_AUTO_PLANE_INDEX*2 + 1;
	  /*Range of indices that can be investigated in the close-packed
	    plane search. E.g. if MAX_AUTO_PLANE_INDEX is 2, then the
	    search index can be -2, -1, 0, 1, or 2, hence 5. This goes
	    for each Miller index of each plane*/
     int packing_test_results[CPP_LIMIT][CPP_LIMIT][CPP_LIMIT];
	  //Stores number of atoms found in planes investigated.
	 int array_x, array_y, array_z;
	  //Array indices for storing results for planes that themselves
	  //have both positive and negative indices.
     
	 int plane_indices[3] = { 0, 0, 0 };
	  //Set of Miller indices.

	 double plane_normal[3];
	  //Stores the normal to a plane being investigated.
	 double temp_loc[3]; 
	  //Temporary atom location storage.
	 double term_sum;
	  //Stores the sum from solving the plane equation.
     
     //Test different planes to see how many atoms are in each plane---
     
     for (int h = -MAX_AUTO_PLANE_INDEX; h <= MAX_AUTO_PLANE_INDEX; ++h)
         {
         array_x = h + MAX_AUTO_PLANE_INDEX;
          //Plane indices go from negative MAX_AUTO_PLANE_INDEX to 
		  //positive MAX_AUTO_PLANE_INDEX in this search.
          //Array indices go alongside that from 0 to CPP_LIMIT - 1.
         for (int k = -MAX_AUTO_PLANE_INDEX; k <= MAX_AUTO_PLANE_INDEX; ++k)
             {
             array_y = k + MAX_AUTO_PLANE_INDEX;
             for (int l = -MAX_AUTO_PLANE_INDEX; l <= MAX_AUTO_PLANE_INDEX; ++l)
                 {
                 array_z = l + MAX_AUTO_PLANE_INDEX;
			     if ( (h == 0) && (k == 0) && (l == 0) )
                  //The (000) plane does not exist. That's just weird
				  //to think about.
                    {
                    continue;     
                    }

                 packing_test_results[array_x][array_y][array_z] = 0;
				  //Initialize packing test results.
				 plane_indices[0] = h;
				 plane_indices[1] = k;
				 plane_indices[2] = l;
                 Determine_Normal_Vector(plane_indices, plane_normal);   
                 for (int a = 0; a < basis_atoms.Size(); ++a)
					 {
					 basis_atoms[a].Get_Atom_Location(temp_loc);
				     term_sum = Solve_Plane_Equation(plane_indices, temp_loc, plane_normal);
                     if (Check_FP_Equality(term_sum, 0.0))  
                      //Found an atom in the plane being considered. A more thorough version
					  //of this function would consider other instances of this atom.
                        {
                        ++packing_test_results[array_x][array_y][array_z];  
                        }                        
                     }
                 }     
             }     
         }
    
     //Order the planes by ranking (most atoms in planes = lowest rank)---
    
     int current_rank = 0; 
     int current_best_closepack[3] = { 0, 0, 0 };
     int max_packing, packing_temp;
     while (current_rank < MAX_FACETING_PLANES)
           {
           max_packing = 0;
           for (int array_x = 0; array_x < CPP_LIMIT; ++array_x)
		    //Find the plane with the most atoms packed in it.
			//(that has not been accounted for).
               {
               for (int array_y = 0; array_y < CPP_LIMIT; ++array_y)
                   {
                   for (int array_z = 0; array_z < CPP_LIMIT; ++array_z)
                       {
                       packing_temp = packing_test_results[array_x][array_y][array_z];
                       if (max_packing <= packing_temp)
                           {
                           current_best_closepack[0] = array_x;
                           current_best_closepack[1] = array_y;
                           current_best_closepack[2] = array_z;  
                           max_packing = packing_test_results[array_x][array_y][array_z];                                         
                           }                 
                       }
                   }     
               }   

		   if (max_packing == 0)
			//No more packing planes to consider. Time to leave.
			   {
			   break;
			   }

           faceting_planes[current_rank][0] = current_best_closepack[0] - MAX_AUTO_PLANE_INDEX;
           faceting_planes[current_rank][1] = current_best_closepack[1] - MAX_AUTO_PLANE_INDEX;
           faceting_planes[current_rank][2] = current_best_closepack[2] - MAX_AUTO_PLANE_INDEX;
		    //Array indices must be shifted to represent the actual Miller indices.
           scaling_factors[current_rank] = 1.0;
		    //Default scaling factor = 1.0.
           packing_test_results[current_best_closepack[0]]
		                       [current_best_closepack[1]][current_best_closepack[2]] = 0;
            //Don't want to consider this plane a second time!
           ++current_rank; 
           num_faceting_planes = current_rank;
           } 
     }