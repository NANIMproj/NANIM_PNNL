
#include "stdafx.h"

#include "DWN_nanop.h"

//Constructors, destructrors, and initialization functions---

nanoparticle::nanoparticle()
      {
      Initialize();
      }

nanoparticle::nanoparticle(double target_diameter, double hole_diameter, int model, 
                           const crystal_system& ze_cryst_syst, bool stuffed,
						   bool neut_part)
      {
      Initialize();
      Set_Properties(target_diameter, hole_diameter, model, 
		             ze_cryst_syst, stuffed, neut_part);
      }

nanoparticle::~nanoparticle()
      {
      }

void nanoparticle::Initialize()
	  {
      particle_atoms.Set_Initial_Allocation(INITIAL_PATOM_COUNT); 
	  model_type = SPH_MODEL;
      effective_diameter = 0.0;
      effective_radius = 0.0;    
      stuff_space = TRUEV;
	  neutralize_particle = FALSEV;
	  }

//Properties setting---

void nanoparticle::Set_Properties(double target_diameter, double pore_diameter,
								  int model, const crystal_system& ze_cryst_syst, 
                                  bool stuffed, bool neut_part)
     {
     effective_diameter = target_diameter;
     effective_radius = target_diameter/2.0;
	 hole_diameter = pore_diameter;
     model_type = model;        
     stuff_space = stuffed;
     cryst_struct = ze_cryst_syst;
	 neutralize_particle = neut_part;
     }

//Nanoparticle construction---

void nanoparticle::Construct_Particle()
	 {
     if (model_type != CON_MODEL)
        {
        Construct_Cryst_Particle();
        }
     else
        {
        Construct_Conti_Particle();                           
        }
	 if (neutralize_particle)
		{
		Neutralize_Particle();
		}
     }

void nanoparticle::Construct_Cryst_Particle()
      {    
      //Get sampling of spatial points to place crystal basis atoms at---

	  double lattice_vecs[3][3];
	  cryst_struct.Get_Vectors(lattice_vecs[0], lattice_vecs[1], lattice_vecs[2]);
	  
	  double largest_vector = cryst_struct.Biggest_Step();
	  double extra_space_parameter = 2.0*largest_vector;
	  if (model_type == CRY_MODEL)
	   //The needed sampling size is harder to pin down with a faceted particle, 
	   //so use a really large sampling volume just to be safe.
	     {
         extra_space_parameter += effective_radius; 
	     }
	  double test_radius = effective_radius + extra_space_parameter;

	  double** points;
	  int point_index = 0;
	  double volume = (4.0/3.0)*PI_CONST*(test_radius*test_radius*test_radius);
	  const double MEM_SCALING_FACTOR = 10000.0/8.0;
	   //Assume that a 8 nm^3 particle would not require more than 10000 lattice points.
	  int max_points = int(MEM_SCALING_FACTOR * (volume + 1.0));
	  Get_Memory(points, max_points, 3);
	  Get_Spherical_Sampling(lattice_vecs[0], lattice_vecs[1], lattice_vecs[2], 
		                     test_radius, points, max_points, point_index);
	    
	  //For faceting nanoparticles, basis positions may be aligned to match 
	  //specific sizes requested by the user. So far, the addition of an odd
	  //number of planes in typical orthorhombics has been figured out there---

	  if (model_type == CRY_MODEL)
	     {
		 if (cryst_struct.Is_OrthoRhombic())
		    {
			double quarter_step;
			for (int a = 0; a < 3; ++a)
			 //Look out for odd # of planes along each dimension---
			    {
				if (cryst_struct.Is_Half_Step_Length(effective_diameter, a))
				 //If requested nanoparticle size = (n + 1/2) * lattice vector,
				 //this may be a request for an odd number of total planes.
				   {
				   cryst_struct.Get_Vector_Mag(quarter_step, a);
				   quarter_step /= 4.0;
				    //If diameter = half step larger in size, radius
				    //is a quarter step larger in size.
				   for (int b = 0; b < point_index; ++b)
					   {
					   points[b][a] -= quarter_step;
				       }
				   }
			    }
		    }
	     }

	  //Add basis atoms to the particle's atomic collection---
  
	  for (int a = 0; a < point_index; ++a)
	     {
         if (model_type == SPH_MODEL)
            {
            cryst_struct.Add_Basis_ToSParticle(particle_atoms, points[a], 
											   hole_diameter/2.0, effective_radius,
											   stuff_space);
               //Add a basis of atoms at that lattice point, not allowing
               //any atoms outside the desired spherical region.
            }
         else if (model_type == CRY_MODEL)
            {
            cryst_struct.Add_Basis_ToFParticle(particle_atoms, points[a],
											   hole_diameter/2.0, effective_radius, 
											   stuff_space);
               //Add a basis of atoms at that lattice point, not allowing
               //any atoms outside the internal volume of the faceted region.   
            }
	     }
	  Zero_Center();
	   //Place nanoparticle at standard position.
	  Free_Memory(points, max_points);
      }

void nanoparticle::Construct_Conti_Particle()
      {
      //Get spatial sampling for the continuous network---

      double lattice_vecs[3][3];
	  cryst_struct.Get_Vectors(lattice_vecs[0], lattice_vecs[1], lattice_vecs[2]);
	  double mesh_pam = lattice_vecs[0][0];
	  double vari_pam = lattice_vecs[0][1];
	   //Get mesh grid (atom-atom spacing) and variance parameters. 

	  double** points;
	  long int point_index = 0;
	  double volume = (4.0/3.0)*PI_CONST*(effective_radius*effective_radius*effective_radius);
	  double mem_scaling_factor = 40000.0/8.0;
	   //Assume that a 8 nm^3 particle would not require more than 40000 atomic locations,
	   //not basis points.
	  const int MAX_TRIES = 10000;
	   //Maximum number of tries in expanding a continuous network of points.
	  int max_points = int(mem_scaling_factor * (volume + 1.0));
	  Get_Memory(points, max_points, 3);
	  Get_Random_SMesh_Sampling(mesh_pam, vari_pam, MAX_TRIES, effective_radius, 
			                    points, max_points, point_index);

	  //Construct the continuous nanoparticle---

      int current_size = 0;
      while (current_size < point_index)
            {
            cryst_struct.Add_Basis_ToPointSet_ImposeSphere(particle_atoms, 
				                              const_cast<const double**>(points), current_size,
											  hole_diameter/2.0, effective_radius,
											  stuff_space);  
            }       
	  Zero_Center();
	   //Place nanoparticle at standard position.
	  Free_Memory(points, max_points);
      }

void nanoparticle::Neutralize_Particle()
	  {
	  double sys_charge = particle_atoms.Return_Total_Charge();
	  Show_Statement("Particle charge before neutralization: ", sys_charge);
	  Zero_Center();

	  //Get first atom-removal test vector---

	  double rand_vec[3];
	  double direc_vec[3] = { 1.0, 1.0, 1.0 };
	  if ( (model_type != CRY_MODEL) || NO_FIRST_ATOM_METHOD)
	   //Crystalline nanoparticle: Remove first atom from
	   //corner.
	     {
		 Get_Rand_Direction(direc_vec, 3);
	     }

	  //Systematically remove ions of the right charge---

	  int test_atom_index;
	  int try_count = 0;
	  int MAX_TRY_COUNT = 5000;
	  while (abs(sys_charge) > 0.3)
	      {
		  Get_NewMag_Vec(direc_vec, rand_vec, effective_radius * 3.0);
		  test_atom_index = particle_atoms.Find_Closest_Atom(rand_vec);
		   //Closest atom to a vector that is outside of the nanoparticle
		   //surface can be linked to a closest atom that is on the surface.
		  ++try_count;
		  if (particle_atoms[test_atom_index].Get_Atom_Charge()/sys_charge > 0.0)
			 {
			 sys_charge -= particle_atoms[test_atom_index].Get_Atom_Charge();
			 particle_atoms.Delete_Atom(test_atom_index);
			 try_count = 0;
			 }

		  if (try_count == MAX_TRY_COUNT)
			 {
			 Show_Warning("PARTICLE NEUTRALIZATION FAILED!");
			 sys_charge = 0.0;
			 }
		  Get_Rand_Direction(direc_vec, 3);
	      }
      }

//Coating functions---

void nanoparticle::Coat_Spherical_Particle(const molecule& coat_mol, double part_head_dist, 
                                           double hh_dist, bool place_on_atom)
     {
     layer_molecule = coat_mol;
      //Record the coating molecule.
     layer_molecule.Standard_Orientation();
	 const double INITIAL_ORIE[3] =  { 0.0, 0.0, -1.0 };
	 Zero_Center();
     
     const double TWO_PI = 2.0*PI_CONST;
     double step_angle = hh_dist / effective_radius;
	  //Effective length of the arc used for molecule-molecule spacings
	  //at the surface, assuming a large surface.
     int num_steps = int(PI_CONST/step_angle); 
      //Number of steps to go along half a long circumference around the surface
	  //of the spherical nanoparticle.
     double breadth = double(num_steps - 1)*step_angle;
      //The exact angular range to be used in positioning the molecules.
     double start_angle = (PI_CONST - breadth)/2.0;
      //Sets the initial angle such that it will create a symmetrica
	  //angle distribution between 0 and PI.
     
     double second_step_angle, second_breadth, second_start_angle;
     int second_num_steps; 
     double current_radius;
     double theta, phi;
	  //In order to cover all points in the monolayer on the spherical
	  //shell, it is necessary to have two angular regions in
	  //the description. One goes from 0 to PI, the other 0 to 2*PI. 

	 int core_size = particle_atoms.Size();
	  //Number of atoms in the core nanoparticle.

     double start_pos[3];
     int nearest_atom_index;
     double nearest_atom_pos[3];
     double final_pos[3];
     int molecule_size = layer_molecule.Size();
	 int molecules_added = 0;
     for (int a = 0; a < num_steps; ++a)
         {
         theta = step_angle*double(a) + start_angle;
         current_radius = effective_radius * sin(theta);
          //This is the radius of the circumference
		  //in a specific slice of the sphere that is
		  //determined by theta.
         
         if (Check_FP_Equality(current_radius, 0.0))
            {
            continue;                                   
            }
         
         second_step_angle = hh_dist / current_radius;
         second_num_steps = int(TWO_PI/second_step_angle);
         second_breadth = double(second_num_steps - 1)*second_step_angle;
         second_start_angle = (TWO_PI - second_breadth)/2.0;
		  //Determine the number of steps to coat one shell of the sphere
		  //at the given value of theta.

         for (int b = 0; b < second_num_steps; ++b)     
             {
             phi = second_step_angle*double(b) + second_start_angle;     
             Get_Cartesian_Coordinates(effective_radius, theta, phi, start_pos);
			  //Determine the Cartesian location corresponding to
			  //the angles theta and phi.
             
             if (place_on_atom)
			  //Surface atom-molecule binding mode.
                {
                nearest_atom_index = particle_atoms.Find_Closest_Atom(0, core_size, start_pos);
                particle_atoms[nearest_atom_index].Get_Atom_Location(nearest_atom_pos);                 
                Set_XYZ(final_pos, nearest_atom_pos);
                }
             else
                {
                Set_XYZ(final_pos, start_pos);                             
                }
             
             Lengthen_Vec(final_pos, part_head_dist);
              //Add the surface-molecule distance to the final vector.
             layer_molecule.Copy_Molecule(particle_atoms, final_pos, INITIAL_ORIE, 
                                          final_pos, FIRST_ATOM);
			 int mol_start = particle_atoms.Size() - molecule_size;
			 int mol_end = particle_atoms.Size();
			 particle_atoms.Monolayer_Overlap_Prevention(core_size, mol_start, mol_start, 
				                                         mol_end, PI_CONST/18.0, 36);
              //Spatial overlap check. If spatial overlap can not be prevented, this
			  //molecule and the one it overlaps with will both be removed.
			 if (mol_end == particle_atoms.Size())
			    {
			    ++molecules_added;
				 //Added one molecule overall.
			    }
			 else
			    {
                --molecules_added;
				 //Removed one molecule overall.
			    } 
             }    
            
         }

	 Show_Statement("Molecule coating of spherical particle added ", molecules_added, " molecules.");
     }

void nanoparticle::Coat_Faceted_Particle(const molecule& coat_mol, double part_head_dist, 
                                         double hh_dist, bool place_on_atom)
     {
     layer_molecule = coat_mol;
	  //Record the coating molecule.
	 layer_molecule.Standard_Orientation();
	 double INITIAL_ORIE[3] = { 0.0, 0.0, -1.0 };
	 Zero_Center();

	 //Get crystal facet information---

     int plane_indices[MAX_FACETING_PLANES][3];
     double plane_intersect[MAX_FACETING_PLANES][3];
     double normal_vecs[MAX_FACETING_PLANES][3];
     double scaling_vals[MAX_FACETING_PLANES];
     double intersect_mag;
     for (int a = 0; a < cryst_struct.Get_Facet_Count(); ++a)
      //Calculate the effective origin in each of the planes that are defining the faceted surface.
         {
         cryst_struct.Get_Faceting_Plane(a, plane_indices[a], scaling_vals[a]);
         cryst_struct.Determine_Normal_Vector(plane_indices[a], normal_vecs[a]);
         intersect_mag = effective_radius*scaling_vals[a];
         Get_NewMag_Vec(normal_vecs[a], plane_intersect[a], intersect_mag);
		 }

	 //Reserve memory for spatial plane sampling---

     int max_1D_points = 1 + int( (effective_diameter)/hh_dist ) * 2;
      //Rough maximum number of molecules that could be added along
	  //one spatial direction along one specific facet.
     int num_points = max_1D_points*max_1D_points; 
	  //Points for a 2D grid on one facet.
	 double** points_to_place;
	 Get_Memory(points_to_place, num_points, 3); 
     
	 int core_size = particle_atoms.Size();
	  //Number of atoms in the core nanoparticle.

	 int nearest_atom_index;
     double term_sum;     
	 bool point_on_facet;
     
     double nearest_atom_pos[3];
     double final_pos[3];
     int molecule_size = layer_molecule.Size();
	 int molecule_count = 0;
     for (int a = 0; a < cryst_struct.Get_Facet_Count(); ++a)
	     {
         Get_Plane_Sampling(plane_indices[a], plane_intersect[a], hh_dist, 
                            points_to_place, num_points);
		  //Obtain a sampling of spatial points in one of the crystal
		  //facets (and just that facet's plane, other facets are not
		  //considered in the point sampling).
         for (int b = 0; b < num_points; ++b)
             {
             point_on_facet = TRUEV; 
             for (int c = 0; c < cryst_struct.Get_Facet_Count(); ++c)
              //If point in the sampling is outside the particle surface, it
              //will fail at least one of the plane equation tests for one 
              //of the other cutting planes.
                 {    
                 if (c == a)
                    {
                    continue;      
                    }
                 
                 term_sum = Solve_Plane_Equation(plane_indices[c], points_to_place[b], plane_intersect[c]);
                    // h(x - x0) + k(y - y0) + l(z - z0) = 0, for points in the plane. 
                 if (term_sum > 0.0)
                    //Values over 0 indicate that the atom is not on the origin-side of the faceting plane.
                    {
                    point_on_facet = FALSEV;  
                    break;       
                    }   
                 }
               
             if (!point_on_facet)
                {
                continue;                    
                }  

             if (place_on_atom)
                {
                nearest_atom_index = particle_atoms.Find_Closest_Atom(0, core_size, points_to_place[b]);    
                particle_atoms[nearest_atom_index].Get_Atom_Location(nearest_atom_pos);                 
                Set_XYZ(final_pos, nearest_atom_pos);
                }
             else
                {
                Set_XYZ(final_pos, points_to_place[b]);                             
                }
             
             Lengthen_Vec(final_pos, part_head_dist);
              //Lengthen the placement vector to account for surface-molecule spacing.
             layer_molecule.Copy_Molecule(particle_atoms, final_pos, INITIAL_ORIE, 
                                                          final_pos, FIRST_ATOM);
			 particle_atoms.Inc_Coll_Size(molecule_size);
			 int mol_start = particle_atoms.Size() - molecule_size;
			 int mol_end = particle_atoms.Size();
			 particle_atoms.Monolayer_Overlap_Prevention(core_size, mol_start, mol_start, 
														 mol_end, PI_CONST/18.0, 36);
			  //Spatial overlap check. If spatial overlap can not be prevented, this
			  //molecule and the one it overlaps with will both be removed.
			 if (mol_end == particle_atoms.Size())
			  //Molecule added.
			    {
                ++molecule_count;
			    }
			 else
			  //Molecule removed.
			    {
                --molecule_count; 
				}
			 }  
         }       
    Free_Memory(points_to_place, num_points);  
	
	Show_Statement("Molecule coating of faceted particle added ", molecule_count, " molecules.");
    }

//Utility functions---

void nanoparticle::Get_Center(double* cen_coors) const
    {
    particle_atoms.Coordinate_Average(cen_coors);              
    }

void nanoparticle::Translate_Particle(const double* shift_coors)
    {
    particle_atoms.Shift_Atomic_Coordinates(shift_coors);      
    }

void nanoparticle::Zero_Center()
	{
    const double ZERO_VEC[3] = { 0.0, 0.0, 0.0 };
	ReAssign_Center(ZERO_VEC);
	}

void nanoparticle::ReAssign_Center(const double* new_center_coors)
    {
    particle_atoms.Translate_To_New_Coordinate_Average(new_center_coors);                                           
    }
                         
void nanoparticle::Rotate(const double* angles)
    {
	double temp_pos[3];
	Get_Center(temp_pos);
	particle_atoms.Rotate_AtomsXYZ(angles, 1, temp_pos);
    }

//Nanoparticle copying---

int nanoparticle::Copy_NanoParticle(atom_collection& atomic_array) const
    {
    atomic_array.Copy_Atomic_Group(particle_atoms);         
    return particle_atoms.Size();                         
    }

int nanoparticle::Copy_NanoParticle(atom_collection& atomic_array,
                                    const double* position) const
    {
    double shift_coors[3];
    double center[3];
    Get_Center(center);
    Set_XYZ(shift_coors, position);
    Sub_XYZ(shift_coors, center);
    atomic_array.Copy_Atomic_Group(particle_atoms, atomic_array.Size(), 0, 
	                               particle_atoms.Size(), shift_coors);   
    return particle_atoms.Size();                 
    }

int nanoparticle::Copy_NanoParticle(atom_collection& atomic_array,
	                                const double* position,
									const double* rot_angles) const
    {
    double shift_coors[3];
    double center[3];
    Get_Center(center);
    Set_XYZ(shift_coors, position);
    Sub_XYZ(shift_coors, center);
    atomic_array.Copy_Atomic_Group(particle_atoms, atomic_array.Size(), 0, 
	                               particle_atoms.Size(), shift_coors, rot_angles);    
    return particle_atoms.Size(); 
    }

//I/O---

void nanoparticle::Print_Atom_Locations() const
    {
    cout << endl << "HERE COMES THE MIGHTY FLOW OF INFORMATION!!!!" << endl;
    particle_atoms.Print_Atoms_Location();
    }

void nanoparticle::Store_Atom_Locations(ofstream& out_file) const
    {
    particle_atoms.Store_Atoms_Location(out_file, TRUEV);
    }
   
void nanoparticle::Save_NanoParticle(const char* save_file) const
    {
    ofstream ze_file(save_file, ios::out);

    if (ze_file.is_open())
      {
      Write_Intro(ze_file, "File Type = NanoParticle Description");
      ze_file << effective_diameter << ' ' << hole_diameter << ' ' 
		      << model_type << ' '  << stuff_space << ' ' 
			  << neutralize_particle << endl;
      cryst_struct.Save_Crystal_System(ze_file);
      ze_file << endl << endl << "Atomic Information List";
      particle_atoms.Atoms_Storage(ze_file);  
	  ze_file.close();
      }
    }
   
void nanoparticle::Load_NanoParticle(const char* load_file)
    {             
    ifstream ze_file(load_file, ios::in);
    if (ze_file.is_open())
      {
      Skip_Phrases(ze_file, 5);
      ze_file >> effective_diameter >> hole_diameter >> model_type 
		      >> stuff_space >> neutralize_particle;
      effective_radius = effective_diameter/2.0;
      cryst_struct.Load_Crystal_System(ze_file);
      Skip_Phrases(ze_file, 3);
      particle_atoms.Atoms_Retrieval(ze_file);
	  ze_file.close();   
      }
    }