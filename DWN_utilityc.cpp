
#include "stdafx.h"

#include "DWN_utilityc.h"

//Coordinate logic (basic)---

double Get_CoorSquare_Sum(double a, double b, double c)
     {
     double coor_square = a*a + b*b + c*c;
     return coor_square;                         
     }

double Get_CoorSquare_Sum(const double* coor)
     {
     double coor_square = coor[0]*coor[0] + coor[1]*coor[1] + coor[2]*coor[2];
     return coor_square;                         
     }

double Get_Distance_FromOrigin(double a, double b, double c)
     {
     double distance = sqrt( Get_CoorSquare_Sum(a, b, c) );
     return distance;
     }  

double Get_Distance_FromOrigin(const double* coor)
	 {
	 double distance= sqrt( Get_CoorSquare_Sum(coor) );
	 return distance;
	 }

double Get_VecMag(const double* vec)
     {
     double distance = Get_Distance_FromOrigin(vec);
     return distance; 
     }
     
void Get_Vec_And_Mag(double* result_vec, double& mag, 
                     const double* start_point, const double* end_point)
     {
     SetSub_XYZ(result_vec, end_point, start_point);
     mag = Get_VecMag(result_vec);                          
     } 

double Get_Dist(const double* coor1, const double* coor2)
     {
     double distances[3];
	 SetSub_XYZ(distances, coor1, coor2);
     double dist_val = Get_VecMag(distances); 
     return dist_val;                     
     }

bool Dist_In_Bounds(const double* coor1, const double* coor2, double min_dist)  
     {
	 bool in_bounds = FALSEV;
     double distances[3];
     distances[0] = abs(coor1[0] - coor2[0]);  
     distances[1] = abs(coor1[1] - coor2[1]); 
     distances[2] = abs(coor1[2] - coor2[2]);

	 if ( (distances[0] < min_dist) && (distances[1] < min_dist) 
		&& (distances[2] < min_dist) )
      //Trying to make function as fast as possible.
	    {
	    if (Get_VecMag(distances) < min_dist)
	       {
           in_bounds = TRUEV;
	       }   
	    }
	 
	 return in_bounds;
     }
     
void Set_XYZ(double* empty_vec, const double* values_vec)
     {
     empty_vec[0] = values_vec[0];
     empty_vec[1] = values_vec[1];
     empty_vec[2] = values_vec[2];                
     }
  
void Add_XYZ(double& x, double& y, double& z, double del_x, double del_y, 
             double del_z)
     {
     x += del_x;
     y += del_y;
     z += del_z;                 
     }  

void Add_XYZ(double* vari, const double* del)
     {
     vari[0] += del[0];
     vari[1] += del[1];
     vari[2] += del[2];           
     }
     
void Sub_XYZ(double& x, double& y, double& z, double del_x, double del_y, 
             double del_z)
     {
     x -= del_x;
     y -= del_y;
     z -= del_z;                
     }  

void Sub_XYZ(double* vari, const double* del)
     {
     vari[0] -= del[0];
     vari[1] -= del[1];
     vari[2] -= del[2];            
     }

void SetAdd_XYZ(double* vari, const double* start_val, const double* add_val)
     {
     Set_XYZ(vari, start_val);
	 Add_XYZ(vari, add_val);
     }

void SetSub_XYZ(double* vari, const double* start_val, const double* sub_val)
     {
     Set_XYZ(vari, start_val);
	 Sub_XYZ(vari, sub_val);
     }

void Divide_XYZ(double& xf, double& yf, double& zf, double xi, double yi, double zi, 
                const double divider)
     {
     xf = xi/divider;   
     yf = yi/divider; 
     zf = zi/divider;                
     }
     
void Divide_XYZ(double* final_val, const double* start_val, double divider)
     {
     final_val[0] = start_val[0]/divider;  
     final_val[1] = start_val[1]/divider; 
     final_val[2] = start_val[2]/divider;        
     }

void Divide_XYZ(double* vec, double divider)
     {
	 Divide_XYZ(vec, vec, divider);
     }

void Multiply_XYZ(double* coor, double multi)
     {
	 coor[0] *= multi;
	 coor[1] *= multi;
	 coor[2] *= multi;
     }
     
void Zero_XYZ(double& x, double& y, double& z)
     {
     x = y = z = 0.0;                  
     }
     
void Zero_XYZ(double* coor)
     {
     coor[0] = coor[1] = coor[2] = 0.0;
     }

void Neg_XYZ(double* coor)
     {
     coor[0] *= -1.0;
     coor[1] *= -1.0;
     coor[2] *= -1.0;                
     }

void Abs_XYZ(double* coor)
     {
     coor[0] = abs(coor[0]);
	 coor[1] = abs(coor[1]);
	 coor[2] = abs(coor[2]);
     }

bool Same_XYZ(const double* coor1, const double* coor2)
     {
     bool same = FALSEV;
     if (Check_FP_Equality(coor1[0], coor2[0])
         && Check_FP_Equality(coor1[1], coor2[1])
         && Check_FP_Equality(coor1[2], coor2[2]))
         {
         same = TRUEV;
         }
     return same;
     }

void Max_Bound_XYZ(double* coor, double bound)
	 {
	 for (int a = 0; a < 3; ++a)
	     {
		 if (coor[a] > bound)
		    {
			coor[a] = bound;
		    }
	     }
     }

void Min_Bound_XYZ(double* coor, double bound)
	 {
	 for (int a = 0; a < 3; ++a)
	     {
		 if (coor[a] < bound)
		    {
			coor[a] = bound;
		    }
	     }
     }

void Max_Bound_XYZ(double coor, double* bounds)
	 {
	 for (int a = 0; a < 3; ++a)
	     {
		 if (coor > bounds[a])
		    {
			coor = bounds[a];
		    }
	     }
     }

void Min_Bound_XYZ(double coor, double* bounds)
	 {
	 for (int a = 0; a < 3; ++a)
	     {
		 if (coor < bounds[a])
		    {
			coor = bounds[a];
		    }
	     }
     }


bool Check_BoundariesMin(const double* coor, const double* bound)
     {
     bool in_bounds = TRUEV;
     for (int a = 0; a < 3; ++a)
         {
         if (Check_FP_LessThan(coor[a], bound[a]))
            {
            in_bounds = FALSEV;    
			break;
            }     
         }
     return in_bounds;                             
     }
     
bool Check_BoundariesMax(const double* coor, const double* bound) 
     {
     bool in_bounds = TRUEV;
     for (int a = 0; a < 3; ++a)
         {
         if (Check_FP_GreaterThan(coor[a], bound[a]))
            {
            in_bounds = FALSEV;         
            break;
            }     
         }
     return in_bounds;                             
     }

bool Check_BoundariesZeroToMax(const double* coor, const double* bound)
     {
     bool in_bounds = TRUEV;
     for (int a = 0; a < 3; ++a)
         {
         if (Check_FP_LessThan(coor[a], 0.0))
            {
            in_bounds = FALSEV;    
            break;    
            } 
         if (Check_FP_GreaterThan(coor[a], bound[a]))
            {
            in_bounds = FALSEV;  
            break;
            }         
         }
     return in_bounds;                                     
     }

bool Check_Boundaries(const double* coor, const double* min, const double* max) 
     {
     bool in_bounds = TRUEV;
     for (int a = 0; a < 3; ++a)
         {
         if (Check_FP_LessThan(coor[a], min[a]))
            {
            in_bounds = FALSEV;    
            break;    
            } 
         if (Check_FP_GreaterThan(coor[a], max[a]))
            {
            in_bounds = FALSEV;  
            break;    
            }         
         }
     return in_bounds;                             
     }

void Convert_RelCoor_To_AbsCoor(double* pos, const double* rel_coor, const double* cell_vecA,
	                            const double* cell_vecB, const double* cell_vecC)
     {
     pos[0] = cell_vecA[0]*rel_coor[0] 
            + cell_vecB[0]*rel_coor[1] 
            + cell_vecC[0]*rel_coor[2];
	  //Obtain x-axis components of absolute coordinates by summing
	  //contributions from all vectors to it.
     pos[1] = cell_vecA[1]*rel_coor[0] 
            + cell_vecB[1]*rel_coor[1] 
            + cell_vecC[1]*rel_coor[2];
	  //Total y-axis component.
     pos[2] = cell_vecA[2]*rel_coor[0] 
            + cell_vecB[2]*rel_coor[1] 
            + cell_vecC[2]*rel_coor[2];    
	  //Total z-axis component.     
     }

//Spherical coordinates---     
     
void Get_Cartesian_Coordinates(double radius, double theta, double phi, 
                               double* coors)
     {
     coors[0] = radius*sin(theta)*cos(phi);
     coors[1] = radius*sin(theta)*sin(phi);
     coors[2] = radius*cos(theta);                                
     }

void Get_Spherical_Coordinates(const double* coors, double& mag, double& theta, 
                               double& phi)
     {
     mag = Get_VecMag(coors);
     theta = acos(coors[2]/mag);
     phi = atan(coors[1]/coors[0]);                             
     }
     
//Perodic boundary coordinate logic---
     
double Get_Dist_OrthoPBC(const double* coor1, const double* coor2, int periodicity, 
						 const double* box_dimensions)
     {
     double distance_coors[3]; 
	 Get_ClosestVec_OrthoPBC(distance_coors, coor1, coor2, periodicity, box_dimensions); 
     double dist_val = Get_VecMag(distance_coors);
     return dist_val;                       
     }

void Get_ClosestVec_OrthoPBC(double* distance_coors, const double* coor1, 
						     const double* coor2, int periodicity, 
							 const double* box_dimensions)
     {
	 SetSub_XYZ(distance_coors, coor2, coor1);
	  //Get the non-PBC vector to start with.
     for (int a = 0; a < periodicity; ++a)
	  //For each dimension, translate along the box dimension
	  //vector until the shortest distance to a mirror image
	  //point is found.
         {
		 while (distance_coors[a] < box_dimensions[a]/2.0)
		  //Vector is too far to the left.
			{
			distance_coors[a] += box_dimensions[a];
			}
		 while (distance_coors[a] > box_dimensions[a]/2.0)
		  //Vector is too far to the right.
			{
			distance_coors[a] -= box_dimensions[a];
			}
         }                      
     }

void Get_PositiveVec_OrthoPBC(double* dist_coors, int periodicity, const double* box_dimensions)
	 {
	 double box_center[3];
	 Set_XYZ(box_center, box_dimensions);
	 Divide_XYZ(box_center, 2.0);
	  //Determine center of the box.
	 Get_ClosestVec_OrthoPBC(dist_coors, box_center, dist_coors, periodicity, box_dimensions);
	  //Get the shortest vector from the box center to the equivalent point.
	 Add_XYZ(dist_coors, box_center);
	  //Finally, add the box center to that vector to get the final vector.
	 }

bool Dist_In_Bounds(const double* coor1, const double* coor2, double min_dist,
	               int periodicity, const double* box_size)  
     {
	 bool in_bounds = FALSEV;
	 double distance = Get_Dist_OrthoPBC(coor1, coor2, periodicity, box_size);
	 if (distance < min_dist)
	       {
           in_bounds = TRUEV;
	       }   
	 return in_bounds;
     }

//Spatial sampling---

void Get_Spatial_Sampling_Pams(const double* lbound, const double* ubound, 
                               const double* step_pam, int dimensions, 
                               double* start_vals, double* end_vals)
     {
     double spread[3];
     int step_counts[3];
     double real_spread[3];
     double waste_space[3];
     for (int a = 0; a < dimensions; ++a)
         {
         spread[a] = ubound[a] - lbound[a];
         step_counts[a] = int( (spread[a] + FP_ERROR_FIX) / step_pam[a]);
         real_spread[a] = double(step_counts[a] - 1) * step_pam[a];
		  //One less step in logic here, due to periodicity of system.
         waste_space[a] = spread[a] - real_spread[a];
         start_vals[a] = lbound[a] + waste_space[a]/2.0;
         end_vals[a] = ( ubound[a] - waste_space[a]/2.0 ) + FP_ERROR_FIX;    
         }              
     }

void Get_Plane_Sampling(const int* ind, const double* origin, double step_size, 
                        double** sample_array, int num_samples)
     {
	 int num_zeros = 0;
	 if (ind[0] == 0) ++num_zeros;
	 if (ind[1] == 0) ++num_zeros;
	 if (ind[2] == 0) ++num_zeros;
	 int grid_factor = SC_GRID;
	 if (num_zeros == 1)
	  //e.g. 110 plane.
	    {
        grid_factor = BCC_GRID;
	    }
	 else if (num_zeros == 0)
	  //e.g. 111 plane.
	    { 
        grid_factor = FCC_GRID;
	    }
	 else if (num_zeros == 3)
	    {
		Show_Statement("A (000) plane? What is that!?"); 
		grid_factor = SC_GRID;
	    }
     Get_Plane_Sampling(ind, origin, step_size, sample_array, num_samples, grid_factor); 
	 }

void Get_Plane_Sampling(const int* ind, const double* origin, double step_size, 
                        double** sample_array, int num_samples, int grid_factor)
     {
	 //Get 2D lattice vectors for sampling---

	 double step_arrowA[2] = {step_size, 0.0};
	 double step_arrowB[2] = {0.0, step_size};
	  //Point sampling stepping values, initialized
	  //to the square grid.

	 if (grid_factor == BCC_GRID)
	  //E.g. 110 plane.
	    {
        step_arrowA[0] *= sqrt(2.0); 
	    }
	 else if (grid_factor == FCC_GRID)
	  //E.g. 111 plane.
	    { 
        step_arrowB[0] = cos(PI_CONST * 2.0/3.0) * step_size;
		step_arrowB[1] = sin(PI_CONST * 2.0/3.0) * step_size;
		}

     //Create an initial plane sampling for the (001) plane---

	 int num_samples_1D = int( sqrt(double(num_samples)) );
	  //Number of points to get along one dimension.
     int sample_count = num_samples_1D*num_samples_1D;
     for (int a = 0; a < num_samples_1D; ++a)
         {
         for (int b = 0; b < num_samples_1D; ++b)
             {
             sample_array[sample_count][0] = double(a)*step_arrowA[0] + double(b)*step_arrowB[0];
             sample_array[sample_count][1] = double(a)*step_arrowA[1] + double(b)*step_arrowB[1];
             sample_array[sample_count][2] = 0.0;
             }     
         }

	 //Determine the average/middle location of all data points---

     double middle_point[3];
      //Middle point of all the above sampling points.
	 double middle_step = double(num_samples_1D - 1)/2.0;
     middle_point[0] = middle_step*step_arrowA[0] + middle_step*step_arrowB[0];
     middle_point[1] = middle_step*step_arrowA[1] + middle_step*step_arrowB[1];
     middle_point[2] = 0.0; 

	 //Rotate all the points such that they collectively define the
	 //requested plane---

	 int start_ind[3] = { 0, 0, 1 };
	 double rot_matrix[9];
     if ( (ind[0] != 0) || (ind[1] != 0) || (ind[2] != 1) )
      //Only rotate if the desired plane is not the (001) plane.
	    { 
        Calc_Cubic_Plane_Rotation_Matrix(start_ind, ind, rot_matrix);   
        for (int a = 0; a < sample_count; ++a)
         //Do the rotation into the desired plane.
           {
           Rotate_Vector(rot_matrix, sample_array[a]); 
           }
        Rotate_Vector(rot_matrix, middle_point);
        }

	 //Rotate the point sampling within the desired plane itself 
	 //to the desired alignment---

	 double row_pointer[3];
	 SetSub_XYZ(row_pointer, sample_array[1], sample_array[0]);
	  //Example point-->point vector that will be used in alignment.
	 double desired_direction[3];
	 double final_direction[3];
	 double rotation_vector[3];
	 double rotation_angle;
     if (grid_factor == SC_GRID)
	 //No alignment necessary.
	    {
        rotation_angle = 0.0;
		Set_XYZ(rotation_vector, origin);
	    }
	 else if (grid_factor == BCC_GRID)
	    {
        desired_direction[0] = desired_direction[1] = desired_direction[2] = 1.0;
		Determine_Closest_Vector(row_pointer, desired_direction, origin, final_direction);
		 //Determine the closest vector in the [111] vector family to the current
		 //point-->point vector that is in the (110) plane.
	    Calculate_Rotation_Parameters(row_pointer, final_direction, 
			                          rotation_vector, rotation_angle);
	     //Determine the rotation matrix required to transform the point-->point vector
		 //to the one desired.
	    }
	 else if (grid_factor == FCC_GRID)
	    {
		desired_direction[0] = desired_direction[1] = 1.0;
		desired_direction[2] = 0.0;
	    Determine_Closest_Vector(row_pointer, desired_direction, origin, final_direction);
	    Calculate_Rotation_Parameters(row_pointer, final_direction, 
		                              rotation_vector, rotation_angle);

		//The following statements correct for shortcomings in the rotation determination.
		//e.g. a faceted nanoparticle will get a symmetrical coating design.

		if (ind[2] < 0)
		   {
		   rotation_angle += 1.0/3.0*PI_CONST;
		   }

		if ((ind[1] + ind[0]) == 0 )
		   {   
		   rotation_angle += 4.0/3.0*PI_CONST;
		   }
	    else
		   {
           rotation_angle -= 4.0/3.0*PI_CONST;
		   }
		}

	 Calc_Rotation_Matrix(rotation_vector, rotation_angle, rot_matrix);
	 for (int a = 0; a < sample_count; ++a)
         //Do the rotation within the plane to get symmetric alignment.
           {
           Rotate_Vector(rot_matrix, sample_array[a]); 
           }
	 Rotate_Vector(rot_matrix, middle_point);

	 //Translate the middle point of the point samples such that it is
	 //located at the requested plane origin---

	 double shift_coors[3];
	 SetSub_XYZ(shift_coors, origin, middle_point);
     for (int a = 0; a < sample_count; ++a)
      //Shift the middle point of the data to the plane origin.
         {
         Add_XYZ(sample_array[a], shift_coors);       	 
		 }
	 }

void Get_OrthoRhombic_Sampling(const double* A_vec, const double* B_vec, const double* C_vec,  
	                           const double* sizes, double** point_array, int max_points, int& point_index)
     {

	 //Get bounds relative to (0, 0, 0) origin---

	 double max_bounds[3];
	 Set_XYZ(max_bounds, sizes);
	 Divide_XYZ(max_bounds, 2.0);

	 //Prepare lattice vector information---
	 
	 double lattice_vecs[3][3];
	 Set_XYZ(lattice_vecs[0], A_vec);
	 Set_XYZ(lattice_vecs[1], B_vec);
	 Set_XYZ(lattice_vecs[2], C_vec);

	 //Find all the points in the sampling---

	 const double MIN_DIST = 0.001;
	  //Minimum distance used to check for a point lying
	  //at the same location as another point.

	 int current_test_pos = 0;
	 point_index = 0;
	 double test_pos[3];
	 double try_pos[3];
	 double abs_try_pos[3];
	 Zero_XYZ(point_array[0]);
	  //Start point propagation at the origin.
	 while (current_test_pos <= point_index)
	  //Finished condition = No more points being added that are to be branched from.
	   {
	   if (point_index >= (max_points - 6))
		//Memory check.
	      {
		  Show_Warning("NOT ENOUGH MEMORY FOR ORTHORHOMBIC POINT SAMPLING!");
		  break;
	      }

	   Set_XYZ(test_pos, point_array[current_test_pos]);

	   //Try moving in all six directions (positive or negative
	   //per lattice vector)---

	   for (int a = 0; a < 3; ++a)
	     {
		 for (int b = 0; b < 2; ++b)
		     {
			 if (b == 0)
			    {
			    SetAdd_XYZ(try_pos, test_pos, lattice_vecs[a]);
			    }
			 else
			    {
                SetSub_XYZ(try_pos, test_pos, lattice_vecs[a]);  
			    }
			 Set_XYZ(abs_try_pos, try_pos);
			 Abs_XYZ(abs_try_pos);
			  //Compare absolute coordinates to the maximum boundary coordinates.
			 if (No_Spatial_Prob(const_cast<const double**>(point_array), point_index + 1, try_pos, MIN_DIST)
			    && Check_BoundariesMax(abs_try_pos, max_bounds))
				{
				++point_index;
				Set_XYZ(point_array[point_index], try_pos);
				}
			 }
	     }

	   ++current_test_pos;
	   }

	 if (COOR_LIST_ORDERING)
	   {
	   Order_Coordinate_List(point_array, point_index);
	   }
     }

void Get_Spherical_Sampling(const double* A_vec, const double* B_vec, const double* C_vec,  
	                        double radius, double** point_array, int max_points, int& point_index)
     {

	 //Prepare lattice vector information---

	 double lattice_vecs[3][3];
	 Set_XYZ(lattice_vecs[0], A_vec);
	 Set_XYZ(lattice_vecs[1], B_vec);
	 Set_XYZ(lattice_vecs[2], C_vec);

	 //Find spherical point sampling---

	 const double MIN_DIST = 0.001;
	  //Minimum distance between placed points, used to check
	  //for points being placed atop each other.

	 int current_test_pos = 0;
	 point_index = 0;
	 double test_pos[3];
	 double try_pos[3];
	 double dist_temp;
	 Zero_XYZ(point_array[0]);
	  //Start point sampling at the origin.
	 while (current_test_pos <= point_index)
	  //Finished condition = No more points being added that are to be branched from.
	   {
	   if (point_index >= (max_points - 6))
		//Memory check.
	      {
		  Show_Warning("NOT ENOUGH MEMORY FOR SPHERICAL POINT SAMPLING!");
		  break;
	      }

	   Set_XYZ(test_pos, point_array[current_test_pos]);

	   //Look for points along each lattice vector
	   //in the positive and negative direction---

	   for (int a = 0; a < 3; ++a)
	     {
	     for (int b = 0; b < 2; ++b)
		   {
		   if (b == 0) 
		       {
			   SetAdd_XYZ(try_pos, test_pos, lattice_vecs[a]);
		       }
		   else 
		       {
			   SetSub_XYZ(try_pos, test_pos, lattice_vecs[a]);
		       }
		   dist_temp = Get_VecMag(try_pos);
	       if (No_Spatial_Prob(const_cast<const double**>(point_array), point_index + 1, try_pos, MIN_DIST)
			   && (Check_FP_LessThan(dist_temp, radius) ) )
		       {
		       ++point_index;
		       Set_XYZ(point_array[point_index], try_pos);
		       }
		   }
	     }
	   ++current_test_pos;
	   }

	 if (COOR_LIST_ORDERING)
	   {
	   Order_Coordinate_List(point_array, point_index);
	   }
	 }

void Expand_Mesh_Grid(double mesh_pam, double vari_pam, long int max_tries, 
	                  double** mesh_grid, double box_size, long int start_index, 
					  long int& point_index)
   //Recursive function used by random mesh generator.
     {
	 double min_pam = mesh_pam * (1.0 - vari_pam);
	 double max_pam = mesh_pam * (1.0 + vari_pam); 
     double test_pos[3];
     double move_coors[3];
     double final_move[3];
	 double point_dist;
     double max_coor = box_size/2.0;
	 double box_min[3] = { -1.0*max_coor, -1.0*max_coor, -1.0*max_coor };
	 double box_max[3] = { max_coor, max_coor, max_coor };
	 double min_dist = min_pam - FP_ERROR_FIX;

	 for (int a = 0; a < max_tries; ++a)
	    {
        Set_XYZ(test_pos, mesh_grid[start_index]);
		move_coors[0] = GetRandDec(-1.0, 1.0);
        move_coors[1] = GetRandDec(-1.0, 1.0);
        move_coors[2] = GetRandDec(-1.0, 1.0);
		point_dist = mesh_pam;
		if (!Check_FP_Equality(min_pam, max_pam))
         //If variance in point separation is not 0....
			{
			point_dist = GetRandDec(min_pam, max_pam);
			}
	    Get_NewMag_Vec(move_coors, final_move, point_dist);
        Add_XYZ(test_pos, final_move);
		if (No_Spatial_Prob(const_cast<const double**>(mesh_grid), point_index, test_pos, min_dist)
			&& Check_Boundaries(test_pos, box_min, box_max) )
		  {
          Set_XYZ(mesh_grid[point_index], test_pos);
          ++point_index; 
		  Expand_Mesh_Grid(mesh_pam, vari_pam, max_tries, mesh_grid, box_size, 
			               point_index - 1, point_index);
		  }
	    }
     }


void Get_Random_OMesh_Sampling(double mesh_pam, double vari_pam, long int max_tries,
	                          const double* sizes, double** mesh_grid, int max_points, 
							  long int& point_index)
     {
	 double max_bounds[3];
	  //Orthorhombic boundaries.
	 Set_XYZ(max_bounds, sizes);
	 Divide_XYZ(max_bounds, 2.0);

	 //Determine mesh grid generation parameters---
		 
	 double min_pam = mesh_pam * (1.0 - vari_pam);
	 double max_pam = mesh_pam * (1.0 + vari_pam); 
	  //Minimum and maximum bond lengths used in creating
	  //the grid.
	 bool exact_pam = Check_FP_Equality(min_pam, max_pam);
	 double min_dist = min_pam - FP_ERROR_FIX;
	  //Minimum distance allowed between points in the sampling.

	 //Get the point sampling---

	 double move_coors[3];
	 double final_move[3];
	 double point_dist;

     int current_test_pos = 0;
	 point_index = 0;
	 double test_pos[3];
	 double try_pos[3];
	 double abs_try_pos[3];
	 Zero_XYZ(mesh_grid[0]);
	  //Start sampling at the origin.
	 while (current_test_pos <= point_index)
	  //Finished condition = No more points being added that are to be branched from.
	   {
	   if (point_index >= (max_points - max_tries))
		//Memory check.
	      {
		  Show_Warning("NOT ENOUGH MEMORY FOR RANDOM ORTHRORHOMBIC POINT SAMPLING!");
		  break;
	      }

	   Set_XYZ(test_pos, mesh_grid[current_test_pos]);
	   for (int a = 0; a < max_tries; ++a)
	      {
	      move_coors[0] = GetRandDec(-1.0, 1.0);
          move_coors[1] = GetRandDec(-1.0, 1.0);
          move_coors[2] = GetRandDec(-1.0, 1.0);
	      point_dist = mesh_pam;
	      if (!exact_pam)
           //If variance in point separation is not 0....
			 {
			 point_dist = GetRandDec(min_pam, max_pam);
		     }
	      Get_NewMag_Vec(move_coors, final_move, point_dist);
		  SetAdd_XYZ(try_pos, test_pos, final_move);
	      Set_XYZ(abs_try_pos, try_pos);
	      Abs_XYZ(abs_try_pos);
	      if (No_Spatial_Prob(const_cast<const double**>(mesh_grid), point_index + 1, try_pos, min_dist)
			  && Check_BoundariesMax(abs_try_pos, max_bounds))
		     {
		     ++point_index;
		     Set_XYZ(mesh_grid[point_index], try_pos);
		     }
	      }
	   ++current_test_pos;
	   }
 
	 if (COOR_LIST_ORDERING)
	   {
	   Order_Coordinate_List(mesh_grid, point_index);
	   }
     }

void Get_Random_SMesh_Sampling(double mesh_pam, double vari_pam, long int max_tries,
	                          double radius, double** mesh_grid, int max_points, 
							  long int& point_index)
     { 
	 //Get random mesh grid generation parameters---

	 double min_pam = mesh_pam * (1.0 - vari_pam);
	 double max_pam = mesh_pam * (1.0 + vari_pam);
	 bool equal_pam = Check_FP_Equality(min_pam, max_pam);
	 double min_dist = min_pam - FP_ERROR_FIX;

	 //Get point sampling---

	 double move_coors[3];
	 double final_move[3];
	 double point_dist;

     int current_test_pos = 0;
	 point_index = 0;
	 double test_pos[3];
	 double try_pos[3];
	 double dist_temp;
	 Zero_XYZ(mesh_grid[0]);
	 while (current_test_pos <= point_index)
	  //Finished condition = No more points being added that are to be branched from.
	   {
	   if (point_index >= (max_points - max_tries))
		//Memory check.
	      {
		  Show_Warning("NOT ENOUGH MEMORY FOR RANDOM SPHERICAL POINT SAMPLING!");
		  break;
	      }

	   Set_XYZ(test_pos, mesh_grid[current_test_pos]);
	   for (int a = 0; a < max_tries; ++a)
	      {
	      move_coors[0] = GetRandDec(-1.0, 1.0);
          move_coors[1] = GetRandDec(-1.0, 1.0);
          move_coors[2] = GetRandDec(-1.0, 1.0);
	      point_dist = mesh_pam;
	      if (!equal_pam)
			 {
			 point_dist = GetRandDec(min_pam, max_pam);
		     }
	      Get_NewMag_Vec(move_coors, final_move, point_dist);
		  SetAdd_XYZ(try_pos, test_pos, final_move);
		  dist_temp = Get_VecMag(try_pos);
	      if (No_Spatial_Prob(const_cast<const double**>(mesh_grid), point_index + 1, try_pos, min_dist)
			  && Check_FP_LessThan(dist_temp, radius))
		     {
		     ++point_index;
		     Set_XYZ(mesh_grid[point_index], try_pos);
		     }
	      }
	   ++current_test_pos;
	   }
 
	 if (COOR_LIST_ORDERING)
	   {
	   Order_Coordinate_List(mesh_grid, point_index);
	   }
     }

void Determine_Closest_Vector(const double* test_vec, const double* vec_family, double* result_vec)
     {
	 double zero_plane[3] = { 0.0, 0.0, 0.0 };
	 Determine_Closest_Vector(test_vec, vec_family, zero_plane, result_vec); 	
     }

void Determine_Closest_Vector(const double* test_vec, const double* vec_family, 
	                          const double* plane_normal, double* result_vec)
     {
     Determine_Closest_Vector(test_vec, vec_family, plane_normal, result_vec, 0);
     }

void Determine_Second_Closest_Vector(const double* test_vec, const double* vec_family, 
	                                 const double* plane_normal, double* result_vec)
     {
     Determine_Closest_Vector(test_vec, vec_family, plane_normal, result_vec, 1);
     }

void Determine_Closest_Vector(const double* test_vec, const double* vec_family, 
	                          const double* plane_normal, double* result_vec,
							  int closeness_factor)
     {
	 double smallest_angle, second_smallest_angle, temp_angle;
     smallest_angle = 2.1*PI_CONST;
	 second_smallest_angle = 2.1*PI_CONST;
	  //Initialize to large angles > 360 degrees.

	 Zero_XYZ(result_vec);
	 double second_result_vec[3];
	  //Second closest vector determined.
	 Zero_XYZ(second_result_vec);

	 double vec_fam_instance[3];
	 int index1, index2;
	 double plane_term;
	 double zero_origin[3] = { 0.0, 0.0, 0.0 };
	 for (int a = 0; a < 3; ++a)
      //Permutate the 48 vectors to consider and find out which one is closest. 
	  //This for loop considers: (abc), (bac), (cba) permutation. Within the loop,
	  //permutation of the type (abc), (acb) is performed. 
	    {

		//Determine the second and third indices, as the value of 'a' sets
	    //the first index (a, b, or c)---

		index1 = a + 1;
		if (index1 > 2)
		  {
          index1 -= 3; 
		  }
		index2 = a + 2;
		if (index2 > 2)
		  {
          index2 -= 3; 
		  }
	    for (int b = 0; b < 2; ++b)
		 //Permutation of type: (abc), (acb).
		  {
	      for (int c = 0; c < 4; ++c)
		   //Permutation of the type: (abc), (-abc), (a-bc), (ab-c).    
		     {
             for (int d = 0; d < 2; ++d)
			  //Permutation of the type: (abc), (-a-b-c).
			   {
			   vec_fam_instance[0] = vec_family[a];
			   if (b == 0)
			     {
		         vec_fam_instance[1] = vec_family[index1];
		         vec_fam_instance[2] = vec_family[index2];
			     }
			   else
			     {
			     vec_fam_instance[1] = vec_family[index2];
		         vec_fam_instance[2] = vec_family[index1];
			     }
		       if (c > 0)
		         {
                 vec_fam_instance[c - 1] *= -1.0;
		         }
		       if (d == 1)
			     {
                 Neg_XYZ(vec_fam_instance);
			     }

			   plane_term = Solve_Plane_Equation(plane_normal, vec_fam_instance, zero_origin);
		       if (Check_FP_Equality(plane_term, 0.0))
			    //Only consider vectors in the plane of interest.
		         { 
		         Get_Angle(test_vec, vec_fam_instance, temp_angle);
		         if (Check_FP_LessThan(temp_angle, smallest_angle))
		           {
				   second_smallest_angle = smallest_angle;
				   smallest_angle = temp_angle;
				   Set_XYZ(second_result_vec, result_vec);
		           Set_XYZ(result_vec, vec_fam_instance);
		           }		
				 else if (Check_FP_LessThan(temp_angle, second_smallest_angle) 
					      && Check_FP_GreaterThan(temp_angle, smallest_angle))
				   {
				   second_smallest_angle = temp_angle;
		           Set_XYZ(second_result_vec, vec_fam_instance);
				   }
		         } 
			   }
		     }
		  }
		}

	 if (closeness_factor == 1)
	   {
       Set_XYZ(result_vec, second_result_vec);
	   }

     }

double Solve_Plane_Equation(const int* indices, const double* point, const double* origin)
 //h*(x-x0) + k*(y-y0) + l(z-z0) = Return value
     {
     double pe_terms[3];
     pe_terms[0] = double(indices[0]) * ( point[0] - origin[0] );
     pe_terms[1] = double(indices[1]) * ( point[1] - origin[1] );
     pe_terms[2] = double(indices[2]) * ( point[2] - origin[2] );
     double term_sum = pe_terms[0] + pe_terms[1] + pe_terms[2];                             
     return term_sum;
     }
 
double Solve_Plane_Equation(const double* normal, const double* point, const double* origin)
     {
     double pe_terms[3];
     pe_terms[0] = normal[0] * ( point[0] - origin[0] );
     pe_terms[1] = normal[1] * ( point[1] - origin[1] );
     pe_terms[2] = normal[2] * ( point[2] - origin[2] );
     double term_sum = pe_terms[0] + pe_terms[1] + pe_terms[2];                             
     return term_sum;                                                 
     }

void Order_Coordinate_List(double** first_list, int num_points)
     { 
	 double** second_list; 	  
     Get_Memory(second_list, num_points, 3);
	 bool* ordered_point = new bool[num_points];

	 for (int a = 0; a < num_points; ++a)
       //Copy over the first list into the new list.
	    {
        ordered_point[a] = FALSEV;  
		Set_XYZ(second_list[a], first_list[a]);
	    }
	 ordered_point[0] = TRUEV;
	  //Keep the first point the same in the new list.

	 double dist_temp, lowest_distance;
	 int lowest_dist_index;
	 const double INIT_DISTANCE = Get_Dist(first_list[0], first_list[1]) * 100000.0;
	 for (int a = 1; a < num_points; ++a)
	  //Fill each point in the ordered array by comparing with the (a - 1) point.	  
	    {
        lowest_distance = INIT_DISTANCE;
		for (int b = 1; b < num_points; ++b)
         //Cycle through list of points.
		    {
            if (!ordered_point[b])
			  {
              dist_temp = Get_Dist(first_list[a - 1], second_list[b]);
			  if (dist_temp < lowest_distance)
			     {
                 lowest_distance = dist_temp;
				 lowest_dist_index = b;
			     }
			  }
		    }
		 Set_XYZ(first_list[a], second_list[lowest_dist_index]);
	     ordered_point[lowest_dist_index] = TRUEV;
	     }
	   
	 delete[] ordered_point;
	 Free_Memory(second_list, num_points);
     }

bool No_Spatial_Prob(const double** points, int num_points, 
	                 const double* test_point, double min_dist)
     {
     bool no_prob = TRUEV;
	 for (int a = 0; a < num_points; ++a)
	    {
        if (Dist_In_Bounds(points[a], test_point, min_dist) )
		  {
          no_prob = FALSEV;
		  break;
		  }
	    }

	 return no_prob;
     }

void Get_Coord_Numbers(const double** points, int num_points, 
	                   int* coor_nums, double min_dist, double max_dist)
	 {
	 Zero_Array(coor_nums, num_points);
	 double dist_temp;
	 bool greater_test, lesser_test;
	 for (int a = 0; a < num_points; ++a)
	     {		 
		 for (int b = (a + 1); b < num_points; ++b)
		     {
			 dist_temp = Get_Dist(points[a], points[b]);
			 greater_test = Check_FP_GreaterOrEqual(dist_temp, min_dist);
			 lesser_test = Check_FP_LessOrEqual(dist_temp, max_dist);
			 if (greater_test && lesser_test)
			    {
				++coor_nums[a];
				++coor_nums[b];
			    }
		     }
	     }
     }

//Rotation logic---

void Take_Cross_Product(const double* vec1, const double* vec2, double* crossP)
     {
     crossP[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
     crossP[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
     crossP[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];                         
     } 

void Take_Dot_Product(const double* vec1, const double* vec2, double& dot_prod)
     {
     double dotP[3];
     dotP[0] = vec1[0]*vec2[0];
     dotP[1] = vec1[1]*vec2[1];
     dotP[2] = vec1[2]*vec2[2];        
     dot_prod = dotP[0] + dotP[1] + dotP[2];                 
     } 

void Get_Angle(const double* vec1, const double* vec2, double& angle)
     {
     double dot_prod;
	 Take_Dot_Product(vec1, vec2, dot_prod);
	 double vec2_mag = Get_VecMag(vec2); 
     double vec1_mag = Get_VecMag(vec1);
     double mag_ratio = dot_prod/(vec2_mag*vec1_mag);
     angle = acos(mag_ratio); 
     }

void Get_Unit_Vec(const double* vec, double* unit_vec)
     {
     double current_mag = Get_VecMag(vec);
	 Divide_XYZ(unit_vec, vec, current_mag);               
     }

void Get_NewMag_Vec(const double* vec, double* new_mag_vec, double new_mag)
     {
     double current_mag = Get_VecMag(vec);
     double scale_factor = new_mag/current_mag;
	 Set_XYZ(new_mag_vec, vec);
	 Multiply_XYZ(new_mag_vec, scale_factor);                   
     } 

void Lengthen_Vec(double* vec, double del_mag)
     {
     double current_mag = Get_VecMag(vec);
     double new_mag = current_mag + del_mag;
     double scale_factor = new_mag/current_mag; 
     vec[0] *= scale_factor;   
     vec[1] *= scale_factor;      
     vec[2] *= scale_factor;                     
     }
     
void Convert_RotAngles_To_RotVecs(double theta, double phi, double* start, double* end)
     {
     start[0] = 0.0;
     start[1] = 0.0;
     start[2] = 1.0;
	  //Set start vector to (0, 0, 1), which has theta and phi
	  //of zero in spherical coordinates.
     Get_Cartesian_Coordinates(1.0, theta, phi, end);     
	  //Translate angles into Cartesian location.
     }

void Calculate_Rotation_Parameters(const double* source_vec, const double* dest_vec, 
                                   double* rot_vec, double& rot_angle)
     {
	 double source_mag = Get_VecMag(source_vec);
	  //Get magnitude of source vector.
	 double dest_mag = Get_VecMag(dest_vec);
	  //Do the same for end vector.
	 double dest_vec_temp[3];
	 Set_XYZ(dest_vec_temp, dest_vec);
	 if (!Check_FP_Equality(source_mag, dest_mag))
	  //Source and destination vectors must be of equal magnitude.
	  //If they aren't, rescale the destination vector.
	    {
        Get_NewMag_Vec(dest_vec, dest_vec_temp, source_mag);
	    }
     
	  if (Same_XYZ(source_vec, dest_vec_temp))
      //No rotation is needed if these two vectors are the same.
        {
        rot_vec[0] = 0.0;
        rot_vec[1] = 0.0;
        rot_vec[2] = 0.0;
        rot_angle = 0.0;
        return;
        }
	 
	 Take_Cross_Product(source_vec, dest_vec_temp, rot_vec); 
      //Get the vector (axis) to use in the rotation of the old vector into the 
      //new vector location.
	 Get_Angle(source_vec, dest_vec_temp, rot_angle);
      //Get the rotation angle (via the dot product).

	 if (Check_FP_Equality(rot_angle, PI_CONST))
	  //Special case: 180 degrees rotation, where the cross product/rotation vector 
	  //is not well defined.
	    { 
		double temp_vec[3] = { 0.0, source_mag, 0.0 };
		 //Make up a vector.
		
		double temp_rot_angle;
		Get_Angle(source_vec, temp_vec, temp_rot_angle);
		if ( (Check_FP_Equality(temp_rot_angle, PI_CONST)) 
		  || (Check_FP_Equality(temp_rot_angle, 0.0)) )
         //Just for safety, make sure the made up vector isn't equivalent in direction
		 //to the source vector. Yeah...I care that much! 
		    {
			temp_vec[0] = source_mag;
			temp_vec[1] = 0.0;
		    }
		Take_Cross_Product(source_vec, temp_vec, rot_vec);
	     //Get a vector that is orthogonal to the source vector for use in rotation.
	    }                                 
     }   

void Calc_Rotation_Matrix(double phi, double theta, double psi, double* rot_matrix)
     {
	 double cos_phi = cos(phi);
	 double sin_phi = sin(phi);
	  //X rotation.

	 double cos_theta = cos(theta);
	 double sin_theta = sin(theta);
	  //Y rotation.

	 double cos_psi = cos(psi);
	 double sin_psi = sin(psi);
	  //Z rotation.

	 //Here comes the rotation matrix---

	 rot_matrix[0] = cos_theta*cos_psi;
	 rot_matrix[1] = -1.0*cos_phi*sin_psi + sin_phi*sin_theta*cos_psi;
	 rot_matrix[2] = sin_phi*sin_psi + cos_phi*sin_theta*cos_psi;
	 rot_matrix[3] = cos_theta*sin_psi;
	 rot_matrix[4] = cos_phi*cos_psi + sin_phi*sin_theta*sin_psi;
	 rot_matrix[5] = -1.0*sin_phi*cos_psi + cos_phi*sin_theta*sin_psi;
	 rot_matrix[6] = -1.0*sin_theta;
	 rot_matrix[7] = sin_phi*cos_theta;
	 rot_matrix[8] = cos_phi*cos_theta;
     }

    
void Calc_Rotation_Matrix(const double* rotation_vector, double angle, 
                          double* rot_matrix)
     {
     if (Check_FP_Equality(angle, 0.0))
      //If no rotation is needed, set the rotation matrix to the unity matrix.
        {
        for (int a = 0; a < 9; ++a)
            {
            if ( (a % 4) == 0)
               {
               rot_matrix[a] = 1.0;
               }
            else
               {
               rot_matrix[a] = 0.0;           
               }     
            }
		return;
        }

     double rot_unit_vec[3];
      //Rotation vector expressed as a unit vector
     Get_Unit_Vec(rotation_vector, rot_unit_vec);
        
     double cos_angle = cos(angle);  
     double sin_angle = sin(angle);
     double one_minus_cos = 1 - cos_angle;
     double ux = rot_unit_vec[0];
     double uy = rot_unit_vec[1];
     double uz = rot_unit_vec[2];
     double uxuy = ux*uy;
     double uxuz = ux*uz;
     double uyuz = uy*uz;
     double ux2 = ux*ux;
     double uy2 = uy*uy;
     double uz2 = uz*uz; 
       
     //Here comes the rotation matrix determination---
     
     rot_matrix[0] = cos_angle + ux2*one_minus_cos; 
     rot_matrix[1] = uxuy*one_minus_cos - uz*sin_angle; 
     rot_matrix[2] = uxuz*one_minus_cos + uy*sin_angle; 
     rot_matrix[3] = uxuy*one_minus_cos + uz*sin_angle; 
     rot_matrix[4] = cos_angle + uy2*one_minus_cos; 
     rot_matrix[5] = uyuz*one_minus_cos - ux*sin_angle; 
     rot_matrix[6] = uxuz*one_minus_cos - uy*sin_angle; 
     rot_matrix[7] = uyuz*one_minus_cos + ux*sin_angle; 
     rot_matrix[8] = cos_angle + uz2*one_minus_cos;   
     
      //Ref: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/                          
     }

void Calc_Rotation_Matrix(const double* start_vec, const double* end_vec, double* rotation_matrix)
     {
     double rotation_vec[3];
     double angle;   
     Calculate_Rotation_Parameters(start_vec, end_vec, rotation_vec, angle);
     Calc_Rotation_Matrix(rotation_vec, angle, rotation_matrix);
     }

void Calc_Cubic_Plane_Rotation_Matrix(const int* start_indices, 
                                      const int* end_indices, double* rot_matrix)
     {
     double old_normal[3];
     double new_normal[3];

	 //Get normal vectors that correspond to the Miller indices---

     for (int a = 0; a < 3; ++a)
         {
         old_normal[a] = 0.0;
         if (start_indices[a] != 0)
            {
            old_normal[a] = 1.0/double(start_indices[a]);                  
            }     
         new_normal[a] = 0.0;
         if (end_indices[a] != 0)
            {
            new_normal[a] = 1.0/double(end_indices[a]);                  
            }
         } 

	 //Calculate the rotation matrix---

     double rot_vec[3];
     double rot_angle;
     Calculate_Rotation_Parameters(old_normal, new_normal, rot_vec, rot_angle);
     Calc_Rotation_Matrix(rot_vec, rot_angle, rot_matrix);                                
     }

void Rotate_Vector(const double* rot_mat, double* vec)
     {
     double x_comp = rot_mat[0]*vec[0] + rot_mat[1]*vec[1] + rot_mat[2]*vec[2];  
     double y_comp = rot_mat[3]*vec[0] + rot_mat[4]*vec[1] + rot_mat[5]*vec[2];  
     double z_comp = rot_mat[6]*vec[0] + rot_mat[7]*vec[1] + rot_mat[8]*vec[2]; 
     vec[0] = x_comp;
     vec[1] = y_comp;
     vec[2] = z_comp;  
     }
     
void Rotate_Vector(const double* rot_mat, double* vec, const double* origin)
     {
	 Sub_XYZ(vec, origin);
     double x_comp = rot_mat[0]*vec[0] + rot_mat[1]*vec[1] + rot_mat[2]*vec[2];  
     double y_comp = rot_mat[3]*vec[0] + rot_mat[4]*vec[1] + rot_mat[5]*vec[2];  
     double z_comp = rot_mat[6]*vec[0] + rot_mat[7]*vec[1] + rot_mat[8]*vec[2]; 
     vec[0] = x_comp;
     vec[1] = y_comp;
     vec[2] = z_comp;
	 Add_XYZ(vec, origin);
     }

 