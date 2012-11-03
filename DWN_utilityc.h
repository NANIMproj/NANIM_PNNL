
#ifndef SPACEY_PING
#define SPACEY_PING

#include <iostream>
#include <cmath>

using namespace std;

#include "DWN_utilityconstants.h"
#include "DWN_utilitym.h"

//Constant enumeration---

const bool PLACE_MODE = 0;
 //Absolute coordinate assignment.
const bool SHIFT_MODE = 1;
 //Relative (translational) coordinate assignment.

const bool CLOCKWISE = 1;
 //Clockwise rotation.
const bool COUNTERCLOCK = 0;
 //Counterclockwise rotation.

//Coordinate logic basics---

double Get_CoorSquare_Sum(double, double, double);
 //Obtain the sum of the squares of the passed values.
double Get_CoorSquare_Sum(const double*);
 //Obtains the sum of the squares of the passed coordinates.

double Get_Distance_FromOrigin(double, double, double);
 //Returns the distance of the passed point from the
 //(0, 0, 0) location.
double Get_Distance_FromOrigin(const double*);
 //Returns the distance of the passed point from the
 //(0, 0, 0) location.
double Get_VecMag(const double*);
 //Returns the magnitude of the vector.

void Get_Vec_And_Mag(double*, double&, const double*, const double*);
 //Returns the vector connecting two points and its magnitude.
double Get_Dist(const double*, const double*);
 //Gets the distance between two points in 3-D space.
bool Dist_In_Bounds(const double*, const double*, double);
 //Returns TRUEV if the two points in space are closer than the passed value.
 
void Set_XYZ(double*, const double*);
 //Sets the passed coordinate set to the values of another set.
void Add_XYZ(double&, double&, double&, double, double, double);
 //Adds with X-Y-Z coordinate sum logic (x_new = x_old + del_x, etc.)
void Add_XYZ(double*, const double*);
 //Adds with X-Y-Z coordinate sum logic (x_new = x_old + del_x, etc).
void Sub_XYZ(double&, double&, double&, double, double, double);
 //Subtracts with X-Y-Z coordinate logic.
void Sub_XYZ(double*, const double*);
 //Subtracts with X-Y-Z coordinate logic.
void SetAdd_XYZ(double*, const double*, const double*);
 //Sets the passed vector to the sum of two vectors.
void SetSub_XYZ(double*, const double*, const double*);
 //Sets the passed vector to the difference between two vectors.
void Divide_XYZ(double&, double&, double&, double, double, double, const double);
 //Divides the passed coordinate set by the passed parameter and stores the
 //result in another coordinate set.
void Divide_XYZ(double*, const double*, double);
 //Divides the passed coordinate set by the passed parameter and stores the
 //result in another coordinate set.
void Divide_XYZ(double*, double);
 //Divides the coordinate set by the passed parameter.
void Multiply_XYZ(double*, double);
 //Multiplies the coordinate set by the passed parameter.
void Zero_XYZ(double&, double&, double&);
 //Sets three variables to zero.
void Zero_XYZ(double*);
 //Sets a coordinate set to zero.
void Neg_XYZ(double*);
 //Multiples all values in the coordinate set by -1.
void Abs_XYZ(double*);
 //Sets all values in the coordinate set to their absolute values.
bool Same_XYZ(const double*, const double*);
 //Returns true if the two coordinate sets are equal.

void Max_Bound_XYZ(double*, double);
 //Enforces a maximum value on each of three coordinates.
void Min_Bound_XYZ(double*, double);
 //Enforces a minimum value on each of three coordinates.
void Max_Bound_XYZ(double, double*);
 //Enforces a set of three maximum values on a coordinate.
void Min_Bound_XYZ(double, double*);
 //Enforces a set of three minimum values on a coordinate.

bool Check_BoundariesMin(const double*, const double*); 
 //Returns TRUEV if coordinates are within certain minimums.
bool Check_BoundariesMax(const double*, const double*); 
 //Returns TRUEV if coordinates are within certain maximums
bool Check_BoundariesZeroToMax(const double*, const double*); 
 //Returns TRUEV if coordinates are within an origin-minimum
 //and certain maximums.
bool Check_Boundaries(const double*, const double*, const double*);
 //Returns TRUEV if coordinates are within certain minimums
 //and maximums.

void Convert_RelCoor_To_AbsCoor(double*, const double*, const double*, 
	                            const double*, const double*);
 //Converts a set of relative coordinates to absolute coordinates via a set of three
 //lattice vectors.

//Spherical coordinate logic---

void Get_Cartesian_Coordinates(double, double, double, double*); 
 //Takes spherical coordinates and returns the corresponding Cartesian coordinates.
void Get_Spherical_Coordinates(const double*, double&, double&, double&);
 //Takes Cartesian coordinates and returns the corresponding spherical coordinates.

//Perodic boundary coordinate logic---

double Get_Dist_OrthoPBC(const double*, const double*, int, const double*);
 //Get the distance between two coordinate sets.
void Get_ClosestVec_OrthoPBC(double*, const double*, const double*, int, const double*);
 //Stores the shortest connecting vector between two points.
void Get_PositiveVec_OrthoPBC(double*, int, const double*);
 //Returns the vector that connects the box origin to a point in the box that 
 //is an equivalent position to the passed point.
bool Dist_In_Bounds(const double*, const double*, double, int, const double*); 
 //Returns TRUEV if the two points in space are closer than the passed value.
 //Applys periodic-boundary conditions.

//Spatial sampling logic---

void Get_Spatial_Sampling_Pams(const double*, const double*, const double*, int, double*, double*);
 //Determines the starting and ending position valules needed to get an 
 //even sampling over the indicated space.
void Get_Plane_Sampling(const int*, const double*, double, double**, int);
 /*Performs a sampling of points in a plane defined by Miller indices and an
   effective origin. Predicts the type of sampling grid (square
   vs. body-centered vs. hexagonal) to use based on Miller indices:
   {100} = square, {110} = rectangular, {111} = hexagonal.*/
void Get_Plane_Sampling(const int*, const double*, double, double**, int, int);
 //Performs a sampling of points in a plane defined by Miller indices, an effective
 //origin, and a requested sampling grid type.
void Get_OrthoRhombic_Sampling(const double*, const double*, const double*,  
	                           const double*, double**, int, int&);
 //Returns an orthorhombic sampling of points in 3D space.
void Get_Spherical_Sampling(const double*, const double*, const double*,  
	                           double, double**, int, int&);
 //Returns a sampling of 3D points inside a sphere.
void Get_Random_OMesh_Sampling(double, double, long int, const double*, double**, int, long int&);
 //Generates a random orthorhombic mesh grid.
void Get_Random_SMesh_Sampling(double, double, long int, double, double**, int, long int&);
 //Generates a random spherical mesh grid.

void Determine_Closest_Vector(const double*, const double*, double*);
 //Takes a Cartesian vector and also a family of vectors (e.g. <100>), and determines 
 //the member of that family that is closest (i.e. makes the smallest angle). 
void Determine_Closest_Vector(const double*, const double*, const double*, double*);
 //Takes a Cartesian vector and also a family of vectors (e.g. <100>), and determines 
 //the member of that family that is closest (i.e. makes the smallest angle) AND located
 //in the indicated plane.
void Determine_Second_Closest_Vector(const double*, const double*, const double*, double*);
 //Takes a Cartesian vector and also a family of vectors (e.g. <100>), and determines 
 //the member of that family that is second closest (i.e. makes the second smallest angle) 
 //AND located in the indicated plane.
void Determine_Closest_Vector(const double*, const double*, const double*, double*, int);
 //Takes a Cartesian vector and also a family of vectors (e.g. <100>), and determines 
 //the member of that family that is n-th closest (i.e. makes the n-th smallest angle) 
 //AND located in the indicated plane; n = 0 (first) or 1 (second) .

double Solve_Plane_Equation(const int*, const double*, const double*);
 //Solves a plane equation in terms of normal vector indices, the plane origin,
 //and a passed location.
double Solve_Plane_Equation(const double*, const double*, const double*);
 //Solves a plane equation in terms of normal vector indices (as FP values), 
 //the plane origin, and a passed location.

void Order_Coordinate_List(double** coor_list, int num_points);
   //Orders a set of points by their location in space.
bool No_Spatial_Prob(const double**, int, const double*, double);
   //Returns TRUEV if the passed point is within a certain distance from any
   //point in a point set.
void Get_Coord_Numbers(const double**, int, int*, double, double);
   //Returns the coordination number for each point in a point set, using
   //coordination sphere logic.

//Rotation logic--- 

void Take_Cross_Product(const double*, const double*, double*);
 //Determines the cross product of two vectors.
void Take_Dot_Product(const double*, const double*, double&);
 //Determines the dot product of two vectors.
void Get_Angle(const double*, const double*, double&);
 //Determines the angle between two vectors.
void Get_Unit_Vec(const double*, double*); 
 //Determines the unit vector of a vector.
void Get_NewMag_Vec(const double*, double*, double);
 //Determines a new vector with differing magnitude from another vector.
void Lengthen_Vec(double*, double);
 //Increases the vector magnitude by the passed amount.
void Convert_RotAngles_To_RotVecs(double, double, double*, double*);
 //Transforms two rotation angles to start vector-->end vector rotation
 //logic.
void Calculate_Rotation_Parameters(const double*, const double*, double*, double&);
 //Returns the rotation vector and angle needed for start vector-->end vector rotation.
void Calc_Rotation_Matrix(double, double, double, double*);
 //Takes three Euler angles in the x-y-z convention and determines the rotation matrix.
void Calc_Rotation_Matrix(const double*, double, double*);
 //Takes rotation vector, rotation angle, and determines the rotation matrix.
void Calc_Rotation_Matrix(const double*, const double*, double*);
 //Takes a start vector-->end vector rotation pair and determines the rotation matrix.
void Calc_Cubic_Plane_Rotation_Matrix(const int*, const int*, double*);
 //Takes integer references for start normal vector-->end normal vector logic
 //and determines the required rotation matrix. Currently only shown to work
 //for Miller indices of cubic systems.
void Rotate_Vector(const double*, double*);
 //Rotates a Cartesian position by the rotation matrix.
void Rotate_Vector(const double*, double*, const double*);
 //Rotates a Cartesian position by the rotation matrix, about an origin.

#endif