
#ifndef WHAT_IS_CONSTANT
#define WHAT_IS_CONSTANT

//Useful constants---

//Conditionals---

const bool FALSEV = 0;
const bool TRUEV = 1;
 //Boolean values.

//Math/basic geometry constants---

const double PI_CONST = 3.141592653589793238;
 //Value of PI.
const double ZERO_LOW_BOUND = -0.0000001;
const double ZERO_UPP_BOUND =  0.0000001;
 //Range of floating-point values that is assumed to be equivalent to zero.
const double FP_ERROR_FIX = 0.0000005;
 //A value to add/subtract from values as to keep the FP math correct when needed.
const double DIST_CONV = 10.0;
 //Converts the value of distance used in this program to Angstroms.

const int SPHERE_SHAPE = 0;
const int BOX_SHAPE = 1;
 //Shapes for geometry fitting.

const int SC_GRID = 0;
const int BCC_GRID = 1;
const int FCC_GRID = 2;
 //Grid factors for sampling points in a plane.

//String processing constants---

const int MAX_I_DIGITS = 9;
const int MAX_POI_DIGITS = 8;
const int MAX_FP_DIGITS = 12;
 //Maximum number of digits in c-string representations of numbers.
 
const int MAX_COMMAND_LENGTH = 100;
 //Expected maximum character length of a command in the scripting language
 //such as "Define Nanostructure."
const int MAX_PHRASE_LENGTH = 100;
 //Expected maximum character length of a phrase in the scripting language
 //such as "Silicon."
const int MAX_FILE_NAME_LENGTH = 100; 
 //Maximum length of a file name created anywhere in this program.
const int EXT_LEN = 5; 
 //Size of character array used to hold a file extension string.

//Molecular geometry constants---

//const int OCTO = 6;
const int TETRA = 4;
const int TRIG = 3;
const int LINEAR = 2;
 //Reference indices for molecular geometries.
const int VERTEX_A = 0;
const int VERTEX_B = 1;
const int VERTEX_C = 2;
const int VERTEX_D = 3;
 //Index references for the vertices within these geometries.

/*const double OCTO_POS[6][3] = { {1.0, 0.0, 0.0} ,
                                {-1.0, 0.0, 0.0} , 
                                {0.0, 1.0, 0.0} ,
                                {0.0, -1.0, 0.0} ,
                                {0.0, 0.0, 1.0} ,
                                {0.0, 0.0, -1.0} ,  */                 
 //One set of points in a unit-vertex octahedron centered around an origin.

const double TETRA_POS[4][3] = { {0.57735, 0.57735, 0.57735} ,
                                 {-0.57735, -0.57735, 0.57735} , 
                                 {-0.57735, 0.57735, -0.57735} ,
                                 {0.57735, -0.57735, -0.57735} };
 //One set of points in a unit-vertex tetrahedron centered around an origin.
const double TRIG_POS[3][3] = { {0.86603, 0.5, 0.0} , 
                                {-0.86603, 0.5, 0.0} , 
                                {0.0, -1.0, 0.0} };
 //Same for trigonal planar shape about an origin.
const double LINEAR_POS[2][3] = { {1.0, 0.0, 0.0} , 
                                  {-1.0, 0.0, 0.0} };
 //Same for linear shape about an origin.

//Force/potential calculation constant values---

const double ELEMENTARY_CHARGE = 1.602176565 * pow(10.0, -19.0);
const double ELECTRIC_CONSTANT = 8.854187817 * pow(10.0, -12.0);
const double INV_4PI_ELEC_CONST = 1.0/(4.0*PI_CONST*ELECTRIC_CONSTANT);
const double ELEM_CHARGE_SQUARE = ELEMENTARY_CHARGE*ELEMENTARY_CHARGE;
const double NM_PER_M = 1.00 * pow(10.0, -9.0);
const double ELEC_POT_CONSTANT_JOULES = INV_4PI_ELEC_CONST*ELEM_CHARGE_SQUARE/NM_PER_M;
const double ELEC_POT_CONSTANT_EV = ELEC_POT_CONSTANT_JOULES/ELEMENTARY_CHARGE;
 
//Output method variables---

const bool COOR_LIST_ORDERING = FALSEV;
 //Indicates if coordinate output should be ordered 
 //such that points should be located near other points
 //that are close to them.

#endif
