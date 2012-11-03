
#ifndef COMM_BEN
#define COMM_BEN

#include <fstream>
#include <ctime>
#include <cstdlib>

using namespace std;

#include "DWN_atom.h"
#include "DWN_atomcoll.h"
#include "DWN_compo.h"
#include "DWN_nanop.h"
#include "DWN_nanos.h"
#include "DWN_utility.h"
#include "DWN_MD.h"
#include "DWN_scope.h"

//File writing---

const bool SIMULATION_POT_LIB = TRUEV;
 //Inclusion of library file lines in simulation files.
 //i.e. "library NANIM.lib" in GULP.

//Files to be used with other programs---

void Write_GULP_File(const atom_collection&, const char*, const char*, bool, bool,
		             bool, bool);
 //Creates a GULP structural optimization file. Uses NPT conditions.
 //Atomic mass information inclusion and core/shell term inclusion
 //are possible.
void Write_GULP_Melt_File(const atom_collection&, const char*, const char*, 
	                      bool, bool, bool, bool, const MD_param&);
 //Writes a GULP file for MD, using a variety of passed MD parameters.

void Write_PDB_File(const atom_collection&, const char*, bool);
 //Writes a PDB file.
 //Reference: http://bmerc-www.bu.edu/needle-doc/latest/atom-format.html#pdb-atom-name-anomalies
void Write_LAMMPS_File(const atom_collection&, const char*, bool);
 //Writes a LAMMPS file.
void Write_XYZ_File(const atom_collection&, const char*, bool);
 //Writes an XYZ atomic coordinate file.
void Write_GSPE_File(const atom_collection&, const char*, bool);
 //Writes a GAUSSIAN '09 input file, for SPE = single point energy calculation.

void Write_QSTEM_File(const atom_collection&, const char*, bool, bool);
 //Writes a QSTEM image simulation file. 
void Write_QSC_File(const char*, const char*, scope_param&, double);
 //Writes a QSC configuration file for QSTEM image simulation. 
void Write_QPBS_File(const char*, const char*);
 //Writes a PBS file for use with QSTEM on the RF research team's supercomputer.
void Write_JEMS_File(const atom_collection&, const char*, bool, bool);
 //Write a JEMS image simulation file.

//Analysis files---

void Write_RAD_File(const atom_collection&, const char*, double);
 //Writes a radial distribution function to file.
void Write_DISTCOUNT_File(const atom_collection&, const char*, double, bool);
 //Writes a distance histogram to file.
void Write_COORCOUNT_File(const atom_collection&, const char*, bool);
 //Writes a coordination number analysis to file, based on variable-size
 //coordination spheres, i.e. determines the coordination number populations for
 //various sizes of the coordination sphere (0 to 2 A, 0 to 3 A, etc.).
void Write_XRD_File(const atom_collection&, const char*, 
					double, bool, double, double, double, double, bool,
					double, double, int);
 //Writes an X-ray diffraction pattern to file, using Debye analysis.
void Write_DWF_File(const atom_collection&, const char*, const char*, bool, double, double);
 /*Writes Debye-Waller Factor analysis (performed for each atom) to file.
   Note: Does not include periodicity effects into spatial coordinate
   analysis, so watch out for outliers.*/

void Interpret_Image(const char*, double, double, int, int*, bool);
 //Analyzes a .img file by outputting a data file for the pixels along with some basic
 //analysis (i.e. maximums, minimums, their averages, etc.). Also performs 
 //an integrated linescan if requested.
 
//File conversion---

void Convert_GULP_To_Collection(atom_collection&, double*, const char*, const atom_collection&, bool);
 //Converts the content of a GULP file into an atomic collection, with optional name
 //preservation, e.g. "C3" atom in GULP file is the "C" atom in this program's user
 //input, but which name should it be assigned in this program's output?
void Convert_GULP_To_Composite(const char*, const char*, const atom_collection&, bool);
 //Converts the content of a GULP file to a composite file.
void Convert_GULP_To_Structure(const char*, const char*, const atom_collection&, bool);
 //Converts the content of a GULP file to a crystal structure file.

void Convert_HIST_To_Collection(atom_collection&, double*, const char*, const atom_collection&, bool,
	                            double, int*, double, bool);
 //Convert .his file into an atomic collection, using an instance (i.e. spatial coordinates
 //at a certain timestep). Patchwork quilting by placing several instances together,
 //side by side, in one output file, is also available.
void Convert_HIST_To_Composite(const char*, const char*, const atom_collection&, 
	                           double, int*, double, bool);
 //Convert .his file to a composite file.
void Convert_HIST_To_Structure(const char*, const char*, const atom_collection&, 
	                           double, int*, double, bool);
 //Converts .his file to a crystal structure file.

void Convert_PDB_To_Collection(atom_collection&, double*, const char*, const atom_collection&, int, bool);
 //Converts a protein database file into an atomic collection, with optional name
 //preservation. Differences in naming between typical PDB files and user
 //input (e.g. "ZN" vs "Zn") are addressed by this function, so the user should not have to worry.
void Convert_PDB_To_Composite(const char*, const char*, const atom_collection&, 
	                          int, bool);
 //Converts a .pdb file to a composite file.
void Convert_PDB_To_Structure(const char*, const char*, const atom_collection&, 
	                          int, bool);
 //Converts a .pdb file to a crystal structure file.

void Convert_Collection_To_Composite(const atom_collection&, const char*);
 //Converts an atomic collection to a composite file.
void Convert_Collection_To_Structure(const atom_collection&, const char*, const double*);
 //Converts an atomic collection to a crystal structure file.

//Partner functions---

void GULP_Ortho_Vector_Output(ofstream&, int, const double*);
 //Outputs the orthorhombic vectors of the simulation box in GULP-file format.
void Load_DW_Factors(double**);
 //Debye-Waller factors loaded, by atomic number, from file.
void Load_Scattering_Coefficients(double**, bool);
 //Scattering factor calculation coefficients loaded, by atomic number, from file.

#endif
