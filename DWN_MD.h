
#ifndef EGGCART_JENSEN
#define EGGCART_JENSEN

#include "DWN_compo.h"

typedef double DW_data;

//Structure for storage of molecular dynamics parameters---

struct MD_param
{
double time_step;
 //Time step used in equation integration, in femtoseconds.
double initial_temp;
 //Initial temperature used in simulation, in K.
double final_temp;
 //Final temperature at last time step, in K.
double initial_temp_time;
 //Time spent simulating at the initial temperature, in picoseconds.
double change_temp_time;
 //Time interval over which initial temp->final temp
 //temperature change is incrementally made, in ps.
double thermo_coeff;
 //Nose-Hoover thermostat coefficient.
double baro_coeff;
 //Parrinello-Rahman barostat coefficient.
double equi_time;
 //Time spent for the initial system equilibration period, in ps.
double prod_time;
 //Time spent in production (averaging of properties of interest), in ps.
double supdate_time;
 //How often the screen (simulation prompt) is updated, in ps.
double fupdate_time;
 //How often the output file is updated, in ps.
bool nvt_run;
 //Indicates if the NVT ensemble is being used.
bool npt_run;
 //Indicates if the NPT ensemble is being used.
bool molecule_commands;
 //Indicates if molecule commands are to be included. These are: 
 //Coulombic subtraction, connectivity (bond) fixing, and internal 
 //kinetic energy removal from molecules.
 

MD_param& operator = (const MD_param&);
 //Assigns all the parameters.
};

//Dedicated class to molecular simulations---

class molecule_sim
   {
   public:

   //Constructors, destructors, and initialization functions---

   molecule_sim();
   molecule_sim(const MD_param&);
    //Gets the simulation parameters.
   ~molecule_sim();

   void Initialize();
    //Sets molecular simulation to default parameters.

   //Properties setting/retrieval---

   void Get_Parameters(MD_param&) const;
   void Set_Parameters(const MD_param&);

   //Molecular dynamics routine---

   void Start_MD(atom_collection&, MD_param&);
    //Starts an MD run with the passed atom collection
    //and MD parameters. 

   //Potential energy calculations---


   //Molecular mechanics---

   //void Rigid_Minimization(atom_collection&, int, int, bool, double);
    //Structural optimization function. Minimizes the energy of the range of atoms
    //indicated, treating those atoms as a rigid body. Only translation and rotation
    //of that body is allowed during energy minimization.

   //Trial moves---

   //void Trial_Move_Translate_Get_Gradient(int, int, double*, double);
    //Determines the local potential gradient for translation, in eV/nm.
   //void Trial_Move_Rotate_Get_Gradient(int, int, double*, double);
    //Determines the local potential gradient for rotation, in eV/radian.

   //Potential calculations---

   //double Calculate_Fast_Total_Potential(bool) const;
    //Gets the total potential energy of the system, as fast as possible.
   //double Calculate_Two_Body_Potential(bool) const;
	//Sums up all two-body potential energies in a system.
   //double Calculate_SR_Two_Body_Potential(bool) const;
	//Sumps up short-range two-pody potentail energies in a system.
	  

   private:

   MD_param mpam;
    //Set of simulation parameters.

   };

#endif