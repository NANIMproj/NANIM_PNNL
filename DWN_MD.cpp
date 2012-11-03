#include "stdafx.h"
#include "DWN_MD.h"

//MD_param structure---

//Assignment operator---

MD_param& MD_param::operator = (const MD_param& mpam)
   {
   time_step = mpam.time_step;
   initial_temp = mpam.initial_temp;
   final_temp = mpam.final_temp;
   initial_temp_time = mpam.initial_temp_time;
   change_temp_time = mpam.change_temp_time;
   thermo_coeff = mpam.thermo_coeff;
   baro_coeff = mpam.baro_coeff;
   equi_time = mpam.equi_time;
   prod_time = mpam.prod_time;
   supdate_time = mpam.supdate_time;
   fupdate_time = mpam.fupdate_time;
   nvt_run = mpam.nvt_run;
   npt_run = mpam.npt_run;
   molecule_commands = mpam.molecule_commands;

   return *this;
    //Proper return for this operation.
   }

//Molecular dynamics class--

//Constructors, destructors and initialization---

molecule_sim::molecule_sim()
   {
   Initialize();
   }

 molecule_sim::molecule_sim(const MD_param& pam)
   {
   mpam = pam;
   }
   
 molecule_sim::~molecule_sim()
   {

   }

void molecule_sim::Initialize()
 //Sets MD simulation parameters to default.
   {
   mpam.time_step = 1.0;
   mpam.initial_temp = 298;
   mpam.final_temp = 298;
    //No dynamic temperature.
   mpam.initial_temp_time = 100.0;
   mpam.change_temp_time = 0.0;
   mpam.thermo_coeff = 0.05;
   mpam.baro_coeff = 0.05;
   mpam.equi_time = 100.0;
   mpam.prod_time = 100.0;
   mpam.supdate_time = 1.0;
   mpam.fupdate_time = 1.0;
   mpam.nvt_run = TRUEV;
   mpam.npt_run = FALSEV;
   mpam.molecule_commands = FALSEV;
   }

//Properties setting/retrieval---

void molecule_sim::Get_Parameters(MD_param& pam) const
   {
   pam = mpam;
   }

void molecule_sim::Set_Parameters(const MD_param& pam)
   {
   mpam = pam;
   }

//General MD Loop---

void Start_MD(atom_collection& atoms, MD_param& parameters)
	{
	int step_count = 0;
	bool simulation_done = FALSEV;

	while (!simulation_done)
		{
		if (step_count == 0)
		 //First time step performed with Euler step.
		   {


		   }

		else
		 //All other time steps performed with Verlet step.
		   {

		   }


		}

	}

