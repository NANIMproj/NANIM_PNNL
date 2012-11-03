
#include "stdafx.h"

#include "DWN_prog.h"

void Start()
      {  
      Show_Statement("WELCOME TO THE PROGRAM! YAY!\n");
      const char FILE_NAME[50] = "CommandFile\\TalkToMe.txt";
	  const char SUPERAUTO_NAME[50] = "CommandFile\\Monologuer.txt";
	  const char MEGAAUTO_NAME[50] = "CommandFile\\Cast.txt";
	  
	  //Main script file setup---

	  ifstream file_reader(FILE_NAME, ios::in);
	  if (!file_reader.is_open())
           {
		   Show_Warning("COULD NOT OPEN MAIN SCRIPT FILE!");
           return;  
	       }
      int file_size = 0;
      char* file_text;
      file_text = new char[MAX_SCRIPT_FILE_SIZE];
      Get_Characterized_File(file_reader, file_text, MAX_SCRIPT_FILE_SIZE, file_size);
        //Dump the script file into a character array for command processing.
	  file_reader.close();
            
	  //Superauto setup---

	  int max_num_runs = 1;
	  int nauto_file_size = 0;
	  char* nauto_file_text = new char[MAX_AUTO_SIZE];
	  char dump_string[20];
	  ifstream nauto_reader(SUPERAUTO_NAME, ios::in);
	  if (nauto_reader.is_open())
	     {
		 nauto_reader >> dump_string >> max_num_runs;
		 if (max_num_runs < 1)
		  //Default = 1 run = No automation.
		    {
			max_num_runs = 1;
		    }
		 Get_Characterized_File(nauto_reader, nauto_file_text, MAX_AUTO_SIZE, 
			                    nauto_file_size);
		 nauto_reader.close();
	     }

	  //Megaauto setup--- (3rd order of automation)

	  int max_string_runs = 1;
	  int sauto_file_size = 0;
      char* sauto_file_text = new char[MAX_AUTO_SIZE];
	  ifstream sauto_reader(MEGAAUTO_NAME, ios::in);
	  if (sauto_reader.is_open())
	     {
		 sauto_reader >> dump_string >> max_string_runs;
		 if (max_string_runs < 1)
		    {
			max_string_runs = 1;
		    }
		 Get_Characterized_File(sauto_reader, sauto_file_text, MAX_AUTO_SIZE,
							    sauto_file_size);
		 sauto_reader.close();
	     }

	  //Process main script commands repeatedly, applying requested changes to
	  //the script via automation---

	  int num_run_count, string_run_count;
	  for (string_run_count = 0; string_run_count < max_string_runs; ++string_run_count)
	     {
		 for (num_run_count = 0; num_run_count < max_num_runs; ++num_run_count)
	       {
		   Script_Reader(file_text, file_size, nauto_file_text, nauto_file_size,
			             sauto_file_text, sauto_file_size, num_run_count, string_run_count);
	       }
	     }

	  //Free memory---

	  delete[] file_text;
	  delete[] nauto_file_text;
	  delete[] sauto_file_text;

	  Show_Statement("\nDone all that reading! Dang I'm tired...need some sleep. Please leave now.");
      system("PAUSE");
      }

void Script_Reader(const char* file_text_script, int file_size,
	               const char* nauto_text_script, int nauto_size,
				   const char* sauto_text_script, int sauto_size,
				   int num_run_count, int text_run_count)
        {
        //Make character-array copies of the changeable scripting files
	    //and apply automation commands---

        char* file_text = new char[MAX_SCRIPT_FILE_SIZE];
		for (int a = 0; a < file_size; ++a)
		    {
		    file_text[a] = file_text_script[a];
		    }
		char* nauto_text = new char[MAX_AUTO_SIZE];
		for (int a = 0; a < nauto_size; ++a)
		    {
		    nauto_text[a] = nauto_text_script[a];
		    }
		
		Show_Statement("\nStarting textual script run # ", text_run_count + 1);
		Show_Statement("\nStarting numerical script run # ", num_run_count + 1);

		Process_SAutomation_Commands(file_text, file_size, nauto_text,
				                     nauto_size, sauto_text_script, 
									 sauto_size, text_run_count);
		Process_NAutomation_Commands(file_text, file_size, nauto_text, 
			                         nauto_size, num_run_count);

		//Following list is of ALL variables used in the reading of the script
		//file. Get ready for a long read.

		 const char* char_pointer;
           //Pointer to current position in the script text.
         double temp_coors[3];
           //Temporary spatial locations.
       
         //Max species counts---
         
         const int MAX_ATOMS = 200;
         const int MAX_CRYST = 30;
         const int MAX_MOLEC = 50;
         const int MAX_NANOS = 20;
         const int MAX_NANOP = 30;
          //Maximum number of definitions in the script text.
          
         //File extensions---
         
		 /*Extensions specific to this program.
		   Note: These are not defined as constant strings,
		   but they are constant. This is the case for
		   most of the following character strings.*/

		 char cryst_ext[EXT_LEN] = "crs";
         char molecule_ext[EXT_LEN] = "mol";
         char nanoparticle_ext[EXT_LEN] = "nnp";
         char nanosurface_ext[EXT_LEN] = "nns";
         char composite_ext[EXT_LEN] = "cop";
		  //Atomic systems.
         char radial_ext[EXT_LEN] = "rad"; 
		 char distcount_ext[EXT_LEN] = "dis";
		 char coorcount_ext[EXT_LEN] = "cor";
		 char analysis_ext[EXT_LEN] = "txt";
		 char xrd_ext[EXT_LEN] = "xrd";
		 char dwf_ext[EXT_LEN] = "dwf";
		  //Analysis files.

		 //Extensions of other programs---

         char GULP_ext[EXT_LEN] = "gin";
		 char RES_ext[EXT_LEN] = "res";
         char GULPDUMP_ext[EXT_LEN] = "gin";
         char GULPTRAJ_ext[EXT_LEN] = "trg";
		 char GRO_ext[EXT_LEN] = "gro";
         char PDB_ext[EXT_LEN] = "pdb";
		 char LAMMPS_ext[EXT_LEN] = "lam";
		 char QSTEM_ext[EXT_LEN] = "cfg";
		 char QCONF_ext[EXT_LEN] = "qsc";
         char JEMS_ext[EXT_LEN] = "jem";
		 char XYZ_ext[EXT_LEN] = "xyz";
		 char PBS_ext[EXT_LEN] = "pbs";
		 char GSPE_ext[EXT_LEN] = "gjf";
		 char FORCE_ext[EXT_LEN] = "frc";
		 char TRJ_ext[EXT_LEN] = "trj";
		 char HIST_ext[EXT_LEN] = "his";
		 char IMG_ext[EXT_LEN] = "img";

		 //Names of files output by this program---

		 const char ANALYSIS_FILE_NAME[24] = "Results\\Listener.txt";
          
		 //File conversion variables---
         
         char convert_signature[14] = "Convert File";
         char* old_file_pointer;
         char* new_file_pointer;
         char temp_old_ext[5];
         char temp_new_ext[5];
		 int model_num_temp;
		 bool name_preservation;
         double time_temp;
		 int patchwork_temp[3];
		 double patchwork_spacing_pam;

         //Atom definition variables---
          
         char atom_signature[13] = "Define Atom";
         atom_collection new_atoms(MAX_ATOMS + 1);
		  //Include extra space for dummy atom.
         char* temp_atom_names[MAX_ATOMS + 1];
         int atom_num_temp; 
         double atom_weight_temp, atom_charge_temp, atom_radius_temp;
         int atom_index = 0;
		  //Number of atoms defined by user.

         //Crystal structure variables---
         
         char crystal_signature[18] = "Define Structure";
         char crystal_load_signature[16] = "Load Structure";
         char latticepam_signature[12] = "LatticePam";
		 char meshpam_signature[9] = "MeshPam";
         char crystal_atom_signature[6] = "Atom";
         char facet_signature[7] = "Facet";
         crystal_system new_crystal_systems[MAX_CRYST];
         char* struct_names[MAX_CRYST];
         double lattice_parameters[9];
		 int num_lattice_parameters_set;
         atom temp_basis_atom;
         bool rel_coor;
         int temp_index;
         char* temp_type_pointer;
         char* name_pointer;
         bool is_special_system;
         int temp_facet_indices[3];
         double temp_scaling_factor;
		 const char* temp_facet_name_pointer;
         int struct_index = 0;

		 //Supercell variables---
         
         char SUPERCELL_signature[11] = "SUPERCELL";
         char super_ext[EXT_LEN] = "sup";
         char* temp_super_name;
         char* temp_cell_struct;
         char temp_super_string[50];
         int temp_super_index;
         int temp_super_sizes[3];
         
         //Molecule variables---
         
         char molecule_signature[17] = "Define Molecule";
         char molecule_load_signature[15] = "Load Molecule";
         char add_tomol_signature[6] = "Add ";
         char water_signature[7] = "Water";
         char alkanethiol_signature[13] = "Alkanethiol";
         char carbonchain_signature[10] = "Carbbone";
         char OH_signature[4] = "OH";
         char SH_signature[4] = "SH";
         char CH2_signature[5] = "CH2";
         char CH3_signature[5] = "CH3";
         char AX_signature[4] = "AX";
         char AXY_signature[5] = "AXY";
         char AXYZ_signature[6] = "AXYZ";
         molecule new_molecules[MAX_MOLEC];
         char* molecule_names[MAX_MOLEC];
         char* temp_tomol_name;
         int temp_atom_index;
         char* temp_atom_name_list[10];
         double temp_bond_lengths[10];
         double temp_start[3];
         double temp_previous_end[3];
         double temp_coor_spaceA[3];
         double temp_coor_spaceB[3];
         int temp_shape_factor;
         int temp_carbon_count;
         int temp_atom_indices[10];
         int molecule_index = 0;
         
         //Nanosurface variables---
         
         char nanosurface_signature[20] = "Define NanoSurface";
         char nanosurface_load_signature[18] = "Load NanoSurface";
         nanosurface new_nanosurfaces[MAX_NANOS];
         char* nanosurface_names[MAX_NANOS];
         int temp_structure_index;
         int temp_indices[3];
		 double temp_angle;
         double temp_sizes[3];
         bool stuffed_surface;
		 bool vary_size_for_neutrality;
		 bool complete_neutralization_needed;
		 double temp_shift_parameters[2];
		 double temp_randomize_parameters[3];
		 double temp_sliceslab_parameters[2];
		 bool staircase_shift_request;
		 bool randomize_request;
		 bool sliceslab_request;
         double temp_s_coat_parameters[3];
         char* temp_s_coat_name;
         int temp_s_coat_index;
         bool surface_coat_type;
         int surface_index = 0;
         
         //Nanoparticle variables---

         char nanoparticle_signature[21] = "Define NanoParticle";
         char nanoparticle_load_signature[19] = "Load NanoParticle";
         nanoparticle new_nanoparticles[MAX_NANOP];
         char* nanoparticle_names[MAX_NANOP];
         int temp_p_structure_index;
         double temp_diameter;
		 double temp_hole_diameter;
         int temp_model;
         bool stuffed_particle;
		 bool neutralize_by_surface_removal;
         double temp_p_coat_parameters[2];
         char* temp_p_coat_name;
         int temp_p_coat_index;
         bool particle_coat_type;
         int particle_index = 0;
         
         //General composite variables---
         
         char gencomp_signature[16] = "Make Composite";
         char gencomp_load_signature[16] = "Load Composite";
         char add_signature[6] = "Add";
		 char sub_signature[6] = "Sub";
		 char bulldoze_signature[6] = "Bul";
		 char invsub_signature[6] = "Inv";
		 char chargebal_signature[11] = "ChargeBal";
         char solvate_signature[9] = "Solvate";
		 char coordex_signature[9] = "Coordex";
         general_composite new_composite;
		 char* temp_file_pointer = NULL;
		 double temp_positions[3];
		 bool temp_mode;
         char* composite_name = NULL;  
         double temp_comp_sizes[3];  
         char* temp_component_name;
         bool temp_add_mode;
         double temp_angles[3];
         int temp_periodicity;
         int search_index;
         int temp_supercell_sizes[3];
		 double temp_void_origin[3];
		 double temp_void_size[3];
         double temp_solvent_box_origin[3];
         double temp_solvent_box_size[3];
		 double temp_chargebal_box_origin[3];
		 double temp_chargebal_box_size[3];
		 double coordex_pams[2];
		 
		 //High-priority (first to be processed) general commands---
         
         char NO_HYDROGEN_IMAGE_signature[12] = "IM_REMOVEH";
         char NO_HYDROGEN_MD_signature[12] = "MD_REMOVEH";
         char NO_HYDROGEN_signature[9] = "NOHYDRO";
		 char ANAL_HYDROGEN_signature[15] = "HYDRO_ANALYZE";
		 char MD_MASS_signature[6] = "MASS";
         bool image_hydrogen, MD_hydrogen, MD_mass, anal_hydrogen;
		 bool super_hydrogen_remove = FALSEV;
		  //Master command. Can only be turned on in code.
		 char CAT_CORESHELL_signature[15] = "CAT_CORESHELL";
		 char ANI_CORESHELL_signature[15] = "ANI_CORESHELL";
		 char CORESHELL_signature[15] = "POL_CORESHELL";
		 char NORAD_signature[7] = "NORAD";
		 char CLEANUP_signature[9] = "CLEANUP";
		 bool ani_coreshell, cat_coreshell, no_rad_calc, clean_up;
         atom_collection temp_collection;

		 //MD parameters---

		 char RUNMD_signature[7] = "RUNMD";
		 char LOADMDPAM_signature[8] = "LOADMD";

		 char MDPAM_signature[15] = "MD Parameters";
		 char TIMESTEP_signature[4] = "dt";
		 char INITIALTEMP_signature[4] = "T0";
		 char FINALTEMP_signature[4] = "Tf";
		 char INITIALTEMPTIME_signature[6] = "T0_t";
		 char FINALTEMPTIME_signature[6] = "Tf_t";
		 char NVT_signature[5] = "NVT";
		 char NPT_signature[5] = "NPT";
		 char THERMOSTAT_signature[12] = "Thermostat";
		 char BAROSTAT_signature[10] = "Barostat";
         char EQUI_signature[6] = "Equi";
		 char PROD_signature[6] = "Prod";
		 char SCREENUPDATE_signature[6] = "Show";
		 char FILEUPDATE_signature[6] = "File";
		 char MOLECULECOMMAND_signature[10] = "Molecule";

		 MD_param mpam;
		  //Stores various MD parameters, from timestep to user update frequency.

		 char RIGIDMIN_signature_DEFAULT[10] = "RIGIDMIN";
		 char RIGIDMIN_signature[10] = "RigidMin";
		 int first_region_added[2];
		 int minimization_region[2];
		 double rigid_precision;
		 
		 //Microscope parameters---

		 char SCOPEPAM_signature[23] = "Microscope Parameters";
		 char SCOPEMODE_signature[6] = "Mode";
		 char Nx_signature[4] = "Nx";
		 char Ny_signature[4] = "Ny";
		 char Nz_signature[4] = "Nz";
		 char Nc_signature[3] = "Nc";
		 char Tx_signature[4] = "Tx";
		 char Ty_signature[4] = "Ty";
		 char Tz_signature[4] = "Tz";
		 char Tc_signature[3] = "Tc";
		 char ScanXMin_signature[10] = "ScanXMin";
		 char ScanXMax_signature[10] = "ScanXMax";
		 char RasterX_signature[10] = "RasterX";
		 char ScanYMin_signature[10] = "ScanYMin";
		 char ScanYMax_signature[10] = "ScanYMax";
		 char RasterY_signature[10] = "RasterY";
		 char ScanMin_signature[9] = "ScanMin";
		 char ScanMax_signature[9] = "ScanMax";
		 char Raster_signature[9] = "Raster";
		 char WindowSizeX_signature[13] = "WindowSizeX";
		 char WindowSizeY_signature[13] = "WindowSizeY";
		 char WindowSize_signature[12] = "WindowSize";
		 char XRes_signature[13] = "ResolutionX";
		 char YRes_signature[13] = "ResolutionY";
		 char Res_signature[12] = "Resolution";
		 char Slice_signature[16] = "SliceThickness";
		 char VOLTAGE_signature[9] = "Voltage";
		 char CURRENT_signature[9] = "Current";
		 char Brightness_signature[12] = "Brightness";
		 char DWELLTIME_signature[7] = "Dwell";
		 char SOURCEANGLE_signature[13] = "SourceAngle";
		 char SourceSize_signature[12] = "SourceSize";
		 char BeamTiltX_signature[11] = "BeamTiltX";
		 char BeamTiltY_signature[11] = "BeamTiltY";
		 char BeamTilt_signature[10] = "BeamTilt";
		 char ConvAngle_signature[11] = "ConvAngle";
		 char Defocus_signature[9] = "Defocus";
		 char DetectorIAngle_signature[16] = "DetectorIAngle";
		 char DetectorOAngle_signature[16] = "DetectorOAngle";
		 char CS3_signature[5] = "Cs3";
		 char CS5_signature[5] = "Cs5";
		 char CC_signature[4] = "Cc";
		 char DE_signature[4] = "dE";
		 char ASTIG2_signature[8] = "2Astig";
		 char ASTIG3_signature[8] = "3Astig";
		 char ASTIGANGLE_signature[7] = "Angle";
		 char SCOPETEMP_signature[6] = "Temp";
		 char TDScount_signature[5] = "TDS";
		 char SliceOutput_signature[10] = "SliceOut";
		 char IDEALIMAGE_signature[7] = "IDEAL";

		 scope_param spam;
		 double res_temp;

         //Low-priority general commands----
         
         char GULP_signature[10] = "GULPOPTI";
         char GULPMELT_signature[10] = "GULPMELT";
		 char GRO_signature[9] = "GROFILE";
		 char PDB_signature[9] = "PDBFILE";
		 char LAMMPS_signature[8] = "LAMMPS";
         char QSTEM_signature[7] = "QSTEM";
		 char QCONF_signature[7] = "QCONF";
		 char QPBS_signature[6] = "QPBS";
         char JEMS_signature[6] = "JEMS";
		 char XYZ_signature[9] = "XYZFILE";
		 char GSPE_signature[6] = "GSPE";
         char FILEDUMP_signature[10] = "FILEDUMP";
         char RADIAL_signature[8] = "RADIAL";
		 char DISTCOUNT_signature[11] = "DISTCOUNT";
		 char COORCOUNT_signature[11] = "COORCOUNT";
		 char ENERGY_signature[12] = "SHOWENERGY";
		 char ELECOMP_signature[10] = "ELECOMP";
		 char ANALYSIS_signature[10] = "ANALYSIS";
		 char CLEANANAL_signature[11] = "CLEANANAL";
		 char STRANGE_signature[9] = "STRANGE";
		 char XRD_signature[9] = "XRDCALC";
		 char IONICSCAT_signature[11] = "IONICSCAT";
		 char TRAJDWF_signature[9] = "TRAJDWF";
		 char IMAGEDWF_signature[9] = "IMAGEDWF";
		 char ELECTROPOT_signature[12] = "ELECTROPOT";
		 char ENERGYPOT_signature[11] = "ENERGYPOT";
         char temp_string[MAX_FILE_NAME_LENGTH];
		 char temp_string2[MAX_FILE_NAME_LENGTH];
         char temp_dump_string[MAX_FILE_NAME_LENGTH];
		 double temp_box_size[3];

		 double energy_value;
		 bool energy_set;

		 const int ELE_STRING_LENGTH = 500;
		 char ele_fit_string[ELE_STRING_LENGTH];
		 bool ele_comp_set;

		 char radial_integration_parameter_signature[10] = "RDF Step";
		 double radial_step;

		 char histogram_parameter_signature[11] = "Histogram";
		 double histo_pam;
		 
		 char strangeness_parameter_signature[9] = "Strange";
		 double strangeness_dist_min, strangeness_dist_max;

		 char cleanup_parameter_signature[7] = "Clean";
		 int cleanup_pam;

		 char XRD_wave_signature[8] = "XRWave";
		 char XRD_min_signature[7] = "XRMin";
		 char XRD_max_signature[7] = "XRMax";
		 char XRD_precision_signature[13] = "XRPrecision";
		 char XRD_DW_signature[6] = "XRDW";
		 char XRD_surface_signature[10] = "XRSurfDW";
		 char XRD_coornum_signature[10] = "XRBulkCN";
		 char XRD_coorcutoff_signature[11] = "XRSurfCut";
		 double xrd_wave, xrd_precision, xrd_min_angle, xrd_max_angle, xrd_dw_factor;
		 double xrd_dw_factor_surface_ratio, xrd_coor_cutoff;
		 int xrd_coor_num;
		 bool xrd_ionic, xrd_surface_effects;
		 char distavg_signature[18] = "Distance Average";
		 double min, max;
		 bool unlike_only, like_only;
		 bool unlike_distavg_set, like_distavg_set, gen_distavg_set;
		 double unlike_distavg, like_distavg, gen_distavg, distavg;
		 double unlike_cooravg, like_cooravg, gen_cooravg, cooravg;

		 bool image_dwf;

		 char shapefit_signature[11] = "Fit Shape";
		 int shape_identity;
		 int coor_num;
		 int shapefit_atom_number_specification;
		 double sphere_fit;
		 double box_fit[3];
		 bool sphere_fit_set;
		 bool box_fit_set;

		 bool strangeness_set;
		 bool is_strange;

		 bool electro_potential_set;
		 double electro_potential;

		 bool energy_potential_set;
		 double energy_potential;

		 bool analysis_request;

		 char dw_coordination_parameter_signature[17] = "DW Coordination";
		 double dw_coordination_min, dw_coordination_max;

		 //File reading variables---

		 char IMAGEREAD_signature[12] = "Read Image";
		 char* temp_image_file_pointer;
		 double gun_size, image_stat_pam;
		 int restriction_pam;
		 bool show_num;
		 int line_scan_params[3];
		 
		 //END VARIABLE DEFINITIONS---

         //Process atom definitions---
         
         for (int index = 0; index < file_size; ++index)
          //Look for atom definitions, which are allowed to be anywhere in the file.
             {
			 if (atom_index == MAX_ATOMS)
				{
				Show_Warning("YOU HAVE TOO MANY ATOMS DEFINED!");
				break;
				}

             char_pointer = &(file_text[index]);
              //Investigate a particular character in the text file.
             if (Command_Check(char_pointer, atom_signature))
              //If true, the function has found a "Define Atom" command in the 
              //file input, which is followed by atomic information.
                {
				MPSNL_CommandCheck_SC(char_pointer, temp_atom_names[atom_index], "Name");
                 //Record the name of the atom.
				MPNL_CommandCheck_SC(char_pointer, "Symbol");
                new_atoms[atom_index].Set_Atom_Name(char_pointer);
                 //Record the symbol of the atom.
				MPNL_CommandCheck_GetPosNum(char_pointer, "AN", atom_num_temp);
                 //Record the atomic number.
				MPNL_CommandCheck_GetFPNum(char_pointer, "AW", atom_weight_temp);
                 //Record the atomic weight.
				MPNL_CommandCheck_GetFPNum(char_pointer, "CH", atom_charge_temp);
                 //Record the atomic charge.
				MPNL_CommandCheck_GetFPNum(char_pointer, "RA", atom_radius_temp);
				atom_radius_temp /= 10.0;
				 //Convert radius to nanometers.
                new_atoms[atom_index].Set_Atom_Properties(atom_num_temp, atom_weight_temp, 
                                                          atom_charge_temp, atom_radius_temp);
                ++new_atoms;
                ++atom_index;
                }
              }

		 Show_Statement("\nDone defining ", atom_index, " atoms!");
		 if (atom_index == 0)
		    {
		    Show_Warning("NO ATOMS ARE DEFINED IN YOUR INPUT FILE!");
		    }

         //Add Dummy atom---

         new_atoms[atom_index].Set_Atom_Name(const_cast<char*>(DUMMY_NAME));
         ++new_atoms;
         ++atom_index;     
                 
         //File conversion commands--- 
         
         for (int index = 0; index < file_size; ++index)
             {
             char_pointer = &(file_text[index]);
             if (Command_Check(char_pointer, convert_signature))
                {
                Move_Pointers_NextLine(char_pointer, old_file_pointer);
                Move_Pointers_NextLine(char_pointer, new_file_pointer);
                Get_Extension(temp_old_ext, old_file_pointer);
                Get_Extension(temp_new_ext, new_file_pointer);            
                 //Got the file names and extensions of those file names. 

				name_preservation = TRUEV;
				if (MPNL_CommandCheck_SC(char_pointer, "PreserveName"))
				      {
					  if (char_pointer[0] == 'N')
					     {
						 name_preservation = FALSEV;
					     }
				      }

                if (Name_Check(temp_old_ext, GULP_ext) && Name_Check(temp_new_ext, composite_ext) )
                 //GULP-->GENERAL COMPOSITE conversion.
                   {
				   Show_Statement("GIN File ", old_file_pointer, " being converted...!");
                   Convert_GULP_To_Composite(old_file_pointer, new_file_pointer, 
					                         new_atoms, name_preservation);                          
                   } 
				else if (Name_Check(temp_old_ext, RES_ext) && Name_Check(temp_new_ext, composite_ext) )
				   {
				   Show_Statement("RES File ", old_file_pointer, " being converted!");
                   Convert_GULP_To_Composite(old_file_pointer, new_file_pointer, 
					                         new_atoms, name_preservation);   
				   }
				 else if (Name_Check(temp_old_ext, GULP_ext) && Name_Check(temp_new_ext, cryst_ext) )
                 //GULP-->GENERAL COMPOSITE conversion.
                   {
				   Show_Statement("GIN File ", old_file_pointer, " being converted!");
                   Convert_GULP_To_Structure(old_file_pointer, new_file_pointer, 
					                         new_atoms, name_preservation);                          
                   } 
				else if (Name_Check(temp_old_ext, RES_ext) && Name_Check(temp_new_ext, cryst_ext) )
				   {
				   Show_Statement("RES File ", old_file_pointer, " being converted!");
                   Convert_GULP_To_Structure(old_file_pointer, new_file_pointer, 
					                         new_atoms, name_preservation);                          
				   }
				else if (Name_Check(temp_old_ext, HIST_ext) && Name_Check(temp_new_ext, composite_ext))
				   {
				   MPNL_CommandCheck_GetFPNum(char_pointer, "Time", time_temp);
				   patchwork_temp[0] = 0;
				   patchwork_temp[1] = 0;
				   patchwork_temp[2] = 0;
				   MPNL_CommandCheck_GetNums(char_pointer, "Quilt", patchwork_temp, 3); 
				   MPNL_CommandCheck_GetFPNum(char_pointer, "Spacing", patchwork_spacing_pam, 0.0);
				   Show_Statement("His file ", old_file_pointer, " being converted ");
				   Show_Statement("Time interval for conversion: ", time_temp);
				   Convert_HIST_To_Composite(old_file_pointer, new_file_pointer, new_atoms, time_temp, 
					                         patchwork_temp, patchwork_spacing_pam, name_preservation);
				   }
				else if (Name_Check(temp_old_ext, HIST_ext) && Name_Check(temp_new_ext, cryst_ext))
				   {
				   MPNL_CommandCheck_GetFPNum(char_pointer, "Time", time_temp);
				   patchwork_temp[0] = 0;
				   patchwork_temp[1] = 0;
				   patchwork_temp[2] = 0;
				   MPNL_CommandCheck_GetNums(char_pointer, "Quilt", patchwork_temp, 3);
				   MPNL_CommandCheck_GetFPNum(char_pointer, "Spacing", patchwork_spacing_pam, 0.0);
				   Show_Statement("His file ", old_file_pointer, " being converted ");
				   Show_Statement("Time interval for conversion: ", time_temp);
				   Convert_HIST_To_Structure(old_file_pointer, new_file_pointer, 
					                         new_atoms, time_temp, patchwork_temp, patchwork_spacing_pam, name_preservation);
				   }
				else if (Name_Check(temp_old_ext, PDB_ext) && Name_Check(temp_new_ext, composite_ext))
				   {
				   model_num_temp = 0;
				   MPNL_CommandCheck_GetNum(char_pointer, "Model", model_num_temp);
				   Show_Statement("PDB file ", old_file_pointer, " being converted.");
				   Convert_PDB_To_Composite(old_file_pointer, new_file_pointer, 
					                        new_atoms, model_num_temp, name_preservation);
				   }
				else if (Name_Check(temp_old_ext, PDB_ext) && Name_Check(temp_new_ext, cryst_ext))
				   {
				   model_num_temp = 0;
				   MPNL_CommandCheck_GetNum(char_pointer, "Model", model_num_temp);
				   Show_Statement("PDB file ", old_file_pointer, " being converted.");
				   Convert_PDB_To_Structure(old_file_pointer, new_file_pointer, 
					                        new_atoms, model_num_temp, name_preservation);
				   }
				else
				   {
				   Show_Warning("FILE EXTENSIONS SPECIFIED FOR CONVERSION ARE NOT FAMILIAR!");
				   }

                                 
                }
             }
             
         //Crystal structure search---
          
         for (int index = 0; index < file_size; ++index)
             {
			 if (struct_index == MAX_CRYST)
				{
				Show_Warning("YOU HAVE TOO MANY CRYSTAL STRUCTURES DEFINED!");
				break;
				}

             char_pointer = &(file_text[index]);
             if (Command_Check(char_pointer, crystal_signature))
               //Found a "Define Structure" command in the file input, which is 
               //followed by crystal structure information.
                { 
				MPSNL_CommandCheck_SC(char_pointer, struct_names[struct_index], "Name");
                 //Get structure label/name.
				MPSNL_CommandCheck_SC(char_pointer, temp_type_pointer, "Struct");
                 //Get structure type (e.g. BCC, ORTHO, etc.) 
                is_special_system = (temp_type_pointer[0] == 'F') || 
                                    (temp_type_pointer[0] == 'B') ||
                                    (temp_type_pointer[0] == 'S');   
				 //Cubic systems (FCC, BCC, SC) are handled differently.
                num_lattice_parameters_set = 0;
				for (int a = 0; a < 9; ++a)
                 //Get lattice parameters.
                    {
					if (MPNL_CommandCheck(char_pointer, latticepam_signature)
						|| MPNL_CommandCheck(char_pointer, meshpam_signature))
					   {
					   lattice_parameters[a] = Get_FP_Number(char_pointer);
						//Get lattice parameter.
					   ++num_lattice_parameters_set;
					   }                             
					else 
					  //If no parameter add command is found, it is time to look for
					  //basis atom additions.
                       {
                       break;                  
                       }
                    }
                
                if (!is_special_system)
                 //If not dealing with special systems (which are handled differently
                 //in user input), define the crystal vectors now.
                   {
                   if ( (temp_type_pointer[0] == 'C') && (temp_type_pointer[1] == 'U'))
                    //CUBIC
                     {
                     new_crystal_systems[struct_index].Set_Cubic_Vectors(lattice_parameters[0]);          
                     }
                   else if ( (temp_type_pointer[0] == 'T') && (temp_type_pointer[1] == 'E') )
                    //TETRAGONAL
                     {
                     new_crystal_systems[struct_index].Set_Tetragonal_Vectors(lattice_parameters[0], 
                                                       lattice_parameters[1]);          
                     }
                   else if ( (temp_type_pointer[0] == 'O') && (temp_type_pointer[1] == 'R') )
                    //ORTHORHOMBIC
                     {
                     new_crystal_systems[struct_index].Set_Orthorhombic_Vectors(lattice_parameters[0], 
                                                       lattice_parameters[1], lattice_parameters[2]);          
                     }   
                   else if ( (temp_type_pointer[0] == 'H') )
                    //HEXAGONAL.
                     {
                     new_crystal_systems[struct_index].Set_Hexagonal_Vectors(lattice_parameters[0],
                                                       lattice_parameters[1]);     
                     }
                   else if ( (temp_type_pointer[0] == 'T') && (temp_type_pointer[1] == 'R'))
                    //TRICLINIC = General starting point of having to define everything.
                    //No symmetry is assumed by this program, but this can be
                    //used to define monoclinic and rhombohedral systems.
                     {
                     new_crystal_systems[struct_index].Set_Vectors(lattice_parameters);          
                     }
                   else if ( (temp_type_pointer[0] == 'A') )
                    //AMORPHOUS. Just set two parameters (grid bond length and variation parameter).
                     {
				     if (num_lattice_parameters_set == 1)
					  //If no variation parameter was specified for amorphous systems, use default.
					   {
					   lattice_parameters[1] = DEFAULT_MESH_VARI;
					   }
                     new_crystal_systems[struct_index].Set_Amor_Vectors(lattice_parameters[0], lattice_parameters[1]);      
                     }
				   else if ( (temp_type_pointer[0] == 'C') && (temp_type_pointer[1] == 'E') )
				    //CELL. Define crystal system with three vector magnitudes and three vector angles,
					//using standard conventions.
				     {
					 new_crystal_systems[struct_index].Set_Vectors(lattice_parameters, &(lattice_parameters[3]));
				     }
                   else
                     {
					 Show_Warning("FOR CRYSTAL SYSTEMS, SPECIFY CUBIC, TETRA, ORTHO, HEX, TRICLINIC, CELL, or AMOR!");                                                                      
                     }
				   }
                while (MPSNL_CommandCheck_GetFPNums(char_pointer, name_pointer, 
					                                crystal_atom_signature, 5, temp_coors, 3))
                 //Get atoms of the crystal basis.
                    {
                    rel_coor = (char_pointer[4] == 'R');
                    if (rel_coor)
                       //"Atom" command for absolute coordinates.
                       //"AtomR" command for relative coordinates.
                       //Hence 'R' character must be skipped
                         {
                         name_pointer = &(name_pointer[1]);          
                         }
                      temp_index = new_atoms.Find_Atom_Name(name_pointer, FALSEV);
                       //Identify the atom requested.
                      if (temp_index == -1)
                         {
						 Show_Warning("ATOM ", name_pointer, " IN CRYSTAL NOT DEFINED!");
                         continue;
                         }
                      temp_basis_atom = new_atoms[temp_index];
                      if (!is_special_system)
                       //FCC, BCC, and SC = Easy to set-up systems for the script
                       //writer where coordinates need not be specified.  
                         {
                         if (rel_coor)
                            {
                            new_crystal_systems[struct_index].Add_Atom_To_Basis_Rel
                            (temp_basis_atom, temp_coors);
                            }
                         else
                            {
                            new_crystal_systems[struct_index].Add_Atom_To_Basis
                            (temp_basis_atom, temp_coors);                  
                            }
                         }               
                    }
                
             while (MPNL_CommandCheck_GetNums(char_pointer, facet_signature, temp_facet_indices, 3))
                   {
                   if (!Is_Text(char_pointer[6]))
				    //Generic specification of faceting planes.
                      {
					  MPNL_CommandCheck_GetFPNum(char_pointer, "Scale", temp_scaling_factor, 1.0);
                      new_crystal_systems[struct_index].Set_Faceting_Plane(temp_facet_indices[0], 
                                                                           temp_facet_indices[1], 
                                                                           temp_facet_indices[2],
                                                                           temp_scaling_factor);
                      }
                   else
				    //Specific set of faceting planes.
					  {
					  temp_facet_name_pointer = &(char_pointer[6]);
					  MPNL_CommandCheck_GetFPNum(char_pointer, "Scale", temp_scaling_factor, 1.0);
					  Set_Faceting_Planes(new_crystal_systems[struct_index], temp_facet_name_pointer, 
										  temp_scaling_factor); 
					  }				   
                   }
                if (is_special_system)
                 //Final definition of the crystal structure for special systems.
                   {
                   if (temp_type_pointer[0] == 'F')
                    //FCC.
                     {
                     new_crystal_systems[struct_index].Set_FCC_System(temp_basis_atom, lattice_parameters[0]);        
                     }
                   else if (temp_type_pointer[0] == 'B')
                    //BCC.
                     {
                     new_crystal_systems[struct_index].Set_BCC_System(temp_basis_atom, lattice_parameters[0]);        
                     }
                   else if (temp_type_pointer[0] == 'S')
                    //SC.
                     {
                     new_crystal_systems[struct_index].Set_SC_System(temp_basis_atom, lattice_parameters[0]);        
                     }
                   }
                ++struct_index;
                }

		     //Crystal structure loading---
             
             if (Command_Check(char_pointer, crystal_load_signature))
                {
				MPSNL_CommandCheck_SC(char_pointer, struct_names[struct_index], "Name");
				MPNL_CommandCheck_SC(char_pointer, "File");
                new_crystal_systems[struct_index].Load_Crystal_System(char_pointer);
                ++struct_index;
                }
             
             }

		 if (struct_index != 0)
			{
			Show_Statement("Done defining ", struct_index, " structures!");
			}

         //Supercell search---
 
         for (int index = 0; index < file_size; ++index)
             {
             char_pointer = &(file_text[index]);
             if (Command_Check(char_pointer, SUPERCELL_signature))
                {
				MPSNL_CommandCheck_SC(char_pointer, temp_super_name, "Name");
				MPSNL_CommandCheck_SC(char_pointer, temp_cell_struct, "Struct");
                temp_super_index = Find_String_Match(temp_cell_struct, struct_names, struct_index);
				MPNL_CommandCheck_GetNums(char_pointer, "Size", temp_super_sizes, 3);
				strcpy(temp_super_string, temp_super_name);
                Make_File_Name(temp_super_string, super_ext);
                new_crystal_systems[temp_super_index].Save_SuperCell(temp_super_string, temp_super_sizes[0], 
                                                                     temp_super_sizes[1], temp_super_sizes[2]);                        
                }
             }
         

		 //Molecule search---

         for (int index = 0; index < file_size; ++index)
             {
			if (molecule_index == MAX_MOLEC)
				{
				Show_Warning("YOU HAVE TOO MANY MOLECULES DEFINED!");
				break;
				}

             char_pointer = &(file_text[index]);
             if (Command_Check(char_pointer, molecule_signature))
				{
				MPSNL_CommandCheck_SC(char_pointer, molecule_names[molecule_index], "Name");
                 //Get molecule label/name.
                while (MPSNL_CommandCheck_GetFPNums(char_pointer, temp_tomol_name, "Add", 4, temp_coors, 3))		
                 //Get first atom/atomic group specified in the file by name and position.
                      {
                      temp_atom_index = new_atoms.Find_Atom_Name(temp_tomol_name, FALSEV);            
                       //Get name of each atom/atomic group. 
                      if (temp_atom_index != -1)
                       //If true, "Add Atom" operation has been requested.
                           {    
                           new_molecules[molecule_index].Add_Atom(new_atoms[temp_atom_index], temp_coors);
                           }   
                      else if (Command_Check(temp_tomol_name, water_signature))
                       //If true, "Add Water" operation has been requested.
                           { 
						   MPNL_CommandCheck_GetFPNums(char_pointer, "Pos", temp_coors, 3);
						   MPNL_CommandCheck_GetNames(char_pointer, "Atoms", 5, temp_atom_name_list, 10);
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 2);
                           new_molecules[molecule_index].Place_Water(temp_coors, 
                                                         new_atoms[temp_atom_indices[0]],
                                                         new_atoms[temp_atom_indices[1]]);
                           }
                      else if (Command_Check(temp_tomol_name, alkanethiol_signature))
                           {
						   MPNL_CommandCheck_GetNum(char_pointer, "Methyl Count", temp_carbon_count);
                           MPNL_CommandCheck_GetNames(char_pointer, "Atoms", 5, temp_atom_name_list, 10);
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 4);
                           new_molecules[molecule_index].Build_Alkanethiol(temp_carbon_count, 
                                                         new_atoms[temp_atom_indices[0]], 
                                                         new_atoms[temp_atom_indices[1]], 
                                                         new_atoms[temp_atom_indices[2]],
														 new_atoms[temp_atom_indices[3]]);
                           } 
                      else if (Command_Check(temp_tomol_name, carbonchain_signature))
                           {
                           MPNL_CommandCheck_GetNum(char_pointer, "Methyl Count", temp_carbon_count);
                           MoleculeAdd_Reader(char_pointer, temp_previous_end, temp_start, temp_atom_name_list);
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 2);
                           new_molecules[molecule_index].Add_Carbon_Chain(temp_previous_end, 
                                                         temp_start, temp_coor_spaceA, temp_coor_spaceB, 
                                                         temp_carbon_count, new_atoms[temp_atom_indices[0]], 
                                                         new_atoms[temp_atom_indices[1]]);              
                           } 
                      else if (Command_Check(temp_tomol_name, OH_signature))
                           {
                           MoleculeAdd_Reader(char_pointer, temp_previous_end, temp_start, temp_atom_name_list); 
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 2); 
                           new_molecules[molecule_index].Add_OH(temp_previous_end, temp_start, 
                                                         new_atoms[temp_atom_indices[0]], 
                                                         new_atoms[temp_atom_indices[1]]);                                 
                           } 
                      else if (Command_Check(temp_tomol_name, SH_signature))
                           {
                           MoleculeAdd_Reader(char_pointer, temp_previous_end, temp_start, temp_atom_name_list); 
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 2); 
                           new_molecules[molecule_index].Add_SH(temp_previous_end, temp_start, 
                                                         new_atoms[temp_atom_indices[0]], 
                                                         new_atoms[temp_atom_indices[1]]);  
                           } 
                      else if (Command_Check(temp_tomol_name, CH2_signature))
                           {
                           MoleculeAdd_Reader(char_pointer, temp_previous_end, temp_start, temp_atom_name_list); 
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 2); 
                           new_molecules[molecule_index].Add_CH2(temp_previous_end, temp_start, 
                                                         new_atoms[temp_atom_indices[0]], 
                                                         new_atoms[temp_atom_indices[1]]);                                   
                           } 
                      else if (Command_Check(temp_tomol_name, CH3_signature))
                           {
                           MoleculeAdd_Reader(char_pointer, temp_previous_end, temp_start, temp_atom_name_list); 
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 2); 
                           new_molecules[molecule_index].Add_CH3(temp_previous_end, temp_start, 
                                                         new_atoms[temp_atom_indices[0]], 
                                                         new_atoms[temp_atom_indices[1]]);                                   
                           } 
                      else if (Command_Check(temp_tomol_name, AX_signature))
                           {
                           MoleculeAdd_Reader(char_pointer, temp_previous_end, temp_start, temp_atom_name_list);      
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 2); 
						   MPNL_CommandCheck_SC(char_pointer, "Shape");
                           Convert_ShapeFactor(char_pointer[0], temp_shape_factor);     
						   MPNL_CommandCheck_GetFPNum(char_pointer, "Bond Length", temp_bond_lengths[0]);                      
                           new_molecules[molecule_index].Add_AX(temp_previous_end, temp_start, 
                                                                new_atoms[temp_atom_indices[0]], 
                                                                new_atoms[temp_atom_indices[1]], 
                                                                temp_bond_lengths[0], temp_shape_factor);                    
                           } 
                      else if (Command_Check(temp_tomol_name, AXY_signature))
                           {
                           MoleculeAdd_Reader(char_pointer, temp_previous_end, temp_start, temp_atom_name_list);      
                           new_atoms.Get_Index_List(temp_atom_name_list, temp_atom_indices, 3);     
                           MPNL_CommandCheck_SC(char_pointer, "Shape");
                           Convert_ShapeFactor(char_pointer[0], temp_shape_factor);      
                           MPNL_CommandCheck_GetFPNum(char_pointer, "Bond Length", temp_bond_lengths[0]); 
						   MPNL_CommandCheck_GetFPNum(char_pointer, "Bond Length", temp_bond_lengths[1]);    
                           new_molecules[molecule_index].Add_AXY(temp_previous_end, temp_start, 
                                                                 new_atoms[temp_atom_indices[0]], 
                                                                 new_atoms[temp_atom_indices[1]], 
                                                                 new_atoms[temp_atom_indices[2]],
                                                                 temp_bond_lengths[0],
                                                                 temp_bond_lengths[1],  temp_shape_factor);                                     
                           } 
                      else if (Command_Check(temp_tomol_name, AXYZ_signature))
                           {
                           MoleculeAdd_Reader(char_pointer, temp_previous_end, temp_start, temp_atom_name_list);      
                           new_atoms.Get_Index_List(temp_atom_name_list, 
                                                    temp_atom_indices, 4);     
                           MPNL_CommandCheck_SC(char_pointer, "Shape");
                           Convert_ShapeFactor(char_pointer[0], temp_shape_factor);   
                           for (int a = 0; a < 3; ++a)
                               {   
                               MPNL_CommandCheck_GetFPNum(char_pointer, "Bond Length", temp_bond_lengths[a]); 
                               }                                  
                           new_molecules[molecule_index].Add_AXYZ(temp_previous_end, temp_start, 
                                                                  new_atoms[temp_atom_indices[0]], 
                                                                  new_atoms[temp_atom_indices[1]], 
                                                                  new_atoms[temp_atom_indices[2]],
                                                                  new_atoms[temp_atom_indices[3]],
                                                                  temp_bond_lengths[0],
                                                                  temp_bond_lengths[1],
                                                                  temp_bond_lengths[2],
                                                                  temp_shape_factor);                                   
                           } 
                      else
                           {
						   Show_Warning("AN ATOM IN ", molecule_names[molecule_index], " NEEDS TO BE DEFINED!");
                           }
                      }
                ++molecule_index;
                }

			 //Molecule loading---
             
             if (Command_Check(char_pointer, molecule_load_signature))
                {
				MPSNL_CommandCheck_SC(char_pointer, molecule_names[molecule_index], "Name");
				MPNL_CommandCheck_SC(char_pointer, "File");
                new_molecules[molecule_index].Load_Molecule(char_pointer);      
                ++molecule_index;
                }
                
             }

		 if (molecule_index != 0)
			{
			Show_Statement("Done defining ", molecule_index, " molecules!");
			}

         //Nanosurface definitions---

         for (int index = 0; index < file_size; ++index)
             {        
			 if (surface_index == MAX_NANOS)
				{
				Show_Warning("YOU HAVE DEFINED TOO MANY NANOSURFACES!");
				break;
				}
             
			 char_pointer = &(file_text[index]); 
             if (Command_Check(char_pointer, nanosurface_signature))
                {
				MPSNL_CommandCheck_SC(char_pointer, nanosurface_names[surface_index], "Name");
                 //Get surface label/name.
				MPNL_CommandCheck_SC(char_pointer, "Struct");
                 //Get surface structural description.
                temp_structure_index = Find_String_Match(char_pointer, struct_names, struct_index);
                 //Determine the structure being requested.
                if (temp_structure_index == -1)
                   {
				   Show_Warning("STRUCTURE ", char_pointer, " NOT DEFINED!");      
                   }          
                for (int a = 0; a < 3; ++a)
                 //Get nanosurface size dimensions.
                  {  
				  MPNL_CommandCheck_GetFPNum(char_pointer, "Size", temp_sizes[a]);
                  }                       
                for (int a = 0; a < 3; ++a)
                 //Get nanosurface terminating face indices.
                  {
				  MPNL_CommandCheck_GetNum(char_pointer, "FaceIndex", temp_indices[a]);
                  }  
				MPNL_CommandCheck_GetFPNum(char_pointer, "Rotate", temp_angle);
				temp_angle *= PI_CONST/180.0;
				MPNL_CommandCheck_SC(char_pointer, "BasisAdd");
				vary_size_for_neutrality = FALSEV;
				complete_neutralization_needed = FALSEV;
				stuffed_surface = TRUEV;
				if (char_pointer[0] == 'P')
				   { 
				   if (MPNL_CommandCheck_SC(char_pointer, "VariableSize"))
					  {
					  if (char_pointer[0] == 'Y')
					   //Surface creation will systematically vary size as to minimize
					   //the total charge of the slab.
						{
						vary_size_for_neutrality = TRUEV;
						}
					  }
				   if (MPNL_CommandCheck_SC(char_pointer, "SurfaceNeutralize"))
					  {
					  if (char_pointer[0] == 'Y')
					   //Surface creation will remove atoms from top surface as to 
					   //bring the charge to zero.
						{
						complete_neutralization_needed = TRUEV;
						}
					  }
				   }
				else if (char_pointer[0] == 'F')
				   {
                   stuffed_surface = FALSEV;
				   }

				staircase_shift_request = FALSEV;
				if (MPNL_CommandCheck_GetFPNums(char_pointer, "Staircase", temp_shift_parameters, 2))
				   {
				   staircase_shift_request = TRUEV;
				   }

				randomize_request = FALSEV;
				if (MPNL_CommandCheck_GetFPNums(char_pointer, "Randomize", temp_randomize_parameters, 3))
				   {
				   randomize_request = TRUEV;
				   }

				sliceslab_request = FALSEV;
				if (MPNL_CommandCheck_GetFPNums(char_pointer, "Slice", temp_sliceslab_parameters, 2))
				   {
				   sliceslab_request = TRUEV;
				   } 
				
				MPSNL_CommandCheck_GetFPNums(char_pointer, temp_s_coat_name, "Coat", 2, temp_s_coat_parameters, 3);
                temp_s_coat_index = Find_String_Match(temp_s_coat_name, molecule_names, molecule_index);
				MPNL_CommandCheck_SC(char_pointer, "Mode");
                if (char_pointer[0] == 'B')
                   {
                   surface_coat_type = BIND;                 
                   }
                else
                   {
                   surface_coat_type = FREE;                  
                   } 
                new_nanosurfaces[surface_index].Set_Properties(new_crystal_systems[temp_structure_index], 
                                                               temp_sizes[0], temp_sizes[1], temp_sizes[2], 
                                                               temp_indices[0], temp_indices[1], temp_indices[2],
                                                               temp_angle, stuffed_surface, vary_size_for_neutrality, 
															   complete_neutralization_needed);
				Show_Statement("Starting construction of surface ", nanosurface_names[surface_index]);
				new_nanosurfaces[surface_index].Construct_Surface();
				if (staircase_shift_request)
				   {
				   Show_Statement("Performing staircase shift on nanosurface.");
				   new_nanosurfaces[surface_index].Staircase_Shift(temp_shift_parameters[0], 
					                               temp_shift_parameters[1]);
				   }
				if (randomize_request)
				   {
				   Show_Statement("Performing randomization of nanosurface.");
				   new_nanosurfaces[surface_index].Random_Shift(temp_randomize_parameters);
				   }
				if (sliceslab_request)
				   {   
				   Show_Statement("Performing slicing of nanosurface.");
				   new_nanosurfaces[surface_index].Slice(temp_sliceslab_parameters[0], 
														 temp_sliceslab_parameters[1]);
				   }

				if (temp_s_coat_index != -1)     
                   { 
                   new_nanosurfaces[surface_index].Coat_Surface(new_molecules[temp_s_coat_index], 
                                                                temp_s_coat_parameters[0],
                                                                temp_s_coat_parameters[1],
                                                                temp_s_coat_parameters[2],
                                                                surface_coat_type);                                                     
                   }      
                ++surface_index; 
                }
                
             //Nanosurface loading---

             if (Command_Check(char_pointer, nanosurface_load_signature))
                {
				MPSNL_CommandCheck_SC(char_pointer, nanosurface_names[surface_index], "Name");
                 //Get surface label/name.
				MPNL_CommandCheck_SC(char_pointer, "File");
				Show_Statement("Loading nanosurface from: ", char_pointer);
                new_nanosurfaces[surface_index].Load_NanoSurface(char_pointer);  
				if (MPNL_CommandCheck_GetFPNums(char_pointer, "Staircase", temp_shift_parameters, 2))
				   {
				   new_nanosurfaces[surface_index].Staircase_Shift(temp_shift_parameters[0], 
					                               temp_shift_parameters[1]);
				   }
				if (MPNL_CommandCheck_GetFPNums(char_pointer, "Randomize", temp_randomize_parameters, 3))
				   {
				   new_nanosurfaces[surface_index].Random_Shift(temp_randomize_parameters);
				   }
				if (MPNL_CommandCheck_GetFPNums(char_pointer, "Slice", temp_sliceslab_parameters, 2))
				   {
				   new_nanosurfaces[surface_index].Slice(temp_sliceslab_parameters[0], 
														 temp_sliceslab_parameters[1]);
				   }
                ++surface_index;
                }
             }

		 if (surface_index != 0)
			{
			Show_Statement("Done defining ", surface_index, " surfaces!");
			}

         //Nanoparticle definitions---  

         for (int index = 0; index < file_size; ++index)
             {
			 if (particle_index == MAX_NANOP)
				{
				Show_Warning("YOU HAVE TOO MANY NANOPARTICLES DEFINED!");
				break;
				}
			 
             char_pointer = &(file_text[index]);
             if (Command_Check(char_pointer, nanoparticle_signature))
                {
				MPSNL_CommandCheck_SC(char_pointer, nanoparticle_names[particle_index], "Name");
                 //Get particle label/name.
				MPNL_CommandCheck_SC(char_pointer, "Struct");
                 //Get particle structural description.
                temp_p_structure_index = Find_String_Match(char_pointer, struct_names, struct_index);
                 //Determine which structure is being referenced.
                if (temp_p_structure_index == -1)
                   {
				   Show_Warning("STRUCTURE ", char_pointer, " NOT DEFINED!");                        
                   }
				MPNL_CommandCheck_GetFPNum(char_pointer, "Radius", temp_diameter);
                 //User input is the radius, so the value must be doubled.
				temp_diameter *= 2.0;
				temp_diameter += FP_ERROR_FIX;
				 //Make sure user can say "I want exactly three lattice constant-width 
				 //based nanoparticle creation" with the outlying atoms not accidentally
				 //cut by the particle creation algorithm.
                 //Get particle diameter.
				MPNL_CommandCheck_GetFPNum(char_pointer, "Minimum", temp_hole_diameter);
				temp_hole_diameter *= 2.0;
				 //Diameter of inner hole within the nanoparticle.
				MPNL_CommandCheck_SC(char_pointer, "Model");
                 //Get nanoparticle model (continuous vs. crystalline) type.
                if (char_pointer[1] == 'O')
                   {
                   temp_model = CON_MODEL;                 
                   }
                else if (char_pointer[1] == 'R')
                   {
                   temp_model = CRY_MODEL;
                   }
                else if (char_pointer[1] == 'P')
                   {
                   temp_model = SPH_MODEL;                      
                   }
				else
				   {
				   Show_Warning("PARTICLE MODEL NOT SPECIFIED CORRECTLY!");
				   }
				MPNL_CommandCheck_SC(char_pointer, "BasisAdd");
				neutralize_by_surface_removal = FALSEV;
				stuffed_particle = TRUEV;
				if (char_pointer[0] == 'P')
				   {
				   if (MPNL_CommandCheck_SC(char_pointer, "NeutralizeParticle"))
					  {
					  if (char_pointer[0] == 'Y')
						{
						neutralize_by_surface_removal = TRUEV;
						}
					  }
				   }
				else if (char_pointer[0] == 'F')
				   {
                   stuffed_particle = FALSEV;
				   }

		        Zero_XYZ(temp_angles);
				MPNL_CommandCheck_GetFPNums(char_pointer, "Rotate", temp_angles, 3);
                 //Phi-theta-psi (x-y-z rotation).
				Divide_XYZ(temp_angles, 180.0/PI_CONST);
				MPSNL_CommandCheck_GetFPNums(char_pointer, temp_p_coat_name, "Coat", 5, 
					                         temp_p_coat_parameters, 2);
                temp_p_coat_index = Find_String_Match(temp_p_coat_name, molecule_names, molecule_index);   
				MPNL_CommandCheck_SC(char_pointer, "Mode");
                if (char_pointer[0] == 'B')
                   {
                   particle_coat_type = BIND;              
                   }
                else
                   {
                   particle_coat_type = FREE;                  
                   }
                new_nanoparticles[particle_index].Set_Properties(temp_diameter, temp_hole_diameter,
												  temp_model, new_crystal_systems[temp_p_structure_index], 
                                                  stuffed_particle, neutralize_by_surface_removal);  
				Show_Statement("Building nanoparticle ", nanoparticle_names[particle_index]);
				new_nanoparticles[particle_index].Construct_Particle();
				if (temp_p_coat_index != -1)
                   {
                   if (temp_model != CRY_MODEL)
                      {
                      new_nanoparticles[particle_index].Coat_Spherical_Particle(
                                                        new_molecules[temp_p_coat_index],
                                                        temp_p_coat_parameters[0],
                                                        temp_p_coat_parameters[1],
                                                        particle_coat_type); 
                      }  
                   else
                      {
                      new_nanoparticles[particle_index].Coat_Faceted_Particle(
                                                        new_molecules[temp_p_coat_index],
                                                        temp_p_coat_parameters[0],
                                                        temp_p_coat_parameters[1],
                                                        particle_coat_type); 
                      }                                                   
                   }          
				new_nanoparticles[particle_index].Rotate(temp_angles);
                ++particle_index; 
                }

			 //Nanoparticle loading---

             if (Command_Check(char_pointer, nanoparticle_load_signature))
                {
				MPSNL_CommandCheck_SC(char_pointer, nanoparticle_names[particle_index], "Name");
				MPNL_CommandCheck_SC(char_pointer, "File");
                new_nanoparticles[particle_index].Load_NanoParticle(char_pointer);      
                ++particle_index;
                }
             
             }

		 if (particle_index != 0)
			{
			Show_Statement("Done defining ", particle_index, " particles!");
			}

         //General composite definition---           

		 first_region_added[0] = 0;
		 first_region_added[1] = 0;
		  //Make a record of the first group of atoms added to the composite, according to
		  //their composite indices.

		 //Loading requests for systems to add to the composite are checked first here---

		 temp_mode = FALSEV;
		 int start_index, end_index;
         for (int index = 0; index < file_size; ++index)
             {   
             char_pointer = &(file_text[index]);
             
             if (Command_Check(char_pointer, gencomp_load_signature))
              //Check for loading first so that user can build composites
              //off of other stored composites.
                {
				MPSNL_CommandCheck_SC(char_pointer, composite_name, "Name");
                 //Get composite label/name.
				MPSNL_CommandCheck_SC(char_pointer, temp_file_pointer, "File");
                 //Get file name.
				MPNL_CommandCheck_GetFPNums(char_pointer, "Pos", temp_positions, 3);
				temp_mode = TRUEV;
				 //Center addition.
				if (MPNL_CommandCheck_SC(char_pointer, "Mode"))
				   {
				   if (char_pointer[0] == 'M')
				    //Minimum coordinate addition.
				     {
				     temp_mode = FALSEV;
				     }
				   }
				Zero_XYZ(temp_angles);
				MPNL_CommandCheck_GetFPNums(char_pointer, "Rotate", temp_angles, 3);
                 //Phi-theta-psi (x-y-z rotation).
				Divide_XYZ(temp_angles, 180.0/PI_CONST);
				start_index = new_composite.Get_Atom_Count();
                new_composite.Load_Composite(temp_file_pointer, temp_positions, 
					                         temp_angles, temp_mode);
				end_index = new_composite.Get_Atom_Count();
				if (first_region_added[1] == 0)
				 //Record rigid minimization bounds.
				    {
				    first_region_added[1] = end_index;
				    }
				if (MPNL_CommandCheck_GetFPNums(char_pointer, coordex_signature, 
					                            coordex_pams, 2))
				 //Coordination number addition to atomic names.
				    {
					new_composite.Coordex_Me(start_index, end_index, coordex_pams[0], coordex_pams[1]);
				    }
                }

			 start_index = new_composite.Get_Atom_Count();
             if (Command_Check(char_pointer, gencomp_signature))    
                 {
				 MPSNL_CommandCheck_SC(char_pointer, composite_name, "Name");
                 //Get composite label/name.
				 MPNL_CommandCheck_GetFPNums(char_pointer, "Size", temp_comp_sizes, 3);    
                 new_composite.Resize(temp_comp_sizes);
                  //Get composite box size parameters.
				 MPNL_CommandCheck_GetPosNum(char_pointer, "Period", temp_periodicity);
                 new_composite.Set_Periodicity(temp_periodicity);
                  //Get the periodicity of the composite box.
                 while (MPSNL_CommandCheck_GetFPNums(char_pointer, temp_component_name, add_signature, 4, temp_coors, 3)
					   || MPSNL_CommandCheck_GetFPNums(char_pointer, temp_component_name, sub_signature, 4, temp_coors, 3) 
					   || MPSNL_CommandCheck_GetFPNums(char_pointer, temp_component_name, invsub_signature, 4, temp_coors, 3)
					   || MPSNL_CommandCheck_GetFPNums(char_pointer, temp_component_name, bulldoze_signature, 4, temp_coors, 3))    
                  //While there are addition/substitution commands to carry out, do exactly that.  
                      {                   
					  if (char_pointer[0] == 'S')
					   //Substitution command.
					     {
						 new_composite.Temp_Set_Sub_Mode();
					     }
					  else if (char_pointer[0] == 'B')
					   //Bulldoze command.
					     {
						 new_composite.Temp_Set_Bulldoze_Mode();
					     }
					  else if (char_pointer[0] =='I')
					    //Inverse substitution command.
					     {
						 new_composite.Temp_Set_InvSub_Mode();
					     }
                      search_index = Find_String_Match(temp_component_name, nanoparticle_names, particle_index);
                       //Nanoparticle add check.
                      if (search_index != -1)
                         {
						 Show_Statement("Adding nanoparticle ", temp_component_name, " to composite.");
						 Zero_XYZ(temp_angles);
						 MPNL_CommandCheck_GetFPNums(char_pointer, "Rotate", temp_angles, 3);
                             //Phi-theta-psi (x-y-z rotation). 
					     Divide_XYZ(temp_angles, 180.0/PI_CONST);
                         new_composite.Add_NanoParticle(new_nanoparticles[search_index], 
							                            temp_coors, temp_angles, TRUEV); 
                         }
                      
                      if (search_index == -1)
                       //Nanosurface add check.
                         {
                         search_index = Find_String_Match(temp_component_name, nanosurface_names, surface_index);   
                         if (search_index != -1)
                            {
							Show_Statement("Adding nanosurface ", temp_component_name, " to composite.");
						    Zero_XYZ(temp_angles);
							MPNL_CommandCheck_GetFPNums(char_pointer, "Rotate", temp_angles, 3);
                             //Phi-theta-psi (x-y-z rotation).
						    Divide_XYZ(temp_angles, 180.0/PI_CONST);
                            new_composite.Add_NanoSurface(new_nanosurfaces[search_index], 
								                          temp_coors, temp_angles, TRUEV, super_hydrogen_remove); 
                            }           
                         }         
                                
                      if (search_index == -1)
                       //Molecule add check.
                         {
                         search_index = Find_String_Match(temp_component_name, molecule_names, molecule_index);
                         if (search_index != -1)
                            {
							Show_Statement("Adding molecule ", temp_component_name, " to composite.");
							MPNL_CommandCheck_SC(char_pointer, "Mode");  
                            if (char_pointer[0] == 'C')
                               {
                               temp_add_mode = MOLECULE_CENTER;                   
                               }
                            else
                               {
                               temp_add_mode = FIRST_ATOM;                   
                               } 
							Zero_XYZ(temp_angles);
							MPNL_CommandCheck_GetFPNums(char_pointer, "Rotate", temp_angles, 3);
                             //Phi-theta-psi (x-y-z rotation).			
					        Divide_XYZ(temp_angles, 180.0/PI_CONST);
                            new_composite.Add_Molecule(new_molecules[search_index], 
                                          temp_coors, temp_angles, 
                                          TRUEV, temp_add_mode);     
                            }                  
                         }
                      
                      if (search_index == -1)
                       //Check for supercell addition.
                         {  
                         search_index = Find_String_Match(temp_component_name, struct_names, struct_index);   
                         if (search_index != -1)
                            {
							Show_Statement("Adding supercell ", temp_component_name, " to composite.");
							MPNL_CommandCheck_GetNums(char_pointer, "Size", temp_supercell_sizes, 3);
                            new_composite.Add_SuperCell(new_crystal_systems[search_index], temp_coors, 
                                                        temp_supercell_sizes[0], temp_supercell_sizes[1],
                                                        temp_supercell_sizes[2], TRUEV);   
                            }       
                         }
                         
                      if (search_index == -1)
                       //Check for single atom addition.  
                         {
                         search_index = new_atoms.Find_Atom_Name(temp_component_name, FALSEV);             
                         if (search_index != -1)
                           {    
						   Show_Statement("Adding atom ", temp_component_name, " to composite.");
                           new_composite.Add_Atom(new_atoms[search_index], temp_coors);
                           }             
                         }
                         
                      if (search_index == -1)
                         {
						 Show_Warning("COMPONENT ", temp_component_name, " IS NOT DEFINED!");
                         }    

					  if (first_region_added[1] == 0)
							{
							first_region_added[1] = new_composite.Get_Atom_Count();
							}
                      } 
                 
				 if (MPNL_CommandCheck_GetFPNums(char_pointer, "Void", temp_void_origin, 3))
				     {
					 MPNL_CommandCheck_GetFPNums(char_pointer, "Size", temp_void_size, 3);
					 Show_Statement("Creating the requested void in the composite.");
					 new_composite.Make_Void(temp_void_origin, temp_void_size, TRUEV);
				     }

				 while (MPSNL_CommandCheck_SC(char_pointer, temp_component_name, chargebal_signature))
				    {
					MPNL_CommandCheck_GetFPNums(char_pointer, "Box Position", temp_chargebal_box_origin, 3);
					MPNL_CommandCheck_GetFPNums(char_pointer, "Box Size", temp_chargebal_box_size, 3);
                    search_index = Find_String_Match(temp_component_name, molecule_names, molecule_index);
					if (search_index == -1)
                     //Molecule not found. Look for atom usage. 
					   {
                       search_index = new_atoms.Find_Atom_Name(temp_component_name, FALSEV);
					   if (search_index != -1)
					     {
						 Zero_XYZ(temp_coors);
					     new_molecules[molecule_index].Add_Atom(new_atoms[search_index], temp_coors);
						 search_index = molecule_index;
						 ++molecule_index;
					     }
					   }
                    if (search_index != -1)
                       {
                       new_composite.Charge_Balance(new_molecules[search_index], 
                                     temp_chargebal_box_origin, temp_chargebal_box_size, TRUEV);               
                       }
                    else
                       {
					   Show_Statement("Note: Charge balance is being performed via ion removal.");
					   new_composite.Charge_Balance(temp_chargebal_box_origin, 
                                     temp_chargebal_box_size);  
                       } 
				    }  

                 while (MPSNL_CommandCheck_GetFPNums(char_pointer, temp_component_name, "Solvate", 8, temp_coors, 3))
                    {
                    search_index = Find_String_Match(temp_component_name, molecule_names, molecule_index);
					MPNL_CommandCheck_GetFPNums(char_pointer, "Box Position", temp_solvent_box_origin, 3);
					MPNL_CommandCheck_GetFPNums(char_pointer, "Box Size", temp_solvent_box_size, 3);
                    if (search_index != -1)
                       {
                       new_composite.FillSpace_With_Solvent(new_molecules[search_index], 
                                     temp_solvent_box_origin, temp_solvent_box_size, temp_coors, TRUEV);     
                       }
                    else
                       {
					   Show_Warning("SOLVATION MOLECULE NOT DEFINED!");        
                       }           
                    }
                 if (MPNL_CommandCheck_GetFPNums(char_pointer, coordex_signature, coordex_pams, 2))
				    {
					new_composite.Coordex_Me(start_index, new_composite.Get_Atom_Count(), 
											 coordex_pams[0], coordex_pams[1]);
				    }
				 }  
             }

		 Show_Statement("Note: Total charge of composite system is: ", new_composite.Get_Charge());

         //General commands search---

		 //Default Settings---

         image_hydrogen = TRUEV;
         MD_hydrogen = TRUEV;  
		 anal_hydrogen = FALSEV;
		 MD_mass = FALSEV;
		 cat_coreshell = FALSEV;
		 ani_coreshell = FALSEV;
		 no_rad_calc = FALSEV;
		 clean_up = FALSEV;
		 histo_pam = 0.005; 
		  //Default spacing in distance histogram.
		 radial_step = 0.001;
		  //Default integration step in RDF calculation.
		 strangeness_dist_min = 0.05;
		 strangeness_dist_max = 0.5;
		  //Default minimum and maximum distances of strangeness = 0.5A/5.0A.
		  //Every atom should have at least one atom within this distance maximum and
		  //no two atoms should be closer than the distance minimum.
		 cleanup_pam = 0;
		 xrd_wave = 0.15418;
		  //Cu K-alpha radiation. Naicker, P. K. et al. J. Phys. Chem. B. 109. 2005. 15243.
		 xrd_ionic = FALSEV;
		  //Don't use ionic scattering factors instead of neutral atom scattering factors.
		 xrd_precision = 0.0015;
		  //Default change in angle per XRD data point calculation, in degrees.
		 xrd_dw_factor = 0.0;
		  //Default DW factor in temperature-effects on XRD patterns of zero.
		  //Note: Equation used is: DW-attenutation of intensity = exp(-2*B*s*s).
		 xrd_min_angle = 3.0 * PI_CONST/180.0;
		 xrd_max_angle = 60.0 * PI_CONST/180.0;
		  //Range of angles considered in XRD_calculation, by default.
		 xrd_dw_factor_surface_ratio = 0.0;
		 xrd_coor_cutoff = 0.0;
		 xrd_coor_num = 0;
		 xrd_surface_effects = FALSEV;
		  //Default: No inclusion of the change in DW factors at an atomic surface.
		 mpam.nvt_run = FALSEV;
		 mpam.npt_run = FALSEV;
		  //NVT vs NPT ensemble for MD simulation.
		 mpam.molecule_commands = FALSEV;
		  //Inclusion of molecule commands in MD simulation.
		 spam.scope_mode = STEM_MODE;
		 spam.cell_x = spam.cell_y = spam.cell_z = 1;
		 spam.tilt_x = spam.tilt_y = spam.tilt_z = 0.0;
		 spam.scan_x_min = spam.scan_y_min = 0.0;
		 spam.scan_x_max = spam.scan_y_max = 6.0;
		 spam.raster_x = spam.raster_y = 60;
		 spam.window_x_size = spam.window_y_size = 2.0; 
		 spam.image_x_res = spam.image_y_res = 0.005;
		  //Microscope parameters: put in nm here, converted to Angstroms below.
		 spam.beam_tilt_x = spam.beam_tilt_y = 0.0;
		 spam.defocus = -4.3;
		  //Defocus of -4.3 nm is the default Scherzer defocus for a Cs3 of 0.005 mm
		  //and voltage of 200 keV.
		 spam.TDS_runs = 0;
		 spam.slice_output_count = 10;
		 image_dwf = FALSEV;
		  //Don't include Debye-Waller factors in image simulation files.

		 //High priority general commands---

         for (int index = 0; index < file_size; ++index)
          //Look for high-priority general commands.
             { 
             char_pointer = &(file_text[index]);
			 Command_Check_Set_BooleanF(char_pointer, NO_HYDROGEN_IMAGE_signature, 
				                        image_hydrogen);
			 Command_Check_Set_BooleanF(char_pointer, NO_HYDROGEN_MD_signature, 
				                        MD_hydrogen); 
             if (Command_Check(char_pointer, NO_HYDROGEN_signature))
                {
                image_hydrogen = FALSEV;
                MD_hydrogen = FALSEV;                             
                }
			 if (MD_hydrogen == FALSEV)
				{
				for (int a = 0; a < atom_index; ++a)
					{
					if (new_atoms[a].Is_Hydrogen())
					 //Prevent the user from making a terrible mistake with charge
					 //sums by removing any assigned hydrogen charge.
						{
						new_atoms[a].Set_Atom_Charge(0.00);
						}
					}
				}
			 Command_Check_Set_BooleanT(char_pointer, ANAL_HYDROGEN_signature, 
				                        anal_hydrogen);
			 Command_Check_Set_BooleanT(char_pointer, MD_MASS_signature, MD_mass);
			 Command_Check_Set_BooleanT(char_pointer, CAT_CORESHELL_signature,
				                        cat_coreshell);
			 Command_Check_Set_BooleanT(char_pointer, ANI_CORESHELL_signature,
				                        ani_coreshell);
			 if (Command_Check(char_pointer, CORESHELL_signature))
				{
				cat_coreshell = TRUEV;
				ani_coreshell = TRUEV;
				}
			 if (Command_Check(char_pointer, histogram_parameter_signature))
			    {
				histo_pam = Get_FP_Number(char_pointer);
			    }
			 if (Command_Check(char_pointer, radial_integration_parameter_signature))
				{
				radial_step = Get_FP_Number(char_pointer);
				}
			 if (Command_Check(char_pointer, strangeness_parameter_signature))
			    {
				MPNL_CommandCheck_GetFPNum(char_pointer, "Min", strangeness_dist_min, 
					                       0.05);
				MPNL_CommandCheck_GetFPNum(char_pointer, "Max", strangeness_dist_max, 
					                       0.5);
			    }
			 Command_Check_Set_BooleanT(char_pointer, NORAD_signature, no_rad_calc);
			 Command_Check_Set_BooleanT(char_pointer, CLEANUP_signature, clean_up);
			 if (Command_Check(char_pointer, cleanup_parameter_signature))
			    {
				cleanup_pam = Get_Number(char_pointer);
			    }
			 if (Command_Check(char_pointer, XRD_wave_signature))
			    {
				xrd_wave = Get_FP_Number(char_pointer);
			    }
			 if (Command_Check(char_pointer, XRD_min_signature))
			    {
			    xrd_min_angle = Get_FP_Number(char_pointer);
			    xrd_min_angle *= PI_CONST/180.0;
			    }
			 if (Command_Check(char_pointer, XRD_max_signature))
			    {
				xrd_max_angle = Get_FP_Number(char_pointer);
				xrd_max_angle *= PI_CONST/180.0;
			    }
			 if (Command_Check(char_pointer, XRD_precision_signature))
			    {
			    xrd_precision = Get_FP_Number(char_pointer);
				xrd_precision *= PI_CONST/180.0;
			    }
			 if (Command_Check(char_pointer, XRD_surface_signature))
			    {
			    xrd_dw_factor_surface_ratio = Get_FP_Number(char_pointer);
				xrd_dw_factor_surface_ratio -= 1.0;
				 //Set relative to change in the DW factor.
				xrd_surface_effects = TRUEV;
			    }
			 if (Command_Check(char_pointer, XRD_coorcutoff_signature))
			    {
				xrd_coor_cutoff = Get_FP_Number(char_pointer);
			    }

			 if (Command_Check(char_pointer, XRD_coornum_signature))
			    {
				xrd_coor_num = Get_Number(char_pointer);
			    }

			 if (Command_Check(char_pointer, XRD_DW_signature))
			    {
				xrd_dw_factor = Get_FP_Number(char_pointer);
			    }

			 Command_Check_Set_BooleanT(char_pointer, IONICSCAT_signature, xrd_ionic);

			 Command_Check_Set_BooleanT(char_pointer, IMAGEDWF_signature, image_dwf);

			 if (Command_Check(char_pointer, LOADMDPAM_signature) 
				 || Command_Check(char_pointer, RUNMD_signature))
			  //Potential loading...work in progress.
			    {
                const char POTENTIAL_FILE[50] = "CommandFile\\Potentials.txt";
			    }

			 if (Command_Check(char_pointer, MDPAM_signature))
			    {
				MPNL_CommandCheck_GetFPNum(char_pointer, TIMESTEP_signature, 
					                       mpam.time_step, 0.2);
				MPNL_CommandCheck_GetFPNum(char_pointer, INITIALTEMP_signature, 
					                       mpam.initial_temp, 298.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, FINALTEMP_signature,
										   mpam.final_temp, mpam.initial_temp);
				MPNL_CommandCheck_GetFPNum(char_pointer, INITIALTEMPTIME_signature,
					                       mpam.initial_temp_time, 0.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, FINALTEMPTIME_signature,
										   mpam.change_temp_time, 0.0);
				if (MPNL_CommandCheck(char_pointer, NVT_signature))
				   {
				   mpam.nvt_run = TRUEV;
				   }
				if (MPNL_CommandCheck(char_pointer, NPT_signature))
				   {
				   mpam.npt_run = TRUEV;
				   }
				MPNL_CommandCheck_GetFPNum(char_pointer, THERMOSTAT_signature, 
					                       mpam.thermo_coeff, 0.05);
				MPNL_CommandCheck_GetFPNum(char_pointer, BAROSTAT_signature, 
					                       mpam.baro_coeff, 0.05);
				MPNL_CommandCheck_GetFPNum(char_pointer, EQUI_signature, 
					                       mpam.equi_time, 100.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, PROD_signature, 
					                       mpam.prod_time, 100.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, SCREENUPDATE_signature, 
					                       mpam.supdate_time, 1.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, FILEUPDATE_signature, 
					                       mpam.fupdate_time, 1.0);
				if (MPNL_CommandCheck(char_pointer, MOLECULECOMMAND_signature))
				   {
				   mpam.molecule_commands = TRUEV;
				   }
			    }

			 if (Command_Check(char_pointer, SCOPEPAM_signature))
			    {
				if (MPNL_CommandCheck_SC(char_pointer, SCOPEMODE_signature)) 
				   {
				   if (char_pointer[0] == 'T')
				      {
					  spam.scope_mode = TEM_MODE;
				      }
				   if (char_pointer[0] == 'C')
				      {
					  spam.scope_mode = CBED_MODE;
				      }
				   }
				if (MPNL_CommandCheck(char_pointer, Nc_signature))
				   {
				   spam.cell_x = Get_Number(char_pointer);
				   spam.cell_y = spam.cell_z = spam.cell_x;
				   }
				if (MPNL_CommandCheck(char_pointer, Nx_signature))
				   {
				   spam.cell_x = Get_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, Ny_signature))
				   {
				   spam.cell_y = Get_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, Nz_signature))
				   {
				   spam.cell_z = Get_Number(char_pointer);
				   }				
			    if (MPNL_CommandCheck(char_pointer, Tc_signature))
				   {
				   spam.tilt_x = Get_FP_Number(char_pointer);
				   spam.tilt_y = spam.tilt_z = spam.tilt_x;
				   }
				if (MPNL_CommandCheck(char_pointer, Tx_signature))
				   {
				   spam.tilt_x = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, Ty_signature))
				   {
				   spam.tilt_y = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, Tz_signature))
				   {
				   spam.tilt_z = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, ScanMin_signature))
				   {
				   spam.scan_x_min = Get_FP_Number(char_pointer);
				   spam.scan_y_min = spam.scan_x_min;
				   }
				if (MPNL_CommandCheck(char_pointer, ScanMax_signature))
				   {
				   spam.scan_x_max = Get_FP_Number(char_pointer);
				   spam.scan_y_max = spam.scan_y_min;
				   }
				if (MPNL_CommandCheck(char_pointer, ScanXMin_signature))
				   {
				   spam.scan_x_min = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, ScanXMax_signature))
				   {
				   spam.scan_x_max = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, ScanYMin_signature))
				   {
				   spam.scan_y_min = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, ScanYMax_signature))
				   {
				   spam.scan_y_max = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, RasterX_signature))
				   {
				   res_temp = Get_FP_Number(char_pointer);
				   spam.raster_x = 
				   int( (spam.scan_x_max - spam.scan_x_min + FP_ERROR_FIX)/res_temp);
				    //Number of points for the beam position to raster across.
				   }
				if (MPNL_CommandCheck(char_pointer, RasterY_signature))
				   {
				   res_temp = Get_FP_Number(char_pointer);
				   spam.raster_y = 
				   int( (spam.scan_y_max - spam.scan_y_min + FP_ERROR_FIX)/res_temp);
				   }
				spam.scan_x_min *= 10.0;
				spam.scan_x_max *= 10.0;
				spam.scan_y_min *= 10.0;
				spam.scan_y_max *= 10.0;
				 //Convert to Angstroms.
				if (MPNL_CommandCheck(char_pointer, Raster_signature))
				   {
				   res_temp = Get_FP_Number(char_pointer);
				   spam.raster_x = 
				   int( (spam.scan_x_max - spam.scan_x_min + FP_ERROR_FIX)/res_temp);
				   spam.raster_y = 
				   int( (spam.scan_y_max - spam.scan_y_min + FP_ERROR_FIX)/res_temp);
				   }
				if (MPNL_CommandCheck(char_pointer, WindowSizeX_signature))
				   {
				   spam.window_x_size = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, WindowSizeY_signature))
				   {
				   spam.window_y_size = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, WindowSize_signature))
				   {
				   spam.window_x_size = Get_FP_Number(char_pointer);
				   spam.window_y_size = spam.window_x_size;
				   }
				if (MPNL_CommandCheck(char_pointer, XRes_signature))
				   {
				   spam.image_x_res = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, YRes_signature))
				   {
				   spam.image_y_res = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, Res_signature))
				   {
				   spam.image_x_res = Get_FP_Number(char_pointer);
				   spam.image_y_res = spam.image_x_res;
				   }
				spam.window_x_size *= 10.0;
				spam.window_y_size *= 10.0;
				spam.image_x_res *= 10.0;
				spam.image_y_res *= 10.0;
				 //Convert to Angstroms.
				MPNL_CommandCheck_GetFPNum(char_pointer, Slice_signature, 
					                       spam.slice_thickness, 0.15);
				spam.slice_thickness *= 10.0;
				MPNL_CommandCheck_GetFPNum(char_pointer, VOLTAGE_signature, 
					                       spam.voltage, 200.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, CURRENT_signature, 
					                       spam.current, 5.0);
			     //pA.
				MPNL_CommandCheck_GetFPNum(char_pointer, Brightness_signature,
										   spam.brightness, 150.0);
				 //Units = 10^8 A/cm^2sr.
				MPNL_CommandCheck_GetFPNum(char_pointer, DWELLTIME_signature,
										   spam.dwell_time, 2.0);
				 //Microseconds, sweet thang.
				MPNL_CommandCheck_GetFPNum(char_pointer, SOURCEANGLE_signature,
					                       spam.source_angle, 0.1);
				 //Source angle, in mrad.
				MPNL_CommandCheck_GetFPNum(char_pointer, SourceSize_signature,
										   spam.source_size, 0.0);
				 //Effective size of source, in nm.
				spam.source_size *= 10.0;
				 //Convert to Angstroms.
				if (MPNL_CommandCheck(char_pointer, BeamTiltX_signature))
				   {
				   spam.beam_tilt_x = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, BeamTiltY_signature))
				   {
				  spam.beam_tilt_y = Get_FP_Number(char_pointer);
				   }
				if (MPNL_CommandCheck(char_pointer, BeamTilt_signature))
				   {
				   spam.beam_tilt_x = Get_FP_Number(char_pointer);
				   spam.beam_tilt_y = spam.beam_tilt_x;
				   }
				MPNL_CommandCheck_GetFPNum(char_pointer, ConvAngle_signature,
					                       spam.convergence_angle, 15.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, DetectorIAngle_signature,
					                       spam.detector_inner_angle, 45.0);
				 //Default = 3.0 * Default convergence angle, for HAADF image.
				MPNL_CommandCheck_GetFPNum(char_pointer, DetectorOAngle_signature,
										   spam.detector_outer_angle, 200.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, DetectorIAngle_signature,
					                       spam.detector_inner_angle2, 0.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, DetectorOAngle_signature,
										   spam.detector_outer_angle2, 15.0);
				 //Default = Default convergence angle, for BF image.
				MPNL_CommandCheck_GetFPNum(char_pointer, CS3_signature, spam.Cs3, 0.005);
				MPNL_CommandCheck_GetFPNum(char_pointer, CS5_signature, spam.Cs5, 0.005);
				if (MPNL_CommandCheck_SC(char_pointer, Defocus_signature))
				   {
				   if (char_pointer[0] == 'S')
				      {
			          scope_sim param_logic(spam);
					  param_logic.Set_Scherzer_Defocus();
					  param_logic.Get_Parameters(spam);
				      }
				   else
				      {
					  spam.defocus = Get_FP_Number(char_pointer);
				      }
				   }
				MPNL_CommandCheck_GetFPNum(char_pointer, CC_signature, spam.Cc, 1.4);
				MPNL_CommandCheck_GetFPNum(char_pointer, DE_signature, spam.dE, 0.5);
				MPNL_CommandCheck_GetFPNum(char_pointer, ASTIG2_signature, 
					                       spam.two_astig, 0.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, ASTIGANGLE_signature,
					                       spam.two_astig_angle, 120.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, ASTIG3_signature,
					                       spam.three_astig, 0.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, ASTIGANGLE_signature,
					                       spam.three_astig_angle, 40.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, SCOPETEMP_signature,
					                       spam.scope_temp, 295.15);
				MPNL_CommandCheck_GetNum(char_pointer, TDScount_signature,
					                     spam.TDS_runs);
				if (MPNL_CommandCheck(char_pointer, SliceOutput_signature))
				   {
				   spam.slice_output_count = Get_Number(char_pointer);
				   }

				if (MPNL_CommandCheck(char_pointer, IDEALIMAGE_signature))
				   {
				   spam.Cs3 = 0.0;
				   spam.Cs5 = 0.0;
				   spam.defocus = 0.0;
				   spam.Cc = 0.0;
				   spam.dE = 0.0;
				   spam.two_astig = 0.0;
				   spam.three_astig = 0.0;
				   }
			    }

			 if (Command_Check(char_pointer, dw_coordination_parameter_signature))
			    {
				MPNL_CommandCheck_GetFPNum(char_pointer, "Min", dw_coordination_min, 0.0);
				MPNL_CommandCheck_GetFPNum(char_pointer, "Max", dw_coordination_max, 0.0);
			    }

			 if (Command_Check(char_pointer, CLEANANAL_signature) &&
				 (num_run_count == 0) && (text_run_count == 0) )
			  //At the start of a new series of tests, clean out the analysis file of 
			  //any previous results it contains.
			    {
				remove(ANALYSIS_FILE_NAME);
			    }

			 if (Command_Check(char_pointer, RIGIDMIN_signature_DEFAULT))
			  //Use default region (first group of atoms added to the composite).
			    {
			    Show_Statement("Performing rigid region energy minimization.");
				new_composite.Rigid_Minimization(first_region_added[0], 
					                             first_region_added[1], 
												 MD_hydrogen, 0.0);
			    }
				
			 if (Command_Check(char_pointer, RIGIDMIN_signature))
			  //User specifies the atomic indices used. 
			    {
				Get_Numbers(char_pointer, minimization_region, 2);
				 //Note: "first_region_added" variable is reused here, but its
				Show_Statement("Performing rigid region energy minimization.");
				--minimization_region[0];
				 //Put indices into computer logic.
				MPNL_CommandCheck_GetFPNum(char_pointer, "Precision", 
					                       rigid_precision, 0.0);
				new_composite.Rigid_Minimization(minimization_region[0], 
					                             minimization_region[1], MD_hydrogen,
												 rigid_precision);
			    }
				 
             }

	     //Clean-Up Routine and I/O Preparation---

         if (composite_name != NULL)
          //Get useful file-writing parameters for the composite system.
          //Note that most file-writing operations require specifically
          //a composite system to be passed.
             {    
			 if (clean_up)
			    {
				Show_Statement("Beginning Clean Up Routine!");
			    new_composite.Clean_Composite(cleanup_pam);
				 //Removes fragmented molecules.
			    }
             new_composite.Get_File_Runners(temp_collection);
			 new_composite.Get_Box_Size(temp_box_size);
             }

		 //I/O and analysis requests---

		 energy_set = FALSEV;
		 box_fit_set = FALSEV;
		 sphere_fit_set = FALSEV;
		 ele_comp_set = FALSEV;
		 unlike_distavg_set = FALSEV;
		 like_distavg_set = FALSEV;
		 gen_distavg_set = FALSEV;
		 strangeness_set = FALSEV;
		 electro_potential_set = FALSEV;
		 energy_potential_set = FALSEV;
		 analysis_request = FALSEV;
         for (int index = 0; index < file_size; ++index)
          //Look for low-priority general commands.
             { 
             char_pointer = &(file_text[index]);

             if (Command_Check(char_pointer, GULP_signature))
              //Write a GULP optimization input file for the defined composite.
                { 
				Make_File_Name(temp_string, composite_name, GULP_ext);
				Make_File_Name(temp_dump_string, composite_name, GULPDUMP_ext);
				Show_Statement("Writing GULP Optimization File!");
                Write_GULP_File(temp_collection, temp_string, temp_dump_string, 
					            MD_hydrogen, MD_mass, cat_coreshell, ani_coreshell);                             
                }
             
             if (Command_Check(char_pointer, GULPMELT_signature))
              //Write a GULP MD input file for the defined composite.
                { 
				Make_File_Name(temp_string, composite_name, GULP_ext);
				Make_File_Name(temp_dump_string, composite_name, RES_ext);
				Show_Statement("Writing GULP Molecular Dyonamics File!");
                Write_GULP_Melt_File(temp_collection, temp_string, temp_dump_string, 
					                 MD_hydrogen, MD_mass, cat_coreshell, ani_coreshell,
								     mpam);                        
                }   

			 if (Command_Check(char_pointer, PDB_signature))
			    {
				Make_File_Name(temp_string, composite_name, PDB_ext);
				Show_Statement("Writing PDB File!");
				Write_PDB_File(temp_collection, temp_string, MD_hydrogen);
			    }

			 if (Command_Check(char_pointer, LAMMPS_signature))
			    {
				Make_File_Name(temp_string, composite_name, LAMMPS_ext);
				Show_Statement("Writing LAMMPS File!");
				Write_LAMMPS_File(temp_collection, temp_string, MD_hydrogen);
			    }
                       
             if (Command_Check(char_pointer, QSTEM_signature))
                {   
				Make_File_Name(temp_string, composite_name, QSTEM_ext);
				Show_Statement("Writing QSTEM File!");
                Write_QSTEM_File(temp_collection, temp_string, 
					             image_hydrogen, image_dwf);                      
                }

			 if (Command_Check(char_pointer, QCONF_signature))
			    {
				Make_File_Name(temp_string, composite_name, QSTEM_ext);
				Make_File_Name(temp_string2, composite_name, QCONF_ext);
				double cell_vec[3];
				temp_collection.Get_Box_Vectors(cell_vec);
				Show_Statement("Writing QSC File!");
                Write_QSC_File(temp_string2, temp_string, spam, cell_vec[2]); 	
			    }
                
			 if (Command_Check(char_pointer, QPBS_signature))
				{
				Make_File_Name(temp_string, composite_name, PBS_ext);
				Show_Statement("Writing QSTEM PBS File!");
                Write_QPBS_File(temp_string, QCONF_ext); 
				}

             if (Command_Check(char_pointer, JEMS_signature))
                { 
				Make_File_Name(temp_string, composite_name, JEMS_ext);
				Show_Statement("Writing JEMS File!");
                Write_JEMS_File(temp_collection, temp_string, image_hydrogen, image_dwf);               
                }
                
			 if (Command_Check(char_pointer, XYZ_signature))
                {  
				Make_File_Name(temp_string, composite_name, XYZ_ext);
				Show_Statement("Writing XYZ File!");
                Write_XYZ_File(temp_collection, temp_string, image_hydrogen);
                }

			 if (Command_Check(char_pointer, GSPE_signature))
				{
				Make_File_Name(temp_string, composite_name, GSPE_ext);
				Show_Statement("Writing GSPE File!");
                Write_GSPE_File(temp_collection, temp_string, MD_hydrogen);
				}

             if (Command_Check(char_pointer, RADIAL_signature) && !no_rad_calc)
                {
				Make_File_Name(temp_string, composite_name, radial_ext);
				Show_Statement("Calculating RDF!");
                Write_RAD_File(temp_collection, temp_string, radial_step);                        
                }

			 if (Command_Check(char_pointer, XRD_signature))
			    {
				Make_File_Name(temp_string, composite_name, xrd_ext);
				Show_Statement("Calculting XRD Pattern!");
				Write_XRD_File(temp_collection, temp_string, xrd_wave, xrd_ionic, xrd_precision,
							   xrd_min_angle, xrd_max_angle, xrd_dw_factor,
							   xrd_surface_effects, xrd_dw_factor_surface_ratio, 
							   xrd_coor_cutoff, xrd_coor_num);
			    }

			 if (Command_Check(char_pointer, TRAJDWF_signature))
			    {
				Make_File_Name(temp_string, composite_name, HIST_ext);
				Make_File_Name(temp_string2, composite_name, dwf_ext);
				Remove_Root(temp_string);
				Show_Statement("Calculating Debye-Waller Factors for the Atoms!");
				Write_DWF_File(temp_collection, temp_string2, temp_string, anal_hydrogen,
							   dw_coordination_min, dw_coordination_max);
			    }

			              
			 if (Command_Check(char_pointer, IMAGEREAD_signature))
			 //.img file analysis
				{
				MPSNL_CommandCheck_SC(char_pointer, temp_image_file_pointer, "File");
				MPNL_CommandCheck_GetFPNum(char_pointer, "Source", gun_size);
				MPNL_CommandCheck_GetFPNum(char_pointer, "Stats", image_stat_pam, -1);
				line_scan_params[2] = 0;
				restriction_pam = 0;
				MPNL_CommandCheck_GetNum(char_pointer, "Restrict", restriction_pam);
				MPNL_CommandCheck_GetNums(char_pointer, "Line", line_scan_params, 3);
				Move_Pointer_NextLine(char_pointer);
				show_num = FALSEV;
			    Command_Check_Set_BooleanT(char_pointer, "ShowNum", 
				                           show_num); 
				Add_File_Extension(temp_image_file_pointer, "img");
				Show_Statement("Analyzing the .img file");
				Interpret_Image(temp_image_file_pointer, gun_size, image_stat_pam, 
					             restriction_pam, line_scan_params, show_num);
				}

			 if (Command_Check(char_pointer, DISTCOUNT_signature))
			    {
				Make_File_Name(temp_string, composite_name, distcount_ext);
				Show_Statement("Writing Distance Histogram File!");
                Write_DISTCOUNT_File(temp_collection, temp_string, histo_pam,
									 anal_hydrogen);
			    }

			 if (Command_Check(char_pointer, COORCOUNT_signature))
			    {
				Make_File_Name(temp_string, composite_name, coorcount_ext);
				Show_Statement("Writing Coordination Histogram File!");
				Write_COORCOUNT_File(temp_collection, temp_string,
									 anal_hydrogen);
			    }

			 if (Command_Check(char_pointer, distavg_signature))
			    {
				MPNL_CommandCheck_GetFPNum(char_pointer, "Min", min);
				MPNL_CommandCheck_GetFPNum(char_pointer, "Max", max);
				unlike_only = FALSEV;
				if (MPNL_CommandCheck_SC(char_pointer, "AB"))
				   {
				   if (char_pointer[0] == 'Y')
				      {
					  unlike_only = TRUEV;
				      }
				   }
				like_only = FALSEV;
				if (MPNL_CommandCheck_SC(char_pointer, "AA"))
				   {
				   if (char_pointer[0] == 'Y')
					  {
					  like_only = TRUEV;
				      }
				   } 
				distavg = temp_collection.Distance_Average(min, max, unlike_only,
							   like_only, anal_hydrogen);
				cooravg = temp_collection.Coordination_Average(min, max, unlike_only,
							   like_only, anal_hydrogen);
				cout << endl << "------DISTANCE AVERAGE YOU REQUESTED IS----"
					 << endl << " From " << min << " to " << max << ": " 
					 << distavg << " nm. Corresponding coordiation number average is: " 
					 << cooravg;
				if (unlike_only)
				   {
				   cout << endl << "Note: This is for unlike (AB) distances only.";
				   unlike_distavg = distavg;
				   unlike_cooravg = cooravg;
				   unlike_distavg_set = TRUEV;
				   }  
				else if (like_only)
				   {
				   cout << endl << "Note: This is for like (AA) distances only.";
				   like_distavg = distavg;
				   like_cooravg = cooravg;
				   like_distavg_set = TRUEV;
				   }
				else if (!like_only && !unlike_only)
				   {
				   gen_distavg = distavg;
				   gen_cooravg = cooravg;
				   gen_distavg_set = TRUEV;
				   }
				}

			 if (Command_Check(char_pointer, shapefit_signature))
			    {
				shape_identity = SPHERE_SHAPE;
				if (MPNL_CommandCheck_SC(char_pointer, "Shape"))
				   {
				   if (char_pointer[0] == 'B')
				      {
					  shape_identity = BOX_SHAPE;
				      }
				   }
				MPNL_CommandCheck_GetFPNum(char_pointer, "NNDist", max);
				unlike_only = FALSEV;
				if (MPNL_CommandCheck_SC(char_pointer, "AB"))
				   {
				   if (char_pointer[0] == 'Y')
				      {
					  unlike_only = TRUEV;
				      }
				   }
				MPNL_CommandCheck_GetNum(char_pointer, "BulkCN", coor_num);
				if (!MPNL_CommandCheck_GetNum(char_pointer, "AtomNum", 
					shapefit_atom_number_specification))
				   {
				   shapefit_atom_number_specification = -1;
				    //Reset to 0 if command is not found. Must be -1.
				   }

				if (shape_identity == SPHERE_SHAPE)
				   {
				   temp_collection.Fit_Sphere(max, unlike_only, coor_num, sphere_fit, 
											  shapefit_atom_number_specification,
											  anal_hydrogen);
				   sphere_fit_set = TRUEV;
				   }
				if (shape_identity == BOX_SHAPE)
				   {
				   temp_collection.Fit_Box(max, unlike_only, coor_num, box_fit, 
										   shapefit_atom_number_specification,
										   anal_hydrogen);
				   box_fit_set = TRUEV;
				   }
				cout << endl << "------"; 
				if (shape_identity == SPHERE_SHAPE) cout << "SPHERE";
				if (shape_identity == BOX_SHAPE) cout << "BOX";
				cout << " FIT YOU REQUESTED IS----";
				if (shape_identity == SPHERE_SHAPE) 
					{
					cout << endl << "Radius is : " << sphere_fit << " nm";
					}
				if (shape_identity == BOX_SHAPE)
				    {
					cout << endl << "Box lengths are: " << box_fit[0] << " " 
						 << box_fit[1] << " " << box_fit[2] << " nm";
				    }
			    }
			 
			 if (Command_Check(char_pointer, ENERGY_signature))
			  //Uses .frc file information.
			    {
                strcpy(temp_string, composite_name);
                Add_File_Extension(temp_string, FORCE_ext);
				char dump_string[30];
				ifstream energy_file(temp_string, ios::in);
				energy_value = 0.0;
				dump_string[0] = '\0';
				if (energy_file.is_open())
				   {
				   energy_file >> dump_string >> energy_value;
				   if (Command_Check(dump_string, "energy"))
					  {
					  energy_set = TRUEV;
					  }
				   cout << endl << "System Energy is: " << energy_value << " eV";
				   energy_file.close();
				   }
				else
				   {
				   Show_Warning("COULD NOT OPEN ENERGY FILE!");
				   }
			    }
			 
			 if (Command_Check(char_pointer, ELECOMP_signature))
				   {
				   temp_collection.Get_Composition(ele_fit_string, ELE_STRING_LENGTH);
				   cout << endl << "Elemental Composition is: " << ele_fit_string;
				   ele_comp_set = TRUEV;
				   }

			 if (Command_Check(char_pointer, STRANGE_signature))
			       {
				   cout << endl << "Determining if there is something strange.";
				   is_strange = temp_collection.Check_Strange_Location(strangeness_dist_min, 
					                            strangeness_dist_max, anal_hydrogen);
				   if (is_strange)
					  {
					  Show_Warning("SOMETHING IS STRANGE ABOUT THIS SYSTEM!");
				      }
				   strangeness_set = TRUEV;
			       }

			 if (Command_Check(char_pointer, ELECTROPOT_signature))
			       {
				   Show_Statement("Calculating the electrostatic potential.");
				   electro_potential = temp_collection.Calculate_ES_Potential(MD_hydrogen);
				   cout.precision(11);
				   cout << endl << "Electrostatic potential is: " 
					    << electro_potential << " eV.";
				   cout.precision(6);
				   electro_potential_set = TRUEV;
			       }

			 if (Command_Check(char_pointer, ENERGYPOT_signature))
			       {
				   /*cout << endl << "Calculating total potential energy.";
				   //energy_potential = temp_collection.Calculate_Fast_Total_Potential(MD_hydrogen);
				   cout.precision(11);
				   //cout << endl << "Total energy potential is: " 
				 	      << energy_potential << " eV.";
				   cout.precision(6);
				   energy_potential_set = TRUEV;*/
			       }

			 if (Command_Check(char_pointer, ANALYSIS_signature))
				   {
				   analysis_request = TRUEV;
				   }
                
             if (Command_Check(char_pointer, FILEDUMP_signature))
              //Make a save file for EVERYTHING defined in the script file.
                { 
			    Show_Statement("Making records of all components used!");
                if (composite_name != NULL)
                   {
                   strcpy(temp_string, composite_name);
                   Make_File_Name(temp_string, composite_ext);
                   new_composite.Save_Composite(temp_string);                  
                   }
                for (int a = 0; a < struct_index; ++a)
                    {
                    strcpy(temp_string, struct_names[a]);
                    Make_File_Name(temp_string, cryst_ext);
                    new_crystal_systems[a].Save_Crystal_System(temp_string);     
                    }
                for (int a = 0; a < molecule_index; ++a)
                    {
                    strcpy(temp_string, molecule_names[a]);
                    Make_File_Name(temp_string, molecule_ext);
                    new_molecules[a].Save_Molecule(temp_string);      
                    }
                for (int a = 0; a < particle_index; ++a)
                    {
                    strcpy(temp_string, nanoparticle_names[a]);
                    Make_File_Name(temp_string, nanoparticle_ext);
                    new_nanoparticles[a].Save_NanoParticle(temp_string); 
                    }        
                for (int a = 0; a < surface_index; ++a)
                    {
                    strcpy(temp_string, nanosurface_names[a]);
                    Make_File_Name(temp_string, nanosurface_ext);
                    new_nanosurfaces[a].Save_NanoSurface(temp_string);     
                    }                
                }
            }
	 
	  if (analysis_request)
		    {
			Show_Statement("Adding information to the analysis file!");
			Analysis_Out(ANALYSIS_FILE_NAME, composite_name, 
				         temp_collection.Real_Size(MD_hydrogen), 
				         ele_comp_set, ele_fit_string,
						 energy_set, energy_value,
						 unlike_distavg_set, unlike_distavg, unlike_cooravg,
						 like_distavg_set, like_distavg, like_cooravg,
						 gen_distavg_set, gen_distavg, gen_cooravg,
						 sphere_fit_set, sphere_fit,
						 box_fit_set, box_fit,
						 strangeness_set, is_strange,
						 electro_potential_set, electro_potential,
						 energy_potential_set, energy_potential);			
	        }
	
	  delete[] file_text;
	  delete[] nauto_text;
      }
      
///Serial reading functions---      
      
void MoleculeAdd_Reader(const char* script_reader, double* previous, double* start, 
                        char** atom_names)
      {
	  MPNL_CommandCheck_GetFPNums(script_reader, "Prev", previous, 3);
	  MPNL_CommandCheck_GetFPNums(script_reader, "Next", start, 3);
	  MPNL_CommandCheck_GetNames(script_reader, "Atoms", 5, atom_names, 10);                       
      }
  
void Convert_ShapeFactor(char second_letter, int& shape_factor)
     {
     if (second_letter == 'E')
        {
        shape_factor = TETRA;                 
        }
     else if (second_letter == 'R')
        {
        shape_factor = TRIG;     
        }
     else if (second_letter == 'I')
        {
        shape_factor = LINEAR;     
        }                         
     }  

void Set_Miller_X_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(1, 0, 0, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, 0, 0, scaling_pam);
	 }

void Set_Miller_Y_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(0, 1, 0, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, -1, 0, scaling_pam);
	 }	

void Set_Miller_Z_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(0, 0, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, 0, -1, scaling_pam);
	 }

void Set_Miller_X_Y_Z_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
	 Set_Miller_X_Facets(ze_cryst, scaling_pam);
	 Set_Miller_Y_Facets(ze_cryst, scaling_pam);
	 Set_Miller_Z_Facets(ze_cryst, scaling_pam);
	 }

void Set_Miller_XY_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(1, 1, 0, scaling_pam);
     ze_cryst.Set_Faceting_Plane(1, -1, 0, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, 1, 0, scaling_pam);
     ze_cryst.Set_Faceting_Plane(-1, -1, 0, scaling_pam);
	 }

void Set_Miller_YZ_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
     ze_cryst.Set_Faceting_Plane(0, 1, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, -1, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, 1, -1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, -1, -1, scaling_pam);	 
	 }

void Set_Miller_XZ_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(1, 0, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(1, 0, -1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, 0, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, 0, -1, scaling_pam);
	 }

void Set_Miller_XY_YZ_XZ_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
	 Set_Miller_XY_Facets(ze_cryst, scaling_pam);
	 Set_Miller_YZ_Facets(ze_cryst, scaling_pam);
	 Set_Miller_XZ_Facets(ze_cryst, scaling_pam);
	 }

void Set_Miller_XYZ_Facets_Set_Zero_Two(crystal_system& ze_cryst, double scaling_pam)
	 {
     ze_cryst.Set_Faceting_Plane(1, 1, 1, scaling_pam);
     ze_cryst.Set_Faceting_Plane(-1, -1, 1, scaling_pam);
     ze_cryst.Set_Faceting_Plane(-1, 1, -1, scaling_pam);
     ze_cryst.Set_Faceting_Plane(1, -1, -1, scaling_pam);  
	 }

void Set_Miller_XYZ_Facets_Set_One_Three(crystal_system& ze_cryst, double scaling_pam)
	 {
     ze_cryst.Set_Faceting_Plane(-1, 1, 1, scaling_pam);
     ze_cryst.Set_Faceting_Plane(1, -1, 1, scaling_pam);
     ze_cryst.Set_Faceting_Plane(1, 1, -1, scaling_pam);
     ze_cryst.Set_Faceting_Plane(-1, -1, -1, scaling_pam);
	 }

void Set_Miller_XYZ_Facets_PosZ(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(1, 1, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, 1, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(1, -1, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, -1, 1, scaling_pam);
	 }

void Set_Miller_XYZ_Facets_NegZ(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(1, 1, -1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, 1, -1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(1, -1, -1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, -1, -1, scaling_pam);
	 }

void Set_Miller_XYZ_Facets(crystal_system& ze_cryst, double scaling_pam)
	 {
	 Set_Miller_XYZ_Facets_Set_Zero_Two(ze_cryst, scaling_pam);
	 Set_Miller_XYZ_Facets_Set_One_Three(ze_cryst, scaling_pam);
	 }

void Set_Miller_XZ_YZ_Facets_PosZ(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(1, 0, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, 0, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, 1, 1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, -1, 1, scaling_pam);
	 }

void Set_Miller_XZ_YZ_Facets_NegZ(crystal_system& ze_cryst, double scaling_pam)
	 {
	 ze_cryst.Set_Faceting_Plane(1, 0, -1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(-1, 0, -1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, 1, -1, scaling_pam);
	 ze_cryst.Set_Faceting_Plane(0, -1, -1, scaling_pam);
	 }

void Set_Miller_Round_Edge_Facets_NegY(crystal_system& ze_cryst, double scaling_pam)
	 {
	 for (int a = -2; a < 3; ++a)
	  //h (index) variation.
		 {
		 for (int b = -2; b < 0; ++b)
		  //k variation.
			 {
			 for (int c = -2; c < 3; ++c)
			  //l variation.
				 {
				 ze_cryst.Set_Faceting_Plane(a, b, c, scaling_pam);
				 }
			 }
		 }
	 }

void Set_Miller_Round_Edge_Facets_NegZ(crystal_system& ze_cryst, double scaling_pam)
	 {
	 for (int a = -2; a < 3; ++a)
	  //h (index) variation.
		 {
		 for (int b = -2; b < 3; ++b)
		  //k variation.
			 {
			 for (int c = -2; c < 0; ++c)
			  //l variation.
				 {
				 ze_cryst.Set_Faceting_Plane(a, b, c, scaling_pam);
				 }
			 }
		 }
	 }

void Set_Faceting_Planes(crystal_system& ze_cryst, const char* facet_set_name, double scaling_parameter)
	 { 
	 const double f1 = 1.0;
	  //Basic scaling factor for faceting = unity.
     if (facet_set_name[0] == 'T')
       //Tetrahedral.
        {
	    Set_Miller_XYZ_Facets_Set_Zero_Two(ze_cryst, f1);
		}
	 else if (facet_set_name[0] == 'O')
	  //Octahedral.
		{
		Set_Miller_XYZ_Facets(ze_cryst, f1);
		}
	 else if ( (facet_set_name[0] == 'C') && (facet_set_name[3] != 'O') )
	   //Cubic.
		{
		Set_Miller_X_Y_Z_Facets(ze_cryst, f1);
		 //{100} family.
		if (facet_set_name[6] == 'R')
		 //CUBIC_ROUNDED system. Not just a basic cube...
		   {
		   const double f2 = 1.13 * scaling_parameter;
		     //Ratio of cutting planes distances determined by
	     	 //Liu's article and the dimensions of the particle seen
			 //experimentally (for scaling parameter = 1.0).
		   Set_Miller_XY_YZ_XZ_Facets(ze_cryst, f2);
			//{110} family.
		   }
	    }
	  else if (facet_set_name[0] == 'C')
	  //Cubo-octohedral.
        {
        const double f2 = 1.2 * scaling_parameter;
		Set_Miller_XYZ_Facets(ze_cryst, f1);
		Set_Miller_X_Y_Z_Facets(ze_cryst, f2);   
        }
	 else if (facet_set_name[0] == 'B')
	  //Block-end cylinder (analogous to a nanosurface with a = b != c).
		{
		const double f2 = 1.0*scaling_parameter;
		 //Short direction (scaling parameter < 1.0), which is
		 //the x and y directions.
		Set_Miller_X_Facets(ze_cryst, f2);
		Set_Miller_Y_Facets(ze_cryst, f2);
		Set_Miller_Z_Facets(ze_cryst, f1);
		 //{100} family.
		if (facet_set_name[6] == 'C')
		 //BLOCK_CHIPPED. Chipped surface at corners of block-end.
		   {
		   const double f3 = 0.5 + (0.63 * scaling_parameter);
		   Set_Miller_XZ_Facets(ze_cryst, f3);
		   Set_Miller_YZ_Facets(ze_cryst, f3);
		   }
		}
	 else if (Command_Check(facet_set_name, "Diamond") )
	  //Diamond.
		{
		const double f2 = f1;
		const double f3 = 2.2;
		Set_Miller_XYZ_Facets_NegZ(ze_cryst, f1);
		Set_Miller_XZ_YZ_Facets_NegZ(ze_cryst, f2);
		 //Base of diamond.
		Set_Miller_X_Y_Z_Facets(ze_cryst, f3);
		Set_Miller_XZ_YZ_Facets_PosZ(ze_cryst, f3);
		Set_Miller_XYZ_Facets_PosZ(ze_cryst, f3);
		 //Top of diamond.
		}
	 else if (Command_Check(facet_set_name, "Disc"))
		{
		const double f2 = 0.2 * scaling_parameter;
		ze_cryst.Set_Faceting_Plane(0, 0, 1, f2);
		ze_cryst.Set_Faceting_Plane(0, 0, -1, f2);
		}
	 else if (facet_set_name[0] == 'D')
	  //Dodecahedral.
		{
		Set_Miller_XY_YZ_XZ_Facets(ze_cryst, f1);
		}
	 else if (facet_set_name[0] == 'P')
	 //(Floating) pod.
		{
		const double f2 = f1 * scaling_parameter;
		ze_cryst.Set_Faceting_Plane(0, 0, 1, f2);
		}
	 else if (Command_Check(facet_set_name, "Arc"))
		{
		const double f2 = 0.3 * scaling_parameter;
		const double f3 = 0.01;
		Set_Miller_Z_Facets(ze_cryst, f2);
		Set_Miller_Round_Edge_Facets_NegY(ze_cryst, f1); 
		Set_Miller_X_Facets(ze_cryst, f1);
		ze_cryst.Set_Faceting_Plane(0, 1, 0, f3);
		}
	 else if (Command_Check(facet_set_name, "Finger"))
		{
		const double f2 = 0.1 * scaling_parameter;
		const double f3 = 0.2;
		Set_Miller_Z_Facets(ze_cryst, f2);
		Set_Miller_Round_Edge_Facets_NegY(ze_cryst, f3); 
		Set_Miller_X_Facets(ze_cryst, f3);
		ze_cryst.Set_Faceting_Plane(0, 1, 0, f1);
		}
	 }

///Automation functioning---

void Process_SAutomation_Commands(char* script_text, int script_size,
	                              char* nauto_text, int nauto_size,
								  const char* sauto_text, int sauto_size,
								  int run_count)
	{
	if (run_count == 0)
	   {
	   return;
	   }

	//Determine the string-replacement automation requests
	//by reading the sauto file---

	const int MAX_REPLACEMENT_COUNT = 1000;
	 //Maximum number of replacement commands.
	int num_replacements = 0;
	char* hunting_pointers[MAX_REPLACEMENT_COUNT];
	char* replacement_pointers[MAX_REPLACEMENT_COUNT];
	const char* sauto_temp;
	int run_temp;
	int rep_count = 0;
	for (int a = 0; a < sauto_size; ++a)
	    {
        sauto_temp = &(sauto_text[a]);
          //Investigate a particular character in the automation text.
        if (Command_Check(sauto_temp, "GEN_REPLACE"))
		   {
		   Move_Pointer_NextLine(sauto_temp);
		   do
		    //Get replacement commands.
		         {
				 Move_Pointers_NextLine(sauto_temp, hunting_pointers[num_replacements]);
				 ++num_replacements;

				 if (num_replacements == MAX_REPLACEMENT_COUNT)
				    {
					Show_Warning("TOO MANY REPLACEMENT COMMANDS IN STRING AUTOMATION!");
					break;
					}
		         }
		   while (!Command_Check(sauto_temp, "RUN_@X@") && 
			      !Command_Check(sauto_temp, "END_@X@"));
		   num_replacements -= 2;
		    //Get rid of strings for last two lines after hunting strings.
		   sauto_temp = &(sauto_text[a]);
		    //Back to start.
		   }
		if (Command_Check(sauto_temp, "RUN_@X@"))
		   {
		   run_temp = Get_Number(sauto_temp);
		   if (run_temp == (run_count + 1))
		    //Found the needed replacement information.
		           {
				   Move_Pointer_NextLine(sauto_temp);
				   do
		             {
				     Move_Pointers_NextLine(sauto_temp, replacement_pointers[rep_count]);
				     ++rep_count;
		             }
		           while (!Command_Check(sauto_temp, "RUN_@X@") && 
					      !Command_Check(sauto_temp, "END_@X@") &&
					       (rep_count != num_replacements) );
				    a = sauto_size;
		            }
		   }
	    }

	if (rep_count != num_replacements)
	   {
	   Show_Warning("NEED TO SPECIFY ALL REPLACEMENTS IN STRING AUTOMATION!");
	   return;
	   }

	//Apply the string automation commands to both number automation and
	//main scripting file---

	char* nauto_temp;
	for (int a = 0; a < nauto_size; ++a)
	    {
		nauto_temp = &(nauto_text[a]);
		for (int b = 0; b < num_replacements; ++b)
		    {
			if (Command_Check(nauto_temp, hunting_pointers[b]))
			 //Found what we are hunting to change, so change it.
			   {
			   String_Replace(nauto_temp, strlen(hunting_pointers[b]), nauto_size - a,
				              replacement_pointers[b], strlen(replacement_pointers[b]));
			   }
		    }
	    }

    char* script_temp;
	for (int a = 0; a < script_size; ++a)
	    {
		script_temp = &(script_text[a]);
		for (int b = 0; b < num_replacements; ++b)
		    {
			if (Command_Check(script_temp, hunting_pointers[b]))
			 //Found what we are hunting to change, so change it.
			   {
			   String_Replace(script_temp, strlen(hunting_pointers[b]), script_size - a,
				              replacement_pointers[b], strlen(replacement_pointers[b]));
			   }
		    }
	    }
	}	

void Process_NAutomation_Commands(char* script_text, int script_size, 
	                              const char* auto_text, int auto_size, 
				                  int run_count)
     {
	 if (run_count == 0)
	    {
		return;
	    }

	 //Determine the number automation commands by reading
	 //the number automation file---

	 char* script_temp; 
	 const char* auto_temp;
	 int num_auto_commands = 0;
     const int MAX_AUTO_COMMANDS = 1000;
	 char* hunting_pointers[MAX_AUTO_COMMANDS];
	 bool is_FP[MAX_AUTO_COMMANDS];
	 bool is_1D[MAX_AUTO_COMMANDS];
	 int integer_changes[MAX_AUTO_COMMANDS][3];
	 double double_changes[MAX_AUTO_COMMANDS][3];

     for (int a = 0; a < auto_size; ++a)
	    {
        auto_temp = &(auto_text[a]);
          //Investigate a particular character in the automation text.
        if (Command_Check(auto_temp, "TEXT_REPLACE_SEQ"))
		   {
		   Move_Pointers_NextLine(auto_temp, hunting_pointers[num_auto_commands]);

		   if (MPNL_CommandCheck_GetFPNums(auto_temp, "FPValsChange", 
			        &(double_changes[num_auto_commands][0]), 3) )
		      {
			  is_FP[num_auto_commands] = TRUEV;
			  is_1D[num_auto_commands] = FALSEV;
		      }
		   else if (MPNL_CommandCheck_GetFPNum(auto_temp, "FPValChange", 
			         double_changes[num_auto_commands][0]) )
		      {
			  is_FP[num_auto_commands] = TRUEV;
			  is_1D[num_auto_commands] = TRUEV;
		      }
		   else if (MPNL_CommandCheck_GetNums(auto_temp, "IValsChange", 
			        &(integer_changes[num_auto_commands][0]), 3) )
		      {
			  is_FP[num_auto_commands] = FALSEV;
			  is_1D[num_auto_commands] = FALSEV;
		      }
		   else if (MPNL_CommandCheck_GetNum(auto_temp, "IValChange", 
			         integer_changes[num_auto_commands][0]) )
		      {
			  is_FP[num_auto_commands] = FALSEV;
			  is_1D[num_auto_commands] = TRUEV;
		      }
		   else
		      {
			  Show_Warning("ONE OF THE SUPER-AUTOMATION COMMANDS IS NOT RECOGNIZED!");
		      }
		   ++num_auto_commands;
		   }
	    }

	 //Apply number automation commands---

	 int temp_integer_changes[3];
	 double temp_double_changes[3];
	 for (int a = 0; a < script_size; ++a)
	    {
	    script_temp = &(script_text[a]);
		for (int b = 0; b < num_auto_commands; ++b)
		  {
		  if (Command_Check(script_temp, hunting_pointers[b]))
		     {
			 if (is_1D[b] && !is_FP[b])
			    {
				temp_integer_changes[0] = run_count * integer_changes[b][0];
				Change_Integer(script_temp, temp_integer_changes[0]);
			    }
			 else if (!is_1D[b] && !is_FP[b])
			    {
				for (int c = 0; c < 3; ++c)
				    {
					temp_integer_changes[c] = run_count * integer_changes[b][c];
				    }
				Change_Integers(script_temp, temp_integer_changes, 3);
			    }
			 else if (is_1D[b] && is_FP[b])
			    {
				temp_double_changes[0] = run_count * double_changes[b][0];
				Change_FP(script_temp, temp_double_changes[0]);
			    }
			 else if (!is_1D[b] && is_FP[b])
			    {
				for (int c = 0; c < 3; ++c)
				    {
					temp_double_changes[c] = run_count * double_changes[b][c];
				    }
				Change_FPS(script_temp, temp_double_changes, 3);
			    }
		     }
		  }
	    }
	 }

///Analysis functions---

void Analysis_Out(const char* analysis_file_name, const char* composite_name, int real_size,
				  bool ele_comp_set, char* ele_fit_string,
				  bool energy_set, double energy_value,
				  bool unlike_distavg_set, double unlike_distavg, double unlike_cooravg,
				  bool like_distavg_set, double like_distavg, double like_cooravg,
				  bool gen_distavg_set, double gen_distavg, double gen_cooravg,
				  bool sphere_fit_set, double sphere_fit,
				  bool box_fit_set, double* box_fit,
				  bool strangeness_set, bool is_strange,
				  bool electro_potential_set, double electro_potential,
				  bool energy_potential_set, double energy_potential)
	 {	 
	 ifstream anal_Ifile(analysis_file_name, ifstream::in);
	 if (!anal_Ifile.is_open())
		{
		Show_Warning("COULD NOT OPEN THE ANALYSIS FILE!");
		return;
		}

	 const int MAX_ELEMENT_LENGTH = 100;
	 const int MAX_ANALYSIS_ELEMENTS = 30;
	 char temp_string[MAX_ELEMENT_LENGTH];
	 char command_list[MAX_ELEMENT_LENGTH*MAX_ANALYSIS_ELEMENTS];
	 bool must_add_labels = TRUEV;
	 bool found_labels = FALSEV;

	 //Check to see if the analysis file already has the list of analysis
	 //elements that is now needed. If so, appending to that list 
	 //will be performed (otherwise, appending to the end of the file)---

	 char space_string[3] = " ";
	 char end_char_string[3] = "|";
	 while  (!anal_Ifile.eof())
		{ 
		anal_Ifile >> temp_string;
		if (Command_Check(temp_string, "Name"))
		 //Get last instance of label listing in the file and check if it
		 //is the same as the labels needed.
		   {
		   command_list[0] = '\0';
		   strcat(command_list, temp_string);
		   found_labels = TRUEV;
		   for (int a = 0; a < MAX_ANALYSIS_ELEMENTS; ++a)
			   {
			   strcat(command_list, space_string);
			   anal_Ifile >> temp_string;
			   strcat(command_list, temp_string);
			   if (Command_Check(temp_string, end_char_string))
				  {
				  a = MAX_ANALYSIS_ELEMENTS;
			      }
		       }
		   }
		}

	 //Create the label string that matches the user request---

	 char label_string[MAX_ELEMENT_LENGTH*MAX_ANALYSIS_ELEMENTS];
	 label_string[0] = '\0';
	 strcat(label_string, "Name Atom_Count ");
     if (ele_comp_set) strcat(label_string, "Elem_Comp ");
	 if (energy_set) strcat(label_string, "Energy ");
	 if (unlike_distavg_set) strcat(label_string, "AB_Dist AB_Coord ");
	 if (like_distavg_set) strcat(label_string, "AA_Dist AA_Coord ");
	 if (gen_distavg_set) strcat(label_string, "Gen_Dist Gen_Coord ");
	 if (sphere_fit_set) strcat(label_string, "Sphere_Fit_Radius ");
	 if (box_fit_set) 
		{
		strcat(label_string, "Box_Fit_X Box_Fit_Y Box_Fit_Z ");
		}
	 if (strangeness_set) strcat(label_string, "Bad_Structure ");
	 if (electro_potential_set) strcat(label_string, "ES_Potential ");
	 if (energy_potential_set) strcat(label_string, "TOT_Potential ");
	 strcat(label_string, end_char_string);

	 //Check if it matches the string list currently in the analysis file---

	 if (found_labels && Command_Check(label_string, command_list))
		{
		must_add_labels = FALSEV;
		}

	 //Writing preparation----

	 ofstream anal_Ofile(analysis_file_name, ofstream::out | ofstream::app);
	  //Open file in appending mode.

	 //Label header addition, if not in the file already---

	 if (must_add_labels)
		{
		anal_Ofile << endl << endl << label_string;
		}

	 //Data addition---

	 anal_Ofile << endl << composite_name << " " << real_size << " ";
	 if (ele_comp_set) anal_Ofile << ele_fit_string << " ";
	 if (energy_set) anal_Ofile << energy_value << " ";
	 if (unlike_distavg_set) anal_Ofile << unlike_distavg << " " << unlike_cooravg << " ";
	 if (like_distavg_set) anal_Ofile << like_distavg << " " << like_cooravg << " ";
	 if (gen_distavg_set) anal_Ofile << gen_distavg << " " << gen_cooravg << " ";
	 if (sphere_fit_set) anal_Ofile << sphere_fit << " ";
	 if (box_fit_set) 
		{
		anal_Ofile << box_fit[0] << " " << box_fit[1] << " " << box_fit[2] << " ";
		}
	 if (strangeness_set)
		{
		if (is_strange)
		   {
		   anal_Ofile << "YES";
		   }
		else
		   {
		   anal_Ofile << "NO";
		   }
		}
	 anal_Ofile.precision(11);
	 if (electro_potential_set) anal_Ofile << electro_potential << " ";
	 if (energy_potential_set) anal_Ofile << energy_potential << " ";
	 anal_Ofile.precision(6);
	 anal_Ofile << space_string << end_char_string;

	 anal_Ifile.close();
	 anal_Ofile.close();
	 }
