
#include "stdafx.h"
#include "DWN_File.h"

//Files to be used with other programs---

void Write_GULP_File(const atom_collection& atoms, const char* Wfile_name, 
                     const char* Dfile_name, bool include_hydro, bool overwrite_mass,
					 bool cat_coreshell, bool ani_coreshell)
	 //GULP structural optimization file.
     {
     ofstream out_file(Wfile_name, ios::out);
     if (out_file.is_open())
        {
		out_file.precision(7);

		//File header/keywords---

        out_file << "conp opti" << endl;

		//Periodicity vectors---

        double cell_temps[3];
        atoms.Get_Box_Vectors(cell_temps);
		Multiply_XYZ(cell_temps, 10.0);
		 //Convert nm to Angstroms.
		GULP_Ortho_Vector_Output(out_file, atoms.Get_Periodicity(), cell_temps);
		 //Output cell vectors.

		//List of atoms---

        out_file << endl << "Cartesian region 1";
		if (!cat_coreshell && !ani_coreshell)
			{
			atoms.Store_Atoms_Location(out_file, include_hydro);
			}
		else
			{
			atoms.Store_Atoms_Location_With_CoreShell(out_file, include_hydro,
				                                     cat_coreshell, ani_coreshell);
			}

		//Inclusion of simulation potential library---

		if (SIMULATION_POT_LIB)
		   {
		   out_file << endl << endl << "library NANIM.lib";
		   }

		//Output files request for dump files and .frc files---

        out_file << endl << endl << "dump " << Dfile_name;     
		char force_name[MAX_FILE_NAME_LENGTH];
		strcpy(force_name, Dfile_name);
	    Change_File_Type(force_name, "frc");
		out_file << endl << endl << "output frc " << force_name;
		out_file.close();
        }                                         
     }
     
void Write_GULP_Melt_File(const atom_collection& atoms, const char* Wfile_name, 
                          const char* Rfile_name, bool include_hydro, bool overwrite_mass,
						  bool cat_coreshell, bool ani_coreshell, const MD_param& mpam)
     {
     ofstream out_file(Wfile_name, ios::out);
     if (out_file.is_open())
        {
		out_file.precision(7);

		//File header/keywords---

        out_file << "md ";
		if (mpam.nvt_run) out_file << "conv ";
		if (mpam.npt_run) out_file << "conp ";
		if (mpam.molecule_commands) out_file << "molecule fix_molecule nointernalke" << endl;

		//Simulation title---

        out_file << "title" << endl << "NANIM Simulation" << endl << "end";

		//Periodicity vectors---

        double cell_temps[3];
        atoms.Get_Box_Vectors(cell_temps);
		Multiply_XYZ(cell_temps, 10.0);
		 //Convert nm to Angstroms.
		GULP_Ortho_Vector_Output(out_file, atoms.Get_Periodicity(), cell_temps);
		 //Output cell vectors.

		//Atomic relative mass for simulation, if requested---

		if (overwrite_mass)
		   {
		   out_file << endl << "element" << endl;
		   atoms.Store_Atom_Types_With_Masses(out_file, include_hydro);
		   out_file << "end";
		   }

		//List of atoms---

        out_file << endl << "Cartesian region 1";
		if (!cat_coreshell && !ani_coreshell)
			{
			atoms.Store_Atoms_Location(out_file, include_hydro);
			}
		else
			{
			atoms.Store_Atoms_Location_With_CoreShell(out_file, include_hydro,
				                                      cat_coreshell, ani_coreshell);
			}

		//Molecular dynamics parameters---

        out_file << endl << endl << "integrator leapfrog verlet" << endl << "ensemble "; 
		 //Integrator for MD equations.
		if (mpam.nvt_run) out_file << "nvt " << mpam.thermo_coeff << endl;
		if (mpam.npt_run) out_file << "npt " << mpam.thermo_coeff << " " << mpam.baro_coeff << endl;
		 //NVT vs. NPT along with pressure/temperature control coefficients.

		//Temperature information, which may include details on how to change temperature
		//over time (e.g. so one simulation file can do an annealing test)---

        out_file << "temperature " << mpam.initial_temp;
		double temp_change = mpam.final_temp - mpam.initial_temp;
		if (!Check_FP_Equality(temp_change, 0.0))
		 //If final temperature matches initial temperature, this can only be the case when
		 //no temperature change is desired. GULP simulation temperature information is for
		 //an initial constant-temperature simulation timeframe followed by changing temperature.
			{
		    double step_count_over_initial_temp = mpam.initial_temp_time/(mpam.time_step/1000.0);
			 //Number of time steps over the initial constant-temperature phase of the simulation.
			double step_count_over_temp_change = mpam.change_temp_time/(mpam.time_step/1000.0);
			 //Number of time steps over the dynamic temperature phase of the simulation.
			double temp_change_per_step = temp_change/step_count_over_temp_change;
			 //Value of the change in temperature in timestep (typically ~0.03 K/fs) for the
			 //dynamic temperature phase of the simulation.
			out_file << " " << temp_change_per_step << " " 
				     << int(step_count_over_temp_change + FP_ERROR_FIX) << " "
					 << int(step_count_over_initial_temp + FP_ERROR_FIX);
		    }

		//Times for equilibrium, production, timesteps, and for results being shown
		//on screen (sample)/in file (write)---

		out_file << endl << "equil " << mpam.equi_time << " ps" << endl
                 << "produ " << mpam.prod_time << " ps" << endl
                 << "timestep " << mpam.time_step << " fs" << endl
                 << "sample " << mpam.supdate_time << " ps" << endl
                 << "write_MD " << mpam.fupdate_time << " ps" << endl;

		//Simulation library inclusion---

		if (SIMULATION_POT_LIB)
		   {
		   out_file << endl << endl << "library NANIM.lib";
		   }

		//Output file requests for dump files and MD history files---
		
        out_file << endl << endl << "dump every 1 " << Rfile_name;
		char history_name[MAX_FILE_NAME_LENGTH];
		strcpy(history_name, Rfile_name);
	    Change_File_Type(history_name, "his");
		out_file << endl << endl << "output history " << history_name;
		out_file.close();
        }                                        
     }

void Write_PDB_File(const atom_collection& atoms, const char* Wfile_name, 
                    bool include_hydro)
	 {
	 char line_feed = char(10);
	  //Correct way to terminate lines in a PDB file.

	 ofstream out_file(Wfile_name, ios::out);
	 if (out_file.is_open())
		{		    
		//Get cell vectors---

		double cell_vecs[3];
		atoms.Get_Box_Vectors(cell_vecs);
		 //Get orthorhombic box sizes.
		Multiply_XYZ(cell_vecs, 10.0);
		 //Convert from nm to Angstroms.

		//Header of PDB file---

		char system_name[MAX_FILE_NAME_LENGTH];
        Copy_With_No_Extension(system_name, Wfile_name);
		Remove_Root(system_name);
         //Use the file name without the root or extension as a title name for the PDB file.
		out_file << "HEADER    DWN FILE GENERATOR PDB FILE";
		out_file << line_feed << "TITLE     " << system_name;
		out_file << line_feed << "REMARK 1   ENJOY THIS PDB FILE";

		//Crystal system information---

		out_file << line_feed << "CRYST1";
		int dimension = 0;
		int period = atoms.Get_Periodicity();
		if (period == 0)
		 //O-D PDB files have a special case for what is placed on the "CRYST" line.
			{
			Show_Statement("Note: You have a 0D-PDB File with according crystal format!");
			}
		for (dimension; dimension < period; ++dimension)
			{
			Right_Justify(out_file, cell_vecs[dimension], 3, 9);
			}
		for (dimension; dimension < 3; ++dimension)
			{
			if (period == 0)
			   {
			   Right_Justify(out_file, 1.0, 3, 9);
				//Default = 1.000 for 0-D case, according to PDB-file format. 
			   }
			else
			   {
			   Right_Justify(out_file, 0.0, 3, 9);
			   }
			}

		//out_file.precision(2);
		out_file << "  90.00  90.00  90.00 P 1           1";
		 //Periodicity vector angles and crystal space group.
		 //Default: Orthorhombic with no inherent symmetry.

		//Model line---

		out_file << line_feed << "MODEL        1";

		//Atom listing---

		double coors[3];
		char atom_name[MAX_ANAME_SIZE];
		int string_length;
		int atom_count = 0;
		int group_tag_index;
		for (int a = 0; a < atoms.Size(); ++a)
			{
			if (!atoms[a].Is_Invisible(include_hydro))
			 //Don't include invisible atoms.
			   {
			   out_file << line_feed << "HETATM";
			    //Leading atom keyword.

			   //Atom index/serial number---

			   ++atom_count;
			   Right_Justify(out_file, atom_count, 5);
				//Right-justify the atom index/serial number.

			   //Atom name---

			   out_file << ' ';
			   atoms[a].Get_Atom_Name(atom_name);
			   string_length = strlen(atom_name);
			   if (!Is_Text(atom_name[1]))
				 //Atomic symbol (i.e. without number) is to be
				 //right-justified within the first two columns
				 //of the atom name area.
			       {
				   out_file << " ";
				   ++string_length;
			       }
			   out_file << atom_name;

			   //Residue name/chain ID/sequence number---

			   Output_Spaces(out_file, 5 - string_length);
			   out_file << "AAA A";
			    //Default residue name and chain ID.
			   group_tag_index = atoms.Get_Group_Tag(a);
			   Right_Justify(out_file, group_tag_index, 4);
			   out_file << "    ";
			    //Use the group tag (i.e. molecular index) as 
				//a "residue" sequence number.

			   //Atomic location---

			   atoms[a].Get_Atom_Location(coors);
			   for (int c = 0; c < 3; ++c)
			       {
				   Right_Justify(out_file, coors[c]*10.0, 3, 8);
			       }

			   //Occupancy and temperature factor---

			   out_file << "  1.00  0.00          ";

			   //Atomic symbol---

			   atoms[a].Get_Atom_Name(atom_name);
			   if (!Is_Text(atom_name[1]))
			    //Right-justify atomic symbol.
			      {
				  out_file << " ";
			      }
			   Remove_Number(atom_name);
			   out_file << atom_name;
			    //Use index-less atom name as the atomic symbol.
			   }
			}

		//ENDMDL line----

		out_file << line_feed << "ENDMDL";
		out_file.close();
		}
	 }

void Write_LAMMPS_File(const atom_collection& atoms, const char* Wfile_name,
                       bool include_hydro)
     {
	 ofstream out_file(Wfile_name, ios::out);
	 if (out_file.is_open())
		{

			
		out_file.close();
		}
	 }

void Write_XYZ_File(const atom_collection& atoms, const char* Wfile_name, bool include_hydro)
     {
     ofstream out_file(Wfile_name, ios::out);                         
     
     if (out_file.is_open())
        {
		//Header---

        out_file << atoms.Real_Size(include_hydro) << endl;
         //Write number of atoms to be included in the XYZ file.
		char system_name[50];
        Copy_With_No_Extension(system_name, Wfile_name);
         //Use the file name without the extension as the system name.
        out_file << system_name << endl; 
         //System name. 
     
		//Atom listing---

        double coors[3];
        char temp_atom_name[MAX_ANAME_SIZE];
		out_file.precision(6);
		out_file << fixed;

		for (int a = 0; a < atoms.Size(); ++a)
            {
			if (!atoms[a].Is_Invisible(include_hydro))
			 //Don't include invisible atoms.
				{
				atoms[a].Get_Atom_Name(temp_atom_name);
				Remove_Number(temp_atom_name);
				 //No indexing of atoms in XYZ file!
				atoms[a].Get_Atom_Location(coors);   
				Multiply_XYZ(coors, 10.0);
				 //Convert to Angstroms.
				out_file << "  " << temp_atom_name << " ";
				Output_Spaces(out_file, 9 - strlen(temp_atom_name));
				for (int a = 0; a < 3; ++a)
					{
					Left_Justify(out_file, coors[a], 6, 16);
					}      
				}
            }
		out_file.close();
        }
     }

void Write_GSPE_File(const atom_collection& atoms, const char* Wfile_name, bool include_hydro)
	 {
	 ofstream out_file(Wfile_name, ios::out);                         
     
     if (out_file.is_open())
        {
		//Header information---

		out_file << "# HF/6-31G(d)" << endl << endl;
		char system_name[50];
        Copy_With_No_Extension(system_name, Wfile_name);
         //Use the file name without the extension as the test name.
        out_file << system_name << " SPE calculation" << endl << endl; 

		//Charge/multiplicity information---

		out_file << int(atoms.Return_Total_Charge() * (1 + FP_ERROR_FIX));
		 //System net charge.
		out_file << " 1" << endl;
		 //Spin multiplicity assumed to be one (singlet).

		//Atomic listing, looks just like the XYZ listing,
		//for the moment---

		out_file.precision(4);
		out_file << fixed;

        double coors[3];
        char temp_atom_name[MAX_ANAME_SIZE];
		for (int a = 0; a < atoms.Size(); ++a)
            {
			if (!atoms[a].Is_Invisible(include_hydro))
			 //Don't consider invisible atoms.
				{
				atoms[a].Get_Atom_Name(temp_atom_name);
				Remove_Number(temp_atom_name);
				 //No indexing of atoms!
				atoms[a].Get_Atom_Location(coors);
				Multiply_XYZ(coors, 10.0);
				 //Convert to Angstroms.
				out_file << temp_atom_name << " ";
				Output_Spaces(out_file, 5 - strlen(temp_atom_name));
				for (int a = 0; a < 3; ++a)
					{
					Left_Justify(out_file, coors[a], 3, 10);
					}      
				}
            }
        out_file.close();
        }
	 }

void Write_QSTEM_File(const atom_collection& atoms, const char* Wfile_name, bool include_hydro, 
	                  bool include_DWF)
     {
     ofstream out_file(Wfile_name, ios::out);                         
     
	 char line_feed = char(10);
	  //Note: Do not use carriage returns with QSTEM!

     if (out_file.is_open())
        {
		out_file.precision(7);

		//Header info of particle count and distance unit---

        out_file << "Number of particles = " << atoms.Real_Size(include_hydro);
         //Number of atoms to be used in the QSTEM image simulation.
        out_file << line_feed << "A = 1.0 Angstrom"; 
         //Unit = Angstroms.

		//QSTEM simulation box vectors---
        
		double cell_vec[3];
		atoms.Get_Box_Vectors(cell_vec);

        for (int a = 0; a < 9; ++a)
         //Output QSTEM cell (which has a periodicity much like a surface
		 //when simulating electron beam interaction).
           {
           out_file << line_feed << "H0(" << 1 + a/3 <<  "," << 1 + a%3 
			 	    << ") = "; 
		   if ( (a % 3) == (a / 3) )
			//Diagonal matrix element.
		      {
		      out_file << cell_vec[a / 3]*10.0; 
		      }
		   else
		      {
			  out_file << '0';
		      }
		   out_file << " A";    
		    //e.g. H0(1,1) = 60 A
		    //H0(1,2) = 0 A, etc.
           } 

		//No velocity line and number of parameters used per atom---
     
        out_file << line_feed << line_feed << ".NO_VELOCITY.";
         //Turns off MD simulation component.

		int entry_count = 3;
		 //Give three Cartesian coordinates per atom.
		if (include_DWF)
		 //Also include a Debye-Waller factor with each atom.
		   {
		   entry_count = 4;
		   }
        out_file << line_feed << line_feed << line_feed << "entry_count = " 
			     << entry_count << line_feed;

		//Atom listing, done with atom name-based categories
		//(i.e. like atoms are grouped together)---
      
        char** atom_names;
		Get_Memory(atom_names, MAX_ATOM_TYPES, MAX_ANAME_SIZE);
	     //List of atom names, one per atom type.
        int num_atom_types;
         //Number of atom names in the system (not exceeding the maximum value).
        atoms.Get_Name_List(atom_names, num_atom_types, MAX_ATOM_TYPES);

		const int MAX_INDEX = 100;
		double** dwf_vals;
		if (include_DWF)
		   {
		   Get_Memory(dwf_vals, MAX_ATOMIC_NUMBER, MAX_INDEX);
		    //Get the Debye-Waller factors from file.
		   Load_DW_Factors(dwf_vals);
		   }

		int atomic_number;
		int atom_index;
        double coors[3];
        bool first_find;
        char temp_atom_name[MAX_ANAME_SIZE];
        for (int name_index = 0; name_index < num_atom_types; ++name_index)
         //Output atomic information list by list, via the atomic symbol.
         //(e.g. list all the Si atoms first, then H atoms, etc.)
            {
            first_find = TRUEV;
            for (int a = 0; a < atoms.Size(); ++a)
                {
                atoms[a].Get_Atom_Name(temp_atom_name);
				atom_index = Get_Number(temp_atom_name);
                if (Name_Check(atom_names[name_index], temp_atom_name))
				 //Found an atom in the category being considered.
                   {
                   if (!atoms[a].Is_Invisible(include_hydro))
					//Don't consider invisible atoms.
                      {
                      if (first_find)
					   //Category header output.
                        {
                        out_file << int(atoms[a].Get_Atom_Rel_Mass() + FP_ERROR_FIX) 
                                 << line_feed << temp_atom_name << line_feed;
                        first_find = FALSEV;
                        }  
                      atoms[a].Get_Atom_Location(coors);                          
                      out_file << coors[0]/cell_vec[0] << " " << coors[1]/cell_vec[1] 
                               << " " << coors[2]/cell_vec[2];
					   //Orthorhombic box logic. QSTEM requires relative coordinates.
					  atomic_number = atoms[a].Get_Atom_Number();
					  if (include_DWF)
					     {
						 out_file << " " << dwf_vals[atomic_number][atom_index];
					     }
					  out_file << line_feed;
                      } 
                   }
                }
            }
        
		if (include_DWF)
		   {
		   Free_Memory(dwf_vals, MAX_ATOMIC_NUMBER);
		   }

		Free_Memory(atom_names, MAX_ATOM_TYPES);
		out_file.close();
        }
     }

void Write_QSC_File(const char* QSCfile_name, const char* CFGfile_name, scope_param& spam,
	                double col_height) 
     {
	 char line_feed = char(10);
	  //Note: Do not use carriage returns with QSTEM!

	 ofstream out_file(QSCfile_name, ios::out);
	 if (!out_file.is_open())
	    {
		Show_Warning("COULD NOT OPEN QSC OUTPUT FILE!");
		return;
	    }

	 out_file.setf(ios::fixed, ios::floatfield);
	 out_file.precision(8);
	  //Set precision settings.

	 char fold_name[40] = "C:\\QSTEM_RUN_FILES\\";
	 //Folder where atomic input file (.cfg) is located, at least by .qsc file's logic.
	  //This can be changed by loading the input file in QSTEM and resaving the .qsc file.
	 char test_file_name[MAX_FILE_NAME_LENGTH];
	 strcpy(test_file_name, CFGfile_name);
	 Remove_Root(test_file_name);
	  //Get name of the .cfg file name used by the .qsc simulation file.
	 char test_folder_name[MAX_FILE_NAME_LENGTH];
	 strcpy(test_folder_name, CFGfile_name);
	 Remove_Extension(test_folder_name);
	 Remove_Root(test_folder_name);
	  //Get full name of folder for storing simulation results.

	 //Header statement---

	 out_file << line_feed << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
	 out_file << line_feed << "% QSTEM configuration file generated by DWN File Maker!";
	 out_file << line_feed << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";

	 //Microscope operating mode---

	 out_file << line_feed << line_feed << "mode: ";
	 if (spam.scope_mode == TEM_MODE)
	    {
		out_file << "TEM";
	    }
	 else if (spam.scope_mode == CBED_MODE)
	    {
		out_file << "CBED";
	    }
	 else 
	    {
		out_file << "STEM";
	    }

	 //"Basic" simulation parameters---

	 out_file << line_feed << "print level: 2"; 
	 out_file << line_feed << "save level:  0";
	  //Output file style for QSTEM results (default = 2 and 0).
	 out_file << line_feed << "filename: \"" << test_file_name << '\"';
	  //File location for .cfg input file.
	 out_file << line_feed << "resolutionX: " << spam.image_x_res;
	 out_file << line_feed << "resolutionY: " << spam.image_y_res;
	  //Point sampling of electron beam/point-to-point resolution of electron beam in simulation.
	  //Typically <= 0.05 A.
	 out_file << line_feed << "NCELLX: " << spam.cell_x;
	 out_file << line_feed << "NCELLY: " << spam.cell_y;
	 out_file << line_feed << "NCELLZ: " << spam.cell_z;
	  //Number of unit cells to periodically repeat when creating simulation box.
	 out_file << line_feed << "v0: " << spam.voltage;
	  //Voltage of electron beam.
	 out_file << line_feed << "tds: ";
	  //Thermal diffuse scattering (thermal movement of atoms).
	 if (spam.TDS_runs != 0)
	    {
	    out_file << "yes";
	    }
	 else
	    {
		out_file << "no";
	    }
	 out_file << line_feed << "temperature: " << spam.scope_temp;
	  //Microscope temperature.

	 //Multi-slice aspects of the calculation---

	 scope_sim param_logic(spam);
	 param_logic.Scale_Thickness(col_height*10.0);
	  //Set slice thickness to be an integer divider of the total column height 
	  //(parameter = per unit cell). Assumes orthorhombic box described in Angstroms.
	  //This function also sets the slice count needed for this thickness.
	 param_logic.Get_Parameters(spam);
	 out_file << line_feed << "slice-thickness: " << spam.slice_thickness;
	  //Slice thickness for potential projection in multislice calculation.
	 out_file << line_feed << "slices: " << spam.slice_count;
	  //Number of slices to span the total column height.
	 out_file << line_feed << "center slices: yes % centering is the way to go";
	  //Try to position slices such that they are close to atomic planes.
	 out_file << line_feed << "slices between outputs: " << spam.slice_output_count;
	  //Interval of slices simulated before a new image is made from the 
	  //electron wavefunction (e.g. make image every 10 slices = 30 Angstroms).

	 out_file << line_feed << line_feed << "xOffset: 0.000";
	 out_file << line_feed << "yOffset: 0.000";
	 out_file << line_feed << "zOffset: 0.000";
	  //Potential offsets.
	 out_file << line_feed << "periodicXY: no";
	 out_file << line_feed << "periodicZ: no";
	  //Periodicity of slices. Assume none.

	 //Scanning window information---

	 if ( (spam.scope_mode == STEM_MODE) || (spam.scope_mode == CBED_MODE) )
	  //Scan window. Describes the rastering of the STEM beam probe
	  //across the specimen.
		{
		out_file << line_feed << line_feed << "%%%%%%%SCAN WINDOW%%%%%%";
		out_file << line_feed << line_feed << "scan_x_start: " << spam.scan_x_min;
		out_file << line_feed << "scan_x_stop: " << spam.scan_x_max;
		out_file << line_feed << line_feed << "scan_x_pixels: " << spam.raster_x;
		out_file << line_feed << "scan_y_start: " << spam.scan_y_min;
		out_file << line_feed << "scan_y_stop: " << spam.scan_y_max;
		out_file << line_feed << "scan_y_pixels: " << spam.raster_y;
		}

	 //Detector information---

	 if (spam.scope_mode == STEM_MODE)
	  //STEM detectors. Supports up to two detectors.
		{
		out_file << line_feed << line_feed << "%%%%%%%%STEM DETECTORS%%%%%%";
		out_file << line_feed << line_feed << "detector: " 
			     << spam.detector_inner_angle << " " 
			     << spam.detector_outer_angle << " detector1 0.000 0.000";
		 //No detector offsets used.
		if (!Check_FP_Equality(spam.detector_outer_angle2, 0.0))
		 //Second detector can be turned off by a zero outer angle value.
		   {
		   out_file << line_feed << "detector: " 
			        << spam.detector_inner_angle2 << " " 
			        << spam.detector_outer_angle2 << " detector2 0.000 0.000";
		   }
	    }

	  //Geometric information---
		
	 out_file << line_feed << line_feed << "%%%%%%GEOMETRIC PROPERTIES%%%%%%%%%";
	 out_file << line_feed << line_feed << "Crystal tilt X: " << spam.tilt_x * PI_CONST/180.0;
	 out_file << line_feed << "Crystal tilt Y: " << spam.tilt_y * PI_CONST/180.0;
	 out_file << line_feed << "Crystal tilt Z: " << spam.tilt_z * PI_CONST/180.0;
	  //Sample tilt, in radians.
	 out_file << line_feed << "Beam tilt X: " << spam.beam_tilt_x;
	 out_file << line_feed << "Beam tilt Y: " << spam.beam_tilt_y;
	  //Beam tilt, in degrees.
	 out_file << line_feed << "Tilt back: no";

	  //Imaging mode properties---

	 out_file << line_feed << line_feed << "%%%%%%IMAGING MODE PARAMETERS%%%%%%%%%";
	 out_file << line_feed << line_feed << "nx: " << param_logic.Probe_Array_Size_X();
	 out_file << line_feed << "ny: " << param_logic.Probe_Array_Size_Y();
	  //Number of points sampled in the virtual electron beam. 
	  //For example, 200 x 200 points in the beam with a spatial resolution of
	  //0.05 A corresponds to a 10 x 10 A size for the virtual probe description.
     out_file << line_feed << "Cs: " << spam.Cs3;
	 out_file << line_feed << "C5: " << spam.Cs5;
	 out_file << line_feed << "Cc: " << spam.Cc;
	 out_file << line_feed << "dV/V: " << param_logic.Rel_Energy_Spread();
	  //Aberrations.
	 out_file << line_feed << "alpha: " << spam.convergence_angle;
	 out_file << line_feed << "defocus: " << spam.defocus;
	  //Beam convergence angle and defocus.
	 out_file << line_feed << "astigmatism: " << spam.two_astig;
	 out_file << line_feed << "astigmatism angle: " << spam.two_astig_angle;
	 out_file << line_feed << "a_33: " << spam.three_astig*10.0;
	 out_file << line_feed << "phi_33: " << spam.three_astig_angle;
	  //More aberrations.
	 out_file << line_feed << line_feed << line_feed 
		      << "Source Size (diameter): " << spam.source_size;
	  //Electron gun size (FWHM).
	 out_file << line_feed << "beam current: " << spam.current;
      //Electrical current in beam.
	 out_file << line_feed << "dwell time: " << spam.dwell_time/1000.0;
	  //Amount of time spent gathering one pixel in a STEM image.
	 out_file << line_feed << "smooth: yes";
	 out_file << line_feed << "gaussian: no";

     //Potential information/Other assorted information---

	 out_file << line_feed << line_feed << "%%%%%%POTENTIAL PARAMETERS%%%%%%%%%%";
	 out_file << line_feed << line_feed << "potential3D: yes";
	 out_file << line_feed << "atom radius: 5.0";
	  //Atom radius used in QSTEM graphical display.
	 out_file << line_feed << "plot V(r)*r: no";
	 out_file << line_feed << "bandlimit f_trans: no";
	 out_file << line_feed << "save potential: no";
	 out_file << line_feed << "saves projected potential: no";
	 out_file << line_feed << "one time integration: yes";
	 out_file << line_feed << "Display Gamma: 0";
	 out_file << line_feed << "Folder: \"" << test_folder_name << '\"';
	  //Output folder for image simulation results.
	 out_file << line_feed << "Runs for averaging: " << spam.TDS_runs;
	  //Number of TDS runs averaged to get the final image.
	 out_file << line_feed << "Structure Factors: WK ";
	 out_file << line_feed << "show Probe: no";
	 out_file << line_feed << "propagation progress interval: 10";
	 out_file << line_feed << "potential progress interval: 1000";
	 out_file << line_feed << "update Web: no";
	  //Do not update my website.
	 out_file << line_feed << "Pendelloesung plot: no";
	 out_file << line_feed << "sequence: 1 1";

	 out_file.close();
     }

void Write_QPBS_File(const char* Wfile_name, const char* qstem_ext)
	 {
	 char line_feed = char(10);
     ofstream out_file(Wfile_name, ios::out);                         
     if (out_file.is_open())
		{
		char system_name[100];
		Copy_With_No_Extension(system_name, Wfile_name);
		Remove_Root(system_name);
         //Use the file name without the extension as the test name.

		char qfile_name[100];
		strcpy(qfile_name, system_name);
		Add_File_Extension(qfile_name, qstem_ext);

		out_file << "#!/bin/bash" << line_feed 
				 << "#PBS -lnodes=1:ppn=1" << line_feed
				 << "#PBS -N " << system_name << line_feed
				 << line_feed
				 << "cd $PBS_O_WORKDIR" << line_feed
				 << line_feed
				 << ". /etc/profile.d/modules.sh" << line_feed
				 << line_feed
				 << "module load qstem" << line_feed
				 << "mpirun -n 1 stem3 " << qfile_name;
		}
	 }
     
void Write_JEMS_File(const atom_collection& atoms, const char* Wfile_name, 
	                 bool include_hydro, bool include_DWF)
     {
     ofstream out_file(Wfile_name, ios::out);                         
     if (out_file.is_open())
        {
		out_file.precision(7);
        
		char system_name[100];
        Copy_With_No_Extension(system_name, Wfile_name);
         //Use the file name without the extension as the test name.
     
		//Header lines---

        out_file <<  "file|/" << Wfile_name << endl;
         //File name header.                  
        out_file << "name|" << system_name << endl;
         //System name.
        out_file << "creator|DWNBRF_FileGenerator" << endl; 
         //Creator name (David Welch-Nigel Browning-Roland Faller File Generator).
        time_t date;
        time (&date); 
        out_file << "date|" << ctime(&date) << endl;      
         //Time stamp.

		//Structural information---

		double cell_vec[3];
		atoms.Get_Box_Vectors(cell_vec);

        out_file << "system|triclinic" << endl;
        out_file << "superCell|true" << endl; 
        out_file << "HMSymbol|1|1|0|0| P 1" << endl;
         //Symmetry information.
        out_file << "rps|0| x , y , z" << endl;
        out_file << "lattice|0|" << cell_vec[0] << endl;
        out_file << "lattice|1|" << cell_vec[1] << endl;
        out_file << "lattice|2|" << cell_vec[2] << endl;
        out_file << "lattice|3|90.0" << endl;
        out_file << "lattice|4|90.0" << endl;
        out_file << "lattice|5|90.0" << endl;
         //Supercell structural information (assumes orthorhombic cell).

		//Loading of DWF factors if needed---

		const int MAX_INDEX = 100;
		double** dwf_vals;
		if (include_DWF)
		   {
		   Get_Memory(dwf_vals, MAX_ATOMIC_NUMBER, MAX_INDEX);
		    //Get the Debye-Waller factors from file.
		   Load_DW_Factors(dwf_vals);
		   }

		//Atom listing---

        double coors[3];
        char temp_atom_name[MAX_ANAME_SIZE];
        int actual_atom_count = 0;
		 //Used for indexing of atoms included in the JEMS file.
		int atomic_number;
		int atom_index;
        for (int a = 0; a < atoms.Size(); ++a)
         //Atom-by-atom information.
          {
          if (!atoms[a].Is_Invisible(include_hydro))
             {
             atoms[a].Get_Atom_Location(coors);
             atoms[a].Get_Atom_Name(temp_atom_name);
             out_file << "atom|" << actual_atom_count << '|' << temp_atom_name 
                      << ",_," << coors[0]/cell_vec[0] << ',' 
					  << coors[1]/cell_vec[1] << ',' << coors[2]/cell_vec[2] 
					  << ',';
			  //Coordinates given relative to the orthorhombic box.
			 out_file.precision(3);
			 out_file << fixed;
			 if (!include_DWF)
			  //Default line.
				{
				out_file << "0.005,1.000,0.100,Def,0" << endl;  
				}
			 else
			  //Include DWFs.
				{
				atomic_number = atoms[a].Get_Atom_Number();
			    atom_index = Get_Number(temp_atom_name);
				out_file << dwf_vals[atomic_number][atom_index] 
				         << ",1.000,0.100,Def,0" << endl;
				}
             ++actual_atom_count;
             }
          }
		if (include_DWF)
		   {
		   Free_Memory(dwf_vals, MAX_ATOMIC_NUMBER);
		   }

		out_file.close();
	    }  
     }

//Analysis files---

void Write_RAD_File(const atom_collection& atoms, const char* Wfile_name, 
                    double integration_interval)
 //g(r) = [ (#atoms/shell interval r1 to r2) / (4*PI*ravg*ravg*(r2-r1)) ] 
 //       / (Npairs/V)
     {
     const double MAX_DISTANCE = 1.0;
	  //Maximum distance analyzed in RDF = 10 Angstroms.
     int num_steps = int(MAX_DISTANCE/integration_interval) + 1;
	  //Number of integrations to perform to get the whole RDF
	  //from 0 to 10 Angstroms.

     char** atom_names;
	 Get_Memory(atom_names, MAX_ATOM_TYPES, MAX_ANAME_SIZE);
      //List of atom names, one per atom types.
     int num_atom_types;
      //Number of atom names in the system.
	 atoms.Get_Name_List(atom_names, num_atom_types, MAX_ATOM_TYPES);
      //Get a list of all names in the atomic collection (e.g. Si, N, etc.).

     char* atomA;
     char* atomB;
      //Temp pointers to the atom name pair being considered in one 
	  //RDF calculation. E.g. O-Si has name pair of O and Si.

     ofstream out_file(Wfile_name, ios::out);                         
     if (out_file.is_open())
        {        
        double* distances_tested = new double[num_steps];
         //X-axis distances to analyze in integration.
        double* RDF_values = new double[num_steps];
         //Y-axis values of the RDF at those x-axis distances.
        for (int a = 0; a < num_steps; ++a)
         //Initialize their values.
           {
           distances_tested[a] = double(a) * integration_interval;  
           RDF_values[a] = 0.0;  
           }
        for (int a = 0; a < num_atom_types; ++a)
            {
            atomA = atom_names[a];
            for (int b = a; b < num_atom_types; ++b)
			 //For each atom name pair that can be considered...
                {
                atomB = atom_names[b];
                atoms.Calculate_RDF(atomA, atomB, RDF_values, integration_interval);
                 //Calculate the RDF for the given atom pair.
                 
                out_file << "Radial distribution function: " << endl;
                out_file << endl << atomA << "-" << atomB; 
				 //Introduce the atom pair.
                for (int a = 0; a < num_steps; ++a)
                 //Output values.
                   {
                   out_file << endl << distances_tested[a] << " " << RDF_values[a];     
                   }
                out_file << endl << endl << endl;
                }
            }   

        delete[] RDF_values;
        delete[] distances_tested;
		Free_Memory(atom_names, MAX_ATOM_TYPES);
		out_file.close();
        }                    
     }

void Write_DISTCOUNT_File(const atom_collection& atoms, const char* Wfile_name,
						  double count_region, bool include_hydro)
     {
	 const double DISTANCE_RANGE_MIN = 0.0;
     const double DISTANCE_RANGE_MAX = 1.0;
      //Length of sampling region, from 0 to 10 A.
	 const double DISTANCE_SPAN = DISTANCE_RANGE_MAX - DISTANCE_RANGE_MIN;
     double bin_size = count_region;
      //Length of counting regions.
     int num_bins = int(DISTANCE_SPAN/count_region) + 1;
     
	 //Fill the bins and output the number of distance elements
	 //in each bin---

     ofstream out_file(Wfile_name, ios::out);                         
     if (out_file.is_open())
        {        
        double* distances_tested = new double[num_bins];
         //X-axis values/bin locations in distance space.
		int* counts = new int[num_bins];
		 //Y-axis values/number of elements in each bin.
        for (int a = 0; a < num_bins; ++a)
         //Initialize values.
           {
           distances_tested[a] = DISTANCE_RANGE_MIN + double(a) * bin_size;  
           }

		out_file << "Distance counts: " << endl;
		for (int a = 0; a < (num_bins - 1); ++a)
		 //Fill each bin.
		   {
		   counts[a] = atoms.Pair_Count(distances_tested[a], 
			                 distances_tested[a + 1], FALSEV, FALSEV, include_hydro);
			//Count the number of distances in the distance region represented by
		    //the bin being considered.
		   out_file << endl << " Distance Range of " << distances_tested[a] << " to " 
					<< distances_tested[a + 1] << " has " << counts[a] << " pairs.";
		   }

		//Additionally, determine the average distance in regions separated with
		//empty bins---

		bool at_zero = TRUEV;
		bool formerly_at_zero;
		double distance_start = 0.0;
		out_file << endl << endl << "Sums over distinguished counting regions: " << endl;
		for (int a = 0; a < (num_bins - 1); ++a)
		 //Get averages within ranges that are surrounded by zero counts/empty bins.
		   {
		   formerly_at_zero = at_zero;
		   at_zero = (counts[a] == 0);
		   if (!at_zero && formerly_at_zero)
			//Encountered the start of a filled-bin range.
		      {
			  distance_start = distances_tested[a];
		      }
		   if (at_zero && !formerly_at_zero)
			//Encountered the end of a filled-bin range. Do the analysis.
		      {
		      out_file << endl << "In the range of " << distance_start << " to " 
				       << distances_tested[a] << " the overall average distance is: "
				       << atoms.Distance_Average(distance_start, distances_tested[a], 
				                                 FALSEV, FALSEV, include_hydro);
		      }
		   }

		out_file.close();
        delete[] counts;
        delete[] distances_tested;
        }                
     }

void Write_COORCOUNT_File(const atom_collection& atoms, const char* Wfile_name,
						  bool include_hydro)
	 {
	 const double DISTANCE_RANGE_MIN = 0.00;
     const double DISTANCE_RANGE_MAX = 1.0 + FP_ERROR_FIX;
      //Length of sampling region.
	 const double DISTANCE_SPAN = DISTANCE_RANGE_MAX - DISTANCE_RANGE_MIN;
     double step_size = 0.01;
      //Step size in coordination shell analysis used for increasing
	  //the coordination sphere cut-off distance.
     int num_steps = int(DISTANCE_SPAN/step_size);
     double radius_cutoff;
	  //Cutoff radius of the coordination sphere.

	 const int MAX_COOR = 1000;
	 int coordination_counts[MAX_COOR];
	  //Bins for different coordination numbers.
	 bool is_first;

     ofstream out_file(Wfile_name, ios::out);                         
     if (out_file.is_open())
        {        
		out_file << "Coordination counts (maximum = 1000): " << endl << endl
				 << "Total atom count: " << atoms.Real_Size(include_hydro) << endl;
		for (int a = 0; a < num_steps; ++a)
		   {
		   radius_cutoff = double(a + 1)*step_size + DISTANCE_RANGE_MIN;
		   atoms.Coordination_Count(DISTANCE_RANGE_MIN, radius_cutoff, 
	                                FALSEV, FALSEV, include_hydro, 
									coordination_counts, MAX_COOR);
		   is_first = TRUEV;
		   for (int b = 1; b < MAX_COOR; ++b)
		       {
			   if (coordination_counts[b] != 0)
				  {
				  if (is_first)
				     {
				     out_file << endl << "For coordination sphere of :"
						  << DISTANCE_RANGE_MIN << " to " << radius_cutoff
						  << ": ";  
					 is_first = FALSEV;
				     }
				  out_file << endl << "CN " << b << " count: " << coordination_counts[b];
			      }
		       }
		   }
		out_file.close();
        }  
	 }


void Write_XRD_File(const atom_collection& atoms, const char* Wfile_name, 
					double wavelength, bool ionic_scat, double precision,
					double min_angle, double max_angle, double dw_factor,
					bool xrd_surface, double surf_fact, double coor_cutoff,
					int coor_num)
     {
	 double angle_range = max_angle - min_angle;
	  //Span of the x-axis of the XRD pattern.
	 int num_points = int((angle_range + FP_ERROR_FIX)/precision);
	  //Number of angles to calculate the XRD intensity for.
	 double* xrd_data = new double[num_points];
	 double angle;
	 const double SCALE = 1.000;
	  //Arbitrary scaling of output values.

	 const double RAD_TO_2THETA = 2.0*180.0/PI_CONST;
	  //Converts radians to 2*degrees.

	 double** scattering_coefficients;
	 Get_Memory(scattering_coefficients, MAX_ATOMIC_NUMBER, 8);
	 Load_Scattering_Coefficients(scattering_coefficients, ionic_scat); 
	  //Scattering coefficients to be used in calculation of atomic scattering factors.

	 atoms.Get_Debye_Intensity(xrd_data, min_angle, precision, num_points, wavelength,
		                       const_cast<const double**>(scattering_coefficients), dw_factor,
							   xrd_surface, surf_fact, coor_cutoff, coor_num);
	  //Calculate the XRD pattern.

	 ofstream out_file(Wfile_name, ios::out);
     if (out_file.is_open())
	    {
	    out_file << "X-RAY DIFFRACTION PATTERN!" << endl << endl;
		angle = min_angle;
		out_file.precision(4);
		out_file << fixed;
		for (int a = 0; a < num_points; ++a)
		 //Output the XRD pattern in terms of angles and corresponding intensities.
		    {
			angle += precision;
			out_file << endl << "2*Angle " << angle*RAD_TO_2THETA << " : " 
				     << xrd_data[a]*SCALE;
		    }

		//Attempt peak analysis to include in the data file---

		const int MAX_PEAKS = 100;
		 //Maximum number of peaks to be found in the XRD pattern.
	    int num_peaks;
		double** peaks; 
		Get_Memory(peaks, MAX_PEAKS, 4);
		Analyze_Peaks(xrd_data, num_points, min_angle + precision, precision, 
			          peaks, num_peaks, MAX_PEAKS);
		 //Characterize the peaks in the XRD pattern.

		for (int a = 0; a < num_peaks; ++a)
		    {
		    out_file << endl << endl << "Peak # " << a;
			 //Peak label.
			out_file << endl << "2*Angle (at Maximum): " << peaks[a][0]*RAD_TO_2THETA;
			 //Angle at peak maximum.
			out_file << endl << "Intensity at Max: " << peaks[a][1]*SCALE;
			 //Peak intensity.
			out_file << endl << "Integrated Intensity: " 
				             << peaks[a][2]*SCALE*RAD_TO_2THETA;
			 //Integrated peak intensity.
			out_file << endl << "Full Width at Half Maximum (2*Angle): " 
				             << peaks[a][3]*RAD_TO_2THETA;
			 //FWHM of the peak.
		    }

		
	    Free_Memory(peaks, MAX_PEAKS);
	    out_file.close();
	    }
	 else
	    {
		Show_Warning("COULD NOT OPEN XRD FILE FOR WRITING!");
	    }

	 Free_Memory(scattering_coefficients, MAX_ATOMIC_NUMBER);
	 delete[] xrd_data;
     }

void Write_DWF_File(const atom_collection& atoms, const char* Wfile_name,
	                const char* hist_name, bool include_hydro,
					double cn_min, double cn_max)
     {
	 ifstream in_file(hist_name, ios::in);
	  //History file to be read.
	 ofstream out_file(Wfile_name, ios::out);
	  //DWF file for analysis.
	 if (in_file.is_open() && out_file.is_open())
	    {
		int atom_count, temp_integer_read;
		 //Number of atoms in system and temp space for integer in input.
		double temp_fp_read;
		 //Temp space for unwanted floating-point values in input.
		char dump_string[100];
		 //Temp space for strings in the input.
		int attempt_count = 0;
		 //Number of attempts in searching for the keyword "timestep"
		 //in the history file, which precedes all other information
		 //for an instance in simulation history.
		const int MAX_ATTEMPTS = 1000;
		 //Maximum number of string reads while searching for "timestep"
		 //keyword before concluding that the file is not in proper format.

		double temp_coors[3];
		double temp_displacement[3];
		 //Temporary spatial components.
		double* rmsd; 
		double* rmsd_x;
		double* rmsd_y;
		double* rmsd_z;
		 //Running sums of spatial deviations (root-mean-square) 
		 //for each atom in the system. Memory allocation comes later.
	
		double** origins;
		 //Average position over time for each atom in the system.
		char** atom_names;
		 //Name of each atom in the system.
		bool first_step = TRUEV;
		 //Boolean that indicates if the timestep being considered is
		 //the first one.
		int num_steps = 0;
		 //Number of timesteps/instances considered.
		double square_displacement_sum;
		 //Temp square displacement sum.

		int dimensionality;
		 //Dimensionality of system, read from history file.
		double cell_vecs[9];
		 //Cell vectors of system.

		//Go to the first time step----

		while (!Command_Check(dump_string, "timestep"))
		     {
			 in_file >> dump_string;
			 ++attempt_count;
			 if (attempt_count > MAX_ATTEMPTS)
			    {
				Show_Warning("HISTROY FILE NOT IN THE RIGHT FORMAT!");
			    break;
			    }
		     }
		attempt_count = 0;

		//Process spatial coordinates of atoms at each time step
		//to determine the average position of those atoms----

		while (Command_Check(dump_string, "timestep"))
		 //Loop through all the time steps.
		     {
			 in_file >> temp_integer_read >> atom_count;
		      //Atom count is after step count in the trajectory file.


			 //Memory allocation and initialization---

			 if (first_step)
			  {
			  rmsd = new double[atom_count];
		      rmsd_x = new double[atom_count];
		      rmsd_y = new double[atom_count];
		      rmsd_z = new double[atom_count];
		      Get_Memory(origins, atom_count, 3);	
		      Get_Memory(atom_names, atom_count, MAX_ANAME_SIZE);
		      for (int a = 0; a < atom_count; ++a)
			     {
			     Zero_XYZ(origins[a]);
			     rmsd[a] = 0.0;
			     rmsd_x[a] = 0.0;
			     rmsd_y[a] = 0.0;
			     rmsd_z[a] = 0.0;
			     atom_names[a][0] = '\0';
			     }
		      first_step = FALSEV;
			   //Indicate that at least one timestep has been found.
		      }

			 in_file >> temp_integer_read >> dimensionality >> temp_fp_read;
			  //Skip past time step details aside from dimensionality
			  //of the simulation cell, which influences the next step
			  //in the read.
			 for (int a = 0; a < dimensionality; ++a)
				{
				in_file >> cell_vecs[a*3] >> cell_vecs[a*3 + 1] 
				        >> cell_vecs[a*3 + 2];
			    }

		     for (int a = 0; a < atom_count; ++a)
		        {
		        in_file >> atom_names[a] >> temp_integer_read >> temp_fp_read >> temp_fp_read;
			     //Get atom name. Skip atom index, atomic mass, and atomic charge.
			    in_file >> temp_coors[0] >> temp_coors[1] >> temp_coors[2];
			     //Get positions.
			    in_file >> temp_fp_read >> temp_fp_read >> temp_fp_read;
			     //Skip velocities.

				Add_XYZ(origins[a], temp_coors);
				 //Add the coordinates to the running position sum of the atom.
			    }

			 ++num_steps;
			  //One more timestep read.
			 in_file >> dump_string;
			  //Get the next string read, which should be the "timestep" command.
		     }

		//Turn the location sum into the location average (divide by N timesteps)---

		for (int a = 0; a < atom_count; ++a)
		    {
			Divide_XYZ(origins[a], origins[a], double(num_steps));
		    }

		//Get the time-average of the spatial-average---

		double average_origin[3];
		Zero_XYZ(average_origin);
		for (int a = 0; a < atom_count; ++a)
	 	    {
			Add_XYZ(average_origin, origins[a]);
		    }
		Divide_XYZ(average_origin, average_origin, double(atom_count));

		//Return to the file start for round two---

		in_file.clear();
		in_file.seekg(0, ios::beg);

		//FIND THE FIRST TIME STEP AGAIN---

		while (!Command_Check(dump_string, "timestep"))
		 //Go to first time step.
		     {
			 in_file >> dump_string;
			 ++attempt_count;
			 if (attempt_count > MAX_ATTEMPTS)
			    {
				Show_Warning("HISTORY FILE NOT IN THE RIGHT FORMAT!");
			    break;
			    }
		     }

		//Process statistical deviations of atoms from their average positions---

		while (Command_Check(dump_string, "timestep"))
		 //Loop through all the time steps.
		     {
			 in_file >> temp_integer_read >> atom_count;
		     in_file >> temp_integer_read >> dimensionality >> temp_fp_read;
			 for (int a = 0; a < dimensionality; ++a)
				{
				in_file >> cell_vecs[a*3] >> cell_vecs[a*3 + 1] 
				        >> cell_vecs[a*3 + 2];
			    }

		     for (int a = 0; a < atom_count; ++a)
		        {
		        in_file >> dump_string >> temp_integer_read >> temp_fp_read >> temp_fp_read;
			     //Skip atom index, atomic mass, and atomic charge.
			    in_file >> temp_coors[0] >> temp_coors[1] >> temp_coors[2];
			     //Get positions.
			    in_file >> temp_fp_read >> temp_fp_read >> temp_fp_read;
			     //Skip velocities.

				//Get running average of square devations of atoms
				//from their average positions---

				SetSub_XYZ(temp_displacement, temp_coors, origins[a]);
				square_displacement_sum = Get_CoorSquare_Sum(temp_displacement);
				rmsd[a] += square_displacement_sum/3.0;
				 //Average 1-D square displacement.
				rmsd_x[a] += temp_displacement[0]*temp_displacement[0];
				rmsd_y[a] += temp_displacement[1]*temp_displacement[1];
				rmsd_z[a] += temp_displacement[2]*temp_displacement[2];
			    }
			 in_file >> dump_string;
			  //Move to next "timestep" command.
		     }

		//Get RMSDs from running square displacement averages---

		for (int a = 0; a < atom_count; ++a)
		     {
		     rmsd[a] = sqrt(rmsd[a]/double(num_steps));
			 rmsd_x[a] = sqrt(rmsd_x[a]/double(num_steps));
			 rmsd_y[a] = sqrt(rmsd_y[a]/double(num_steps));
			 rmsd_z[a] = sqrt(rmsd_z[a]/double(num_steps));
		     } 

		//Coordination number analysis---

		int* cn_vals;
		bool do_cn_calc = TRUEV;
		cn_min *= 10.0;
		cn_max *= 10.0;
		 //Convert nm to Angstroms.
		if (Check_FP_Equality(cn_max, 0.0))
		 //Coordination analysis is turned off by the assignment
		 //of a zero value to the coordination cut-off distance.
		    {
			do_cn_calc = FALSEV;
		    }
		if (!first_step && do_cn_calc)
		    {
			cn_vals = new int[atom_count];
			Get_Coord_Numbers(const_cast<const double**>(origins), atom_count, 
				              cn_vals, cn_min, cn_max);
		    }

		//OUTPUT THE INFORMATION---

		out_file << "Debye Waller Parameters Per Atom From History" 
			     << endl << endl
				 << "Atom_Index Atom_Name RMSD_1D RMSD_X RMSD_Y RMSD_Z "
				 << "B_Factor Dist_From_Origin";
		if (do_cn_calc)
		   {
		   out_file << " CN_from_" << cn_min << "_to_" << cn_max;
		   }
		out_file << endl << endl;
		double distance_temp;
		if (!first_step)
		   {
		   for (int a = 0; a < atom_count; ++a)
		      {
			  distance_temp = Get_Dist(average_origin, origins[a]);
			  out_file << a << " " << atom_names[a] << " " << rmsd[a] << " " 
			           << rmsd_x[a] << " " << rmsd_y[a] << " " << rmsd_z[a] << " "
				       << 8*PI_CONST*PI_CONST*rmsd[a]*rmsd[a] 
					   << " " << distance_temp;
			  if (do_cn_calc)
			     {
				 out_file << " " << cn_vals[a];
			     }
			  out_file << endl;
		      }

		   //Free memory----
		   
		   if (do_cn_calc)
		      {
		      delete[] cn_vals;
		      }
		   Free_Memory(atom_names, atom_count);
		   Free_Memory(origins, atom_count);
		   delete[] rmsd_z;
		   delete[] rmsd_y;
		   delete[] rmsd_x;
		   delete[] rmsd;
		   }
	    }
	 else
	    {
		Show_Warning("FILE COULD NOT BE OPENED FOR WRITING DWF ANALYSIS!");
	    }
	 
	 in_file.close();
	 out_file.close();
	 }

int compare (const void * a, const void * b)
 //Used by quick sort in the Interpret_Image(...) function for sorting image pixel contrast values during analysis.
{
return ( *(int*)a - *(int*)b );
}

void Interpret_Image(const char* file_pointer, double source_size, double stat_pam, 
					 int restriction_pam, int* line_scan_params, bool show_nums)
	 {
	 ifstream image_file(file_pointer, ios::in | ios::binary);
	 int header_size, param_size, comment_size, Ny, Nx, complex_flag, data_size, version_num;

	 if (!image_file.is_open())
		{
		Show_Warning("COULD NOT OPEN IMAGE FILE FOR READING!");
		return;
		}

	 //Read the header information---

	 image_file.read(reinterpret_cast<char*>(&header_size), 4);
	 image_file.read(reinterpret_cast<char*>(&param_size), 4);
	 param_size += 3;
	  //Thickness, resolution values not included in sum in file.
	 image_file.read(reinterpret_cast<char*>(&comment_size), 4);
	 image_file.read(reinterpret_cast<char*>(&Ny), 4);
	 image_file.read(reinterpret_cast<char*>(&Nx), 4);
	 image_file.read(reinterpret_cast<char*>(&complex_flag), 4);
	 image_file.read(reinterpret_cast<char*>(&data_size), 4);
	 image_file.read(reinterpret_cast<char*>(&version_num), 4);

	 //Skip past parameter space---

	 double* param_read = new double[param_size];
	 image_file.read(reinterpret_cast<char*>(param_read), param_size*8);

	 //Read the comment space---

	 char* comment_read = new char[comment_size + 1];
	 image_file.read(comment_read, comment_size);
	 comment_read[comment_size] = '\0';

	 //Read the data---

	 float* data_buffer = new float[Ny*Nx];
	 for (int a = 0; a < Ny*Nx; ++a)
		 {
		 image_file.read(reinterpret_cast<char*>(&(data_buffer[a])), 4);
		 }
	 if (!image_file)
		{
		Show_Warning("IMAGE FILE READING FAILED!");
		}

	 //Output analysis---

	 char anal_name[MAX_FILE_NAME_LENGTH];
	 anal_name[0] = '\0';
	 strcpy(anal_name, file_pointer);
	 Change_File_Type(anal_name, "txt");
	 ofstream analysis_file(anal_name, ios::out);
	 if (!analysis_file.is_open())
		{
		Show_Warning("COULD NOT OPEN ANALYSIS FILE DURING IMAGE INTERPRETATION!");
		return;
		}	


	 //Statistical analysis---

	 float* data_buffer_temp;
	 if (stat_pam != -1)
		{
	    double pixel_average = 0.0;
	    int number_handled = 0;
		data_buffer_temp = new float[Nx*Ny];
		int temp_x_pixel, temp_y_pixel;
		bool low_x_test, high_x_test, low_y_test, high_y_test;
		for (int a = 0; a < Nx*Ny; ++a)
			{
			temp_x_pixel = (a/Ny);
			temp_y_pixel = (a%Ny);
			low_x_test = (temp_x_pixel >= restriction_pam);
			high_x_test = (temp_x_pixel < (Nx - restriction_pam) );
			low_y_test = (temp_y_pixel >= restriction_pam);
			high_y_test = (temp_y_pixel < (Ny - restriction_pam) );
			if (low_x_test && low_y_test && high_x_test && high_y_test)
			    {
				data_buffer_temp[number_handled] = data_buffer[a];
				pixel_average += data_buffer[a];
				++number_handled;
			    }
			}
		pixel_average /= float(number_handled);
		Show_Statement("Number of pixels analyzed in general statistics: ", number_handled);
		qsort (data_buffer_temp, number_handled, sizeof(float), compare);

		if (stat_pam >= number_handled)
		   {
		   stat_pam = number_handled - 1;
		   }

		analysis_file << endl << "General Statistics: (Restriction Pam = " << restriction_pam << ')' << endl;

		double average_max = 0.0;
		analysis_file << endl << endl << "Maximums:";

		for (int a = (number_handled - 1); a > (number_handled - 1) - stat_pam; --a)
			{
		    analysis_file << ' ' << data_buffer_temp[a];
			average_max += data_buffer_temp[a];
			}
		analysis_file << " Average: " << average_max/double(stat_pam);

		double average_min = 0.0;
		analysis_file << endl << endl << "Minimums:";
		for (int a = 0; a < stat_pam; ++a)
			{
		    analysis_file << ' ' << data_buffer_temp[a];
			average_min += data_buffer_temp[a];
			}
		analysis_file << " Average: " << average_min/double(stat_pam);

		analysis_file << endl << endl << "Overall Average: " 
					  << pixel_average << endl;
		}

	 //Linescan---

	 if (line_scan_params[2] > 0)
		{
		int del_x = 0;
		int del_y = 0;
		int integrate_x = 0;
		int integrate_y = 0;
		if (line_scan_params[1] < line_scan_params[0])
         //Y-axis scan.
		   {
		   del_y = 1;
		   integrate_x = 1;
		   }
		else if (line_scan_params[1] > line_scan_params[0])
	     //X-axis scan.
		   {
		   del_x = 1;
		   integrate_y = 1;
		   }
		else
	     //Diagonal scan for param 1 = param 0.
		   {
		   del_x = 1;
		   del_y = 1;
		   integrate_x = 1;
		   integrate_y = -1;
		   }

		int cur_x = line_scan_params[0];
		int cur_y = line_scan_params[1];
		 //Start at requested position.
		int new_x, new_y;
		int data_index;
		double point_value;

		int integrated_scan_shift = 1 + (line_scan_params[2])/2;
		int start_shift = 1 - integrated_scan_shift;
		int end_shift = line_scan_params[2] - integrated_scan_shift;
		 //Pixel range used for integrated line scan, relative to the
		 //central line of pixels considered.
		analysis_file << endl << "Line scan: (Integration pam = " << line_scan_params[2] << ')' << endl;
		while ( (cur_x < Nx) && (cur_y < Ny) )
			  {
			  point_value = 0.0;
			  for (int a = start_shift; a <= end_shift; ++a)
				  {
				  new_y = cur_y + a*integrate_y;
				  new_x = cur_x + a*integrate_x;
				  if ( (new_x >= 0) && (new_x < Nx) && (new_y >= 0) && (new_y < Ny) )
					 {
				     data_index = new_y + new_x*Ny;
					 point_value += data_buffer[data_index];
					 }
				  }
			  analysis_file << endl << cur_x*param_read[1] << " " << cur_y*param_read[2] << " : " 
							<< point_value/line_scan_params[2];
			  cur_x += del_x;
			  cur_y += del_y;
			  }
		Show_Statement("Done performing line-scan.");
		}

	 //Show numbers---

	 if (show_nums)
		{
	    analysis_file << endl << endl << "List of pixels: " << endl;
		for (int a = 0; a < Nx*Ny; ++a)
			{
			analysis_file << endl << (a/Ny)*param_read[1] << " " << (a%Ny)*param_read[2] << " : "
						  << data_buffer[a];
			}
		}

     if (stat_pam != -1)
		 {
		 delete[] data_buffer_temp;
		 }

	 analysis_file.close();
	 delete[] data_buffer;
	 delete[] comment_read;
	 delete[] param_read;
	 }

//File conversion---     

//GULP conversion---

void Convert_GULP_To_Atom_Collection(atom_collection& atom_list, double* vectors, const char* gulp_file_name, 
	                                 const atom_collection& atom_types, bool preserve_name)
     {

	 //Commands of interest within the GULP file---

     char atom_region_signature_ABS[6] = "cart";
	  //Cartesian absolute coordinates.
     char atom_region_signature_REL[6] = "frac";
	  //Relative coordinates.
	 char atom_core_signature[6] = "core";
	 char atom_shel_signature[6] = "shel";
	  //Core-shell terms.
     char vectors1_signature[6] = "pvec";
	  //1-D cell.
     char vectors2_signature[6] = "svec";
	  //2-D cell.
     char vectors3_signature[6] = "vect";
	  //3-D cell.

	 //File reading booleans---

     bool atom_abs_on = FALSEV; 
	  //Absolute coordinates.
     bool atom_rel_on = FALSEV;
	  //Relative coordinates.
     bool done = FALSEV;
	  //Done file reading.
     
	 //Lattice vectors---
	 
	 atom_list.Set_Periodicity(0);
	 Zero_Array(vectors, 9);

	 //Reading temporaries for atoms---

	 atom temp_atom;
     int temp_name_index;
	  //Temporary index.
     char temp_atom_name[MAX_ANAME_SIZE];
	  //Temporary atom name.
     char temp_coreshell_desc[8];
	  //Temporary core/shell description.
     double temp_coors[3], temp_coorsB[3];
	  //Temporary spatial coordinates.


	 char dump_string[100];
	  //Temporary storage space for reading.
     ifstream gulp_reader(gulp_file_name, ios::in);
     if (gulp_reader.is_open())
        {
        while (!gulp_reader.eof() && !done)
         //Get lattice vectors first and then loop through atom information.
              {
              gulp_reader >> dump_string;
              if (Command_Check(dump_string, vectors1_signature))
                 {
                 gulp_reader >> vectors[0];
				 atom_list.Set_Periodicity(1);
                 }
              
              if (Command_Check(dump_string, vectors2_signature))
                 {
				 gulp_reader >> vectors[0] >> vectors[1]
				             >> vectors[3] >> vectors[4];	 
				 atom_list.Set_Periodicity(2);
                 }         
              
              if (Command_Check(dump_string, vectors3_signature))
                 { 
                 for (int a = 0; a < 9; ++a)
                     {
                     gulp_reader >> vectors[a];
                     } 
				 atom_list.Set_Periodicity(3);
                 }

              if (Command_Check(dump_string, atom_region_signature_ABS))
               //Found the list of atoms and their coordinates.
                 {
                 atom_abs_on = TRUEV;
                 } 
              if (Command_Check(dump_string, atom_region_signature_REL))
                 {
                 atom_rel_on = TRUEV;
                 } 
              
              if (atom_abs_on || atom_rel_on)
			   //Loop through atom listing.
                 {
				 for (int a = 0; a < 9; ++a)
			       //Convert lattice vectors to nm.
					{
					vectors[a] /= 10.0;
					}

                 while (!gulp_reader.eof() && !done)
                       {  
                       gulp_reader >> temp_atom_name >> temp_coreshell_desc;
					    //Get atom name and core/shell description.
					   
					   Read_XYZ(gulp_reader, temp_coors, 1.0);
					    //Get Cartesian coordinates.
					   if (atom_abs_on)
					      {
						  Divide_XYZ(temp_coors, 10.0);
						   //Convert to nm.
					      }
					   else
					      {
						  Convert_RelCoor_To_AbsCoor(temp_coorsB, temp_coors, &vectors[0], 
							                         &vectors[3], &vectors[6]);
						  Set_XYZ(temp_coors, temp_coorsB);
					      }

					   gulp_reader.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					     //Ignore charge, occupancy, and radius, if present.
					   if (Command_Check(temp_coreshell_desc, atom_shel_signature))
						 //Don't need to store shell information, only interested in atom cores.
					      {
						  continue;
					      }

					   //Done reading. Set the corresponding atom and add it to the atomic collection---

                       temp_name_index = atom_types.Find_Atom_Name(temp_atom_name, FALSEV);              
                       if (temp_name_index != -1)
                        //Found a defined atom.
                          {
						  temp_atom = atom_types[temp_name_index];
						  if (preserve_name)
						     {
						     temp_atom.Set_Atom_Name(temp_atom_name);
						     }
                          }         
                       else
					    //Quit if an unfamiliar name is encountered.
                          {
                          done = TRUEV;
                          continue;
                          }   
                       temp_atom.Set_Atom_Location(temp_coors);
                       atom_list.Add_Atom(temp_atom);
                       }             
                 }
              }
		gulp_reader.close();
	    }
     }

void Convert_GULP_To_Composite(const char* gulp_file_name, const char* comp_file_name, 
                               const atom_collection& atom_types, bool preserve_name)
     {
	 atom_collection atom_list(INIT_CNUM_ATOMS);
	 double vectors[9];
     Convert_GULP_To_Atom_Collection(atom_list, vectors, gulp_file_name, atom_types, preserve_name);
     Show_Statement("Converted a GULP file with: ", atom_list.Size(), " atoms");
	 Convert_Collection_To_Composite(atom_list, comp_file_name);      
     }

void Convert_GULP_To_Structure(const char* gulp_file_name, const char* cryst_file_name, 
                               const atom_collection& atom_types, bool preserve_name)
     {
	 atom_collection atom_list(INIT_CNUM_ATOMS);
	 double vectors[9];
     Convert_GULP_To_Atom_Collection(atom_list, vectors, gulp_file_name, atom_types, preserve_name);
     Show_Statement("Converted a GULP file with: ", atom_list.Size(), " atoms");
	 Convert_Collection_To_Structure(atom_list, cryst_file_name, vectors);
     }

//History conversion---

void Convert_HIST_To_Collection(atom_collection& atom_list, double* vectors, const char* hist_file_name,
	                            const atom_collection& atom_types, double sim_time,
								int* quilt_pams, double quilt_spacing_factor, bool preserve_name)
     {

	 //Reading signatures and booleans---

	 char timestep_signature[10] = "timestep";
	 bool reached_coor_list = FALSEV;
	  //Boolean for reaching list of atomic coordinates.
     bool done = FALSEV;
	  //Boolean for completion of file reading.

	 //Reading data temporaries---

	 int step_count;
	 int dump_integer;
     char dump_string[500];
     
	 //Vector parameters---

	 int dimensionality;
	 double temp_vectors[9];
	 Zero_Array(temp_vectors, 9);
	  //Lattice vectors for new crystal system.

	 //Reading temporaries for atom---
    
	 int atom_index;
	 atom temp_atom;
	 int atom_count, temp_name_index;
     char temp_atom_name[MAX_ANAME_SIZE];
     double temp_coors[3];
	 double temp_quilt_shift[3];
     double temp_charge, temp_mass, timestep_size, actual_time, temp_vel;

	 //Simulation time and patchworking variables---

	 double target_sim_time = sim_time;
	  //Simulation timestep at which to take spatial coordinates.
	 int patch_index = 0;
	  //Index for patchwork quilting.
	 int max_patch_count = quilt_pams[0]*quilt_pams[1]*quilt_pams[2];
	  //Maximum number of patches to make.
	 if (max_patch_count == 0)
	  //Default setting = only use one instance.
	    {
	    max_patch_count = 1;
	    }
	 double quilt_pos[3];
	 Zero_XYZ(quilt_pos);
	 double lattice_step_temp[3];
	  //Lattice stepping for making the patchwork quilt.

     ifstream hist_reader(hist_file_name, ios::in);
     if (hist_reader.is_open())
        {
        while (!hist_reader.eof() && (patch_index < max_patch_count))
		 //Get lattice vectors and atomic coordinates associated
		 //with particular time steps.
              {
			  hist_reader >> dump_string;
			  if (Command_Check(dump_string, "END_HIST"))
			   //Exit clause.
			     {
				 break;
			     }
			  if (Command_Check(dump_string, timestep_signature))
			     { 
				 hist_reader >> step_count >> atom_count >> dump_integer 
					         >> dimensionality >> timestep_size;
				 actual_time = double(step_count)*timestep_size ;
				  //
				 if (actual_time > (target_sim_time - FP_ERROR_FIX))
				  //Found the instance of interest in the history file.
				    {
					reached_coor_list = TRUEV;
				    }
			     }
              
			  if (reached_coor_list)
                 {
				 atom_list.Set_Periodicity(dimensionality);
				 if (dimensionality == 1)
				    {
				    hist_reader >> temp_vectors[0];
				    }
				 else if (dimensionality == 2)
				    {
				    hist_reader >> temp_vectors[0] >> temp_vectors[1] 
				                >> temp_vectors[3] >> temp_vectors[4];
				    }
				 else if (dimensionality == 3)
				    {
				    for (int a = 0; a < 9; ++a)
					   {
					   hist_reader >> temp_vectors[a];
					   }
				    }
				 for (int a = 0; a < 9; ++a)
				    //Convert to nm.
				     {
					 temp_vectors[a] /= 10.0;
				     }

				 if (patch_index == 0)
				  //For patchwork quilt logic, use the first set of unit cell
				  //vectors for quilting.
				     {
					 for (int a = 0; a < 9; ++a)
					     {
					     vectors[a] = temp_vectors[a];
					     }
				     }

                 for (int a = 0; a < atom_count; ++a)
				  //Read in atomic information.
                     {  
                     hist_reader >> temp_atom_name >> atom_index 
					             >> temp_mass >> temp_charge;
				     --atom_index;
					  //Get atom index from 0 to (max index - 1).
					   
					 Read_XYZ(hist_reader, temp_coors, 0.100);
					  //Read in coordinates in nm.

					 //Done reading about the atom. Set new atom properties---

                     temp_name_index = atom_types.Find_Atom_Name(temp_atom_name, FALSEV);              
                     if (temp_name_index != -1)
                        //Found a defined atom.
                          {
						  temp_atom = atom_types[temp_name_index];
						  if (preserve_name)
						     {
						     temp_atom.Set_Atom_Name(temp_atom_name);
						     }
                          }         
                       else
                          {
                          done = TRUEV;
                          continue;
                          }   
					   Add_XYZ(temp_coors, quilt_pos);
					   temp_atom.Set_Atom_Location(temp_coors);
					   atom_list.Add_Atom(temp_atom);
					   hist_reader >> temp_vel >> temp_vel >> temp_vel;
					    //Skip past velocity information.
                       }   
				 reached_coor_list = FALSEV;
				 done = FALSEV;

				 //Prepare next piece of the patchwork quilt---

				 target_sim_time += sim_time;
				 ++patch_index;
				 lattice_step_temp[0] = double( patch_index % quilt_pams[0] );
				 lattice_step_temp[1] = double( (patch_index / quilt_pams[0]) % quilt_pams[1] );
				 lattice_step_temp[2] = double( (patch_index / (quilt_pams[0]*quilt_pams[1])) 
					                  % quilt_pams[2] );
				 Set_XYZ(temp_quilt_shift, lattice_step_temp);
				 Multiply_XYZ(temp_quilt_shift, quilt_spacing_factor);
                 Convert_RelCoor_To_AbsCoor(quilt_pos, lattice_step_temp, &vectors[0], &vectors[3], &vectors[6]);
				 Add_XYZ(quilt_pos, temp_quilt_shift);
                 }
              }
	    hist_reader.close();
		}
     }

void Convert_HIST_To_Composite(const char* gulp_file_name, const char* comp_file_name, 
                               const atom_collection& atom_types, double sim_time,
							   int* quilt_pams, double spacing_pam, bool preserve_name)
     {
	 atom_collection atom_list(INIT_CNUM_ATOMS);
	 double vectors[9];
	 Convert_HIST_To_Collection(atom_list, vectors, gulp_file_name, atom_types, sim_time, 
		                        quilt_pams, spacing_pam, preserve_name);
	 Show_Statement("Converted a HIST file with: ", atom_list.Size(), " atoms");
     Convert_Collection_To_Composite(atom_list, comp_file_name);        
     }

void Convert_HIST_To_Structure(const char* gulp_file_name, const char* cryst_file_name, 
                               const atom_collection& atom_types, double sim_time,
							   int* quilt_pams, double spacing_pam, bool preserve_name)
     {
	 atom_collection atom_list(INIT_CNUM_ATOMS);
	 double vectors[9];
	 Convert_HIST_To_Collection(atom_list, vectors, gulp_file_name, atom_types, sim_time, 
		                        quilt_pams, spacing_pam, preserve_name);
	 Show_Statement("Converted a HIST file with: ", atom_list.Size(), " atoms");
     Convert_Collection_To_Structure(atom_list, cryst_file_name, vectors);    
     }

//PDB Conversion---

void Convert_PDB_To_Collection(atom_collection& atom_list, double* vectors, const char* pdb_file_name, 
	                           const atom_collection& atom_types, int model_num, bool preserve_name)
     {

	 //Reading temporaries---

	 char dump_character;
	 char dump_string[100];
	 int read_length;
	 bool done = FALSEV;

     //Crystal information---

	 int cryst_num;
	 double vecm[3];
	 double veca[3];
	 bool is_periodic = FALSEV;

	 //PDB model locator---

	 int temp_model_num;
	 bool found_model = FALSEV;
	 if (model_num == 0)
	  //Default setting. Don't need to look for model commands.
	    {
		model_num = 1;
		found_model = TRUEV;
	    }

	 //Atom information---
	
	 atom temp_atom;
	 int atom_index;
	 char temp_atom_name[MAX_ANAME_SIZE + 10];
	 char temp_atom_symbol[MAX_ANAME_SIZE];
	 double temp_coors[3];
	 double temp_dwf;
	 double temp_occupancy;
	 int temp_atom_index;

	 ifstream pdb_reader(pdb_file_name, ios::in);
	 if (!pdb_reader.is_open())
	    {
		Show_Warning("COULD NOT OPEN PDB FILE!");
		return;
	    }

	 while (!pdb_reader.eof())
	    {
	    pdb_reader >> dump_string;

		if (Command_Check(dump_string, "CAVEAT"))
		   {
		   Show_Warning("YOUR PDB FILE HAS A CAVEAT DECLARED IN IT!");
		   }

		if (Command_Check(dump_string, "CRYST"))
		   {
		   cryst_num = int(dump_string[5] - '0');
		   Read_XYZ(pdb_reader, vecm, 0.100);
		   Read_XYZ(pdb_reader, veca, PI_CONST/180.0);
			 //Get cell vector magnitude and angles.
		   if (!Check_FP_Equality(vecm[0], 0.1) 
			   || !Check_FP_Equality(vecm[1], 0.1)
			   || !Check_FP_Equality(vecm[2], 0.1) )
			 //PDB file specifies no unit cell by the condition that all
			 //cell vector magnitudes are 1.0 A and all angles are 90.0 degrees.
			 //If all magnitudes are 1.0 A, this conditional fails.
		      {
			  is_periodic = TRUEV;
		      }
		   }

		if (Command_Check(dump_string, "MODEL"))
		 //Sometimes multiple models of an atomic system are included in
		 //the PDB file, so find the requested model, which is usually
		 //the first one.
		   {
		   pdb_reader >> temp_model_num;
		   if (temp_model_num == model_num)
		      {
			  found_model = TRUEV;
		      }
		   else
		      {
			  found_model = FALSEV;
		      }
		   }

		if ( found_model &&
			( Command_Check(dump_string, "ATOM") || 
			  Command_Check(dump_string, "HETATM") ) )
		 //Found atom information to convert.
		   {

		   //Atom index information---

		   if ( (dump_string[0] == 'H') && (strlen(dump_string) > 6)  )
		    //Have to watch out for e.g. "HETATM13000" case where no space
			//separates atom command and atom index.
		      {
		   	  dump_string[5] = ' ';
              atom_index = Get_Number(dump_string); 
		      }
		   else
		      {
			  pdb_reader >> atom_index;
		      }

		   //Atom name--

		   pdb_reader >> temp_atom_name;
		   read_length = strlen(temp_atom_name);
		   if (read_length > 4)
			//Due to PDB file formatting, its possible for the components of the atom
			//description to blend together with no spaces between them. A real pain
			//in the butt for typical C++ reading logic.
		      {
			  if (Is_Number(temp_atom_name[2]))
			   //Find first number, which should be the last character of the
			   //atom name if components are blending together.
			     {
				 temp_atom_name[3] = '\0';
			     }
			  else
			     {
			     temp_atom_name[4] = '\0';
			     }
		      }

		   pdb_reader >> noskipws; 
		    //Time to move past character by character (including white space).
		   for (int a = 0; a < (16 - read_length); ++a)
			//Move past residue information to the spatial coordinates.
		      {
			  pdb_reader >> dump_character;
		      }
		    pdb_reader >> skipws;

		   //Get atomic coordinates and skip occupancy and Debye-Waller Factor.
		   //Keep atomic symbol on hand for name searching---

		   Read_XYZ(pdb_reader, temp_coors, 0.100);
		   pdb_reader >> temp_occupancy >> temp_dwf >> temp_atom_symbol;

		   if (Name_Check(temp_atom_symbol, "CA") )
			//Special case to look out for: alpha-carbon vs calcium atom. If 
			//the atomic symbol is "CA", it is definitely a calcium atom. 
		      {
			  temp_atom_name[1] = 'a';
			   //Convert string "CA" to "Ca," for calcium.
			  }

		   //Done reading. Define and add the atom---

		   temp_atom_index = atom_types.Find_Atom_Name(temp_atom_name, FALSEV);  
		   if (temp_atom_index == -1)
			//Try searching user-defined atoms' names with atomic symbol.
		      {
			  temp_atom_index = atom_types.Find_Atom_Name(temp_atom_symbol, TRUEV);
		      }

           if (temp_atom_index != -1)
              //Found a defined atom.
              {
              temp_atom = atom_types[temp_atom_index];  
			  if (preserve_name)
			     {
			     temp_atom.Set_Atom_Name(temp_atom_name);
			     }
              }         
           else
              {
			  Show_Warning("ATOM ", temp_atom_name, " CAN NOT BE READ!");
              continue;
              } 
		   temp_atom.Set_Atom_Location(temp_coors);		   
           atom_list.Add_Atom(temp_atom);
		   }

		if (Command_Check(dump_string, "END"))
	     //End of file reached.
		   {
		   break;
		   }
	    }
	 pdb_reader.close();
     }

void Convert_PDB_To_Composite(const char* pdb_file_name, const char* comp_file_name, 
	                          const atom_collection& atom_types, int model_num,
							  bool preserve_name)
     {
     atom_collection atom_list(INIT_CNUM_ATOMS);
	 double vectors[9];
	 Convert_PDB_To_Collection(atom_list, vectors, pdb_file_name, atom_types, 
		                       model_num, preserve_name);
	 Show_Statement("Converted a PDB file with: ", atom_list.Size(), " atoms");
	 Convert_Collection_To_Composite(atom_list, comp_file_name);
      }  

void Convert_PDB_To_Structure(const char* pdb_file_name, const char* cryst_file_name, 
	                          const atom_collection& atom_types, int model_num, 
							  bool preserve_name)
	 {
     atom_collection atom_list(INIT_CNUM_ATOMS);
	 double vectors[9];
	 Convert_PDB_To_Collection(atom_list, vectors, pdb_file_name, atom_types, 
		                       model_num, preserve_name);
	 Show_Statement("Converted a PDB file with: ", atom_list.Size(), " atoms");
	 Convert_Collection_To_Structure(atom_list, cryst_file_name, vectors);
     }

//Collection conversion---

void Convert_Collection_To_Composite(const atom_collection& atom_list, const char* comp_file_name)
     {
     general_composite temp_composite;
	 const double ZERO_VEC[3] = { 0.0, 0.0, 0.0 };
	 temp_composite.Resize(ZERO_VEC);
     temp_composite.Add_Atom_Collection(atom_list, ZERO_VEC, TRUEV);   
	 //double box_vec[3];
	 int periodicity;
	 //atom_list.Get_Box_Vectors(box_vec);
	 periodicity = atom_list.Get_Periodicity();
	 //temp_composite.Resize(box_vec);
	 temp_composite.Set_Periodicity(periodicity);
     temp_composite.Save_Composite(comp_file_name); 
     }

void Convert_Collection_To_Structure(const atom_collection& atom_list, const char* cryst_file_name, 
	                                 const double* vectors)
     {
     crystal_system new_cryst_struct;
	 new_cryst_struct.Set_Vectors(vectors);
	 double temp_coors[3];
	 for (int a = 0; a < atom_list.Size(); ++a)
		 {
	     atom_list[a].Get_Atom_Location(temp_coors);
		 new_cryst_struct.Add_Atom_To_Basis(atom_list[a], temp_coors);
		 }
	 new_cryst_struct.Save_Crystal_System(cryst_file_name);
     }

//Partner functions---

void GULP_Ortho_Vector_Output(ofstream& writer, int periodicity, const double* cell_vec)
	  {	  
	  //0-D = no periodicity vector output.

      if (periodicity == 1)
	   //1-D (e.g. polymer) simulation.
         {
         writer << "pvectors" << endl << cell_vec[0] << endl;                        
         }
      if (periodicity == 2)
	   //2-D (e.g. surface) simulation.
         {
         writer << "svectors" << endl 
                << cell_vec[0] << " 0.0" << endl
                << "0.0 " << cell_vec[1] << endl;
         }
      else if (periodicity == 3)
		//3-D (e.g. bulk lattice) simulation.
         {
         writer << "vectors" << endl
                << cell_vec[0] << " 0.0 0.0" << endl
                << "0.0 " << cell_vec[1] << " 0.0" << endl
                << "0.0 0.0 " << cell_vec[2];      
         }
	  }


void Load_DW_Factors(double** factors)
     {
	 const char DWF_FILE[40] = "CommandFile\\DW_Factors.txt"; 
	 char name_string[40] = "\0";
	 int atom_num, atom_index;
	 double coefficient;

	 const int MAX_INDEX = 100;
	 for (int a = 0; a < MAX_ATOMIC_NUMBER; ++a)
		 {
		 Zero_Array(factors[a], MAX_INDEX);
		 }

	 ifstream in_file(DWF_FILE, ios::in);
	 if (in_file.is_open())
	    {
		in_file >> name_string >> name_string;
		while (1)
		      {
			  in_file >> name_string >> atom_num >> coefficient;
			   //Grab atom name, atomic number, and DWF factor.
			  if (Command_Check(name_string, "END_DWF") || in_file.eof())
			     {
				 break;
			     }
			  atom_index = Get_Number(name_string);
			  factors[atom_num][atom_index] = coefficient;
		      }
		in_file.close();
	    }
	 else
	    {
		Show_Warning("COULD NOT FIND DWF FILE!");
	    }
     }


void Load_Scattering_Coefficients(double** coeffs, bool ionic_scat)
     {
	 const char SCAT_FILE[40] = "CommandFile\\Scattering_Factors.txt"; 
	 char name_string[40];
	 name_string[0] = '\0';
	 int atom_num;
	 double coefficient;

	 for (int a = 0; a < MAX_ATOMIC_NUMBER; ++a)
	  //Initialize the coefficients.
	     {
		 Zero_Array(coeffs[a], 8);
	     }

	 /*For the following logic, the scattering coefficient list in file 
	   is done atom by atom, with some atoms being followed by another line
	   with scattering coefficients for the ionic form. Hence, if ionic
	   scattering is turned on, that second line, if present, is the one
	   wanted. */

	 ifstream in_file(SCAT_FILE, ios::in);
	 bool take_info;
	 if (in_file.is_open())
	    {
		in_file >> name_string >> name_string;
		while (1)
		      {
			  in_file >> name_string >> atom_num;
			  if (Command_Check(name_string, "END_SCAT") || in_file.eof())
			     {
				 break;
			     }

			  take_info = FALSEV;
			  if (Check_FP_Equality(coeffs[atom_num][0], 0.0) || ionic_scat)
			   //This ensures that either neutral or ionic scattering factors
			   //are taken from the file, depending on requested settings.
			     {
				 take_info = TRUEV;
			     }
			  for (int a = 0; a < 8; ++a)
			   //Get the scattering coefficients.
			      {
				  in_file >> coefficient;
				  if (take_info)
				     {
                     coeffs[atom_num][a] = coefficient;
					 if (a % 2)
					    {
					    coeffs[atom_num][a] *= 0.01;
						 //Change b coefficients to be in units of nm^2.
					    }
				     }
			      }

		      }
		in_file.close();
	    }
	 else
	    {
		Show_Warning("COULD NOT FIND SCATTERING FACTORS FILE!");
	    }
     }