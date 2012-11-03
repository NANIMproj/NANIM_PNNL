
#include "stdafx.h"
#include "DWN_scope.h"

//Structure for parameters of microscope simulation--

//Assignment operator---

scope_param& scope_param::operator = (const scope_param& spam)
     {
	 scope_mode = spam.scope_mode;
	 cell_x = spam.cell_x;
	 cell_y = spam.cell_y;
	 cell_z = spam.cell_z;
	 tilt_x = spam.tilt_x;
	 tilt_y = spam.tilt_y;
	 tilt_z = spam.tilt_z;
	 //scan_x_min = spam.scan_x_min;
	 scan_y_min = spam.scan_y_min;
	 scan_x_max = spam.scan_x_max;
	 scan_y_max = spam.scan_y_max;
	 raster_x = spam.raster_x;
	 raster_y = spam.raster_y;
	 window_x_size = spam.window_x_size;
	 window_y_size = spam.window_y_size;
	 image_x_res = spam.image_x_res;
	 image_y_res = spam.image_y_res;
	 slice_thickness = spam.slice_thickness;
	 slice_count = spam.slice_count;
	 voltage = spam.voltage;
	 current = spam.current;
	 brightness = spam.brightness;
	 dwell_time = spam.dwell_time;
	 source_angle = spam.source_angle;
	 source_size = spam.source_size;
	 convergence_angle = spam.convergence_angle;
	 beam_tilt_x = spam.beam_tilt_x;
	 beam_tilt_y = spam.beam_tilt_y;
	 defocus = spam.defocus;
	 detector_inner_angle = spam.detector_inner_angle;
	 detector_outer_angle = spam.detector_outer_angle;
	 detector_inner_angle2 = spam.detector_inner_angle2;
	 detector_outer_angle2 = spam.detector_outer_angle2;
	 Cs3 = spam.Cs3;
	 Cs5 = spam.Cs5;
	 Cc = spam.Cc;
	 dE = spam.dE;
	 two_astig = spam.two_astig;
	 three_astig = spam.three_astig;
	 two_astig_angle = spam.two_astig_angle;
	 three_astig_angle = spam.three_astig_angle;
	 scope_temp = spam.scope_temp;
	 TDS_runs = spam.TDS_runs;
	 slice_output_count = spam.slice_output_count;

	 return *this;
	  //Proper return for this operation.
     }

//Microscope simulation class---

//Constructors, destructors, and initialization functions---

scope_sim::scope_sim()
     {
	 Initialize();
	  //Set microscope parameters to default state.
	 }

scope_sim::scope_sim(const scope_param& pam)
     {
     spam = pam;
     }

scope_sim::~scope_sim()
     {

     }

void scope_sim::Initialize()
 //Parameters from initial test simulations of PbS nanoparticles used.
     {
	 spam.scope_mode = STEM_MODE;
	 spam.cell_x = spam.cell_y = spam.cell_z = 1;
	 spam.tilt_x = spam.tilt_y = spam.tilt_z = 0.0;
     spam.scan_x_min = spam.scan_y_min = 0.0;
	 spam.scan_x_max = spam.scan_y_max = 60.0;
	 spam.raster_x = spam.raster_y = 120;
	 spam.window_x_size = spam.window_y_size = 10.0;
	 spam.image_x_res = spam.image_y_res = 0.05;
	 spam.slice_thickness = 2.967;
	 spam.slice_count = 50;
     spam.voltage = 200.0;
	 spam.brightness = 150.0;
	 spam.dwell_time = 2.0;
	 spam.source_angle = 10.0;
	 spam.source_size = 1.0;
	  //FWHM effective gun size.
	 spam.convergence_angle = 30.0;
	 spam.beam_tilt_x = spam.beam_tilt_y = 0.0;
	 spam.defocus = -4.3;
	 spam.detector_inner_angle = 100.0;
	 spam.detector_outer_angle = 170.0;
	 spam.detector_inner_angle2 = 0.0;
	 spam.detector_outer_angle2 = 20.0;
	 spam.Cs3 = 0.005;
	 spam.Cs5 = 3.2;
	 spam.Cc = 1.1;
	 spam.dE = 0.6;
	 spam.two_astig = 10.0;
	 spam.two_astig_angle = 0.0;
	 spam.three_astig = 0.0;
	 spam.three_astig_angle = 0.0;
	 spam.scope_temp = 295.15;
	 spam.TDS_runs = 15;
	 spam.slice_output_count = 20; 
     }

//Properties setting/retrieval---

void scope_sim::Get_Parameters(scope_param& pam) const
     {
     pam = spam;  
     }

void scope_sim::Set_Parameters(const scope_param& pam)
     {
	 spam = pam;
     }

//Properties calculations---

double scope_sim::Calculate_Wavelength()
     {
	 double scope_calc_temp1 = 1.0 + (1000.0*spam.voltage/(2.0*511000.0));
	 double scope_calc_temp2 = 2.0*9.109*pow(10.0, -31.0)*
							   1.602*pow(10.0, -19.0)*
					           1000*spam.voltage*scope_calc_temp1;
	 double scope_calc_temp3 = (6.626 * pow(10.0, -34.0))/sqrt(scope_calc_temp2);
	 double wavelength = scope_calc_temp3*pow(10.0, 9.0);
	 return wavelength;
     }

double scope_sim::Rel_Energy_Spread()
     {
	 return (spam.dE/(1000.0*spam.voltage));
     }


void scope_sim::Set_Scherzer_Defocus()
     {
	 double wavelength = Calculate_Wavelength();
	 double scherzer_defocus = -1.2*sqrt(spam.Cs3*pow(10.0, 6.0)*wavelength);
	 spam.defocus = scherzer_defocus;
     }
	 
void scope_sim::Scale_Thickness(double z_height)
     {
	 double z_col = z_height * double(spam.cell_z);
	  //Get total z column length.
	 double slice_ratio = z_col / spam.slice_thickness;
	 int slice_temp = int(slice_ratio - FP_ERROR_FIX) + 1;
	  //Number of slices needed to span the whole column,
	  //without changing slice thickness.
	 spam.slice_thickness = z_col/double(slice_temp);
	  //Reduced slice thickness such that the slices exactly
	  //span the column.
	 spam.slice_count = int ( (z_col + FP_ERROR_FIX) / spam.slice_thickness);
	  //Number of slices needed to span the whole column, with new slice thickness.
     }


int scope_sim::Probe_Array_Size_X()
     {
	 return int( (spam.window_x_size + FP_ERROR_FIX) / spam.image_x_res);
     }

int scope_sim::Probe_Array_Size_Y()
     {
	 return int( (spam.window_y_size + FP_ERROR_FIX) / spam.image_y_res);
     }
