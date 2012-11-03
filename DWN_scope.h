
#ifndef MATT_VITAMINS
#define MATT_VITAMINS

//Microscope file constants---

const int STEM_MODE = 0;
 //Scanning transmission electron microscopy.
const int TEM_MODE = 1;
 //Transmission electron microscopy.
const int CBED_MODE = 2;
 //Convergent-beam electron diffraction.

//Structure for parameters of microscope simulation--

struct scope_param
   {
  
   //General parameters---

   int scope_mode;
    //STEM, TEM, or CBED mode?
   int cell_x, cell_y, cell_z;
    //Repetition of input structure cell to creat simulation box,
    //which is in turn periodic in the x and y directions.
   double tilt_x, tilt_y, tilt_z;
    //Sample tilting, in degrees.

   //Rastering/beam positioning parameters---

   double scan_x_min, scan_x_max;
    //Scan start and end positions for probe rastering, along x axis, in Angstroms. 
   double scan_y_min, scan_y_max;
    //Scan start and end positions for probe rastering, along y axis, in Angstroms.
   int raster_x, raster_y;
    //Number of pixels (probe positions) along x and y directions.

   //Beam representation in simualtion---

   double window_x_size, window_y_size;
    /*Electron beam size in simulation in the x and y directions, in Angstroms. 
      Beyond these sizes, any part of the beam (i.e. from beam spreading) that would be
      present in real life electron microscopy will be neglected in simulation.
      e.g. window size = 20 A --> This means to simulate the beam only out to 10 A from 
	  the beam center.*/
   double image_x_res, image_y_res;
    /*Probe resolution in the x and y directions, in Angstroms. This is essentially the 
	  spatial sampling in the simulation of the electron beam much like timestep is 
	  temporal sampling in MD. In other words, the beam can not be infinitely represented 
	  in simulation so this is the spacing between points of the electron beam that are sampled.*/

   //Multislice parameters---

   double slice_thickness;
    //Slice thickness in multislice calculation, in Angstroms.
   int slice_count;
    //Number of slices in multislice calculation.

   //Microscope parameters---

   double voltage, current, brightness, dwell_time;
    //Electron probe features. Units: keV, pA, 10^8 A/cm^2sr, microsecond.
   double source_angle, source_size, convergence_angle, beam_tilt_x, beam_tilt_y, defocus;
    //Pre-specimen spatial probe features. Units: Degrees or Angstroms.
    //Source size = effective gun size = FWHM of gun emission current.
   double detector_inner_angle, detector_outer_angle, detector_inner_angle2, detector_outer_angle2;
    //Detector(s) annular ranges for non-TEM imaging.
   double Cs3, Cs5, Cc, dE, two_astig, three_astig, two_astig_angle, three_astig_angle;
    //Aberration/lens coefficients. Units: mm, mm, mm, eV, nm, nm, degrees, degrees.
   double scope_temp;
    //Temperature of microscope/specimen in Kelvin.
   int TDS_runs;
    //Number of frozen phonon configurations to be sampled and finally averaged.
   int slice_output_count;
	//Interval of slices simulated before a new image is made from the 
	//electron wavefunction (e.g. make image every 10 slices = 30 Angstroms).

   scope_param& operator = (const scope_param&);
    //Sets all properties of one equal to the other.
   };

//Formal class for microscope simulation, which relies heavily on the above structure---

class scope_sim
   {
   public:

   //Constructors, destructors, and initialization functions---

   scope_sim();
   scope_sim(const scope_param&);
    //Initialize the group of microscope/simulation parameters.
   ~scope_sim();

   void Initialize();
    //Sets microscope simulation to default parameters.

   //Properties setting/retrieval---

   void Get_Parameters(scope_param&) const;
    //Sets all parameters of the passed microscope parameters
    //to those of this class's set of parameters.
   void Set_Parameters(const scope_param&);
    //Sets all parameters of this class's set of parameters
    //to those passed to this function.

   //Properties calculations---

   double Calculate_Wavelength();
    //Calculates relativistic wavelength (in nm) corresponding to the electron beam voltage.
   double Rel_Energy_Spread();
    //Gets energy spread (in eV/V) that is relative to the beam voltage = dE/(1000.0*voltage).
   
   void Set_Scherzer_Defocus();
    //Sets the defocus to the Scherzer defocus value (in nm) that best balances the C3 spherical
    //aberration.
   void Scale_Thickness(double);
    //Changes the multislice thickness such that the thickness evenly
    //divides the overall height of the column (in Angstroms) through which the beam
    //passes. Thus sets the multislice thickness and slice count parameters.

   int Probe_Array_Size_X();
   int Probe_Array_Size_Y();
    //Calculates the size of the probe array = beam window size / beam spatial resolution.

   private:

   scope_param spam;
    //Parameters of the electron imaging simulation.
   };

#endif