
#include "stdafx.h"

#include "DWN_atom.h"

//Constructors and destructors---

atom::atom()
 //Basic constructor.
       {
       Set_Atom_Name(DUMMY_NAME);
       Zero_XYZ(atom_location);
	    //Default to the origin.
	   Set_Atom_Properties(0, 0.0, 0.0, -5.0);
	    /*Prevent dummy atom from having any default
	      spatial effects (i.e. other atoms should not 
	      not even see dummy atom). This is best ensured by
	      having a large negative radius value.*/
	   group_tag = -1;
	    //Initialize group tag to -1 (not in a molecule).
       }

atom::atom(const char* name, int num, double rel_mass, double charge, double radius)
       {
       Set_Atom_Name(name);
	   Zero_XYZ(atom_location);
	    //Default to the origin.
       Set_Atom_Properties(num, rel_mass, charge, radius); 
	   group_tag = -1;
	    //Initialize group tag to -1 (not in a molecule).
       }

atom::atom(const char* name, const double* loc, int num, 
		   double rel_mass, double charge, double radius)
       {
       Set_Atom_Specifics(name, loc);
       Set_Atom_Properties(num, rel_mass, charge, radius);    
	   group_tag = -1;
	    //Initialize group tag to -1 (not in a molecule).
       }


atom::atom(const char* name, const double* loc, int num, double rel_mass, 
	       double charge, double radius, int tag)
       {
       Set_Atom_Specifics(name, loc);
       Set_Atom_Properties(num, rel_mass, charge, radius);
	   Set_Atom_Group_Tag(tag);
       }

atom::atom(const atom& ze_atom)
	   {
	   *this = ze_atom;
	   }

atom::~atom()
 //Destructor.
       {
                  
       }

//Properties settings (non-spatial)---

void atom::Set_Atom_Specifics(const char* name, const double* loc)
       {
       Set_Atom_Name(name);
       Set_Atom_Location(loc);
       }

void atom::Set_Atom_Properties(int num, double rel_mass, double charge, double radius)
       {
	   Set_Atom_Number(num);
       Set_Atom_Rel_Mass(rel_mass);
	   Set_Atom_Charge(charge);
	   Set_Atom_Radius(radius);
       }

void atom::Set_Atom_Name(const char* name)
       {
       strcpy(atom_name, name);
       }

void atom::Set_Atom_Name(const char* name, bool only_first_capital)
 //E.g. If bool = TRUEV, Name of "CG1" (.pdb) becomes "C."
       {
	   if (only_first_capital == FALSEV)
	      {
		  Set_Atom_Name(name);
	      }
	   else
	      {
          String_Copy_StopAt_NotInitial_Captial(atom_name, name);
	      }
	   }

void atom::Index_Atom_Name(int index)
       {
	   Remove_Number(atom_name);
	    //Make sure no index is currently present.
	   Add_Number(atom_name, index);
       }

void atom::Set_Atom_Number(int at_num)
       {
	   atomic_number = at_num;
	   if (at_num >= MAX_ATOMIC_NUMBER)
		//Atomic numbers must be confied to the periodic table.
	      {
		  Show_Warning("ASSIGNED ATOMIC NUMBER IS TOO LARGE!");
	      }
       }

void atom::Set_Atom_Rel_Mass(double mass)
       {
	   rel_atomic_mass = mass;
       }

void atom::Set_Atom_Charge(double charge)
	   {
	   atomic_charge = charge;
	   }

void atom::Set_Atom_Radius(double radius)
       {
       atomic_radius = radius;
       }

void atom::Set_Atom_Group_Tag(int tag)
	   {
	   group_tag = tag;
	   }

//Properties setting (spatial)---
      
void atom::Set_Atom_Location(double x, double y, double z)
       {
       atom_location[0] = x;
       atom_location[1] = y;
       atom_location[2] = z;
       }

void atom::Set_Atom_Location(const double* loc)
       {
	   Set_XYZ(atom_location, loc);                  
       }

void atom::AddTo_Atom_Location(double del_x, double del_y, double del_z)      
       {
       atom_location[0] += del_x;
       atom_location[1] += del_y;
       atom_location[2] += del_z;
       }
       
void atom::AddTo_Atom_Location(const double* del)      
       {
	   Add_XYZ(atom_location, del);
       }

void atom::SubFrom_Atom_Location(double del_x, double del_y, double del_z)
       {
       atom_location[0] -= del_x;
       atom_location[1] -= del_y;
       atom_location[2] -= del_z;                  
       }

void atom::SubFrom_Atom_Location(const double* del)
       {
	   Sub_XYZ(atom_location, del);
       }

//Information retriveal---

void atom::Get_Atom_Name(char* ze_name) const
      {
      strcpy(ze_name, atom_name);                     
      }

bool atom::Same_Name(const atom& ze_atom) const
      {
      return ( strcmp(atom_name, ze_atom.atom_name) == 0 );                    
      }

int atom::Get_Atom_Number() const
      {
      return atomic_number;                     
      }

double atom::Get_Atom_Rel_Mass() const
      {
      return rel_atomic_mass;                        
      }

double atom::Get_Atom_Charge() const
      {
      return atomic_charge;
      }

double atom::Get_Atom_Radius() const
      {
      return atomic_radius;
      }

int atom::Get_Atom_Group_Tag() const
	  {
	  return group_tag;
	  }

void atom::Get_Atom_Location(double& x, double& y, double& z) const
       {
       x = atom_location[0];
       y = atom_location[1];
       z = atom_location[2];
       }

void atom::Get_Atom_Location(double* loc) const
       {
	   Set_XYZ(loc, atom_location);
       }

//Distance logic---

double atom::Get_Distance_To_Origin() const
       {
       return (Get_VecMag(atom_location));                  
       }

double atom::Get_Distance(const atom& ze_atom) const
       {
       double dist_val = Get_Dist(atom_location, ze_atom.atom_location);
       return dist_val;                        
       }

double atom::Get_Distance(const atom& ze_atom, int periodicity, 
                          const double* box_size) const
       {
       double dist_val = Get_Dist_OrthoPBC(atom_location, ze_atom.atom_location, 
                                           periodicity, box_size);   
       return dist_val;                      
       }

double atom::Get_Distance(const double* spatial_coors) const
       {
	   double dist_val = Get_Dist(atom_location, spatial_coors);
	   return dist_val;
       }

double atom::Get_Distance(const double* spatial_coors, 
	                      int periodicity, const double* box_size) const
       {
	   double dist_val = Get_Dist_OrthoPBC(atom_location, spatial_coors, 
		                                   periodicity, box_size);
	   return dist_val;
       }

bool atom::Covalent_Overlap(const atom& ze_atom) const
       {
	   double radii_sum = RADII_SCALING_FACTOR * 
		                  (atomic_radius + ze_atom.atomic_radius);
	   bool bond = Dist_In_Bounds(atom_location, ze_atom.atom_location, radii_sum);
	   return bond;
       }

bool atom::Covalent_Overlap(const atom& ze_atom, int periodicity, 
                            const double* box_size) const
	   {
	   double radii_sum = RADII_SCALING_FACTOR * 
		                     (atomic_radius + ze_atom.atomic_radius);
	   bool bond = Dist_In_Bounds(atom_location, ze_atom.atom_location, radii_sum,
		                          periodicity, box_size);
	   return bond;
       }

bool atom::Same_Molecule(const atom& ze_atom) const
	   {
	   return ( (group_tag == ze_atom.group_tag) && (group_tag != -1) );
	   }

//Specific identity functions---

bool atom::Is_Hydrogen() const
       {
       bool is_hydro = (atomic_number == 1);
       return is_hydro;
       }
        
bool atom::Is_Dummy_Atom() const
       {
       bool is_dummy = Name_Check(atom_name, const_cast<char*>(DUMMY_NAME));
       return is_dummy;             
       }
       
bool atom::Is_Invisible(bool include_hydrogens) const
       {
       bool invisi = FALSEV;
       if (!include_hydrogens && Is_Hydrogen())
          {
          invisi = TRUEV;                   
          } 
       else if (Is_Dummy_Atom())
          {
          invisi = TRUEV;                 
          }              
       return invisi;    
       }

//Operator overloads---

void atom::operator = (const atom& ze_atom)
       {
       strcpy(atom_name, ze_atom.atom_name);
	   Set_XYZ(atom_location, ze_atom.atom_location);
	   Set_Atom_Properties(ze_atom.atomic_number, ze_atom.rel_atomic_mass, 
		                   ze_atom.atomic_charge, ze_atom.atomic_radius);
	   Set_Atom_Group_Tag(ze_atom.group_tag);
       }

bool atom::operator == (const atom& ze_atom)
 //Does comparison operation by location alone!
       {
       bool same_location = Same_XYZ(atom_location, ze_atom.atom_location); 
       return same_location;    
       }

//Information output---

void atom::Print_Atom_Location() const
       {
       cout << endl << atom_name << " ";
	   Show_XYZ(atom_location, DIST_CONV);                       
       }
       
void atom::Store_Atom_Location(ofstream& out_file) const
       {
       out_file << endl << atom_name << " ";
	   Write_XYZ(out_file, atom_location, DIST_CONV);                            
       }

void atom::Store_Atom_Location_With_Charge(ofstream& out_file) const
       {
       out_file << endl << atom_name << " ";
	   Write_XYZ(out_file, atom_location, DIST_CONV);
       out_file << " " << atomic_charge; 
       }

void atom::Store_Atom_Location_With_CoreShell(ofstream& out_file) const
       {
       out_file << endl << atom_name << " core ";
	   Write_XYZ(out_file, atom_location, DIST_CONV);
	   out_file << endl << atom_name << " shel ";
	   Write_XYZ(out_file, atom_location, DIST_CONV);
       }

void atom::Store_Atom_Info(ofstream& Wfile) const
     {
     Wfile << endl << atom_name << " "; 
	 Write_XYZ(Wfile, atom_location, 1.0);
     Wfile << " " << atomic_number << " " << rel_atomic_mass 
           << " " << atomic_charge << " " << atomic_radius 
		   << " " << group_tag;
     }
       
void atom::Load_Atom_Info(ifstream& Rfile)
     {
     Rfile >> atom_name; 
	 Read_XYZ(Rfile, atom_location, 1.0);
     Rfile >> atomic_number >> rel_atomic_mass >> atomic_charge 
		   >> atomic_radius >> group_tag;  
     }