
//MASTER/SPECIAL COMMANDS TEXT FILE, e.g. for super hydrogen removal---

//Review alignment logic in plane sampling function---

//Neighbor lists being used to get index list for what is relevant to
//a spatial analysis of any type. This could be used in combination with
//quilting to obtain averaging---

//LAMMPS file writing---

//Rotation logic cubic? Also, more logic for the rotation for plane-->plane
//that makes it aware of what arrangements are at in the final plane. Where
//does intraplanar configuration end up?

//Charge balancing checks for molecule removal---

//Component addition modes/spatial analysis in general in composite class
//for molecule fragmentation...Removing odd tricks?---

//NPT solvent molecule addition--adding certain numbers of molecules.

//Average XRD pattern from history file---

//Reading in history file for XRD calculation---

//Dehydration command

//Atom removal during charge balancing by coordination number

//Neighbor lists---

//Make cubic cuts with density remaining the same, e.g. for amorphous simulations---

//More detailed molecule cleaning function that analyzes sp3, etc. geometry of molecules---

//Alignment function between multiple groups---

//Adding functionality so that multiple composites can be loaded at the same time---

//Fit_Box function having more details
//i.e. no hydrogen inclusion and specific atom analysis.

//COORDINATE AVERAGE function is not periodic---
//Second function with that ability.

//Covalent and vdW radii...maybe

//GULP File reading---Conversion of "cell" commands in addition
//to "vectors" using function in crystal structure that takes
//3 lengths and 3 angles.

//Make Find_Atom_Name() more efficient with its name scanning process.

//Make box/ near-neighbor list for atom collection to make
//inverse substitution faster

//Specific strangeness parameters...other forms of strangeness...
//DONE PARTIALLY---

//Water cutting with a box creation in nanosurface definition
//Chunks missing in nanosurface definition of water box???



//Inverse substitution---removal of just oxygens and not hydrogens too??

//Bond angle distribution/full coordination number analysis 
//in ANALYSIS capabilities. RADIAL DISTRIBUTION PART.

//Add user parameter for Strangeness(...) distance criterion.

//Reposition corner function in atom collection
//And use it on Fit_Box(...) function at the start
//of the function so as to 100% prevent math mistakes.

//Orientation in plane sampling
//Other point placements for surface monolayers?
//c(4 x 2) Placement and positioning over fcc holes....
//Use of unit cell in monolayer placement?

//Spatial overlap with radii for solvent placement

//Correct other particle labels---

//Memory check, small head-head spacings lead to bad memory access

//More molecule building tools----

//Initial orientation for monolayers on surfaces (changeable)

//BCC Grid/110 faces not tested in plane sampling. Also, ability to vary: plane for sampling, lattice for points in that plane, AND rotation angle applied to the final set

//Read over the code
//Extensive testing...also memory check
//Also, convert "int" to "long int" where appropriate, also consider other data types.
//Make use of break commands
//Compartmentalize code further.


//Have supercell as an optional surface construction mode.

//1D and 2D systems with crystal structure (still use same system)

//Extending composite system to complex dimensions

//Also, split coating functions into space samplers in the utility files. Applies to surface and solvent too.



//Definitely going to be last---

//More complicated amorphous structure generation (also, 1D, 2D...)---

//Periodicity for multislice method and nanosurface cutting
 //Should probably test the significance of this first before getting worked up over it. 

//Use trajectory files to make movies--- 
 
//Potential declaration help-- 

//Other factors in atom list for QSTEM such as DW factor 

//Indexing atoms by molecule groups...
 
//Surface slab creation, try to speed up big systems. Speeding up {100} planes for cubic systems should be quick. Recursion algorithms!


//General statistics file output by program---
//DONE---

//Small_Dist(...) having PBC version.
//DONE---

//Solvent addition using covalent radii instead of distance parameter
//DONE---

//Remove reverse ordering in GULP Files. Then change reference atoms used in
//Rigid_Minimization(...) function to atoms 0 and 1 for ease of use.
//DONE---

//List of atoms with charges not equal to zero for potential/force calculation---
 //DONE-----

//Spherical/cubic voids in composite
 //DONE----

//Scripting language---combine moving pointers with command checks somehow
 //JUDGED UNNECESSARY

//MD enabler??? Position/velocity initializer???
 //JUDGED UNNECESSARY


//Function that makes sure all atomic names have a unique number added to
//the end, for input files like "Au1 ... Au2 ... Au3 ... etc.
/// DONNNNNNNEEEEE----------

//Basis atom addition to crystal structure should be able to be 
//relative to crystal vectors.
// DONNNNNNEEEEE----------

//Atom function -- Construct a list of atomic names from a collection
//for easy access.
/// DONNNNNEEEEE-----

//Relative coordinate input in conversion file.
//Convert .gin to .comp file TODAY
 //DONEE------
 
//Numbers in scripting names, be able to say "Struct Si3N4 0.0 0.0 0.0"
//DONEEE----

//Cu is intrepreted by carbon (C) in atom identification.
//DONEEE----

//Add individual atoms to the composite
//DONEEE----

//Check loading files for all system types (intro strings???)
//Make function to skip text (e.g. skip 5 words).
//DONEE-----

//Capitalize constants
//DONE----

//Make sure memory initialization/deletion is good in all systems
//DONE-----

 //Improve load-then-add functionality of general composite
 //DONE---

//Read birthday wishes on FB
//DONE WITH MUCH LOVE---

//Keep solvent filler out of the box for surface. Solvent area vs system area?
 //DONE----

//Size violation check to 3D
 //Changed to an indepedent parameter. May need FUTURE consideration---
 
//Separate utility aspects 
 //DONE----

//Also: Rotation of MgO 110?? Bonds at edges???
 //CHECK DONE---

//Also: RDF periodicity inclusion
//Periodicity utility files
//DONE---

//Big surfaces of a-Si3N4 made program fail. WHYYYY?
//Takes a lot of time to calculate and watch the maximum atom counts in program---

//RDF values seem strange around 2a
 //DONE---
 
//Molecule building tools---
 //DONE---

//Hydrogens turned off
 //DONE---

//(Keep water hydrogens turned on for MD)--- 
//Have blank "X" atom option in molecule building that allows certain atoms to 
//just not be added during molecule construction.
 //DONE WITH DUMMY ATOM XXXX----
 
 //Specify coating functions in terms of BIND TO ATOM vs IGNORING the atomic
//structure of the surface 
 //DONE----
 
//Constant-ify the code
 //DONE AFTER ~2 to 3 HOURS. HOLY CRAP! WHY DO I DO THIS TO MYSELF? 

//Check nanoparticle coating.
 //DONE. KEEP IN MIND WITH CONTINUOUS MODEL A SULFUR ATOM CAN BE VERY CLOSE
 //TO MULTIPLE GOLD ATOMS.

//Crystalline nanop model
 //DONE. Rather quickly too...so EAT THAT!
 
 //Even distribution of coating function over surface slab like what was done
//for particles per Josh's advice
 //DONE

//Different add modes for molecules--FIRST ATOM and MOLECULE CENTER  
//DONE

//Atom collection
//DONE.

//Rigid regions being outputted second (reverse order?) 
//DONE.

//Some scripting readers (such as Get_Num()) fail if there is extra characters
//or white space at the end of a line. This is weak point.
//DONE.

//Rotation of single molecule placement (see composite) 
//DONE.

//Random close-packed spheres. Research further crystalline nanoparticles and
//observed structure.
//DONE. A pretty good random algorithm employed =) 

//Memory limits
//DONE. When a program gives the weirdest of errors, there is a good chance
//of memory mismanagement/corruption somewhere. Don't forget this again!

//Make stuffed attribute overwrittable by text file
//DONE---

//Output weights for MD files
//DONE---

//Parameter for a-Au mesh?
//Found one paper. May be sufficient to use crystalline bond length
//Perhaps best to reproduce experimental density if possible.
//DONE FOR NOW---

//Composite system---add ions for total charge balance
//DONE----

//Coating surfaces and non-faceted NP. Use new spatial overlap function.
//DONE---

//Molecule rotation during solvent placement an addition??? Atom collection rotation?
//PARTIALLY DONE--- SOLVENT PLACEMENT DOES NOT NEED RANDOM PLACEMENT SO FAR

//2D memory allocation/deallocation functions
//DONE---

//Check faceted coating--
//DONE---

//180 degree rotation matrix
//DONE---

//Rotation algorithm--standard treatment for unequal vectors passed to it---
//DONE--

//Finish setting SIO orientation and see if it helps---
//DONE---

//Symmetry of coating...secondary rotation to get orientation complete.
//DONE---

//PNNL Online training
//Register the C++ Express 
//DONE---

//Continue symmetry analysis
//DONE---

//Check on that weird inverse order thing for plane rotation function calls in cryst system.
//DONE---

//Check program logic to make sure other parts of lines don't trigger unwanted results
//MPNL that goes to the first character of the next line with the char pointer
//Or use of temps
//Perhaps MPSNL functions should only move the first pointer so far...
//DONE---

//AddXXX command that adds a group of atoms while removing all in the area of that group
//DONE---

//Shells -- both for speedup and for making just a surface layer for simulation
//DONE---

//Charge balance by removing charged species (and solvent approach)
//DONE---

//Systematically vary size parameters until charge_neutrality is achieved.
//DONE---

//Also look for neutrality by shifting the locations of the point sampling.
//Or actual cleavage algorithm...
//DONE---

//Core-shell terms
//DONE---
