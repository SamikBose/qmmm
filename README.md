# qmmm
The codes and scripts related to QMMM calculations are placed here. This will mostly cover two sections:
i) To split up a MD generated structure file (preferentially in pdb or gro format) in QM and MM region. The default structure is rhodopsin, however, this code can be implemented for any biological system with slight modificantions.
ii) Generation of single point QM calculation input file. Right now the default input file is for Qchem 5.0 but can be modified for any software.

Important things to note:
i) H-atom is used as link atom.
ii) The nearest atom of MM region is removed to take care of the spurious over-polarization problem. The charge of that atom is distributed within three nearest MM atoms (standard procedure).

Prerequisite files: 
1. Topology file containing charge and atomtypes of the molecules of the system.
2. Structure file in pdb or gro format.
3. python compiler.

How does the codes work?
1. Generation of the txt file: 
Run the top_modify.py file. It will take the originial topology file as input (e.g. rhodopsin_charmm36.top). It will generate a modified topology file containg only the atomtypes and the charges. Be careful with the counters in the for loops and change them according to your requirement. The number of lines in the modified topology file (e.g. rhodopsin_charmm36_EM_mod.top) should match with the number lines in the structure file (em.gro). Then paste the two files (structure file and the modified topology file i.e. em.gro and rhodopsin_charmm36_EM_mod.top) from the command line in the linux terminal and store it to a txt file (i.e. em_topology_structure.txt in the example). This will be used as the input for the qmmm job set-up code.
2. Running the qmmm single point calculation code:
It takes the txt file as input. The cut-off radii of the water and ions are to be set at the beginning of the code. The output file of this code is a qmmm single point calculation of eom-ee-ccsd targetting first two excited states. The keywords of the qmmm calculations can be changed/modified at the line number 24 of the code. PLease go through the code and try to understand it line-by-line. All important lines have comments above or with them. For example the motivation of this qmmm code is to fragment the qm and mm region of a rhodopsin in water. It has the retinal-LYS complex as the photo-active chromophore site, which is required to be treated by the qm methods. Hence, in the code you will see 255LYR as the marker of the QM region.
