# qmmm
The codes and scripts related to QMMM calculations will be placed here. This will mostly cover two sections:
i) To split up a MD generated structure file (preferentially in pdb or gro format) in QM and MM region. The default structure is rhodopsin, however, this code can be implemented for any biological system with slight modificantions.
ii) Generation of single point QM calculation input file. Right now the default input file is for Qchem 5.0 but can be modified for any software.

Important things to note:
i) H-atom is used as link atom.
ii) The nearest atom of MM region is removed to take care of the spurious over-polarization problem. The charge of that atom is distributed within three nearest MM atoms (standard procedure).
