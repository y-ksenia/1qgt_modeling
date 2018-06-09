
#include "pdbio.h"

void doKabsch(PDB* refPDB, PDB* initialPDB, PDB* finalPDB, bool (*atomsComparator)(PDBAtom atom1, PDBAtom atom2));

