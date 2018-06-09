/*
 * main.cpp
 *
 *  Created on: 11.03.2014
 *      Author: zhmurov
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "pdbio.h"
#include "kabsch.h"

bool atomsComparator(PDBAtom atom1, PDBAtom atom2){
	if(strcmp(atom1.segment, atom2.segment) == 0){
		if(strncmp(atom1.name, "CA", 2) == 0 && strncmp(atom2.name, "CA", 2) == 0){
			if(atom1.resid == atom2.resid){
				return true;
			}
		}
	}
	return false;
}


int main(int argc, char *argv[]){

	FILE* vmdscript = fopen("Chains/loadall.vmd", "w");

	PDB refpdb;
	readPDB("2g33_full.vdb", &refpdb);
	int i;
	int segid;
	for(i = 0; i < refpdb.atomCount; i++){
		segid = (60*i)/refpdb.atomCount + 1;
		sprintf(refpdb.atoms[i].segment, "%c%d", refpdb.atoms[i].chain, segid);
	}

	writePDB("2g33_full_segname.vdb", &refpdb);

	PDB pdb;
	PDB rotpdb;
	char filename[1024];
	
	readPDB("D_full.pdb", &pdb);
	readPDB("D_full.pdb", &rotpdb);

	for(segid = 1; segid <= 60; segid++){

		for(i = 0; i < pdb.atomCount; i++){
			pdb.atoms[i].chain = 'D';
			sprintf(pdb.atoms[i].segment, "%c%d", pdb.atoms[i].chain, segid);
			sprintf(rotpdb.atoms[i].segment, "%c%d", rotpdb.atoms[i].chain, segid);
		}
		writePDB("temp.pdb", &refpdb);
	
		doKabsch(&refpdb, &pdb, &rotpdb, &atomsComparator);
	
		sprintf(filename, "Chains/D_full_%d.pdb", segid);
		writePDB(filename, &rotpdb);
		printf(vmdscript, "mol new %s\n", filename);
	}

	fclose(vmdscript);	
}

