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
	char filename[1024] = "A_full.pdb";
	readPDB(filename, &pdb);
	PDB fullpdb;
	fullpdb.atoms = (PDBAtom*)calloc(pdb.atomCount*500, sizeof(PDBAtom));
	fullpdb.atomCount = 0;
	char chain;

	for(segid = 1; segid <= 60; segid++){
		for(chain = 'A'; chain <= 'D'; chain++){
			sprintf(filename, "%c_full.pdb", chain);
			readPDB(filename, &pdb);
			readPDB(filename, &rotpdb);
			for(i = 0; i < pdb.atomCount; i++){
				pdb.atoms[i].chain = chain;
				rotpdb.atoms[i].chain = chain;
				sprintf(pdb.atoms[i].segment, "%c%d", pdb.atoms[i].chain, segid);
				sprintf(rotpdb.atoms[i].segment, "%c%d", rotpdb.atoms[i].chain, segid);
			}
		
			
			doKabsch(&refpdb, &pdb, &rotpdb, &atomsComparator);

			for(i = 0; i < rotpdb.atomCount; i++){
				memcpy(&fullpdb.atoms[fullpdb.atomCount], &rotpdb.atoms[i], sizeof(PDBAtom));
				fullpdb.atomCount++;
			}
			sprintf(filename, "Chains/%c_full_%d.pdb", chain, segid);
			writePDB(filename, &rotpdb);
			fprintf(vmdscript, "mol new %s\n", filename);
		}
	}
	sprintf(filename, "Chains/Full.pdb");
	writePDB(filename, &fullpdb);
	
	fclose(vmdscript);
}

