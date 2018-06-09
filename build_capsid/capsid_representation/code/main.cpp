#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "pdbio.h"

int main(int argc, char *argv[]){
	PDB initialPDB;
	PDB finalPDB;
	readPDB("1qgt.pdb", &initialPDB);
	readPDB("1qgt.pdb", &finalPDB);

	FILE* vmdscript = fopen("load_capsomers.vmd", "w");

	double cmx = 0;
	double cmy = 0;
	double cmz = 0;
	char lastchain = 'A';
	char outfilename[1024];
	int counter = 1;
	finalPDB.atomCount = 0;
	int i;
	for(i = 0; i < initialPDB.atomCount; i++){
		if(lastchain != initialPDB.atoms[i].chain || i == initialPDB.atomCount - 1){
			cmx /= finalPDB.atomCount;
			cmy /= finalPDB.atomCount;
			cmz /= finalPDB.atomCount;
			//save
			sprintf(outfilename, "%c%d.pdb", lastchain, counter);
			writePDB(outfilename, &finalPDB);
			fprintf(vmdscript, "mol new %s\n", outfilename);

			fprintf(vmdscript, "mol modcolor 0 top Chain\n");

			fprintf(vmdscript, "set mol%d [atomselect top all]\n", counter);
			fprintf(vmdscript, "$mol%d moveby {%f %f %f}\n", counter, -cmx, -cmy, -cmz); //so that STRIDE will work fine
			fprintf(vmdscript, "mol modstyle 0 top NewCartoon\n");
			fprintf(vmdscript, "$mol%d moveby {%f %f %f}\n", counter, cmx, cmy, cmz);
			fprintf(vmdscript, "mol modmaterial 0 top AOChalky\n");

			/*fprintf(vmdscript, "mol addrep top\n");
			fprintf(vmdscript, "mol modstyle 1 top QuickSurf\n");
			fprintf(vmdscript, "mol modcolor 1 top ColorID 8\n");
			fprintf(vmdscript, "mol modmaterial 1 top Glass1\n");*/

			cmx = 0;
			cmy = 0;
			cmz = 0;
			
			finalPDB.atomCount = 0;
	
			lastchain = initialPDB.atoms[i].chain; 

			counter++;
		}
		cmx += initialPDB.atoms[i].x;
		cmy += initialPDB.atoms[i].y;
		cmz += initialPDB.atoms[i].z;
		
		memcpy(&finalPDB.atoms[finalPDB.atomCount], &initialPDB.atoms[i], sizeof(PDBAtom));
		finalPDB.atomCount ++;
		
	}
}

