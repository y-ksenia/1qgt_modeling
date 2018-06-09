#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pdbio.h"

/*
 * pdbio.c
 *
 *  Created on: Nov 9, 2008
 *      Author: zhmurov
 */

#define buf_size 80
//#define DEBUGPDBIO

/*
 * Private methods
 */
void parseAtomLine(PDB* pdbData, char* line, int currentAtom);
void parseSSBondLine(PDB* pdbData, char* line, int currentSSBond);
void parseSymmetryLine(PDB* pdbData, char* line);
void parseMatrixLine(PDB* pdbData, char* line);
void parseCrystalLine(PDB* pdbData, char* line);


/*
 * Parses data from PDB (Protein Data Bank) file format into PDB object.
 * Only ATOM and SSBOND entries are considered, some fields are ignored
 * (marked befor corresponding method).
 *
 * For detailed information about PDB file format, please visit
 * http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Parameters:
 * 		filename: name of a pdb file to parse
 * 		pdbData: pointer of an object to save data into
 */

void readPDB(const char* filename, PDB* pdbData){
	printf("Reading %s.\n", filename);
	int ss_count = 0, atoms_count = 0, symmetry_count = 0, matrix_count = 0;
	char buffer[buf_size];
	FILE* file = fopen(filename, "r");
	if ( file != NULL ){
		while(fgets(buffer, buf_size, file) != NULL){
			if(strncmp(buffer,"SSBOND",6) == 0){
				ss_count++;
			}
			if(strncmp(buffer, "ATOM", 4) == 0 || strncmp(buffer, "HETATM", 6) == 0){
				atoms_count++;
			}
			if(strncmp(buffer, "REMARK 290   SMTRY", 18) == 0){
				symmetry_count++;
			}
			if(strncmp(buffer, "MTRIX", 5) == 0){
				matrix_count++;
			}
		}
		symmetry_count /= 3;
		matrix_count /= 3;
		printf("Found: %d atoms, %d S-S bonds, %d symmetry entries, %d rotational matrices\n",
				atoms_count, ss_count, symmetry_count, matrix_count);

		pdbData->atomCount = atoms_count;
		pdbData->ssCount = ss_count;
		pdbData->symmetryCount = symmetry_count;
		pdbData->matrixCount = matrix_count;
		pdbData->atoms = (PDBAtom*)malloc(atoms_count*sizeof(PDBAtom));
		pdbData->ssbonds = (PDBSSBond*)malloc(ss_count*sizeof(PDBSSBond));
		pdbData->symmetries = (PDBSymmetry*)malloc(symmetry_count*sizeof(PDBSymmetry));
		pdbData->matrices = (PDBMatrix*)malloc(matrix_count*sizeof(PDBMatrix));

		int current_atom = 0;
		int current_ss = 0;

		rewind(file);

		while(fgets(buffer, buf_size, file) != NULL){
			//char* pch = strtok(buffer, " ");
			//printf("%s\n", buffer);
			if(strncmp(buffer,"SSBOND",6) == 0){
				parseSSBondLine(pdbData, buffer, current_ss);
				current_ss++;
			}
			if(strncmp(buffer, "ATOM", 4) == 0 || strncmp(buffer, "HETATM", 6) == 0){
				parseAtomLine(pdbData, buffer, current_atom);
				current_atom ++;
			}
			if(strncmp(&buffer[0], "REMARK 290   SMTRY", 18) == 0){
				//printf("%s\n", buffer);
				parseSymmetryLine(pdbData, buffer);
			}
			if(strncmp(buffer, "MTRIX", 5) == 0){
				parseMatrixLine(pdbData, buffer);
			}
			if(strncmp(buffer, "CRYST", 5) == 0){
				parseCrystalLine(pdbData, buffer);
			}
		}
		int i;
		/*for(i = 0; i < pdbData->symmetryCount; i++){
			printf("================\n");
			printf("%d:\n", pdbData->symmetries[i].index);
			printf("%f\t%f\t%f\t%f\n", pdbData->symmetries[i].r11, pdbData->symmetries[i].r12, pdbData->symmetries[i].r13, pdbData->symmetries[i].tx);
			printf("%f\t%f\t%f\t%f\n", pdbData->symmetries[i].r21, pdbData->symmetries[i].r22, pdbData->symmetries[i].r23, pdbData->symmetries[i].ty);
			printf("%f\t%f\t%f\t%f\n", pdbData->symmetries[i].r31, pdbData->symmetries[i].r32, pdbData->symmetries[i].r33, pdbData->symmetries[i].tz);
			printf("================\n");
		}*/
		/*for(i = 0; i < pdbData->matrixCount; i++){
			printf("================\n");
			printf("%d:\n", pdbData->matrices[i].index);
			printf("%f\t%f\t%f\t%f\n", pdbData->matrices[i].M11, pdbData->matrices[i].M12, pdbData->matrices[i].M13, pdbData->matrices[i].V1);
			printf("%f\t%f\t%f\t%f\n", pdbData->matrices[i].M21, pdbData->matrices[i].M22, pdbData->matrices[i].M23, pdbData->matrices[i].V2);
			printf("%f\t%f\t%f\t%f\n", pdbData->matrices[i].M31, pdbData->matrices[i].M32, pdbData->matrices[i].M33, pdbData->matrices[i].V3);
			printf("================\n");
		}*/
	printf("Done reading '%s'.\n", filename);
	fclose(file);
	} else {
		perror(filename);
		exit(0);
	}
}

/*
 * Parses single line of 'ATOM' entry from pdb file.
 * ATOM entry format in PDB:
 *
 * COLUMNS      DATA TYPE        FIELD      DEFINITION
 * ------------------------------------------------------
 *  1 -  6      Record name      "ATOM    "
 *  7 - 11      Integer          id		    Atom serial number.
 * 13 - 16      Atom             name       Atom name.
 * 17           Character        altLoc     Alternate location indicator.
 * 18 - 20      Residue name     resName    Residue name.
 * 22           Character        chainID    Chain identifier.
 * 23 - 26      Integer          resSeq     Residue sequence number.
 * 27           AChar            iCode      Code for insertion of residues. (Ignored)
 * 31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angstroms
 * 39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angstroms
 * 47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angstroms
 * 55 - 60      Real(6.2)        occupancy  Occupancy.
 * 61 - 66      Real(6.2)        tempFactor Temperature factor.
 * 77 - 78      LString(2)       element    Element symbol, right-justified. (Ignored)
 * 79 - 80      LString(2)       charge     Charge on the atom. (Ignored)
 *
 * For more details visit http://www.wwpdb.org/documentation/format23/sect9.html
 *
 *
 * Method parameters:
 * 		pdbData:	pointer to pdbData object to read into
 * 		line: line from the file to parse
 * 		currentSSBond: indicates location in array of ssbonds
 * 			where new bond should be saved
 */

void parseAtomLine(PDB* pdbData, char* line, int currentAtom){
	char id[6];
	char atomName[5];
	char resName[4], chain, altLoc;
	char resid[5];
	char x[9], y[9], z[9];
	char occupancy[7];
	char beta[7];
	char segment[5];
	strncpy(id, &line[6], 5);
	id[5] = '\0';
	strncpy(atomName, &line[12], 4);
	atomName[4] = '\0';
	altLoc = line[16];
	strncpy(resName, &line[17], 3);
	resName[3] = '\0';
	strncpy(&chain, &line[21], 1);
	strncpy(resid, &line[22], 4);
	resid[4] = '\0';
	strncpy(x, &line[30], 8);
	x[8] = '\0';
	strncpy(y, &line[38], 8);
	y[8] = '\0';
	strncpy(z, &line[46], 8);
	z[8] = '\0';
	strncpy(occupancy, &line[54], 6);
	occupancy[6] = '\0';
	strncpy(beta, &line[60], 6);
	beta[6] = '\0';
	strncpy(segment, &line[72], 4);
	segment[4] = '\0';
	int i;
	for(i = 0; i < 5; i++){
		if(segment[i] == '\n' || segment [i] == '\r'){
			segment[i] = ' ';
		}
	}
	strcpy(pdbData->atoms[currentAtom].name, strtok(atomName, " "));
	pdbData->atoms[currentAtom].altLoc = altLoc;
	pdbData->atoms[currentAtom].name[4] = 0;
	pdbData->atoms[currentAtom].chain = chain;
	pdbData->atoms[currentAtom].resid = atoi(resid);
	strcpy(pdbData->atoms[currentAtom].resName, resName);
	pdbData->atoms[currentAtom].resName[3] = 0;
	pdbData->atoms[currentAtom].id = atoi(id);
	pdbData->atoms[currentAtom].x = atof(x);
	pdbData->atoms[currentAtom].y = atof(y);
	pdbData->atoms[currentAtom].z = atof(z);
	pdbData->atoms[currentAtom].occupancy = atof(occupancy);
	pdbData->atoms[currentAtom].beta = atof(beta);
	strcpy(pdbData->atoms[currentAtom].segment, segment);
#ifdef DEBUGPDBIO
	printAtom(pdbData->atoms[currentAtom]);
#endif
}

/*
 * Parses single line of 'SSBOND' entry from pdb file.
 * SSBOND entry format in PDB:
 *
 * COLUMNS        DATA TYPE       FIELD         DEFINITION
 * -------------------------------------------------------------------
 *  1 -  6        Record name     "SSBOND"
 *  8 - 10        Integer         serNum       Serial number.
 * 12 - 14        LString(3)      "CYS"        Residue name.
 * 16             Character       chainID1     Chain identifier.
 * 18 - 21        Integer         seqNum1      Residue sequence number.
 * 22             AChar           icode1       Insertion code. (Ignored)
 * 26 - 28        LString(3)      "CYS"        Residue name.
 * 30             Character       chainID2     Chain identifier.
 * 32 - 35        Integer         seqNum2      Residue sequence number.
 * 36             AChar           icode2       Insertion code. (Ignored)
 * 60 - 65        SymOP           sym1         Symmetry oper for 1st resid (Ignored)
 * 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid (Ignored)
 *
 * For more details visit http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Method parameters:
 * 		pdbData:	pointer to pdbData object to read into
 * 		line: line from the file to parse
 * 		currentSSBond: indicates location in array of ssbonds
 * 			where new bond should be saved
 */
void parseSSBondLine(PDB* pdbData, char* line, int currentSSBond){
	char chain;
	char resid[4];
	strncpy(&chain, &line[15], 1);
	strncpy(resid, &line[17], 4);
	pdbData->ssbonds[currentSSBond].chain1 = chain;
	pdbData->ssbonds[currentSSBond].resid1 = atoi(resid);
	strncpy(&chain, &line[29], 1);
	strncpy(resid, &line[31], 4);
	pdbData->ssbonds[currentSSBond].chain2 = chain;
	pdbData->ssbonds[currentSSBond].resid2 = atoi(resid);
#ifdef DEBUG
	printf("SS #%d: %c%d - %c%d\n", current_ss, pdbData->ssbonds[current_ss].chain1, pdbData->ssbonds[current_ss].resid1,
			pdbData->ssbonds[current_ss].chain2, pdbData->ssbonds[current_ss].resid2);
#endif
}


void parseSymmetryLine(PDB* pdbData, char* line){
	char index[4], r1[9], r2[9], r3[9], t[13];
	index[3] = '\0';
	r1[8] = '\0';
	r2[8] = '\0';
	r3[8] = '\0';
	t[12] = '\0';

	int lineid = atoi(&line[18]);
	//printf("lineid = %d\n", lineid);

	strncpy(index, &line[20], 3);

	strncpy(r1, &line[24], 8);
	strncpy(r2, &line[34], 8);
	strncpy(r3, &line[44], 8);
	strncpy(t, &line[54], 12);
	//printf("'%s'\t'%s'\t'%s'\t'%s'\t'%s'\n", index, r1, r2, r3, t);
	//printf("'%s'\n", line);
	int currentSymmetry = atoi(index) - 1;
	if(lineid == 1){
		pdbData->symmetries[currentSymmetry].index = atoi(index);
		pdbData->symmetries[currentSymmetry].r11 = atof(r1);
		pdbData->symmetries[currentSymmetry].r12 = atof(r2);
		pdbData->symmetries[currentSymmetry].r13 = atof(r3);
		pdbData->symmetries[currentSymmetry].tx = atof(t);
	} else if(lineid == 2){
		pdbData->symmetries[currentSymmetry].index = atoi(index);
		pdbData->symmetries[currentSymmetry].r21 = atof(r1);
		pdbData->symmetries[currentSymmetry].r22 = atof(r2);
		pdbData->symmetries[currentSymmetry].r23 = atof(r3);
		pdbData->symmetries[currentSymmetry].ty = atof(t);
	} else if(lineid == 3){
		pdbData->symmetries[currentSymmetry].index = atoi(index);
		pdbData->symmetries[currentSymmetry].r31 = atof(r1);
		pdbData->symmetries[currentSymmetry].r32 = atof(r2);
		pdbData->symmetries[currentSymmetry].r33 = atof(r3);
		pdbData->symmetries[currentSymmetry].tz = atof(t);
	}
}

void parseMatrixLine(PDB* pdbData, char* line){
	char index[4], Mn1[11], Mn2[11], Mn3[11], Vn[11];
	index[3] = '\0';
	Mn1[10] = '\0';
	Mn2[10] = '\0';
	Mn3[10] = '\0';
	Vn[10] = '\0';

	int lineid = atoi(&line[5]);
	//printf("lineid = %d\n", lineid);

	strncpy(index, &line[7], 3);

	strncpy(Mn1, &line[10], 10);
	strncpy(Mn2, &line[20], 10);
	strncpy(Mn3, &line[30], 10);
	strncpy(Vn, &line[45], 10);
	/*printf("'%s'\t'%s'\t'%s'\t'%s'\t'%s'\n", index, Mn1, Mn2, Mn3, Vn);
	printf("'%s'\n", line);*/
	int currentMatrix = atoi(index) - 1;
	if(lineid == 1){
		pdbData->matrices[currentMatrix].index = atoi(index);
		pdbData->matrices[currentMatrix].M11 = atof(Mn1);
		pdbData->matrices[currentMatrix].M12 = atof(Mn2);
		pdbData->matrices[currentMatrix].M13 = atof(Mn3);
		pdbData->matrices[currentMatrix].V1 = atof(Vn);
	} else if(lineid == 2){
		pdbData->matrices[currentMatrix].index = atoi(index);
		pdbData->matrices[currentMatrix].M21 = atof(Mn1);
		pdbData->matrices[currentMatrix].M22 = atof(Mn2);
		pdbData->matrices[currentMatrix].M23 = atof(Mn3);
		pdbData->matrices[currentMatrix].V2 = atof(Vn);
	} else if(lineid == 3){
		pdbData->matrices[currentMatrix].index = atoi(index);
		pdbData->matrices[currentMatrix].M31 = atof(Mn1);
		pdbData->matrices[currentMatrix].M32 = atof(Mn2);
		pdbData->matrices[currentMatrix].M33 = atof(Mn3);
		pdbData->matrices[currentMatrix].V3 = atof(Vn);
	}
}

void parseCrystalLine(PDB* pdbData, char* line){
	strncpy(pdbData->crystal.line, line, 80);
	char a[9], b[9], c[9], alpha[7], beta[7], gamma[7];
	strncpy(a, &line[6], 9);
	strncpy(b, &line[15], 9);
	strncpy(c, &line[24], 9);
	strncpy(alpha, &line[33], 7);
	strncpy(beta, &line[40], 7);
	strncpy(gamma, &line[47], 7);
	pdbData->crystal.a = atof(a);
	pdbData->crystal.b = atof(b);
	pdbData->crystal.c = atof(c);
	pdbData->crystal.alpha = atof(alpha);
	pdbData->crystal.beta = atof(beta);
	pdbData->crystal.gamma = atof(gamma);
	printf("Crystal: %f\t%f\t%f\t%f\t%f\t%f\n", pdbData->crystal.a, pdbData->crystal.b, pdbData->crystal.c,
			pdbData->crystal.alpha, pdbData->crystal.beta, pdbData->crystal.gamma);
}

/*
 * Saves data from PDB object into a PDB file format
 *
 * For detailed information about PDB file format, please visit
 * http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Parameters:
 * 		filename: name of a file to save into (will be overwritten/created)
 * 		pdbData: pointer of an object to get data from
 */
void writePDB(const char* filename, PDB* pdbData){
	FILE* file = fopen(filename, "w");
	int i;
	//fprintf(file, "%s\n", pdbData->crystal.line);
/*
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000
REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   2  0.000000  1.000000  0.000000       47.80000
REMARK 290   SMTRY3   2  0.000000  0.000000 -1.000000        0.00000
*/
	/*for(i = 0; i < pdbData->symmetryCount; i++){
		PDBSymmetry sym = pdbData->symmetries[i];
		fprintf(file, "REMARK 290   SMTRY%d %3d %9.6f %9.6f %9.6f      %9.5f\n", 1, sym.index, sym.r11, sym.r12, sym.r13, sym.tx);
		fprintf(file, "REMARK 290   SMTRY%d %3d %9.6f %9.6f %9.6f      %9.5f\n", 2, sym.index, sym.r21, sym.r22, sym.r23, sym.ty);
		fprintf(file, "REMARK 290   SMTRY%d %3d %9.6f %9.6f %9.6f      %9.5f\n", 3, sym.index, sym.r31, sym.r32, sym.r33, sym.tz);
	}*/
	/*for(i = 0; i < pdbData->ssCount; i++){
		fprintf(file, "SSBOND %3d CYS %c %4d    CYS %c %4d\n",
								i + 1,
								pdbData->ssbonds[i].chain1,
								pdbData->ssbonds[i].resid1,
								pdbData->ssbonds[i].chain2,
								pdbData->ssbonds[i].resid2);
	}*/
	for(i = 0; i < pdbData->atomCount; i++){
		fprintf(file, "ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %s       \n",
								pdbData->atoms[i].id,
								pdbData->atoms[i].name,
								pdbData->atoms[i].altLoc,
								pdbData->atoms[i].resName,
								pdbData->atoms[i].chain,
								pdbData->atoms[i].resid,
								pdbData->atoms[i].x,
								pdbData->atoms[i].y,
								pdbData->atoms[i].z,
								pdbData->atoms[i].occupancy,
								pdbData->atoms[i].beta,
								pdbData->atoms[i].segment);
	}
	fprintf(file, "END");
	printf("Done writing '%s'.\n", filename);
	fclose(file);
}

/*
 * Prints PDBAtom object in a PDB format.
 * Parameters:
 * 		atomData: PDBAtom object to print
 */
void printAtom(PDBAtom atomData){
	printf("ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      PROT\n",
			atomData.id,
			atomData.name,
			atomData.altLoc,
			atomData.resName,
			atomData.chain,
			atomData.resid,
			atomData.x,
			atomData.y,
			atomData.z,
			atomData.occupancy,
			atomData.beta);
}

/*
 * Prints PDBAtom object in a PDB format into a file.
 * Parameters:
 * 		atomData: PDBAtom object to print
 */
void printAtomToFile(FILE* file, PDBAtom atomData){
	fprintf(file, "ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %c\n",
			atomData.id,
			atomData.name,
			atomData.altLoc,
			atomData.resName,
			' ',
			atomData.resid,
			atomData.x,
			atomData.y,
			atomData.z,
			atomData.occupancy,
			atomData.beta,
			atomData.chain);
}
