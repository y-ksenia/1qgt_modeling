/*
 * main.c
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "topio.h"
#include "psfio.h"
#include "config_reader.h"
#include "header.h"

extern void readPDB(char* pdb_filename, PDB* pdbdata);
extern void savePDB(char* filename);

float pairsThreshold;
float rLimitBond;
int covalentRepulsion;
PDB pdbdata;

int checkCovalent(int i, int j);
int checkNative(int i, int j);
int checkPairs(int i, int j);
float getDistance(Atom a1, Atom a2);

int main(int argc, char* argv[]){
	parseFile(argv[1]);
	char pdb_filename[100];
	getMaskedParameter(pdb_filename, "structure");
	readPDB(pdb_filename, &pdbdata);
	pairsThreshold = getFloatParameter("pairs_threshold");
	rLimitBond = getFloatParameter("R_limit_bond");
	covalentRepulsion = getYesNoParameter("covalentRepulsion");
	float eh;
	float eh_inter;
	float eh_intra;
	char ehString[30];
	getParameter(ehString, "eh");
	char tempFilename[100];
	if(strcmp(ehString, "O") == 0){
		printf("Taking eh's values from occupancy column.\n");
		eh = -1.0;
	} else if(strcmp(ehString, "B") == 0){
		printf("Taking eh's values from beta column.\n");
		eh = -2.0;
	} else if(strcmp(ehString, "inter-intra") == 0){
		eh = -3.0;
		eh_inter = getFloatParameter("eh_inter");
		eh_intra = getFloatParameter("eh_intra");
		printf("eh's values are %f for inter-chain interactions,%f for intra-chain.\n",
				eh_inter, eh_intra);
	} else {
		eh = atof(ehString);
		if(eh == 0){
			printf("WARNING: Value of eh in a configuration file is not valid. Should be:\n "
					" - Positive integer, or\n"
					" - 'O' to take values from occupancy column, or\n"
					" - 'B' to take from beta column.\n"
					"eh is set to 0.\n");
		}
	}
	printf("Creating SOP-model...\n");

	sop.aminoCount = 0;
	sop.bondCount = 0;
	sop.nativeCount = 0;
	sop.pairCount = 0;

	int i, j, b;

	for(i = 0; i < pdbdata.atomCount; i++){
		if(strcmp(pdbdata.atoms[i].name, "CA") == 0){
			sop.aminoCount ++;
		}
	}
	printf("Found %d aminos.\n", sop.aminoCount);
	sop.aminos = (Atom*)calloc(sop.aminoCount, sizeof(Atom));
	j = 0;
	for(i = 0; i < pdbdata.atomCount; i++){
		if(strcmp(pdbdata.atoms[i].name, "CA") == 0){
			memcpy(&sop.aminos[j], &pdbdata.atoms[i], sizeof(Atom));
			j++;
		}
	}
/*
	for(i = 0; i < 10; i++){
		printf("%d %s %s\n", i, sop.aminos[i].name, sop.aminos[i].segName);
	}

	exit(0);
*/
	for(i = 0; i < sop.aminoCount; i++){
		for(j = i + 1; j < sop.aminoCount; j++){
			if(checkCovalent(i, j)){
				sop.bondCount ++;
			}
		}
	}
	printf("Found %d covalent bonds.\n", sop.bondCount);
	sop.bonds = (CovalentBond*)calloc(sop.bondCount, sizeof(CovalentBond));
	sop.bondCount = 0;
	for(i = 0; i < sop.aminoCount; i++){
		for(j = i + 1; j < sop.aminoCount; j++){
			if(checkCovalent(i, j)){
				sop.bonds[sop.bondCount].i = i;
				sop.bonds[sop.bondCount].j = j;
				sop.bonds[sop.bondCount].r0 = getDistance(sop.aminos[i], sop.aminos[j]);
				sop.bondCount ++;
			}
		}
	}

	for(i = 0; i < sop.aminoCount; i++){
		for(j = i + 1; j < sop.aminoCount; j++){
			if(checkNative(i, j)){
				sop.nativeCount ++;
			}
		}
	}
	printf("Found %d native contacts.\n", sop.nativeCount);
	sop.natives = (NativeContact*)calloc(sop.nativeCount, sizeof(NativeContact));
	sop.nativeCount = 0;
	getMaskedParameter(tempFilename, "intra");
	FILE* intra = fopen(tempFilename, "w");
	getMaskedParameter(tempFilename, "inter");
	FILE* inter = fopen(tempFilename, "w");
	int interCount = 0;
	int intraCount = 0;
	for(i = 0; i < sop.aminoCount; i++){
		for(j = i + 1; j < sop.aminoCount; j++){
			if(checkNative(i, j)){
				sop.natives[sop.nativeCount].i = i;
				sop.natives[sop.nativeCount].j = j;
				sop.natives[sop.nativeCount].r0 = getDistance(sop.aminos[i], sop.aminos[j]);

				if(eh > 0){
					sop.natives[sop.nativeCount].eh = eh;
				} else
				if(eh == -1.0){
					sop.natives[sop.nativeCount].eh = sqrtf(sop.aminos[i].occupancy*sop.aminos[j].occupancy);
				} else
				if(eh == -2.0){
					sop.natives[sop.nativeCount].eh = sqrtf(sop.aminos[i].beta*sop.aminos[j].beta);
				} else
				if(eh == -3.0){
					//if(abs(i - j) < 262 && sop.aminos[i].chain == sop.aminos[j].chain){
					if( strcmp(sop.aminos[i].segName, sop.aminos[j].segName) == 0 ){
						fprintf(intra, "%f\n", sop.natives[sop.nativeCount].r0);
						sop.natives[sop.nativeCount].eh = eh_intra;
						intraCount ++;
					} else {
						fprintf(inter, "%f\n", sop.natives[sop.nativeCount].r0);
						sop.natives[sop.nativeCount].eh = eh_inter;
						interCount ++;
					}
				}
				sop.nativeCount ++;
			}
		}
	}
	fclose(intra);
	fclose(inter);
	printf("Intra-chain contacts: %d.\n", intraCount);
	printf("Inter-chain contacts: %d.\n", interCount);
/*
	for(i = 0; i < sop.aminoCount; i++){
		for(j = i + 1; j < sop.aminoCount; j++){
			if(checkPairs(i, j)){
				sop.pairCount ++;
			}
		}
	}
	printf("Found %d possible pairs.\n", sop.pairCount);
	sop.pairs = (PossiblePair*)calloc(sop.pairCount, sizeof(PossiblePair));
	sop.pairCount = 0;
	for(i = 0; i < sop.aminoCount; i++){
		for(j = i + 1; j < sop.aminoCount; j++){
			if(checkPairs(i, j)){
				sop.pairs[sop.pairCount].i = i;
				sop.pairs[sop.pairCount].j = j;
				sop.pairCount ++;
			}
		}
	}
*/
	printf("Computing center of mass..\n");
	float rx = 0.0;
	float ry = 0.0;
	float rz = 0.0;
	for(i = 0; i < sop.aminoCount; i++){
		rx += sop.aminos[i].x;
		ry += sop.aminos[i].y;
		rz += sop.aminos[i].z;
	}
	rx /= sop.aminoCount;
	ry /= sop.aminoCount;
	rz /= sop.aminoCount;
	printf("CM coordinates: (%f, %f, %f)\n", rx, ry, rz);

	printf("Saving radial distribution..\n");

	getMaskedParameter(tempFilename, "rdist");
	FILE* rdistFile = fopen(tempFilename, "w");
	float x, y, z, r;
	for(i = 0; i < sop.aminoCount; i++){
		x = sop.aminos[i].x - rx;
		y = sop.aminos[i].y - ry;
		z = sop.aminos[i].z - rz;
		r = sqrtf(x*x + y*y + z*z);
		fprintf(rdistFile, "%f\n", r);
	}
	fclose(rdistFile);

	char top_filename[1024];
	getMaskedParameter(top_filename, "topology");
	saveTOP(top_filename);
	char coord_filename[1024];
	getMaskedParameter(coord_filename, "coordinates");
	savePDB(coord_filename);

	PSF cgpsf;
	cgpsf.natom = sop.aminoCount;
	cgpsf.atoms = (PSFAtom*)calloc(cgpsf.natom, sizeof(PSFAtom));
	cgpsf.nbond = 0;
	cgpsf.nbond = sop.bondCount;
	cgpsf.bonds = (PSFBond*)calloc(cgpsf.nbond, sizeof(PSFBond));
	cgpsf.ncmap = 0;
	cgpsf.nimphi = 0;
	cgpsf.nnb = 0;
	cgpsf.nphi = 0;
	cgpsf.ntheta = 0;


	for(i = 0; i < sop.aminoCount; i++){
		cgpsf.atoms[i].id = i+1;
		sprintf(cgpsf.atoms[i].name, "%s", sop.aminos[i].name);
		sprintf(cgpsf.atoms[i].type, "%s", sop.aminos[i].name);
		sprintf(cgpsf.atoms[i].segment, "%s", sop.aminos[i].segName);

		cgpsf.atoms[i].m = sop.aminos[i].beta;
		cgpsf.atoms[i].q = sop.aminos[i].occupancy;


		cgpsf.atoms[i].resid = sop.aminos[i].resid;
		sprintf(cgpsf.atoms[i].resName, "%s", sop.aminos[i].resName);
	}

	for(b = 0; b < sop.bondCount; b++){
		cgpsf.bonds[b].i = sop.bonds[b].i + 1;
		cgpsf.bonds[b].j = sop.bonds[b].j + 1;
	}

	char psf_filename[1024];
	getMaskedParameter(psf_filename, "psf");
	writePSF(psf_filename, &cgpsf);


	return 0;
}

int checkCovalent(int i, int j){
	if((sop.aminos[i].resid == sop.aminos[j].resid + 1 && sop.aminos[i].chain == sop.aminos[j].chain)
			|| (sop.aminos[j].resid == sop.aminos[i].resid + 1 && sop.aminos[j].chain == sop.aminos[i].chain)){
		if(abs(i - j) < 100){
			return 1;
		}
	}
	// For SS-bonds
	int k;

	if (sop.aminos[i].resid == 61 && sop.aminos[j].resid == 61){
		double dx = sop.aminos[i].x - sop.aminos[j].x;
 		double dy = sop.aminos[i].y - sop.aminos[j].y;
		double dz = sop.aminos[i].z - sop.aminos[j].z;
		double dist = sqrt(dx * dx + dy * dy + dz * dz);
		if (dist < 8){
			return 1; 
		}
	}
	/*
	for(k = 0; k < pdbdata.ssCount; k++){
		if((sop.aminos[i].resid == pdbdata.ssbonds[k].resid1) && (sop.aminos[i].chain == pdbdata.ssbonds[k].chain1) &&
				(sop.aminos[j].resid == pdbdata.ssbonds[k].resid2) && (sop.aminos[j].chain == pdbdata.ssbonds[k].chain2)){
			return 1;
		}
		if((sop.aminos[i].resid == pdbdata.ssbonds[k].resid2) && (sop.aminos[i].chain == pdbdata.ssbonds[k].chain2) &&
						(sop.aminos[j].resid == pdbdata.ssbonds[k].resid1) && (sop.aminos[j].chain == pdbdata.ssbonds[k].chain1)){
			return 1;
		}
	}
	*/
	/*for(k = 0; k < pdbdata.ssCount; k++){
		if((sop.aminos[i].resid == pdbdata.ssbonds[k].resid1) && ( strncmp(sop.aminos[i].segName, pdbdata.ssbonds[k].segment1, 3) == 0) &&
				(sop.aminos[j].resid == pdbdata.ssbonds[k].resid2) && ( strncmp(sop.aminos[j].segName, pdbdata.ssbonds[k].segment2, 3) == 0) ){
			printf("%d CYS %s%d --- %d CYS %s%d\n", i, sop.aminos[i].segName, sop.aminos[i].resid,
													j, sop.aminos[j].segName, sop.aminos[j].resid);
			return 1;
		}
		if((sop.aminos[i].resid == pdbdata.ssbonds[k].resid2) && ( strncmp(sop.aminos[i].segName, pdbdata.ssbonds[k].segment2, 3) == 0) &&
						(sop.aminos[j].resid == pdbdata.ssbonds[k].resid1) && ( strncmp(sop.aminos[j].segName, pdbdata.ssbonds[k].segment1, 3) == 0 )){
			printf("%d CYS %s%d --- %d CYS %s%d\n", i, sop.aminos[i].segName, sop.aminos[i].resid,
													j, sop.aminos[j].segName, sop.aminos[j].resid);
			return 1;
		}
	}
	*/
	return 0;
}

int checkNative(int i, int j){
	if(sop.aminos[i].resid >= 143 || sop.aminos[j].resid >= 143){
		return 0;
	}
	if(!checkCovalent(i, j) && getDistance(sop.aminos[i], sop.aminos[j]) < rLimitBond){
		return 1;
	}
	return 0;
}

int checkPairs(int i, int j){
	if(checkCovalent(i, j) && covalentRepulsion == 0){
		return 0;
	}
	if(checkNative(i, j)){
		return 0;
	}
	if(getDistance(sop.aminos[i], sop.aminos[j]) > pairsThreshold){
		return 0;
	}
	if(i == j) {
		return 0;
	}
	return 1;
}

float getDistance(Atom a1, Atom a2){
	float dx, dy, dz;
	dx = a1.x - a2.x;
	dy = a1.y - a2.y;
	dz = a1.z - a2.z;
	return sqrtf(dx*dx + dy*dy + dz*dz);
}
