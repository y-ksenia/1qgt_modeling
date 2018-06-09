#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "dcdio.h"

DCD dcdIn;
DCD dcdOut;

typedef struct {
	double x;
	double y;
	double z;
} double3;

int parseDCD(const char* dcdInFilename, const char* dcdOutFilename){
	int i;
	double a;
	double curaver;
	
	dcdOpenRead(&dcdIn, dcdInFilename);
	dcdReadHeader(&dcdIn);
	dcdIn.frame.X = (float*)calloc(dcdIn.header.N, sizeof(float));
	dcdIn.frame.Y = (float*)calloc(dcdIn.header.N, sizeof(float));
	dcdIn.frame.Z = (float*)calloc(dcdIn.header.N, sizeof(float));
	
	createDCD(&dcdOut, dcdIn.header.N, dcdIn.header.nfile, dcdIn.header.npriv, dcdIn.header.delta, dcdIn.header.nsavc, dcdIn.hasUC, dcdIn.uc.a, dcdIn.uc.b, dcdIn.uc.c);
	dcdOpenWrite(&dcdOut, dcdOutFilename);
	dcdWriteHeader(dcdOut);
	dcdOut.frame.X = (float*)calloc(dcdOut.header.N, sizeof(float));
	dcdOut.frame.Y = (float*)calloc(dcdOut.header.N, sizeof(float));
	dcdOut.frame.Z = (float*)calloc(dcdOut.header.N, sizeof(float));

	int frame = 0;
	while(dcdReadFrame(&dcdIn) != -1){
		printf("%d:\t", frame);
		int i;
		for(i = 0; i < dcdIn.header.N; i++){
			dcdOut.frame.X[i] = dcdIn.frame.X[i];
			dcdOut.frame.Y[i] = dcdIn.frame.Y[i];
			dcdOut.frame.Z[i] = dcdIn.frame.Z[i];
		}
		dcdWriteFrame(dcdOut);
		frame ++;
	}
	

}

int main(int argc, char* argv[]){
	//parseDCD("/home/yagafarova/dcd/structures/dna_mM50_push_1.dcd", "/home/yagafarova/dcd/structures/dna_mM50_push_1.corrected.dcd");
	parseDCD("/home/yagafarova/dcd/structures/dna_mM50_push_5.dcd", "/home/yagafarova/dcd/structures/dna_mM50_push_5.corrected.dcd");

}

