#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "dcdio.h"
#include "topio.h"

#define STRIDE 1

DCD dcd;
TOPData top;


int ParseDCD(const char* name, char* dcdFilename, const char* topfilename, int stride){	
	readTOP(topfilename, &top);
	int i;
	double a;
	dcdOpenRead(&dcd, dcdFilename);

	dcdReadHeader(&dcd);
	dcd.frame.X = (float*)calloc(top.atomCount, sizeof(float));
	dcd.frame.Y = (float*)calloc(top.atomCount, sizeof(float));
	dcd.frame.Z = (float*)calloc(top.atomCount, sizeof(float));
	long long int frame = 0;
	int* idx_na = (int*)calloc(top.atomCount, sizeof(int));
	int count_na = 0;
	int* idx_cl = (int*)calloc(top.atomCount, sizeof(int));
	int count_cl = 0;
	for(int i = 0; i < top.atomCount; i++){
		if(atoi(top.atoms[i].type) == 2){
			idx_na[count_na] = i;
			count_na++;
		}
		if(atoi(top.atoms[i].type) == 3){
			idx_cl[count_cl] = i;
			count_cl++;
		}
	}
	

	char datNA[1024];
	char datCL[1024];
	
	
	while(dcdReadFrame(&dcd) != -1){
		
		sprintf(datNA, "../DNA/dsDNA/%sbp/50mM/radius_vectors/Na.%lld.dat", name, frame);
		sprintf(datCL, "../DNA/dsDNA/%sbp/50mM/radius_vectors/Cl.%lld.dat", name, frame);
			
		FILE* out2 = fopen(datNA, "w");
		FILE* out3 = fopen(datCL, "w");

		printf("%lld\n", frame);

		for(int i = 0; i < count_na; i++){
			double x = dcd.frame.X[idx_na[i]];
			double y = dcd.frame.Y[idx_na[i]];
			double z = dcd.frame.Z[idx_na[i]];
	
			double r2 = x*x + y*y + z*z;
			double r = sqrtf(r2);
			fprintf(out2, "%f\n", r);
		}
		for(int i = 0; i < count_cl; i++){
			double x = dcd.frame.X[idx_cl[i]];
			double y = dcd.frame.Y[idx_cl[i]];
			double z = dcd.frame.Z[idx_cl[i]];
	
			double r2 = x*x + y*y + z*z;
			double r = sqrtf(r2);
			fprintf(out3, "%f\n", r);
		}
		
		fclose(out2);
		fclose(out3);
		
		frame ++;
	}
}



int main(int argc, char *argv[]){
	char* name = argv[1];
	//long long int max_frame = atoi(argv[2]);
	int stride = atoi(argv[2]);
	char dcdfilename[1024];
	// sprintf(dcdfilename, "%s/%s_push.dcd", name, name);
	sprintf(dcdfilename, "../DNA/dsDNA/%sbp/50mM/push_from_30nm/dsDNA_%sbp_mM50_push_from_30nm.dcd", name, name);
	
	char topfilename[1024];
	// sprintf(topfilename, "%s/%s.top", name, name);
	sprintf(topfilename, "../DNA/dsDNA/%sbp/50mM/topology/dsDNA_%sbp_mM50.top", name, name);
	
	ParseDCD(name, dcdfilename, topfilename, stride);
}
