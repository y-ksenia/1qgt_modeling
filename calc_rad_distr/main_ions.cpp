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


int ParseDCD(const char* bp, const char* mM, char* dcdFilebp, const char* topfilebp, int stride){	
	readTOP(topfilebp, &top);
	int i;
	double a;
	dcdOpenRead(&dcd, dcdFilebp);

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
		
		sprintf(datNA, "../DNA/dsDNA/%sbp/%smM/push_from_30nm/radius_vectors/Na.%lld.dat", bp, mM, frame);
		sprintf(datCL, "../DNA/dsDNA/%sbp/%smM/push_from_30nm/radius_vectors/Cl.%lld.dat", bp, mM, frame);
			
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
	char* bp = argv[1];
	char* mM = argv[2];
	//long long int max_frame = atoi(argv[2]);
	int stride = atoi(argv[3]);
	char dcdfilebp[1024];
	// sprintf(dcdfilebp, "%s/%s_push.dcd", bp, bp);
	sprintf(dcdfilebp, "../DNA/dsDNA/%sbp/%smM/push_from_30nm/dsDNA_%sbp_%smM_push.dcd", bp, mM, bp, mM);
	
	char topfilebp[1024];
	// sprintf(topfilebp, "%s/%s.top", bp, bp);
	sprintf(topfilebp, "../DNA/dsDNA/%sbp/%smM/topology/dsDNA_%sbp_%smM.top", bp, mM, bp, mM);
	
	ParseDCD(bp, mM, dcdfilebp, topfilebp, stride);
}
