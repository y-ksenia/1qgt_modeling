#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "dcdio.h"
#include "xyzio.h"

DCD dcd;
XYZ xyz;


int ParseDCD(const char* bp, const char* mM, int run, char* dcdFilename, const char* xyzfilename, int stride){

    readXYZ(xyzfilename, &xyz);
    printf("Finish reading xyz\n");
    int i;
    double a;
    dcdOpenRead(&dcd, dcdFilename);

    dcdReadHeader(&dcd);
    dcd.frame.X = (float*)calloc(xyz.atomCount, sizeof(float));
    dcd.frame.Y = (float*)calloc(xyz.atomCount, sizeof(float));
    dcd.frame.Z = (float*)calloc(xyz.atomCount, sizeof(float));
    printf("Finish reading dcd\n");
    long long int frame = 0;
    // int* idx_dna = (int*)calloc(xyz.atomCount, sizeof(int));
    int count_dna = xyz.atomCount;


    char datDNA[1024];


    while(dcdReadFrame(&dcd) != -1){
    //for(long long int frame = 0; frame < max_frame; frame = frame + stride){
        if (frame % stride == 0) {

            sprintf(datDNA, "../DNA/dsDNA/%sbp/%smM/push_100k_steps/radius_vectors/run_%d/DNA.%lld.dat", bp, mM, run, frame);
            // printf("../DNA/dsDNA/%sbp/%smM/ppush_100k_steps/radius_vectors/run_%d/DNA.%lld.dat", bp, mM, run, frame);

            FILE* out1 = fopen(datDNA, "w");
            printf("Run_%d\t %lld\n", run, frame);

            for(int i = 0; i < count_dna; i++){
                double x = dcd.frame.X[i];
                double y = dcd.frame.Y[i];
                double z = dcd.frame.Z[i];

                double r2 = x*x + y*y + z*z;
                double r = sqrtf(r2);
                fprintf(out1, "%f\n", r);
            }

            fclose(out1);
        }

        frame ++;

    }
}



int main(int argc, char *argv[]){
    char* bp = argv[1];
    char* mM = argv[2];
    //long long int max_frame = atoi(argv[2]);
    int stride = atoi(argv[3]);
    char dcdfilename[1024];
    char xyzfilename[1024];
    for (int i = 1; i < 6; i++) {
	    sprintf(dcdfilename, "../DNA/dsDNA/%sbp/%smM/push_100k_steps/central_line/run_%d/dsDNA_%sbp_%smM.dcd", bp, mM, i, bp, mM);
    	sprintf(xyzfilename, "../DNA/dsDNA/%sbp/%smM/push_100k_steps/central_line/run_%d/dsDNA_%sbp_%smM.xyz", bp, mM, i, bp, mM);
    	ParseDCD(bp, mM, i, dcdfilename, xyzfilename, stride);
    }
}
