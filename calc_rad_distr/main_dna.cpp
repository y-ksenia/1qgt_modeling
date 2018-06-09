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


int ParseDCD(const char* name, char* dcdFilename, const char* xyzfilename, int stride){

    readXYZ(xyzfilename, &xyz);
    printf("Finish reading xyz");
    int i;
    double a;
    dcdOpenRead(&dcd, dcdFilename);

    dcdReadHeader(&dcd);
    dcd.frame.X = (float*)calloc(xyz.atomCount, sizeof(float));
    dcd.frame.Y = (float*)calloc(xyz.atomCount, sizeof(float));
    dcd.frame.Z = (float*)calloc(xyz.atomCount, sizeof(float));
    printf("Finish reading dcd");
    long long int frame = 0;
    // int* idx_dna = (int*)calloc(xyz.atomCount, sizeof(int));
    int count_dna = xyz.atomCount;


    char datDNA[1024];


    while(dcdReadFrame(&dcd) != -1){
    //for(long long int frame = 0; frame < max_frame; frame = frame + stride){
        if (frame % stride == 0) {

            sprintf(datDNA, "../DNA/dsDNA/%sbp/50mM/radius_vectors/DNA.%lld.dat", name, frame);

            FILE* out1 = fopen(datDNA, "w");
            printf("%lld\n", frame);

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
    char* name = argv[1];
    //long long int max_frame = atoi(argv[2]);
    int stride = atoi(argv[2]);
    char dcdfilename[1024];
    sprintf(dcdfilename, "../DNA/dsDNA/%sbp/50mM/central_line/dsDNA_%sbp_mM50_central_line.dcd", name, name);

    char xyzfilename[1024];
    sprintf(xyzfilename, "../DNA/dsDNA/%sbp/50mM/central_line/dsDNA_%sbp_mM50_central_line.xyz", name, name);

    ParseDCD(name, dcdfilename, xyzfilename, stride);
}
