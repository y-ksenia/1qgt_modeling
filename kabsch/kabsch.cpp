/*
 * kabsch.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: zhmurov
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "xyzio.h"
#include "kabsch.h"

typedef struct {
	double x;
	double y;
	double z;
} Coord;

Coord cmInitial;
Coord cmRef;
Coord vector;

int atomRefCount;
int* initialRef;
int* refRef;

double U1[3][3];
double U2[3][3];

void computeGeometricalCenters(XYZ* initialXYZ, XYZ* refXYZ);
void computeRotationalMatrices(XYZ* initialXYZ, XYZ* refXYZ);

double computeDeterminant(double A[3][3]);
void computeInverseMatrix(double A[3][3], double Ainv[3][3]);
void computeMatrixSqrt(double A[3][3], double Asqrt[3][3]);
void computeTransposeMatrix(double A[3][3], double AT[3][3]);
void computeMatrixProduct(double A[3][3], double B[3][3], double C[3][3]);
void printMatrix(const char* title, double A[3][3]);

double doKabsch(XYZ* refXYZ, XYZ* initialXYZ, XYZ* finalXYZ, int offset, int N){

	int i, j;

	//atomRefCount = 2*N/4;
	atomRefCount = N/4;

	initialRef = (int*)calloc(atomRefCount, sizeof(int));
	refRef = (int*)calloc(atomRefCount, sizeof(int));
	for(i = 0; i < N/4; i++){
		initialRef[i] = N/4 + i;
		refRef[i] = i + offset;		
	}
	/*for(i = 0; i < N/4; i++){
		initialRef[i+N/4] = N/4 + i + N;
		refRef[i+N/4] = offset + i + N;		
	}*/


	// Computing geometrical centers
	printf("Geometrical centers:\n");
	computeGeometricalCenters(initialXYZ, refXYZ);
	printf("Initial: (%5.2f, %5.2f, %5.2f)\n", cmInitial.x, cmInitial.y, cmInitial.z);
	printf("Reference: (%5.2f, %5.2f, %5.2f)\n", cmRef.x, cmRef.y, cmRef.z);
	// Translation vector
	vector.x = cmRef.x - cmInitial.x;
	vector.y = cmRef.y - cmInitial.y;
	vector.z = cmRef.z - cmInitial.z;
	printf("Moving initial protein along (%5.2f, %5.2f, %5.2f).\n",
			vector.x, vector.y, vector.z);

	// Creating new PDB data with shifted coordinates
	for(i = 0; i < initialXYZ->atomCount; i++){
		finalXYZ->atoms[i].x += vector.x;
		finalXYZ->atoms[i].y += vector.y;
		finalXYZ->atoms[i].z += vector.z;

	}
	computeRotationalMatrices(initialXYZ, refXYZ);

	vector.x = cmRef.x - cmInitial.x;
	vector.y = cmRef.y - cmInitial.y;
	vector.z = cmRef.z - cmInitial.z;

	for(i = 0; i < initialXYZ->atomCount; i++){
		double x = finalXYZ->atoms[i].x - cmRef.x;
		double y = finalXYZ->atoms[i].y - cmRef.y;
		double z = finalXYZ->atoms[i].z - cmRef.z;
		finalXYZ->atoms[i].x = U1[0][0]*x + U1[1][0]*y + U1[2][0]*z;
		finalXYZ->atoms[i].y = U1[0][1]*x + U1[1][1]*y + U1[2][1]*z;
		finalXYZ->atoms[i].z = U1[0][2]*x + U1[1][2]*y + U1[2][2]*z;
		finalXYZ->atoms[i].x += cmRef.x;
		finalXYZ->atoms[i].y += cmRef.y;
		finalXYZ->atoms[i].z += cmRef.z;
	}
	double rmsd = 0.0;
	for(i = 0; i < 2*N; i++){
		double dx = finalXYZ->atoms[i].x - refXYZ->atoms[i].x;
		double dy = finalXYZ->atoms[i].y - refXYZ->atoms[i].y;
		double dz = finalXYZ->atoms[i].z - refXYZ->atoms[i].z;
		rmsd += dx*dx + dy*dy + dz*dz;
	}
	rmsd = rmsd/(2*N);
	rmsd = sqrt(rmsd);
	printf("%d: %f\n", offset, rmsd);
	return rmsd;
}

int main(int argc, char *argv[]){
	int offset;
	double minrmsd = 10000000000.0;
	int minoffset;
	int N = 4000;
	double rmsd;
	XYZ refXYZ;
	XYZ initialXYZ;
	XYZ finalXYZ;


	
	
	for(offset = 0; offset < 3*N/4; offset++){
		readXYZ("structures/NH_1kk.10.5nm.xyz", &refXYZ);
		readXYZ("structures/NH_10kk.10.5nm.xyz", &initialXYZ);
		readXYZ("structures/NH_10kk.10.5nm.xyz", &finalXYZ);
		rmsd = doKabsch(&refXYZ, &initialXYZ, &finalXYZ, offset, N);
		if(rmsd < minrmsd){
			minoffset = offset;
			minrmsd = rmsd;
		}
	}
	rmsd = doKabsch(&refXYZ, &initialXYZ, &finalXYZ, minoffset, N);
	writeXYZ("rotated.xyz", &finalXYZ);

}

void computeGeometricalCenters(XYZ* initialXYZ, XYZ* refXYZ){
	int a, i, j;
	cmInitial.x = 0.0;
	cmInitial.y = 0.0;
	cmInitial.z = 0.0;
	cmRef.x = 0.0;
	cmRef.y = 0.0;
	cmRef.z = 0.0;
	for(a = 0; a < atomRefCount; a++){
		i = initialRef[a];
		j = refRef[a];
		XYZAtom atom1 = initialXYZ->atoms[i];
		XYZAtom atom2 = refXYZ->atoms[j];
		cmInitial.x += atom1.x;
		cmInitial.y += atom1.y;
		cmInitial.z += atom1.z;
		cmRef.x += atom2.x;
		cmRef.y += atom2.y;
		cmRef.z += atom2.z;
	}
	cmInitial.x /= atomRefCount;
	cmInitial.y /= atomRefCount;
	cmInitial.z /= atomRefCount;
	cmRef.x /= atomRefCount;
	cmRef.y /= atomRefCount;
	cmRef.z /= atomRefCount;
}

void computeRotationalMatrices(XYZ* initialXYZ, XYZ* refXYZ){
	int a, i, j, k;
	// Computing matrices for Kabsch rotations
	double** coords1 = (double**)malloc(atomRefCount*3*sizeof(double));
	for(i = 0; i < atomRefCount; i++){
		coords1[i] = (double*)calloc(3, sizeof(double));
	}
	double** coords2 = (double**)malloc(atomRefCount*3*sizeof(double));
	for(i = 0; i < atomRefCount; i++){
		coords2[i] = (double*)calloc(3, sizeof(double));
	}
	for(a = 0; a < atomRefCount; a++){
		i = initialRef[a];
		j = refRef[a];
		XYZAtom atom1 = refXYZ->atoms[j];
		XYZAtom atom2 = initialXYZ->atoms[i];
		coords1[a][0] = atom1.x - cmRef.x;
		coords1[a][1] = atom1.y - cmRef.y;
		coords1[a][2] = atom1.z - cmRef.z;
		coords2[a][0] = atom2.x - cmRef.x;
		coords2[a][1] = atom2.y - cmRef.y;
		coords2[a][2] = atom2.z - cmRef.z;
	}
	double C[3][3];
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			C[i][j] = 0;
			for(k = 0; k < atomRefCount; k++){
				C[i][j] += coords1[k][i]*coords2[k][j];
			}
		}
	}

	printMatrix("C", C);
	double det = computeDeterminant(C);
	printf("Determinant: %f\n", det);
	double Cinv[3][3];
	computeInverseMatrix(C, Cinv);
	printMatrix("Cinv", Cinv);
	double I[3][3];
	computeMatrixProduct(Cinv, C, I);
	printMatrix("C*Cinv", I);
	double CT[3][3];
	computeTransposeMatrix(C, CT);
	printMatrix("CT", CT);
	double CTC[3][3];
	computeMatrixProduct(C, CT, CTC);
	printMatrix("CT*C", CTC);
	double CTCsqrt[3][3];
	computeMatrixSqrt(CTC, CTCsqrt);
	printMatrix("(CT*C)^1/2", CTCsqrt);
	double CTCsqrtSqr[3][3];
	computeMatrixProduct(CTCsqrt, CTCsqrt, CTCsqrtSqr);
	printMatrix("((CT*C)^1/2)^2", CTCsqrtSqr);
	computeMatrixProduct(Cinv, CTCsqrt, U1);
	printMatrix("U1", U1);
	double U1T[3][3];
	computeTransposeMatrix(U1, U1T);
	double U1U1T[3][3];
	computeMatrixProduct(U1T, U1, U1U1T);
	printMatrix("U1*U1T", U1U1T);
	det = computeDeterminant(U1);
	printf("det = %f\n", det);
}

double computeDeterminant(double A[3][3]){
	double det;
	det = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2]) - A[1][0]*(A[0][1]*A[2][2] - A[2][1]*A[0][2]) + A[2][0]*(A[0][1]*A[1][2] - A[1][1]*A[0][2]);
	return det;
}
void computeInverseMatrix(double A[3][3], double Ainv[3][3]){
	double det = computeDeterminant(A);
	Ainv[0][0] = (A[1][1]*A[2][2] - A[2][1]*A[1][2])/det;
	Ainv[0][1] = (A[2][1]*A[0][2] - A[0][1]*A[2][2])/det;
	Ainv[0][2] = (A[0][1]*A[1][2] - A[1][1]*A[0][2])/det;

	Ainv[1][0] = (A[2][0]*A[1][2] - A[1][0]*A[2][2])/det;
	Ainv[1][1] = (A[0][0]*A[2][2] - A[2][0]*A[0][2])/det;
	Ainv[1][2] = (A[1][0]*A[0][2] - A[0][0]*A[1][2])/det;

	Ainv[2][0] = (A[1][0]*A[2][1] - A[2][0]*A[1][1])/det;
	Ainv[2][1] = (A[2][0]*A[0][1] - A[0][0]*A[2][1])/det;
	Ainv[2][2] = (A[0][0]*A[1][1] - A[1][0]*A[0][1])/det;

}
void computeMatrixSqrt(double A[3][3], double Asqrt[3][3]){
	double Y[3][3];
	double Z[3][3];
	double Yinv[3][3];
	double Zinv[3][3];
	int i, j, k;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			Y[i][j] = A[i][j];
			if(i == j){
				Z[i][j] = 1.0;
			} else {
				Z[i][j] = 0.0;
			}
		}
	}
	for(k = 0; k < 100; k++){
		computeInverseMatrix(Y, Yinv);
		computeInverseMatrix(Z, Zinv);
		for(i = 0; i < 3; i++){
			for(j = 0; j < 3; j++){
				Y[i][j] = 0.5*(Y[i][j] + Zinv[i][j]);
				Z[i][j] = 0.5*(Z[i][j] + Yinv[i][j]);
			}
		}
	}
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			Asqrt[i][j] = Y[i][j];
		}
	}
}

void computeTransposeMatrix(double A[3][3], double AT[3][3]){
	int i, j;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			AT[i][j] = A[j][i];
		}
	}
}

void computeMatrixProduct(double A[3][3], double B[3][3], double C[3][3]){
	int i, j, k;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			C[i][j] = 0;
			for(k = 0; k < 3; k++){
				C[i][j] += A[k][j] * B[i][k];
			}
		}
	}
}

void printMatrix(const char* title, double A[3][3]){
	printf("%s:\n", title);
	int i, j;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			printf("%f\t", A[i][j]);
		}
		printf("\n");
	}
}



