/*
 * pdbio.h
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */

#ifndef PDBIO_H_
#define PDBIO_H_
#include <stdio.h>

/*
 * Structures
 */

typedef struct {

  int id;
  char   name[5], chain, resName[4], altLoc;
  int    resid;
  double x, y, z;

  char segment[5];

  double occupancy;
  double beta;

} PDBAtom;

typedef struct {

	int resid1;
	char chain1;

	int resid2;
	char chain2;

} PDBSSBond;
/*
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000
REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   2  0.000000  1.000000  0.000000       47.80000
REMARK 290   SMTRY3   2  0.000000  0.000000 -1.000000        0.00000
*/

typedef struct {
	int index;
	double r11, r12, r13, r21, r22, r23, r31, r32, r33, tx, ty, tz;
} PDBSymmetry;;

typedef struct {
	int index;
	double M11, M12, M13, M21, M22, M23, M31, M32, M33, V1, V2, V3;
} PDBMatrix;
/*
CRYST1   37.530   68.420   47.670  90.00 105.36  90.00 P 1 21 1      2 
*/
typedef struct {
	char line[80];
	double a, b, c, alpha, beta, gamma;
} PDBCrystal;

typedef struct {
	int atomCount;
	PDBAtom* atoms;
	int ssCount;
	PDBSSBond* ssbonds;
	int symmetryCount;
	PDBSymmetry* symmetries;
	int matrixCount;
	PDBMatrix* matrices;
	PDBCrystal crystal;
} PDB;

/*
 * Public methods
 */
void readPDB(const char* filename, PDB* pdbData);
void writePDB(const char* filename, PDB* pdbData);
void printAtom(PDBAtom atomData);
void printAtomToFile(FILE* file, PDBAtom atomData);


#endif /* PDBIO_H_ */
