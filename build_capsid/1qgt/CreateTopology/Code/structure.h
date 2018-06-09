#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

typedef struct {

  int id;
  char   name[5], chain, resName[4], segName[3];
  int    resid;
  double x, y, z;

  double occupancy;
  double beta;

  int ssbond_aminoid; //TO BE REMOVED

} Atom;

typedef struct {

	int resid1;
	char chain1;
	char segment1[3];

	int resid2;
	char chain2;
	char segment2[3];

} SSBond;

typedef struct {
	int atomCount;
	Atom* atoms;
	int ssCount;
	SSBond* ssbonds;
} PDB;

typedef struct {
	int i;
	int j;
	float r0;
} CovalentBond;

typedef struct {
	int i;
	int j;
	float r0;
	float eh;
} NativeContact;

typedef struct {
	int i;
	int j;
} PossiblePair;

typedef struct {
	int aminoCount;
	int bondCount;
	int nativeCount;
	int pairCount;
	Atom* aminos;
	CovalentBond* bonds;
	NativeContact* natives;
	PossiblePair* pairs;
} SOP;

#endif
