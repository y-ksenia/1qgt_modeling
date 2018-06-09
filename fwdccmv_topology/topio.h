/*
 * topio.h
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */
#ifndef _TOPIO_H_
#define _TOPIO_H_

typedef struct {
	int i;
	int j;
	int func;
	float c0;
	float c1;
	float c2;
	float c3;
} TOPPair;


extern void loadTOP(char* filename);
extern void saveTOP(char* filename);

#endif
