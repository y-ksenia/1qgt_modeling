/*
 * config_reader.h
 *
 *  Created on: Feb 27, 2009
 *      Author: zhmurov
 */

extern void parseFile(char* filename);

extern int getParameter(char* paramValue, const char* paramName);
extern int getIntegerParameter(const char* paramName);
extern long long int getLongIntegerParameter(const char* paramName);
extern float getFloatParameter(const char* paramName);
extern int getYesNoParameter(const char* paramName);

extern int getMaskedParameter(char* result, const char* paramName);
