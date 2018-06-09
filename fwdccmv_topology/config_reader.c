/*
 * config_reader.c
 *
 *  Created on: Jan 23, 2009
 *      Author: zhmurov
 */
#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define buf_size 1024
#define name_length 1024
#define value_length 1024

int paramCount;
char** paramNames;
char** paramValues;

void parseFile(char* filename);

int getParameter(char* paramValue, const char* paramName);
int getIntegerParameter(const char* paramName);
long long int getLongIntegerParameter(const char* paramName);
float getFloatParameter(const char* paramName);
int getYesNoParameter(const char* paramName);

int getMaskedParameter(char* result, const char* paramName);
void applyMask(char* result, const char* parameterMask);
void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);

void parseFile(char* filename){
	FILE* file = fopen(filename, "r");
	if(file != NULL){
		printf("Parsing '%s' parameters file...\n", filename);
	} else {
		printf("ERROR: Parameters file '%s' can not be found. Program will exit.\n", filename);
		exit(0);
	}
	paramCount = 0;
	char buffer[buf_size];
	while(fgets(buffer, buf_size, file) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			//printf("%s\n", buffer);
			paramCount++;
		}
	}
	paramNames = (char**)malloc(name_length*paramCount*sizeof(char));
	paramValues = (char**)malloc(value_length*paramCount*sizeof(char));
	int i;
	for(i = 0; i < paramCount; i++){
		paramNames[i] = (char*)malloc(name_length*sizeof(char));
		paramValues[i] = (char*)malloc(value_length*sizeof(char));
	}
	rewind(file);
	i = 0;
	while(fgets(buffer, buf_size, file) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			char* pch = strtok(buffer, " \t");
			strcpy(paramNames[i], pch);
			pch = strtok(NULL, " \t\n");
			strcpy(paramValues[i], pch);
#ifdef DEBUG
			printf("%s\t%s\n", paramNames[i], paramValues[i]);
#endif
			i++;
		}
	}
	fclose(file);
	printf("Done.\n");
}

int getParameter(char* paramValue, const char* paramName){
	int i;
	for(i = 0; i < paramCount; i++){
		if(strcmp(paramName, paramNames[i]) == 0){
			strcpy(paramValue, paramValues[i]);
#ifdef DEBUG
			printf("'%s' = '%s'\n", paramName, paramValue);
#endif
			return 0;
		}
	}
	printf("WARNING: Parameter '%s' have not been found.\n", paramName);
	strcpy(paramValue, "0");
	return -1;
}

int getIntegerParameter(const char* paramName){
	char paramValue[value_length];
	int error = getMaskedParameter(paramValue, paramName);
	int result = atoi(paramValue);
	if(result == 0 && strcmp(paramValue, "0") != 0){
		printf("WARNING: Wrong value of %s in a configuration file ('%s'). Should be integer. Assigning zero value.\n", paramName, paramValue);
		return 0;
	}
	if(error != 0){
		return 0;
	}
	return result;
}

long long int getLongIntegerParameter(const char* paramName){
	char paramValue[value_length];
	int error = getMaskedParameter(paramValue, paramName);
	long long int result = atol(paramValue);
	if(result == 0 && strcmp(paramValue, "0") != 0){
		printf("WARNING: Wrong value of %s in a configuration file ('%s'). Should be integer. Assigning zero value.\n", paramName, paramValue);
		return 0;
	}
	if(error != 0){
		return 0;
	}
	return result;
}

float getFloatParameter(const char* paramName){
	char paramValue[value_length];
	int error = getMaskedParameter(paramValue, paramName);
	float result = atof(paramValue);
	if(result == 0.0 && strcmp(paramValue, "0") != 0){
		printf("WARNING: Wrong value of %s in a configuration file ('%s'). Should be float. Assigning zero value.\n", paramName, paramValue);
		return 0.0;
	}
	if(error != 0){
		return 0.0;
	}
	return result;
}

int getYesNoParameter(const char* paramName){
	char paramValue[value_length];
	int error = getMaskedParameter(paramValue, paramName);
	if(error != 0){
		return 0;
	}
	if(strcmp(paramValue, "YES") == 0 || strcmp(paramValue, "Yes") == 0 || strcmp(paramValue, "yes") == 0
			|| strcmp(paramValue, "Y") == 0 || strcmp(paramValue, "y") == 0
			|| strcmp(paramValue, "ON") == 0 || strcmp(paramValue, "On") == 0 || strcmp(paramValue, "on") == 0
			|| strcmp(paramValue, "TRUE") == 0 || strcmp(paramValue, "True") == 0 || strcmp(paramValue, "true") == 0){
		return 1;
	} else {
		return 0;
	}
}

int getMaskedParameter(char* result, const char* paramName){
	char parameterMask[value_length];
	int error = getParameter(parameterMask, paramName);
	if(error == 0){
		applyMask(result, parameterMask);
	}
	return error;
}

void applyMask(char* result, const char* parameterMask){
	char tempstring[1024];
	strcpy(tempstring, parameterMask);
	int i;
	for(i = 0; i < paramCount; i++){
		char paramName[name_length];
		sprintf(paramName, "<%s>", paramNames[i]);
		replaceString(result, tempstring, paramValues[i], paramName);
		strcpy(tempstring, result);
	}
}

void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace){
	//printf("Looking for %s in %s.\n", stringToReplace, initialString);
	int i;
	int len1, len2;
	for(i = 0; i < strlen(resultString); i++){
		resultString[i] = 0;
	}
	//printf("Result string: %s\n", resultString);
	if(strstr(initialString, stringToReplace) != NULL){
		len1 = strlen(initialString) - strlen(strstr(initialString, stringToReplace));
		//printf("len1: %d\n", len1);
		strncpy(resultString, initialString, len1);
		//printf("Result string: %s\n", resultString);
		strncpy(&resultString[len1], replacementString, strlen(replacementString));
		len2 = len1 + strlen(stringToReplace);
		len1 += strlen(replacementString);
		//printf("Result string: %s\n", resultString);

		strcpy(&resultString[len1], &initialString[len2]);
		strcat(resultString, "\0");
		//printf("Found. Result: %s\n", resultString);
	} else {
		strcpy(resultString, initialString);
	}
}
