#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "xyzio.h"
#include "topio.h"

typedef struct {
	float x;
	float y;
	float z;
} float3; 

int diff_count_CL_NA(float a){
	if(a > 0){
		return (int)a;
	}else{
		return -(int)a;
	}
}

// Функция, генерирующая случайное действительное число от min до max
double random(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}


//На входе - имя входных данных argv[1], имя для выходных данных argv[2], размер коробки argv[3] в А; необходимая концентрация в мМ argv[4]
int main(int argc, char *argv[]){

	char* input_name = argv[1]; //input the name for top and xyz
	char* output_name = argv[2]; 
	double box_size = atoi(argv[3]);
	int mM = atoi(argv[4]);
 	
 	float cutoff_ion = 5;
 	float cutoff = 5;
	srand(time(NULL));

	char inTOP[1024];
	char inXYZ[1024];
	sprintf(inTOP, "structures/%s.top", input_name);
	sprintf(inXYZ, "structures/%s.xyz", input_name);

 	TOPData top;
	readTOP(inTOP, &top);
	XYZ xyz;
	readXYZ(inXYZ, &xyz);

	char outTOP[1024];
	char outXYZ[1024];
	sprintf(outTOP, "structures/%s.top", output_name);
	sprintf(outXYZ, "structures/%s.xyz", output_name);

	FILE* datout = fopen(outTOP, "w");
	FILE* datout2 = fopen(outXYZ, "w");
	
	int countTOP = 0;
	for (int i = 0; i < top.atomCount; i++){
		if(atoi(top.atoms[i].type) == 1){
			countTOP++;
		}
	}

	int ion_count = 2 * (int) (mM * 6.02 * pow(box_size,3) * pow(10,7)) + countTOP;

	float d_p_charge = 0.0;
	for (int i = 0; i < countTOP; i++){
		d_p_charge += top.atoms[i].charge;
	}

	//printf("d_p_charge = %f\n", d_p_charge);

	int count_Cl, count_Na;
	if(int(d_p_charge) % 2 != 0){
		ion_count--;
	}
	if(d_p_charge > 0){
		count_Na = (ion_count - diff_count_CL_NA(d_p_charge))/2;
		count_Cl = (ion_count + diff_count_CL_NA(d_p_charge))/2;
	}else{
		count_Na = (ion_count + diff_count_CL_NA(d_p_charge))/2;
		count_Cl = (ion_count - diff_count_CL_NA(d_p_charge))/2;		
	}


	//XYZ output FILE

	fprintf(datout2, "%d\n", countTOP + ion_count);
	fprintf(datout2, "Created by kir_min\n");
	int count_xyz = 1;
	for (int i = 0; i < countTOP; i++){
		fprintf(datout2, "%s\t%f\t%f\t%f\n", "C", xyz.atoms[i].x, xyz.atoms[i].y, xyz.atoms[i].z);
		count_xyz++;
	}

	float x,y,z;
	float dx, dy, dz, dr;
	int count = 0;
	bool while_flag = true;
	bool flag;
	float3* ions;
	ions = (float3*)calloc(ion_count, sizeof(float3));
	while (while_flag){
		//x = random(-200.0, 200.0);
		//y = random(-200.0, 200.0);
		//z = random(-200.0, 200.0);
		x = random(-box_size/2, box_size/2);
		y = random(-box_size/2, box_size/2);
		z = random(-box_size/2, box_size/2);
		for (int i = 0; i < countTOP; i++){
			dx = xyz.atoms[i].x - x;
			dy = xyz.atoms[i].y - y;
			dz = xyz.atoms[i].z - z;
			dr = sqrt(dx*dx + dy*dy + dz*dz);
			if(dr > cutoff){
				flag = true;
			}else{
				flag = false;
				break;
			}
		}
		if(flag){
			for (int i = 0; i < count; i++){
				dx = ions[i].x - x;
				dy = ions[i].y - y;
				dz = ions[i].z - z;
				dr = sqrt(dx*dx + dy*dy + dz*dz);
				if(dr <= cutoff_ion){
					flag = false;
				}
			}
			if(flag){
				ions[count].x = x;
				ions[count].y = y;
				ions[count].z = z;
				printf("count = %d\n", count);
				count++;
			}
		}
		if(count == ion_count ){//-1
			while_flag = false;
		}
	}
	printf("ion_count_haha = %d\n", ion_count);
	for (int i = 0; i < ion_count; i++){
		fprintf(datout2, "%s\t%f\t%f\t%f\n", "C", ions[i].x, ions[i].y, ions[i].z);
		count_xyz++;
	}
	printf("Na = %d\tCl = %d\n", count_Na, count_Cl);
	printf("Done writting xyz\n");

	
	//TOP output FILE
	int type, resid;
	char resName[10];
	char name[10];
	char chain[1];
	float charge;
	float mass;


	//ATOMS
	fprintf(datout, "Created by kir_min\n");
	fprintf(datout, "\n");
	fprintf(datout, "[ atoms ]\n");
	fprintf(datout, ";\tnr\ttype\tresnr\tresidue\tatom\tcgnr\tcharge\tmass\n");
	int count_top = 1;
	for (int i = 0; i < countTOP; i++){
		fprintf(datout, "\t");
		fprintf(datout, "%d\t%d\t%d\t%s\t%s\t%c\t%.2f\t%.3f\n", count_top, atoi(top.atoms[i].type), top.atoms[i].resid, top.atoms[i].resName, top.atoms[i].name,
			top.atoms[i].chain, top.atoms[i].charge, top.atoms[i].mass);
		count_top++;
		
	}
	for (int i = 0; i < count_Na; i++){
		type = 2;
		resid = 0;
		sprintf(resName, "ION");
		//printf("%s\n", resName);
		strcpy(name, "NA");
		strcpy(chain, "Y");
		charge = 1.0;
		mass = 22.989769;
		fprintf(datout, "\t");
		fprintf(datout, "%d\t%d\t%d\t%s\t%s\t%s\t%.2f\t%.3f\n", count_top, type, resid, "ION", name, chain, charge, mass);
		//printf("%d\t%d\t%d\t'%s'\t%s\t%s\t%.2f\t%.3f\n", count_top, type, resid, "ION", name, chain, charge, mass);
		count_top++;
		
	}
	for (int i = 0; i < count_Cl; i++){
		type = 3;
		resid = 0;
		strcpy(resName, "ION");
		strcpy(name, "CL");
		strcpy(chain, "Z");
		charge = -1.0;
		mass = 35.453;
		fprintf(datout, "\t");
		fprintf(datout, "%d\t%d\t%d\t%s\t%s\t%s\t%.2f\t%.3f\n", count_top, type, resid, "ION", name, chain, charge, mass);
		count_top++;
		
	}
	//BONDS
	fprintf(datout, "\n");
	fprintf(datout, "[ bonds ]\n");
	fprintf(datout, ";\tai\taj\tfunc\tc0\tc1\tc2\tc3\n");

	for (int i = 0; i < top.bondCount; i++){
		fprintf(datout, "\t");
		fprintf(datout, "%d\t%d\t%d\t%f\n", top.bonds[i].i, top.bonds[i].j, top.bonds[i].func, top.bonds[i].c0);
		
	}
	//PAIRS
	/*fprintf(datout, "\n");
	fprintf(datout, "[ pairs ]\n");
	fprintf(datout, ";\tai\taj\tfunc\tc0\tc1\tc2\tc3\n");
	for (int i = 0; i < top.pairsCount; i++){
		fprintf(datout, "\t");
		fprintf(datout, "%d\t%d\t%d\t%f\t%f\n", top.pairs[i].i, top.pairs[i].j, top.pairs[i].func, top.pairs[i].c0, top.pairs[i].c1);
	};*/
	//ANGLE
	fprintf(datout, "\n");
	fprintf(datout, "[ angles ]\n");
	fprintf(datout, ";\tai\taj\tak\tfunc\tc0\tc1\tc2\tc3\n");
	for (int i = 0; i < top.angleCount; i++){
		fprintf(datout, "\t");
		fprintf(datout, "%d\t%d\t%d\t%d\t%.0f\n", top.angles[i].i, top.angles[i].j, top.angles[i].k, top.angles[i].func, top.angles[i].c0);
	}
	//EXCLUSIONS
	fprintf(datout, "\n");
	fprintf(datout, "[ exclusions ]\n");
	fprintf(datout, ";\tai\taj\n");
	for (int i = 0; i < top.exclusionCount; i++){
		fprintf(datout, "\t");
		fprintf(datout, "%d\t%d\t%d\n", top.exclusions[i].i, top.exclusions[i].j, top.exclusions[i].func);
	}
	fclose(datout);
	printf("Done writting top\n");
	fclose(datout2);
 	return 0;
 }
