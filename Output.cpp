/*
This file output results to .csv files. 
It will output to two .csv files, one contains the HMax, Epsilon_s, M0, and Epsilon_0.5 data of all functions
The other contains the H and M values for all epsilon values and all functions, in order to plot curves.
*/
#include <iostream>
#include <fstream>
#include "Output.h"
#include "Data.h"

using namespace std;

extern Data data;

Output::Output()
{}

Output::~Output()
{}

void Output::ResultOutput(int funcnum)
{
	
	FILE* FourFeatureOutput;
	char FinalFileName[100];
	sprintf_s(FinalFileName,"ResultOutput//FourFeatureValues.csv");
	//getchar();
	fstream finalfoutclear(FinalFileName,ios::out);
	finalfoutclear.close();
	errno_t err=fopen_s(&FourFeatureOutput,FinalFileName,"a+");
	if (err==0)
	{
		fprintf(FourFeatureOutput,", Hmax, Epsilon s, M0, Epsilon 0.5 \n");
		for (int i=0;i<data.MaxFuncTested;i++)
		{
			fprintf(FourFeatureOutput, "F%d, %f, %f, %f, %f \n",data.StartFunction+i, data.HMaxVector[i], data.Epsilon_sVector[i], data.M0Vector[i], data.Epsilon_05Vector[i]);
		}
	}
	fclose(FourFeatureOutput);
}
