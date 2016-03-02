/*-------------------------------------------------------
Started on 2016.2.27 by Mike Yi SUN

This is a program reading from benchmark functions and making samplings on each functions 
to calculate their information content data. The data will be outputed to .csv files first.
May also add functions to predict function classifications based on Information Content data.

Sampling methods should involve Monte Carlo and Latin Hypercube

5 tracks:
1. Main framework & data structures
2. Functions for sampling
3. Calculating information contents based on samples
4. Output data to files
5. Prediction functions
-------------------------------------------------------*/

#include <iostream>
#include <time.h>
#include <string>
#include "Sampling.h"
#include "Data.h"
#include "Info_Content_Cal.h"

using namespace std;
using std::string;

int main(int __argc, char *__argv[])
{
	srand(time(NULL));

	extern Data data;
	DataOperator dataoperator;
	Sampling sampling;
	Info_Content_Cal info_con_cal;


	data.SamplingMethod = "LHD"; //method of sampling, Latin Hypercube Design (LHD) or Random (RNG)
	data.SequencingMethod = "NN"; //if use LHD, what is the method to sequence samples, Nearest Neighbor (NN) or Random (RNG)
	data.D = 2;
	data.NumSample = 8;
	data.EpsilonMin = 0.00001;
	data.EpsilonMax = 1000000000000000;
	data.NumEpsilon = 1000;
	int StartFunction = 1;//the first function that will be tested
	int MaxFunctionTested = 1;//how many function is not being tested, 1 stands for 1 function.
	int FunctionCount = 0;

	dataoperator.DataIni();

	while(FunctionCount<MaxFunctionTested)
	{
		data.Func_Num = StartFunction+FunctionCount;//assign the current function being tested to data.Func_Num

		//update the variable ranges
		for (int d=0;d<data.D;d++)
		{
			data.Min[d]=-5;
			data.Max[d]=5;
		}

		sampling.SamplingProcess();
		//update objective function values
		for (int i=0;i<data.NumSample;i++)
		{
			data.ObjectiveSequence[i] = i*i - 3*i;
		}
		info_con_cal.Info_Con_Process();

		FunctionCount++;
	}

	dataoperator.DataRelease();

	printf_s("\nprogram end \n");
	getchar();
	return 0;
}
