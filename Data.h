//This is the header file for Data
#include <string>
using std::string;

#ifndef _Data_H
#define _Data_H

struct Data//This structure contains all global variables used for data structuring
{
	int D;//solution dimension
	int NumSample;//Number of samples. Also the number of sections on each dimension.
	string SamplingMethod;//the sampling method to be used
	string SequencingMethod;//the sequencing method to be used

	int Func_Num; //the number of function being tested
	double *Min;
	double *Max;//min and max value of each variable, need to be reassigned for each function

	double EpsilonMax;//The max value of control variable that determines the level of acceptance of ruggidness
	double EpsilonMin;//The min value of control variable that determines the level of acceptance of ruggidness
	int NumEpsilon;//The number of epsilon values we would like to test within the EpsilonMin and EpsilonMax range

	double **SampleSequence;//to record the samples in sequence.
	double *ObjectiveSequence;//to record the objective values of all samples
};

class DataOperator
{
public:
	DataOperator();
	virtual~DataOperator();

	void DataIni();
	void DataRelease();
};

#endif