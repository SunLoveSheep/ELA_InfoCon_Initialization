#include "Data.h"

Data data;

DataOperator::DataOperator(){
}

DataOperator::~DataOperator(){
}

void DataOperator::DataIni()
{
	data.Min=new double[data.D];
	data.Max=new double[data.D];
	data.ObjectiveSequence = new double[data.NumSample];
	data.SampleSequence = new double*[data.NumSample];
	for (int i=0;i<data.NumSample;i++)
	{
		data.SampleSequence[i] = new double[data.D];
		for (int d=0;d<data.D;d++)
		{
			data.SampleSequence[i][d]=0;
		}
	}
}

void DataOperator::DataRelease()
{
	delete []data.Min;
	delete []data.Max;
	delete []data.ObjectiveSequence;
	for (int i=0;i<data.NumSample;i++)
	{
		delete []data.SampleSequence[i];
	}
	delete []data.SampleSequence;
}